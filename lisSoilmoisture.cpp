/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#include "lisemqt.h"
#include "model.h"

#define Aavg(a,b)  (0.5*(a+b))
#define Havg(a,b)  (2.0/(1.0/a+1.0/b))
#define Savg(a,b)  sqrt(a * b)
#define tol2 0.2
#define tol1 0.01
/*
void TWorld::calcNewNodalValues(long i_, double Hnew[], double Wnew[], double C1[], double C2[])
{
    SOIL_LIST s;
    s = crSoil[i_];

    for(int j = 0; j < nNodes; j++) {
        double Hx = 0.5*(Hnew[j] + s.h[j]);

        //check is H negative or not!, note hb is negative
        if (Hx < s.hb[j])
            s.K[j] = s.Ks[j]*pow(s.hb[j]/Hx, 2.0+3.0*s.lambda[j]);
        else
            s.K[j] = s.Ks[j];

        if (Hnew[j] < s.hb[j]) {
            Wnew[j] = s.thetar[j] + pow(s.hb[j]/Hnew[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
        } else {
            Wnew[j] = s.pore[j];
        }

        if (j < nNodes-1) {

            // C = dtheta/dh and theta_e(h) = (hb/h)^lambda
            // https://www.derivative-calculator.net/
            // a+b(c/x)^d => -[bd(c/x)^d]/x

            if (fabs(Hnew[j]-s.h[j]) <= 3*tol2)
                C1[j] = -1.0/Hx *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(s.hb[j]/Hx, s.lambda[j]) * s.dz[j]/s.dts;
            else
                C1[j] = (Wnew[j]-s.theta[j])/(Hnew[j]-s.h[j]) * s.dz[j]/s.dts;

            C2[j] = 0;

       }
    }
}


void TWorld::solveFiniteElement(long i_,double *influx, double Hnew[], double C1[], double C2[])
{
    SOIL_LIST s;
    s = crSoil[i_];
    double ANE, FNN1, FNN2, A1, A2;

   // double H[nNodes];
    double K1[nNodes];
    double K2[nNodes];
    double A[nNodes];
    double D[nNodes];
    double F[nNodes];
    double S[nNodes];
    int nN = nNodes-1;

    for(int j = 0; j < nNodes; j++) {
        S[j] = 0;
    }

    // K1 and K2 are avg between node and upper and lower node
    for(int j = 1; j < nNodes; j++)
        K1[j] = 0.5*(s.K[j-1]+s.K[j]);
    for(int j = 0; j < nNodes-1; j++)
        K2[j] = 0.5*(s.K[j]+s.K[j+1]);
    K1[0] = K2[0];
    K2[nN] = K1[nN];

    // sort of Simpson trapezium solution
    F[0] = (2.0*C1[0] + C1[1])/6.0;
    A[0] = -K1[0];
    D[0] = -A[0]+F[0];
    for(int j = 1; j < nNodes; j++) {
        F[j] = (C1[j-1] + 4.0*C1[j] + C1[j+1])/6.0;
        A[j] = -K2[j];
        D[j] = K1[j] + K2[j] + F[j];
    }
    ANE = A[nN-1];
    F[nN] = (C1[nN-1] + 2.0*C1[nN])/6.0;
    D[nN] = -A[nN-1] + F[nN];

    FNN1 = F[nN];
    F[0] = F[0]*s.h[0] - s.dz[0]*(K1[0] - (2*S[0]+S[1])/6);
    for(int j = 1; j < nNodes-1; j++) {
        F[j] = F[j]*s.h[j]
               + 0.5*(s.dz[j]+s.dz[j-1]) * K1[j]
               - 0.5*(s.dz[j]+s.dz[j+1]) * K2[j]
               - s.dz[j] * (S[j-1]+4*S[j]+S[j+1])/6;
    }
    F[nN] = F[nN]*s.h[nN] + s.dz[nN] * (K2[nN] - (S[nN-1]+2*S[nN])/6);

    // upper boundary condition
    if (s.ponded){
        Hnew[0] = s.h[0] + *influx*s.dts;
        *influx = D[0]*Hnew[0] - F[0];
        F[0] = Hnew[0];
        A1 = A[0];
        A[0] = 0;
        D[0] = 1;
        F[1] = F[1] - A1 * Hnew[0];
    } else {
       F[0] = F[0] + *influx;
    }

    // lower boundary condition
    D[nN] = -A[nN-1];
    FNN2 = F[nN];
    F[nN] = 0;

 // calc Hnew with Gaussian elimination and back substitution
    for(int j = 1; j < nNodes; j++) {
        A2 = A[j-1]/D[j-1];
        D[j] = D[j] - A2 * A[j-1];
        F[j] = F[j] - A2 * F[j-1];
    }
    Hnew[nN] = F[nN]/D[nN];
    for (int j = nNodes-2; j >= 0; j--) {
        Hnew[j] = (F[j] - A[j]*Hnew[j+1])/D[j];
    }

//    for(int j = 1; j < nNodes; j++) {
//        //if (lowerBC = 2) and (z[i] >= GWlevel) then break;
//        Hnew[j] = std::min(Hnew[j],0.0);//s.hb[j]);
//    }

  // calc boundary fluxes

    if (s.ponded && Hnew[1] < 0) {
        *influx = *influx + A1 * Hnew[1];
        // de influx is nu gelijk aan:
        //K1*((hnew1-hnew2)/dz+1) + 0.5*dé/dh * (hnew1-h1) }
    }

    s.drain =0.5*(s.Ks[nN]*pow(s.hb[nN]/s.h[nN], 2.0+3.0*s.lambda[nN]) +
                  s.Ks[nN]*pow(s.hb[nN]/Hnew[nN], 2.0+3.0*s.lambda[nN]));

}

*/
void TWorld::cell_Soilwater(long i_)
{
    SOIL_LIST s;
    s = crSoil[i_];

    double dtmin = 0.01;
    bool stopit = false;
    double Hold[nNodes];
    double Hnew[nNodes];
    double Wnew[nNodes];
    double C1[nNodes];
   // double C2[nNodes];
    double WH0;

    int c = s.c;
    int r = s.r;
    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc; // flood in kin wave
    //WH0 = 0.1;

    double influx = WH0/_dt;// m/s
    double Infact = 0;

    s.dtsum = 0;
   // s.dts = _dt/10;


//    if (WH0 > 0) {
//        Hnew[0] = WH0;
//        s.h[0] = WH0;
//        Wnew[0] = s.pore[0];
//    }

//    double moist = s.dz[1]*(s.theta[1]-s.thetar[1]);
//    moist = moist + WH0;
//    if (moist/s.dz[1]+s.thetar[1] > s.pore[1]) {
//        s.ponded = true;
//        s.theta[0] = s.pore[0];
//        s.theta[1] = s.pore[1];
//        s.hn[0] = 0;
//        s.hn[1] = 0;
//    } else {
//        s.theta[1] = s.thetar[1] + moist/s.dz[1];
//        s.hn[1] = s.hb[1]/pow((s.theta[1]-s.thetar[1])/(s.pore[1]-s.thetar[1]), 1.0/s.lambda[1]);
//        s.hn[1] = s.hn[0];
//        s.theta[0] = s.theta[1];
//        influx = 0;
//    }

    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.hn[j];
        //Wnew[j] = s.theta[j];
    }
    double cnt = 0;
    // outer loop timestep lisem
    do {
        //calcSinkterm(i_); // ETa

        int NIT = 0;
        int NITMAX = 12;

        // iteration per gridcell
        do {
            s.ponded = Hnew[0] >= 0;

            for(int j = 0; j < nNodes; j++)
                Hold[j] = Hnew[j];

            // Brooks Corey for H,K,theta and C
            //calcNewNodalValues(i_, Hnew, Wnew, C1, C2);
            // Galerkin 3-diagonal scheme and back substitution
            //solveFiniteElement(i_, &influx, Hnew, C1, C2);

            // values of K, theta and C for H
            for(int j = 0; j < nNodes; j++) {
                double Hx = 0.5*(Hnew[j] + s.h[j]);

                if (Hx < s.hb[j])
                    s.K[j] = s.Ks[j]*pow(s.hb[j]/Hx, 2.0+3.0*s.lambda[j]);
                else
                    s.K[j] = s.Ks[j];

                if (Hnew[j] < s.hb[j])
                    Wnew[j] = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hnew[j], s.lambda[j]);
                else
                    Wnew[j] = s.pore[j];

                if (fabs(Hnew[j]-s.h[j]) <= 3*tol2)
                    C1[j] = -1.0/Hx *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(s.hb[j]/Hx, s.lambda[j]) * s.dz[j]/s.dts;
                else
                    C1[j] = (Wnew[j]-s.theta[j])/(Hnew[j]-s.h[j]) * s.dz[j]/s.dts;

            }

            double ANE, FNN1, FNN2, A1, A2;
            double K1[nNodes];
            double K2[nNodes];
            double A[nNodes];
            double D[nNodes];
            double F[nNodes];
            double S[nNodes];
            int nN = nNodes-1;

            for(int j = 0; j < nNodes; j++) {
                S[j] = 0;
            }

            // K1 and K2 are avg between node and upper and lower node
            for(int j = 1; j < nNodes; j++)
                K1[j] = Aavg(s.K[j-1],s.K[j]);
            for(int j = 0; j < nNodes-1; j++)
                K2[j] = Aavg(s.K[j],s.K[j+1]);
            K1[0] = K2[0];
            K2[nN] = K1[nN];

            // sort of Simpson trapezium solution
            F[0] = (2.0*C1[0] + C1[1])/6.0;
            A[0] = -K1[0];
            D[0] = -A[0]+F[0];
            for(int j = 1; j < nNodes-1; j++) {
                F[j] = (C1[j-1] + 4.0*C1[j] + C1[j+1])/6.0; // dtheta/dh * dz/dt
                A[j] = -K2[j];
                D[j] = K1[j] + K2[j] + F[j];
            }
            //free drainabe bc
            ANE = A[nN-1];
            F[nN] = (C1[nN-1] + 2.0*C1[nN])/6.0;
            D[nN] = -A[nN-1] + F[nN];

            FNN1 = F[nN];
            F[0] = F[0]*s.h[0] - s.dz[0]*(K1[0] - (2*S[0]+S[1])/6);
            for(int j = 1; j < nNodes-1; j++) {
                F[j] = F[j]*s.h[j]
                       + 0.5*(s.dz[j]+s.dz[j-1]) * K1[j]
                       - 0.5*(s.dz[j]+s.dz[j+1]) * K2[j]
                       - s.dz[j] * (S[j-1]+4*S[j]+S[j+1])/6;
            }
            F[nN] = F[nN]*s.h[nN] + s.dz[nN] * (K2[nN] - (S[nN-1]+2*S[nN])/6);

            // upper boundary condition
            if (s.ponded){
                Hnew[0] = s.h[0] + (influx-Infact)*s.dts;
                Infact = D[0]*Hnew[0] - F[0];
                F[0] = Hnew[0];
                A1 = A[0];
                A[0] = 0;
                D[0] = 1;
                F[1] = F[1] - A1 * Hnew[0];
            } else {
                Infact = influx;
                F[0] = F[0] + Infact; // direct adding, units??
            }

            // lower boundary condition
            D[nN] = -A[nN-1];
            FNN2 = F[nN];
            F[nN] = 0;

            // calc Hnew with Gaussian elimination and back substitution
            for(int j = 1; j < nNodes; j++) {
                A2 = A[j-1]/D[j-1];
                D[j] = D[j] - A2 * A[j-1];
                F[j] = F[j] - A2 * F[j-1];
            }
            Hnew[nN] = F[nN]/D[nN];
            for (int j = nNodes-2; j >= 0; j--) {
                Hnew[j] = (F[j] - A[j]*Hnew[j+1])/D[j];
            }

            for(int j = 1; j < nNodes; j++)
                Hnew[j] = std::max(0.0, Hnew[j]);
            // node o can be > 0 ?

            // calc boundary fluxes
            if (s.ponded && Hnew[1] < 0) {
                Infact = Infact + A1 * Hnew[1];
                // de influx is nu gelijk aan:
                //K1*((hnew1-hnew2)/dz+1) + 0.5*dé/dh * (hnew1-h1) }
            }

            double Drain =0.5*(s.Ks[nN]*pow(s.hb[nN]/s.h[nN], 2.0+3.0*s.lambda[nN]) +
                               s.Ks[nN]*pow(s.hb[nN]/Hnew[nN], 2.0+3.0*s.lambda[nN]));
            if (Hnew[nN] < 0)
                s.drain = 0.5*(Drain + FNN2 - Hnew[nN]*FNN1);
            else
                s.drain = 0.5*(Drain + FNN2);

            s.ponded = Hnew[0] >= 0;

            stopit = true;
            double tol = 0;
            for(int j = 0; j < nNodes; j++) {
                tol = tol1*fabs(Hold[j]) + tol2;
                if (s.ponded)
                    tol /= 2.0;
                if (fabs(Hnew[j] - Hold[j]) > tol) {
                    stopit = false;
                    break;
                }
            }

            NIT++;

            if (!stopit && NIT > NITMAX) {
                // restore and do again
                s.dts /= 2.0;
                s.dtsum -= s.dts;
                for(int j = 0; j < nNodes; j++)
                    Hnew[j] = 0.5*(s.h[j]+Hnew[j]);
                stopit = true;
            }

        } while(!stopit);
        // end SoilMoisture in PAS code

        //next: if not NextRainPulse then AdjustTime;
        for(int j = 0; j < nNodes; j++) {
            s.h[j] = Hnew[j];
            s.hn[j] = Hnew[j];
            if (s.h[j] < s.hb[j])
                s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
            else
                s.theta[j] = s.pore[j];
            Wnew[j] = s.theta[j]; // not needed?
        }


        double dtlast = s.dts;
        double factor = 0.5 + 1/sqrt((double)NIT);
        s.dts = s.dts * factor;
        s.dts = std::max(s.dts,dtmin);
        s.dts = std::min(s.dts,_dt-s.dtsum);
        factor = s.dts/dtlast;

        s.dtsum += s.dts;
        cnt += 1.0;

//        for(int j = 0; j < nNodes; j++) {
//            double dH = (Hnew[j] - s.h[j])*factor;
//            s.h[j] = Hnew[j];
//            Hnew[j] = std::min(0.0, Hnew[j] + dH); // ???????
//            s.hn[j] = Hnew[j];

//            if (s.h[j] < s.hb[j])
//                s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
//            else
//                s.theta[j] = s.pore[j];
//        }

        if (i_==3738) {
            QString S;
            for(int j = 0; j < nNodes; j++) {
                S = S + QString(" %1").arg(Hnew[j]);
            }
            qDebug() << NIT << s.dtsum << Hnew[0] << Hnew[1] << influx;
        }

    } while(s.dtsum < _dt);

    s.dts = _dt/cnt;

    if (FloodDomain->Drc == 0) {
        Infact = std::min(Infact, WH->Drc);
        WH->Drc -= Infact; //runoff in kinwave or dyn wave
    } else {
        Infact = std::min(Infact, hmx->Drc);
        hmx->Drc -= Infact; // flood in kin wave
    }

    InfilVol->Drc = Infact*SoilWidthDX->Drc*DX->Drc;

    ThetaI1a->Drc = (s.theta[0]);// + s.theta[1] + s.theta[2] + s.theta[3])/4;

//     if (i_ == 1000)
//        qDebug() << (time - BeginTime)/60 << dInf;
//??????????
//    if (s.dts < dtmin)
//        s.dts = dtorg;

    // calcCorrections(i_);

}



