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

void TWorld::calcNewNodalValues(long i_, double *Hnew, double *K, double *C1)
{
}


void TWorld::solveFiniteElement(long i_,double *Hnew, double *K, double *C1)
{


}


void TWorld::cell_Soilwater(long i_)
{
    SOIL_LIST s;
    s = crSoil[i_];

    double dtmin = 0.01;
    bool stopit = false;
    int c = s.c;
    int r = s.r;

    bool freeDrainage = !SwitchImpermeable;

    s.drain = 0;

    double ANE, FNN1, FNN2, A1, A2;

    double *Hold = new double [nNodes];
    double *Hnew = new double [nNodes];

    double *C1 = new double [nNodes];
    double *S  = new double [nNodes];
    double *K  = new double [nNodes];

    double *K1 = new double[nNodes];
    double *K2 = new double[nNodes];
    double *A  = new double[nNodes];
    double *D  = new double[nNodes];
    double *F  = new double[nNodes];

    double WH0;
    double cnt = 0;
    int nN = nNodes-1;

    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc; // flood in kin wave

    s.InfPot = WH0/_dt;// m/s
    s.Infact = 0;
    s.dtsum = 0;

    #pragma omp parallel for num_threads(userCores)
    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.h[j];
        S[j] = 0;
    }
        //calcSinkterm(i_); // ETa
    // outer loop timestep lisem
    do {


        int NIT = 0;
        int NITMAX = 12;

        // iteration per gridcell
        do {
            s.ponded = Hnew[0] >= 0;

            #pragma omp parallel for num_threads(userCores)
            for(int j = 0; j < nNodes; j++)
                Hold[j] = Hnew[j];

            // Brooks Corey for H,K,theta and C
            //calcNewNodalValues(i_, Hnew, K, C1);
            // Galerkin 3-diagonal scheme and back substitution
            //solveFiniteElement(i_, Hnew, K, C1);

            // values of K, theta and C for H
            #pragma omp parallel for num_threads(userCores)
            for(int j = 0; j < nNodes; j++) {
                double Hx = 0.5*(Hnew[j] + s.h[j]);

                if (Hx < s.hb[j])
                    K[j] = s.Ks[j]*pow(s.hb[j]/Hx, 2.0+3.0*s.lambda[j]);
                else
                    K[j] = s.Ks[j];

                double Wnew;
                if (Hnew[j] < s.hb[j])
                    Wnew = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hnew[j], s.lambda[j]);
                else
                    Wnew = s.pore[j];

                if (fabs(Hnew[j]-s.h[j]) <= 3*tol2) {
                    if (Hx < s.hb[j])
                        C1[j] = s.theta[j] * 1e-6/s.pore[j] -1.0/Hx *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(s.hb[j]/Hx, s.lambda[j]) * s.dz[j]/s.dts;
                    else
                        C1[j] = s.theta[j] * 1e-6/s.pore[j];
                } else
                    C1[j] = (Wnew-s.theta[j])/(Hnew[j]-s.h[j]) * s.dz[j]/s.dts;
            }

            // K1 and K2 are avg between node and upper and lower node
            #pragma omp parallel for num_threads(userCores)
            for(int j = 1; j < nNodes; j++)
                K1[j] = Aavg(K[j-1],K[j]);

            #pragma omp parallel for num_threads(userCores)
            for(int j = 0; j < nNodes-1; j++)
                K2[j] = Aavg(K[j],K[j+1]);
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
            ANE = A[nN-1];
            A[nN] = -K2[nN];
            F[nN] = (C1[nN-1] + 2.0*C1[nN])/6.0;
            D[nN] = -A[nN-1] + F[nN];

          //  if (i_ == 3738) qDebug() << "a" << Hnew[0] << F[0] << A[0] << D[0];

            FNN1 = F[nN];
            // include gravity term (i.e. k only) and sink term
            F[0] = F[0]*s.h[0] - s.dz[0]*(K1[0] - (2*S[0]+S[1])/6);
            for(int j = 1; j < nNodes-1; j++) {
                F[j] = F[j]*s.h[j]
                       + 0.5*(s.dz[j]+s.dz[j-1]) * K1[j]
                       - 0.5*(s.dz[j]+s.dz[j+1]) * K2[j]
                       - s.dz[j] * (S[j-1]+4*S[j]+S[j+1])/6;
            }
            F[nN] = F[nN]*s.h[nN] + s.dz[nN] * (K2[nN] - (S[nN-1]+2*S[nN])/6);

            //if (i_ == 3738) qDebug() << "b" << Hnew[0] << F[0] << A[0] << D[0];

            // upper boundary condition
            if (s.ponded){
                Hnew[0] = s.h[0] + (s.InfPot-s.Infact)*s.dts;
                s.Infact = std::min(s.InfPot, D[0]*Hnew[0] - F[0]);
                F[0] = Hnew[0];
                A1 = A[0];
                A[0] = 0;
                D[0] = 1;
                F[1] = F[1] - A1 * Hnew[0];
            } else {
                s.Infact = s.InfPot;
                F[0] = F[0] + s.Infact; // direct adding, units??
            }

            // lower boundary condition
            if (freeDrainage) {
                D[nN] = -A[nN-1];
                FNN2 = F[nN];
                F[nN] = 0;
            } else {
                // fixed flux (0)
                s.drain = 0;
                F[nN] = F[nN] - s.drain;
            }

            // calc Hnew with Gaussian elimination and back substitution
            for(int j = 1; j < nNodes; j++) {
                A2 = A[j-1]/D[j-1];
                D[j] = D[j] - A2 * A[j-1];
                F[j] = F[j] - A2 * F[j-1];
            }
            Hnew[nN] = F[nN]/D[nN];

            //if (i_ == 3738) qDebug() << nN << Hnew[nN] << F[nN] << A[nN] << D[nN];
            for (int j = nNodes-2; j >= 0; j--) {
                Hnew[j] = (F[j] - A[j]*Hnew[j+1])/D[j];
               // if (i_ == 3738) qDebug() << j << Hnew[j] << F[j] << A[j] << D[j];
            }

            //if (i_ == 3738) qDebug() << "c" << Hnew[0] << Hnew[1] << F[0] << A[0] << D[0];

          //  #pragma omp parallel for num_threads(userCores)
            //for(int j = 1; j < nNodes; j++)
              //  Hnew[j] = std::min(0.0, Hnew[j]);
            // node 0 can be > 0 ?

            // calc boundary fluxes

            if (s.ponded && Hnew[1] < 0) {
                s.Infact = s.Infact + A1 * Hnew[1];
            }

            if (freeDrainage)
                s.drain = 0.5*(K[nN] + FNN2 - std::min(Hnew[nN],0.0)*FNN1);

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

        #pragma omp parallel for num_threads(userCores)
        for(int j = 0; j < nNodes; j++) {
            s.h[j] = Hnew[j];
            if (s.h[j] < s.hb[j])
                s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
            else
                s.theta[j] = s.pore[j];
        }

        double factor = 0.5 + 1/sqrt((double)NIT);
        s.dts = s.dts * factor;
        s.dts = std::max(s.dts,dtmin);
        s.dts = std::min(s.dts,_dt-s.dtsum);

        s.dtsum += s.dts;
        cnt += 1.0;
    } while(s.dtsum < _dt);

    s.dts = _dt/(cnt);

    if (FloodDomain->Drc == 0) {
        s.Infact = std::min(s.Infact, WH->Drc);
        WH->Drc -= s.Infact; //runoff in kinwave or dyn wave
    } else {
        s.Infact = std::min(s.Infact, hmx->Drc);
        hmx->Drc -= s.Infact; // flood in kin wave
    }

    InfilVol->Drc = s.Infact*SoilWidthDX->Drc*DX->Drc;

    ThetaI1a->Drc = s.theta[0];
    for (int j = 1; j < nN1_; j++)
        ThetaI1a->Drc += s.theta[j];
    ThetaI1a->Drc /= (double)nN1_;

    ThetaI2a->Drc = s.theta[nN1_];
    for (int j = nN1_+1; j < nNodes; j++)
        ThetaI2a->Drc += s.theta[j];
    ThetaI2a->Drc /= (double)nN2_;

    Perc->Drc = s.drain*_dt;

    double gw = 0;
    for (int j = nN; j > 0; j--) {
        if (s.h[j] >= 0)
            gw = gw + s.dz[j];
    }
    Lw->Drc = gw;
//    if (i_ == 3738) qDebug() << cnt << s.dts << s.drain << gw << Hnew[nN];
    if (i_==3738) {
        QString S;
        for(int j = 0; j < nNodes; j++) {
            S = S + QString(" %1").arg(s.h[j]);
        }
        //qDebug() << s.Infact << S;
    }

    delete(Hnew);
    delete(Hold);
    delete(C1);
    delete(S);
    delete(K);
    delete(K1);
    delete(K2);
    delete(A );
    delete(D );
    delete(F );
}



