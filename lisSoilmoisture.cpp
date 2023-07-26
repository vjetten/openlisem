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


void TWorld::solveFiniteElement(long i_,double influx, double *Hnew, , double *Wnew, double *Cnew)
{
    SOIL_LIST s;
    s = crSoil[i_];
    double ANE, FNN1, FNN2, A1, A2;

    double H[nNodes];
    double K1[nNodes];
    double K2[nNodes];
    double A[nNodes];
    double D[nNodes];
    double F[nNodes];
    double S[nNodes];

    for(int j = 0; j < nNodes; j++) {
        S[j] = 0;
        H[j] = s.h[j];
    }

    for(int j = 1; j < nNodes; j++)
        K1[j] = 0.5*(s.K[j-1]+s.K[j])/s.dz[j];
    for(int j = 0; j < nNodes-1; j++)
        K2[j] = 0.5*(s.K[j]+s.K[j+1])/s.dz[j];
    K1[0] = K2[0];
    K2[nNodes] = K1[nNodes];

    F[0] = (2.0*C1[0] + C2[0])/6.0;
    A[0] = -K1[0];
    D[0] = -A[0]+F[0];
    for(int j = 1; j < nNodes; j++) {
        F[j] = (C1[j-1] + 2.0*C2[j-1] + 2.0*C1[j] + C2[j])/6.0;
        A[j] = -K2[j];
        D[j] = K1[j] + K2[j] + F[j];
    }
    ANE = A[nNodes-1];
    F[nNodes] = (C1[nNodes-1] + 2.0*C2[nNodes-1])/6.0;
    D[nNodes] = -A[nNodes] + F[nNodes];

    FNN1 = F[nNodes];
    F[0] = F[0]*s.h[0] - s.dz[0] * (K1[0] - (2*S[0]+S[1])/6);
    for(int j = 1; j < nNodes; j++) {
        F[j] = F[j]*s.h[j] + 0.5*(s.dz[j]+s.dz[j-1]) * K1[j]
                           - 0.5*(s.dz[j]+s.dz[j+1]) * K2[j]
                           - s.dz[j] * (S[j-1]+4*S[j]+S[j+1])/6;
    }
    F[nNodes] = F[nNodes]*H[nNodes] + s.dz[nNodes] * (K2[nNodes] - (S[nNodes-1]+2*S[nNodes])/6);

    // upper boundary condition
    if (s.ponded){
        Hnew[0] = H[0] + s.dts * (-influx);
        influx = D[0]*Hnew[0] - F[0];
        F[0] = Hnew[0];
        A1 = A[0];
        A[0] = 0;
        D[0] = 1;
        F[1] = F[1] - A1 * Hnew[0];
    } else {
       influx = throughfall;
       F[0] = F[0] + influx;
    }

    // lower boundary condition
//    Case (lowerBC):
//        0 : { // drainage
//        D[NNodes] = -A[NNodes-1];
//        FNN2 = F[NNodes];
//        F[NNodes] = 0;
//    }
//    1 : { //fixed potential in lower node
//        Drain = F[NNodes] - D[NNodes]*Hnew[NNodes];
//        F[NNodes] = Hnew[NNodes];
//        D[NNodes] = 1.0;
//        F[NNodes-1] = F[NNodes-1] - ANE*Hnew[NNodes];
//        A[NNodes-1] = 0;
//    }
//    2 : { //fixed flux: drain spec by user, e.g. when groundwater
//        Drain = GroundwaterFlux/1440;
//        F[NNodes] = F[NNodes] - Drain;
//    }

 // calc Hnew with Gaussian elimination and back substitution
    for(int j = 1; j < nNodes; j++) {
        A2 = A[j-1]/D[j-1];
        D[j] = D[j] - A2 * A[j-1];
        F[j] = F[j] - A2 * F[j-1];
    }
    Hnew[nNodes] = F[nNodes]/D[nNodes];

    for (int j = nNodes-1; j >= 0; j--) {
        Hnew[j] = (F[j] - A[j]*Hnew[j+1])/D[j];
//        if (Hnew[j] > s.hb[j])
//            Hnew[j] = s.hb[j];
        // Hnew cannot become 0?
    }

    for(int j = 1; j < nNodes; j++) {
        //if (lowerBC = 2) and (z[i] >= GWlevel) then break;

    }

  // calc boundary fluxes

  if (s.ponded && Hnew[1] < 0) {
    influx = influx + A1 * Hnew[1];
    // de influx is nu gelijk aan:
    //K1*((hnew1-hnew2)/dz+1) + 0.5*dé/dh * (hnew1-h1) }
  }

//    Case LowerBC of
//    // free drainage
//    0 : {
//            Drain =0.5*(fKWH(mat[nNodes], 2, H[nNodes], false) +
//                             fKWH(mat[nNodes], 2, Hnew[nNodes], false));
//    if (Hnew[nNodes] < 0)
//       Drain = 0.5*(drain + FNN2 - Hnew[nNodes]*FNN1)
//    else
//       Drain = 0.5*(drain + FNN2);
//    }
//    // fixed potential
//    1 : Drain = Drain - ANE*Hnew[nNodes-1];
//}


}

void TWorld::calcNewNodalValues(long i_, double *Hnew, double *Wnew, double Cnew)
{
    SOIL_LIST s;
    s = crSoil[i_];


    for(int j = 0; j < nNodes; j++) {
        double Hx = 0.5*(Hnew[j] + s.h[j]);

        //check is H negative or not!, note hb is negative
        if (Hx < s.hb[j])
            s.K[j] = s.Ks[j]*pow(Hx/s.hb[j], 2.0+3.0*s.lambda[j]);
        else
            s.K[j] = s.Ks[j];

        if (Hnew[j] < s.hb[j]) {
            Wnew[j] = s.thetar[j] + pow(s.hb[j]/Hnew[j], s.lambda[j]); * (s.pore[j]-s.thetar[j]);
            Cnew[j] = (s.pore[j]-s.thetar[j])*s.lambda[j]*pow(s.hb[j]/Hnew[j], s.lambda[j]-1);
        } else {
            Wnew[j] = s.pore[j];
            Cnew[j] = Wnew[j]*1e-6/s.pore[j];
        }
        Cnew[j] *= s.dz[j]/s.dts; //????????????? check
    }
}

void TWorld::cell_Soilwater(long i_)
{
    SOIL_LIST s;
    s = crSoil[i_];

    double dtmin = 0.1;
    bool stopit = false;
    double Hold[nNodes];
    double Hnew[nNodes];
    double Wnew[nNodes];
    double Cnew[nNodes];

    s.dtsum = 0;

    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.h[j];
        Wnew[j] = s.theta[j]; //????
        Cnew[j] =
    }

    // outer loop timestep lisem
    do {
        //calcSinkterm(i_);

        int NIT = 0;
        int NITMAX = 12;

        // iteration per gridcell
        do {
            s.ponded = Hnew[0] > 0;

            for(int j = 0; j < nNodes; j++)
                Hold[j] = Hnew[j];

            calcNewNodalValues(i_, Hnew, Wnew);

            solveFiniteElement(i_, influx, Hnew, Wnew, Cnew);

            s.ponded = Hnew[0] > 0;

            stopit = true;
            double tol = 0;
            for(int j = 0; j < nNodes; j++) {
                tol = 0.01*fabs(Hold[j]) + 0.2;
                if (s.ponded)tol /= 2.0;
                if (fabs(Hnew[j] - Hold[j]) > tol) {
                    stopit = false;
                    break;
                }
            }
            NIT++;

            if (!stopit && NIT > NITMAX) {
                s.dts /= 2.0;
                s.dtsum -= s.dts;
                for(int j = 0; j < nNodes; j++)
                    Hnew[j] = 0.5*(s.h[j]+Hnew[j]);
            }

        } while(stopit || s.dts < dtmin);

        s.dtsum += s.dts;
        double dtlast = s.dts;
        double factor =0.5 + 1/sqrt((double)NIT);
        s.dts = s.dts * factor;
        s.dts = std::min(s.dts,_dt);
        s.dts = std::max(s.dts,dtmin);
        factor = s.dts/dtlast;

        for(int j = 0; j < nNodes; j++) {
            double dH =(Hnew[j] - s.h[j])*factor;
            s.h[j] = Hnew[j];
            Hnew[j] = Hnew[j] + dH;
        }
        //move(Wnew^, W^, (NN+1) * sizeof(double));

    } while(s.dtsum < _dt);

//??????????
//    if (s.dts < dtmin)
//        s.dts = dtorg;

    // calcCorrections(i_);

}



