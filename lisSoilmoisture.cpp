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


//DO NOT USE OMP, MUCH SLOWER
// Calculates one dimensional soil water balance based on
// Richards equation, using a fully implicit, mass lumped,
// Galerkin scheme Finite Element Method.
// Based on FORTRAN program WORM by Van Genuchten (1987),
// Research report no.121 of USDA Salinity Laboratory and
// modified by J.B.Kool and G.H. de Rooij (1989) LUW.
//
// Adapted for rainfall, evapotranspiration and groundwater level
// change. Conductivity and retention curves according to
// Brooks Corey functions (1966).
// H, Hnew : pressure head begin and end of timestep dt
// K1, K2 : hydraulic conductivity upper and lower layer
// C1     : differential water capacity, 2 dif. methods
// F, A, D : Galerkin scheme tri-diagonal matrix coefficients
// S : sink term = ET


void TWorld::calcSinkterm(long i_, double *S)
{
    SOIL_LIST s = crSoil[i_];
    int r = s.r;
    int c = s.c;
    s.ponded = s.h[0] >= 0;

    double ETafactor = 1;
    double Ld = 12;

    // SwitchDailyET = true;
    double day = trunc(time/86400.0);
    double hour = std::min(24.0,std::max(0.0, time/3600.0-day*24.0));
    if (SwitchDailyET) {
        double declination = -23.45 * M_PI/180.0 * cos(2*M_PI*(day+10)/365.0);
        Ld = 24.0/M_PI*(acos(-tan(declination)*tan(latitude/180.0*M_PI)));  // daylength in hour
        if (std::isnan(Ld)) Ld = 12.0;
        ETafactor = std::max(0.,sin((-0.5-hour/Ld)*M_PI)) / Ld*_dt/3600.0*M_PI*0.5;
        //<= this ensures that the sum of all steps in a day amounts to the daily ET, regardless of _dt
    }

//    if (Rain->Drc* 3600000.0/_dt > rainfallETa_threshold) {
//        ETa->Drc = 0;
//        ETp->Drc = 0;
//    }

    if (i_ == 3738) qDebug() << day << hour << ETp->Drc;
    if (ETp->Drc*ETafactor > 0 && Rain->Drc* 3600000.0/_dt < rainfallETa_threshold) {
        double AreaSoil = SoilWidthDX->Drc * DX->Drc;
        double Cover_ = Cover->Drc;
        double ETp_ = ETp->Drc * ETafactor;
        double tot = 0;
        double etanet = 0;

        ETpCum->Drc += ETp_;
/*
        //  interception decrease, drying out canopy
        double CStor_  = CStor->Drc;
        if (CStor_ > 0) {
            double ETa_int = ETp_;

            ETa_int = std::min(ETa_int, CStor_);
            CStor_ = CStor_- ETa_int;

            RainCum->Drc = std::max(0.0, RainCum->Drc-ETa_int);
            if (CStor_ < 1e-6)
                RainCum->Drc = 0;

            etanet = std::max(0.0, ETp_ - ETa_int);
            Interc->Drc = Cover_ * CStor_ * CHAdjDX->Drc; //????
            IntercETa->Drc += Cover_ * ETa_int * CHAdjDX->Drc;
            CStor->Drc = CStor_;
        }

        if (SwitchLitter) {
            double CvL = Litter->Drc;
            double LCS = LCStor->Drc;

            double ETa_int = std::min(etanet, LCS);
            etanet = std::max(0.0, ETp_ - ETa_int);
            LCStor->Drc = LCS- ETa_int;
            IntercETa->Drc += CvL * ETa_int * AreaSoil;
            LInterc->Drc =  CvL * LCS * AreaSoil;
        }

        if (SwitchHouses)
        {
            double CvH = HouseCover->Drc;
            double HS = HStor->Drc;

            double ETa_int = std::min(etanet, HS);
            etanet = std::max(0.0, ETp_ - ETa_int);
            HStor->Drc = HS - ETa_int;
            IntercETa->Drc += CvH * ETa_int * AreaSoil;
            double roofsurface = (_dx * DX->Drc * CvH); // m2
            IntercHouse->Drc =  roofsurface * HS;
        }
*/
        // if no water on the surface
        if (!s.ponded) {
            //transpiration under Cover
            for (int j = 0; j < nNodes; j++) {
                if (s.rootz[j] = 0)
                    break;

                // van genuchten H50 = -3.5 m
                double f = 1.0/(1.0+pow(s.h[j]/-3.5,1.5));
                if (s.h[j] > s.hb[j]) f = 0; // saturation
                if (s.h[j] < -16.0) f = 0; // wilting point -16000 cm
                S[j] =  ETp_ * Cover_ * f * s.rootz[j]; // etanet * Cover_ *
            }
            // surface evaporation (1-Cover)
            S[0] += (s.theta[0]-s.thetar[0])/(s.pore[0]-s.thetar[0])*ETp_*(1-Cover_);

                if (i_==3738) {
                    QString ss;
                    for(int j = 0; j < nNodes; j++) {
                        ss = ss + QString(" %1").arg(S[j]);
                    }
                    qDebug() << Cover_ << ETp_ << ss;
                }

            for (int j = 0; j < nNodes; j++)
                tot += S[0];
        } else {
            // ponded
            for (int j = 0; j < nNodes; j++)
                S[j] = 0;

            double ETa_pond = ETp_;
            if (hmxWH->Drc > ETp_) {
                double WHRunoff_ = WHrunoff->Drc;
                ETa_pond = std::min(ETa_pond, WHRunoff_);
                WHRunoff_ = WHRunoff_ - ETa_pond;
                WHroad->Drc = WHRunoff_;
                WH->Drc = WHRunoff_ + WHstore->Drc;
                WHrunoff->Drc = WHRunoff_;
                WaterVolall->Drc = CHAdjDX->Drc * (WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;

                tot = tot + ETa_pond;
            }
        }
        ETa->Drc = tot;
        ETaCum->Drc += tot;
    }
}

void TWorld::cell_Soilwater(long i_)
{
    SOIL_LIST s;
    s = crSoil[i_];

    double dtmin = 0.01;
    int NITMAX = 12;
    bool stopit = false;
    int c = s.c;
    int r = s.r;    
    double WH0;
    double cnt = 0;
    int nN = nNodes-1;
    bool freeDrainage = !SwitchImpermeable;

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

    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc; // flood in kin wave

    s.InfPot = WH0/_dt;// m/s
    s.Infact = 0;
    s.dtsum = 0;
    s.drain = 0;

    //#pragma omp parallel for num_threads(userCores)
    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.h[j];
        S[j] = 0;
    }

    calcSinkterm(i_, S); // ETa

    // outer loop timestep lisem
    do {
        int NIT = 0;

        // iteration per gridcell
        do {
            s.ponded = Hnew[0] >= 0;

            //#pragma omp parallel for num_threads(userCores)
            for(int j = 0; j < nNodes; j++)
                Hold[j] = Hnew[j];

            // Brooks Corey for H,K,theta and C
            //calcNewNodalValues(i_, Hnew, K, C1);
            // Galerkin 3-diagonal scheme and back substitution
            //solveFiniteElement(i_, Hnew, K, C1);

            // values of K, theta and C for H
            //#pragma omp parallel for num_threads(userCores)
            for(int j = 0; j < nNodes; j++) {
                double Hx = std::min(s.hb[j], 0.5*(Hnew[j] + s.h[j]));

                K[j] = s.Ks[j]*pow(s.hb[j]/Hx, 2.0+3.0*s.lambda[j]);

//                if (Hx < s.hb[j])
//                    K[j] = s.Ks[j]*pow(s.hb[j]/Hx, 2.0+3.0*s.lambda[j]);
//                else
//                    K[j] = s.Ks[j];

                double Wnew;
                Wnew = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/std::min(Hnew[j],s.hb[j]), s.lambda[j]);
//                if (Hnew[j] < s.hb[j])
//                    Wnew = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hnew[j], s.lambda[j]);
//                else
//                    Wnew = s.pore[j];

                if (fabs(Hnew[j]-s.h[j]) <= 3*tol2) {
                    if (Hx < s.hb[j])
                        C1[j] = s.theta[j] * 1e-6/s.pore[j] - 1.0/Hx *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(s.hb[j]/Hx, s.lambda[j]);
                    else
                        C1[j] = s.theta[j] * 1e-6/s.pore[j];
                } else
                    C1[j] = (Wnew-s.theta[j])/(Hnew[j]-s.h[j]);
                C1[j] *= s.dz[j]/s.dts;
            }

            // K1 and K2 are avg between node and upper and lower node
            //#pragma omp parallel for num_threads(userCores)
            for(int j = 1; j < nNodes; j++)
                K1[j] = Aavg(K[j-1],K[j]);

            //#pragma omp parallel for num_threads(userCores)
            for(int j = 0; j < nNodes-1; j++)
                K2[j] = Aavg(K[j],K[j+1]);
            K1[0] = K2[0];
            K2[nN] = K1[nN];

            // sort of Simpson trapezium solution
            F[0] = (2.0*C1[0] + C1[1])/6.0;
            A[0] = -K1[0];
            D[0] = -A[0]+F[0];
            //#pragma omp parallel for num_threads(userCores)
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
            //#pragma omp parallel for num_threads(userCores)
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
            //#pragma omp parallel for num_threads(userCores)
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

          //  //#pragma omp parallel for num_threads(userCores)
//            for(int j = 1; j < nNodes; j++)
//                Hnew[j] = std::min(0.0, Hnew[j]);
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
                //#pragma omp parallel for num_threads(userCores)
                for(int j = 0; j < nNodes; j++)
                    Hnew[j] = 0.5*(s.h[j]+Hnew[j]);
                stopit = true;
            }

        } while(!stopit);
        // end SoilMoisture in PAS code

        //#pragma omp parallel for num_threads(userCores)
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
    //#pragma omp parallel for num_threads(userCores)
    for (int j = 1; j < nN1_; j++)
        ThetaI1a->Drc += s.theta[j];
    ThetaI1a->Drc /= (double)nN1_;

    ThetaI2a->Drc = s.theta[nN1_];
    //#pragma omp parallel for num_threads(userCores)
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
//    if (i_==3738) {
//        QString S;
//        for(int j = 0; j < nNodes; j++) {
//            S = S + QString(" %1").arg(s.h[j]);
//        }
//        //qDebug() << s.Infact << S;
//    }

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



