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

#define tol2 0.2
#define tol1 0.01
#define cell (_nrRows/2*_nrCols)+_nrCols/3
#define show(a,b,c) if (i_==cell) qDebug() << a << b << c

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

double TWorld::calculateDayLength(double latitude, int dayNumber)
{
    const double degreesToRadians = PI / 180.0;

    // Convert latitude from degrees to radians
    latitude *= degreesToRadians;

    // Earth's axial tilt in degrees

    const double axialTiltRadians = 23.44 * degreesToRadians;

    // Day angle in radians
    double dayAngle = 2 * PI * (dayNumber - 1) / 365;

    // Calculate the declination angle in radians
    double declination = asin(sin(axialTiltRadians) * sin(dayAngle));

    // Calculate the hour angle at sunrise and sunset in radians
    double hourAngle = acos(-tan(latitude) * tan(declination));

    // Calculate day length in hours
    double dayLength = (2.0 * hourAngle) * (180.0 / PI) / 15.0;

    return dayLength;
}


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
        //Ld = calculateDayLength(day, latitude);
        Ld = (2.0*acos(-tan(latitude*0.01745329) * tan(asin(0.397789 * sin(0.017214*(day-1)))))) * 3.8197186;
        ETafactor = std::max(0.0,sin((-0.5-hour/Ld)*PI)) / Ld*_dt/3600.0*PI*0.5;
        //<= this ensures that the sum of all steps in a day amounts to the daily ET, regardless of _dt
    }

//    if (i_ == 3738) qDebug() << day << hour << Ld;

    if (ETp->Drc*ETafactor > 0 && Rain->Drc* 3600000.0/_dt < rainfallETa_threshold) {
        double Cover_ = Cover->Drc;
        double ETp_ = ETp->Drc * ETafactor;
        double tot = 0;
        double etanet = ETp_;

        ETpCum->Drc += ETp_;

        //  interception decrease, drying out canopy, litter, roofs
        double CStor_  = CStor->Drc;
        if (CStor_ > 0) {
            double ETa_int = ETp_;

            ETa_int = std::min(ETa_int, CStor_);
            CStor_ = CStor_- ETa_int;

            RainCum->Drc = std::max(0.0, RainCum->Drc-ETa_int);
            if (CStor_ < 1e-6)
                RainCum->Drc = 0;

            etanet = std::max(0.0, etanet - Cover_*ETa_int);
            Interc->Drc = Cover_ * CStor_ * CHAdjDX->Drc;
            IntercETa->Drc += Cover_ * ETa_int * CHAdjDX->Drc;
            CStor->Drc = CStor_;
        }

        if (SwitchLitter) {
            double CvL = Litter->Drc;
            double LCS = LCStor->Drc;
            if (CvL > 0 && LCS > 0) {
                double ETa_int = std::min(etanet, LCS);
                etanet = std::max(0.0, etanet - CvL*ETa_int);
                LCStor->Drc = LCS - ETa_int;
                IntercETa->Drc += CvL * ETa_int * CHAdjDX->Drc;
                LInterc->Drc =  CvL * LCS * CHAdjDX->Drc;
            }
        }

        if (SwitchHouses)
        {
            double CvH = HouseCover->Drc;
            double HS = HStor->Drc;
            if (CvH > 0 && HS > 0) {
                double ETa_int = std::min(etanet, HS);
                etanet = std::max(0.0, etanet - CvH * ETa_int);
                HStor->Drc = HS - ETa_int;
                IntercETa->Drc += CvH * ETa_int * CHAdjDX->Drc;
                //double roofsurface = (_dx * DX->Drc * CvH); // m2
                IntercHouse->Drc =  (_dx * DX->Drc * CvH) * HS;
            }
        }
        //transpiration under Cover from rootzone
        for (int j = 0; j < nNodes; j++) {
            // van genuchten H50 = -3.5 m
            double f = 1.0/(1.0+pow(s.h[j]/-3.5,1.5));
            if (s.h[j] > s.hb[j]) f = 0; // saturation
            if (s.h[j] < -16.0) f = 0; // wilting point -16000 cm
            S[j] =  etanet * Cover_ * f * s.rootz[j];
        }

        // add surface evaporation (1-Cover) to top node
        if (hmxWH->Drc < 1e-6)
           S[0] += (s.theta[0]-s.thetar[0])/(s.pore[0]-s.thetar[0])*etanet*(1-Cover_);

        for (int j = 0; j < nNodes; j++) {
            tot += S[j];
        }


        // ponded, evap only outside cover
        double ETa_pond = etanet*(1-Cover_);
        if (hmxWH->Drc > ETa_pond) {
            ETa_pond = std::min(ETa_pond, WH->Drc);
            WH->Drc -= ETa_pond;
            WHrunoff->Drc = std::max(0.0, WH->Drc - WHstore->Drc);
            WHroad->Drc = WHrunoff->Drc;
            WHrunoff->Drc = WHrunoff->Drc;
            WaterVolall->Drc = CHAdjDX->Drc * (WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;
            tot = tot + ETa_pond;
        }

        ETa->Drc = tot;
        ETaCum->Drc += tot;

        if (i_== cell) {
            QString ss;
            for(int j = 0; j < nNodes; j++) {
                ss = ss + QString(" %1").arg(S[j]);
            }
         //   qDebug() << tot << etanet << ss;
        }

    }
}

void TWorld::cell_Soilwater(long i_)
{
    SOIL_LIST s;
    s = crSoil[i_];

    double dtmin = 0.01*_dt;
    double dtmax = SoilWBdtfactor*_dt;
    int NITMAX = 12;
    bool stopit = false;
    int c = s.c;
    int r = s.r;    
    double WH0;
    double cnt = 0;
    int nN = nNodes-1;
    bool freeDrainage = !SwitchImpermeable;

    double ANE, FNN1, FNN2, A1, A2;

    double *Hold = new double[nNodes];
    double *Hnew = new double[nNodes];
    double *Wnew = new double[nNodes];
    double *C1 = new double[nNodes];
    double *S  = new double[nNodes];
    double *K  = new double[nNodes];
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
    s.Infact = WH0/_dt;
    s.dtsum = 0;
    s.drain = 0;

    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.h[j];
        Hold[j] = s.h[j];
        S[j] = 0;
    }

    if (SwitchIncludeET)
        calcSinkterm(i_, S); // ETa

    s.Ks[0] = Ksateff->Drc/_dt;
    s.pore[0] = Poreeff->Drc;
    // include crusting compaction etc in top layer;


    s.dts=0.5*SoilWBdtfactor*_dt;

    // outer loop timestep lisem
    do {
        int NIT = 0;
        cnt = 0;

        // iteration
        do {

            //======= Brooks Corey for H,K,theta and C
            for(int j = 0; j < nNodes; j++) {
                double Hx = std::min(s.hb[j], 0.5*(Hnew[j] + Hold[j]));

                K[j] = s.Ks[j]*pow(std::min(1.0,s.hb[j]/Hx), 2.0+3.0*s.lambda[j]);

                if (Hnew[j] <= s.hb[j])
                    Wnew[j] = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hnew[j], s.lambda[j]);
                else
                    Wnew[j] = s.pore[j];
                double W = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/(Hnew[j]-0.01), s.lambda[j]);
                C1[j] = (Wnew[j]-W)/0.01;

//                if (fabs(Hnew[j]-Hold[j]) <= 3*tol2) {
//                    if (Hx < s.hb[j])
//                        C1[j] = s.theta[j] * 1e-6/s.pore[j]
//                                -1.0/Hx *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(s.hb[j]/Hx, s.lambda[j]);
//                    else {
//                        //double f = std::max(1e-6,(0-Hx)/(0-s.hb[j]));
//                        C1[j] = s.theta[j] * 1e-6/s.pore[j]; //????
//                    }
//                } else
//                    C1[j] = (Wnew[j]-s.theta[j])/(Hnew[j]-Hold[j]);

                C1[j] *= s.dz[j]/s.dts;
             }

            // K1 and K2 are avg between node and upper and lower node
            for(int j = 1; j < nNodes; j++) {
                switch (KavgType) {
                    case 0: K1[j] = Aavg(K[j-1],K[j]);
                    case 1: K1[j] = Savg(K[j-1],K[j]);
                    case 2: K1[j] = Havg(K[j-1],K[j],s.dz[j-1],s.dz[j]);
                    case 3: K1[j] = Mavg(K[j-1],K[j]);
                }
            }
            for(int j = 0; j < nNodes-1; j++) {
                switch (KavgType) {
                    case 0: K2[j] = Aavg(K[j],K[j+1]);
                    case 1: K2[j] = Savg(K[j],K[j+1]);
                    case 2: K2[j] = Havg(K[j],K[j+1],s.dz[j],s.dz[j+1]);
                    case 3: K2[j] = Mavg(K[j],K[j+1]);
                }
            }
            K1[0] = K2[0];
            K2[nN] = K1[nN];

            s.ponded = Hnew[0] > 0;
//            if (!s.ponded) {
//               // calculate available space in profile in cm: (pore-theta)*dz
//               double space = 0;
//               for(int j = 0; j < nNodes; j++)
//                   space += (s.pore[j] - Wnew[j]) * s.dz[j];

//               s.ponded = WH0 > space;
//            }

            s.Infact = -s.Ks[0]*(WH0-std::min(s.hb[0],Hnew[0]))/s.dz[0] + K[0];
            if (s.Infact <= s.InfPot)
               s.ponded = true;

            //======== Galerkin 3-diagonal scheme and back substitution

            // sort of Simpson trapezium solution
            F[0] = (2.0*C1[0] + C1[1])/6.0;
            A[0] = -K1[0];
            D[0] = -A[0]+F[0];
            for(int j = 1; j < nNodes-1; j++) {
                F[j] = (C1[j-1] + 4.0*C1[j] + C1[j+1])/6.0; // dtheta/dh * dz/dt
                A[j] = -K2[j];
                D[j] = (K1[j] + K2[j]) + F[j];
            }
            ANE = A[nN-1];
            //A[nN] = -K2[nN];
            F[nN] = (C1[nN-1] + 2.0*C1[nN])/6.0;
            D[nN] = -A[nN-1] + F[nN];
            FNN1 = F[nN];
            // include gravity term (i.e. k only) and sink term S[]
            F[0] = F[0]*Hold[0] - s.dz[0]*(K1[0] - (2*S[0]+S[1])/6);
            for(int j = 1; j < nNodes-1; j++) {
                F[j] = F[j]*Hold[j]
                       + 0.5*(s.dz[j]+s.dz[j-1]) * K1[j]
                       - 0.5*(s.dz[j]+s.dz[j+1]) * K2[j]
                       - s.dz[j] * (S[j-1]+4*S[j]+S[j+1])/6;
            }
            F[nN] = F[nN]*Hold[nN] + 0.5*s.dz[nN]*(K2[nN] - (S[nN-1]+2*S[nN])/3);

            // upper boundary condition
            if (s.ponded){
                Hnew[0] = Hold[0] +  s.Infact*s.dts;
                // in worm staat hiet de potential aan het oppervlak dus WH0, alleen als [0] de oppervlakte is
                s.Infact = D[0]*Hnew[0] - F[0];

              //  F[0] = Hnew[0];  //??? F is a flux, hoe kun je H erbij optellen?
                F[0] = F[0] + s.Infact;
                A1 = A[0];
                A[0] = 0;
                D[0] = 1;
                F[1] = F[1] - A1 * Hnew[0];
            } else {
                s.Infact = s.InfPot;
                F[0] = F[0] + s.Infact; // direct adding, units?? F is a flux
            }

            // lower boundary condition
            if (freeDrainage) {
                D[nN] = -A[nN-1];
                FNN2 = F[nN];
                F[nN] = 0;
            } else {
                // fixed flux (for instance 0)
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
            for (int j = nNodes-2; j >= 0; j--) {
                Hnew[j] = (F[j] - A[j]*Hnew[j+1])/D[j];
            }

            for(int j = 0; j < nNodes; j++)
                Hnew[j] = std::min(0.0, Hnew[j]);
            // CHECK THIS? GW can be positive

            //======== calc boundary fluxes

            //??????????????????
//            if (s.ponded && Hnew[1] < 0) {
//                s.Infact = s.Infact + A1 * Hnew[1];
//            }
            if (freeDrainage)
                s.drain = 0.5*(K[nN] + FNN2 - std::min(Hnew[nN],0.0)*FNN1);

            s.ponded = Hnew[0] > 0;

            // stop if Hnew and Hold are close
            stopit = true;
            for(int j = 0; j < nNodes; j++) {
                double tol = tol1*fabs(Hold[j]) + tol2;
                if (s.ponded) tol /= 2;
                if (fabs(Hnew[j] - Hold[j]) > tol) {
                    stopit = false;
                    break;
                }
            }
            NIT++;

            if (!stopit && NIT > NITMAX) {
                // restore and do again
                s.dts /= 2.0;
                s.dts = std::max(s.dts,dtmin);
                for(int j = 1; j < nNodes; j++)  // do not include surface node
                  //  Hnew[j] = Hold[j];
                    Hnew[j] = 0.5*(Hold[j]+Hnew[j]);
                stopit = true;
            }

            for(int j = 0; j < nNodes; j++) {
                Hold[j] = Hnew[j];
            }


           //if (NIT > NITMAX) stopit = true;
        } while(!stopit);
        // end SoilMoisture in PAS code


        // change timestep for next iteration
//        double factor = 0.5 + 1/sqrt((double)NIT);
//        s.dts = s.dts * factor;
//        s.dts = std::max(s.dts,dtmin);
//        s.dts = std::min(s.dts,_dt-s.dtsum);
//        s.dtsum += s.dts;

        // SWATRE method is faster!
        double dt = dtmax;
        for(int j = 0; j < nNodes; j++) {
           double mdih = 0.2 + 0.02*fabs(Hnew[j]);
           double dih  = fabs(Hnew[j] - Hold[j]);
           if (dih > 0.10)
              dt = dt*mdih/dih;
        }
        s.dts = std::min(dt,dtmax);
        s.dts = std::max(dt,dtmin);
        s.dts = std::min(s.dts,_dt-s.dtsum);
        s.dtsum += s.dts;

        if (i_ == cell) qDebug() << NIT << s.dts << s.dtsum << _dt;
    } while(s.dtsum < _dt);

    for(int j = 0; j < nNodes; j++) {
        s.h[j] = Hnew[j];
        if (s.h[j] < s.hb[j])
            s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
        else
            s.theta[j] = s.pore[j];
    }

    double in = s.Infact*_dt;
    if (FloodDomain->Drc == 0) {
        in = std::min(WH->Drc , in);
        WH->Drc -= in; //runoff in kinwave or dyn wave
    } else {
        in = std::min(hmx->Drc , in);
        hmx->Drc -= in; // flood in kin wave
    }
    s.Infact = in/_dt;

    InfilVol->Drc = s.Infact*_dt*SoilWidthDX->Drc*DX->Drc;

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
    if (i_ == cell) qDebug() << cnt << WH0 << s.InfPot << s.Infact  << s.ponded;
    if (i_==cell) {
        QString S;
        QString S1;
        for(int j = 0; j < nNodes; j++) {
            S = S + QString(" %1").arg(s.h[j]);
            S1 = S1 + QString(" %1").arg(s.theta[j]);
        }
        qDebug() << S;
        qDebug() << S1;
    }

    delete[] Hold;
    delete[] Hnew;
    delete[] Wnew;
    delete[] C1;
    delete[] S;
    delete[] K;
    delete[] K1;
    delete[] K2;
    delete[] A ;
    delete[] D ;
    delete[] F ;
}


void TWorld::cell_SWATRECalc(long i_)
{
    SOIL_LIST s;
    s = crSoil[i_];

    double dtmin = 0.01;
    int c = s.c;
    int r = s.r;

    double WH0;
    int nN = nNodes-1;
    double qbot = 0;
    double qtop = s.InfPot;
    double alpha;
    bool fltsat = false;

    double *Hold = new double [nNodes];
    double *Hnew = new double [nNodes];
    double *K = new double [nNodes];
    double *C1 = new double [nNodes];
    double *S  = new double [nNodes];
    double *kavg  = new double [nNodes];
    double *thoma  = new double[nNodes];
    double *thomb  = new double[nNodes];
    double *thomc  = new double[nNodes];
    double *thomf  = new double[nNodes];
    double *beta  = new double[nNodes];
    double *disnod = new double[nNodes];
    double *Wnew = new double[nNodes];
    double *z = new double[nNodes];


    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc; // flood in kin wave

    s.InfPot = WH0/_dt;// m/s   //qtop
    s.Infact = WH0/_dt;
    s.dtsum = 0;

    s.Ks[0] = Ksateff->Drc/_dt;
    s.pore[0] = Poreeff->Drc;
    // include crusting compaction etc in top layer;

    //#pragma omp parallel for num_threads(userCores)
    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.h[j];
        Hold[j] = s.h[j];
        S[j] = 0;
    }

//    if (SwitchIncludeET)
//        calcSinkterm(i_, S); // ETa

    // outer loop timestep lisem
    do {

            //======= Brooks Corey for H,K,theta and C
            for(int j = 0; j < nNodes; j++) {
                disnod[j] = s.dz[j];

            double Hx = Hnew[j] > s.hb[j] ? s.hb[j] : Hnew[j];
                    //std::min(Hnew[j],s.hb[j]);
                K[j] = s.Ks[j]*pow(s.hb[j]/Hx, 2.0+3.0*s.lambda[j]);

                Wnew[j] = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hx, s.lambda[j]);
                double W = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/(Hx-0.01), s.lambda[j]);
                C1[j] = (Wnew[j]-W)/0.01;
              //  C1[j] = -1.0/Hx *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(Hx/Hnew[j], s.lambda[j]);
           }

            for (int j = 1; j < nNodes; j++) {
                switch (KavgType) {
                    case 0: kavg[j] = Aavg(K[j],K[j-1]);
                    case 1: kavg[j] = Savg(K[j],K[j-1]);
                    case 2: kavg[j] = Havg(K[j],K[j-1],1.0,1.0);
                    case 3: kavg[j] = Mavg(K[j],K[j-1]);
                }
            }
            kavg[0] = kavg[1];

            //--- boundary conditions ---//

        //----- TOP -----//
        qtop = -WH0/s.dts;
        //----- BOTTOM -----//
        // bottom is 0 or copy of flux of last 2 layers
        if (SwitchImpermeable)
           qbot = 0;
        else
           qbot = kavg[nN]*(Hnew[nN]-Hnew[nN-1])/disnod[nN] - kavg[nN];

        // 1st check flux aginst max flux
        double qmax = Savg( s.Ks[0],K[0])*(WH0-Hnew[0]) / disnod[0] - kavg[0];

        if (i_ == cell)
           qDebug() << "qm" << qtop << qmax << WH0 << K[0] << Hnew[0];

        // maximum possible flux, compare to real top flux available
        s.ponded = (qtop < qmax); // both negative, why?

        if (!s.ponded) {
           // calculate available space in profile in cm: (pore-theta)*dz
           double space = 0;
           for(int j = 0; j < nNodes; j++)
               space += (s.pore[j] - Wnew[j]) * s.dz[j];

           s.ponded = abs(-qtop * s.dts) > space;
        }
        fltsat = true;
        for(int j = 0; j < nNodes; j++)
            if (Hnew[j] < 0) fltsat = false;

        for(int j = 0; j < nNodes; j++) {
            Hold[j] = Hnew[j];
            s.theta[j] = Wnew[j]; // thetaprev? used where
        }

        //====== headcalc
        //First node : 0 (include boundary cond. qtop or pond)
        if ( s.ponded ) //|| (fltsat && (qtop <= qbot)) )
        {
            thoma[0] = -s.dts * kavg[0];
            thomc[0] = -s.dts * kavg[1];
            thomb[0] = -thomc[0] + C1[0] + s.dts*kavg[0];
            thomf[0] = C1[0]*s.h[0] +s.dts/(-s.dz[0]) * (kavg[0] - kavg[1]) + s.dts*kavg[0]*WH0;
        }
        else
        {
            s.ponded = false;
            thoma[0] = -s.dts * kavg[0];
            thomc[0] = -s.dts * kavg[1];
            thomb[0] = -thomc[0] + C1[0];
            thomf[0] = C1[0]*s.h[0] +
                       s.dts/(-s.dz[0]) * (- qtop - kavg[1]);
        }

        //Intermediate nodes: i = 1 to n-2
        for (int j = 1; j < nNodes-1; j++)
        {
            thoma[j] = -s.dts*kavg[j]; // /s.dz[j]/disnod[j];
            thomc[j] = -s.dts*kavg[j+1]; // /s.dz[j]/disnod[j+1];
            thomb[j] = -thoma[j] - thomc[j] + C1[j];
            thomf[j] = C1[j]*Hnew[j] + s.dts/(-s.dz[j])*(kavg[j]-kavg[j+1]);
        }
        // last node
        thoma[nN] = -s.dts*kavg[nN]; // /s.dz[nN]/disnod[nN];
        thomb[nN]= -thoma[nN] + C1[nN];
        thomf[nN] = C1[nN]*Hnew[nN] + s.dts/(-s.dz[nN])*(kavg[nN]+qbot);


        //Gaussian elimination and backsubstitution h - first time
        alpha = thomb[0];
        Hnew[0] = thomf[0] / alpha;
        for (int j = 1; j < nN; j++) {
            beta[j] = thomc[j-1] / alpha;
            alpha = thomb[j] - thoma[j] * beta[j];
            Hnew[j] = (thomf[j] - thoma[j] * Hnew[j-1]) / alpha;
        }
        for (int j = (nN-1); j >= 0; j--)
            Hnew[j] -= beta[j+1] * Hnew[j+1];

        //correct tridiagonal matrix
        for (int j = 0; j < nNodes; j++) {
            double Hx = std::min(Hnew[j],s.hb[j]);

        Wnew[j] = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hx, s.lambda[j]);
        double W = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/(Hx-0.01), s.lambda[j]);
        C1[j] = (Wnew[j]-W)/0.01;
            double dimocaNew = -1.0/Hx *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(Hx/Hnew[j], s.lambda[j]);
            thomb[j] = thomb[j] - C1[j] + dimocaNew;
            thomf[j] = thomf[j] - C1[j]*Hold[j] + dimocaNew*Hnew[j]
                       - Wnew[j] + s.theta[j];
        }

        // back substitution
        alpha = thomb[0];
        Hnew[0] = thomf[0] / alpha;
        for (int j = 1; j < nNodes; j++) {
            beta[j] = thomc[j-1] / alpha;
            alpha = thomb[j] - thoma[j] * beta[j];
            Hnew[j] = (thomf[j] - thoma[j] * Hnew[j-1]) / alpha;
        }
        for (int j = (nN-1); j >= 0; j--)
            Hnew[j] -= beta[j+1] * Hnew[j+1];

        for (int j = 0; j < nNodes; j++)
            if(Hnew[j] > 0) Hnew[j] = 0;

        //====== headcalc

        if (i_==cell) {
            QString S;
            QString S1;
            QString S2;
            QString S3;
            for(int j = 0; j < nNodes; j++) {
               S = S + QString(" %1").arg(thoma[j]);
               S1 = S1 + QString(" %1").arg(thomb[j]);
               S2 = S2+ QString(" %1").arg(Hnew[j]);
               S3 = S3 + QString(" %1").arg(C1[j]);
            }
        //   qDebug() << S;
        //   qDebug() << S1;
         //   qDebug() << S2;
         //   qDebug() << S3;
        }
        // determine new boundary fluxes
        if (SwitchImpermeable)
           qbot = 0;
        else
           qbot = kavg[nN]*(Hnew[nN]-Hnew[nN-1])/disnod[nN] - kavg[nN];

        if ( s.ponded || (fltsat && (qtop < qbot)) )
           qtop = -kavg[0] * ((Hnew[0] - WH0)/disnod[0] + 1);
        // adjust top flux

        WH0 += qtop*s.dts;
        // decrease pond with top flux
        s.Infact -= qtop;
        // this is the last flux?

        double dt = s.dts;
        double precParam = 5.0;
        double accur1 = 0.3 - 0.02 * precParam;
        double accur2 = 0.03 - 0.002 * precParam;

        for(int j = 0; j < nNodes; j++)
        {
           double mdih = accur1 + accur2 * std::max(1.0, fabs(Hnew[j]));
           double dih  = fabs(Hnew[j] - Hold[j]);
           if (dih > 0.10)
              dt = std::min(dt, s.dts*mdih/dih);
        }
        s.dts = std::max(dt,dtmin);
        s.dts = std::min(s.dts,_dt-s.dtsum);
        s.dtsum += s.dts;


       // if (i_ == 3738) qDebug() << cnt << NIT << s.dts << s.dtsum << _dt;
    } while(s.dtsum < _dt);

    if (FloodDomain->Drc == 0) {
//        s.Infact = std::min(s.Infact, WH->Drc);
//        WH->Drc -= s.Infact; //runoff in kinwave or dyn wave
        WH->Drc = WH0; //runoff in kinwave or dyn wave
    } else {
        //s.Infact = std::min(s.Infact, hmx->Drc);
        hmx->Drc = WH0;//s.Infact; // flood in kin wave
    }

    InfilVol->Drc = -s.Infact*_dt*SoilWidthDX->Drc*DX->Drc;

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


    s.drain = -qbot;
    Perc->Drc = s.drain*_dt;

    double gw = 0;
    for (int j = nN; j > 0; j--) {
        if (s.h[j] >= 0)
            gw = gw + s.dz[j];
    }
    Lw->Drc = gw;

    if (i_ == cell) qDebug() << WH0 << s.InfPot << s.Infact  << s.ponded;
    if (i_==cell) {
        QString S;
        QString S1;
        for(int j = 0; j < nNodes; j++) {
            S = S + QString(" %1").arg(s.h[j]);
            S1 = S1 + QString(" %1").arg(s.theta[j]);
        }
        qDebug() << S;
        qDebug() << S1;
    }


    delete[] Hnew;
    delete[] Hold;
    delete[] C1;
    delete[] S;
    delete[] kavg;
    delete[] thoma;
    delete[] thomb;
    delete[] thomc;
    delete[] thomf;
    delete[] Wnew;
    delete[] disnod;
    delete[] z;
}



