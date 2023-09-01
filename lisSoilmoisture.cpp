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

//DO NOT USE OMP, MUCH SLOWER
// Calculates one dimensional soil water balance based on
// Richards equation, using a fully implicit, mass lumped,
// Galerkin scheme Finite Element Method.
// Based on FORTRAN program WORM-V by Van Genuchten (1987),
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

// i: 1 = theta from h, 2 k from h, 3 C from h, 4 h from theta



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

        if (r==_nrRows/2 && c == _nrCols/2) {
            QString ss;
            for(int j = 0; j < nNodes; j++) {
                ss = ss + QString(" %1").arg(S[j]);
            }
            //   qDebug() << tot << etanet << ss;
        }

    }
}


void TWorld:: VanGenuchten(SOIL_LIST s, double Hnew[], double K[], double C1[], bool analytical)
{
    for(int j = 0; j < nNodes; j++) {
        double Se, Kr;
        double m = 1-1/s.vg_n[j];
        double Hx = std::min(Hnew[j]*100, s.hb[j]*100);

        Se = std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hx), s.vg_n[j]), -m);

        if (Se == 1.0)
            Kr = 1.0;
        else
            if (Se < 0.04)
                Kr = qSqrt(Se) * std::pow(m * std::pow(Se,1/m), 2.0);
            else
                Kr = qSqrt(Se) * std::pow(1-std::pow(1-std::pow(Se,1/m),m), 2.0);
        K[j] = s.Ks[j]*Kr;

        double W = s.thetar[j]+(s.pore[j]-s.thetar[j])*Se;
        if (!analytical) {
            double W1 = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hx-0.01), s.vg_n[j]), -m);
            C1[j] = (W-W1)/0.01;
        } else {
            if (Hnew[j] < s.hb[j])
                C1[j] = (W-s.thetar[j])*(s.vg_n[j]-1)*s.vg_alpha[j]*std::pow(fabs(s.vg_alpha[j]*fabs(Hx)),s.vg_n[j]-1)/
                            (1+std::pow(abs(s.vg_alpha[j]*fabs(Hx)),s.vg_n[j])) + W*1e-6/s.pore[j];
            else
                C1[j] = W*1e-6/s.pore[j];
        }
      //  if (s.r==_nrRows/2 && s.c == _nrCols/2)
         //   qDebug()<<j << Hnew[j] << C1[j];
    }
}

void TWorld::BrooksCorey(SOIL_LIST s, double Hnew[], double K[], double C1[], bool analytical)
{
    for(int j = 0; j < nNodes; j++) {
        double Hx = Hnew[j];
        if (Hx < s.hb[j])
            K[j] = s.Ks[j]*pow(std::min(1.0,s.hb[j]/Hx), 2.0+3.0*s.lambda[j]);
        else
            K[j] = s.Ks[j];

        // differential moisture capacity dtheta/dh (tangent of pF curve)
        if (!analytical) {
            double Hx = std::min(Hx, s.hb[j]);
            double Wnew = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hnew[j], s.lambda[j]);
            double W = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/(Hnew[j]-0.01), s.lambda[j]);
            C1[j] = (Wnew-W)/0.01;
        } else {
           if (Hnew[j] < s.hb[j])
                C1[j] = s.theta[j] * 1e-6/s.pore[j]
                        -1.0/Hnew[j] *(s.pore[j]-s.thetar[j])*s.lambda[j]*pow(s.hb[j]/Hnew[j], s.lambda[j]);
            else
                C1[j] = s.theta[j] * 1e-6/s.pore[j];
        }
    }
}


void TWorld::cell_Soilwater(long i_)
{
    SOIL_LIST s;
    s = crSoil[i_];

    double dtmin = 0.01*_dt;
    double dtmax = SoilWBdtfactor*_dt;
    s.dts = dtmax;
    int NITMAX = 12;
    bool stopit = false;
    int c = s.c;
    int r = s.r;
    double WH0, WH1;
    int nN = nNodes-1;
    bool freeDrainage = !SwitchImpermeable;

    double ANE, FNN1, FNN2, A1, A2;
    double qmax = 0;

    double *Hold = new double[nNodes];
    double *Hnew = new double[nNodes];
    double *C1 = new double[nNodes];
    double *S  = new double[nNodes];
    double *K  = new double[nNodes];
    double *K1 = new double[nNodes];
    double *K2 = new double[nNodes];
    double *A  = new double[nNodes];
    double *D  = new double[nNodes];
    double *F  = new double[nNodes];

    // find node above groundwater
    int GWnode = nN;
    if (SwitchGWflow) {
        if (GWWH->Drc > 0) {
            int j = 0;
            for (j = 1; j < nNodes; j++)
            {
                if (s.z[j] >= s.SD - GWWH->Drc)
                    break;
            }
            GWnode = j-1;
        }
    }
    nN = GWnode;

    //    if (r==_nrRows/2 && c == _nrCols/2)
    //        qDebug() << GWnode << GWWH->Drc << s.z[GWnode];

    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc;// flood in kin wave
    WH1 = WH0;
    s.InfPot = WH0/_dt;  // m/s
    s.Infact = WH0/_dt;
    s.dtsum = 0;
    s.drain = 0;

    // fill Hnew and Hold
    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.h[j];
        Hold[j] = Hnew[j];
        S[j] = 0;
    }


    if (SwitchGWflow) {
        for(int j = GWnode+1; j < nNodes; j++) {
            Hnew[j] = s.z[j]-(s.SD-GWWH->Drc);
            Hold[j] = Hnew[j];
            // h = 0 at depth SD2-GWWH
            // h = GWWH at depth SD2
            // h = ? at depth s.z
            //TODO:if GW lowers and nodes become negative!
        }
    }


    // ET calculation goes into sink term S
    if (SwitchIncludeET)
        calcSinkterm(i_, S); // ETa

    // initial dts
    s.dts=dtmax;
    int cnt = 0;

    // outer loop timestep lisem
    do {
        int NIT = 0;
        s.InfPot = WH1/s.dts;  // m/s

        // iteration for Hnew
        do {
            SwitchVanGenuchten = true;//false;
            SwitchBrooksCorey = false;//true;

            if(SwitchVanGenuchten)
                VanGenuchten(s, Hnew, K, C1, false);
            if(SwitchBrooksCorey)
                BrooksCorey(s, Hnew, K, C1, false);

            for(int j = 0; j < nNodes; j++) {

//                if (r == _nrRows/2 && c == _nrCols/2)
//                    qDebug() <<"vg" << Hnew[j] << K[j] << C1[j];

                C1[j] *= s.dz[j]/s.dts;

                if (SwitchGWflow && j >= GWnode)
                    K[j] *= GW_recharge;
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

            //======== Galerkin 3-diagonal scheme and back substitution

            A[0] = -K1[0];
            D[0] = K1[0]+ 0.5*C1[0];
            F[0] = 0.5*C1[0]*Hnew[0] - s.dz[0]*(K1[0] + (2*S[0]+S[1])/6);
            for(int j = 1; j < nN; j++) {
                A[j] = -K2[j];
                D[j] = K1[j] + K2[j] + C1[j];
                F[j] = C1[j]*Hnew[j]
                       + 0.5*(s.dz[j]+s.dz[j-1]) * K1[j]
                       - 0.5*(s.dz[j]+s.dz[j+1]) * K2[j]
                       - s.dz[j] * (S[j-1]+4*S[j]+S[j+1])/6;
            }
            A[nN] = -K2[nN];
            D[nN] = K2[nN] + 0.5*C1[nN];
            F[nN] = 0.5*C1[nN]*Hnew[nN] + s.dz[nN]*(K2[nN] - (S[nN-1]+2*S[nN])/6);


            s.ponded = Hnew[0] > 0;

            // check for ponding and first estimate of infil with darcy
            // with conductivity between Ksat and first node average K1
            qmax = Savg(s.Ks[0],K[0])*((WH1-Hnew[0])/s.dz[0] - 1);
            //qmax = -s.Ks[0]*((Hnew[0]-WH1)/s.dz[0] + 1);
            if (s.InfPot > 0 && fabs(qmax) < fabs(s.InfPot))
                s.ponded = true;

//            if (r == _nrRows/2 && c == _nrCols/2)
//                qDebug() <<"pond" << WH0 << WH1 << s.InfPot << s.Infact << qmax;

            // upper boundary condition
            if (s.ponded){

          //     Hnew[0] = Hold[0] + s.dts*(s.InfPot-qmax); //s.dts*std::min(s.InfPot,qmax);//
          //     Hnew[0] = Hold[0] + s.dts*std::min(s.InfPot,qmax);
                Hnew[0] = WH1;
                s.Infact = fabs(qmax);//D[0]*Hnew[0] - F[0];
                F[0] = F[0] + s.Infact;
                F[0] = Hnew[0];

                A[0] = 0;//-s.Ks[0];
                D[0] = 1;
                F[1] = F[1] + K1[0]*Hnew[0];
            } else {

                s.Infact = s.InfPot; // infil is net rainfall
                F[0] = F[0] + s.Infact;
                            if (r == _nrRows/2 && c == _nrCols/2)
                                qDebug() <<"pond" << F[0];
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
            for (int j = nN-1; j >= 0; j--)
                Hnew[j] = (F[j] - A[j]*Hnew[j+1])/D[j];


        //    for(int j = 1; j < nNodes; j++)
        //        Hnew[j] = std::min(0.0, Hnew[j]);
            // not necessary, and surface can be + so not for the top node anyway!

            //======== calc boundary fluxes
            if (s.ponded && Hnew[1] < 0) {
                s.Infact = s.Infact + A1 * Hnew[1];
            }
            if (freeDrainage) {
                s.drain = 0.5*(K[nN] + FNN2 - Hnew[nN]*0.5*C1[nN]);
            }

            stopit = true;
            bool iterate = false;
            if (iterate) {
               // s.ponded = Hnew[0] > 0;

                //stop if Hnew and Hold are close
                // do not include top node?
                for(int j = 0; j <= nN; j++) {
                    if (Hnew[j] < 0) {
                        double tol = tol1*fabs(Hold[j]) + tol2;
                        //if (s.ponded) tol /= 2;
                        if (fabs(Hnew[j] - Hold[j]) > tol) {
                            stopit = false;
                            break;
                        }
                    }
                }

                NIT++;

                if (!stopit && NIT > NITMAX) {
                    // try again with smaller dts
                    s.dts /= 2.0;
                    s.dts = std::max(s.dts,dtmin);
                    for(int j = 1; j < nNodes; j++)
                    //    Hnew[j] = Hold[j];
                      Hnew[j] = 0.5*(Hold[j]+Hnew[j]);
                    stopit = true;
                }
            }

            for(int j = 0; j < nNodes; j++)
                Hold[j] = Hnew[j];

        } while(!stopit);


        // change timestep for next iteration
        //        double factor = 0.5 + 1/sqrt((double)NIT);
        //        s.dts = s.dts * factor;
        //        s.dts = std::min(s.dts,dtmax);
        //        s.dts = std::max(s.dts,dtmin);
        //        s.dts = std::min(s.dts,_dt-s.dtsum);
        //        s.dtsum += s.dts;

        // SWATRE method is faster!
        double dt = dtmax;
        for(int j = 1; j < nNodes; j++) {
            double mdih = tol2 + tol1*fabs(Hnew[j]);
            double dih  = fabs(Hnew[j] - Hold[j]);
            if (dih > 0.10)
                dt = dt*mdih/dih;
        }
        s.dts = std::min(dt,dtmax);
        s.dts = std::max(dt,dtmin);
        s.dts = std::min(s.dts,_dt-s.dtsum);
        s.dtsum += s.dts;


        if (WH1 > 0) {
            WH1 = WH1 - s.Infact*s.dts;
            WH1 = std::max(0.0, WH1);
        }
        s.InfPot = WH1/s.dts;

        //if (r == _nrRows/2 && c == _nrCols/2)qDebug() << NIT << s.dts << s.dtsum << _dt;
        cnt++;
    } while(s.dtsum < _dt);

    for(int j = 0; j < nNodes; j++) {
        s.h[j] = Hnew[j];
//        if (s.h[j] < s.hb[j])
//            s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
//        else
//            s.theta[j] = s.pore[j];
        if (Hnew[j] < 0.0)
            s.theta[j] = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hnew[j]*100), s.vg_n[j]), -(1-1/s.vg_n[j]));
        else
            s.theta[j] = s.pore[j];

    }

    //    double infil = s.Infact*_dt;
    //    if (FloodDomain->Drc == 0) {
    //        infil = std::min(WH->Drc , infil);
    //        WH->Drc -= infil; //runoff in kinwave or dyn wave
    //    } else {
    //        infil = std::min(hmx->Drc , infil);
    //        hmx->Drc -= infil; // flood in kin wave
    //    }
    //    s.Infact = infil/_dt;

    if (FloodDomain->Drc == 0) {
        WH->Drc = WH1; //runoff in kinwave or dyn wave
    } else {
        hmx->Drc = WH1; // flood in kin wave
    }
    s.Infact = (WH0-WH1)/_dt;

    InfilVol->Drc = s.Infact*_dt*SoilWidthDX->Drc*DX->Drc;

    ThetaI1a->Drc = s.theta[1];
    for (int j = 2; j <= nN1_; j++)
        ThetaI1a->Drc += s.theta[j];
    ThetaI1a->Drc /= (double)nN1_;

    if (SwitchTwoLayer || SwitchThreeLayer) {
        ThetaI2a->Drc = s.theta[nN1_+1];
        for (int j = nN1_+2; j < nN1_+nN2_+1; j++)
            ThetaI2a->Drc += s.theta[j];
        ThetaI2a->Drc /= (double)nN2_;
    }
    if (SwitchThreeLayer) {
        ThetaI3a->Drc = s.theta[nN1_+nN2_+1];
        for (int j = nN2_+nN1_+2; j < nN1_+nN2_+nN3_+1; j++)
            ThetaI3a->Drc += s.theta[j];
        ThetaI3a->Drc /= (double)nN3_;
    }

    Perc->Drc = s.drain*_dt;

    if (r == _nrRows/2 && c == _nrCols/2) {
        qDebug() << WH0 << WH1 << s.Infact*_dt  << s.drain;
        QString S;
        QString S1;
        for(int j = 0; j < nNodes; j++) {
            S = S + QString(" %1").arg(s.h[j]);
            //S1 = S1 + QString(" %1").arg(s.theta[j]);
        }
        qDebug() << cnt << S;
        //qDebug() << S1;
    }

    delete[] Hold;
    delete[] Hnew;
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


    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc; // flood in kin wave

    double WH1 = WH0;
    s.dtsum = 0;
    s.dts = SoilWBdtfactor*_dt;

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

        VanGenuchten(s, Hnew, K, C1, false);

        if (r == _nrRows/2 && c == _nrCols/2) {
            QString S;
            for(int j = 0; j < nNodes; j++) {
                S = S + QString(" %1").arg(Hnew[j]);
            }
            qDebug() << "zero" << S;
        }
        for (int j = 1; j < nNodes; j++) {
        //    if (r==_nrRows/2 && c == _nrCols/2)
           //     qDebug()<<j << Hnew[j] << C1[j];

            switch (KavgType) {
                case 0: kavg[j] = Aavg(K[j],K[j-1]);
                case 1: kavg[j] = Savg(K[j],K[j-1]);
                case 2: kavg[j] = Havg(K[j],K[j-1],1.0,1.0);
                case 3: kavg[j] = Mavg(K[j],K[j-1]);
            }
        }
        kavg[0] = kavg[1];

        for(int j = 0; j < nNodes; j++)
            Hold[j] = Hnew[j];

        //--- boundary conditions ---//

        //----- TOP -----//
        qtop = -WH1/s.dts;

        //----- BOTTOM -----//
        // bottom is 0 or copy of flux of last 2 layers
        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = kavg[nN]*(Hnew[nN]-Hnew[nN-1])/s.dz[nN] - kavg[nN];

        // 1st check flux aginst max flux
        double qmax = Savg(s.Ks[0],K[0])*(WH1-Hnew[0]) / s.dz[0] - kavg[0];
       // double qmax = -Savg(s.Ks[0],K[1])*((WH1-Hnew[1]) / s.dz[0] + 1);

        // maximum possible flux, compare to real top flux available
        if (fabs(qtop) > fabs(qmax))
            s.ponded = true;

      //  if (r==_nrRows/2 && c == _nrCols/2)
       //     qDebug() << "top" << Hnew[0];//s.ponded << qtop << qmax << WH1 << K[0];



        if ( s.ponded ) //|| (fltsat && (qtop <= qbot)) )
        {
            thomc[0] = -s.dts * kavg[1]/ s.dz[0]; // / s.dz[1];
            thomb[0] = -thomc[0] + C1[0] + s.dts*kavg[0]/s.dz[0]; // / s.dz[0];
            thomf[0] = C1[0]*Hnew[0]
                       - s.dts/s.dz[0]*(kavg[0] - kavg[1])
                       + s.dts*kavg[0]*WH1/s.dz[0]; // / s.dz[0];
        }
        else
        {
            s.ponded = false;
            thomc[0] = -s.dts * kavg[1]/s.dz[0]; // / s.dz[1];
            thomb[0] = -thomc[0] + C1[0];
            thomf[0] = C1[0]*Hnew[0] - s.dts/s.dz[0] * (fabs(qtop) - kavg[1]);
        }

        //Intermediate nodes: i = 1 to n-2
        for (int j = 1; j < nN; j++) {

            thoma[j] = -s.dts*kavg[j]/s.dz[j]; //s.dz[j];
                thomc[j] = -s.dts*kavg[j+1]/s.dz[j]; //s.dz[j+1];
            thomb[j] = -thoma[j] - thomc[j] + C1[j]; //D
            thomf[j] = C1[j]*Hnew[j] - s.dts/s.dz[j]*(kavg[j]-kavg[j+1]); //F
        }
        // last node
        thoma[nN] = -s.dts*kavg[nN]/s.dz[nN]; //s.dz[nN];
        thomb[nN]= -thoma[nN] + C1[nN];
        thomf[nN] = C1[nN]*Hnew[nN] - s.dts/s.dz[nN]*(kavg[nN]+qbot); //qbot = negative


        //Gaussian elimination and backsubstitution h - first time
        alpha = thomb[0];
        Hnew[0] = thomf[0] / alpha;
        for (int j = 1; j < nNodes; j++) {
            beta[j] = thomc[j-1] / alpha;
            alpha = thomb[j] - thoma[j] * beta[j];
            Hnew[j] = (thomf[j] - thoma[j] * Hnew[j-1]) / alpha;
        }
        for (int j = (nN-1); j >= 0; j--)
            Hnew[j] -= beta[j+1] * Hnew[j+1];

        if (r == _nrRows/2 && c == _nrCols/2) {
            QString S;
            for(int j = 0; j < nNodes; j++) {
                S = S + QString(" %1").arg(Hnew[j]);
            }
            qDebug() << "first" << S;
        }

        //correct tridiagonal matrix
        for (int j = 0; j < nNodes; j++) {
            //double Hx = std::min(Hnew[j],s.hb[j]);

            double Se;
            double m = 1-1/s.vg_n[j];
            double Hx = std::min(Hnew[j]*100, 0.0);

            Se = std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hx), s.vg_n[j]), -m);

            double W = s.thetar[j]+(s.pore[j]-s.thetar[j])*Se;
            double W1 = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hx-0.01), s.vg_n[j]), -m);
            double Cnew = (W-W1)/0.01;

            thomb[j] = thomb[j] - C1[j] + Cnew;
            thomf[j] = thomf[j] - C1[j]*Hold[j] + Cnew*Hnew[j]- W + s.theta[j];
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

        if (r == _nrRows/2 && c == _nrCols/2) {
            QString S;
            for(int j = 0; j < nNodes; j++) {
                S = S + QString(" %1").arg(Hnew[j]);
            }
            qDebug() << "2nd" << S;
        }




        for (int j = 1; j < nNodes; j++)
            if(Hnew[j] > 0) Hnew[j] = 0;

        //====== headcalc

        // determine new boundary fluxes
        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = kavg[nN]*(Hnew[nN]-Hnew[nN-1])/s.dz[nN] - kavg[nN];

//        if ( s.ponded )//|| (fltsat && (qtop < qbot)) )
//            qtop = -kavg[0] * ((Hnew[0] - WH1)/s.dz[0] + 1);
        // adjust top flux

        WH1 += qtop*s.dts;
        WH1 = std::max(0.0,WH1);
        // decrease pond with top flux

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

    } while(s.dtsum < _dt);

    for(int j = 0; j < nNodes; j++) {
        s.h[j] = Hnew[j];
//        if (s.h[j] < s.hb[j])
//            s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
//        else
//            s.theta[j] = s.pore[j];

        if (Hnew[j] < 0.0)
            s.theta[j] = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hnew[j]*100), s.vg_n[j]), -(1-1/s.vg_n[j]));
        else
            s.theta[j] = s.pore[j];
    }

    if (FloodDomain->Drc == 0) {
        WH->Drc = WH1; //runoff in kinwave or dyn wave
    } else {
        hmx->Drc = WH1;//s.Infact; // flood in kin wave
    }
    s.Infact = (WH0-WH1)/_dt;
    s.drain = -qbot;

    InfilVol->Drc = s.Infact*_dt*SoilWidthDX->Drc*DX->Drc;

    ThetaI1a->Drc = s.theta[1];
    for (int j = 2; j <= nN1_; j++)
        ThetaI1a->Drc += s.theta[j];
    ThetaI1a->Drc /= (double)nN1_;

    if (SwitchTwoLayer || SwitchThreeLayer) {
        ThetaI2a->Drc = s.theta[nN1_+1];
        for (int j = nN1_+2; j < nN1_+nN2_+1; j++)
            ThetaI2a->Drc += s.theta[j];
        ThetaI2a->Drc /= (double)nN2_;
    }
    if (SwitchThreeLayer) {
        ThetaI3a->Drc = s.theta[nN1_+nN2_+1];
        for (int j = nN2_+nN1_+2; j < nN1_+nN2_+nN3_+1; j++)
            ThetaI3a->Drc += s.theta[j];
        ThetaI3a->Drc /= (double)nN3_;
    }

    Perc->Drc = s.drain*_dt;

    if (r == _nrRows/2 && c == _nrCols/2) {
      //  qDebug() << WH0 << WH1 << s.Infact  << s.drain;
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
   // delete[] Wnew;
   // delete[] disnod;
   // delete[] z;
}



