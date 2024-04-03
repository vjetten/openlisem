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


double TWorld::calcSinkterm(long i_, double WH, double *S)
{
    SOIL_LIST s = crSoil[i_];
    int r = s.r;
    int c = s.c;
    //    if (i_ == 3738) qDebug() << day << hour << Ld;

    if (ETp->Drc*ETafactor >0 && Rain->Drc* 3600000.0/_dt > rainfallETa_threshold) {

        double AreaSoil = SoilWidthDX->Drc * DX->Drc;
        double Cover_ = Cover->Drc;
        double ETp_ = ETp->Drc * ETafactor; // potential ETp
        double tot = 0;
        double etanet = ETp_;
        double ETpshade = ETp_*(1-Cover_)+0.15*ETp_*(Cover_);
       //double eta = 0;

        ETpCum->Drc += ETp_;

       //  interception decrease, drying out canopy
        double CStor_  = CStor->Drc;
        if (CStor_ > 0) {
            double ETa_int = ETp_;

            ETa_int = std::min(ETa_int, CStor_);
            CStor_ = CStor_- ETa_int;

            RainCum->Drc = std::max(0.0, RainCum->Drc-ETa_int);
            if (CStor_ < 1e-5)
               RainCum->Drc = 0;

            // restart the cumulative process when CStor is dried out

            Interc->Drc = Cover_ * CStor_ * AreaSoil;
            IntercETa->Drc += Cover_ * ETa_int * AreaSoil;
            CStor->Drc = CStor_;
        }

        if (SwitchHouses)
        {
            double CvH = HouseCover->Drc;
            double HS = HStor->Drc;

            double ETa_int = std::min(ETp_, HS);
            HStor->Drc = HS - ETa_int;
            IntercETa->Drc += CvH * ETa_int * AreaSoil;
            double roofsurface = (_dx * DX->Drc * CvH); // m2
            IntercHouse->Drc =  roofsurface * HS;
        }

        // on ground level energy is shared
        if (SwitchLitter) {
            double CvL = Litter->Drc;
            double LCS = LCStor->Drc;

            double ETa_int = std::min(ETpshade, LCS);
            LCStor->Drc = LCS - ETa_int;
            IntercETa->Drc += CvL * ETa_int * CHAdjDX->Drc;
            LInterc->Drc =  CvL * LCS * CHAdjDX->Drc;
        }

        //transpiration under Cover from rootzone
        for (int j = 0; j < nNodes; j++) {
            // van genuchten H50 = -3.5 m
            double f = 1.0/(1.0+pow(s.h[j]/-3.5,1.5));
            if (s.h[j] > s.hb[j]) f = 0; // saturation
            if (s.h[j] < -16.0) f = 0; // wilting point -16000 cm
            etanet = ETp_* Cover_;
            S[j] =  etanet * f * s.rootz[j];
        }

        // add surface evaporation (1-Cover) to top node
        if (!s.ponded) {
           etanet = ETp_*(1-Cover_);
           if (s.h[0] > -16)
               S[0] += (s.theta[0]-s.thetar[0])/(s.pore[0]-s.thetar[0])*etanet;
        } else {
            WH = WH - ETpshade;
            WH = std::max(0.0, WH);
        }

        for (int j = 0; j < nNodes; j++) {
            tot += S[j];
        }

        ETa->Drc = tot;
        ETaCum->Drc += tot;


        // if (r==_nrRows/2 && c == _nrCols/2) {
        //     QString ss;
        //     for(int j = 0; j < nNodes; j++) {
        //         ss = ss + QString(" %1").arg(S[j]);
        //     }
        //     qDebug() << "S" << tot << etanet << ss;
        // }

    }
    return (WH);
}


void TWorld:: VanGenuchten(SOIL_LIST s, double Hnew[], double K[], double C1[], bool analytical)
{
    for(int j = 0; j < nNodes; j++) {
        double Se, Kr;
        double m = 1-1/s.vg_n[j];
        double Hx = std::min(Hnew[j], 0.0);

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
            double W1 = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hx-0.1), s.vg_n[j]), -m);
            C1[j] = (W-W1)/0.01;
        } else {
            if (Hnew[j] < s.hb[j])
                C1[j] = (W-s.thetar[j])*(s.vg_n[j]-1)*s.vg_alpha[j]*std::pow(fabs(s.vg_alpha[j]*fabs(Hx)),s.vg_n[j]-1)/
                            (1+std::pow(abs(s.vg_alpha[j]*fabs(Hx)),s.vg_n[j])) + W*1e-6/s.pore[j];
            else
                C1[j] = W*1e-6/s.pore[j];
        }
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
            double Wnew = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hx, s.lambda[j]);
            double W = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/(Hx-0.01), s.lambda[j]);
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

void TWorld::getThetafromH(int j, SOIL_LIST s)
{
    if (SwitchBrooksCorey) {
        if (s.h[j] < s.hb[j])
            s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
    } else {
        if (s.h[j] < 0.0)
            s.theta[j] = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(s.h[j]), s.vg_n[j]), -(1-1/s.vg_n[j]));
        else
            s.theta[j] = s.pore[j];
    }
}

void TWorld::getHfromTheta(int j, SOIL_LIST s)
{
    double se = (s.theta[j] - s.thetar[j])/(s.pore[j]-s.thetar[j]);
    if (SwitchBrooksCorey) {
        double hh = std::pow(se, (1.0/s.lambda[j]));
        s.h.replace(j,s.hb[j]/hh);
    } else {
        double n = s.vg_n[j];
        double m = 1-1/n;
        s.h.replace(j, -std::pow((std::pow(1/se,1/m)-1),1/n)/s.vg_alpha[j]);
        // kPa to m water
    }
}

void TWorld::cell_Soilwater(long i_)
{
    SOIL_LIST s = crSoil[i_];

    double dtmin = 0.01*_dt;
    double dtmax = std::min(_dt, SoilWBdtfactor);
    s.dts = dtmax;
    int NITMAX = 12;
    bool stopit = false;
    int c = s.c;
    int r = s.r;
    double WH0, WH1;
    int nN = nNodes-1;
    bool freeDrainage = !SwitchImpermeable;

    double *H = new double[nNodes];
    double *Hnew = new double[nNodes];
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
        WH0 = hmx->Drc;// flood in kin wave

    s.InfPot = 0;  // m/s
    s.Infact = 0;
    s.dtsum = 0;
    s.drain = 0;

/*
    if (SwitchGWflow) {
        double dif = std::max(0.0, SoilDepth2init->Drc - GWWH->Drc);
        //distance GW to surface
        if (dif < s.z[1]) {
            double store = dif*(s.pore[1]-s.theta[1]);
            if (store > WH1) {
                WH1 = 0;
                store = WH0-store;
                double m = (s.theta[1]-s.thetar[1])*dif + WH0;
                s.theta[1] = m/dif + s.thetar[1];
                getHfromTheta(1, s);
            } else {
                WH1 -= store;
                store = 0;
                s.theta[1] = s.pore[1];
                s.h[1] = 0;
            }

            if (FloodDomain->Drc == 0) {
                WH->Drc = WH1;
            } else {
                hmx->Drc = WH1;
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
            return;

        }
    }


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
*/

    // fill Hnew and Hold
    for(int j = 0; j < nNodes; j++) {
        H[j] = s.h[j];
        Hnew[j] = s.h[j];
        S[j] = 0;
    }

    // if (SwitchGWflow) {
    //     for(int j = GWnode+1; j < nNodes; j++) {
    //         Hnew[j] = s.z[j]-(s.SD-GWWH->Drc);
    //         Hold[j] = Hnew[j];
    //         // h = 0 at depth SD2-GWWH
    //         // h = GWWH at depth SD2
    //         // h = ? at depth s.z
    //         //TODO:if GW lowers and nodes become negative!
    //     }
    // }


    // ET calculation goes into sink term S
    // return reduced WH if ponding
    WH1 = WH0;
    if (SwitchIncludeET)
        WH1 = calcSinkterm(i_, WH0, S);

    int cnt = 0;

    // outer loop timestep lisem
    do {
        int NIT = 0; // iteration steps, max 12
        s.InfPot = WH1/_dt;  //potential infil is all water at the surface m/s
        double ANE, FFN, FNN2, A1;

        // iteration for Hnew
        do {

            // SwitchVanGenuchten = false;
            // SwitchBrooksCorey = true;
            //======== Get nodal values
            if(SwitchVanGenuchten)
                VanGenuchten(s, H, K, C1, false); //????Hnew???c
            if(SwitchBrooksCorey)
                BrooksCorey(s, H, K, C1, false);
            for(int j = 0; j < nNodes; j++) {
                C1[j] *= s.dz[j]/s.dts;

            //     if (SwitchGWflow && j >= GWnode)
            //         K[j] *= GW_recharge;
            }
            if (r == _nrRows/2 && c == _nrCols/2) {
                 QString S;
                 QString S1;
                 for(int j = 0; j < nNodes; j++) {
                     S = S + QString(" %1").arg(s.h[j]*100);
                    S1 = S1 + QString(" %1").arg(s.theta[j]);
                 }
                qDebug() << time/86400 << S;
                qDebug() << S1;
            }
            // K1 and K2 are avg between node and upper and lower node
            for(int j = 1; j < nNodes; j++) {
                switch (KavgType) {
                case 0: K1[j] = Aavg(K[j-1],K[j]); break;
                case 1: K1[j] = Savg(K[j-1],K[j]); break;
                case 2: K1[j] = Havg(K[j-1],K[j],s.dz[j-1],s.dz[j]); break;
                case 3: K1[j] = Mavg(K[j-1],K[j]); break;
                }
            }
            for(int j = 0; j < nNodes-1; j++) {
                switch (KavgType) {
                case 0: K2[j] = Aavg(K[j],K[j+1]); break;
                case 1: K2[j] = Savg(K[j],K[j+1]); break;
                case 2: K2[j] = Havg(K[j],K[j+1],s.dz[j],s.dz[j+1]); break;
                case 3: K2[j] = Mavg(K[j],K[j+1]); break;
                }
            }
            K1[0] = K2[0];
            K2[nN] = K1[nN];

            //======== Galerkin 3-diagonal scheme and back substitution
            // checked against WORM-V !!!

            // C is differential moisture capacity dtheta/dh. That means that F has the unit of m?
            // but is multiplied here with dz/dt !?
            F[0] = (2*C1[0] +C1[1])/6.0;
            A[0] = -K1[0];
            D[0] = -A[0] + F[0];
            for(int j = 1; j < nN; j++) {
                F[j] = (C1[j-1] + 4*C1[j] + C1[j+1])/6.0;
                A[j] = -K2[j];
                D[j] = K1[j] + K2[j] + F[j];
            }
            ANE = A[nN];
            F[nN] = (C1[nN] +2*C1[nN])/6.0;
            D[nN] = -A[nN-1] + F[nN];
            FFN = F[nN];

            //include gravity term (i.e. k only) and sink term
            //NOTE in WORM for first term originally -K1 + S1 but in paper
            //  by Van Genuchten -K1 -S1 !!!
            F[0] = F[0]*H[0] - 0.5*s.dz[0]*(K1[0]-(2.0*S[0]+S[1])/3.0);
            for(int j = 1; j < nN; j++) {
                F[j] = F[j]*H[j] + 0.5*(s.dz[j]+s.dz[j-1]) * K1[j]
                                 - 0.5*(s.dz[j]+s.dz[j+1]) * K2[j]
                                 - s.dz[j] * (S[j-1]+4*S[j]+S[j+1])/6;
            }
            F[nN] = F[nN]*H[nN] + 0.5*s.dz[nN]*(K2[nN]-(S[nN-1]+2.0*S[nN])/3.0);

            //======== Boundary conditions

            // check for ponding
            s.ponded = Hnew[0] > 0; //or H[0] > 0 ???? check

            // NOT necessary, SWATRE code:
            // double qmax = -0.5*(s.Ks[0]+K[0])*((Hnew[0]-WH1)/s.dz[0] + 1);
            // if (fabs(qmax) < fabs(s.InfPot))
            //     s.ponded = true;
            // if (r == _nrRows/2 && c == _nrCols/2)
            // qDebug() << "qmax" << s.ponded << qmax << s.InfPot;

            // upper boundary condition
            if (s.ponded){
                Hnew[0] = H[0] + WH1; //SOAP: throughfall -influx = the water layer on the surface, but what is influx? prev timestep?
                s.Infact = D[0] * H[0] - F[0];
                F[0] = Hnew[0];
                A1 = A[0];
                A[0] = 0;
                D[0] = 1;
                F[1] = F[1] - A1*Hnew[0]; //was +
                // THIS IS ALSO THE CONDITION UNDER GW, SO D=1!, A = 0, F = H???
            } else {
                s.Infact = s.InfPot;
                F[0] = F[0] + s.InfPot; // in WORM this is a flux
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
                double A2 = A[j-1]/D[j-1];
                D[j] = D[j] - A2 * A[j-1];
                F[j] = F[j] - A2 * F[j-1];
            }
            Hnew[nN] = F[nN]/D[nN];
            for (int j = nN-1; j >= 0; j--)
                Hnew[j] = (F[j] - A[j]*Hnew[j+1])/D[j];


//            for(int j = 1; j < nNodes; j++)
//                Hnew[j] = std::min(0.0, Hnew[j]);
            // not necessary? and surface can be + so not for the top node anyway!

            //======== calc boundary fluxes

            if (s.ponded && Hnew[1] < 0) {
               s.Infact = s.Infact + A1 * Hnew[1];
            }

            //CHECK!!!!!!!!!
            if (freeDrainage) {
                if (Hnew[nN] < 0)
                    s.drain = 0.5*(K[nN] + FNN2 - Hnew[nN]*0.5*C1[nN]);
                else
                    s.drain = 0.5*(K[nN] + FNN2);
            }

            s.ponded = Hnew[0] > 0;

            stopit = true;
            //stop if Hnew and Hold are close
            // do not include top node?
            for(int j = 1; j <= nN; j++) {
                if (Hnew[j] < 0) {
                    double tol = tol1*fabs(H[j]) + tol2;
                    if (s.ponded) tol /= 2;
                    if (fabs(Hnew[j] - H[j]) > tol) {
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
                  Hnew[j] = 0.5*(H[j]+Hnew[j]);
                stopit = true;
            }

            //?????????? do this or not?
            for(int j = 0; j < nNodes; j++)
                H[j] = Hnew[j];

        } while(!stopit);
            // iteration Hnew


        //========= change timestep for next iteration

        //        double factor = 0.5 + 1/sqrt((double)NIT);
        //        s.dts = s.dts * factor;

        // SWATRE method is faster!
        double dt = dtmax;
        for(int j = 1; j < nNodes; j++) {
            double mdih = tol2 + tol1*fabs(Hnew[j]);
            double dih  = fabs(Hnew[j] - H[j]);
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

       // if (r == _nrRows/2 && c == _nrCols/2) qDebug() << "iteration" <<  NIT << s.dts << s.dtsum << _dt;
        cnt++;

    } while(s.dtsum < _dt);
        // outer _dt loop

    for(int j = 0; j < nNodes; j++) {
        s.h[j] = Hnew[j];
        getThetafromH(j, s);
    }

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

    // if (r == _nrRows/2 && c == _nrCols/2) {
    //    // qDebug() << "timeloop" << WH0 - WH1 << s.Infact*3600000  << s.h[0] << s.h[1];
    //      QString S;
    //      QString S1;
    //      for(int j = 0; j < nNodes; j++) {
    //          S = S + QString(" %1").arg(s.h[j]*100);
    //         S1 = S1 + QString(" %1").arg(s.theta[j]);
    //      }
    //     qDebug() << time/86400 << S;
    //     qDebug() << S1;
    // }

    // put the results back
    crSoil[i_] = s;

    delete[] H;
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
    double *disnod  = new double[nNodes];


    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc; // flood in kin wave

    double WH1 = WH0;
    s.dtsum = 0;
    s.dts = SoilWBdtfactor;//*_dt;

    //#pragma omp parallel for num_threads(userCores)
    for(int j = 0; j < nNodes; j++) {
        Hnew[j] = s.h[j];
        Hold[j] = s.h[j];
        S[j] = 0;
        disnod[j] = j < nN ? s.z[j+1]-s.z[j] : s.dz[j];
    }

    //    if (SwitchIncludeET)
    //        calcSinkterm(i_, S); // ETa

    // outer loop timestep lisem
    do {

        SwitchVanGenuchten = true;
        //SwitchBrooksCorey = false;//true;
        bool anl = false;

        if(SwitchVanGenuchten)
            VanGenuchten(s, Hnew, K, C1, anl);
        else
        //if(SwitchBrooksCorey)
            BrooksCorey(s, Hnew, K, C1, anl);

//        if (r == _nrRows/2 && c == _nrCols/2) {
//            QString S;
//            for(int j = 0; j < nNodes; j++) {
//                S = S + QString(" %1").arg(Hnew[j]);
//            }
//            qDebug() << "zero" << S;
//        }
        for (int j = 1; j < nNodes; j++) {
        //    if (r==_nrRows/2 && c == _nrCols/2)
           //     qDebug()<<j << Hnew[j] << C1[j];

            switch (KavgType) {
                case 0: kavg[j] = Aavg(K[j],K[j-1]);break;
                case 1: kavg[j] = Savg(K[j],K[j-1]);break;
                case 2: kavg[j] = Havg(K[j],K[j-1],1.0,1.0);break;
                case 3: kavg[j] = Mavg(K[j],K[j-1]);break;
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

        // 1st check flux aginst max Darcy flux
        double qmax = -Savg(s.Ks[0],K[0])*((Hnew[0]-WH1) / s.dz[0] + 1);
        //qtop = -kavg[0] * ((h[0] - pond)/DistNode(p)[0] + 1);

        s.ponded = false;
        // maximum possible flux, compare to real top flux available
        if (fabs(qtop) > fabs(qmax))
            s.ponded = true;

//        double store = 0;
//        for (int j = 1; j < nNodes; j++) {
//            store = store + (s.pore[j] - s.theta[j])*s.dz[j];
//        }
//        if (WH1 > store)
//            s.ponded = true;

      //  if (r==_nrRows/2 && c == _nrCols/2)
       //     qDebug() << "top" << Hnew[0];//s.ponded << qtop << qmax << WH1 << K[0];


        if ( s.ponded ) //|| (fltsat && (qtop <= qbot)) )
        {
            thomc[0] = -s.dts * kavg[1]/ (s.dz[0]/disnod[1]);
            thomb[0] = -thomc[0] + C1[0] + s.dts*kavg[0]/(disnod[0]/s.dz[0]);
            thomf[0] = C1[0]*Hnew[0] - s.dts/s.dz[0]*(kavg[0]-kavg[1])
                       + s.dts*kavg[0]*WH1/(disnod[0]/s.dz[0]);
        }
        else
        {
            thomc[0] = -s.dts * kavg[1]/(s.dz[0]/disnod[1]);
            thomb[0] = -thomc[0] + C1[0];
            thomf[0] = C1[0]*Hnew[0] - s.dts/s.dz[0] * (fabs(qtop) - kavg[1]);
        }

        //Intermediate nodes: i = 1 to n-2
        for (int j = 1; j < nN; j++) {

            thoma[j] = -s.dts*kavg[j]/(s.dz[j]/disnod[j]);
            thomc[j] = -s.dts*kavg[j+1]/(s.dz[j]/disnod[j+1]);
            thomb[j] = -thoma[j] - thomc[j] + C1[j]; //D
            thomf[j] = C1[j]*Hnew[j] - s.dts/s.dz[j]*(kavg[j]-kavg[j+1]); //F
        }
        // last node
        thoma[nN] = -s.dts*kavg[nN]/(s.dz[nN]/disnod[nN]);
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

//        if (r == _nrRows/2 && c == _nrCols/2) {
//            QString S;
//            for(int j = 0; j < nNodes; j++) {
//                S = S + QString(" %1").arg(Hnew[j]);
//            }
//            qDebug() << "first" << S;
//        }

        //correct tridiagonal matrix
        for (int j = 0; j < nNodes; j++) {
            double Cnew = 1e-6;
            double Wnew = s.theta[j];
            double Wold = s.theta[j];
            if (SwitchVanGenuchten) {
                double Hx = std::min(Hnew[j], 0.0);
                double m = 1-1/s.vg_n[j];
                double Se = std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hx), s.vg_n[j]), -m);
                Wnew = s.thetar[j]+(s.pore[j]-s.thetar[j])*Se;
                double Wnew1 = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hx-0.01), s.vg_n[j]), -m);
                Cnew = (Wnew-Wnew1)/0.01;

                Wold = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hold[j]), s.vg_n[j]), -m);

            } else {
                double Hx = std::min(Hnew[j], s.hb[j]);
                Wnew = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/Hx, s.lambda[j]);
                double Wnew1 = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(s.hb[j]/(Hx-0.01), s.lambda[j]);
                Cnew = (Wnew-Wnew1)/0.01;
                Wold = s.thetar[j] + (s.pore[j]-s.thetar[j])*pow(std::min(1.0,s.hb[j]/Hold[j]), s.lambda[j]);
            }

            thomb[j] = thomb[j] - C1[j] + Cnew;
            thomf[j] = thomf[j] - C1[j]*Hold[j] + Cnew*Hnew[j]- Wnew + Wold;
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

//        if (r == _nrRows/2 && c == _nrCols/2) {
//            QString S;
//            for(int j = 0; j < nNodes; j++) {
//                S = S + QString(" %1").arg(Hnew[j]);
//            }
//            qDebug() << "2nd" << S;
//        }

        for (int j = 1; j < nNodes; j++)
            if(Hnew[j] > 0) Hnew[j] = 0;

        //====== headcalc

        // determine new boundary fluxes
        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = kavg[nN]*(Hnew[nN]-Hnew[nN-1])/disnod[nN] - kavg[nN];

        if ( s.ponded )       //|| (fltsat && (qtop < qbot)) )
            qtop = qmax;//-kavg[0] * ((Hnew[0] - WH1)/disnod[0] + 1);
        // okay: qtop is incoming water (-WH1/s.dts) unless ponded, then Darcy flux

        WH1 -= fabs(qtop)*s.dts;
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
        if (SwitchVanGenuchten) {
            if (Hnew[j] < -0.01)
                s.theta[j] = s.thetar[j]+(s.pore[j]-s.thetar[j])*std::pow(1+std::pow(s.vg_alpha[j]*fabs(Hnew[j]), s.vg_n[j]), -(1-1/s.vg_n[j]));
            else
                s.theta[j] = s.pore[j];
        } else {
            if (s.h[j] < s.hb[j])
                s.theta[j] = s.thetar[j] + pow(s.hb[j]/s.h[j], s.lambda[j])*(s.pore[j]-s.thetar[j]);
            else
                s.theta[j] = s.pore[j];
        }

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
        qDebug() << WH0 << WH1 << s.Infact;
        QString S;
        QString S1;
        for(int j = 0; j < nNodes; j++) {
            S = S + QString(" %1").arg(s.h[j]);
            S1 = S1 + QString(" %1").arg(s.theta[j]);
        }
       qDebug() << S;
     //    qDebug() << S1;
    }

    // put the results back
    crSoil[i_] = s;

    delete[] Hnew;
    delete[] Hold;
    delete[] C1;
    delete[] S;
    delete[] kavg;
    delete[] thoma;
    delete[] thomb;
    delete[] thomc;
    delete[] thomf;
    delete[] disnod;   
}



/*
void TWorld::cell_SoilwaterExpl(long i_)
{
    SOIL_LIST s = crSoil[i_];

    double dtmin = 0.01*_dt;
    double dtmax = SoilWBdtfactor*_dt;
    s.dts = _dt;//dtmax;
    int c = s.c;
    int r = s.r;
    double WH0, WH1;
    bool freeDrainage = !SwitchImpermeable;

    double *S  = new double[nNodes];
    double *K  = new double[nNodes];
    double *C1  = new double[nNodes];
    double *Ka = new double[nNodes];
    double *Q  = new double[nNodes];
    double *Hnew  = new double[nNodes];

    // SwitchVanGenuchten = false;
    // SwitchBrooksCorey = !SwitchVanGenuchten;

    if (FloodDomain->Drc == 0)
        WH0 = WH->Drc; //runoff in kinwave or dyn wave
    else
        WH0 = hmx->Drc;// flood in kin wave

    WH1 = WH0;

    s.InfPot = 0;  // m/s
    s.Infact = 0;
    s.dtsum = 0;
    s.drain = 0;

    // fill Hnew and Hold
    for(int j = 0; j < nNodes; j++) {
        S[j] = 0;
        Hnew[j] = s.h[j];
    }

    //get K and C for current H
    if(SwitchVanGenuchten)
        VanGenuchten(s, Hnew, K, C1, false);
    if(SwitchBrooksCorey)
        BrooksCorey(s, Hnew, K, C1, false);

    // ET calculation goes into sink term S
    // if (SwitchIncludeET)
    //     calcSinkterm(i_, S); // ETa

    // initial dts
    s.dts=dtmax;
    int cnt = 0;

    // K1 and K2 are avg between node and upper and lower node
    for(int j = 1; j < nNodes; j++) {
        switch (KavgType) {
        case 0: Ka[j] = Aavg(K[j-1],K[j]); break;
        case 1: Ka[j] = Savg(K[j-1],K[j]); break;
        case 2: Ka[j] = Havg(K[j-1],K[j],s.dz[j-1],s.dz[j]); break;
        case 3: Ka[j] = Mavg(K[j-1],K[j]); break;
        }
    }
    Ka[0] = Ka[1];
    s.InfPot = WH1/s.dts;

    // explicit state and fluxes
    // ponding?
    // Q=-K[h]*((h1-h2)/(z1-z2)+1) but dz is z2-z1 so reversed and positive
    // therefore in this case it should be Q=-K[h]*((h2-h1)/dz+1), where h <= 0 or on the surface can be positive
    if (WH1 > 0) {
        Q[0] = -s.Ks[0]*((s.h[0]-WH1)/s.dz[0] + 1);
        Q[0] = std::max(Q[0], -s.InfPot);
    } else {
        Q[0] = -Ka[0]*((s.h[0]-0)/s.dz[0] + 1); //??? is there something better
    }
    // if (r == _nrRows/2 && c == _nrCols/2)
    //     qDebug() << WH1 << Q[0] << s.h[0] << s.theta[0] << s.dz[0] << Ka[0] << K[0] << K[1];

    for(int j = 1; j < nNodes-1; j++) {
        Q[j] = -Ka[j]*((s.h[j]-s.h[j-1])/s.dz[j] + 1);
    }

    for(int j = 0; j < nNodes; j++) {
        getThetafromH(j,s);
        double moist = s.theta[j] * s.dz[j];
        double moistmax = s.pore[j] * s.dz[j];
        double dm = moistmax-moist;
        if ((Q[j]-Q[j+1])*s.dts > dm) {
            Q[j] = (dm+Q[j+1]*s.dts)/s.dts;
        }
        moist = moist + (Q[j]-Q[j+1])*s.dts;
        s.theta[j] = moist/s.dz[j];
        s.theta[j] = std::min(s.theta[j],s.pore[j]);

        getHfromTheta(j,s);
    }

    s.drain = Q[nNodes-1];
    if (WH1 > 0)
        WH1 = std::max(0.0, WH1 + Q[0]*s.dts);

    // outer loop timestep lisem

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

    // if (r == _nrRows/2 && c == _nrCols/2) {
    //  //   qDebug() << WH0 << WH1 << Q[0];
    //     QString S;

    //     // for(int j = 0; j < nNodes; j++) {
    //     //     S = S + QString(" %1").arg(s.h[j]);
    //     //     //S1 = S1 + QString(" %1").arg(s.theta[j]);
    //     // }
    //    // qDebug() <<  S;
    // }

    // put the results back
    crSoil[i_] = s;

    delete[] Hnew;
    delete[] C1;
    delete[] S;
    delete[] Ka;
    delete[] K;
    delete[] Q;

}
*/
