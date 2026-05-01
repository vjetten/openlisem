/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/


#include "model.h"


//---------------------------------------------------------------------------
/**
 * @fn double TWorld::MaxConcentration(double watvol, double sedvol)
 * @brief Calculates concentration with a maximum of MAXCONC, changes sed vol and deposisition
 *
 * @param watvol : the watervolume
 * @param sedvol : the sediment mass
 * @return sediment concentration (kg/m3)
 * @see MAXCONC
 *
 */
double TWorld::MaxConcentration(double watvol, double sedvol)
{
    double conc = 0;//MAXCONC;//0;
    if (watvol > 1e-6) {
        conc = qMin(sedvol/watvol, MAXCONC);
    }
    return conc;
}

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::GetSV(double d)
 * @brief get settling velocity of sediment with grain size d
 *
 * @param d : the grain size (in micrometer)
 * @return The settling velocity
 */
double TWorld::GetSV(double d )
{
    if (SwitchSV == 2) {
        double dm = d / 1e6;
        double ds = dm * pow(1.65*GRAV/1e-12,0.333333);
        return SVCHCalibration*1e-6/dm*(ds*ds*ds)*pow(38.1+0.93*pow(ds,12.0/7.0), -7.0/8.0);
        //    // zhiyao et al, 2008
    } else {
        if(d < 100) {
            return SVCHCalibration*2*(2650.0-1000.0)*GRAV*pow(d/2000000.0, 2)/(9*0.001);
            //Stokes range settling velocity
        } else {
            double dm = d/1000.0;
            return SVCHCalibration*10.0 *sqrt(1.0 + 0.01 *(1.65* GRAV *dm*dm*dm )-1.0)/dm;
            //Settling velocity by Zanke (1977)
        }
    }

}

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::calcTCSuspended(int r,int c, int method,double U, int type)
 * @brief Calculates suspended layer transport capacity
 *
 * Calculates suspended load sediment transport capacity.
 * Based on govers, Van Rijn (simplified or full)
 * or Wu wang & Jia (for multiclass sediment).
 *
 * @param _d : The grain class (only needed when grain size distribution is used)
 * @param method : the TC method used
 * @param U : velicty can be channel of overland or flood
 * @param type : channel (0) or flood (1) or overland (2)
 */

double TWorld::calcTCSuspended(int r,int c, int method, double h, double w, double U, int type)
{
    double R=0, hs=0, S = 0, man = 0.01;
    double d50m;

    if (type == SUSPchannel) {
        // river
        d50m = D50CH->Drc/1000000.0;
        hs = ChannelSSDepth->Drc;
        S = ChannelGrad->Drc;
        R = (w*h)/(2*h+w);
        man = ChannelN->Drc;
    } else
        if (type == SUSPflood) {
            // flood
            d50m = D50->Drc/1000000.0;
            hs = SSDepthFlood->Drc;
            S = Grad->Drc;
            //R = (w*h)/(2*h+w); //?
            R = h;
            man = N->Drc;
        } else
            if (type == SUSPrunoff) {
                // kin wave
                hs = WHrunoff->Drc;
                S = Grad->Drc;
                R = WHrunoff->Drc;
                man = N->Drc;
            }

    //when water height is insignificant, transport capacity is zero
    //this is necessary since some of the used equations have strange behaviour
    //for these water heights or velocities. (h and v outside of valid range)
    if(h < MIN_HEIGHT || hs < MIN_HEIGHT)
        return 0;
    if(U < MIN_FLUX)
        return 0;

    double ps = 2650.0;
    double pw = 1000.0;
    double tc = 0;

    if(method == FSHAIRSINEROSE)
    {
        double om =  U*S;
        double omcr = 0.004;
        tc =  d50m/SettlingVelocitySS->Drc* 0.013/GRAV * 1.650 * qMax(0.0, om - omcr)/h ;
        //    m/ (m/s)* kg* m/s /m  dimensionless?

    } else
        if(method == FSGOVERS)
        {
            //### Calc transport capacity
            double uc = 100.0*U*S; //in cm/s  in this formula assuming grad is SINE
            double ucr = 0.4;   // critical unit streampower in cm/s
            double cg = cgovers->Drc;//pow((d50m+5)/0.32, -0.6);
            double dg = dgovers->Drc;//pow((d50m+5)/300, 0.25);
            tc = ps * cg * pow(qMax(0.0, uc-ucr), dg); // kg/m3

        } else
            if(method == FSRIJN)
            {
                //https://www.leovanrijn-sediment.com/papers/Formulaesandtransport.pdf
                //double kinvis = 1e-6;
                //double Ds = d50m * pow(1.65*GRAV/1e-12,1.0/3.0); //dimensionless sed size
                //cr,suspension = 0.3/(1+ D*) + 0.1 [1-exp(-0.05D*)
                // critical shields parameter
                //double cs = 0.3/(1+Ds)+0.1*(1-exp(-0.05*Ds));
                // Ucritical, suspension= 5.75 [log(12h/(6D50))] [cr,suspension (s-1) g D50]0.5

                // double ucr;
               // ucr = 5.75*(log10(12*hs/(6.0*d50m))*qSqrt(cs*1.650*GRAV*d50m); //1.650 = s-1
//                if( d50m < 0.0005) // 500 mu, so always the first one!
//                    ucr = 0.19 * pow(d50m, 0.1) * log10(2.0* h/d50m); //p5
//                else
//                    ucr = 8.5  * pow(d50m, 0.6) * log10(2.0* h/d50m);
//                // set in motion critical U
                double Ds = d50m * 25296; //pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0)); // let op: D* = 25296*D50m! R2 = 1
              //  double shields = 0.3/(1+Ds)+0.1*(1-exp(-0.05*DS));
              //  double ucr = 5.75*log10(2*h/d50m)*sqrt(shields*1.65*GRAV*d50m);

                double ucr = UcrCHCalibration*2.8*pow(h/d50m,0.1)*sqrt(1.65*GRAV*d50m);
//                // suspension critical U

               //page 5

                double me = qMax(0.0,U - ucr)/sqrt(GRAV * d50m * 1.65);
                //p15 mobility parameter

                double qs = 0.03 * ps*U*d50m * me*me * pow(Ds, -0.6); // kg/s/m
                // van rijn 2007?, p 17, eq 6.4
                ChannelQsr->Drc = qs;
                tc =  qs/ (U * h); //kg/s/m / (m2/s) =  kg/m3   => WH or WHs

            }else
                if(method == FSRIJNFULL)
                {
                    /*
                 * D50	d50m	D*
                10	0.00001	0.25
                30	0.00003	0.76
                50	0.00005	1.26
                100	0.0001	2.53
                300	0.0003	7.59
                500	0.0005	12.65
                1000	0.001	25.30
                2000	0.002	50.59
                */
                    //van Rijn full (1984) following page 1632
                    double kinvis = 1e-6;
                    double ds = d50m * pow(1.65*GRAV/(kinvis*kinvis),(1.0/3.0)); //
                    //double chezy = 18 * log10(4 * R/d90m);
                    double chezy = 1/man*pow(R,1/6);
                    double uc = U * sqrt(GRAV)/chezy;

                    //shields functions
                    double uscr = 0.055;
                    if(ds <150 && ds >= 20)
                        uscr = 0.013*pow(ds,0.29);
                    if(ds < 20 && ds >= 10)
                        uscr = 0.04*pow(ds,-0.10);
                    if(ds < 10 && ds >= 4)
                        uscr = 0.14*pow(ds,-0.64);
                    if(ds <4)
                        uscr = 0.24*pow(ds,-1);
                    uscr = UcrCHCalibration*sqrt(uscr * 1.65*GRAV * d50m);

                    double T = qMax(((uc*uc)/(uscr*uscr) - 1),0.0);  //transport stage parameter
                    double bsv = sqrt(GRAV * h * S); // bed shear velocity
                    double a = 0.1;  // half of the bedform height in m
                    double ca = 0.015 * (d50m/a) * pow(T,1.5)/pow(ds,0.3); //eq 38 reference concentration
                    double sv = SettlingVelocitySS->Drc;//GetSV(D50->Drc);

                    double beta = qMin(1.0 + 2.0*(sv/bsv)*(sv/bsv),5.0);
                    double powcb = 0.75; // not clear, between 0.4 and 1
                    double phi = 2.5 * pow(sv/bsv,0.8) * powcb;
                    double Z = sv/(beta*bsv*0.41); //suspension parameter, to do with upward turbulent versus gravity
                    double Zs = Z + phi;
                    double ad = 0.1; // else F is not valid!
                    double F = (pow(ad,Zs) - pow(ad,1.2))/(pow(1.0-ad,Zs)* (1.2 - Zs));
                    double qs =  F * U * h * ca;
                    tc = ps * qs/ (U * h);
                } else
                    if(method == FENGELUND)
                    {
                        //https://www.hec.usace.army.mil/confluence/rasdocs/rassed1d/1d-sediment-transport-technical-reference-manual/computing-transport-capacity/sediment-transport-potential/engelund-hansen
                        double U2 = U*U;
                        double C = 1/man*pow(R,1/6);  //Chezy
                        //double Tb = pw*U2*GRAV/(C*C);
                        //double shields = Tb/((ps-pw)*GRAV*d50m);
                        double shields = U2/(C*C*1.65*d50m);
                        double qs = U2*pow(shields, 1.5)*sqrt(d50m/(GRAV*1.65));

                        //http://ponce.sdsu.edu/onlineengelundhansen.php
                        //double shields1 = h*S/(1.65*d50m);
                        //double qs = 0.001*(U2/(2*GRAV*S*h))*pow(shields1,2.5)*dw*sqrt(1.65*GRAV*d50m*d50m*d50m);
                        //0.001 is ton to kg
                        // gives almost the same value
                        ChannelQsr->Drc = qs;

                        tc =  qs/ (U * h); //kg/s/m / (m2/s) =  kg/m3   => WH or WHs
                } else
                        if(method == FSWUWANGJIA) {

                        /*
                        // NOT USED, FOR MULTIPLE GRAINSIZES

                    double phk = 0;
                    double pek = 0;
                    double sv = settlingvelocities.at(_d);
                    double gd = graindiameters.at(_d)/1000000.0;
                    if (type == 0) {
                        FOR_GRAIN_CLASSES
                        {
                            //LET OP : RW_D and W_D !!!!
                            phk += RW_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                            pek += RW_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
                        }
                    } else {
                        FOR_GRAIN_CLASSES
                        {
                            //LET OP : RW_D and W_D !!!!
                            phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                            pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
                        }
                    }

                    double ppk = 1;
                    ppk = pow(phk/pek,0.6);
                    if(pek == 0 )
                    {
                        return 0;
                    }

                    double css = 0.03* (ps - pw) * (gd) * ppk;

                    double qs = 0.0000262 *pow(qMax(( pw * 0.01 * h * GRAV * S /css) - 1.0, 0.0)* U/(sqrt(sv)),2.2);
                    qs = qs * 1 * sqrt((ps/pw - 1)*GRAV*pow(gd,3.0));

                    tc = ps * qs/ (U * h);
*/
                }
    return qMax(qMin(tc,MAXCONC ),0.0);
}
//--------------------------------------------------------------------------
/**
 * @fn double TWorld::calcTCBedload(int r,int c, int _d, int method, bool river)
 * @brief Calculates suspended layer transport capacity
 *
 * Calculates suspended load sediment transport capacity.
 * Based on govers, Van Rijn (simplified or full)
 * or Wu wang & Jia (for multiclass sediment).
 *
 * @param _d : The grain class (only needed when grain size distribution is used)
 * @param method : the TC method used
 */
double TWorld::calcTCBedload(int r,int c, int method, double h, double w, double U, int type)
{
    double R,  hb, n, S;

    if (type == 0) {
        //    h = ChannelWH->Drc;
        hb = ChannelBLDepth->Drc;
        n = qMax(0.001, ChannelN->Drc);
        S = ChannelGrad->Drc;
        //w = ChannelWidth->Drc;
        R = (w*h)/(2*h+w);
    } else
        if (type == 1) {
            //   h = hmx->Drc;
            hb = BLDepthFlood->Drc;
            n = qMax(0.001, N->Drc);
            S = Grad->Drc;
           // w = ChannelAdj->Drc*rillfactor;
            R = (w*h)/(2*h+w);
        }

    //when water height is insignificant, transport capacity is zero
    //this is necessary since some of the used equations have strange behaviour
    //for these water heights or velocities. (h and v outside of valid range)
    if(h < MIN_HEIGHT || hb < MIN_HEIGHT)
        return 0;
    if(U < MIN_FLUX)
        return 0;

    double ps = 2650.0; //2400.0;
    double pw = 1000.0;
    double d50m = (D50->Drc/1000000.0);
    double d90m = (D90->Drc/1000000.0);
    if (type == 0) {
        d90m = D90CH->Drc/1000000.0;
        d50m = D50CH->Drc/1000000.0;
    }

    double tc = 0;

    if(method == FSRIJN)
    {
        //Van rijn simplified (2007?)
        double ucr;
        if( d50m < 0.0005)
            ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0*R/d90m);
        else
            ucr  = 8.5 * pow(d50m, 0.6) * log10(4.0*R/d90m);

        double me = qMax((U - ucr)/(sqrt(GRAV * d50m * ((ps/pw) - 1.0))),0.0);
        //        double qs = 0.005 * ps * U * h * pow(d50m/h,1.2) * pow(me, 2.4);
        double qs = 0.015 * ps*U*h * pow(d50m/h,1.2) * pow(me, 1.5); //eq 6.2
        // in kg/m/s /(m2/s) = kg/m3
        tc =  qs/ (U * hb);

    }else if(method == FSRIJNFULL)
    {
        //van Rijn full (1984)  see page 1450 1984_JHE_VanRijn_a.pdf
        double kinvis = 1e-6;

        double _dm = d90m; //d50m; interpretation -> assume all bedload particles are d90? Van RIjn deals mostly with sand

        double ds = _dm * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
        double chezy = 18 * log(4 * h/d90m);  // h or hb or radius?
        double us = sqrt(GRAV) * U/chezy;

        // shield equations, full
        double uscr = 0.055;
        if(ds < 150 && ds >= 20)
            uscr = 0.013*pow(ds,0.29);
        if(ds < 20 && ds >= 10)
            uscr = 0.04*pow(ds,-0.10);
        if(ds < 10 && ds >= 4)
            uscr = 0.14*pow(ds,-0.64);
        if(ds <4)
            uscr = 0.24*pow(ds,-1);
        uscr = sqrt(uscr * (ps/pw - 1)*GRAV * _dm);  // effective bed shear velocity

        double T = qMax((us*us)/(uscr*uscr) - 1,0.0); // transport stage parameter
        double qs = 0.053 * (pow(T,2.1)/pow(ds,0.3)) * sqrt((ps/pw -1)*GRAV)*_dm*sqrt(_dm); // eq 22
        tc = ps * qs/ (U * hb);

    }else if(method == FSWUWANGJIA)
    {
        /*
        double na = (pow(graindiameters.at(_d)/100000.0,(1.0/6.0))/20.0)/n;
        double phk = 0;
        double pek = 0;
        if (type == 0) {
            FOR_GRAIN_CLASSES
            {
                //LET OP : RW_D and W_D !!!!
                phk += RW_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += RW_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }
        } else {
            FOR_GRAIN_CLASSES
            {
                //LET OP : RW_D and W_D !!!!
                phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }
        }
        double ppk = 1;
        ppk = pow(phk/pek,0.6);

        if(pek == 0 )
            return 0;

        double R = w*h/(2*h+w);
        double css = 0.03* (ps - pw) * (graindiameters.at(_d)/1000000.0) * ppk;

        double qs = 0.0053 *pow(qMax(pow(na,1.5)*((pw * R * GRAV * 0.1 * S/css)) - 1.0, 0.0),2.2);
        qs = qs * 1 * sqrt((ps/pw - 1)*GRAV*pow(graindiameters.at(_d)/1000000.0,3.0));

        tc = ps * qs/ (U * hb);
*/
    }

    return qMax(qMin(tc,MAXCONCBL),0.0);
}

//---------------------------------------------------------------------------
// NOT USED FOR NOW
/**
 * @fn double TWorld::DetachMaterial(int r,int c, int d,bool channel, bool flood,bool bl,double detachment)
 * @brief Calculates real detachment from potential detachment.
 *
 * This cell uses the real time calculated effective erosion coefficient
 * to calculate actual erosion. When the soil layer has less sediment left then
 * the potential erosion, an analytical solution is used to converge
 * both sediment in the soil layer and sediment in tranport to a
 * stable value. When both the channel and flood parameter are false,
 * overland flow detachment is assumed.
 *
 * @param r : Row nr of the cell
 * @param c : Column nr of the cell
 * @param d : Grain diameter class
 * @param Channel : Channel detachment?
 * @param flood : Flood detachment?
 * @param bl : Bed Load detachment?
 * @param detachment : Potential detachment
 * @return Actual detachment
 */

double TWorld::DetachMaterial(int r,int c, int d,bool channel, bool flood,bool bl,double detachment)
{
    /*
     * NOTE: the actual detachment is immediately taken
     * from the soil layers. It is therefore assumed that when using
     * this function, the actual detachment is added to the
     * sediment in flow. THIS IS REQUIRED! to maintain mass
     * balance
     */
    // when there is no usage of material depth, there is no deposited layer
    // actual erosion can then be calculated using the original erosion efficiency coefficient
    if(!SwitchUseMaterialDepth)
    {
        if(channel)
            return detachment *= ChannelY->Drc;
        else
            return detachment *= Y->Drc;
    }

    //first check if it is channel detachment
    if(channel)
    {
        //calculate depth of deposited layer
        double depdepth = qMax((RStorageDep->Drc / (BulkDens))/(ChannelWidth->Drc * DX->Drc),0.0);

        //linear decrease in influence from lower soil layer
        //from 0 to MixingDepth, with fac1 for bottom layer, fac2 for top layer
        double fac1 = 1.0;
        if(RSedimentMixingDepth->Drc > MIN_HEIGHT)
            fac1 = qMax(0.0,1.0 - depdepth/RSedimentMixingDepth->Drc);
        double fac2 = 1-fac1;

        //new erosion coefficient bases on soil layer mixinfactors
        double newY = ChannelY->Drc * fac1 + fac2 * 1.0;

        //multiply potential detachment by erosion coefficient
        detachment = detachment *newY;

        //to remove small rounding errors that lead to negative values
        detachment = qMax(detachment,0.0);

        //check wat we can detache from the top and bottom layer of present material
        double dleft = detachment;
        double deptake = 0;
        double mattake = 0;
        detachment = 0;

        //take from the total storage
        deptake = qMin(dleft,RStorageDep->Drc);
        RStorageDep->Drc -= deptake;

        //add to the detachment what we have taken from the first soil layer
        detachment += deptake;

        //if the deposited layer is empty
        //use erosion efficiency of bottom layer again
        if(newY > 0)
            dleft *= ChannelY->Drc/newY;
        else
            dleft = 0;

        //bottom soil layer can be infinite
        if(!((RStorage->Drc) < -1))
        {
            //take from the total storage
            mattake = qMin(dleft,RStorage->Drc);
            RStorage->Drc -= mattake;
            //add to the detachment what we have taken from the second soil layer
            detachment += mattake;
        } else {
            //all left potential detachment is added to detachment (infinite soil layer)
            detachment += dleft;
        }

        //finally, return the total detachment
        return qMax(0.0,detachment);

    } else
        if(flood) {
            double depdepth = qMax((StorageDep->Drc / (BulkDens))/(_dx * DX->Drc),0.0);

            //linear decrease in influence from lower soil layer
            //from 0 to MixingDepth, with fac1 for bottom layer, fac2 for top layer
            double fac1 = 1.0;
           // if(SedimentMixingDepth->Drc > MIN_HEIGHT)
                fac1 = qMax(0.0,1.0 - depdepth/SedimentMixingDepth->Drc);
            double fac2 = 1-fac1;
            //new erosion coefficient bases on soil layer mixinfactors
            double newY = Y->Drc * fac1 + fac2 * 1.0;

            //multiply potential detachment by erosion coefficient
            detachment = detachment *newY;

            //to remove small rounding errors that lead to negative values
            detachment = qMax(detachment,0.0);

            //check wat we can detach from the top and bottom layer of present material
            double dleft = detachment;
            double deptake = 0;
            double mattake = 0;
            detachment = 0;

            //take from the total storage
            deptake = qMin(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;

            //add to the detachment what we have taken from the first soil layer
            detachment += deptake;

            //if the deposited layer is empty
            //use erosion efficiency of bottom layer again
            if(newY > 0)
                dleft *= Y->Drc/newY;
            else
                dleft = 0;

            if(!((Storage->Drc) < -1))
            {
                //take from the total storage
                mattake = qMin(dleft,Storage->Drc);
                Storage->Drc -= mattake;
                //add to the detachment what we have taken from the second soil layer
                detachment += mattake;
            } else {
                //all left potential detachment is added to detachment (infinite soil layer)
                detachment += dleft;
            }

            //finally, return the total detachment
            return qMax(0.0,detachment);

            //if it is neither flood nor channel detachment, overland flow is assumed
        } else {
            //calculate depth of deposited layer
            double depdepth = qMax((StorageDep->Drc / (BulkDens))/(_dx * DX->Drc),0.0);

            //linear decrease in influence from lower soil layer
            //from 0 to MixingDepth, with fac1 for bottom layer, fac2 for top layer
            double fac1 = qMax(0.0,1.0 - depdepth/SedimentMixingDepth->Drc);
            double fac2 = 1 - fac1;
            if(SedimentMixingDepth->Drc < MIN_HEIGHT)
            {
                fac1 = 1;
                fac2 = 0;
            }
            //new erosion coefficient bases on soil layer mixinfactors
            double newY = Y->Drc * fac1 + fac2 * 1.0;
            //multiply potential detachment by erosion coefficient

            detachment = detachment *newY;

            //to remove small rounding errors that lead to negative values
            detachment = qMax(detachment,0.0);
            //check wat we can detache from the top and bottom layer of present material
            double dleft = detachment;
            double deptake = 0;
            double mattake = 0;
            detachment = 0;

            //take from the total storage
            deptake = qMin(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;

            //add to the detachment what we have taken from the first soil layer
            detachment += deptake;
            //if the deposited layer is empty
            //use erosion efficiency of bottom layer again
            if(newY > 0)
                dleft *= Y->Drc/newY;
            else
                dleft = 0;

            //bottom soil layer can be infinite
            if(!((Storage->Drc) < -1))
            {
                //take from the total storage
                mattake = qMin(dleft,Storage->Drc);
                Storage->Drc -= mattake;
                //add to the detachment what we have taken from the second soil layer
                detachment += mattake;
            } else {
                //all left potential detachment is added to detachment (infinite soil layer)
                detachment += dleft;
            }

            //finally, return the total detachment
            return qMax(0.0,detachment);
        }
}
