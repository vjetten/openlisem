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
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

#include "lisemqt.h"
#include "model.h"

// redistribution of soilwater after infiltration
// out[put is new Lw and new Thetaeff and ThetaI2
void TWorld::cell_Redistribution1(int r, int c)
{
    double Percolation, pore, theta, thetar, theta_E;
    double SoilDep2, pore2, theta2, thetar2, FC, FC2;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    pore = Poreeff->Drc;
    thetar = 0.025 * pore;
    theta = Thetaeff->Drc;
    FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;

    int ww = 0;

    Percolation = 0;

    double factor = 0.1;
    if (SwitchTwoLayer) {
        pore2 = ThetaS2->Drc;
        thetar2 = 0.025 * pore2;
        theta2 = ThetaI2->Drc;
        SoilDep2 = SoilDepth2->Drc;
        FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore2;
    }

    if (Lw_ < 0.001)
       return;
    // nothing to redistribute

    if (fact->Drc > 0)
        return;
    // no redistribution while there is infiltration

    if (SwitchImpermeable) {
        if (SwitchTwoLayer && Lw_ > SoilDepth2->Drc-0.01)
            return;
        if (!SwitchTwoLayer && Lw_ > SoilDepth1->Drc-0.01)
            return;
    }

//    if (r == 200 && c == 200)
//        qDebug() << "a" << fact->Drc << theta << theta2 << Lw_;
//if (Lw_ < 0) {
//    qDebug() << "a" << theta << theta2 << Lw_;
//}
    if(SwitchTwoLayer) {
        if (Lw_ < SoilDep1) {

            double store = (SoilDep1 - Lw_) * (pore-theta); // space in SD1 under Lw_

            if (store == 0) {
                Lw_ = SoilDep1;
                ww = -1;
            } else {
                theta_E = (theta-thetar)/(pore-thetar);
                Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc); // m/timestep
              //  Percolation = sqrt(Percolation * (Ksateff->Drc*_dt/3600000.0));
                //flux across boundary is ks+ke/2, average after Swatre
                double moistw = Lw_ * (pore-thetar); //available sat moisture above Lw_
                double dm = (pore-FC)*Lw_ * factor;
                Percolation = std::min(dm, Percolation);

                moistw -= Percolation; // decrease moistw with percolation, can be 0

                Lw_ = moistw/(pore-thetar); // new Lw_

                if (Percolation <= store) {
                    // if percolation fits in store layer 1
                    double m1 = theta*(SoilDep1-Lw_) + Percolation;
                    theta = m1/(SoilDep1-Lw_);

                    if (theta >= pore) {
                        theta = pore;
                        Lw_= SoilDep1;
                    }
                    ww = 1;
                } else {

                    //double Perc1 = store;
                    //double Perc2 = std::max(0.0,Percolation - Perc1);
                    // some spills over in layer 2
                    // Lw_ is in layer 1
                    double m1 = theta*(SoilDep1-Lw_);
                    double m2 = theta2*(SoilDep2-SoilDep1);
                    double Perc1 = m1/(m1+m2)*Percolation;
                    double Perc2 = m2/(m1+m2)*Percolation;

                    m1 = theta*(SoilDep1-Lw_) + Perc1;
                    theta = m1/(SoilDep1-Lw_);

                    if (theta >= pore) {
                        theta = pore;
                        Lw_ = SoilDep1;
                    }

                    m2 = theta2*(SoilDep2-SoilDep1) + Perc2;
                    theta2 = m2/(SoilDep2-SoilDep1);

                    if (theta2 >= pore2) {
                        theta2 = pore2;
                        Lw_ = SoilDep2;
                    }
                    ww = 2;
                }
            }
        } else {
            //Lw_ > SoilDep1

            // get the potential flux in m
            theta_E = (theta2-thetar2)/(pore2-thetar2);
            Percolation = Ksat2->Drc*_dt/3600000.0 * pow(theta_E, bca2->Drc); // m/timestep
           //Percolation = sqrt(Percolation * (Ksat2->Drc*_dt/3600000.0));
            //flux across boundary is (ks+ke)/2, average after Swatre

            double moist1 = SoilDep1 * (pore - thetar); // moisture above wettingfront in layer 1
            double moist2 = (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2
            //available moisture
            double dm1 = (pore - FC)*SoilDep1*factor;
            double dm2 = (pore2 - FC2)*(Lw_-SoilDep1)*factor;
            //available moisture to move

            Percolation = std::min(Percolation, dm1+dm2);
            //cannot have more flux than available water

            double perc2 = std::min(Percolation, dm2); // part taken from layer 2
            double perc1 = std::max(0.0, Percolation-perc2); // part taken from layer 1, can be 0

            moist2 -= perc2;
            double m2 = (SoilDep2-Lw_)*theta2 + perc2;
            Lw_ = moist2/(pore2-thetar2) + SoilDep1;
            theta2 = m2/(SoilDep2-Lw_);
            theta = pore;

            if (theta2 >= pore2) {
                Lw_ = SoilDep2;
                theta2 = pore2;
            }

            ww = 3;

            if (perc1 > 0) {
                // opercolation moves into layer 1
                moist1 -= perc1;

                // sat layer 1 now becomes free until FC
                Lw_ = moist1/(pore-thetar);
                theta = FC;// (FC)/(SoilDep1-Lw_);

//                if (theta >= pore) {
//                    theta=pore;
//                    Lw_ = SoilDep1;
//                }
                ww = 4;
            }

        }
    }
//    if (r == 200 && c == 200)
//        qDebug() << "b" << ww << theta << theta2 << Lw_;

    Thetaeff->Drc = theta;
    ThetaI2->Drc = theta2;
    Lw->Drc = Lw_;//std::max(0.0, Lw_);

}
//---------------------------------------------------------------------------
// redistribution of water from the wetting front to the layer below. Lw decreases and thetai increases
// percolation is done after this
void TWorld::cell_Redistribution(int r, int c)
{
    double Percolation, pore, theta, thetar, theta_E;
    double SoilDep2, pore2, theta2, thetar2, moistw, m1, m2;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    pore = Poreeff->Drc;
    thetar = 0.025 * pore;
    theta = Thetaeff->Drc;
    Percolation = 0;

    double factor = 0.2;

    //double FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;

    if (Lw_ == 0)
       return;
    // nothing to redistribute

    if (WH->Drc >0 || hmx->Drc > 0)
        return;
    // no redistribution while there is water on the surface

    if (SwitchImpermeable) {
        if (SwitchTwoLayer && Lw_ > SoilDepth2->Drc-0.01)
            return;
        if (!SwitchTwoLayer && Lw_ > SoilDepth1->Drc-0.01)
            return;
    }

    if (SwitchTwoLayer) {

        pore2 = ThetaS2->Drc;
        thetar2 = 0.025 * pore2;
        theta2 = ThetaI2->Drc;
        SoilDep2 = SoilDepth2->Drc;
    }

    if (SwitchTwoLayer && Lw_ < SoilDep1) {
        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc); // m/timestep
       // Percolation = sqrt(Percolation * (Ksateff->Drc*_dt/3600000.0));
        //flux across boundary is ks+ke/2, average after Swatre

        moistw = Lw_ * (pore-thetar); //available sat moisture above Lw_
        Percolation = std::min(moistw*factor, Percolation);

        m1 = theta * (SoilDep1-Lw_);
        m2 = theta * (SoilDep2-SoilDep1);

     }

     if (SwitchTwoLayer && Lw_ >= SoilDep1) {

        // get the potential flux in m
        theta_E = (theta2-thetar2)/(pore2-thetar2);
        Percolation = Ksat2->Drc*_dt/3600000.0 * pow(theta_E, bca2->Drc); // m/timestep
      //  Percolation = sqrt(Percolation * (Ksat2->Drc*_dt/3600000.0));
        //flux across boundary is (ks+ke)/2, average after Swatre

        //get the available moisture
        moistw = SoilDep1 * (pore - thetar); // moisture above wettingfront in layer 1
        moistw += (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2

        m1 = 0;
        m2 = theta2 * (SoilDep2-Lw_);
     }

     Percolation = std::min(Percolation, moistw*factor);
        //cannot have more flux than total water*factor

     Lw_ = moistw/(0.5*((pore-thetar)+(pore2-thetar2)));

     if (Lw_ >= SoilDep2) {
        m2 += Percolation;
        theta2 = m2/(SoilDep2-Lw_);
     } else {
        double p1 = (SoilDep1-Lw_)/(SoilDep2) * Percolation;
        double p2 = (1-(SoilDep1-Lw_)/(SoilDep2)) * Percolation;

        m1 += p1;
        m2 += p2;

        theta = m1/(SoilDep1-Lw_);
        theta2 = m2/(SoilDep2-SoilDep1);
     }

    Thetaeff->Drc = theta;
    ThetaI2->Drc = theta2;
    Lw->Drc = Lw_;
}

//---------------------------------------------------------------------------
// percolation from the bottom of the soil profile
double TWorld::cell_Percolation(int r, int c, double factor)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    if(SwitchTwoLayer) {
        pore = ThetaS2->Drc;
        thetar = 0.025 * pore;
        theta = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;
        double ksat = factor*Ksat2->Drc*_dt/3600000;

        if(theta > thetar) {
            // percolation in m per timestep
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, bca2->Drc);

            if (Lw->Drc > SoilDepth1->Drc)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // if Wet Fr still in first layer percolation only make 2nd drier

            double moisture = std::max(0.0, dL*(theta-thetar));
            if (Lw_ < SoilDep2-0.001) {
                // decrease thetaeff because of percolation
                Percolation = std::min(Percolation, moisture*0.5);
                moisture -= Percolation;
                theta = moisture/dL + thetar;
                theta = std::max(thetar, theta);
            } else {
                // wetting front = soildepth2, dL = 0, moisture = 0
                // assume theta goes back to 0.7 pore and decrease the wetting fornt
                Percolation = sqrt(ksat * Percolation);  //k = sqrt(ks*k)
                double FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore;
                theta = FC2;
                double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - Percolation/(pore - theta));
                Percolation = (Lwo-Lw_)*(pore-FC2);
            }            
            ThetaI2->Drc = theta;
            Lw->Drc = Lw_;
            //Perc->Drc = Percolation;
            return(Percolation);
        }
    } else {
        // one layer
        double pore = Poreeff->Drc;
        thetar = 0.025 * pore;
        double theta = Thetaeff->Drc;
        double ksat = factor*Ksateff->Drc*_dt/3600000;

        if (theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, bca1->Drc);
            if (Lw_ < SoilDep1) {
                // wetting front has not reached bottom, make soil drier
                // decrease thetaeff because of percolation
                dL = SoilDep1 - Lw_;
                double moisture = dL*(theta-thetar);
                Percolation = std::min(Percolation, moisture);
                moisture -= Percolation;
                theta = moisture/dL + thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume tehta goes back to 0.7 pore and decrease the wetting fornt
                double FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;
                theta = FC;
                double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - Percolation/(pore - theta));
                Percolation = (Lwo-Lw_)*(pore-FC);
            }
            Thetaeff->Drc = theta;
//            Perc->Drc = Percolation;
            Lw->Drc = Lw_;
            return(Percolation);
        }
    }
    return(0);
}

//---------------------------------------------------------------------------
void TWorld::cell_Interception(int r, int c)
{
    // all variables are in m
    double Cv = Cover->Drc;
    double Rainc_ = Rainc->Drc;
    double RainCum_ = RainCum->Drc;
    double AreaSoil = SoilWidthDX->Drc * DX->Drc;
    double RainNet_ = Rainc_;

    if (Cv > 0)
    {
        double CS = CStor->Drc;
        //actual canopy storage in m
        double Smax = CanopyStorage->Drc;
        //max canopy storage in m

        if (Smax > 0) {
            CS = Smax*(1-exp(-kLAI->Drc*RainCum_/Smax));
        }

        LeafDrain->Drc = std::max(0.0, Cv*(Rainc_ - (CS - CStor->Drc)));
        // diff between new and old strage is subtracted from rainfall
        // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!
        // note: cover already implicit in LAI and Smax, part falling on LAI is cover*rainfall

//        if(r==100&&c==200)
//            qDebug()<< CS << Smax << RainCum_ << exp(-kLAI->Drc*RainCum_/Smax);
        CStor->Drc = CS;
        // put new storage back in map
        Interc->Drc =  Cv * CS * AreaSoil;
        // Interc->Drc =  Cv * CS * _dx * DX->Drc;
        // WHY: cvover already takes care of this, trees can be above a road or channel

        RainNet_ = LeafDrain->Drc + (1-Cv)*Rainc_;
        // net rainfall is direct rainfall + drainage
        // rainfall that falls on the soil, used in infiltration
    }

    if (SwitchLitter) {
        double CvL = Litter->Drc;
        if (hmx->Drc == 0 && WH->Drc == 0 && CvL > 0 && RainNet_ > 0)
        {

            double Smax = LitterSmax/1000.0;
            // assume simply that the cover linearly scales between 0 and LtterSmax of storage

            double LCS = LCStor->Drc;
            //actual canopy storage in m

            LCS = std::min(LCS + RainNet_, Smax);
            // add water to the storage, not more than max

            double drain = std::max(0.0, CvL*(RainNet_ - (LCS - LCStor->Drc)));
            // diff between new and old storage is subtracted from leafdrip

            LCStor->Drc = LCS;
            // put new storage back in map for next dt

            LInterc->Drc =  CvL * LCS * AreaSoil;
            // only on soil surface, not channels or roads, in m3

            RainNet_ = drain + (1-CvL)*RainNet_;
            //recalc
        }
    }

    // all variables are in m
    if (SwitchHouses)
    {
        double CvH = HouseCover->Drc;
        if (CvH > 0 &&  RainNet_ > 0)
        {
            //house on cell in m2
            double HS = HStor->Drc;
            //actual roof storage in m

            double Hmax = RoofStore->Drc;
            //max roof storage in m

            HS = std::min(HS + RainNet_, Hmax);

            double housedrain = std::max(0.0, CvH * (RainNet_ - (HS - HStor->Drc)));
            // overflow in m3/m2 of house

            HStor->Drc = HS;
            // put new storage back in maps in m

            double roofsurface = (AreaSoil * CvH); // m2
            // double roofsurface = (_dx * DX->Drc * CvH); // m2
            // user should assure housecover is correct with respect to channel and roads
            IntercHouse->Drc =  roofsurface * HS;
            // total interception in m3,exclude roads, channels

            RainNet_ = housedrain + (1-CvH)*RainNet_;
            // net rainfall is direct rainfall + drainage
            // rainfall that falls on the soil, used in infiltration

            // filling raindrums with surplus drainage from roofs
            // drum is recalculated to m based on roof surface
            double DS = 0;
            if (SwitchRaindrum && DrumStore->Drc > 0)
            {
                double Dmax = DrumStore->Drc/roofsurface;
                //max drum storage in m as if roof storage is more

                DS = DStor->Drc;
                //actual drum storage in m

                DS = std::min(DS + RainNet_, Dmax);
                // fill tank to max
                double drumdrain = std::max(0.0, HouseCover->Drc * (RainNet_ - (DS - DStor->Drc)));

                DStor->Drc = DS;
                // put new drum storage back in maps in m3

                IntercHouse->Drc += roofsurface * DS;
                // total interception in m3,exclude roads, channels

                RainNet_ = drumdrain + (1-CvH)*RainNet_;
            }
        }
    }

    RainNet->Drc = RainNet_;

}

//---------------------------------------------------------------------------

