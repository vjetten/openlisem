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

void TWorld::cell_Redistribution1(int r, int c)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    Percolation = 0;


    if(SwitchTwoLayer) {
        theta = (Thetaeff->Drc + ThetaI2->Drc)/2.0;
        pore = (Poreeff->Drc + ThetaS2->Drc)/2.0;
        double ksa = 0.5*(Ksateff->Drc+Ksat2->Drc)*_dt/3600000.0;

        theta_E = theta/pore;//(theta-thetar)/(pore-thetar);
        Percolation = ksa * pow(theta_E, 0.5*(bca1->Drc+bca2->Drc)); // m/timestep
        Percolation = 0.5*(Percolation + ksa);
        // assumed average percolation based on all parameters

        double SoilDep2 = SoilDepth2->Drc;

        dL = std::min(Lw_, Percolation/(pore-theta));
        Percolation = dL*(pore-theta);
        Lw_ = std::max(0.0, Lw_-dL);

        double m1 = (SoilDep1 - Lw_) * theta; // moisture above wettingfront
        double m2 = (SoilDep2 - SoilDep1) * ThetaI2->Drc; // moisture above wettingfront

        m1 += m1/(m1+m2)*Percolation;
        m2 += m2/(m1+m2)*Percolation;

        Thetaeff->Drc = m1/(SoilDep1-Lw_);
        ThetaI2->Drc = m2/(SoilDep2-SoilDep1);
        Lw->Drc = Lw_;
    }


}
//---------------------------------------------------------------------------
// redistribution of water from the wetting front to the layer below. Lw decreases and thetai increases
// percolation is done after this
void TWorld::cell_Redistribution(int r, int c)
{
    double Percolation, pore, theta, thetar, theta_E;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    pore = Poreeff->Drc;
    thetar = 0.025 * pore;
    theta = Thetaeff->Drc;
    Percolation = 0;
    double factor = 0.2;
    double dth = (pore-thetar) * factor;
    double FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;

    if (Lw_ == 0)
       return;

    if (SwitchImpermeable) {
        if (SwitchTwoLayer && Lw_ > SoilDepth2->Drc-0.01)
            return;
        if (!SwitchTwoLayer && Lw_ > SoilDepth1->Drc-0.01)
            return;
    }

    if(SwitchTwoLayer && Lw_ < SoilDep1) {

        double SoilDep2 = SoilDepth2->Drc;

        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc); // m/timestep
        Percolation = sqrt(Percolation * (Ksateff->Drc*_dt/3600000.0));
        //flux across boundary is ks+ke/2, average after Swatre

        double moistw = Lw_ * dth; //available moisture
        Percolation = std::min(moistw, Percolation);
        moistw -= Percolation; // decrease moistw with percolation, can be 0

        double Lwo = Lw_;
        Lw_ = moistw/dth;
        // factor*moisture above wettingfront, redistribution water

        // add Percolation to total unsat part
        double m0 = (Lwo-Lw_)*dth;        // moisture in part that was evacuated
        double m1 = (SoilDep1 - Lwo) * theta; // moisture below org wettingfront
        double m2 = (SoilDep2 - SoilDep1) * ThetaI2->Drc; // moisture below wettingfront in Layer 2

        m0 += m0/(m0+m1+m2)*Percolation;
        m1 += m1/(m0+m1+m2)*Percolation;
        m2 += m2/(m0+m1+m2)*Percolation;

        Thetaeff->Drc = (m0+m1)/(SoilDep1-Lw_);
        ThetaI2->Drc = m2/(SoilDep2-SoilDep1);

if (Thetaeff->Drc > Poreeff->Drc) qDebug() << Thetaeff->Drc;

        Lw->Drc = Lw_;
    }

    if(SwitchTwoLayer && Lw_ >= SoilDep1) {
        // and wetting front is in second layer

        double pore2 = ThetaS2->Drc;
        double thetar2 = 0.025 * pore2;
        double theta2 = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;
        double dth2 = (pore2-thetar2) * factor;
        double FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore2;
        //double TT2 = (pore2-FC2)*(SoilDep2-Lw_)/Ksat2->Drc;

        // get the potential flux in m
        theta_E = (theta2-thetar2)/(pore2-thetar2);
        Percolation = Ksat2->Drc*_dt/3600000.0 * pow(theta_E, bca2->Drc); // m/timestep
        Percolation = sqrt(Percolation * (Ksat2->Drc*_dt/3600000.0));
        //flux across boundary is (ks+ke)/2, average after Swatre

        //get the available moisture*factor above the wet front in m
        // factor is to avoid that the moisture goes back to thetar suddenly
        double moist1 = SoilDep1 * dth; // moisture above wettingfront in layer 1
        double moist2 = (Lw_-SoilDep1) * dth2; // moisture above wettingfront in layer 2

        Percolation = std::min(Percolation, moist1+moist2);
        //cannot have more flux than total water

        double perc2 = std::min(Percolation, moist2); // part taken from layer 2
        double perc1 = std::max(0.0, Percolation-perc2); // part taken from layer 1, can be 0

        // drain moist2 first, can be 0
        moist2 -= perc2;
        // Lw_ is now at soildep1 or lower
        double Lwo = Lw_;
        Lw_ = moist2/dth2 + SoilDep1;
        double m0 = (Lwo-Lw_)*dth2;        // moisture in part that was evacuated
        double m2 = (SoilDep2 - Lwo) * theta2; //old moisture below in L2        
        theta2 = (m2+m0+perc2)/(SoilDep2-Lw_);  //total new moisture below new Lw
        theta2 = std::min(theta2, pore2);
        // if redistribution moves Lw into layer 1
        if (perc1 > 0) {
            // moist2 is 0 and Lw_ is soildep1
            moist1 -= perc1;
            Lwo = Lw_;
            Lw_ = moist1/dth; // new Lw_
            m0 = (Lwo-Lw_)*dth;        // moisture in part that was veacuated
            theta = (m0+perc1)/(SoilDep1-Lw_);
            theta = std::min(pore,theta);
        }

        Thetaeff->Drc = theta;
        ThetaI2->Drc = theta2;
    }
    Lw->Drc = Lw_;

}

//---------------------------------------------------------------------------
// percolation from the bottom of the soil profile
void TWorld::cell_Percolation(int r, int c)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    Percolation = 0;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    if(SwitchTwoLayer) {
        pore = ThetaS2->Drc;
        thetar = 0.025 * pore;
        theta_E = 0;
        theta = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;

      //  if(theta > thetar) {
            theta_E = theta/pore; //(theta-thetar)/(pore-thetar);
            //Percolation = Ksat2->Drc*_dt/3600000 * pow(theta_E, bca2->Drc);
            // percolation in m per timestep

            if (Lw->Drc > SoilDepth1->Drc)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // if Wet Fr still in first layer percolation only make 2nd drier

            double moisture = std::max(0.0, dL*(theta-thetar));

            if (Lw_ < SoilDep2-0.01) {
                // decrease thetaeff because of percolation
                Percolation = Ksat2->Drc*_dt/3600000 * pow(theta_E, bca2->Drc);
                Percolation = std::min(Percolation, moisture);
                moisture -= Percolation;
                theta = moisture/dL;//+thetar;
            } else {

                // wetting front = soildepth2, dL = 0, moisture = 0
                // assume theta goes back to 0.7 pore and decrease the wetting fornt
                Percolation = Ksat2->Drc*_dt/3600000;
                double FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore;
                theta = FC2;//(pore - thetar);
                Lw_ -= std::max(0.0, Percolation/(pore - theta));
            }
            ThetaI2->Drc = theta;
        //}
    } else {
        // one layer
        double pore = Poreeff->Drc;
        thetar = 0.025 * pore;
        double theta = Thetaeff->Drc;

        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc);

        if (Percolation > 0) {
            if (Lw_ < SoilDep1) {
                // wetting front has not reached bottom, make soil drier
                // decrease thetaeff because of percolation
                dL = SoilDep1 - Lw_;
                double moisture = dL*(theta-thetar);
                moisture -= Percolation;
                theta = moisture/dL+thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume tehta goes back to 0.7 pore and decrease the wetting fornt
                double FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;
                theta = FC;//(pore - thetar);
                //theta = 0.7*(pore - thetar);
                Lw_ -= std::max(0.0, Percolation/(pore - theta));
            }
            Thetaeff->Drc = theta;
        }
    }

    Lw->Drc = Lw_;
    Perc->Drc = Percolation;
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
        //else
        //  CS = 0;

        LeafDrain->Drc = std::max(0.0, Cv*(Rainc_ - (CS - CStor->Drc)));
        // diff between new and old strage is subtracted from rainfall
        // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!
        // note: cover already implicit in LAI and Smax, part falling on LAI is cover*rainfall

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

