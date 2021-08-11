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

//---------------------------------------------------------------------------
void TWorld::cell_Redistribution(int r, int c)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    Percolation = 0;

    if(!SwitchTwoLayer || Lw_ < SoilDep1) {
        // ======= one layer ============
        // or wetting front is not yet in second layer
        double pore = Poreeff->Drc;
        thetar = 0.025 * pore;
        double theta = Thetaeff->Drc;

        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc); // m/timestep
        Percolation = 0.5*(Percolation + (Ksateff->Drc*_dt/3600000.0));
        //flux across boundary is ks+ke/2, average after Swatre

        if (Percolation > 0) {
            // only redristribution if there is a wetting front
            if (Lw_ > 0) {
                double moistw = Lw_ * (pore-thetar); // moisture above wettingfront
                double Perc1 = std::min(moistw, Percolation);
                moistw -= Perc1;
                Lw_ = moistw/(pore-thetar);
                // decrease Lw with Percolation and add to unsat part
                dL = std::max(0.0, SoilDep1 - Lw_);
                double moisture = dL*(theta-thetar);
                moisture += Perc1;
                theta = std::min(pore, moisture/dL+thetar);

                Lw->Drc = Lw_;
                Thetaeff->Drc = theta;
            }
        } //Percolation > 0
    } // one layer

    if(SwitchTwoLayer && Lw_ > SoilDep1) {
        // ======= two layers ============
        // and wetting front is in second layer

        pore = ThetaS2->Drc;
        thetar = 0.025 * pore;
        theta = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;

        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = Ksat2->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc); // m/timestep
        Percolation = 0.5*(Percolation + (Ksat2->Drc*_dt/3600000.0));
        //flux across boundary is ks+ke/2, average after Swatre

        if (Percolation > 0) {
            // only redristribution if there is a wetting front
            double moistw = (Lw_-SoilDep1) * (pore-thetar); // moisture above wettingfront
            double Perc1 = std::min(moistw, Percolation);
            moistw -= Perc1;
            Lw_ = moistw/(pore-thetar) + SoilDep1;
            // decrease Lw with Percolation and add to unsat part
            dL = std::max(0.0, SoilDep2 - Lw_); // depth below wetting fornt
            double moisture = dL*(theta-thetar);
            moisture += Perc1;
            theta = std::min(pore, moisture/dL+thetar);

            Lw->Drc = Lw_;
            ThetaI2->Drc = theta;
        } //Percolation > 0
    } // two layer
}

//---------------------------------------------------------------------------
void TWorld::cell_Percolation(int r, int c)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    Percolation = 0;
    if(SwitchTwoLayer) {
        pore = ThetaS2->Drc;
        thetar = 0.025 * pore;
        theta_E = 0;
        theta = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;

        if(theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = Ksat2->Drc * pow(theta_E, bca2->Drc);
            // percolation in m

            if (Lw->Drc > SoilDepth1->Drc)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // if Wet Fr still in first layer percolation only make 2nd drier

            double moisture = dL*(theta-thetar);

            if (moisture > Percolation) {
                // decrease thetaeff because of percolation
                moisture -= Percolation;
                theta = moisture/dL+thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume theta goes back to 0.7 pore and decrease the wetting fornt
                theta = 0.7*(pore - thetar);
                Lw_ -= std::max(0.0, Percolation/(pore - theta));
            }
            ThetaI2->Drc = theta;
        }
    } else {
        // ======= one layer ============
        double pore = Poreeff->Drc;
        thetar = 0.025 * pore;
        double theta = Thetaeff->Drc;

        if(theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc); // m/timestep
        }

        if (Percolation > 0) {
            if (Lw_ > 0) {
                double moistw = Lw_ * (pore-thetar); // moisture above wettingfront
                double Perc1 = std::min(moistw, Percolation);
                moistw -= Perc1;
                Lw_ = moistw/(pore-thetar);
                // decrease Lw with Percolation and add to unsat part
                dL = std::max(0.0, SoilDep1 - Lw_);
                double moisture = dL*(theta-thetar);
                moisture += Perc1 - Percolation;
                moisture = std::max(0.0,moisture);
                theta = moisture/dL+thetar;
            }



            dL = std::max(0.0, SoilDep1 - Lw_);
            if (dL > 0) {
                double moistw = Lw_ * (pore-thetar);
                double moisture = dL*(theta-thetar);
                if (moistw > Percolation)
                    moistw -= Percolation;
                    Lw_ = moistw/(pore-thetar);

                }


            }


            double moisture = dL*(theta-thetar);
            if (moisture > Percolation) {
                // wetting front has not reached bottom, make soil drier
                // decrease thetaeff because of percolation
                moisture -= Percolation;
                theta = moisture/dL+thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume tehta goes back to half pore and decrease the wetting fornt
                theta = 0.7*(pore - thetar);
                Lw_ -= std::max(0.0, Percolation/(pore - theta));
            }
            Thetaeff->Drc = theta;
        } //Percolation > 0
    }

    if (Percolation > 0) {
        double moisture = dL*(theta-thetar);
        if (moisture > Percolation) {
            // wetting front has not reached bottom, make soil drier
            // decrease thetaeff because of percolation
            moisture -= Percolation;
            theta = moisture/dL+thetar;
        } else {
            // wetting front = soildepth1, dL = 0, moisture = 0
            // assume tehta goes back to half pore and decrease the wetting fornt
            theta = 0.7*(pore - thetar);
            Lw_ -= std::max(0.0, Percolation/(pore - theta));
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

