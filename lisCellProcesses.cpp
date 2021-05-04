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

        if(theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = Ksat2->Drc * pow(theta_E, bca->Drc);
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
        // one layer
        double pore = Poreeff->Drc;
        thetar = 0.025 * pore;
        double theta = Thetaeff->Drc;

        if(theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca->Drc);
        }

        if (Percolation > 0) {
            dL = std::max(0.0, SoilDep1 - Lw_);
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
        }
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



void TWorld::do_CellProcesses()
{
      RainfallMap();         // get rainfall from table or mpas
      SnowmeltMap();         // get snowmelt

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        if (SwitchRainfall) {
//            double timeminprev = (time-_dt) / 60; //prev time in minutes
//            int  rainplace;
//            double tt = 3600000.0;

//            for (rainplace = 0; rainplace < nrRainfallseries; rainplace++)
//                if (timeminprev < RainfallSeriesM[rainplace].time)
//                    break;

//            if (RainfallSeriesM[rainplace].isMap)
//            {
//                auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(RainfallSeriesM[rainplace].name)));
//                if (pcr::isMV(_M->Drc))
//                {
//                    QString sr, sc;
//                    sr.setNum(r); sc.setNum(c);
//                    ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+RainfallSeriesM[rainplace].name;
//                    throw 1;
//                }
//                else
//                    Rain->Drc = _M->Drc *_dt/tt;
//            }
//            else
//            {
//                Rain->Drc = RainfallSeriesM[rainplace].intensity[(int) RainZone->Drc-1]*_dt/tt;
//            }

//            if (!rainStarted) {
//                if(Rain->Drc > 0)
//                    rainStarted = true;
//            }
//            if (rainStarted && RainstartTime == -1)
//                RainstartTime = time;

//            Rainc->Drc = Rain->Drc * _dx/DX->Drc;
//            // correction for slope dx/DX, water spreads out over larger area
//            RainCumFlat->Drc += Rain->Drc;
//            // cumulative rainfall
//            RainCum->Drc += Rainc->Drc;
//            // cumulative rainfall corrected for slope, used in interception
//            RainNet->Drc = Rainc->Drc;
//            // net rainfall in case of interception

            if (Rainc->Drc > 0)
                cell_Interception(r,c);
        }

        //========================

        double FloodDomain_ = FloodDomain->Drc;
        if (FloodDomain_ == 0) {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
        } else {
            hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        }
        double RW = RoadWidthDX->Drc;
        if ((SwitchRoadsystem || SwitchHardsurface) && RW > 0) {
            WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
        }

        //===== INFILTRATION =========

        if (InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE) {

            cell_InfilMethods(r, c);
            if (!SwitchImpermeable)
                cell_Percolation(r, c);

        }
        if (InfilMethod == INFIL_SWATRE) {
            cell_InfilSwatre(r, c);
        }

        //===== SURFACE STORAGE =====

        double SW = SoilWidthDX->Drc;
        double WHr = WHroad->Drc;
        double WHs; //WHstore
        double WH_ = WH->Drc;

        //### surface storage on rough surfaces
        WHs = std::min(WH_, MDS->Drc*(1-exp(-1.875*(WH_/std::max(0.01,0.01*RR->Drc)))));
        // non-linear release fo water from depression storage
        // resembles curves from GIS surface tests, unpublished

        double FW = std::min(ChannelAdj->Drc, SW + RW);
        // calculate flowwidth by fpa*surface + road, excludes channel already

        WHrunoff->Drc = ((WH_ - WHs)*SW + WHr*RW)/FW;
        FlowWidth->Drc = FW;

        WaterVolall->Drc = DX->Drc*(WH_*SW + WHr*RW);
        WHstore->Drc = WHs;

        //========================

        if (SwitchErosion) {
            double wh = FloodDomain_ == 0 ? WH->Drc : hmx->Drc;
            cell_SplashDetachment(r,c,wh);
        }

    }}
report(*Lw, "lw");

}

//---------------------------------------------------------------------------
void TWorld::cell_ETa(int r, int c)
{
    //NOTE: all in meter
    double tot = 0;
    double ETp = 5.0/(12.0*3600.0) * 0.001 *_dt;  //5m/day in m per dt
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SwitchTwoLayer ? SoilDepth2->Drc-SoilDepth1->Drc : 0.0;
    double thetai1 = ThetaI1->Drc;
    double thetai2 = SwitchTwoLayer ? ThetaI2->Drc : 0.0;
    double thetas1 = ThetaS1->Drc;
    double thetas2 = SwitchTwoLayer ? ThetaS2->Drc : 0.0;
    double cover = Cover->Drc;

    // soil ET from top 15 cm where there is no cover and no ponding
    double ETa_soil = 0;
    ETa_soil = hmx->Drc > 0 || WH->Drc > 0 ? 0.0 : thetai1/ThetaS1->Drc * ETp;
    ETa_soil = thetai1 < 0.05 ? 0.0 : ETa_soil;
    ETa_soil *= 1-cover;
    if (ETa_soil > 0) {
        double moist = thetai1 * std::min(150.0,SoilDep1); // assume soil evap fomr top 15 cm
        ETa_soil = std::min(ETa_soil, moist);
        moist = moist - ETa_soil;
        double moist1 = thetai1 * std::max(0.0, SoilDep1-150.0);
        ThetaI1->Drc = (moist1+moist)/SoilDep1;
        tot = tot + ETa_soil;
    }

    // interception ET where there is cover
    double cstor = CStor->Drc;
    double ETa_int = cover * ETp;
    ETa_int = std::min(cstor, ETa_int);
    cstor = cstor-ETa_int;
    Interc->Drc =  cover * cstor * SoilWidthDX->Drc * DX->Drc;
    CStor->Drc = cstor;
    tot = tot + ETa_int;

    //transpiration form the entire soil where cover
    thetai1 = ThetaI1->Drc;
    double ETa_Tr = ETa_int == 0 ? cover * ETp : 0.0;
    double moist = thetai1*SoilDep1 + thetai2 * SoilDep2;
    double the = (thetai1+thetai2)/(thetas1+thetas2);
    double ETfactor = 1.0-(1.0/(1.0+pow(the/((thetas1+thetas2)*0.5), 5.0)));
    ETa_Tr = std::min(moist*0.95, ETa_Tr);
    ETa_Tr = ETa_Tr * ETfactor;

    moist = moist - ETa_Tr;
    double moist1 = moist*SoilDep1/(SoilDep1+SoilDep2);
    double moist2 = moist*SoilDep2/(SoilDep1+SoilDep2);
    ThetaI1->Drc = moist1/SoilDep1;
    if (SwitchTwoLayer)
        ThetaI2->Drc = moist2/SoilDep2;
    tot = tot + ETa_Tr;

    // surface waters
    double wh = WH->Drc;
    double whr = WHrunoff->Drc;
    double hm = hmx->Drc;
    double ETa_pond = hm > 0 || wh > 0 ? ETp : 0.0;
    if (hm > 0) {
        ETa_pond = std::min(ETa_pond, hmx->Drc);
        hm = hm-ETa_pond;
        hmx->Drc = hm;
    } else {
        ETa_pond = std::min(ETa_pond, wh);
        wh = wh-ETa_pond;
        if (wh < WHstore->Drc) {
            WHstore->Drc = wh;
            whr = 0;
        } else {
            whr = wh - WHstore->Drc;
        }
        WaterVolall->Drc = whr*CHAdjDX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        WHroad->Drc = whr;
        WH->Drc = wh;
        WHrunoff->Drc = whr;
    }
    tot = tot + ETa_pond;

    ETa->Drc = tot;
    ETaCum->Drc += tot;
}
