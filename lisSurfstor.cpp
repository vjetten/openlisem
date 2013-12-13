
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisSurfstor.cpp
  \brief calculate surface storage and flow width

functions: \n
- void TWorld::GridCell(void) \n
- void TWorld::addRainfallWH(void) \n
- void TWorld::SurfaceStorage(void)\n
 */


#include "model.h"
#define tiny 1e-8



//---------------------------------------------------------------------------
void TWorld::GridCell(void)
{
    FOR_ROW_COL_MV
    {
        if (BufferID->Drc > 0)
            RoadWidthDX->Drc = 0;
        //VJ 100609 cannot have a road with a buffer, to complicated

        if (RoadWidthDX->Drc + HouseWidthDX->Drc > _dx)
            HouseWidthDX->Drc = _dx-RoadWidthDX->Drc;
        // road takes priority

        if (SwitchIncludeChannel)
        {
            if (RoadWidthDX->Drc + ChannelWidthUpDX->Drc > _dx)
                RoadWidthDX->Drc = ChannelAdj->Drc;
            if (HouseWidthDX->Drc + ChannelWidthUpDX->Drc > _dx)
                HouseWidthDX->Drc = ChannelAdj->Drc;
        }
        // channel takes priority

        /** wheeltracks are not implemented yet */
        //      WheelWidthDX->Drc = 0;
        //      if (SwitchWheelPresent)
        //         WheelWidthDX->Drc = max(0,_dx-RoadWidthDX->Drc-ChannelWidthUpDX->Drc)/_dx * WheelWidth->Drc;
        // adjust wheelwidth in cells with other surfaces
        /** TODO is wheelwidth needed or just an extra map? */

        SoilWidthDX->Drc = max(0, _dx - ChannelWidthUpDX->Drc
                               - RoadWidthDX->Drc);
        //- HouseWidthDX->Drc);

        //      SoilWidthDX->Drc = max(0, _dx - ChannelWidthUpDX->Drc
        //                             - GullyWidthDX->Drc
        //                             - RoadWidthDX->Drc
        //                             - WheelWidthDX->Drc);
    }
}

//---------------------------------------------------------------------------
/// Adds new rainfall afterinterception to runoff water nheight or flood waterheight
void TWorld::addRainfallWH(void)
{
    FOR_ROW_COL_MV
    {
        if (FloodDomain->Drc > 0)
            hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        else
        {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
            // add net to water rainfall on soil surface (in m)

            if (SwitchBuffers && !SwitchSedtrap)
                if(BufferID->Drc > 0 && BufferVol->Drc > 0)
                {
                    WH->Drc = 0;
                    BufferVol->Drc  += (Rainc->Drc + Snowmeltc->Drc) * DX->Drc * _dx;
                }
            // buffers and not full yet (buffervol > 0) then add rainflal to buffers and set WH to zero
            // not for sed traps, behave normally

            if (GrassFraction->Drc > 0)
                WHGrass->Drc += RainNet->Drc + Snowmeltc->Drc;
            // net rainfall on grass strips, infil is calculated separately for grassstrips

            if (RoadWidthDX->Drc > 0)
                WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
            // assume no interception and infiltration on roads, gross rainfall
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::SurfaceStorage(void)
{
    FOR_ROW_COL_MV
    {
        double RRm = 0.01*RR->Drc; // assume RR in cm convert to m
        double wh = WH->Drc, whflow = 0;
        double SDS;
        double mds = MDS->Drc;  // mds is in meters
        double WaterVolrunoff;

        //### surface storage
        SDS = 0.1*mds;
        // arbitrary minimum depression storage is 10% of max depr storage, in m

        if (mds > 0)
            whflow = (wh-SDS) * (1-exp(-1000*wh*(wh-SDS)/(mds-SDS)));
        //could be: whflow = (wh-SDS) * (1-exp(-wh/mds));
        // non-linear release fo water from depression storage
        // resemles curves from GIS surface tests, unpublished
        else
            whflow = wh;

        // whflow = max(0, wh-SDS);
        // subtract surface storage and calc water available for runoff, in m
        // assumed on soilsurface because there is the roughness

        WHstore->Drc = wh - whflow;
        // average water stored on flowwidth and not available for flow, in m

        //houses no surf storage
        if (SwitchHouses)
        {
            WHstore->Drc *= (1-HouseCover->Drc);
            whflow = wh - WHstore->Drc;
        }

        WaterVolrunoff = DX->Drc*( whflow*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
        // runoff volume available for flow, surface + road
        WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
        // all water in the cell incl storage

        //### the trick: use ponded area for flowwidth
        if (RRm == 0)
            fpa->Drc = 1;
        else
            fpa->Drc = 1-exp(-1.875*(wh/RRm));
        // fraction ponded area of a gridcell

        if (SwitchChannelFlood)
            if (FloodDomain->Drc > 0)
                fpa->Drc = 1;

        FlowWidth->Drc = max(0.01*_dx, fpa->Drc*SoilWidthDX->Drc + RoadWidthDX->Drc);
        // calculate flowwidth by fpa*surface + road, excludes channel already

        if (GrassFraction->Drc > 0)
            FlowWidth->Drc = GrassWidthDX->Drc + (1-GrassFraction->Drc)*FlowWidth->Drc;
        // assume grassstrip spreads water over entire width

        //Houses
        if (SwitchHouses)
            FlowWidth->Drc = (1-0.5*HouseCover->Drc)*FlowWidth->Drc;
        // assume house severely restricts flow width, 0.5 is arbitrary
        // cannot be zero flowwidth in 100% house pixel because watwer would not go anywhere

        if (FlowWidth->Drc > 0)
            WHrunoff->Drc = WaterVolrunoff/(DX->Drc*FlowWidth->Drc);
        else
            WHrunoff->Drc = 0;
        // average WHrunoff from soil surface + roads, because kin wave can only do one discharge
        // this now takes care of ponded area, so water height is adjusted

        //WHrunoffCum->Drc += WHrunoff->Drc * 1000;
        // cumulative runoff for output maps, in mm
        // make no sense? runoff isa water layer independent of the timestep,
        // you cannot acumulate it like this.

    }
}
//---------------------------------------------------------------------------

// TRIAL TO PUT ALL PROCESSES IN ONE LOOP
void TWorld::Allprocs(void)
{
    double timeminprev = (time-_dt) / 60; //prev time in minutes
    int  rainplace;
    double tt = 3600000.0;

    for (rainplace = 0; rainplace < nrRainfallseries; rainplace++)
        if (timeminprev < RainfallSeriesM[rainplace].time)
            break;
    // find rainfall time interval

    InfilVol->fill(0);
    fact->fill(0);
    fpot->fill(0);


    FOR_ROW_COL_MV
    {
        if (BufferID->Drc > 0)
            RoadWidthDX->Drc = 0;
        //VJ 100609 cannot have a road with a buffer, to complicated

        if (RoadWidthDX->Drc + HouseWidthDX->Drc > _dx)
            HouseWidthDX->Drc = _dx-RoadWidthDX->Drc;
        // road takes priority

        if (SwitchIncludeChannel)
        {
            if (RoadWidthDX->Drc + ChannelWidthUpDX->Drc > _dx)
                RoadWidthDX->Drc = ChannelAdj->Drc;
            if (HouseWidthDX->Drc + ChannelWidthUpDX->Drc > _dx)
                HouseWidthDX->Drc = ChannelAdj->Drc;
        }
        // channel takes priority

        SoilWidthDX->Drc = max(0, _dx - ChannelWidthUpDX->Drc
                               - RoadWidthDX->Drc);


        //rainfall

        if (SwitchRainfall)
        {

            Rain->Drc = RainfallSeriesM[rainplace].intensity[(int) RainZone->Drc-1]*_dt/tt;
            // Rain in m per timestep from mm/h, rtecord nr corresponds map nID value -1
            Rainc->Drc = Rain->Drc * _dx/DX->Drc;
            // correction for slope dx/DX, water spreads out over larger area

            //TODO: weighted average if dt larger than table dt

            RainCum->Drc += Rainc->Drc;
            // cumulative rainfall corrected for slope, used in interception
            RainNet->Drc = Rainc->Drc;
        }

        //interception
        double CS = CStor->Drc;
        //actual canopy storage in m
        double Smax = CanopyStorage->Drc;
        //max canopy storage in m
        double LAIv;
        if (SwitchInterceptionLAI)
            LAIv = LAI->Drc;
        else
            LAIv = (log(1-Cover->Drc)/-0.4)/max(0.9,Cover->Drc);
        //Smax is based on LAI and LAI is the average of a gridcell, already including the cover
        // a low cover means a low LAI means little interception
        // avoid division by 0

        if (SwitchBuffers && !SwitchSedtrap)
            if(BufferID->Drc > 0)
                Smax = 0;
        // no interception with buffers, but sedtrap can have interception

        if (SwitchHardsurface)
            Smax *= (1-HardSurface->Drc);
        //VJ 110111 no interception on hard surfaces

        if (PlantHeight->Drc < WH->Drc)
        {
            Smax = 0;
            CS = 0;
        }
        //VJ no interception when water level is heigher than plants


        if (Smax > 0)
        {
            //double k = 1-exp(-CanopyOpeness*LAIv);
            //VJ !!!!!!2013 05 18 BUG ! was k=exp(-coLAI) MUST BE 1-exp
            double k = exp(-CanopyOpeness*LAIv);
            //VJ 131010 NOT !!!!! a dense canopy has a low openess factor, so little direct throughfall and high CS

            CS = Smax*(1-exp(-k*RainCum->Drc/Smax));
            //      CS = Smax*(1-exp(-0.0653*LAIv*RainCum->Drc/Smax));
            //VJ 110209 direct use of openess, astons value quite open for eucalypt.
            //A good guess is using the cover LAI relation
            //and interpreting cover as openess: k = exp(-0.45*LAI) BUT this is 0.065!!!
        }
        else
            CS = 0;
        // 0.0653 is canopy openess, based on Aston (1979), based on Merriam (1960/1973), De Jong & Jetten 2003
        // is not the same as canopy cover. it also deals with how easy rainfall drips through the canopy
        // possible to use equation from Ahston but for very open Eucalypt

        CS = max(0, CS * (1-StemflowFraction));
        //VJ 110206 decrease storage with stemflow fraction!

        LeafDrain->Drc = max(0, Cover->Drc*(Rainc->Drc - (CS - CStor->Drc)));
        // diff between new and old strage is subtracted from rainfall
        // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!
        // note: cover already implicit in LAI and Smax, part falling on LAI is cover*rainfall

        CStor->Drc = CS;
        // put new storage back in map
        Interc->Drc =  Cover->Drc * CS * SoilWidthDX->Drc * DX->Drc; //*
        // only on soil surface, not channels or roads, in m3
        // cover already implicit in CS, Smax

        RainNet->Drc = LeafDrain->Drc + (1-Cover->Drc)*Rainc->Drc;
        // net rainfall is direct rainfall + drainage
        // rainfall that falls on the soil, used in infiltration

        if (SwitchHouses)
        {
            if (HouseCover->Drc > 0)
            {
                double HS, DS;
                //actual roof storage in m
                double Hmax = RoofStore->Drc;
                //max roof storage in m
                // HouseCover->Drc = qMin(HouseCover->Drc, 0.95);

                double Dmax =0;
                if (SwitchRaindrum)
                    if (HouseCover->Drc > 0)
                        Dmax = DrumStore->Drc/(CellArea->Drc*HouseCover->Drc);
                //max drum storage in m

                double housedrain = 0;
                //overflow in m

                if (Hmax > 0)
                {
                    double k = 1.0;
                    // speed of filling of roof storage, very quickly
                    HS = Hmax*(1-exp(-k*RainCum->Drc/Hmax));
                    //roof storage in m
                }
                else
                    HS = 0;

                if (Dmax > 0)
                {
                    double k = 0.05;
                    // speed of filling of raindrum near house, slower
                    DS = Dmax*(1-exp(-k*RainCum->Drc/Dmax));
                    //drum storage in m

                }
                else
                    DS = 0;

                housedrain = max(0, (RainNet->Drc - (HS - HStor->Drc) - (DS - DStor->Drc)));
                // diff between new and old strage is subtracted from rainfall
                // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!

                HStor->Drc = HS;
                DStor->Drc = DS;
                // put new storage back in maps in m and m

                IntercHouse->Drc = HouseCover->Drc * (HS+DS) * SoilWidthDX->Drc * DX->Drc;//_dx*_dx
                // total interception in m3

                RainNet->Drc = HouseCover->Drc*housedrain + (1-HouseCover->Drc)*RainNet->Drc;
                // net rainfall is direct rainfall + drainage
                // rainfall that falls on the soil, used in infiltration
            }
            else
                IntercHouse->Drc = 0;
        }

        if (FloodDomain->Drc > 0)
            hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        else
        {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
            // add net to water rainfall on soil surface (in m)

            if (SwitchBuffers && !SwitchSedtrap)
                if(BufferID->Drc > 0 && BufferVol->Drc > 0)
                {
                    WH->Drc = 0;
                    BufferVol->Drc  += (Rainc->Drc + Snowmeltc->Drc) * DX->Drc * _dx;
                }
            // buffers and not full yet (buffervol > 0) then add rainflal to buffers and set WH to zero
            // not for sed traps, behave normally

            if (GrassFraction->Drc > 0)
                WHGrass->Drc += RainNet->Drc + Snowmeltc->Drc;
            // net rainfall on grass strips, infil is calculated separately for grassstrips

            if (RoadWidthDX->Drc > 0)
                WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
            // assume no interception and infiltration on roads, gross rainfall
        }

//infiltration

            InfilVol->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

            // potential water volume on surface before infil

            // calculate effective ksat for various situations
            if (InfilMethod != INFIL_SWATRE  && InfilMethod != INFIL_NONE)
            {
                Ksateff->Drc = Ksat1->Drc;
                if (SwitchInfilCrust || SwitchInfilCompact)
                    Ksateff->Drc = Ksat1->Drc*(1-CompactFraction->Drc-CrustFraction->Drc) +
                            KsatCrust->Drc*CrustFraction->Drc + KsatCompact->Drc*CompactFraction->Drc;
                // when not switched on fractions and ksat are 0
                // avg ksat of "normal" surface with crusting and compaction fraction
                // adjust effective infil for crusting and compaction
                //VJ 110106 adapted this calculation

                if (SwitchHardsurface && HardSurface->Drc > 0)
                    Ksateff->Drc = (1-HardSurface->Drc)*Ksateff->Drc;// =  0;
                //VJ 110111 no infiltration on hard surfaces

                //houses
                if (SwitchHouses)
                    Ksateff->Drc = Ksateff->Drc * (1-HouseCover->Drc);
                //VJ decrease ksat for celss with houses

                if (GrassFraction->Drc > 0)
                    Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;

                Ksateff->Drc *= ksatCalibration;
                // apply runfile/iface calibration factor

                if (SwitchBuffers && !SwitchSedtrap)
                    if(BufferID->Drc > 0)
                        Ksateff->Drc = 0;
                //VJ 1000608 no infil in buffers, , but sedtrap can have infil

            double fact1 = 0;
            double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
            double fwh = WH->Drc; // in m, WH is old WH + net rainfall
            double Psi = Psi1->Drc;
            double space = max(ThetaS1->Drc-ThetaI1->Drc, tiny);

            if (SoilDepth1->Drc <= tiny)
            {
                fpot->Drc = 0;
                fact->Drc = 0;
                continue;
            }
            // outcrops etc: no infil in this cell

            if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
            {
                Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
                // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
                Psi = Psi2->Drc;
                space = max(ThetaS2->Drc-ThetaI2->Drc, tiny);
            }
            // two layers

            // calculate potential infiltration fpot in m
            switch (InfilMethod)
            {
            case INFIL_KSAT : fpot->Drc = Ks; break;
            case INFIL_GREENAMPT :
            case INFIL_GREENAMPT2 :
                fpot->Drc = Ks*(1.0+(Psi+fwh)/(L1->Drc+L2->Drc)); break;
            case INFIL_SMITH :
            case INFIL_SMITH2 :
                double B = (fwh + Psi)*space;
                double Cdexp = exp(Fcum->Drc/B);
                fpot->Drc = Ks*Cdexp/(Cdexp-1);
                break;
            }

            fact1 = min(fpot->Drc, fwh);
            // actual infil in m, cannot have more infil than water on the surface

            fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc, &FFull->Drc);
            // adjust fact and increase L1 and L2, for twolayer, impermeable etc

                if (RoadWidthDX->Drc == _dx)
                {
                    fact->Drc = 0;
                    fpot->Drc = 0;
                }
                // to make sure WH is not decreasd when all is road

                WH->Drc -= fact->Drc;
                if (WH->Drc < 0) // in case of rounding of errors
                {
                    fact->Drc += WH->Drc;
                    WH->Drc = 0;
                }
                // subtract fact->Drc from WH, cannot be more than WH

                Fcum->Drc += fact->Drc;
                // cumulative infil in m

            FSurplus->Drc = min(0, fact->Drc - fpot->Drc);
            // negative surplus of infiltration in m for kinematic wave in m

            if (FFull->Drc == 1)
                FSurplus->Drc = 0;
            //VJ 101216 if soil full and impermeable: no surplus and no extra infil in kin wave

            InfilVol->Drc -= DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
            // infil volume is WH before - water after
        }


            double RRm = 0.01*RR->Drc; // assume RR in cm convert to m
            double wh = WH->Drc, whflow = 0;
            double SDS;
            double mds = MDS->Drc;  // mds is in meters
            double WaterVolrunoff;

            //### surface storage
            SDS = 0.1*mds;
            // arbitrary minimum depression storage is 10% of max depr storage, in m

            if (mds > 0)
                whflow = (wh-SDS) * (1-exp(-1000*wh*(wh-SDS)/(mds-SDS)));
            //could be: whflow = (wh-SDS) * (1-exp(-wh/mds));
            // non-linear release fo water from depression storage
            // resemles curves from GIS surface tests, unpublished
            else
                whflow = wh;

            // whflow = max(0, wh-SDS);
            // subtract surface storage and calc water available for runoff, in m
            // assumed on soilsurface because there is the roughness

            WHstore->Drc = wh - whflow;
            // average water stored on flowwidth and not available for flow, in m

            //houses no surf storage
            if (SwitchHouses)
            {
                WHstore->Drc *= (1-HouseCover->Drc);
                whflow = wh - WHstore->Drc;
            }

            WaterVolrunoff = DX->Drc*( whflow*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
            // runoff volume available for flow, surface + road
            WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
            // all water in the cell incl storage

            //### the trick: use ponded area for flowwidth
            if (RRm == 0)
                fpa->Drc = 1;
            else
                fpa->Drc = 1-exp(-1.875*(wh/RRm));
            // fraction ponded area of a gridcell

            if (SwitchChannelFlood)
                if (FloodDomain->Drc > 0)
                    fpa->Drc = 1;

            FlowWidth->Drc = max(0.01*_dx, fpa->Drc*SoilWidthDX->Drc + RoadWidthDX->Drc);
            // calculate flowwidth by fpa*surface + road, excludes channel already

            if (GrassFraction->Drc > 0)
                FlowWidth->Drc = GrassWidthDX->Drc + (1-GrassFraction->Drc)*FlowWidth->Drc;
            // assume grassstrip spreads water over entire width

            //Houses
            if (SwitchHouses)
                FlowWidth->Drc = (1-0.5*HouseCover->Drc)*FlowWidth->Drc;
            // assume house severely restricts flow width, 0.5 is arbitrary
            // cannot be zero flowwidth in 100% house pixel because watwer would not go anywhere

            if (FlowWidth->Drc > 0)
                WHrunoff->Drc = WaterVolrunoff/(DX->Drc*FlowWidth->Drc);
            else
                WHrunoff->Drc = 0;
            // average WHrunoff from soil surface + roads, because kin wave can only do one discharge
            // this now takes care of ponded area, so water height is adjusted

            //WHrunoffCum->Drc += WHrunoff->Drc * 1000;
            // cumulative runoff for output maps, in mm
            // make no sense? runoff isa water layer independent of the timestep,
            // you cannot acumulate it like this.

            double Perim;
            const double beta = 0.6;
            const double _23 = 2.0/3.0;
            double beta1 = 1/beta;
            //double kinvisc = 1.1e-6; // 15 degrees celcius water
            double NN = N->Drc;

            if (SwitchChannelFlood)
                NN = N->Drc * qExp(mixing_coefficient*hmx->Drc);
            // slow down water in flood zone
            //NOTE CALCULATE REYNOLDS AND SEE IF TURBULENCE MAKES COEFFICIENT HIGHER
            //LAMINAR MEANS NON-MIXING SO NN = N
            //    tma->Drc = hmx->Drc * UVflood->Drc/1.1e-6;

            // avg WH from soil surface and roads, over width FlowWidth
            Perim = 2*WHrunoff->Drc+FlowWidth->Drc;

            if (Perim > 0)
                R->Drc = WHrunoff->Drc*FlowWidth->Drc/Perim;
            else
                R->Drc = 0;

            Alpha->Drc = pow(NN/sqrt(Grad->Drc) * pow(Perim, _23),beta);

            if (Alpha->Drc > 0)
                Q->Drc = pow((FlowWidth->Drc*WHrunoff->Drc)/Alpha->Drc, beta1);
            else
                Q->Drc = 0;

            V->Drc = pow(R->Drc, _23)*sqrt(Grad->Drc)/NN;


    }//ROWCOL

}
