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
 \file lisChannelflow.cpp
 \brief Channel hydrology and sediment detachment and movement processes.

functions: \n
- void TWorld::CalcVelDischChannel(void) calculate Velocity, alpha and Q in the channel \n
- void TWorld::ChannelFlow(void) calculate channelflow, channelheight, do kinematic wave \n
*/

#include "model.h"

//---------------------------------------------------------------------------
// V, alpha and Q in the channel
void TWorld::CalcVelDischChannel(void)
{
   /*
    dw      FW      dw
   \  |            |  /
    \ |         wh | /
     \|____________|/
  */
   FOR_ROW_COL_MV_CH
   {
      double Perim, Radius, Area;
      const double beta = 0.6;
      const double _23 = 2.0/3.0;
      double beta1 = 1/beta;
      double wh = ChannelWH->Drc;
      double FW = ChannelWidth->Drc;
      double grad = sqrt(ChannelGrad->Drc);
      double dw = 0.5*(ChannelWidthUpDX->Drc - FW); // extra width when non-rectamgular

      if (dw > 0)
      {
         //			Perim = FW + 2*sqrt(wh*wh + dw*dw);
         Perim = FW + 2*wh/cos(atan(ChannelSide->Drc));
         Area = FW*wh + wh*dw;// NOT*2 BECAUSE TRIANGLE !!;
      }
      else
      {
         Perim = FW + 2*wh;
         Area = FW*wh;
      }
      //Perim = ChannelWidth->Drc + 2*wh/cos(atan(ChannelSide->Drc));
      // cos atanb more expensive than sqrt ?
      //Area = ChannelWidth->Drc*wh + wh*(ChannelWidthUpDX->Drc - ChannelWidth->Drc);

      ChannelPerimeter->Drc = Perim;
      //VJ 110109 needed for channel infil

      if (Perim > 0)
         Radius = Area/Perim;
      else
         Radius = 0;

      ChannelAlpha->Drc = pow(ChannelN->Drc/grad * powl(Perim, _23),beta);

      if (ChannelAlpha->Drc > 0)
         ChannelQ->Drc = pow(Area/ChannelAlpha->Drc, beta1);
      else
         ChannelQ->Drc = 0;

      ChannelV->Drc = pow(Radius, _23)*grad/ChannelN->Drc;
   }
//   else // WHY THIS??? RELIC !!!
//   {
//      ChannelAlpha->Drc = 0;
//      ChannelQ->Drc = 0;
//      ChannelV->Drc = 0;
//      ChannelPerimeter->Drc = 0;
//   }
}
//---------------------------------------------------------------------------
//! calc channelflow, channelheight, kin wave
void TWorld::ChannelFlow(void)
{
   if (!SwitchIncludeChannel)
      return;

   // calculate new channel WH , WidthUp and Volume
   FOR_ROW_COL_MV_CH
   {
      /*---- Water ----*/

      ChannelQsn->Drc = 0;
      Channelq->Drc = 0;
      ChannelWH->Drc = 0;

      ChannelWaterVol->Drc += RunoffVolinToChannel->Drc;
      // add inflow to channel
      ChannelWaterVol->Drc += Rainc->Drc*ChannelWidthUpDX->Drc*DX->Drc;
      // add rainfall in m3, no interception, rainfall so do not use ChannelDX

      if (SwitchBuffers && ChannelBufferVol->Drc > 0)
      {
         ChannelBufferVol->Drc -= ChannelWaterVol->Drc;
         ChannelWaterVol->Drc = 0;
         // add inflow from slopes and rainfall to buffer
      }

      if (ChannelSide->Drc == 0)// && ChannelWidth->Drc > 0)// rectangular channel
      {
         ChannelWidthUpDX->Drc = ChannelWidth->Drc;
         ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
      }
      else  // non-rectangular
      {
         /*
   ABC fornula
    dw      w       dw
   \  |            |  /
    \ |          h |a/  <= tan(a) is channelside = tan angle of side wall
     \|____________|/
   area = vol/dx = h*w + h*dw
   dw = h*tan(a)
   vol = w*h + dw*h = w*h + tan(a)*h*h
   tan(a)h^2 + wh - vol/dx = 0
     aa (h2)   +   bb(h) +  cc
*/
         double aa = ChannelSide->Drc;  //=tan(a)
         double bb = ChannelWidth->Drc; //=w
         double cc = -ChannelWaterVol->Drc/ChannelDX->Drc; //=vol/DX

         ChannelWH->Drc = (-bb + sqrt(bb*bb - 4*aa*cc))/(2*aa);
         if (ChannelWH->Drc < 0)
         {
            ErrorString = QString("channel water height is negative at row %1, col %2").arg(r).arg(c);
            throw 1;
         }
         ChannelWidthUpDX->Drc = ChannelWidth->Drc + 2*ChannelSide->Drc*ChannelWH->Drc;
      }

      if (ChannelWidthUpDX->Drc > _dx)
      {
         ErrorString = QString("channel width > dx at row %1, col %2").arg(r).arg(c);
         throw 1;
      }

      ChannelWidthUpDX->Drc = min(0.9*_dx, ChannelWidthUpDX->Drc);
      // new channel width with new WH, goniometric, side is top angle tan, 1 is 45 degr
      // cannot be more than 0.9*_dx

      if (RoadWidthDX->Drc > 0)
         ChannelWidthUpDX->Drc = min(0.9*_dx-RoadWidthDX->Drc, ChannelWidthUpDX->Drc);
      // channel cannot be wider than _dx-road
      /* TODO zit al in gridcell, nodig hier? */

      if (SwitchChannelInfil)
      {
         Channelq->Drc =  -(ChannelKsat->Drc *  ChannelPerimeter->Drc/3600000.0);
         //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
      }
      // NOTE: for buffers channelksat = 0

   }

   CalcVelDischChannel();
   // alpha, V and Q from Manning

   /*---- Sediment ----*/

   if (SwitchErosion)
   {
      ChannelFlowDetachment();

      FOR_ROW_COL_MV_CH
      {
         ChannelQs->Drc = ChannelQ->Drc * ChannelConc->Drc;
      }
   }

   ChannelQn->setMV();

   if (useSorted)
   {
      KinematicSorted(lddlistch, lddlistchnr, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, ChannelDX,
                      ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed);
   }
   else
   {
      FOR_ROW_COL_MV_CH
      {
         if (LDDChannel->Drc == 5)
         {
            Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, ChannelDX,
                      ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed, SubsMaps);
         }
      }
   }

   ChannelQn->cover(LDD, 0); // avoid missing values around channel for adding to Qn for output
   ChannelQs->cover(LDD, 0);

   FOR_ROW_COL_MV_CH
   {
      double ChannelArea = ChannelAlpha->Drc*pow(ChannelQn->Drc, 0.6);
      // in buffers ChannelQn = 0;

      ChannelWH->Drc = ChannelArea/((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2);
      // water height is not used except for output i.e. watervolume is cycled

      double diff = Channelq->Drc*_dt + ChannelWaterVol->Drc - (ChannelArea * ChannelDX->Drc) - ChannelQn->Drc*_dt;
      //difference between fluxes and store in and out of channel cell

      if (SwitchBuffers && ChannelBufferVol->Drc > 0)
      {
         //qDebug()<< ChannelBufferVol->Drc << Channelq->Drc*_dt << ChannelWaterVol->Drc << (ChannelArea * ChannelDX->Drc) << ChannelQn->Drc*_dt<< diff;
      }
      else
         if (SwitchChannelInfil)
            InfilVolKinWave->Drc += diff;
      //VJ 110111 add channel infil to infil for mass balance

      ChannelWaterVol->Drc = ChannelArea * ChannelDX->Drc;
      // total water vol after kin wave in m3, going to the next timestep
      // in a buffer ChannelArea = 0 so channelvolume is also 0

      if (SwitchErosion)
      {
         ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc, ChannelDep->Drc);
         // correct for very high concentrations, max 850 g/l
      }
   }
}
//---------------------------------------------------------------------------
void TWorld::ChannelFlood(void)
{
   /*
   if (!SwitchIncludeChannel)
      return;

   double DEMmin = 1e20;
   int flag = 0;

   FOR_ROW_COL_MV
   {
      DEMmin = DEM->Drc < DEMmin ? DEM->Drc : DEMmin;
      tm->Drc = 0;
      tma->Drc = 0;
      tmb->Drc = 0;
   }

   FOR_ROW_COL_MV_CH
   {
      if (ChannelWH->Drc - ChannelHeight->Drc > 0)
      {
         tma->Drc = ChannelWH->Drc - ChannelHeight->Drc;
         tm->Drc = 1;
      }
      // water above channel
   }
   tma->report("whbef");

   //   if (flag == 0)
   //    return;
   //no overflow
   for( int j = 0; j < 10; j++)
   {
      FOR_ROW_COL_MV
      {
         double dHdL[5];
         double H[5], h[5], z[5];
         int dc[5] = {0, 0, +1, 0, -1};
         int dr[5] = {0, -1, 0, 1, 0};
         double dQ = 0;
         double flow = 0.1;

         for (int i = 0; i < 5; i++)
         {
            int rr, cc;
            rr = r+dr[i];
            cc = c+dc[i];

            if (INSIDE(rr, cc) && !IS_MV_REAL8(&LDD->Data[rr][cc]))
            {
               z[i] = (DEM->Data[rr][cc]);// - DEMmin + 1);
               // gravity potential
               h[i] = tma->Data[rr][cc];
               // water potential
               H[i] = z[i] + h[i];
               // hydraulic potential
               dHdL[i] = sin(atan((H[i]-H[0])/_dx));
               // difference Darcy in EW and NS directions

               if (i > 0)
                  dQ = dQ + flow * _dx * h[i]*dHdL[i];
            }

         }

         tmb->Drc += dQ;

         h[0] = h[0] + dQ/(_dx*_dx);

         tma->Drc = h[0];
      }
   }
   FOR_ROW_COL_MV_CH
   {
      if (tma->Drc > 0)
      {
         ChannelWH->Drc = tma->Drc + ChannelHeight->Drc;
         ChannelWaterVol->Drc = ChannelWH->Drc * ChannelWidth->Drc * ChannelDX->Drc;

      }
      // water above channel
   }

   tmb->report("dq");
   tma->report("whaft");
*/
}

//# There are 4 steps:
//# 1) gravity based gradual GW flow
//# 2) Mass balance correction because there is no iteration
//# 3) fast GW flow through cracks reaching the river in one timestep
//# 4) Force GW to stay below the surface

//### (1) gravity based groundwater flow ###

//# Gravity based flow: potential differences
//# between GW surface in NS and EW directions
//# Total flow in/out cell is:
//# Darcy : Q = q*A = K sin(a)*(h*dx)
//# potential diff in a direction is diff hydr head over a distance along the slope dH/dL
//# distance along slope is: angle of flow is dH/dx  dL = dx*
//# A is wet cross section of the flow (cell width dx * water height h)
//# dQ = Ks * SUM [ sin(a)*h*dx ] of 4 neighbour cells in EW and NS direction
//# dQ is in m3/timestep and is added to the central cell
//# Converted to height by division by the cell area

//#NOTE: everything is converted to meters and meters/day

//dx = celllength();

//ldd2 = ldd(2*mask); #south, row + 1
//ldd4 = ldd(4*mask); #west, col - 1
//ldd6 = ldd(6*mask); #east, col + 1
//ldd8 = ldd(8*mask); #north, row - 1
//# set up basic directions for groundwater movement

//z = DEM - mapminimum(DEM) - soildepth/1000 + mapmaximum(soildepth/1000)*2;
//# z = bedrock surface in m above a fixed datum
//# max soildepth * 2 added to make z positive (not really needed)
//# gravity potential equals dem of bedrock, soildepth is in mmm, convert to m

//##!!!!!!!!!!!!!!!!!!!!!
//rise = Perc/theta_s;
//GWDepth = GWDepth + rise;
//# GW rise because of incoming percolation in mm spread out in pore space
//# this should be rise = Perc/(theta_s-theta) but that gives problems with
//# extreme fluctuations
//##!!!!!!!!!!!!!!!!!!!!!

//h = GWDepth/1000;
//# GWdepth in mm, h = matric potential in m
//SumGWbefore = maptotal(h);
//# sum GW before movement is done, needed for mass balance corrcetions

//H = h + z;
//# total hydraulic potential in m (i.e. terrain height and water level)

//dHdL2 = sin(atan((upstream(ldd2, H)-H)/dx));
//dHdL4 = sin(atan((upstream(ldd4, H)-H)/dx));
//dHdL6 = sin(atan((upstream(ldd6, H)-H)/dx));
//dHdL8 = sin(atan((upstream(ldd8, H)-H)/dx));
//# Darcy: Q = ks*dH/dL =
//# sin = dH/dL! => sine of potential differences between central cell
//# dH/dx = tan so atan(dH/dx) is angle
//# and cells in 4 directions EW ans NS (in m)

//h2 = max(h,upstream(ldd2, h));
//h4 = max(h,upstream(ldd4, h));
//h6 = max(h,upstream(ldd6, h));
//h8 = max(h,upstream(ldd8, h));
//# water height in 4 directions for cross section
//# maximum of central water height or in 4 directions

//dQ = GWksat_cal * Ksat/1000 * dx * (h2*dHdL2 + h8*dHdL8 + h4*dHdL4 + h6*dHdL6);
//# sum of all fluxes in m3/day, ksat in m/day, divide by 1000
//# H2,4,6,8 is the water height in 4 directions, dx is cell width,
//# so h*dx is cross section of flow A. Q is then Ksat*A*sin(dH/dL)
//h = h + dQ/(dx*dx);
//# add in/out flow to the cell in m

//h = max(0, min(h, soildepth/1000));
//# h must between 0 and soildepth: 0 < h < soildepth in meters
//# careful: potential source fo mass balance errors by cutting of

//### (2) GW mass balance correction instead of iteration ###

//SumGWafter = maptotal(h);
//# sum GW after the movement
//errorh = (SumGWbefore - SumGWafter)*mask;
//#total mass balance error in GW depth (before - after)
//wetcells = maptotal(scalar(h gt 0))*mask;
//# calc which cells have GW
//h = h + if(h gt 0, errorh/wetcells, 0);
//# smooth out the error over all wet cells
//h = max(0, min(h, soildepth/1000));
//# correct again: h must between 0 and soildepth: 0 < h < soildepth
//SumGWafter = maptotal(h);

//report error${1}_${2}_${3}_${4}_${5}_${6}.tss = 1000*(SumGWbefore - SumGWafter)/maptotal(wetcells);
//# average mass balance cell error in mm

//GWDepth = h*1000;
//# convert from m back to mm for comparison with the other fluxesin the model

//### (3) fast GW flow through cracks ###

//GWloss = if (GWDepth gt GWksatCrack_level, GWksatCrack_cal * GWDepth, 0);
//# subtract a loss when GWDepth is above a threshold

//GWloss = min(GWloss, GWDepth);
//# can not be more than there is
//GWDepth = GWDepth - GWloss;
//# lower the GW with GWloss each timestep

//GWloss = GWloss * theta_s;
//# GWloss is a fraction of the depth, that means the flux
//# to the river is GWloss*theta_s

//### (4) force GW below the margin zone ###

//GWDepth = if(outcrop, 0, GWDepth);
//# no groundwater where outcrops (should not occur anyway)
//##!!!!!!!!!!!!!!!!!!!!!!!!!!!
//GWSurp =  max(0, GWDepth - (soildepth-margin));
//GWSurp = if(outcrop, 0, GWSurp);
//GWDepth = GWDepth - GWsurp;
//# correct all groundwater entering the "margin" zone
//# to avoid situations where the soil is saturated to the surface
//# a saturated soil means unsatdepth = 0 which stops all fluxes!

//##!!!!!!!!!!!!!!!!!!!!!!!!!!!
