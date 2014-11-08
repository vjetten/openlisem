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
 \file lisChannelflood.cpp
 \brief Channel flood using a various solutions of St Venant equations: \n
        and more stable 1st and 2nd order st Venant following the fullSWOF2D code (univ Orleans)\n
        called before ChannelFlow(), takes old channel overflow height and spreads it out, puts new channelWH \n
        back into channel before kin wave of channel is done in ChannelFlow()
        
functions: \n
- void TWorld::ChannelOverflow(void) Mixing of flood and overflow in channel cells, source of overflow
- void TWorld::ChannelFlood(void) Calculate channelflood height maps (hmx, QFlood, UVFlood) and FloodDomain
*/

#include "lisemqt.h"
#include "model.h"
#include "global.h"

#define MINHMX 0.0
//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
void TWorld::ChannelOverflow(void)
{
  FOR_ROW_COL_MV_CH
  {
    if (ChannelDepth->Drc > 0 && ChannelMaxQ->Drc == 0 && LDD->Drc != 5)
      {
        double fc = ChannelWidthUpDX->Drc/_dx;

        // fraction reaching the channel
        double levee = ChannelLevee->Drc;
        double chdepth = ChannelDepth->Drc + levee;


        if (hmx->Drc < MINHMX+levee && ChannelWH->Drc < chdepth+MINHMX)
          continue;


        double whlevel = (ChannelWH->Drc - chdepth)*fc + _max(0.0, hmx->Drc-levee)*(1-fc);
        // new water level = weighed values of channel surplus level + hmx, levee is counted as barrier
        // can be negative if channelwh is below channel depth and low hmx level

        //if average water level is positive, water redistributes instantaneously and
        // hmx and channel wh are equal
        // normally this goes by a wave velocity probably because water level is flat V = sqrt(gh)
        if (whlevel > 0)
          {
            hmx->Drc = _min(hmx->Drc, levee);
            // cutoff hmx at levee but can be smaller
            hmx->Drc += whlevel;
            ChannelWH->Drc = whlevel + chdepth;

          }
        else
          {
            //qDebug() << nr << r << c << chdepth << whlevel << ChannelWH->Drc << ChannelWH->Drc+hmx->Drc/fc << hmx->Drc <<  fc;
            ChannelWH->Drc = ChannelWH->Drc+_min(0.0,hmx->Drc-levee-MINHMX)/fc;
            hmx->Drc = _min(MINHMX+levee, hmx->Drc);
          }
      }
    // ChannelWaterVol->Drc = ChannelWH->Drc * (ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2.0 * ChannelDX->Drc;
    // recalc chan volume for looping
    // only necessary if flooding is done after kin wave channel, gives instability here
  } // channel cells
}

//---------------------------------------------------------------------------
// correct mass balance
double TWorld::correctMassBalance(double sum1, CTMap *M, double minV)
{
  double sum2 = 0;
  double n = 0;
  FOR_ROW_COL_MV
  {
    if(M->Drc > 0)
      {
        sum2 += M->Drc;
        if(M->Drc > minV)
          n += 1;
      }
  }
  // toal and cells active for M

  double dh = (n > 0 ? (sum1 - sum2)/n : 0);
  FOR_ROW_COL_MV
  {
    if(M->Drc > minV)
      {
        M->Drc += dh;
        M->Drc = _max(M->Drc , 0.0);
      }
  }
  return dh;
}
//---------------------------------------------------------------------------
void TWorld::FloodSpuriousValues()
{
  tm->fill(0);
  FOR_ROW_COL_MV
  {
    if (hmx->Drc > F_extremeHeight)
      {
        tm->Drc = hmx->getWindowAverage(r, c);
      }
  }

  FOR_ROW_COL_MV
  {
    if ((hmx->Drc > F_extremeHeight*2) || (hmx->Drc > F_extremeHeight && hmx->Drc > tm->Drc + F_extremeDiff))
      {
        double htmp = hmx->Drc;
        hmx->Drc = _min( tm->Drc, _min(hmx->Drc, Hmx->Drc));
        qDebug() << hmx->Drc << Hmx->Drc << tm->Drc << htmp << r << c ;
      }
  }

}
//---------------------------------------------------------------------------
void TWorld::FloodBoundary()
{
//  tm->copy(hmx);
  FOR_ROW_COL_MV
  {
    if (FloodEdge->Drc > 0 && hmx->Drc > 0)
      {
   //     qDebug() << Qflood->Drc*_dt/(DX->Drc*ChannelAdj->Drc);
        hmx->Drc = _max(0, hmx->Drc - Qflood->Drc*_dt/(DX->Drc*ChannelAdj->Drc));
        floodBoundaryTot += Qflood->Drc*_dt;
      }
  }
   // tm->calcMap(hmx, SUB);
  //  tm->report("diffh");

  //volumme weg is sum diff * cell size
}
//---------------------------------------------------------------------------
void TWorld::FloodMaxandTiming()
{
  // floodwater volume and max flood map
  FloodWaterVol->fill(0);
  FOR_ROW_COL_MV
  {
    FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;

    floodHmxMax->Drc = _max(floodHmxMax->Drc, hmx->Drc);
    if (hmx->Drc > minReportFloodHeight)
      timeflood->Drc += _dt/60;
    // for output
  }

  floodVolTotMax = 0;
  floodAreaMax = 0;
  double area = _dx*_dx;
  FOR_ROW_COL_MV
  {
    if (floodHmxMax->Drc > minReportFloodHeight)
      {
        floodVolTotMax += floodHmxMax->Drc*area;
        floodAreaMax += area;
      }
  }

  FOR_ROW_COL_MV
  {
    if (hmx->Drc > minReportFloodHeight && FloodTimeStart->Drc == 0)
      {
        //            FloodTimeStart->Drc = (time - RainpeakTime)/60;
        FloodTimeStart->Drc = (time - RainstartTime)/60;
        // time since first pixel received rainfall
      }
  }
}
//---------------------------------------------------------------------------
// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
  if (!SwitchIncludeChannel)
    return;
  if (!SwitchChannelFlood)
    return;

  ChannelOverflow();
  // mix overflow water and flood water in channel cells

  double sumh_t = hmx->mapTotal();
  double dtflood = 0;

  startFlood = false;
  FOR_ROW_COL_MV
  {
    if (hmx->Drc > 0)
      {
        startFlood = true;
        break;
      }
  }

  if (SwitchFloodSWOForder2)
    {
      dtflood = fullSWOF2Do2(hmx, Uflood, Vflood, DEM);

      FOR_ROW_COL_MV
      {
        UVflood->Drc = 0.5*(fabs(Uflood->Drc)+fabs(Vflood->Drc));

        Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
      }
    }
  else
    if (SwitchFloodSWOForder1)
      {
        dtflood = fullSWOF2Do1(hmx, Uflood, Vflood, DEM);
        FOR_ROW_COL_MV
        {
          UVflood->Drc = 0.5*(Uflood->Drc+Vflood->Drc);
          Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
        }
      }
    else
      if (SwitchFloodExplicit)
        {
          dtflood = floodExplicit();
        }

  FloodSpuriousValues();


  correctMassBalance(sumh_t, hmx, 1e-12);
  // correct mass balance


  ChannelOverflow();
  // mix overflow water and flood water in channel cells

  //new flood domain
  double cells = 0;
  sumh_t = 0;
  FOR_ROW_COL_MV
  {
    if (hmx->Drc > 0)
      {
        FloodDomain->Drc = 1;
        cells += 1.0;
        sumh_t += hmx->Drc;
      }
    else
      FloodDomain->Drc = 0;
  }


  Hmx->copy(hmx);

  FloodMaxandTiming();

  FloodBoundary();

  //double avgh = (cells > 0 ? (sumh_t)/cells : 0);
  double area = cells*_dx*_dx;
  //    debug(QString("Flooding (dt %1 sec, n %2): avg h%3 m, area %4 m2").arg(dtflood,6,'f',3).arg(iter_n,4).arg(dh ,6,'e',1).arg(area,8,'f',1));
  debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(cells));
  // some error reporting


}
