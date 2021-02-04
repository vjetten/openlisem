
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisSnowmelt.cpp
  \brief Get snowmelt adta and make a map. Snowmelt is like rainfall, melt intensities

functions: \n
- void TWorld::SnowmeltMap(void) \n
 */

#include <memory>
#include "io.h"
#include "model.h"


//---------------------------------------------------------------------------
/**
 snowmelt intensity read is that reported with the next line: example
 0 0\n
 5 2.3   ->from 0 to 5 minutes intensity is 2.3\n
 7.5 4.5 ->from 5 to 7.5 minutes intensity is 4.5\n
 etc. */
void TWorld::SnowmeltMap(void)
{
   double timeminprev = (time-_dt) / 60; //prev time in minutes
   int  place;
   double tt = 3600000.0;

   if (!SwitchSnowmelt)
      return;

   for (place = 0; place < nrSnowmeltseries; place++)
      if (timeminprev < SnowmeltSeriesM[place].time)
         break;

   if (SnowmeltSeriesM[place].isMap)
   {
      auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(
          SnowmeltSeriesM[place].name)));

      #pragma omp parallel for num_threads(userCores)
      FOR_ROW_COL_MV_L {
          if (pcr::isMV(_M->Drc)) {
              QString sr, sc;
              sr.setNum(r); sc.setNum(c);
              ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+SnowmeltSeriesM[place].name;
              throw 1;
          }
          else
          Snowmelt->Drc = _M->Drc *_dt/tt;
      }}
   }
   else
   {
#pragma omp parallel for num_threads(userCores)
      FOR_ROW_COL_MV_L {
         Snowmelt->Drc = SnowmeltSeriesM[place].intensity[(int) SnowmeltZone->Drc-1]*_dt/tt;
         // Rain in m per timestep from mm/h, rtecord nr corresponds map nID value -1
      }}
   }
#pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_L
   {
       Snowmeltc->Drc = Snowmelt->Drc * _dx/DX->Drc;
       // correction for slope dx/DX, water spreads out over larger area
       SnowmeltCum->Drc += Snowmeltc->Drc;
       // cumulative rainfall corrected for slope, used in interception
       //RainNet->Drc = Rainc->Drc;
       // net rainfall in case of interception

   }}
}
//---------------------------------------------------------------------------
