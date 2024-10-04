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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
\file swatinit.cpp
\brief SWATRE: initialize soil profile with inithead maps data and clean up after run

functions:
- SOIL_MODEL * TWorld::InitSwatre(cTMap *profileMap); \n
- void TWorld::CloseSwatre(SOIL_MODEL *s); \n
*/

#include "lerror.h"
#include "model.h"

//--------------------------------------------------------------------------------
SOIL_MODEL *TWorld::InitSwatre(cTMap *profileMap)
{
   SOIL_MODEL *s = (SOIL_MODEL *)malloc(sizeof(SOIL_MODEL));
   /* TODO check if this needs freeing when error */

   int  i, n, nrNodes  = ((zone == nullptr) ? -1 : zone->nrNodes);
   int nodeDataIncr = nrNodes+1;
   long nrCells = _nrCols*_nrRows;

   s->minDt = swatreDT;
   s->pixel = new PIXEL_INFO[nrCells];
//Fill(*tmb,0);
   // set initial values
   for (i = 0; i < nrCells; i++)
   {
      s->pixel[i].profile = nullptr;
      s->pixel[i].h = new REAL8[nodeDataIncr];
      s->pixel[i].theta = new REAL8[nodeDataIncr];
      s->pixel[i].k = new REAL8[nodeDataIncr];
      for (n = 0; n < nrNodes; n++) {
         s->pixel[i].h[n] = -1e10;
         s->pixel[i].theta[n] = 0.01;
         s->pixel[i].k[n] = 0;
      }

      s->pixel[i].nrNodes = nrNodes;  //set to 1 for output of a pixel
      s->pixel[i].dumpHid = 0;  //set to 1 for output of a pixel
      s->pixel[i].tiledrain = 0;
      s->pixel[i].tilenode = -1;
   //   s->pixel[i].repellency = 0;// no repellency as init
      // set tiledrain to 0, and tiledepth to -1 (above surface)
      s->pixel[i].currDt = swatreDT;
   }

   // give each pixel a profile 
   FOR_ROW_COL_MV
   {
      s->pixel[r*_nrCols+c].profile = profileList[swatreProfileNr.indexOf((int)profileMap->Drc)];
            //ProfileNr(profileMap->Drc);
      // profileNr throws an error if profile nr not found

      // if (SwitchWaterRepellency)
      //    s->pixel[r*_nrCols+c].repellency = (int)RepellencyCell->Drc;

      if(SwitchDumpH || SwitchDumpTheta || SwitchDumpK) {
          s->pixel[r*_nrCols+c].dumpHid = SwatreOutput->Drc;
      }

   }

   // fill the inithead structure of each pixel and set tiledrain depth if any
   for (n = 0; n < nrNodes; n++)
   {
      QString fname = QString("%1.%2").arg(initheadName).arg(n+1, 3, 10, QLatin1Char('0'));
      // make inithead.001 to .00n name

      cTMap *inith = ReadMap(LDD,fname);
      // get inithead information

      FOR_ROW_COL_MV
      {
         s->pixel[r*_nrCols+c].h[n] = inith->data[r][c];//*psiCalibration;

         // find depth of tilenode
         if (SwitchIncludeTile) {
             if (!pcr::isMV(TileDepth->Drc) && TileDepth->Drc > 0)
             {
                 // NOTE depth is in m while node info is in cm, so *100
                 // endComp is the depth at the bottom of the compartment, so the tile is <= endcomp
                 if (s->pixel[r*_nrCols+c].profile->zone->endComp[n] > TileDepth->Drc*100)
                     s->pixel[r*_nrCols+c].tilenode = n-1;
             }
         }
      }
   }
   return(s);
}
//--------------------------------------------------------------------------------
/// soil model instance to be freed
void TWorld::CloseSwatre(SOIL_MODEL *s)
{
   if (s == nullptr)
      return;

   swatreProfileDef.clear();
   swatreProfileNr.clear();

   for (int i = 0; i < _nrCols*_nrRows; i++)
      delete[] s->pixel[i].h;

   free(s->pixel);
   free(s);
   s = nullptr;
}
//--------------------------------------------------------------------------------
