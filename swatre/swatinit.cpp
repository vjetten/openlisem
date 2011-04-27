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
\file swatinit.cpp
\brief SWATRE: initialize soil profile with inithead maps data and clean up after run

functions:
- SOIL_MODEL * TWorld::InitSwatre(TMMap *profileMap, QString initHeadMaps, double minDt); \n
- void TWorld::CloseSwatre(SOIL_MODEL *s); \n
*/



#include "csf.h"
#include "error.h"
#include "model.h"

//--------------------------------------------------------------------------------
SOIL_MODEL *TWorld::InitSwatre(
   TMMap *profileMap,
   QString initHeadMaps,
   TMMap *tiledepthMap,
   double minDt)
{
	SOIL_MODEL *s = (SOIL_MODEL *)malloc(sizeof(SOIL_MODEL));

   /* TODO check if this needs freeing when error */

   int  i, n, nrNodes  = ((zone == NULL) ? -1 : zone->nrNodes);
	int nodeDataIncr = nrNodes+1;
   long nrCells = _nrCols*_nrRows;

	if (nrNodes == -1)
	{
		Error("SWATRE: can't call \'initswatre\' before \'swatre input\'");
      //		return(NULL);
	}
	s->minDt = minDt;
   s->pixel = new PIXEL_INFO[nrCells];
   
   // obsolete
   //	s->nrCells = nrCells;
   //	s->precision = precision;
   //	s->geometric = geom;
   //	s->swatreBottomClosed = bottomClosed;
   //	s->calibrationfactor = calibration;

	for (i = 0; i < nrCells; i++)
	{
		s->pixel[i].h = new REAL8[nodeDataIncr];
		for (n = 0; n < nrNodes; n++)
			s->pixel[i].h[n] = -1e10;
      s->pixel[i].dumpHid = 0;
      if ( i == (_nrRows/2*_nrCols + 0.5*_nrCols))
         s->pixel[i].dumpHid = 1;
      //SetMemMV(&s->pixel[i].h,nodeDataIncr,CR_REAL8);

      s->pixel[i].tiledrain = 0;
      s->pixel[i].tilenode = -1;
      // set tiledrain to 0, and tiledepth to -1 (above surface)


   }
	for (n = 0; n < nrNodes; n++)
	{
		QString fname = QString("%1.%2").arg(initHeadMaps)
            .arg(n+1, 3, 10, QLatin1Char('0'));
		// make inithead.001 to .00n name

		TMMap *inith = ReadMap(LDD,fname);
		// get inithead information

		FOR_ROW_COL_MV
		{
			if (profileMap->Drc == -1 || ProfileNr(profileMap->Drc) == NULL)
			{
				Error(QString("SWATRE: profile nr '%1' is missing").arg(profileMap->Drc));
				return (NULL);
			}

         s->pixel[r*_nrCols+c].profile = ProfileNr(profileMap->Drc);

         // find depth of tilenode
         if (!IS_MV_REAL8(&tiledepthMap->Drc) && tiledepthMap->Drc > 0)
         {
            // NOTE depth is in m while node info is in cm, so *100
            // endComp is the depth at the bottom of the compartment, so the tile is <= endcomp
            if (s->pixel[r*_nrCols+c].profile->zone->endComp[n] >= tiledepthMap->Drc*100)
               s->pixel[r*_nrCols+c].tilenode = n;
         }
         s->pixel[r*_nrCols+c].currDt = minDt;
         s->pixel[r*_nrCols+c].h[n] = inith->Data[r][c];
		}
	}
	return(s);
}
//--------------------------------------------------------------------------------
/// soil model instance to be freed
void TWorld::CloseSwatre(SOIL_MODEL *s)
{
   if (s == NULL)
      return;

   for (int i = 0; i < _nrCols*_nrRows; i++)
		delete[] s->pixel[i].h;

	free(s->pixel);
	free(s);
	s = NULL;
}
//--------------------------------------------------------------------------------
