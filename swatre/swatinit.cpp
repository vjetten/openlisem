/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/
//#include "swatremisc.h"
#include "csf.h"
#include "error.h"
#include "model.h"

//--------------------------------------------------------------------------------
SOIL_MODEL *TWorld::InitSwatre(
		TMMap *profileMap,
		QString initHeadMaps,    /* init head name */
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
		return(NULL);
	}
//	s->nrCells = nrCells;
	s->minDt = minDt;
   s->pixel = new PIXEL_INFO[nrCells];
   
//	s->precision = precision;
//	s->geometric = geom;
//	s->swatreBottomClosed = bottomClosed;
//	s->calibrationfactor = calibration;

	for (i = 0; i < nrCells; i++)
	{
		s->pixel[i].h = new REAL8[nodeDataIncr];
		for (n = 0; n < nrNodes; n++)
			s->pixel[i].h[n] = -1e10;
		//SetMemMV(&s->pixel[i].h,nodeDataIncr,CR_REAL8);
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
         s->pixel[r*_nrCols+c].currDt = minDt;
         s->pixel[r*_nrCols+c].h[n] = inith->Data[r][c];
		}
	}
	return(s);
}
//--------------------------------------------------------------------------------
/* soil model instance to be freed */
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
