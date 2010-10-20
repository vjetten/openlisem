/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/
#include "misc.h"
#include "csf.h"
#include "swatre_p.h"
#include "swatre_g.h"
#include "swat_inp.h"
#include "error.h"
#include "model.h"

//--------------------------------------------------------------------------------
SOIL_MODEL *TWorld::InitSwatre(
		TMMap *profileMap,
		QString initHeadMaps,    /* init head name */
		double minDt,         /* minimum timestep, is also initial timestep */
		double precision,     /* precision factor to determine if if timestep can be altered */
		double calibration,
		bool geom,
		bool bottomClosed)
{
	SOIL_MODEL *s = (SOIL_MODEL *)malloc(sizeof(SOIL_MODEL));
	//TODO check if this needs freeing when error
	int  i, n, nrNodes  = NrZoneNodes();
	int nodeDataIncr = nrNodes+1;
	long nrCells = nrCols*nrRows;

	if (nrNodes == -1)
	{
		Error("SWATRE: can't call \'initswatre\' before \'swatre input\'");
		return(NULL);
	}
	s->nrCells = nrCells;
	s->precision = precision;
	s->minDt = minDt;
	s->geometric = geom;
	s->swatreBottomClosed = bottomClosed;
	s->calibrationfactor = calibration;
	s->pixel = new PIXEL_INFO[nrCells];

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

			s->pixel[r*nrCols+c].profile = ProfileNr(profileMap->Drc);
			s->pixel[r*nrCols+c].currDt = minDt;
			s->pixel[r*nrCols+c].h[n] = inith->Data[r][c];
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
    
	for (int i = 0; i < s->nrCells; i++)
		delete[] s->pixel[i].h;

	free(s->pixel);
	free(s);
	s = NULL;
}
//--------------------------------------------------------------------------------
