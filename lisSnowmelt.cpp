
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
  \file lisSnowmelt.cpp
  \brief Get snowmelt adta and make a map. Snowmelt is like rainfall, melt intensities

functions: \n
- void TWorld::GetSnowmeltData(void) \n
- void TWorld::SnowmeltMap(void) \n
 */

#include "model.h"

//---------------------------------------------------------------------------
/// read snowmelt file and put the dta in SnowmeltSeries
/*
void TWorld::GetSnowmeltData(void)
{
	QFile fff(snowmeltFileName);
	QFileInfo fi(snowmeltFileName);
	QString S;
	bool ok;
	int j = 0;

	if (!fi.exists())
	{
		ErrorString = "Snowmelt file not found: " + snowmeltFileName;
		throw 1;
	}

	nrSnowmeltstations = 0;
	nrSnowmeltseries = 0;

	fff.open(QIODevice::ReadOnly | QIODevice::Text);
	S = fff.readLine();
	// read the header
	while (S.isEmpty())
		S = fff.readLine();
	QStringList SL = S.split(QRegExp("\\s+")); //<== white spave character as split
	nrSnowmeltstations = SL[SL.size()-2].toInt(&ok, 10);


	if (!ok)
	{
		ErrorString = "Cannot read nr rainfall stations in header rainfall file";
		throw 1;
	}

	for(int r=0; r < nrSnowmeltstations+1; r++)
		S = fff.readLine();
	// read column headers

	while (!fff.atEnd())
	{
		S = fff.readLine();
		if (!S.trimmed().isEmpty())
			nrSnowmeltseries++;
	}
	// count rainfall records, skip empty lines
	if (nrSnowmeltseries <= 1)
	{
		ErrorString = "Snowmelt records <= 1!";
		throw 1;
	}

	nrSnowmeltseries++;
	SnowmeltSeries = new double*[nrSnowmeltseries];
	for(int r=0; r < nrSnowmeltseries; r++)
		SnowmeltSeries[r] = new double[nrSnowmeltstations+1];
	// make structure to contain Snowmelt
	//SnowmeltSeries is matrix with rows is data and 1st column is time, other columns are stations

	fff.close();
	// close file and start again

	fff.open(QIODevice::ReadOnly | QIODevice::Text);
	S = fff.readLine();
	for (int i = 0; i < nrSnowmeltstations; i++)
		S = fff.readLine();
	// read header data
	while (!fff.atEnd())
	{
		S = fff.readLine();
		if (S.trimmed().isEmpty()) continue;

		QStringList SL = S.split(QRegExp("\\s+"));
		if (SL.size()-1 < nrSnowmeltstations)
		{
			QString ss;
			ErrorString = "Snowmelt: Nr stations specified in header = " + ss.setNum(nrSnowmeltstations);
			ErrorString += ", nr columns available = "+ ss.setNum(SL.size()-1);
			throw 1;
		}
		SnowmeltSeries[j][0] = SL[0].toDouble();
		// time in min
		for (int i = 1; i < nrSnowmeltstations+1; i++)
			SnowmeltSeries[j][i] = SL[i].toDouble();
      // snowmelt intensities
		j++;
	}
	SnowmeltSeries[nrSnowmeltseries-1][0] = 1e20;
	for (int i = 1; i < nrSnowmeltstations+1; i++)
		SnowmeltSeries[nrSnowmeltseries-1][i] = 0;
	// end series with 0 value and extreme time
	fff.close();
}
*/
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
      CTMap *_M = new CTMap();
      _M->PathName = SnowmeltSeriesM[place].name;
      bool res = _M->LoadFromFile();
      if (!res)
      {
         ErrorString = "Cannot find Snowmelt map " +_M->PathName;
         throw 1;
      }

      for (int r = 0; r < _nrRows; r++)
         for (int c = 0; c < _nrCols; c++)
            if (!IS_MV_REAL8(&LDD->Drc) && IS_MV_REAL8(&_M->Drc))
            {
               QString sr, sc;
               sr.setNum(r); sc.setNum(c);
               ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+SnowmeltSeriesM[place].name;
               throw 1;
            }
      FOR_ROW_COL_MV
      {
         Snowmelt->Drc = _M->Drc *_dt/tt;
         Snowmeltc->Drc = Snowmelt->Drc * _dx/DX->Drc;
         SnowmeltCum->Drc += Snowmeltc->Drc;
      }

      _M->KillMap();
   }
   else
   {
      FOR_ROW_COL_MV
      {
         Snowmelt->Drc = SnowmeltSeriesM[place].intensity[(int) SnowmeltZone->Drc-1]*_dt/tt;
         // Rain in m per timestep from mm/h, rtecord nr corresponds map nID value -1
         Snowmeltc->Drc = Snowmelt->Drc * _dx/DX->Drc;
         // correction for slope dx/DX, water spreads out over larger area

         //TODO: weighted average if dt larger than table dt

         SnowmeltCum->Drc += Snowmeltc->Drc;
      }
   }

}
/*
void TWorld::SnowmeltMap(void)
{
	double timemin = time / 60;  //time in minutes
	double timeminp = (time-_dt) / 60; //prev time in minutes
	int place;
   double tt = 3600000.0;

	if (!SwitchSnowmelt)
		return;

//	for (placep = 0; placep < nrSnowmeltseries; placep++)
//		if (timeminp < SnowmeltSeries[placep][0])
//			break;
	for (place = 0; place < nrSnowmeltseries; place++)
		if (timemin < SnowmeltSeries[place][0])
			break;


	FOR_ROW_COL_MV
	{
      int col = (int) SnowmeltZone->Drc;

      Snowmelt->Drc = SnowmeltSeries[place][col]/3600000 * _dt;
      // Snowmelt in m per timestep
      Snowmeltc->Drc = Snowmelt->Drc * _dx/DX->Drc;
      // DO NOT correct for slope dx/DX, snow is already on the surface
      // TODO: CORRECT OR NOT

      // TODO: weighted average if dt larger than table dt

      SnowmeltCum->Drc += Snowmeltc->Drc;
	}
}
*/
//---------------------------------------------------------------------------
