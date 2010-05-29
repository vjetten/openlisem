/*---------------------------------------------------------------------------
project: openLISEM
name: lisRainintc.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisRainintc.cpp:
- read rainfall file
- transfer rainfall to a map
- calc interception
---------------------------------------------------------------------------*/

#include "model.h"

//---------------------------------------------------------------------------
// read rainfall files of different types
void TWorld::GetRainfallData(void)
{
	QFile fff(rainFileName);
	QFileInfo fi(rainFileName);
	QString S;
	bool ok;
	int j = 0;

	if (!fi.exists())
	{
		ErrorString = "Rainfall file not found: " + rainFileName;
		throw 1;
	}

	nrstations = 0;
	nrrainfallseries = 0;

	fff.open(QIODevice::ReadOnly | QIODevice::Text);
	S = fff.readLine();
	// read the header
	if (S.contains("RUU CSF TIMESERIE", Qt::CaseInsensitive)) // file is old lisem file
	{
		QStringList SL = S.split(QRegExp("\\s+"));
		nrstations = SL[SL.size()-2].toInt(&ok, 10);
	}

	else // file is PCRaster timeseries
	{
		S = fff.readLine();
		nrstations = S.toInt(&ok, 10);
	}
	if (!ok)
	{
		ErrorString = "Cannot read nr rainfall stations in header rainfall file";
		throw 1;
	}

	for(int r=0; r < nrstations+1; r++)
		S = fff.readLine();
	// read column headers

	while (!fff.atEnd())
	{
		S = fff.readLine();
		if (!S.trimmed().isEmpty())
			nrrainfallseries++;
	}
	// count rainfall records, skip empty lines
	if (nrrainfallseries <= 1)
	{
		ErrorString = "rainfall records <= 1!";
		throw 1;
	}

	nrrainfallseries++;
	RainfallSeries = new double*[nrrainfallseries];
	for(int r=0; r < nrrainfallseries; r++)
		RainfallSeries[r] = new double[nrstations+1];
	// make structure to contain rainfall
	//RainfallSeries is matrix with rows is data and 1st column is time, other columns are stations

	fff.close();
	// close file and start again

	fff.open(QIODevice::ReadOnly | QIODevice::Text);
	S = fff.readLine();
	for (int i; i < nrstations; i++)
		S = fff.readLine();
	// read header data
	while (!fff.atEnd())
	{
		S = fff.readLine();
		if (S.trimmed().isEmpty()) continue;

		QStringList SL = S.split(QRegExp("\\s+"));
		if (SL.size()-1 < nrstations)
		{
			QString ss;
			ErrorString = "Rainfall: Nr stations specified in header = " + ss.setNum(nrstations);
			ErrorString += ", nr columns available = "+ ss.setNum(SL.size()-1);
			throw 1;
		}
		RainfallSeries[j][0] = SL[0].toDouble();
		// time in min
		for (int i = 1; i < nrstations+1; i++)
			RainfallSeries[j][i] = SL[i].toDouble();
		// rainfall intensities
		j++;
	}
	RainfallSeries[nrrainfallseries-1][0] = 1e20;
	for (int i = 1; i < nrstations+1; i++)
		RainfallSeries[nrrainfallseries-1][i] = 0;
	// end series with 0 value and extreme time
	fff.close();
}
//---------------------------------------------------------------------------
// rainfall intensity read is that reported with the next line: example
// 0 0
// 5 2.3   ->from 0 to 5 minutes intensity is 2.3
// 7.5 4.5 ->from 5 to 7.5 minutes intensity is 4.5
// etc.
void TWorld::Rainfall(void)
{
	double timemin = time / 60;  //time in minutes
	double timeminp = (time-_dt) / 60; //prev time in minutes
	int placep, place;

	for (placep = 0; placep < nrrainfallseries; placep++)
		if (timeminp < RainfallSeries[placep][0])
			break;
	for (place = 0; place < nrrainfallseries; place++)
		if (timemin < RainfallSeries[place][0])
			break;

	FOR_ROW_COL_MV
	{
			int col = (int) RainZone->Drc;

			Rain->Drc = RainfallSeries[place][col]/3600000 * _dt;
			// Rain in m per timestep
			Rainc->Drc = Rain->Drc * _dx/DX->Drc;
			// correction for slope dx/DX, water spreads out over larger area

			//TODO: weighted average if dt larger than table dt

			RainCum->Drc += Rainc->Drc;
			// cumulative rainfall corrected for slope, used in interception
	}
}
//---------------------------------------------------------------------------
//- interception seen as rigid storage SMax filling up and overflowing
//- overflow flux is identical to rainfall flux in intensity
//- SMax is the storage of the plants inside the gridcell, not the average storage of the gridcell
// so if a single tree inside a cell has an SMax of 2mm even if it covers 10%, the Smax of that cell is 2
// - therefore the same goes for LAI: the LAI of the plants inside the gridcell
// this is also easier to observe. The LAI from a satellite image is the average LAI of a cell, must be divided by Cover
void TWorld::Interception(void)
{
	// all variables are in m

	FOR_ROW_COL_MV
	{
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

		if (Smax > 0)
			CS = Smax*(1-exp(-0.0653*LAIv*RainCum->Drc/Smax));
		else
			CS = 0;
		// 0.0653 is canopy openess, based on Aston (1979), based on Merriam (1960/1973), De Jong & Jetten 2003
		// is not the same as canopy cover. it also deals with how easy rainfall drips through the canopy
		//possible to use equation from Ahston but for very open Eucalypt

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
	}
}
//---------------------------------------------------------------------------



