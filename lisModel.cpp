/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisModel.cpp:
- creating and destroying TWorld
- mass balance calculations
- link to simple interface output
- main model time loop calling process functions
- run and stop thread functions
---------------------------------------------------------------------------*/

#include <QtGui>
#include "ifacebasic.h"
#include "model.h"
#include "global.h"


//---------------------------------------------------------------------------
TWorld::TWorld(QObject *parent) :
QThread(parent)
{
	moveToThread(this);
}
//---------------------------------------------------------------------------
TWorld::~TWorld()
{
}
//---------------------------------------------------------------------------
void TWorld::DEBUGs(QString SSS)
{
	// work around to show a QString when debugging
	char ptr[128];
	strncpy(ptr, SSS.toAscii().constData(), 126);
	int i = 1;
}
//---------------------------------------------------------------------------
void TWorld::MassBalance(void)
{
	// WATER in m3

	double areafac = 1000/CatchmentArea; // conversion from m3 total to mm area average
//	double rainfall = Rain->MapTotal() * CatchmentArea;

	RainAvgmm = Rain->MapTotal()/nrCells * 1000; // avg mm over whole area
	RainTotmm = RainTotmm + RainAvgmm;

	tm->copy(Rain);
   tm->calc(DX, MUL);
   tm->calcV(_dx, MUL);
   double rainfall = tm->MapTotal(); // in m3
   RainTot += rainfall; // in m3

	double oldpeak = Rainpeak;
	Rainpeak = max(Rainpeak, rainfall);
	if (oldpeak < Rainpeak)
		RainpeakTime = time;

	IntercTot = Interc->MapTotal();
	IntercTotmm = IntercTot*areafac;

	InfilTot += InfilVol->MapTotal() + InfilVolKinWave->MapTotal(); //m3
	InfilKWTot += InfilVolKinWave->MapTotal();
	InfilTotmm = InfilTot*areafac;

	tm->calc2(WHstore, FlowWidth, MUL); //m2
	tm->calcV(_dx, MUL); //m3
	SurfStorTot = tm->MapTotal();

	WaterVolTot = WaterVolall->MapTotal();//m3
	WaterVolTotmm = WaterVolTot*areafac;

	Qtot += Qoutflow->MapTotal();
	// sum all outflow m3 for all timesteps, is already mult by dt
	Qtotmm = Qtot*areafac;

	// NOTE peak time is detected in lisOverlandlfow.cpp

	FOR_ROW_COL_MV
	{
		TotalWatervol->Drc = WaterVolall->Drc;
	}
	// SEDIMENT in kg
	if (SwitchErosion)
	{
		DetTotSplash += DETSplash->MapTotal();
		DetTotFlow += DETFlow->MapTotal();
		DepTot += DEP->MapTotal();
		DetTot += DETSplash->MapTotal() + DETFlow->MapTotal();
		SoilLossTot += Qsoutflow->MapTotal();
		SedVolTot = SedVol->MapTotal();

		FOR_ROW_COL_MV
		{
			TotalSedvol->Drc = SedVol->Drc;
			TotalConc->Drc = (TotalWatervol->Drc > 0? TotalSedvol->Drc/TotalWatervol->Drc : 0);
		}
	}

	// needed for output conc
	// Channel
	if (SwitchIncludeChannel)
	{
		WaterVolTot += ChannelWaterVol->MapTotal(); //m3
		Qtot += ChannelQoutflow->MapTotal();

		FOR_ROW_COL_MV
		{
			TotalWatervol->Drc += ChannelWaterVol->Drc;
		}

		if (SwitchErosion)
		{
			ChannelDetTot += ChannelDetFlow->MapTotal();
			ChannelDepTot += ChannelDep->MapTotal();
			ChannelSedTot = ChannelSedVol->MapTotal();

			SoilLossTot += ChannelQsoutflow->MapTotal();

			FOR_ROW_COL_MV
			{
				TotalSedvol->Drc += ChannelSedVol->Drc;
				TotalConc->Drc = (TotalWatervol->Drc > 0? TotalSedvol->Drc/TotalWatervol->Drc : 0);
				//TODO add gully, wheeltracks etc
			}
		}
	}

	// Mass Balance
	if (RainTot > 0)
		MB = (RainTot - IntercTot - InfilTot - WaterVolTot - Qtot)/RainTot*100;

	if (SwitchErosion && DetTot > 0)
		MBs = (DetTot + ChannelDetTot - SoilLossTot - SedVolTot - ChannelSedTot + DepTot + ChannelDepTot)/DetTot*100;
}
//---------------------------------------------------------------------------
// fill output structure to talk to interface
// op is declared in ifacebasic
// report to screen, hydrographs and maps
void TWorld::Output()
{
	runstep++;

	op.dx = _dx;
	op.SwitchErosion = SwitchErosion;
	op.SwitchIncludeChannel = SwitchIncludeChannel;
	op.MB = MB;
	op.runstep = runstep;
	op.maxstep = (int) ((EndTime-BeginTime)/_dt);
	op.EndTime = EndTime/60.0;
	op.CatchmentArea = CatchmentArea;

	op.RainTotmm=RainTotmm;
	op.WaterVolTotmm=WaterVolTotmm;
	op.Qtotmm=Qtotmm;
	op.Qtot=Qtot;
	op.Qpeak=Qpeak;
	op.QpeakTime=QpeakTime/60;
	op.InfilTotmm=InfilTotmm;
	op.IntercTotmm=IntercTotmm;
	op.InfilKWTotmm=InfilKWTot; // infil part in kin wave not used
	op.RunoffFraction = (RainTotmm > 0 ? Qtotmm/RainTotmm : 0);

	op.MBs = MBs;
	op.DetTotSplash=DetTotSplash*0.001;
	op.DetTotFlow=DetTotFlow*0.001;
	op.DepTot=DepTot*0.001;
	op.SedVolTot=SedVolTot*0.001;

	op.ChannelDetTot=ChannelDetTot*0.001;
	op.ChannelDepTot=ChannelDepTot*0.001;
	op.ChannelSedTot=ChannelSedTot*0.001;

	op.SoilLossTot=SoilLossTot*0.001;

	op.t = time_ms.elapsed()*0.001/60.0;
	op.time = time/60;
	op.maxtime = op.t/runstep * op.maxstep;

	emit show();

	if (runstep > 0 && runstep % printinterval == 0)
		printstep++;
	// printstep determines map reporting frequency

	ReportTimeseriesNew();
	// report hydrographs ande swdigraophs at all points in outpoint.map

	ReportTotalsNew();
	// report totals to a text file

	ReportMaps();
	// report all maps and mapseries

}
//---------------------------------------------------------------------------
// the actual model with the main loop
void TWorld::DoModel()
{

	time_ms.start();
	// get time to calc run length
	temprunname = op.runfilename;

	try
	{
		emit debug("reading and initializing data");

		IntializeOptions();
		// set all to 0 and false
		GetRunFile();
		ParseInputData();
		// get and parse runfile

		InitMapList();
		// map structure to destroy data automatically
		GetInputData();
		IntializeData();
		GetRainfallData();
		// get all input data and create and initialize all maps and variables

		BeginTime = getvaluedouble("Begin time") * 60;
		EndTime = getvaluedouble("End time") * 60;
		_dt = getvaluedouble("Timestep");
		//time vraiables in sec

		DEBUG(QString("running"));

		runstep = 0;
		printstep = 0;
		for (time = BeginTime; time < EndTime; time += _dt)
		{
			mutex.lock();
			if(stopRequested) {emit debug("User interrupt...");}
			if(stopRequested) break;
			mutex.unlock();
			// check if user wants to quit

			GridCell();          // set channel widths, flowwidths road widths etc
			Rainfall();
			Interception();
			Infiltration();
			SoilWater();
			SurfaceStorage();
			SplashDetachment();
			FlowDetachment();
			OverlandFlow();
			ChannelFlow();
			MassBalance();

			Output();
		}

		DestroyData();  // destroy all maps automatically

		emit done("finished");
	}
	catch(...)  // if an error occurred
	{
		DestroyData();
		emit done("ERROR STOP: "+ErrorString);
	}
}
//---------------------------------------------------------------------------
void TWorld::run()
{
	QTimer::singleShot(0, this, SLOT(DoModel()));
	exec();
}
//---------------------------------------------------------------------------
void TWorld::stop()
{
	QMutexLocker locker(&mutex);
	stopRequested = true;
}
//---------------------------------------------------------------------------
