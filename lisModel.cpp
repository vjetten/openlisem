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
#include "lisemqt.h"
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
//	int i = 1;
}
//---------------------------------------------------------------------------
// fill output structure to talk to interface
// op is declared in ifacebasic
// report to screen, hydrographs and maps
void TWorld::OutputUI()
{
	if (runstep > 0 && runstep % printinterval == 0)
		printstep++;
	// printstep determines map reporting frequency
	runstep++;

	op.dx = _dx;
	op.MB = MB;
	op.runstep = runstep;
	op.maxstep = (int) ((EndTime-BeginTime)/_dt);
	op.EndTime = EndTime/60.0;
//	op.BeginTime = BeginTime/60.0;

	op.CatchmentArea = CatchmentArea;

	op.RainTotmm=RainTotmm;
	op.WaterVolTotmm=WaterVolTotmm-SurfStoremm;
	op.Qtotmm=Qtotmm;
	op.Qtot=QtotOutlet;
	op.Qpeak=Qpeak;
	op.QpeakTime=QpeakTime/60;
	op.RainpeakTime=RainpeakTime/60;
	op.InfilTotmm=InfilTotmm;
	op.SurfStormm=SurfStoremm;
	op.IntercTotmm=IntercTotmm;
	op.InfilKWTotmm=InfilKWTot; // infil part in kin wave not used
	op.RunoffFraction = (RainTotmm > 0 ? Qtotmm/RainTotmm : 0);

	op.MBs = MBs;
	op.DetTotSplash=DetSplashTot*0.001;
	op.DetTotFlow=DetFlowTot*0.001;
	op.DepTot=DepTot*0.001;
	op.SedTot=SedTot*0.001;

	op.ChannelDetTot=ChannelDetTot*0.001;
	op.ChannelDepTot=ChannelDepTot*0.001;
	op.ChannelSedTot=ChannelSedTot*0.001;

	op.SoilLossTot=SoilLossTotOutlet*0.001;

	op.t = time_ms.elapsed()*0.001/60.0;
	op.time = time/60;
	op.maxtime = op.t/runstep * op.maxstep;

	op.P = RainAvgmm* 3600/_dt;
	op.Q = Qoutput->DrcOutlet;
	op.Qs = Qsoutput->DrcOutlet;
	op.C = TotalConc->DrcOutlet;

	op.BufferVolTot = BufferVolTot;
	op.BufferSedTot = BufferSedTot*0.001; //ton

	emit show();
	// send the op structure with data to function Showit in the interface
	// in file LisUIModel

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
	temprunname = QString(op.LisemDir+"openlisemtmp.run");

	try
	{
		DEBUG("reading and initializing data");

		IntializeOptions();
		// set all to 0 and false
		InitMapList();
		// map structure to destroy data automatically

		DEBUG("GetRunFile()");
		GetRunFile();
		DEBUG("ParseInputData()");
		ParseInputData();
		// get and parse runfile

		DEBUG("GetInputData()");
		GetInputData();
		DEBUG("IntializeData()");
		IntializeData();
		DEBUG("GetRainfallData()");
		GetRainfallData();
		if (SwitchSnowmelt)
		{
			DEBUG("GetSnowmeltData()");
			GetSnowmeltData();
		}
		// get all input data and create and initialize all maps and variables

		BeginTime = getvaluedouble("Begin time") * 60;
		EndTime = getvaluedouble("End time") * 60;
		_dt = getvaluedouble("Timestep");
		//time vraiables in sec

		runstep = 0; // NOTE runstep is used to initialize graph!
		printstep = 1;

		DEBUG("Running...");

		for (time = BeginTime; time < EndTime; time += _dt)
		{
			mutex.lock();
			if(stopRequested) DEBUG("User interrupt...");
			if(stopRequested) break;
			mutex.unlock();

			mutex.lock();
			if (waitRequested) DEBUG("User pause...");
			if (waitRequested) condition.wait(&mutex);
			mutex.unlock();
			// check if user wants to quit

			GridCell();          // set channel widths, flowwidths road widths etc
			RainfallMap();
			SnowmeltMap();
			Interception();
			Infiltration();
			//SoilWater();
			SurfaceStorage();
		   CalcVelDisch();

			SplashDetachment();
			FlowDetachment();

			ToChannel();    // fraction going into channel
			OverlandFlow(); // slope kin wave
			ChannelFlow();  // channel erosion and kin wave

			Totals();
			MassBalance();

			OutputUI();
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
