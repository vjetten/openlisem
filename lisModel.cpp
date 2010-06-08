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
//#include "ifacebasic.h"
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
	int i = 1;
}
//---------------------------------------------------------------------------
// totals for screen and file output and mass balance
void TWorld::Totals(void)
{
	//***** WATER *****//

	RainAvgmm = Rain->MapTotal()*1000/nrCells;
	RainTotmm = RainTotmm + RainAvgmm;
	// avg area rainfall in mm

	tm->calc2V(Rain, (_dx*_dx), MUL); //in m3
	double rainfall = tm->MapTotal(); // in m3
	RainTot += rainfall; // in m3

	double oldpeak = Rainpeak;
	Rainpeak = max(Rainpeak, rainfall);
	if (oldpeak < Rainpeak)
		RainpeakTime = time;

	IntercTot = Interc->MapTotal();
	IntercTotmm = IntercTot*1000/(_dx*_dx*nrCells);//CatchmentArea;
	// interception in mm and m3

	InfilTot += InfilVol->MapTotal() + InfilVolKinWave->MapTotal(); //m3
	InfilKWTot += InfilVolKinWave->MapTotal(); // not really used, available for output when needed
	InfilTotmm = max(0,InfilTot*1000/(_dx*_dx*nrCells));//CatchmentArea;
	// infiltration mm and m3

	tm->calc2V(WHstore, 1000, MUL); //mm
	SurfStoremm = tm->MapTotal()/nrCells;
	// surface storage CHECK THIS
	// does not go to MB, is already in tot water vol

	WaterVolTot = WaterVolall->MapTotal();//m3
	WaterVolTotmm = WaterVolTot*1000/(_dx*_dx*nrCells);//CatchmentArea; //mm
	// water on the surface in runoff in m3 and mm
	//NOTE: surface storage is already in here so does not need to be accounted for in MB

	Qtot += Qoutflow->MapTotal();
	// sum outflow m3 for all timesteps for all pits, is already mult by dt
	// needed for mass balance
	Qtotmm = Qtot*1000/(_dx*_dx*nrCells);//CatchmentArea;
	// in mm for screen output
/*
	FOR_ROW_COL_MV
	{
		if (Outlet->Drc == 1)
			QtotOutlet += Qoutflow->Drc;
		// for screen output, total main outlet in m3
		TotalWatervol->Drc = WaterVolall->Drc;
		// for sed conc calc output
	}
	*/
	QtotOutlet += Qoutflow->DrcOutlet;
	// for screen output, total main outlet in m3
	TotalWatervol->copy(WaterVolall);
	// for sed conc calc output


	if (SwitchIncludeChannel)
	{
		WaterVolTot += ChannelWaterVol->MapTotal(); //m3
		// add channel vol to total
		WaterVolTotmm = WaterVolTot*1000/(_dx*_dx*nrCells);//CatchmentArea; //mm
		// recalc in mm for screen output

		Qtot += ChannelQoutflow->MapTotal();
		// add channel outflow (in m3) to total for all pits
		Qtotmm = Qtot*1000/(_dx*_dx*nrCells);//CatchmentArea;
		// recalc in mm for screen output

		QtotOutlet += ChannelQoutflow->DrcOutlet;
		// add channel outflow (in m3) to total for main outlet
		TotalWatervol->calc(ChannelWaterVol,ADD);
		// add channel volume to total for sed conc calc

		/*
		FOR_ROW_COL_MV
		{
			if (Outlet->Drc == 1)
				QtotOutlet += ChannelQoutflow->Drc;
			// add channel outflow (in m3) to total for main outlet

			TotalWatervol->Drc += ChannelWaterVol->Drc;
			// add channel volume to total for sed conc calc
		}
		*/
	}

	if (SwitchBuffers)
	{
		BufferVolTot = BufferVol->MapTotal(); // in m3
		if (SwitchIncludeChannel)
			BufferVolTot += ChannelBufferVol->MapTotal();
		//sum up all volume remaining in all buffers (so the non-water!)
		BufferVolTot = BufferVolTotInit - BufferVolTot;
		//subtract this from the initial volume to get the water in the buffers
	}

	// output fluxes for reporting
	FOR_ROW_COL_MV
	{
		Qoutput->Drc = 1000*(Qn->Drc + ChannelQn->Drc); // in l/s
	}

	oldpeak = Qpeak;
	Qpeak = max(Qpeak, Qoutput->DrcOutlet);
	if (oldpeak < Qpeak)
		QpeakTime = time;
	// peak flow and peak time calculation, based on sum channel and runoff

	//***** SEDIMENT *****//

	if (SwitchErosion)
	{
		DetSplashTot += DETSplash->MapTotal();
		DetFlowTot += DETFlow->MapTotal();
		DepTot += DEP->MapTotal();
		DetTot += DETSplash->MapTotal() + DETFlow->MapTotal();
		SedTot = Sed->MapTotal();
		// all in kg/cell

		SoilLossTot += Qsoutflow->MapTotal();
		// sum all sed in all pits (in kg), needed for mass balance

		SoilLossTotOutlet += Qsoutflow->DrcOutlet;
		// for screen output, total main outlet sed loss
		TotalSed->copy(Sed);
		// for sed conc

		if (SwitchIncludeChannel)
		{
			ChannelDetTot += ChannelDetFlow->MapTotal();
			ChannelDepTot += ChannelDep->MapTotal();
			ChannelSedTot = ChannelSed->MapTotal();

			SoilLossTot += ChannelQsoutflow->MapTotal();
			// add sed outflow for all pits to total soil loss

			SoilLossTotOutlet += ChannelQsoutflow->DrcOutlet;
				// add channel outflow (in kg) to total for main outlet

			TotalSed->calc(ChannelSed, ADD);
				// for sed conc file output
		}

		FOR_ROW_COL_MV
			TotalConc->Drc = (TotalWatervol->Drc > 0? TotalSed->Drc/TotalWatervol->Drc : 0);
		// for file output

		if (SwitchBuffers)
		{
			BufferSedTot = BufferSed->MapTotal();
			if (SwitchIncludeChannel)
				BufferSedTot += ChannelBufferSed->MapTotal();
			BufferSedTot = BufferSedTotInit - BufferSedTot;
		}
		//TODO add gully, wheeltracks etc

		// spatial totals for output
		FOR_ROW_COL_MV
		{
			Qsoutput->Drc = Qsn->Drc + ChannelQsn->Drc;  // sum channel and OF sed output in kg/s

			TotalDetMap->Drc += DETSplash->Drc + DETFlow->Drc;
			TotalDepMap->Drc += DEP->Drc;
			if (SwitchIncludeChannel)
			{
				TotalDetMap->Drc += ChannelDetFlow->Drc;
				TotalDepMap->Drc += ChannelDep->Drc;
			}
			TotalSoillossMap->Drc = TotalDetMap->Drc + TotalDepMap->Drc;
		}
	}
}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
	// Mass Balance water
	if (RainTot > 0)
		MB = (RainTot - IntercTot - InfilTot - WaterVolTot  - BufferVolTot- Qtot)/RainTot*100;

	// Mass Balance sediment
	if (SwitchErosion && DetTot > 0)
		MBs = (DetTot + ChannelDetTot - SoilLossTot - SedTot - ChannelSedTot + DepTot + ChannelDepTot - BufferSedTot)/DetTot*100;
}
//---------------------------------------------------------------------------
// fill output structure to talk to interface
// op is declared in ifacebasic
// report to screen, hydrographs and maps
void TWorld::Output()
{
	if (runstep > 0 && runstep % printinterval == 0)
		printstep++;
	// printstep determines map reporting frequency
	runstep++;

	op.dx = _dx;
	op.SwitchErosion = SwitchErosion;
	op.SwitchIncludeChannel = SwitchIncludeChannel;
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

		GetRunFile();
		ParseInputData();
		// get and parse runfile
		GetInputData();
		IntializeData();
		GetRainfallData();
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
			if(stopRequested) {emit debug("User interrupt...");}
			if(stopRequested) break;
			mutex.unlock();
			// check if user wants to quit

			GridCell();          // set channel widths, flowwidths road widths etc
			Rainfall();
			Interception();
			Infiltration();
			//SoilWater();
			SurfaceStorage();
			SplashDetachment();
			FlowDetachment();

			ToChannel();    // fraction going into channel
			OverlandFlow(); // slope kin wave
			ChannelFlow();  // channel erosion and kin wave

			Totals();
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
