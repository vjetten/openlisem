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
void TWorld::MassBalance(void)
{
// WATER in m3
     double areafac = 1000/CatchmentArea;
     tm->copy(Rain);
     tm->calc(DX, MUL);
     tm->calcV(_dx, MUL);
     double rainfall = tm->MapTotal(); // in m3
     RainTot += rainfall; // in m3
     RainTotmm = RainTot*areafac;

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

// SEDIMENT in kg
     if (SwitchErosion)
     {
        DetTotSplash += DETSplash->MapTotal();
        DetTotFlow += DETFlow->MapTotal();
        DepTot += DEP->MapTotal();
        DetTot += DETSplash->MapTotal() + DETFlow->MapTotal();
        SoilLossTot += Qsoutflow->MapTotal();
        SedVolTot = SedVol->MapTotal();
     }

// Channel
     if (SwitchIncludeChannel)
     {
       WaterVolTot += ChannelWaterVol->MapTotal(); //m3
       Qtot += ChannelQoutflow->MapTotal();

       if (SwitchErosion)
       {
          ChannelDetTot += ChannelDetFlow->MapTotal();
          ChannelDepTot += ChannelDep->MapTotal();
          ChannelSedTot = ChannelSedVol->MapTotal();

          SoilLossTot += ChannelQsoutflow->MapTotal();
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
	op.Qpeak=Qpeak*1000;
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

    if (runstep > 0 && runstep % 4 == 0)
    	printstep++;

	DEBUGv(printstep);
    // for display and write timeseries maps to disk

//	ReportTimeseries();
   	ReportTimeseriesNew();
	// report hydrographs ande swdigraophs at all points in outpoint.map

//	ReportTotals();
    // report totals to a text file
	ReportTotalsNew();


	ReportMaps();
	// report all maps and mapseries

}
//---------------------------------------------------------------------------
// the actual model with the main loop
void TWorld::DoModel()
{
    /*

	for (int i = 0; i < 1000; i++)
    {
    	DEBUGv(i);
    	WMap();
        mutex.lock();
        if(stopRequested) {emit debug("User interrupt...");}
        if(stopRequested) break;
        mutex.unlock();
    }
    */
  time_ms.start();
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

/*
void WMap()
{
	MAP *m, *out;
	REAL4 *mapData;
	m = Mopen("D:\\data\\jantiene\\Trier_17Oct03\\grad_tor_20m.map",M_READ);
	if (m == NULL)
		throw 1;
	UINT4 nrCols = RgetNrCols(m);
	UINT4 nrRows = RgetNrRows(m);


	(void)RuseAs(m, CR_REAL4);
	out = Rcreate("D:\\data\\jantiene\\Trier_17Oct03\\res\\try.map",RgetNrRows(m), RgetNrCols(m), CR_REAL4, VS_SCALAR,
			MgetProjection(m), RgetX0(m), RgetY0(m), RgetAngle(m), RgetCellSize(m));

	(void)RuseAs(out, CR_REAL4);
	mapData = (REAL4 *) Rmalloc(out, RgetNrCols(m));
	for (int c = 0; c < nrCols; c++)
		mapData[c] = 0;
	for(UINT4 r=0; r < nrRows; r++)
	{
		RgetRow(m, r, mapData);
	//	RputSomeCells(out,r*nrCols, nrCols, mapData);
			if (RputRow(out, r, mapData) != RgetNrCols(m))
			throw 1;
	}

	free(mapData);
	Mclose(out);
	Mclose(m);
}
*/
