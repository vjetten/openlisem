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
     double areafac = 1000/(DX->MapTotal()*_dx);
     tm->copy(Rain);
     //RainTotmm += tm->MapTotal() *1000/nrCells; // avg in mm
     tm->calc(DX, MUL);
     tm->calcV(_dx, MUL);
     RainTot += tm->MapTotal(); // in m3
     RainTotmm = RainTot*areafac;

     IntercTot = Interc->MapTotal();
     IntercTotmm = IntercTot*areafac;

     InfilTot += InfilVol->MapTotal() + InfilVolKinWave->MapTotal(); //m3
     InfilKWTot += InfilVolKinWave->MapTotal();
     InfilTotmm = InfilTot*areafac;

     tm->calc2(WHstore, FlowWidth, MUL); //m2
     tm->calcV(_dx, MUL); //m3
     SurfStorTot = tm->MapTotal();

     WaterVolTot = WaterVol->MapTotal();//m3
     WaterVolTotmm = WaterVolTot*areafac;

     Qtot += Qoutflow->MapTotal();
     // sum all outflow m3 for all timesteps, is already mult by dt
     Qtotmm = Qtot*areafac;

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
          DetTot += ChannelDetFlow->MapTotal();
          DepTot += ChannelDep->MapTotal();
          SoilLossTot += ChannelQsoutflow->MapTotal();
          SedVolTot += ChannelSedVol->MapTotal();
       }
     }

// Mass Balance
     if (RainTot > 0)
       MB = (RainTot - IntercTot - InfilTot - WaterVolTot - Qtot)/RainTot*100;

     if (SwitchErosion && DetTot > 0)
       MBs = (DetTot - SoilLossTot - SedVolTot + DepTot)/DetTot*100;
}
//---------------------------------------------------------------------------
void TWorld::Output()
{
        runstep++;
        // for display and write timeseries maps to disk

	op.MB = MB;
	op.runstep = runstep;
	op.maxstep = (int) ((EndTime-BeginTime)/_dt);

	op.RainTot=RainTotmm;
	op.WaterVolTot=WaterVolTotmm;
	op.Qtotmm=Qtotmm;
	op.Qtot=Qtot;
	op.Qpeak=Qpeak;
	op.InfilTot=InfilTotmm;
	op.IntercTot=IntercTotmm;
	op.InfilKWTot=InfilKWTot;

	op.MBs = MBs;
	op.DetTotSplash=DetTotSplash*0.001;
	op.DetTotFlow=DetTotFlow*0.001;
	op.SoilLossTot=SoilLossTot*0.001;
	op.SedVolTot=SedVolTot*0.001;
	op.DepTot=DepTot*0.001;
	op.t = time_ms.elapsed()*0.001/60.0;
	op.time = time/60;
	op.maxtime = op.t/runstep * op.maxstep;

	emit show(runstep);

	ReportTimeseries();

	//fpot->report("fpot",runstep);
	//fact->report("fact",runstep);
	//   RainCum->report("rainc", runstep);
    // Fcum->report("fcum", runstep);
    // WHstore->report("sstor", runstep);
    //  Qn->report("Qn", runstep);
}//---------------------------------------------------------------------------
void TWorld::DoModel()
{
  time_ms.start();
  temprunname = op.runfilename;

  try
  {
     emit debug("reading data");

     IntializeOptions(); // all switches to false, clear names
     GetRunFile();
     QString sss;
     sss.setNum(nrnamelist);
     DEBUG(sss+ " variables read from runfile");
     ParseInputData();

    // SwitchIncludeChannel = op.SwitchIncludeChannel;// from interface
    // SwitchErosion = op.SwitchErosion;

     InitMapList();
     GetInputData();
     IntializeData();
     GetRainfallData();
     BeginTime = getvaluedouble("Begin time") * 60; //read min and convert to sec
     EndTime = getvaluedouble("End time") * 60;
     _dt = getvaluedouble("Timestep");

     emit debug("running");

     runstep = 0;
     for (time = BeginTime; time < EndTime; time += _dt)
     {
       mutex.lock();
       if(stopRequested) break;
       mutex.unlock();

       GridCell();
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

    DestroyData();

    emit done("finished");

  }
  catch(int i)
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
