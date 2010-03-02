//---------------------------------------------------------------------------
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
// WATER

     tm->copy(Rain);
     tm->calc(DX, MUL);
     tm->calcV(_dx, MUL);
     RainTot += tm->MapTotal();

     IntercTot = Interc->MapTotal();

     InfilTot += InfilVol->MapTotal() + InfilVolKinWave->MapTotal();
     InfilKWTot += InfilVolKinWave->MapTotal();

     tm->calc2(WHstore, FlowWidth, MUL); //m2
     tm->calcV(_dx, MUL); //m3
     SurfStorTot = tm->MapTotal();

     WaterVolTot = WaterVol->MapTotal();
     //m3

     Qtot += Qoutflow->MapTotal();
     // sum all outflow m3/s for all timesteps

// SEDIMENT
     if (SwitchErosion)
     {
        DetTotSplash += DETSplash->MapTotal();
        DetTotFlow += DETFlow->MapTotal();// + (SwitchIncludeChannel ? ChannelDetFlow->MapTotal(): 0);
        DepTot += DEP->MapTotal();// + (SwitchIncludeChannel ? ChannelDep->MapTotal() : 0);
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

	op.RainTot=RainTot;
	op.WaterVolTot=WaterVolTot;
	op.Qtot=Qtot;
	op.InfilTot=InfilTot;
	op.IntercTot=IntercTot;
	op.InfilKWTot=InfilKWTot;

	op.MBs = MBs;
	op.DetTotSplash=DetTotSplash*0.001;
	op.DetTotFlow=DetTotFlow*0.001;
	op.SoilLossTot=SoilLossTot*0.001;
	op.SedVolTot=SedVolTot*0.001;
	op.DepTot=DepTot*0.001;
	op.t = time_ms.elapsed()*0.001;
	op.time = time/60;

	emit show(runstep);

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
     ParseInputData();

     SwitchIncludeChannel = op.SwitchIncludeChannel;// from interface
     SwitchErosion = op.SwitchErosion;

     InitMapList();
     GetInputData();
     IntializeData();
     GetRainfallData();
     BeginTime = getvaluedouble("Begin time") * 60; //read min and convert to sec
     EndTime = getvaluedouble("End time") * 60;
     _dt = getvaluedouble("Timestep");

     emit debug("running");
     DEBUG(resultDir);
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
