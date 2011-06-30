
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
  \file lisModel.cpp
  \brief Central model file with the main loop. From here all processes are called.

  functions: \n
  - void TWorld::OutputUI() fill output structure 'op' with results to talk to the interface.\n
  - void TWorld::DoModel() the main model function with the timeloop. It is a 'slot' linked to a signal.\n
  - void TWorld::run() Run is called from the interface to activate DoModel() \n
  - void TWorld::stop() Stops the loop on user request.\n

*/


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
/** fill output structure 'op' with results to talk to the interface:
    report to screen, hydrographs and maps */
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
    
	op.RainTotmm=RainTotmm + SnowTotmm;
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
    
   op.P = (RainAvgmm + SnowAvgmm)*3600/_dt;
   op.Q = Qoutput->DrcPlot; //Outlet;
   op.Qs = Qsoutput->DrcPlot; //Outlet;
   op.C = TotalConc->DrcPlot; //Outlet;
   op.Qtile = TileQoutflow->DrcPlot; //Outlet;
   // VJ 110630 show hydrograph for selected output point
    
	op.BufferVolTot = BufferVolTot;
	op.BufferSedTot = BufferSedTot*0.001; //ton
    
	emit show();
   // send the op structure with data to function worldShow the interface
	// in file LisUIModel
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
      DEBUG("ParseRunfileData()");
      ParseRunfileData();
		// get and parse runfile
        
		DEBUG("GetInputData()");
		GetInputData();
		DEBUG("IntializeData()");
		IntializeData();
		if (SwitchRainfall)
		{
            DEBUG("GetRainfallData()");
            GetRainfallData();
        }
		if (SwitchSnowmelt)
		{
			DEBUG("GetSnowmeltData()");
			GetSnowmeltData();
		}
		// get all input data and create and initialize all maps and variables
        
      CountLandunits();
      //VJ 110110 for output totals per landunit

		BeginTime = getvaluedouble("Begin time") * 60;
		EndTime = getvaluedouble("End time") * 60;
		_dt = getvaluedouble("Timestep");
		op.BeginTime = BeginTime/60; // for graph drawing
		op.EndTime = EndTime/60;// for graph drawing
		//time vraiables in sec
        
		runstep = 0; // NOTE runstep is used to initialize graph!
		printstep = 1;

      // VJ 110630 show hydrograph for selected output point
      FOR_ROW_COL_MV
      {
         if (op.outputpointnr == PointMap->Drc)
         {
            r_plot = r;
            c_plot = c;
            op.outputpointdata = QString("Hydrograph point %1 [row %2; col %3]").arg(op.outputpointnr).arg(r).arg(c);
            if( op.outputpointnr == 1)
               op.outputpointdata = QString("Hydrograph main outlet");
         }
      }

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
         ToTiledrain();    // fraction going into tiledrain directly from surface

         OverlandFlow(); // slope kin wave
			ChannelFlow();  // channel erosion and kin wave
         TileFlow();     // tile drain flow kin wave

			Totals();
			MassBalance();
            
			OutputUI();

         reportAll(); //VJ 110114 now separate
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
