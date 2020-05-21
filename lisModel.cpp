
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
  - void TWorld::DoModel() the main model function with the timeloop. It is a 'slot' linked to a signal.\n
  - void TWorld::run() Run is called from the interface to activate DoModel() \n
  - void TWorld::stop() Stops the loop on user request.\n

*/

#include <QtGui>
#include "lisemqt.h"
//#include "model.h"
#include "global.h"

/*
 * 1D hydrology
 * runoff 1D kin or 2D diffusive wave -> liserosion + flooding 2Ddyn + swofsediment
 * runoff/flood 2D dyn wave: swof + swof sediment
 */

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
void TWorld::saveMBerror2file(bool doError, bool start)
{
    if (doError && start) {
        //create error file
        QFile efout(resultDir+errorFileName);
        efout.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream eout(&efout);
        eout << "#water mass balance error (%)\n";
        eout << "2\n";
        eout << "run step\n";
        eout << "error\n";
        //  eout << "runtime\n";
        efout.flush();
        efout.close();

        if (SwitchErosion) {
            QFile esfout(resultDir+errorSedFileName);
            esfout.open(QIODevice::WriteOnly | QIODevice::Text);
            QTextStream esout(&esfout);
            esout << "#sediment mass balance error (%)\n";
            esout << "2\n";
            esout << "run step\n";
            esout << "MBs error\n";
            esfout.flush();
            esfout.close();
        }
    }


    if (doError) {
        QFile efout(resultDir+errorFileName);
        efout.open(QIODevice::Append | QIODevice::Text);
        QTextStream eout(&efout);
        eout << " " << runstep << " " << MB << /*" " << op.t <<*/ "\n";
        efout.flush();
        efout.close();

        if (SwitchErosion) {
            QFile esfout(resultDir+errorSedFileName);
            esfout.open(QIODevice::Append | QIODevice::Text);
            QTextStream esout(&esfout);
            esout << " " << runstep << " " << MBs <<  "\n";
            esfout.flush();
            esfout.close();
        }
    }

}

//---------------------------------------------------------------------------
// the actual model with the main loop
void TWorld::DoModel()
{
    if (!noInterface)
        temprunname = QString(op.LisemDir+"openlisemtmp.run");
    else
        temprunname = op.runfilename;

    timestampRun = QDateTime().currentDateTime().toString("yy.MM.dd-hh.mm");

    mapFormat = "PCRaster";

    //  QString DT = QDateTime().currentDateTime().toString("hh.mm-yy.MM.dd");
    errorFileName = QString(resultDir + "error"+ timestampRun +".txt");
    errorSedFileName = QString(resultDir + "errorsed"+ timestampRun +".txt");
    time_ms.start();
    // get time to calc run length


    //TODO: check material depth
    //TODO: check grainsize classes
    //TODO: check grass strips, just inthe input maps or special characteristics
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

        BeginTime = getvaluedouble("Begin time") * 60;
        EndTime = getvaluedouble("End time") * 60;
        _dt = getvaluedouble("Timestep");
        op.BeginTime = BeginTime/60; // for graph drawing
        op.EndTime = EndTime/60;// for graph drawing
        //VJ get time here else combomaps goes wrong for rainfall intensity

        //time vraiables in sec
        DEBUG("GetInputData()");
        GetInputData();
        DEBUG("IntializeData()");
        IntializeData();

        DEBUG("GetComboMaps()");
        GetComboMaps();

        DEBUG("setupDisplayMaps()");
        setupDisplayMaps();
        // reset all display output maps for new job
        // must be done after Initialize Data because then we know how large the map is

        if (SwitchRainfall)
        {
            DEBUG("GetRainfallData()");
            GetRainfallDataM(rainFileName, true);
        }
        if (SwitchSnowmelt)
        {
            DEBUG("GetSnowmeltData()");
            GetRainfallDataM(snowmeltFileName, false);
        }
        // get all input data and create and initialize all maps and variables

        CountLandunits();
        //VJ 110110 for output totals per landunit

        runstep = 0; //  runstep is used to initialize graph!
        printstep = 1; // printstep determines report frquency

        DEBUG("setupHydrographData()");
        setupHydrographData();

        saveMBerror2file(SwitchDumpMassBallance, true);

        InfilEffectiveKsat();
        // calc effective ksat from all surfaces once

        //start multithreading threadpool
        ThreadPool = new LisemThreadPool();
        ThreadPool->InitThreads(this);
        ThreadPool->SetMaskInitial(DEM);
        copy(*CoreMask, *ThreadPool->CoreMask);
        ThreadPool->StartReportThread(this);

        //create a function object referring to the cellprocesses wrapper
        wrapCellProcesses1D = std::bind((&TWorld::CellProcesses),this,std::placeholders::_1);
        fcompute2 = std::bind((&TWorld::CellProcesses2),this,std::placeholders::_1);

        DEBUG("Running...");

        for (time = BeginTime; time < EndTime; time += _dt)
        {
            if (runstep > 0 && runstep % printinterval == 0)
                printstep++;
            runstep++;

            if (noInterface && !noOutput)
            {
                int maxstep = static_cast <int>((EndTime-BeginTime)/_dt) ;
                qDebug() << runstep << maxstep << time/60 ;
            }

            mutex.lock();
            if(stopRequested) DEBUG("User interrupt... finishing time step");
            mutex.unlock();

            mutex.lock();
            if (waitRequested) DEBUG("User pause...");
            if (waitRequested) condition.wait(&mutex);
            mutex.unlock();
            // check if user wants to quit or pause

            //DEBUG(QString("Running timestep %1").arg((this->time - this->BeginTime)/_dt));

            //these functions read files, so they can not be multithreaded well
            RainfallMap();         // get rainfall from table or mpas
            SnowmeltMap();         // get snowmelt

            //do cell specific stuff, hydrology and splash detachment, threaded
            ThreadPool->RunCellCompute(wrapCellProcesses1D);
            ThreadPool->WaitForAll();

            ToChannel();           // overland flow water and sed flux going into or out of channel, in channel cells

            ToTiledrain();         // fraction going into tiledrain directly from surface

            OverlandFlow(); // overland flow 1D (non threaded), 2Ddiff or 2Ddyn (threaded), if 2Ddyn then also SWOFsediment!

            ChannelFlood(); // st venant channel 2D flooding from channel, only for kyn wave

            // flow detachment
     //        ThreadPool->RunCellCompute(fcompute2);
     //        ThreadPool->WaitForAll();

            OrderedProcesses();  //do ordered solutions such as channel LDD etc., non threaded

            //wait for the report thread that was started in the previous timestep
         //   ThreadPool->WaitForReportThread();

            Totals();            // calculate all totals and cumulative values

            MassBalance();       // check water and sed mass balance
            OutputUI();          // fill the "op" structure for screen output

            saveMBerror2file(SwitchDumpMassBallance, false);

//            std::function<void(int)> freport = std::bind((&TWorld::Wrapper_ReportAll),this,std::placeholders::_1);
//            ThreadPool->RunReportFunction(freport);
            mapFormat = op.format;

            ReportTimeseriesNew();
            // report hydrographs ande sedigraphs at all points in outpoint.map

            ReportTotalsNew();
            // report totals to a text file

            ReportMaps();
            // report all maps and mapseries

            ReportLandunits();
            // reportc stats per landunit class

            ChannelFloodStatistics();
            // report buildings submerged in flood level classes in 5cm intervals




            if (!noInterface)
                emit show();
            // send the op structure with data to function worldShow in LisUIModel.cpp

            if(stopRequested)
                time = EndTime;
        }

        //close the threads
        ThreadPool->WaitForReportThread();
        ThreadPool->WaitForAll();
        ThreadPool->Close();

        DestroyData();  // destroy all maps automatically
        DEBUG("Data destroyed");

        if (!noInterface)
            emit done("finished");
        else
        {
            if (!noOutput)
                qDebug() << "Done";

            QApplication::quit();
        }


    }
    catch(...)  // if an error occurred
    {
        DestroyData();

        if (!noInterface)
            emit done("ERROR STOP: "+ErrorString);
        else
        {
            qDebug() << "ERROR STOP: "+ErrorString;
            QApplication::quit();
        }
    }
}

void TWorld::CellProcesses(int thread)
{
    SetFlowBarriers(thread);     // update the presence of flow barriers
    GridCell(thread);            // set channel widths, flowwidths road widths etc

    Interception(thread);        // vegetation interception
    InterceptionLitter(thread);  // litter interception
    InterceptionHouses(thread);  // urban interception

    addRainfallWH(thread);       // adds rainfall to runoff water height or flood water height

    Infiltration(thread);        // infil of overland flow/flood water, decrease WH

    SoilWater(thread);           // simple soil water balance, percolation from lower boundary
    SurfaceStorage(thread);      // surface storage and flow width, split WH in WHrunoff and WHstore

    //doETa(thread);

    CalcVelDisch(thread);        // overland flow velocity, discharge and alpha for erosion

    SplashDetachment(thread);    // splash detachment

    FlowDetachment(thread);      // flow detachment, V used is from calcveldis for diff and kin, but not dynamic

    //Pestmobilisation();         // experimental
}

void TWorld::CellProcesses2(int thread) // obsolete for now
{
    CalcVelDisch(thread);        // overland flow velocity, discharge and alpha for erosion

    FlowDetachment(thread);      // flow detachment
}


void TWorld::OrderedProcesses()
{
    ChannelAddBaseandRainNT();  //NT is non threaded

    ChannelWaterHeightFromVolumeNT();

    CalcVelDischChannelNT(); // alpha, V and Q from Manning

    ChannelFlowDetachment();  //detachment, deposition for SS and BL

    ChannelFlow();         // channel erosion and kin wave

    TileFlow();          // storm/tile drain flow kin wave

    StormDrainFlow();

}

/* non threaded version for reference:
            SetFlowBarriers();     // update the presence of flow barriers
            GridCell();            // set channel widths, flowwidths road widths etc
            RainfallMap();         // get rainfall from table or mpas
            SnowmeltMap();         // get snowmelt
            Interception();        // vegetation interception
            InterceptionLitter();  // litter interception
            InterceptionHouses();  // urban interception
            addRainfallWH();       // adds rainfall to runoff water height or flood water height
            Infiltration();        // infil of overland flow water, decrease WH
            SoilWater();           // simple soil water balance, percolation from lower boundary
            SurfaceStorage();      // surface storage and flow width, split WH in WHrunoff and WHstore
            CalcVelDisch();        // overland flow velocity, discharge and alpha for erosion
            SplashDetachment();    // splash detachment
            ToFlood();             // overland flow water added to flood (not in channel cells)
            ToChannel();           // water and sed flux going into channel in channel cells
            ToTiledrain();         // fraction going into tiledrain directly from surface
            OverlandFlow();        // overland flow wave for water and sed
            CalcVelDisch();        // overland flow velocity, discharge and alpha for erosion
            FlowDetachment();      // flow detachment
            ChannelWaterHeight();  // add rainfall and runoff to channel and get channel WH from volume
            ChannelFlood();        // st venant channel 2D flooding from channel
            CalcVelDischChannel(); // alpha, V and Q from Manning
            ChannelFlow();         // channel erosion and kin wave
            TileFlow();            // tile drain flow kin wave
            Totals();              // calculate all totals and cumulative values
            MassBalance();         // check water and sed mass balance
            */
