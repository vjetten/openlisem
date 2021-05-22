/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
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
    if (!op.doBatchmode)
        temprunname = QString(op.LisemDir+"openlisemtmp.run");
    else
        temprunname = op.runfilename;

    mapFormat = "PCRaster";

    errorFileName = QString(resultDir + "error-"+ op.timeStartRun +".txt");
    errorSedFileName = QString(resultDir + "errorsed-"+ op.timeStartRun +".txt");
    time_ms.start();
    // get time to calc run length
    startTime=omp_get_wtime()/60.0;

    //qDebug() << SwitchErosion << Switch2DDiagonalFlow << F_pitValue;

    try
    {
        DEBUG("reading and initializing data");

        IntializeOptions(); // reset all options

        InitMapList();
        // map structure to destroy data automatically

        DEBUG("GetRunFile()");
        GetRunFile();
        DEBUG("ParseRunfileData()");
        ParseRunfileData();
        // get and parse runfile

        //BeginTime = getTimefromString(bt)*60; // in seconds!
        //EndTime = getTimefromString(et)*60;
        double btd = getvaluedouble("Begin time day");
        double btm = getvaluedouble("Begin time");
        double etd = getvaluedouble("End time day");
        double etm = getvaluedouble("End time");
        BeginTime = (btd*1440+btm)*60; //in sec
        EndTime = (etd*1440+etm)*60;   //in sec
           // qDebug() << "time" << btd << btm << etd << etm;
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

        SwitchSnowmelt = false;
        if (SwitchRainfall)
        {
            DEBUG("Get Rainfall Data Information");
            if (SwitchRainfallSatellite) {
                GetSpatialMeteoData(rainSatFileName, 0);
            }
            else
                GetRainfallData(rainFileName);
          //  op.maxRainaxis = getmaxRainfall();
        }
        if (SwitchIncludeET)
        {
            DEBUG("Get ET Data Information");
            if (SwitchETSatellite)
                GetSpatialMeteoData(ETSatFileName, 1);
            else
                GetRainfallData(ETFileName);
        }
        if (SwitchSnowmelt)
        {
            DEBUG("Get Snowmelt Data Information");
            if (SwitchSnowmeltSatellite)
                GetSpatialMeteoData(SnowmeltSatFileName, 2);
            else
                GetSnowmeltData(snowmeltFileName);
        }


        if (SwitchChannelInflow)
        {
            DEBUG("GetDischargeData()");
            GetDischargeData(dischargeinFileName);
        }

        // get all input data and create and initialize all maps and variables

        CountLandunits();
        //VJ 110110 for output totals per landunit

        runstep = 0; //  runstep is used to initialize graph!
        printstep = 1; // printstep determines report frquency

        DEBUG("setupHydrographData()");
        setupHydrographData();

        bool saveMBerror = false;
        saveMBerror2file(saveMBerror, true);

        InfilEffectiveKsat();  // calc effective ksat from all surfaces once
        SetFlowBarriers();     // update the presence of flow barriers, static for now, unless breakthrough
        GridCell();            // static for now

        DEBUG("Running...");
        for (time = BeginTime; time < EndTime; time += _dt)
        {

            if (runstep > 0 && runstep % printinterval == 0)
                printstep++;
            runstep++;

            if(stopRequested) {
                mutex.lock();
                DEBUG("User interrupt... finishing time step");
                mutex.unlock();
            }

            if (waitRequested) {
                mutex.lock();
                DEBUG("User pause...");
                condition.wait(&mutex);
                mutex.unlock();
            }
            // check if user wants to quit or pause

            CellProcesses();
         //   do_CellProcesses();

            ToTiledrain();         // fraction going into tiledrain directly from surface

            OverlandFlow();        // overland flow 1D (non threaded), 2Ddyn (threaded), if 2Ddyn then also SWOFsediment!

            OrderedProcesses();    // do ordered LDD solutions channel, tiles, drains, non threaded

            Totals();            // calculate all totals and cumulative values

            MassBalance();       // check water and sed mass balance

            OutputUI();          // fill the "op" structure for screen output and calc some output maps

            reportAll();

            emit show(noInterface); // send the 'op' structure with data to function worldShow in LisUIModel.cpp

            saveMBerror2file(saveMBerror, false);

            if(stopRequested)
                time = EndTime;
        }

        if (SwitchEndRun)
            ReportMaps();

        DestroyData();  // destroy all maps automatically        
        DEBUG("Data destroyed");

        emit done("finished");

        //???????????
        if(batchmode)
            QApplication::quit();
    }
    catch(...)  // if an error occurred
    {
        DestroyData();

//        if (!noInterface)
            emit done("ERROR STOP: "+ErrorString);
//        else
//        {
//            qDebug() << "ERROR STOP: "+ErrorString;
//          if(batchmode)
//                QApplication::quit();
//        }
    }
}

void TWorld::CellProcesses()
{
    if (SwitchRainfallSatellite)
        GetRainfallSatMap();         // get rainfall from table
    else
        GetRainfallMap();         // get rainfall from maps

    if (SwitchIncludeET) {
        if (SwitchETSatellite)
            GetETSatMap();
    }


    GetSnowmeltMap();         // get snowmelt

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (Rainc->Drc > 0)
            cell_Interception(r,c);
        if (FloodDomain->Drc > 0) {
            hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        } else {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
            // add net to water rainfall on soil surface (in m)

            //    if (SwitchGrassStrip && GrassWidthDX->Drc > 0)
            //        WHGrass->Drc += RainNet->Drc + Snowmeltc->Drc;
            // net rainfall on grass strips, infil is calculated separately for grassstrips
        }
        if (RoadWidthHSDX->Drc > 0)
            WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;

        if (InfilMethod == INFIL_SWATRE) {
           cell_InfilSwatre(r, c);
        }
        else
        if (InfilMethod != INFIL_NONE) {
           cell_InfilMethods(r, c);
           if (!SwitchImpermeable)
               cell_Percolation(r, c);

        }

        if (SwitchIncludeET) {
            cell_ETa(r, c);
        }

        double wh = WH->Drc;
        double SW = SoilWidthDX->Drc;
        double RW = RoadWidthHSDX->Drc;
        double WHr = WHroad->Drc;
        double WHs;
        //### surface storage on rough surfaces
        WHs = std::min(wh, MDS->Drc*(1-exp(-1.875*(wh/std::max(0.01,0.01*RR->Drc)))));
        // non-linear release fo water from depression storage
        // resemles curves from GIS surface tests, unpublished
        double FW = std::min(ChannelAdj->Drc, SW + RW);
        // calculate flowwidth by fpa*surface + road, excludes channel already

        WHrunoff->Drc = ((wh - WHs)*SW + WHr*RW)/FW;
        FlowWidth->Drc = FW;

        WaterVolall->Drc = DX->Drc*(wh*SW + WHr*RW);
        // all water in the cell incl storage
        WHstore->Drc = WHs;

        if (SwitchErosion) {
            double wh = FloodDomain->Drc == 0 ? WH->Drc : hmx->Drc;
            cell_SplashDetachment(r,c,wh);
        }


    }}
//    Interception();        // vegetation interception
//    InterceptionLitter();  // litter interception
//    InterceptionHouses();  // urban interception

//    addRainfallWH();       // adds rainfall to runoff water height or flood water height

//    Infiltration();        // infil of overland flow/flood water, decrease WH

//    SoilWater();           // simple soil water balance, percolation from lower boundary

//    SurfaceStorage();      // surface storage and flow width, split WH in WHrunoff and WHstore

//    doETa();

//    Pestmobilisation();         // experimental

//    SplashDetachment();    // splash detachment

}


// these are all non-threaded
void TWorld::OrderedProcesses()
{
    SwitchChannelKinWave = true;// set to false for experimental swof in channel

    ChannelAddBaseandRain();  // add baseflow o, subtract infil, add rainfall

    ChannelFlow();            //channel kin wave for water and sediment

    ChannelFlowDetachmentNew();  //detachment, deposition for SS and BL

    TileFlow();          // tile drain flow kin wave

    StormDrainFlow();    // storm drain flow kin wave

}

