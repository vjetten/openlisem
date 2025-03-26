/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
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
void TWorld::saveMBerror2file( bool start) //bool doError,
{
    if (start) {
        //create error file
        QFile efout(resultDir+errorFileName);
        efout.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream eout(&efout);
        eout << "#mass balance error (%)\n";

        efout.flush();
        efout.close();
    }


 //   if (doError) {
        QFile efout(resultDir+errorFileName);
        efout.open(QIODevice::Append | QIODevice::Text);
        QTextStream eout(&efout);
        eout << " " << runstep << "," << MB << "," << (SwitchErosion ? MBs : 0.0) << "\n";
        efout.flush();
        efout.close();
  //  }

}
//---------------------------------------------------------------------------
// the actual model with the main loop
void TWorld::DoModel()
{
    if (!op.doBatchmode)
        temprunname = QString(op.userAppDir+"openlisemtmp.run");
    else
        temprunname = op.runfilename;

    mapFormat = "PCRaster";

    errorFileName = QString(resultDir + "error-"+ op.timeStartRun +".csv");
    //errorSedFileName = QString(resultDir + "errorsed-"+ op.timeStartRun +".txt");
    time_ms.start();
    // get time to calc run length
    startTime=omp_get_wtime()/60.0;

    ETafactorTot = 0;


    try
    {
        DestroyData();

        DEBUG("reading and initializing data");

        IntializeOptions(); // reset all options

        DEBUG("GetRunFile()");
        GetRunFile();
        DEBUG("Parse Runfile");
        ParseRunfileData();
        // get and parse runfile

        QString S = resultDir + QFileInfo(op.runfilename).fileName();
        QFile::copy(op.runfilename, S);


        //time vraiables in sec
        // double btd = getvaluedouble("Begin time day");
        // double btm = getvaluedouble("Begin time");
        // double etd = getvaluedouble("End time day");
        // double etm = getvaluedouble("End time");
        double btd, etd, btm, etm;
        QString beginTimeString = getvaluestring("Begin Time");
        QString endTimeString = getvaluestring("End Time");
        bool dayOk, minuteOk;
        QStringList parts = beginTimeString.split(":");
        if (parts.size() == 2) {
            btd = parts[0].toInt(&dayOk);
            btm = parts[1].toInt(&minuteOk);
        }

        parts = endTimeString.split(":");
        if (parts.size() == 2) {
            etd = parts[0].toInt(&dayOk);
            etm = parts[1].toInt(&minuteOk);
        }

        btd -= 1.0; // because day 1, minute 10 is in fact minute 10 in the first day
        etd -= 1.0;

        if (SwitchEventbased) {
            DEBUG("Day in start and end time is ignored.");
        }

        _dt = getvaluedouble("Timestep");

        if (SwitchEventbased) {
            BeginTime = (btm)*60; //for running in sec
            EndTime = (etm)*60;   //in sec
            op.BeginTime = BeginTime/60; // for graph drawing in min
            op.EndTime = EndTime/60;
        } else {
            BeginTime = (btd*1440+btm)*60; //for running in sec
            EndTime = (etd*1440+etm)*60;   //in sec
            op.BeginTime = BeginTime/60;// for graph drawing in min
            op.EndTime = EndTime/60;
        }

        //get all maps
        DEBUG("Get Input Maps");
        GetInputData();
        DEBUG("Intialize Database");
        IntializeData();

        // MC - no_ui probalbly this can be skipped for noInterface??
        setupDisplayMaps();
        // reset all display output maps for new job
        // must be done after Initialize Data because then we know how large the map is
        // clear() calls the destruction of all elements in the sturcture

        if (SwitchRainfall)
        {
            RainfallSeries.clear();
            RainfallSeriesMaps.clear();
            raintime.clear();
            DEBUG("Get Rainfall Data");
            if (SwitchRainfallSatellite) {
                GetSpatialMeteoData(rainSatFileName, 0);
            } else {
                GetRainfallStationData(rainFileName);
            }
        }

        if (SwitchIncludeET)
        {
            ETSeries.clear();
            ETSeriesMaps.clear();
            ETtime.clear();

            DEBUG("Get EvapoTranspiaration Data");
            if (SwitchETSatellite) {
                GetSpatialMeteoData(ETSatFileName, 1);
            } else {
                GetETStationData(ETFileName);
            }
        }

        SwitchSnowmelt = false;
        // if (SwitchSnowmelt)

        // }

        if (SwitchDischargeUser)
        {
            DischargeSeries.clear();
            dischargetime.clear();

            DEBUG("GetUserDischargeData()");
            GetUserDischargeData(dischargeinFileName);
        }

        if (SwitchWaveUser)
        {
            WHSeries.clear();
            WHtime.clear();

            DEBUG("GetWHboundaryData()");
            GetWHboundaryData(WaveinFileName);
        }

        // get all input data and create and initialize all maps and variables

        CountLandunits();
        //VJ 110110 for output totals per landunit

        runstep = 0; //  runstep is used to initialize graph!
        printstep = 1; // printstep determines report frequency in report()

        DEBUG("setupHydrographData()");
        setupHydrographData(); // reset hydrograph display

        //bool saveMBerror = true;
        //saveMBerror2file(true); //saveMBerror,

        SetFlowBarriers();     // update the presence of flow barriers, static for now, unless breakthrough
        GridCell();            // static for now

        _dt_user = _dt;

        DEBUG(" ");

        GetComboMaps(); // moved to outside timeloop!

        InfilEffectiveKsat();

        // ---- THE TIME LOOP ----
        for (time = BeginTime; time < EndTime; time += _dt)
        {            
            // printstep determines report frequency in #define report(...)
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

            GetInputTimeseries(); // get rainfall, ET, snowmelt, discharge

            InfilDynamicCrusting(); // if crusting recalc Ksateff and Poreff becuase of crusting effect

            HydrologyProcesses();  // hydrological processes in one loop, incl splash

            OverlandFlow(); // overland flow 1D (non threaded), 2Ddyn (threaded), if 2Ddyn then also SWOFsediment!

            // these are all non-threaded
            ChannelFlowandErosion();    // do ordered LDD solutions channel, tiles, drains, non threaded

            TileFlow();          // tile drain flow kin wave
                                 // storm drain flow kin wave
            //StormDrainFlow();

            TotalsHydro();       // calculate all totals and cumulative values
            TotalsFlow();
            TotalsSediment();

            MassBalance();       // check water and sed mass balance

            reportToUI();          // fill the "op" structure for screen output and calc some output maps

            reportToFile();         // report hydrograohs, totals, maps etc to files

            emit show(noInterface); // send the 'op' structure with data to function worldShow in LisUIModel.cpp

            //saveMBerror2file(false); //saveMBerror

            if(stopRequested)
                time = EndTime;

            // show progress in console without GUI
            if (op.doBatchmode) {
                int x;
                x = std::round((op.t / op.maxtime) * 100) ;
                printf("\rprogress: %d %%                     ", x);
                // or use qDebug()
            }
             // MC - maybe not the most sophisticated solution but noInterface works again
        }

        if (SwitchEndRun)
            ReportMaps();

        emit done("finished");

        if (op.doBatchmode)
        {
            // delete all maps
            qDeleteAll(maplistCTMap.begin(),maplistCTMap.end());
            maplistCTMap.clear();

            //delete swatre 3D soil layer structure if exists
            if (initSwatreStructure)
                FreeSwatreInfo();

            qDebug() << "\nfinished after "<< op.maxtime << "minutes\n";
            if (noInterface)
                QCoreApplication::quit();
            else
                QApplication::quit();
            // close the world model
        }
    }
    catch(...)  // if an error occurred
    {
        emit done("ERROR STOP: "+ErrorString);
        if (op.doBatchmode) {qDebug() << "ERROR STOP "<< ErrorString;
            if (noInterface)
                QCoreApplication::quit();
            else
                QApplication::quit();
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::GetInputTimeseries()
{
    // get meteo data
    if(SwitchRainfall) {
        if (SwitchRainfallSatellite)
            GetRainfallMapfromSat(time);         // get rainfall from maps
        else
            GetRainfallMapfromStations(time);  // get rainfall from stations
    }

    if (SwitchIncludeET) {
        if (SwitchETSatellite)
            GetETSatMap(time); // get rainfall from maps
        else
            GetETMapfromStations(time);   // get rainfall from stations
    }

    if (SwitchDischargeUser) {
        GetDischargeMapfromStations(time);
    }

    if (SwitchWaveUser) {
        GetWHboundaryMap(time);
    }

//    if (SwitchSnowmelt) {
//        if (SwitchSnowmeltSatellite)
//            ; //TODO snowmelt satellite
//        else
//            GetSnowmeltMap(time);  // get snowmelt from stations
//    }

}
//---------------------------------------------------------------------------
// all hydrologuical processes in one big parallel loop for speed
void TWorld::HydrologyProcesses()
{
   // double soiltot1 = SoilWaterMass();

    if (SwitchIncludeET) {
        if (SwitchDailyET)
            ETafactor = getETaFactor(); // based on daylength if daily values, converts from m/day to m/timestep directly
        else
            ETafactor = 1.0;   // if not ETfactor can be 1,.0 because ET is already in m/timestep
    }

    // above ground
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (SwitchInterception)
            cell_Interception(r,c);
        // all interception on plants, houses, litter
        // result is rainnet (and leafdrip for erosion)

        if (SwitchIncludeET)
            cell_ETa(r,c);
        // interception and soil surface evap, also ET from Green and Ampt, not SWATRE

        // floododmain is used if kinwave + overflow to separate WH runoiff from 2D hmx flood
        if (FloodDomain->Drc > 0) {
            hmx->Drc += RainNet->Drc;// + Snowmeltc->Drc; // only used in kin wave plus flood from channel, hmx is flood water
        } else {
            WH->Drc += RainNet->Drc;// + Snowmeltc->Drc;  // used in 2D flow and kin wave
        }

        // incoming wave at boundary
        if (SwitchWaveUser) {
            WHboundRain->Drc += RainNet->Drc;
            if (WHboundarea->Drc > 0) {
                // WHbound is the forced water level in area with value '1', plus cum rainfall
                WH->Drc = WHbound->Drc + WHboundRain->Drc;
            }
        }
    }}

    if (SwitchInfiltration) {
        // non SWATRE infiltration, redistribution and percolation
        if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_SOAP) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                cell_InfilMethods(r, c);

                if (SwitchTwoLayer) {
                    cell_Redistribution2(r, c);
                    //cell_Channelinfow2(r, c);
                } else {
                    cell_Redistribution1(r, c);
                    //cell_Channelinfow1(r, c);
                }

                if (!SwitchImpermeable)
                    Perc->Drc = cell_Percolation(r, c, 1.0);
           }}
        }

        // SWATRE infiltration
        if (InfilMethod == INFIL_SWATRE) {
            InfilSwatre();
        }
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        cell_SurfaceStorage(r, c);
        //calc surf storage and total watervol and WHrunoff
    }}

    if (SwitchErosion) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            cell_SplashDetachment(r, c);
            // if (SwitchSlopeStability)
            //     cell_SlopeStability(r, c);
        }}
    }
    //MoistureContent();
    // double soiltot2 = SoilWaterMass();
    // if (InfilMethod != INFIL_SOAP)
    //     SoilMoistDiff = soiltot2 - soiltot1;

}
//---------------------------------------------------------------------------


