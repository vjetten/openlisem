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
TWorld::TWorld(QObject *parent) : QObject(parent)
{
}
//---------------------------------------------------------------------------
TWorld::~TWorld()
{
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
        if (!efout.open(QIODevice::WriteOnly | QIODevice::Text))
            return;
        QTextStream eout(&efout);
        eout << "#mass balance error (%)\n";

        efout.flush();
        efout.close();
    }


 //   if (doError) {
        QFile efout(resultDir+errorFileName);
        if (!efout.open(QIODevice::Append | QIODevice::Text))
            return;
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
    QTextStream consoleout(stdout); // for info with -ni batch mode
    if (noInterface) {
        consoleout << "\nrunning OpenLISEM with:" << op.runfilename << "\n\n";
        consoleout.flush();
    }

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

        // impose cmdline params
        if (op.calhydro.size() > 0 || op.calflow.size() > 0 || op.caleros.size() > 0) {
            QStringList calstrings;
            calstrings << "Smax calibration"
                       << "RR calibration"
                       << "Ksat calibration"
                       << "Ksat2 calibration"
                       << "Ksat3 calibration"
                       << "Theta calibration"
                       << "Psi calibration"

                       << "N calibration"
                       << "Channel N calibration"
                       << "Channel Ksat calibration"
                       << "Boundary water level calibration"
                       << "Culvert size calibration"

                       << "Aggregate stability calibration"
                       << "Cohesion calibration"
                       << "Grain Size calibration D50"
                       << "Grain Size calibration D90"
                       << "Cohesion Channel calibration";

            for (int j = 0; j < nrrunnamelist; j++) {
                QString p1 = runnamelist[j].name;

                if (op.calhydro.size() > 0) {
                    for (int i = 0; i < 7; i++) {
                        if (p1.compare(calstrings[i])==0) {
                            runnamelist[j].value = op.calhydro[i];
                            qDebug() << p1 << calstrings[i] << op.calhydro[i];
                        }
                    }
                }
                if (op.calflow.size() > 0) {
                    for (int i = 0; i < 5; i++) {
                        if (p1.compare(calstrings[i+7])==0) {
                            runnamelist[j].value = op.calflow[i];
                            qDebug() << p1 << calstrings[i+7] << op.calflow[i];
                        }
                    }
                }
                if (op.caleros.size() > 0) {
                    for (int i = 0; i < 5; i++) {
                        if (p1.compare(calstrings[i+12])==0) {
                            runnamelist[j].value = op.caleros[i];
                            qDebug() << p1 << calstrings[i+12] << op.caleros[i];
                        }
                    }
                }
            }
        }

        QString S = resultDir + QFileInfo(op.runfilename).fileName();
        QFile::copy(op.runfilename, S);

        // QSaveFile file(resultDir + op.explanation + ".txt");
        // if (!file.open(QIODevice::WriteOnly))
        //     return false;
        // return file.commit();

        //time vraiables in sec
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
            BeginTime = ((btd)*1440+btm)*60; //for running in sec
            EndTime = ((etd)*1440+etm)*60;   //in sec
            op.BeginTime = BeginTime/60;// for graph drawing in min
            op.EndTime = EndTime/60;
        }
        if (EndTime < BeginTime + 60) {
            ErrorString = "End time must be > Begin time + 1 minute.";
            throw 1;
        }

        //get all maps
        DEBUG("Get Input Maps");
        GetInputData();
        DEBUG("Intialize Database");
        IntializeData();

        // MC - no_ui probalbly this can be skipped for noInterface??
        /// TODO
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
        //for output totals per landunit

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
            if (SwitchIncludeET && SwitchDailyET) {
                // increase timestep between rainfall
                if (longdt > _dt) {
                    bool do_longdt = true;
                    FOR_ROW_COL_MV_L {
                        if (/*WH->Drc > 1e-3 &&*/ V->Drc > 5e-5 ) {
                            do_longdt = false;
                            break;
                        }
                    }}
                    if (do_longdt)
                        _dt = longdt;
                    else
                        _dt =_dt_user;
                }
            }
            savemaptodisk = false;
            // printstep determines report frequency in #define report(...)
            if (runstep > 0 && runstep % printinterval == 0) {
                savemaptodisk = true;
                printstep++;
            }
            if (SwitchReportMapsEnd)
                savemaptodisk = false;

            runstep++;

            if(stopRequested) {
                mutex.lock();
                DEBUG("User interrupt... finishing time step");
                time = EndTime;
                mutex.unlock();
            }

            if (waitRequested) {
                mutex.lock();
                DEBUG("User pause...");
                mu_condition.wait(&mutex);
                mutex.unlock();
            }
            // check if user wants to quit or pause

            GetInputTimeseries(); // get rainfall, ET, snowmelt, discharge

            GetETparameters();

            InfilDynamicCrusting(); // if crusting recalc Ksateff and Poreff becuase of crusting effect

            HydrologyProcesses();  // hydrological processes in one loop, incl splash

            ToTiledrain();  // fraction going into tiledrain directly from surface

            OverlandFlow(); // overland flow 1D (non threaded), 2Ddyn (threaded), if 2Ddyn then also SWOFsediment!

            // if (SwitchIncludeChannel) {
            //     ChannelRainandInfil();  // subtract infil, retention,  add rainfall
            //     ChannelBaseflow();              // add stationary and GW baseflow if selected
            // }

            // these are all non-threaded
            ChannelFlowandErosion();    // do ordered LDD solutions channel, tiles, drains, non threaded

            TileFlow();          // tile drain flow kin wave

            TotalsHydro();       // calculate all totals and cumulative values
            TotalsFlow();
            TotalsSediment();

            MassBalance();       // check water and sed mass balance

            reportToUI();        // fill the "op" structure for screen and file output and calc some COMBO output maps

            reportToFile();      // report hydrograhs, totals, maps etc to files
            // reporting to file is done in nthe same thread, mdoes not need a mutex lock

            // because showing is done outside the Thread in the GUI, a mutex.lock() is needed
            // mu_condition gives a wakeAll() signal at the end of the display in showWorld()
            if (!noInterface) {
                emit show(); // send the 'op' structure with data to function worldShow in LisUIModel.cpp
                mutex.lock();
                mu_condition.wait(&mutex);   // Wait for GUI to finish drawing
                mutex.unlock();
            }

            //saveMBerror2file(false); //saveMBerror

            // show progress in console without GUI
            if (op.doBatchmode && noInterface) {
                int x = 0;
                x = std::round(op.t/op.maxtime * 100) ;
                consoleout << "\rprogress: " << QString("step %1        %2 %   end time %3").arg(runstep).arg(x, -3).arg(op.maxtime) << "        ";
                consoleout.flush();

                // THIS SHOULD ALSO WORK IN LINUX ???
            }
        } // TIME LOOP

        if(SwitchReportMapsEnd) {
            ReportMaps();
            ReportMapSeries();
            if(SwitchDumphead) {
                for (int i = 0; i < SwatreSoilModel->pixel[0].profile->zone->nrNodes; i++) {

                    QString dig = QString("%1").arg(i+1, 3, 10, QLatin1Char('0'));
                    QString hname = QString("head0000.") + dig;
                    QString tname = QString("theta000.") + dig;

                    #pragma omp parallel for num_threads(userCores)
                    FOR_ROW_COL_MV_L {
                        if (ProfileID->Drc <= 0 || fractionImperm->Drc > 0.999) {
                            tma->Drc = 0;
                            tmb->Drc = 0;
                        } else {
                            tma->Drc =  qMin(0.0, SwatreSoilModel->pixel[i_].h[i]);
                            tmb->Drc = FindValue(tma->Drc, SwatreSoilModel->pixel[i_].profile->horizon[i], H_COL, THETA_COL);
                        }
                    }}
                    report(*tma, hname);
                    report(*tmb, tname);
                }
            }        }

        if (!noInterface) {
            // wrap up and close the thread
            emit done("Finished");
        }

        if (op.doBatchmode) {
            if (!noInterface) {
                mutex.lock();
                emit ScreenShot();
                mu_condition.wait(&mutex);   // Wait for GUI to finish drawing
                mutex.unlock();
            }
            // delete all maps
            qDeleteAll(maplistCTMap.begin(),maplistCTMap.end());
            maplistCTMap.clear();

            // //delete swatre 3D soil layer structure if exists
            if (initSwatreStructure)
                FreeSwatreInfo();

            if (noInterface) {
                consoleout << "\n\n Finished after "<< op.maxtime << "minutes";
                consoleout.flush();
                //QCoreApplication::quit();
                // no longer used because app.exec() is not called, just let it exit
            } else {
                QApplication::quit();
            }
            // close the world model
        }
    }
    catch(...)  // if an error occurred
    {
        if (!noInterface) {
            emit done("ERROR STOP: "+ErrorString);
        }
        if (op.doBatchmode) {
            if (noInterface) {
                consoleout << "ERROR STOP "<< ErrorString;
                consoleout.flush();
                #ifdef Q_OS_WIN
                system("pause"); // waits for a key press
                #endif
               // QCoreApplication::quit();
                // no longer used because app.exec() is not called, just let it exit
            } else {
                QApplication::quit();
            }
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
        hmxWH->Drc = hmx->Drc+WH->Drc;

        // incoming wave at boundary
        if (SwitchWaveUser) {
            WHboundRain->Drc += RainNet->Drc;
            if (WHboundarea->Drc > 0) {
                // WHbound is the forced water level in area with value '1', plus cum rainfall
                WH->Drc = WHbound->Drc + WHboundRain->Drc;
                // assume there is no hmx in WHboundarea
            }
        }
    }}

    if (SwitchInfiltration) {
        // non SWATRE infiltration, redistribution and percolation
        if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_SOAP) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                cell_InfilMethods(r, c);

                if (SwitchThreeLayer) {
                    cell_RedistributionUnsat(r, c);
                    cell_Redistribution3(r, c);
                } else {
                    if (SwitchTwoLayer) {
                        cell_RedistributionUnsat(r, c);
                        cell_Redistribution2(r, c);
                        cell_Tiledrain2(r,c);
                        //cell_Channelinfow2(r, c);
                    } else {
                        cell_Redistribution1(r, c);
                        cell_Tiledrain1(r,c);
                        //cell_Channelinfow1(r, c);
                    }
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
        cell_SplashDetachment();
            // if (SwitchSlopeStability)
            //     cell_SlopeStability(r, c);
    }
    //MoistureContent();
    // double soiltot2 = SoilWaterMass();
    // if (InfilMethod != INFIL_SOAP)
    //     SoilMoistDiff = soiltot2 - soiltot1;

}
//---------------------------------------------------------------------------


