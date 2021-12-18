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

    loc = QLocale::system(); // current locale
    loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
    QLocale::setDefault(loc); // set as default
    //QLocale somel;

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

    try
    {
       // DEBUG("reading and initializing data");

        IntializeOptions(); // reset all options

        InitMapList();
        // map structure to destroy data automatically

   //     DEBUG("GetRunFile()");
        GetRunFile();
        DEBUG("Parse Runfile");
        ParseRunfileData();
        // get and parse runfile

        //BeginTime = getTimefromString(bt)*60; // in seconds!
        //EndTime = getTimefromString(et)*60;
        double btd = getvaluedouble("Begin time day");
        double btm = getvaluedouble("Begin time");
        double etd = getvaluedouble("End time day");
        double etm = getvaluedouble("End time");
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
            BeginTime = (btd*1440+btm)*60; //for eunning in sec
            EndTime = (etd*1440+etm)*60;   //in sec
            op.BeginTime = BeginTime/60;// for graph drawing in min
            op.EndTime = EndTime/60;
        }
        //VJ get time here else combomaps goes wrong for rainfall intensity

        //time vraiables in sec
//        DEBUG("Get Input Data");
        GetInputData();
        DEBUG("Intialize Input Data()");
        IntializeData();

    //    DEBUG("setupDisplayMaps()");
        setupDisplayMaps();
        // reset all display output maps for new job
        // must be done after Initialize Data because then we know how large the map is

        SwitchSnowmelt = false;
        if (SwitchRainfall)
        {
           // qDebug() << rainSatFileName << rainFileName;
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
                   //     qDebug() << ETSatFileName << ETFileName;
            DEBUG("Get ET Data Information");
            if (SwitchETSatellite)
                GetSpatialMeteoData(ETSatFileName, 1);
            else
                GetETData(ETFileName);
        }
        if (SwitchSnowmelt)
        {
            DEBUG("Get Snowmelt Data Information");
            if (SwitchSnowmeltSatellite)
                GetSpatialMeteoData(snowmeltSatFileName, 2);
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

      //  DEBUG("setupHydrographData()");
        setupHydrographData();

        bool saveMBerror = false;
        saveMBerror2file(saveMBerror, true);

        InfilEffectiveKsat();  // calc effective ksat from all surfaces once
        SetFlowBarriers();     // update the presence of flow barriers, static for now, unless breakthrough
        GridCell();            // static for now

        _dt_user = _dt;

        DEBUG("Running...");

        GetComboMaps(); // moved to outside timeloop!

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

            GetInputTimeseries(); // get rainfall, ET, snowmelt, discharge

            HydrologyProcesses();  // hydrological processes in one loop, incl splash

            OverlandFlow(); // overland flow 1D (non threaded), 2Ddyn (threaded), if 2Ddyn then also SWOFsediment!

            ChannelandTileflow();    // do ordered LDD solutions channel, tiles, drains, non threaded

            Totals();            // calculate all totals and cumulative values

            MassBalance();       // check water and sed mass balance

            OutputUI();          // fill the "op" structure for screen output and calc some output maps

            reportAll();         // report maps and files to screen and disk

            emit show(noInterface); // send the 'op' structure with data to function worldShow in LisUIModel.cpp

            saveMBerror2file(saveMBerror, false);

            if(stopRequested)
                time = EndTime;                       
        }

        if (SwitchEndRun)
            ReportMaps();

        DEBUG("Free data structure memory");
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
void TWorld::GetInputTimeseries()
{
    // get meteo data
    if (SwitchRainfallSatellite)
        GetRainfallMap();         // get rainfall from maps
    else
        GetRainfallMapfromStations();         // get rainfall from stations

    if (SwitchIncludeET) {
        if (SwitchETSatellite)
            GetETSatMap(); // get rainfall from maps
        else
            GetETMap();   // get rainfall from stations
    }

    if (SwitchSnowmelt) {
        if (SwitchSnowmeltSatellite)
            ; //TODO snowmelt satellite
        else
            GetSnowmeltMap();  // get snowmelt from stations
    }

}
//---------------------------------------------------------------------------
void TWorld::HydrologyProcesses()
{

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        if (Rainc->Drc > 0)
            cell_Interception(r,c);
        // all interception on plants, houses, litter
        // result is rainnet (and leafdrip for erosion)

        if (FloodDomain->Drc > 0) {            
            hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        } else {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
        }
        // add net to water rainfall on soil surface (in m)
        // when kin wave and flooded hmx exists else always WH

        if (RoadWidthHSDX->Drc > 0)
            WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
        // tarred roads have no interception ?

        // infiltration by SWATRE of G&A+percolation
        if (InfilMethod == INFIL_SWATRE) {
           cell_InfilSwatre(r, c);
        }
        else
        if (InfilMethod != INFIL_NONE) {
           cell_InfilMethods(r, c);

           cell_Redistribution(r, c);

           if (!SwitchImpermeable && !SwitchChannelBaseflow)
               Perc->Drc = cell_Percolation1(r, c, 1.0);
            // if baseflow this already does percollation, no need to do this twice
        }

        // deposit all sediment still in flow when infiltration causes WH to become minimum
        if(SwitchErosion) {
            if(SwitchKinematic2D == K2D_METHOD_DYN) {
                if(WH->Drc < 1e-6) {
                    DepFlood->Drc -= SSFlood->Drc;
                    SSFlood->Drc = 0;
                    SSCFlood->Drc = 0;
                    SSTCFlood->Drc = 0;
                    if (SwitchUse2Phase) {
                        DepFlood->Drc -= BLFlood->Drc;
                        BLFlood->Drc = 0;
                        BLCFlood->Drc = 0;
                        BLTCFlood->Drc = 0;
                    }
                    Conc->Drc = 0; // after dynwave conc is sum of SS and BL!
                }
            } else {
                if (FloodDomain->Drc > 0) {
                    if(hmx->Drc < 1e-6) {
                        DepFlood->Drc -= SSFlood->Drc;
                        SSFlood->Drc = 0;
                        SSCFlood->Drc = 0;
                    }

                } else {
                    if(WH->Drc < 1e-6) {
                        DEP->Drc -= Sed->Drc;
                        Sed->Drc = 0;
                        Conc->Drc = 0;
                    }
                }
            }
        }

        cell_SurfaceStorage(r, c);

        if (SwitchErosion) {
            double wh = FloodDomain->Drc == 0 ? WH->Drc : hmx->Drc;
            cell_SplashDetachment(r,c,wh);
        }
    }}

    if (SwitchIncludeET)
        doETa();
    // ETa is subtracted from canopy, soil water surfaces
    // divided over 12 hours in a day with sine curve

    avgTheta();

}
//---------------------------------------------------------------------------
// these are all non-threaded
void TWorld::ChannelandTileflow()
{
    SwitchChannelKinWave = true;// set to false for experimental swof in channel

    ChannelRainandInfil();   // subtract infil, add rainfall

    ChannelBaseflow();       // calculate baseflow

    if (SwitchChannelKinwaveDt) {
        if (_dt_user > _dtCHkin) {
            double n = _dt_user/_dtCHkin;
            _dt = _dt_user/n;
        }
    }

    for (double t = 0; t < _dt_user; t+=_dt) {
        ChannelFlow();            //channel kin wave for water and sediment
    }
    _dt=_dt_user;
    //ChannelFillDam();

    ChannelFlowDetachmentNew();  //detachment, deposition for SS and BL

    TileFlow();          // tile drain flow kin wave

    StormDrainFlow();    // storm drain flow kin wave
}

