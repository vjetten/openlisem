
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
// the actual model with the main loop
void TWorld::DoModel()
{
    if (!noInterface)
        temprunname = QString(op.LisemDir+"openlisemtmp.run");
    else
        temprunname = op.runfilename;

  //  QString DT = QDateTime().currentDateTime().toString("hh.mm-yy.MM.dd");
    errorFileName = QString(resultDir + "error"+QDateTime().currentDateTime().toString("yy.MM.dd-hh.mm")+".txt");

    time_ms.start();
    // get time to calc run length


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

        BeginTime = getvaluedouble("Begin time") * 60;
        EndTime = getvaluedouble("End time") * 60;
        _dt = getvaluedouble("Timestep");
        op.BeginTime = BeginTime/60; // for graph drawing
        op.EndTime = EndTime/60;// for graph drawing
        //time vraiables in sec

        runstep = 0; // NOTE runstep is used to initialize graph!
        printstep = 1; // printstep determines report frquency

        // VJ 110630 show hydrograph for selected output point
        bool found = false;
        if (op.outputpointnr > 1)
        {
            FOR_ROW_COL_MV
            {
                if (op.outputpointnr == PointMap->Drc)
                {
                    r_plot = r;
                    c_plot = c;
                    op.outputpointdata = QString("point %1 [row %2; col %3]").arg(op.outputpointnr).arg(r).arg(c);
                    found = true;
                }
            }
            if (!found)
            {
                ErrorString = QString("Point %1 for hydrograph outputpoint 2 not found").arg(op.outputpointnr);
                throw 1;
            }
        }
        else
            op.outputpointdata = QString("Main Outlet");

        QFile efout(resultDir+errorFileName);
        efout.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream eout(&efout);
        eout << "#error tryout\n";
        eout << "2\n";
        eout << "time\n";
        eout << "error\n";
        efout.flush();
        efout.close();

        InfilEffectiveKsat();
        // calc effective ksat from all surfaces once

        DEBUG("Running...");

        for (time = BeginTime; time < EndTime; time += _dt)
        {
            if (runstep > 0 && runstep % printinterval == 0)
                printstep++;
            runstep++;

            if (noInterface && !noOutput)
            {
                int maxstep = (int)(EndTime-BeginTime)/_dt ;
                qDebug() << runstep << maxstep << time/60 ;
            }

            mutex.lock();
            if(stopRequested) DEBUG("User interrupt...");
            if(stopRequested) break;
            mutex.unlock();
            mutex.lock();
            if (waitRequested) DEBUG("User pause...");
            if (waitRequested) condition.wait(&mutex);
            mutex.unlock();
            // check if user wants to quit or pause

            GridCell();            // set channel widths, flowwidths road widths etc
            RainfallMap();         // get rainfall from table or mpas
            SnowmeltMap();         // get snowmelt

            Interception();        // vegetation interception
            InterceptionHouses();  // urban interception

            addRainfallWH();       // adds rainfall to runoff water height or flood water height

            Infiltration();        // infil of overland flow water, decrease WH
            InfiltrationFloodNew();// infil in flooded area, decrease hmx

            SoilWater();           // simple soil water balance, percolation from lower boundary
            SurfaceStorage();      // surface storage and flow width, split WH in WHrunoff and WHstore

          //  RainfallToFlood();       // converts rainfall on flat areas to flood instead of runoff
            // OBSOLETE with runoff in 2D

            CalcVelDisch();        // overland flow velocity, discharge and alpha for erosion

            SplashDetachment();    // splash detachment
            FlowDetachment();      // flow detachment

            //Pestmobilisation();  // experimental

            ToFlood();             // overland flow water added to flood (not in channel cells)
            ToChannel();           // water and sed flux going into channel in channel cells
            ToTiledrain();         // fraction going into tiledrain directly from surface

            CalcVelDisch();        // recalc overland flow velocity, discharge and alpha because of extractions

            OverlandFlowNew();     // overland flow kin wave for water and sed

            ChannelWaterHeight();  // add rainfall and runoff to channel and get channel WH from volume

            ChannelFlood();        // st venant channel 2D flooding from channel

            CalcVelDischChannel(); // alpha, V and Q from Manning

            ChannelFlow();         // channel erosion and kin wave

            TileFlow();          // tile drain flow kin wave

            //SumSedimentClasses();  //sum the induvidual sediment class layers

            Totals();            // calculate all totals and cumulative values
            MassBalance();       // check water and sed mass balance

            QFile efout(resultDir+errorFileName);
            efout.open(QIODevice::Append | QIODevice::Text);
            QTextStream eout(&efout);
            eout << " " << runstep << " " << MB << "\n";
            efout.flush();
            efout.close();

            reportAll();          // report all maps and timeseries

            OutputUI();          // fill the "op" structure for screen output
            // show after report calc is done

            if (!noInterface)
                emit show();
            // send the op structure with data to function worldShow in LisUIModel.cpp
        }

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
//---------------------------------------------------------------------------
// Make the maps to bedrawn in the interface as a copy in the op starcture
// reason is that all pointers are destroyed after the run so when lisem finishes
// the information on the output screen points to an empty pointer
// by copying the info remains available
// initialize maps for output to screen
// must be done after Initialize Data because then we know how large the map is
void TWorld::setupDisplayMaps()
{
    if (op.baseMap != 0)
    {
        delete op.baseMap;
        delete op.baseMapDEM;
        delete op.channelMap;
        delete op.roadMap;
        delete op.houseMap;
    }

    op.baseMap = new cTMap();
    op.baseMapDEM = new cTMap();
    op.channelMap = new cTMap();
    op.roadMap = new cTMap();
    op.houseMap = new cTMap();

    op.baseMap->MakeMap(LDD, 0);
    op.baseMapDEM->MakeMap(LDD, 0);
    op.channelMap->MakeMap(LDD, 0);
    op.roadMap->MakeMap(LDD, 0);
    op.houseMap->MakeMap(LDD, 0);


}
//---------------------------------------------------------------------------
void TWorld::GetComboMaps()
{
    //combo box maps
    ClearComboMaps();

    QList<double> Colormap;
    QList<QString> Colors;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.0005);
    Colormap.append(0.01);
    Colormap.append(0.05);
    Colormap.append(0.1);
    Colormap.append(0.9);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#8C8CFF");
    Colors.append("#8080FF");
    Colors.append("#4040ff");
    Colors.append("#0000FF");
    Colors.append("#00006F");
    Colors.append("#FF0000");
    Colors.append("#FF3300");
    AddComboMap(0,"Overland Discharge","l/s",Q,Colormap,Colors,true,false,1000.0);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"Channel Discharge","l/s",ChannelQ,Colormap,Colors,true,false,1000.0);
        AddComboMap(0,"Total Discharge","l/s",COMBO_QOFCH,Colormap,Colors,true,false,1000.0);
    }

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.25);
    Colormap.append(0.75);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#00FF00");
    Colors.append("#FFFF00");
    Colors.append("#FF0000");
    Colors.append("#A60000");
    AddComboMap(0,"Overland Flow Velocity","m/s",V,Colormap,Colors,false,false,1.0);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"Channel Velocity","m/s",ChannelV,Colormap,Colors,false,false,1.0);
    }

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.5);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#5477ff");
    Colors.append("#0023b1");
    Colors.append("#001462");

    AddComboMap(0,"Overland Flow Height","m",WHrunoff,Colormap,Colors,false,false,1.0);

    if(InfilMethod != INFIL_NONE)
    {
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.1);
        Colormap.append(0.4);
        Colormap.append(1.0);

        Colors.clear();
        Colors.append("#f6f666");
        Colors.append("#FFFF55");
        Colors.append("#8080FF");
        Colors.append("#0000AA");
        AddComboMap(0,"Infiltration","mm",InfilmmCum,Colormap,Colors,false,false,1.0);
    }

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.5);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#8888FF");
    Colors.append("#0000FF");
    Colors.append("#FF0000");

    AddComboMap(0,"Rainfall Cumulative","mm",RainCumFlat,Colormap,Colors,false,false,1.0);
    AddComboMap(0,"Rainfall Intensity","mm/h",Rain,Colormap,Colors,false,false,1.0);

    if(SwitchChannelFlood)
    {
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#5477ff");
        Colors.append("#0023b1");
        Colors.append("#001462");
        AddComboMap(0,"Flood Height","m",hmx,Colormap,Colors,false,false,1.0);

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.75);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#00FF00");
        Colors.append("#FFFF00");
        Colors.append("#FF0000");
        Colors.append("#A60000");
        AddComboMap(0,"Flood Velocity","m/s",UVflood,Colormap,Colors,false,false,1.0);

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#A60000");
        Colors.append("#FF0000");
        Colors.append("#FFFF00");
        Colors.append("#00FF00");
        Colors.append("#007300");

        AddComboMap(0,"Flood Start Time","min",floodTimeStart,Colormap,Colors,false,false,1.0);

    }
    if(SwitchErosion)
    {
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.3);
        Colormap.append(0.5);
        Colormap.append(0.70);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#457A60");
        Colors.append("#96B547");
        Colors.append("#FFFFFF");
        Colors.append("#FFFF00");
        Colors.append("#FF0000");

        QString unit = "kg/cell";
        double factor = 1.0;
        if(ErosionUnits == 2)
        {
            factor = 1.0/(_dx*_dx);
            unit = "kg/m2";
        }else if (ErosionUnits == 0)
        {
            factor = 10.0/(_dx*_dx);
            unit = "t/ha";
        }
        AddComboMap(1,"Total Soil Loss",unit,TotalSoillossMap,Colormap,Colors,false,true,factor);


        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(0.9);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#FAFAD2");
        Colors.append("#FFFF66");
        Colors.append("#808000");
        Colors.append("#804000");

        AddComboMap(1,"Sediment Load","kg/m2",COMBO_SS,Colormap,Colors,false,false,1.0/(_dx*_dx));

        if(SwitchUseGrainSizeDistribution)
        {
            Colormap.clear();
            Colormap.append(0.0);
            Colormap.append(0.5);
            Colormap.append(0.9);
            Colormap.append(1.0);
            Colors.clear();
            Colors.append("#FAFAD2");
            Colors.append("#FFFF66");
            Colors.append("#808000");
            Colors.append("#804000");

            FOR_GRAIN_CLASSES
            {
                AddComboMap(1,"Overland S.L. Grain Class " + QString::number(d),"kg/m2",Sed_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx));
            }
            if(SwitchIncludeChannel)
            {
                if(SwitchUse2Layer)
                {
                    FOR_GRAIN_CLASSES
                    {
                        AddComboMap(1,"Channel BL S.L. Grain Class " + QString::number(d),"kg/m2",RBL_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx));
                    }
                    FOR_GRAIN_CLASSES
                    {
                        AddComboMap(1,"Channel SS S.L. Grain Class " + QString::number(d),"kg/m2",RSS_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx));
                    }
                }else
                {
                    FOR_GRAIN_CLASSES
                    {
                        AddComboMap(1,"Channel S.L. Grain Class " + QString::number(d),"kg/m2",RBL_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx));
                    }
                }
            }
            if(SwitchChannelFlood)
            {
                FOR_GRAIN_CLASSES
                {
                    AddComboMap(1,"Flood BL S.L. Grain Class " + QString::number(d),"kg/m2",BL_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx));
                }
                FOR_GRAIN_CLASSES
                {
                    AddComboMap(1,"Flood SS S.L. Grain Class " + QString::number(d),"kg/m2",SS_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx));
                }
            }

        }



    }


}
//---------------------------------------------------------------------------
void TWorld::ClearComboMaps()
{

    for(int i =op.ComboMapsSafe.length() - 1; i >-1 ; i--)
    {
        delete op.ComboMapsSafe.at(i);
    }
    op.ComboMapsSafe.clear();

    op.ComboLists.clear();
    op.ComboMaps.clear();
    op.ComboMapsSafe.clear();
    op.ComboColorMap.clear();
    op.ComboColors.clear();
    op.ComboLogoritmic.clear();
    op.ComboSymColor.clear();
    op.ComboMapNames.clear();
    op.ComboUnits.clear();
    op.ComboScaling.clear();

    op.comboboxset = false;
}
//---------------------------------------------------------------------------
void TWorld::AddComboMap(int listn, QString name, QString unit,cTMap * map,QList<double> ColorMap, QList<QString> Colors, bool log,bool symcol, double scale)
{
    op.ComboLists.append(listn);
    op.ComboMaps.append(map);
    op.ComboMapsSafe.append(new cTMap());
    op.ComboMapsSafe.at(op.ComboMapsSafe.length()-1)->MakeMap(LDD,0.0);
    op.ComboColorMap.append(ColorMap);
    op.ComboColors.append(Colors);
    op.ComboLogoritmic.append(log);
    op.ComboSymColor.append(symcol);
    op.ComboMapNames.append(name);
    op.ComboUnits.append(unit);
    op.ComboScaling.append(scale);

    op.comboboxset = false;
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
