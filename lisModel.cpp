
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

        DEBUG("setupHydrographData()");
        setupHydrographData();

        //create error file
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

        ThreadPool->StartReportThread(this);

        DEBUG("Running...");

        bool first = true;

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



            DEBUG(QString("Running timestep %1").arg((this->time - this->BeginTime)/_dt));


            ////START CALCULATIONS
            ///Function wrappers for multithreading calls

                //CELL PROCESSES
                CellProcessWrapper();

                //DYNAMIC PROCESSES
                DynamicProcessWrapper();

            ////END CALCULATIONS


            ThreadPool->WaitForReportThread();
            mutex.lock();
            if(stopRequested) DEBUG("User interrupt...");
            if(stopRequested) break;
            mutex.unlock();
            mutex.lock();
            if (waitRequested) DEBUG("User pause...");
            if (waitRequested) condition.wait(&mutex);
            mutex.unlock();
            // check if user wants to quit or pause

            Totals();            // calculate all totals and cumulative values

            MassBalance();       // check water and sed mass balance

            QFile efout(resultDir+errorFileName);
            efout.open(QIODevice::Append | QIODevice::Text);
            QTextStream eout(&efout);
            eout << "timestep " << runstep << " mass balance : " << MB << " average DT : " << UF_DTAverage << "\n";
            efout.flush();
            efout.close();

            OutputUI();          // fill the "op" structure for screen output

            DEBUG("Report to files");

            reportAll();          // report all maps and timeseries

            DEBUG("Report to interface");


            if (!noInterface)
                emit show();
            // send the op structure with data to function worldShow in LisUIModel.cpp
        }

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

void TWorld::CellProcessWrapper()
{

    //these proceses can not be done multithreaded yet

    RainfallMap();         // get rainfall from table or mpas
    SnowmeltMap();         // get snowmelt

    //run the created function on seperate threads
    ThreadPool->RunCellCompute(fcompute);
    //now wait till all threads are done!
    ThreadPool->WaitForAll();
}

void TWorld::CellProcesses(int thread)
{

    //based on the value of thread, a part of the map is used!

    GridCell(thread);            // set channel widths, flowwidths road widths etc

    Interception(thread);        // vegetation interception
    InterceptionLitter(thread);  // litter interception
    InterceptionHouses(thread);  // urban interception

    addRainfallWH(thread);       // adds rainfall to runoff water height or flood water height

    Infiltration(thread);        // infil of overland flow water, decrease WH

    SoilWater(thread);           // simple soil water balance, percolation from lower boundary
    SurfaceStorage(thread);      // surface storage and flow width, split WH in WHrunoff and WHstore

    SplashDetachment(thread);    // splash detachment


}

void TWorld::DynamicProcessWrapper()
{


    //set input for unified flow model
    //put in the multithreaded cell processes?
    UF_SetInput();

    ////Functions below are not yet multithreaded
    //does multithreading itself
    SlopeStability();      // slope stability calculations

    //slope failure must be before flow calculations. Slope failure calculations can be used in entrainment/deposition
    SlopeFailure();        // slope failure, transfers solids and liquids to unified flow equations

    Seismic();

    //does multithreading itself
    UnifiedFlow();      	//Unified flow method

    //set output of unified flow equations for lisem
    UF_SetOutput();

}
