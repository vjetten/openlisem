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
 \file LisUIModel.cpp
 \brief manage the model runs

  This file maintains the link between the interface and the LISEM model
  with data and graph displayed on the screen
  it creates the model "world" W and activates and kills the model thread
 */


#include "lisemqt.h"
#include "model.h"
#include "global.h"

//---------------------------------------------------------------------------
/** Run the model:
Save the current inyerface as a temporari run file, read by the model
Make the model world and run it
*/
void lisemqt::runmodel()
{

    if(W)
    {
        if (W->waitRequested)
            pausemodel();
        return;
    }


    //NOTE op.runfilename is set in function openRunFile()
    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this,"openLISEM",QString("Load a runfile first!"));
        return;
    }

    label_runfilename->setText(QFileInfo(op.runfilename).fileName());
    /* TODO if run from commandline this name must exist */

    savefile(QString(op.LisemDir+"openlisemtmp.run"));
    // save the current settings as a runfile that is read by the model
    // in savefile(string) the runfile is updated with all user options and map names

    tabWidget->setCurrentIndex(2);
    tabWidget_out->setCurrentIndex(0);
    //switch to output screen

    startplot = true;
    stopplot = false;

    // initialize output graphs and maps
    initPlot();

    initMapPlot();

    initOutputData();


    initOP();
    // reset op structure

    showOutputData();
    // wipe the result screen after op reset

    //=======================================================================================//
    W = new TWorld();
    // make the model world !!!

    connect(W, SIGNAL(show(void)),this, SLOT(worldShow(void)),Qt::BlockingQueuedConnection);//Qt::QueuedConnection);
    connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
    connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
    // connect emitted signals from the model thread to the interface routines that handle them

    W->stopRequested = false;
    // stoprequested is used to stop the thread with the interface
    W->waitRequested = false;
    // waitrequested is used to pause the thread with the interface, only on windows machines!
    W->noInterface = false;
    W->noOutput = false;
    W->batchmode = false;
    // run without Qt interface on openlisemtmp.run only

    first3d = true;

    W->start();
    // start the model thread, executes W->run()
    //=======================================================================================//

}
//---------------------------------------------------------------------------
void lisemqt::pausemodel()
{
    if(W)
    {
        W->waitRequested = !W->waitRequested;
        if (!W->waitRequested)
        {
            runAct->setChecked(true);
            stopAct->setChecked(false);
            pauseAct->setChecked(false);
            label_debug->setText("User continue...");
            W->condition.wakeAll();
        }
        else
        {
            stopAct->setChecked(false);
            runAct->setChecked(false);
            pauseAct->setChecked(true);
        }
    }
    else
    {
        stopAct->setChecked(false);
        runAct->setChecked(false);
        pauseAct->setChecked(false);
    }
}
//---------------------------------------------------------------------------
void lisemqt::stopmodel()
{
    creator->DestroyWorldFromLisem();

    if(W)
        W->stopRequested = true;
    if(W->waitRequested)
    {
        W->condition.wakeAll();
    }
}
//---------------------------------------------------------------------------
void lisemqt::worldShow()
{
    progressBar->setMaximum(op.maxstep);
    progressBar->setValue(op.runstep);

    startPlots(); // called once using bool startplot

    showOutputData(); // show output data for all and point x

    SetTextHydrographs(); // show text hydrograph data

    GetPlotData(); // get the plot data from the output structure

    showPlot(); // show main plot for point X

  //  showSmallPlot(); // show small plot next map for point X

    // draw maps
    showBaseMap(); // show shaded relief base map, only once, set startplot to false

    showChannelMap(); // show channel map

    showRoadMap(); // show road map

    showHouseMap(); // show building structures map

    startplot = false;

    showMap(); // show map


    if(first3d)
    {
        creator->CreateWorldFromLisem();
        first3d = false;
    }else
    {
        creator->UpdateWorldFromLisem();
    }


    if (doShootScreens)
       shootScreen();
}
//---------------------------------------------------------------------------
void lisemqt::worldDone(const QString &results)
{
    label_debug->setText(results);
    if (results.contains("ERROR"))
        QMessageBox::critical(this,QString("openLISEM"), results, QMessageBox::Ok );

    tabWidget->setCurrentIndex(2);
    shootScreen();

    // arrive here after model emits done signal
    if (W)
    {
        delete W;
        W = NULL;
    }
    //free the world instance


    stopplot = true;

    // free the map plot discharge bdata
    QFile::remove(QString(op.LisemDir+"openlisemtmp.run"));

    // delete the temp run file
    qDebug() << QString(op.LisemDir+"openlisemtmp.run")<< "deleted";

    stopAct->setChecked(false);
    runAct->setChecked(false);
    pauseAct->setChecked(false);

    if (doBatchmode)
        close();
}
//---------------------------------------------------------------------------
// this function is linked to the debug signal emitted from the model world
void lisemqt::worldDebug(const QString &results)
{
    label_debug->setText(results);
    // show messages from the World model on the screen
}
//---------------------------------------------------------------------------
void lisemqt::initOP()
{

    op.ComboMaps.clear();
    op.ComboMapsSafe.clear();
    op.ComboColorMap.clear();
    op.ComboColors.clear();
    op.ComboLogaritmic.clear();
    op.ComboSymColor.clear();
    op.ComboMapNames.clear();
    op.ComboUnits.clear();
    op.ComboScaling.clear();

    op.comboboxset = false;

    op.graindiameters.clear();
    op.baseMap = NULL;
    op.baseMapDEM = NULL;
    op.channelMap = NULL;
    op.roadMap = NULL;
    op.houseMap = NULL;
    op.vegcover = NULL;
    op.vegheight = NULL;
    op.randomroughness = NULL;
    op.gl_flow_height = NULL;
    op.gl_flow_u = NULL;
    op.gl_flow_v = NULL;
    op.gl_flow_c = NULL;

    op.runstep = 0;
    op.printstep = 0;
    op.maxstep = 0;
    op.CatchmentArea = 0;
    op.dx = 0;
    op.t = 0;
    op.time = 0;
    op.maxtime = 0;
    op.EndTime = 0;
    op.BeginTime = 0;

    op.MB = 0;
    op.Qtot = 0;
    op.RainpeakTime = 0;
    op.RunoffFraction = 0;
    op.FloodTotMax = 0;
    op.FloodAreaMax = 0;

    op.BaseFlowtot = 0;
    op.LitterStorageTot = 0;
    op.WaterVolTotchannelmm = 0;

    op.Pmm = 0;
    op.Qtotmm = 0;
    op.IntercTotmm = 0;
    op.IntercHouseTotmm = 0;
    op.WaterVolTotmm = 0;
    op.InfilTotmm = 0;
    op.RainTotmm = 0;
    op.SurfStormm = 0;
    op.InfilKWTotmm = 0;

    op.WHflood = 0;
    op.Qflood = 0;

    op.MBs = 0;
    op.DetTot = 0;
    op.DetTotSplash = 0;
    op.DetTotFlow = 0;
    op.DepTot = 0;
    op.SoilLossTot = 0;
    op.SedTot = 0;

    op.ChannelVolTot = 0;
    op.ChannelSedTot = 0;
    op.ChannelDepTot = 0;
    op.ChannelDetTot = 0;
    op.ChannelWH = 0;
    op.FloodSed = 0;
    op.FloodDepTot = 0;
    op.FloodDetTot = 0;

    op.volFloodmm = 0;

    op.F_solution = 1;
    op.F_scheme = 1;
    op.F_fluxLimiter = 1;
    op.F_replaceV = 1;
    op.F_maxVelocity = 99.0,
    op.F_extremeHeight = 5.0,
    op.F_extremeDiff = 2.0;
    op.F_courant = 0.2;
    op.F_courant_diffusive = 0.2;
    op.F_Maxiter = 200;
    op.F_SigmaDiffusion = 1;

}
