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

    setFloodOP(true);

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
    if(W)
        W->stopRequested = true;
}
//---------------------------------------------------------------------------
void lisemqt::worldShow()
{
    progressBar->setMaximum(op.maxstep);
    progressBar->setValue(op.runstep);

    startPlots(); // called once using bool startplot

    showOutputData(); // show output data for all and point x

    showPlot(); // show main plot for point X

    showSmallPlot(); // show small plot next map for point X

    // draw maps
    showBaseMap(); // show shaded relief base map, only once, set startplot to false

    showChannelMap(); // show channel map

    showRoadMap(); // show road map

    showHouseMap(); // show building structures map

    startplot = false;

    showMap(); // show map

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

    killPlot();
    // clear() plot data

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
    op.DrawMap = NULL;
    op.DrawMap1 = NULL;
    op.DrawMap2 = NULL;
    op.DrawMap3 = NULL;
    op.DrawMap4 = NULL;
    op.DrawMap5 = NULL;
    op.DrawMap6 = NULL;
    op.DrawMap7 = NULL;
    op.baseMap = NULL;
    op.baseMapDEM = NULL;
    op.channelMap = NULL;
    op.roadMap = NULL;
    op.houseMap = NULL;

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
    op.Q = 0;
    op.Qtot = 0;
    op.Qpeak = 0;
    op.RainpeakTime = 0;
    op.QpeakTime = 0;
    op.RunoffFraction = 0;
    op.FloodTotMax = 0;
    op.FloodAreaMax = 0;

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
    op.QPlot = 0;
    op.QtotPlot = 0;
    op.SoilLossTotPlot = 0;
    op.QpeakPlot = 0;

    op.MBs = 0;
    op.Qs = 0;
    op.C = 0;
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

    op.BufferVolTot = 0;
    op.BufferSedTot = 0;

    op.volFloodmm = 0;

    op.F_solution = 1;
    op.F_scheme = 1;
    op.F_fluxLimiter = 1;
    op.F_replaceV = 1;
    op.F_maxVelocity = 99.0,
    op.F_extremeHeight = 5.0,
    op.F_extremeDiff = 2.0;
    op.F_courant = 0.2;
    op.F_Maxiter = 200;

}
//---------------------------------------------------------------------------
void lisemqt::setFloodOP(bool)
{
    op.F_solution = E_floodSolution->text().toInt();
    op.F_fluxLimiter = E_FloodFluxLimiter->text().toInt();
    op.F_scheme = E_FloodScheme->text().toInt();
    op.F_replaceV = (E_FloodReplaceVcheck->isChecked() ? 1:0);
    op.F_maxVelocity = E_FloodMaxVelocity->text().toDouble(),
    op.F_extremeHeight = E_FloodExtremeHeight->text().toDouble(),
    op.F_extremeDiff = E_FloodExtremeDiff->text().toDouble();
    op.F_courant = E_courantFactor->text().toDouble();
    op.F_Maxiter = E_FloodMaxIter->text().toInt();


//    if (p1.compare("Flooding courant factor")==0)        namelist[j].value = E_courantFactor->text();
//    //  if (p1.compare("Flooding SWOF csf factor")==0)   namelist[j].value = E_cflFactor->text();
//    if (p1.compare("Flooding SWOF scheme")==0)           namelist[j].value = E_FloodScheme->text();
//    if (p1.compare("Flooding SWOF flux limiter")==0)     namelist[j].value = ;
//    if (p1.compare("Flooding SWOF Reconstruction")==0)   namelist[j].value = E_FloodReconstruction->text();
//    if (p1.compare("Include levees")==0)                 namelist[j].value.setNum((int)checkLevees->isChecked());
//    if (p1.compare("Minimum reported flood height")==0)  namelist[j].value = E_floodMinHeight->text();
//    if (p1.compare("Flooding mixing coefficient")==0)    namelist[j].value = E_mixingFactor->text();
//    if (p1.compare("Flooding runoff partitioning")==0)   namelist[j].value = E_runoffPartitioning->text();
//    if (p1.compare("Flood initial level map")==0)        namelist[j].value.setNum((int)checkFloodInitial->isChecked());
//    if (p1.compare("Flood limit max velocity")==0)       namelist[j].value = E_FloodReplaceV->text();
//    if (p1.compare("Flood max velocity threshold")==0)   namelist[j].value = E_FloodMaxVelocity->text();
//    if (p1.compare("Flood extreme value height")==0)     namelist[j].value = E_FloodExtremeHeight->text();
//    if (p1.compare("Flood extreme value difference")==0) namelist[j].value = E_FloodExtremeDiff->text();
}
