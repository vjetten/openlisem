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
//#include "model.h"
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

   // label_runfilename->setText(QFileInfo(op.runfilename).fileName());
    /* TODO if run from commandline this name must exist */

    savefile(QString(op.LisemDir+"openlisemtmp.run"));
    // save the current settings as a runfile that is read by the model
    // in savefile(string) the runfile is updated with all user options and map names

    tabWidget->setCurrentIndex(2);
    tabWidget_out->setCurrentIndex(0);
    //switch to output screen

    checkBoxComboMaps->setChecked(true);
    checkBoxComboMaps2->setChecked(false);

    checkMapFlowBarriers->setChecked(false);
    //transparencyBarrier->setEnabled(checkFlowBarriers->isChecked());
    //checkMapFlowBarriers->setEnabled(checkFlowBarriers->isChecked());
    transparencyBarrier->setEnabled(checksatImage->isChecked());
    checkMapFlowBarriers->setEnabled(checksatImage->isChecked());

    checkMapChannels->setChecked(false);
    transparencyChannel->setEnabled(checkIncludeChannel->isChecked());
    checkMapChannels->setEnabled(checkIncludeChannel->isChecked());

    checkMapBuildings->setChecked(false);
    transparencyHouse->setEnabled(checkHouses->isChecked());
    checkMapBuildings->setEnabled(checkHouses->isChecked());

    checkMapRoads->setChecked(false);
    transparencyRoad->setEnabled(checkRoadsystem->isChecked());
    checkMapRoads->setEnabled(checkRoadsystem->isChecked());

    startplot = true;
    stopplot = false;

    // initialize output graphs and maps
    initPlot();

    initMapPlot();

    initOutputData();

    initOP();
    // reset op structure

    showOutputData();

    //=======================================================================================//

    W = new TWorld();
    // make the model world !!!

    connect(W, SIGNAL(show(void)),this, SLOT(worldShow(void)),Qt::BlockingQueuedConnection);//Qt::QueuedConnection);
    connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
    connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
    // connect emitted signals from the model thread to the interface routines that handle them

    WhasStopped = false;
    W->stopRequested = false;
    // stoprequested is used to stop the thread with the interface
    W->waitRequested = false;
    // waitrequested is used to pause the thread with the interface, only on windows machines!
    W->noInterface = false;
    W->noOutput = false;
    W->batchmode = false;
    // run without Qt interface on openlisemtmp.run only

    //setFloodOP();

    W->start();
    // start the model thread, executes W->run()

    E_runFileList->setEnabled(false);
    label_1->setEnabled(false);
    toolButton_fileOpen->setEnabled(false);
    toolButton_deleteRun->setEnabled(false);

    op.timeStartRun =  QDateTime().currentDateTime().toString("yyMMdd-hhmm");//"hh.mm-yy.MM.dd");

    screenShotDir = CheckDir(E_ResultDir->text(), true) + "screens" + op.timeStartRun + "/";
    QDir(E_ResultDir->text()).mkpath(screenShotDir);
    qDebug() << screenShotDir;
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
    if(W) {
        W->stopRequested = true;
        WhasStopped = true;
    }
}
//---------------------------------------------------------------------------
void lisemqt::worldShow()
{
    progressBar->setMaximum(op.maxstep);
    progressBar->setValue(op.runstep);

    startPlots(); // called once using bool startplot

    showOutputData(); // show output data for all and point x


    GetPlotData(); // get the plot data from the output structure

    SetTextHydrographs(); // show text hydrograph data

    showPlot(); // show main plot for point X

    showBaseMap(); // show shaded relief base map, only once, set startplot to false

    showChannelMap(); // show channel map

    showRoadMap(); // show road map

    showHouseMap(); // show building structures map

    showFlowBarriersMap();

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
        //https://stackoverflow.com/questions/31442006/properly-delete-qthread
        W->quit();
        W->wait();
       // W->deleteLater();
        delete W;
        W = nullptr;
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

    E_runFileList->setEnabled(true);
    label_1->setEnabled(true);
    toolButton_fileOpen->setEnabled(true);
    toolButton_deleteRun->setEnabled(true);

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
/*
struct output{
    int runstep;
    int printstep;
    int maxstep;

    QList<int> OutletIndices;
    QList<int> OutletLocationX;
    QList<int> OutletLocationY;
    QList<QList<double>*> OutletQ;
    QList<QList<double>*> OutletQs;
    QList<QList<double>*> OutletC;
    QList<QList<double>*> OutletChannelWH;
    QList<double> OutletQpeak;
    QList<double> OutletQpeaktime;
    QList<double> OutletQtot;
    QList<double> OutletQstot;
    double timestep;

    int
    F_solution,
    F_scheme,
    F_SSMethod,
    F_BLMethod,
    F_SigmaDiffusion,
    F_fluxLimiter,
    double
    F_MinTimestepFlood,
    F_courant,
    F_courant_diffusive;

    double CatchmentArea, dx, t,time, maxtime, EndTime, BeginTime;

    double
    // water
    MB, Qtot,  Qtile, Qtiletot, RunoffFraction, RainpeakTime,
    Qtotmm,  IntercTotmm, IntercHouseTotmm, WaterVolTotmm,InfilTotmm,
    RainTotmm, SurfStormm, InfilKWTotmm, Pmm, BaseFlowtotmm,LitterStorageTotmm,WaterVolTotchannelmm,
    floodBoundaryTot, floodBoundarySedTot,
    // channel
    ChannelVolTotmm, ChannelSedTot, ChannelDepTot, ChannelDetTot, ChannelWH,
    // flood
    FloodTotMax, FloodAreaMax, WHflood, Qflood, volFloodmm,
    FloodDetTot, FloodDepTot, FloodSedTot,
    // sediment
    MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedTot;


    cTMap *baseMap;
    cTMap *baseMapDEM;
    cTMap *channelMap;
    cTMap *roadMap;
    cTMap *houseMap;
    cTMap *flowbarriersMap;
    cTRGBMap *Image;

    QList<double> graindiameters;

    //combox selection of drawn map
    QList<int> ComboLists;
    QList<cTMap *> ComboMaps;
    QList<cTMap *> ComboMapsSafe;
    QList<QList<double>> ComboColorMap;
    QList<QList<QString>> ComboColors;
    QList<bool> ComboLogaritmic;
    QList<bool> ComboSymColor;
    QStringList ComboMapNames;
    QStringList ComboUnits;
    QList<double> ComboScaling;
    QList<double> userMinV;
    QList<double> userMaxV;
    QList<double> comboStep;

    bool comboboxset;
    bool has_image;

    QString runfilename;
    QString LisemDir;
    QString format;

    bool doBatchmode;
};

*/
    op.OutletIndices.clear();
    op.OutletLocationX.clear();
    op.OutletLocationY.clear();
    op.OutletQ.clear();
    op.OutletQs.clear();
    op.OutletC.clear();
    op.OutletQpeak.clear();
    op.OutletQpeaktime.clear();
    op.OutletChannelWH.clear();
    op.OutletQtot.clear();
    op.OutletQstot.clear();
    op.has_image = false;
    op.Image = nullptr;

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
    op.baseMap = nullptr;
    op.baseMapDEM = nullptr;
    op.channelMap = nullptr;
    op.roadMap = nullptr;
    op.houseMap = nullptr;
    op.flowbarriersMap = nullptr;
    op.Image = nullptr;

    op.runstep = 0;
    op.printstep = 0;
    op.maxstep = 0;
    op.CatchmentArea = 1.0;
    op.dx = 0;
    op.t = 0;
    op.time = 0;
    op.maxtime = 0;
    op.EndTime = 0;
    op.BeginTime = 0;
    op.MB = 0;
    op.Qtot = 0;
    op.Qtile = 0;
    op.Qtiletot = 0;
    op.RainpeakTime = 0;
    op.RunoffFraction = 0;
    op.FloodTotMax = 0;
    op.FloodAreaMax = 0;
    op.BaseFlowtotmm = 0;
    op.LitterStorageTotmm = 0;
    op.WaterVolTotchannelmm = 0;
    op.Pmm = 0;
    op.Qtotmm = 0;
    op.IntercTotmm = 0;
    op.IntercHouseTotmm = 0;
    op.WaterVolTotmm = 0;
    op.StormDrainTotmm = 0;
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
    op.ChannelVolTotmm = 0;
    op.ChannelSedTot = 0;
    op.ChannelDepTot = 0;
    op.ChannelDetTot = 0;
    op.ChannelWH = 0;
    op.FloodSedTot = 0;
    op.FloodDepTot = 0;
    op.FloodDetTot = 0;
    op.volFloodmm = 0;

//    op.F_solution = 1;
//    op.F_scheme = 1;
//    op.F_fluxLimiter = 1;
//    op.F_MinTimestepFlood = 0.001;
//    op.F_courant = 0.2;
//    op.F_courant_diffusive = 0.2;
//    op.F_Maxiter = 200;
//    op.F_SigmaDiffusion = 1;

    op.format = "PCRaster";

}
//---------------------------------------------------------------------------
void lisemqt::setFloodOP()
{
//    op.F_solution = E_floodSolution->text().toInt();
//    op.F_fluxLimiter = E_FloodFluxLimiter->text().toInt();
//    op.F_scheme = 1;//E_FloodScheme->text().toInt();
//    op.F_SSMethod = E_SSMethod->text().toInt();
//    op.F_SigmaDiffusion = E_SigmaDiffusion->text().toDouble();
//    op.F_BLMethod = E_BLMethod->text().toInt();
//    //op.F_MinTimestepFlood = E_TimestepMinFlood->text().toDouble();
//    op.F_courant = E_courantFactor->text().toDouble();
//    op.F_courant_diffusive = E_courantFactorSed->text().toDouble();
//    op.F_Maxiter = E_FloodMaxIter->text().toInt();
}
