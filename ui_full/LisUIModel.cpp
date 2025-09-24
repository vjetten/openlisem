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
 \file LisUIModel.cpp
 \brief manage the model runs

  This file maintains the link between the interface and the LISEM model
  with data and graph displayed on the screen
  it creates the model "world" W and activates and kills the model thread
 */


#include "lisemqt.h"
#include "global.h"

#define Clear2D(List) qDeleteAll(List.begin(), List.end());List.clear();

void lisemqt::ClearOP()
{
    Clear2D(op.OutletQ);
    Clear2D(op.Wavein);
    Clear2D(op.OutletQs);
    Clear2D(op.OutletC);
    Clear2D(op.OutletChannelWH);

    op.OutletIndices.clear();
    op.OutletLocationX.clear();
    op.OutletLocationY.clear();
    op.OutletQpeak.clear();
    op.OutletQpeaktime.clear();
    op.OutletQtot.clear();
    op.OutletQstot.clear();
    op.Pmm.clear();
    op.Time.clear();
    op.Qtile.clear();
    op.EndPointX.clear();
    op.EndPointY.clear();
    op.ObsPointX.clear();
    op.ObsPointY.clear();
    op.lddch_.clear();

    delete op.baseMap;
    delete op.baseMapDEM;
    delete op.channelMap;
    delete op.outletMap;
    delete op.roadMap;
    delete op.houseMap;
    delete op.hardsurfaceMap;
    delete op.Image;
}

void lisemqt::deleteWStructures()
{
    // destroy ALL maps
    qDeleteAll(W->maplistCTMap.begin(),W->maplistCTMap.end());
    W->maplistCTMap.clear();
    // destroy all network structures
    W->cr_.clear();
    W->crch_.clear();
    W->crlinkedldd_.clear();
    W->crlinkedlddch_.clear();
    W->crldd5_.clear();
    W->crlddch5_.clear();
    W->crout_.clear();
    W->dcr_.clear();
    W->crtile_.clear();

    // delete PGraph;
    // delete QGraph;
    // delete QsGraph;
    // delete CGraph;
    // delete QtileGraph;
    // delete QbGraph;

    QVector <double> zero;
    zero.clear();
    PGraph->setSamples(zero,zero);
    QGraph->setSamples(zero,zero);

    QsGraph->setSamples(zero,zero);
    CGraph->setSamples(zero,zero);

    QtileGraph->setSamples(zero,zero);
    QbGraph->setSamples(zero,zero);

    HPlot->replot();

    ClearOP(); // clear most of the op structure

    // destroy swatre structures
    if (W->initSwatreStructure) {
        W->FreeSwatreInfo(); // free horizon structures, this calls also closeswatre
    }

    // drawing riuvers on screen structures
    Xa.clear();
    Ya.clear();
    Xc.clear();
    Yc.clear();
    op.ObsPointX.clear();
    op.ObsPointY.clear();
    op.EndPointX.clear();
    op.EndPointY.clear();
}
//---------------------------------------------------------------------------
/** Run the model:
Save the current interface as a temporary run file, read by the model
Make the model world and run it
*/
void lisemqt::runmodel()
{
    //NOTE op.runfilename is set in function openRunFile()
    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this,"openLISEM",QString("Load a runfile first!"));
        return;
    }

    // connect emitted signals from the model thread to the interface routines that handle them

    // if the model has stopped and a new run is requested, clear the datastructures
    // until that time the user can look at the old results
    // we do that at the start of a new run and not at the end of a run,
    //because the user wants to still switch maps in the interface after the ruin
    if (stoprun && W) {
        deleteWStructures();
    }

    startplot = true; // user has pressed run, used only to initiatte screen stop, after that set to false!
    stoprun = false; // user has not stopped the run

    label_debug->text().clear();

    lastOptionSceen = tabWidgetOptions->currentIndex();

    showOutputDataZero();

    updateModelData(); // read temporary run file
    QFile f(QString(op.userAppDir+"openlisemtmp.run"));
    if (f.exists())
        f.remove();

    savefile(QString(op.userAppDir+"openlisemtmp.run"));
    // save the current settings as a runfile that is read by the model
    // in savefile(string) the runfile is updated with all user options and map names

    tabWidget->setCurrentIndex(2);
    tabWidget_out->setCurrentIndex(0);
    //switch to output screen

    checkBoxComboMaps->setChecked(true);
    checkBoxComboMaps2->setChecked(false);

    checkMapImage->setChecked(false);
    //transparencyImage->setEnabled(checksatImage->isChecked());
    checkMapImage->setEnabled(checksatImage->isChecked());

    checkMapChannels->setChecked(false);
    checkMapChannels->setEnabled(checkIncludeChannel->isChecked());

    if (checkInfrastructure->isChecked()) {
         checkMapBuildings->setChecked(checkHouses->isChecked());
         checkMapRoads->setChecked(checkRoadsystem->isChecked());
         checkMapHardSurface->setChecked(checkHardsurface->isChecked());
        transparencyHardSurface->setValue(200);
        transparencyRoad->setValue(200);
    }

    showInfoAct->setChecked(true);
    setOutputInfo(true); // show the cursor over the map

    // initialize output graphs
    initPlot();

    initMapPlot();

    initOP();
    // reset op structure

    showOutputData();

    //=======================================================================================//
    // create the world and the thread it runs in

    if (W) {
        W->deleteLater(); // delete after it is finished
    }

    W = new TWorld();
    connect(W, &TWorld::show, this, &lisemqt::worldShow);
    connect(W, &TWorld::done, this, &lisemqt::worldDone);
    connect(W, &TWorld::debug, this, &lisemqt::worldDebug);
    connect(W, &TWorld::timedb, this, &lisemqt::worldDebug);
    //connections to trigger messages and model stop from the interface
    // e.g. if the world emits done, the worldDone is called to stop the model

    // dealing with digit separator comma or dot
    W->loc = QLocale::system(); // current locale
    W->loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
    QLocale::setDefault(W->loc);

    // make a thread to run the world in
    worldThread = new QThread();
    W->moveToThread(worldThread);

    connect(worldThread, &QThread::started, W, &TWorld::DoModel);
    connect(W, &TWorld::done, worldThread, &QThread::quit);
    connect(worldThread, &QThread::finished, worldThread, &QThread::deleteLater); // dlete later means these are automatically deleted when the thread finishes

    W->showInfo = true;

    //WhasStopped = false;
    W->stopRequested = false;
    // stoprequested is used to stop the thread with the interface
    W->waitRequested = false;
    // waitrequested is used to pause the thread with the interface, only on windows machines!
    W->noInterface = false; // batchmode if true show nothing
    W->noOutput = false;// if false then show something on screen
    W->batchmode = false;
    // run without Qt interface on original runfile only

    op.timeStartRun = QDateTime().currentDateTime().toString("yyMMdd-hhmm");
    if (op.explanation != "empty" ) {
        op.timeStartRun = op.explanation;
        checkAddDatetime->setChecked(true);
    }


    if (checkAddDatetime->isChecked()) {
        screenShotDir = E_ResultDir->text() + QString("res"+op.timeStartRun+"/");
        QDir(screenShotDir).mkpath(QString("screens/"));
        screenShotDir = screenShotDir + QString("screens/");
    } else {
        screenShotDir = E_ResultDir->text();
        QDir(screenShotDir).mkpath(QString("screens"+op.timeStartRun+"/"));
        screenShotDir = screenShotDir + QString("screens"+op.timeStartRun+"/");
    }

    // take a screenshot of all option widgets
    tabWidget->setCurrentIndex(0);
    for (int i = 0; i < tabWidgetOptions->count(); i++) {
        tabWidgetOptions->setCurrentIndex(i);
        shootSingleScreen(1);
    }
    tabWidget->setCurrentIndex(2);
    //switch to output screen

    worldThread->start();
    // start the model thread, executes W->run()

    E_runFileList->setEnabled(false);
    checkDoErosion->setEnabled(false);
    label_1->setEnabled(false);
    toolButton_fileOpen->setEnabled(false);
    toolButton_deleteRun->setEnabled(false);

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
          //  label_debug->setText("User continue...");
            W->mu_condition.wakeOne();//wakeAll();
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
// linked to stop button in interface, the current loop is finished before
// the model thread is really stopped, this infact sets time to endtime
// after the loop is finished worldDone is called to end the thread
void lisemqt::stopmodel()
{
    if(W) {
        W->stopRequested = true;
    }
}
//---------------------------------------------------------------------------
void lisemqt::worldShow()
{
    progressBar->setMaximum(op.maxstep);
    progressBar->setValue(op.runstep);

    startPlots(); // called once using bool startplot

    showOutputData(); // show output data of totals as minimumfeedback

    if (!W->noOutput) {
        showPlot(); // show main plot for point X

        showBaseMap(); // show shaded relief base map, only once, set startplot to false

        getOutletMap();

        showChannelVectorNew(); // make channel vectors once

        showRoadMap(); // show road map

        showHouseMap(); // show building structures map

        showHardSurfaceMap(); // show parking lots etc

        showImageMap();

        startplot = false; //if not set to false all the above are done eahc time

        showMap(); // show map with selected data
        // the op structure uses POINTERS to maps. These maps are being used in the thread loop
        // so the action must be locked by mutex, to ensure only one trhead can access the data

        if (doShootScreens)
            shootMultipleScreens();
    }

    //qDebug() << "GUI thread waking up model thread at" << QTime::currentTime();
    W->mutex.lock();
    W->mu_condition.wakeAll();
    W->mutex.unlock();
}
//---------------------------------------------------------------------------
void lisemqt::worldDone(const QString &results)
{
    label_debug->setText(results);
    if (results.contains("ERROR"))
        QMessageBox::critical(this,QString("openLISEM"), results, QMessageBox::Ok );

    tabWidgetOptions->setCurrentIndex(lastOptionSceen);

    tabWidget->setCurrentIndex(2);
    tabWidget_out->setCurrentIndex(0);
    shootSingleScreen(0);
    tabWidget_out->setCurrentIndex(1);
    shootSingleScreen(0);

    stoprun = true;
    startplot = false;

    // free the map plot discharge bdata
    if (QFileInfo(QString(op.userAppDir+"openlisemtmp.run")).exists())
        QFile::remove(QString(op.userAppDir+"openlisemtmp.run"));

    stopAct->setChecked(false);
    runAct->setChecked(false);
    pauseAct->setChecked(false);

    E_runFileList->setEnabled(true);
    checkDoErosion->setEnabled(true);
    label_1->setEnabled(true);
    toolButton_fileOpen->setEnabled(true);
    toolButton_deleteRun->setEnabled(true);

    // not sure if this is needed?

    // if (op.doBatchmode) {
    //     close();
    // }
}
//---------------------------------------------------------------------------
// this function is linked to the debug signal emitted from the model world
void lisemqt::worldDebug(const QString &results)
{
    QString sss = results; //label_debug->text() + results + " - ";
    label_debug->setText(sss);
    // show messages from the World model on the screen
}
//---------------------------------------------------------------------------
void lisemqt::initOP()
{
    op.Pmm.clear();
    op.Time.clear();
    op.Qtile.clear();
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

    op.ComboMaps.clear();
    op.ComboColorMap.clear();
    op.ComboColors.clear();
    op.ComboLogaritmic.clear();
    op.ComboSymColor.clear();
    op.ComboMapNames.clear();
    op.ComboUnits.clear();
    op.ComboScaling.clear();

    op.comboboxset = false;

//the maps are pointers to the real maps, not copies
    op.baseMap = nullptr;
    op.baseMapDEM = nullptr;
    op.channelMap = nullptr;
    op.outletMap = nullptr;
    op.roadMap = nullptr;
    op.houseMap = nullptr;
    op.hardsurfaceMap = nullptr;
    op.Image = nullptr;

    op.EndPointX.clear();
    op.EndPointY.clear();
    op.ObsPointX.clear();
    op.ObsPointY.clear();

    op.runstep = 0;
    op.printstep = 0;
    op.maxstep = 0;
    op.CatchmentArea = 1.0;
    op._dx = 1.0;
    op._llx = 0;
    op._lly = 0;
    op._nrRows = 10;
    op._nrCols = 10;
    op.t = 0;
    op.time = 0;
    op.maxtime = 0;
    op.EndTime = 0;
    op.BeginTime = 0;
    op.MB = 0;
    op.Qtot = 0;
    //op.Qtile = 0;
    op.Qtiletot = 0;
    op.RainpeakTime = 0;
    op.RunoffFraction = 0;
    op.FloodTotMax = 0;
    op.FloodAreaMax = 0;
    op.BaseFlowTotmm = 0;
    op.IntercLitterTotmm = 0;
    op.Qtotmm = 0;
    op.IntercTotmm = 0;
    op.IntercHouseTotmm = 0;
    op.WaterVolTotmm = 0;
    op.StormDrainTotmm = 0;
    op.InfilTotmm = 0;
    op.RainTotmm = 0;
    op.ETaTotmm = 0;
    op.GWlevel = 0;
    op.Theta1 = 0;
    op.Theta2 = 0;
    op.SurfStormm = 0;
    op.InfilKWTotmm = 0;
    op.WHflood = 0;
    op.MBs = 0;
    op.DetTot = 0;
    op.DetTotSplash = 0;
    op.DetTotFlow = 0;
    op.DepTot = 0;
    op.SoilLossTot = 0;
    op.SedTot = 0;
    op.ChannelVolTotmm = 0;
    op.RetentionVolTot = 0;
    op.RetentionVolTotmm = 0;
    op.ChannelSedTot = 0;
    op.ChannelDepTot = 0;
    op.ChannelDetTot = 0;
    op.ChannelWH = 0;
    op.FloodSedTot = 0;
    op.FloodDepTot = 0;
    op.FloodDetTot = 0;
    op.FloodVolmm = 0;
    op.format = "PCRaster";

}
//---------------------------------------------------------------------------

