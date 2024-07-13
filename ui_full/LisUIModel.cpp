/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
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

#define Clear(list) while (list.count()) {if (!list.isEmpty()) list.removeLast();}
#define Clear2D(list) for(int i = 0 ; i < list.count(); i++) {\
            while (list.at(i)->count()) {if (!list.at(i)->isEmpty())list.at(i)->removeLast();}\
            while (list.count()) {if (!list.isEmpty()) list.removeLast();}}



void lisemqt::ClearOP()
{
    /*
    QList<int> OutletIndices;
    QList<int> OutletLocationX;
    QList<int> OutletLocationY;
    QList<QVector<double>*> OutletQ;
    QList<QVector<double>*> Wavein;
    QList<QVector<double>*> OutletQs;  //current kg/s
    QList<QVector<double>*> OutletC;   // avg concetration
    QList<QVector<double>*> OutletChannelWH;
    QVector<double> OutletQpeak;
    QVector<double> OutletQpeaktime;
    QVector<double> OutletQtot;
    QVector<double> OutletQstot;  // sum in kg
    QVector<double> Pmm;
    QVector<double> Time;
    QVector <double> Qtile;
    QVector <double> EndPointX;
    QVector <double> EndPointY;
    QVector <double> ObsPointX;
    QVector <double> ObsPointY;
    QVector <LDD_COORIN> lddch_;

    // map pointers for display
    cTMap *baseMap;
    cTMap *baseMapDEM;
    cTMap *channelMap;
    cTMap *outletMap;
    cTMap *roadMap;
    cTMap *houseMap;
    cTMap *hardsurfaceMap;
    cTRGBMap *Image;

    QList<int> ComboLists;
    QList<cTMap *> ComboMaps;
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
    */

    Clear2D(op.OutletQ);
    Clear2D(op.Wavein);
    Clear2D(op.OutletQs);
    Clear2D(op.OutletC);
    Clear2D(op.OutletChannelWH);

    Clear(op.OutletIndices);
    Clear(op.OutletLocationX);
    Clear(op.OutletLocationY);
    Clear(op.OutletQpeak);
    Clear(op.OutletQpeaktime);
    Clear(op.OutletQtot);
    Clear(op.OutletQstot);
    Clear(op.Pmm);
    Clear(op.Time);
    Clear(op.Qtile);
    Clear(op.EndPointX);
    Clear(op.EndPointY);
    Clear(op.ObsPointX);
    Clear(op.ObsPointY);
    Clear(op.lddch_);
    // Clear(op.ComboLists);
    // qDeleteAll(op.ComboMaps.begin(), op.ComboMaps.end());
    // op.ComboMapNames.clear();
    // op.ComboUnits.clear();
    // Clear(op.ComboLogaritmic);
    // Clear(op.ComboSymColor);
    // Clear(op.ComboScaling);
    // Clear(op.userMinV);
    // Clear(op.userMaxV);
    // Clear(op.comboStep);
    // Clear(op.comboStep);
    delete op.baseMap;
    delete op.baseMapDEM;
    delete op.channelMap;
    delete op.outletMap;
    delete op.roadMap;
    delete op.houseMap;
    delete op.hardsurfaceMap;
    delete op.Image;
}

//---------------------------------------------------------------------------
/** Run the model:
Save the current interface as a temporary run file, read by the model
Make the model world and run it
*/
void lisemqt::runmodel()
{

    if (W)
    {
        if (W->waitRequested) {
            pausemodel();
            qDebug() << "pauze";
            return;
        }
    }

    // if the model has stopped and a new run is requested, clear the datastructures
    // until that time the user can look at the old results
    if (stoprun && W) {
        // destroy ALL maps
        qDeleteAll(W->maplistCTMap.begin(),W->maplistCTMap.end());
        W->maplistCTMap.clear();
        // destroy all networlk structures
        Clear(W->cr_);
        Clear(W->crch_);
        Clear(W->crlinkedldd_);
        Clear(W->crlinkedlddch_);
        Clear(W->crldd5_);
        Clear(W->crlddch5_);
        Clear(W->crout_);
        Clear(W->dcr_);
        Clear(W->crtile_);

        QVector <double> zero;
        zero.clear();
        PGraph->setSamples(zero,zero);
        QGraph->setSamples(zero,zero);
        QbGraph->setSamples(zero,zero);
        QsGraph->setSamples(zero,zero);
        CGraph->setSamples(zero,zero);
        QtileGraph->setSamples(zero,zero);
        HPlot->replot();

        ClearOP(); // clear most of the op structure

        // destroy swatre data structures
        if (W->InfilMethod == INFIL_SWATRE && W->initSwatreStructure)
        {
            W->FreeSwatreInfo();
            if (W->SwatreSoilModel)
                W->CloseSwatre(W->SwatreSoilModel);
            if (W->SwatreSoilModelCrust)
                W->CloseSwatre(W->SwatreSoilModelCrust);
            if (W->SwatreSoilModelCompact)
                W->CloseSwatre(W->SwatreSoilModelCompact);
            if (W->SwatreSoilModelGrass)
                W->CloseSwatre(W->SwatreSoilModelGrass);
        }
    }

    startplot = true; // user has pressed run, used only to initiatte screen stop, after that set to false!
    stoprun = false; // user has not stopped the run

    label_debug->text().clear();

    //NOTE op.runfilename is set in function openRunFile()
    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this,"openLISEM",QString("Load a runfile first!"));
        return;
    }

    lastOptionSceen = tabWidgetOptions->currentIndex();

    showOutputDataZero();

    updateModelData();
    QFile f(QString(op.LisemDir+"openlisemtmp.run"));
    if (f.exists())
        f.remove();

    savefile(QString(op.LisemDir+"openlisemtmp.run"));
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

    checkMapBuildings->setChecked(false);
    checkMapBuildings->setEnabled(checkHouses->isChecked());

    checkMapRoads->setChecked(false);
    checkMapRoads->setEnabled(checkRoadsystem->isChecked());

    checkMapHardSurface->setChecked(false);
    checkMapHardSurface->setEnabled(checkHardsurface->isChecked());
    transparencyHardSurface->setEnabled(checkHouses->isChecked() || checkRoadsystem->isChecked() || checkHardsurface->isChecked());
    //transparencyHardSurface->setVisible(false);

    sedgroup->setVisible(checkDoErosion->isChecked());
    tabWidget_totout->setTabEnabled(1,checkDoErosion->isChecked() );

    showInfoAct->setChecked(true);
    setOutputInfo(true); // show the cursor over the map

    // initialize output graphs
    initPlot();

    initMapPlot();

    initOutputData();

    initOP();
    // reset op structure

    showOutputData();

    //=======================================================================================//

    // moved to lisemqt, nneeds to be done only once
    // W = new TWorld();
    // // make the model world !!!
    // connect(W, SIGNAL(show(bool)),this, SLOT(worldShow(bool)),Qt::BlockingQueuedConnection);
    // connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
    // connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
    // connect(W, SIGNAL(timedb(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
    // // connect emitted signals from the model thread to the interface routines that handle them

    W->showInfo = true;

    //WhasStopped = false;
    W->stopRequested = false;
    // stoprequested is used to stop the thread with the interface
    W->waitRequested = false;
    // waitrequested is used to pause the thread with the interface, only on windows machines!
    W->noInterface = true; // if true then show something
    W->noOutput = false;
    W->batchmode = false;
    // run without Qt interface on original runfile only

    op.timeStartRun = QDateTime().currentDateTime().toString("yyMMdd-hhmm");

    if (checkAddDatetime->isChecked()) {
        screenShotDir = E_ResultDir->text() + QString("res"+op.timeStartRun+"/");
        QDir(screenShotDir).mkpath(QString("screens/"));
        screenShotDir = screenShotDir + QString("screens/");
    } else {
        screenShotDir = E_ResultDir->text();
        QDir(screenShotDir).mkpath(QString("screens"+op.timeStartRun+"/"));
        screenShotDir = screenShotDir + QString("screens"+op.timeStartRun+"/");
    }
    //qDebug() << screenShotDir;

    // take a screenshot of all option widgets
    tabWidget->setCurrentIndex(0);
    for (int i = 0; i < 9; i++) {
        tabWidgetOptions->setCurrentIndex(i);
        shootSingleScreen(1);
    }
    tabWidget->setCurrentIndex(2);
    //switch to output screen

    W->start();
    // start the model thread, executes W->run()

    E_runFileList->setEnabled(false);
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
void lisemqt::worldShow(bool showall)
{
    progressBar->setMaximum(op.maxstep);
    progressBar->setValue(op.runstep);

    startPlots(); // called once using bool startplot

    showOutputData(); // show output data for all and point x

    if (!showall)
        return;

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

    if (doShootScreens)
        shootMultipleScreens();
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


    // arrive here after model emits done signal
     if (W) {
        W->quit();
        W->wait();
    }
    stoprun = true;
    startplot = false;

    // free the map plot discharge bdata
    if (QFileInfo(QString(op.LisemDir+"openlisemtmp.run")).exists())
        QFile::remove(QString(op.LisemDir+"openlisemtmp.run"));

    // delete the temp run file
    //qDebug() << QString(op.LisemDir+"openlisemtmp.run")<< "deleted";

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
    QString sss = results; //label_debug->text() + results + " - ";
    label_debug->setText(sss);
    // show messages from the World model on the screen
}
//---------------------------------------------------------------------------
// this function is linked to the debug signal emitted from the model world
void lisemqt::worldTimedb(const QString &results)
{
   // label_realtime->setText(results);
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
   // op.OutletQb.clear();
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
    op.ComboColorMap.clear();
    op.ComboColors.clear();
    op.ComboLogaritmic.clear();
    op.ComboSymColor.clear();
    op.ComboMapNames.clear();
    op.ComboUnits.clear();
    op.ComboScaling.clear();

    op.comboboxset = false;

    // delete op.baseMap;
    // delete op.baseMapDEM;
    // delete op.channelMap;
    // delete op.outletMap;
    // delete op.roadMap;
    // delete op.houseMap;
    // delete op.hardsurfaceMap;
    // delete op.Image;
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
    //op.WaterVolTotchannelmm = 0;
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
    op.format = "PCRaster";

}
//---------------------------------------------------------------------------

