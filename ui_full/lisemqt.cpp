/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2022  Victor Jetten
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
  \file lisemqt.cpp
  \brief Main UI functions

  NOTE:
  namelist is a structure contining an exact copy of the runfile (incl empty lines etc)
  DEFmaps is an array of stringlists of the maps as they appear in the map tree structure
  maplist is a list of input maps and descriptions for the interface map tree structure,
  created and maintained automatically (DEFmaps is the interface version of maplist)

  in namelist and maplist "name" is the description and "value" is the map name (filename)

-# namelist and DEFmaps are first filled with default values (in LisUIDefaultNames.cpp)
-# maplist is created in fillMapnames()
-# runfile is read and parsed, namelist is adapted with runfile choices (in LisUIrunfile.cpp)
-# user makes changes in options and mapnames in interface
-# when run is pressed: the namelist is updated with new choices and saved to a tmp runfile for the model

update of the runfile before running:
-# in savefile (in lisemqt.cpp) the namelist is updated with all the new variables: updateModelData()
-# updateModelData() (in LisUIrunfile.cpp) puts all new choices in the namelist so that is can be saved to the runfile

*/

#include "lisemqt.h"
#include "model.h"
#include "global.h"
#include <omp.h>

output op;
// declaration of variable structure between model and interface.
// All model results are put in this structure and sent from the model
// to the interface each timestep, defined in LisUIoutput.h


//--------------------------------------------------------------------
lisemqt::lisemqt(QWidget *parent, bool doBatch, QString runname)
    : QMainWindow(parent)
{
    setupUi(this);
    // set up interface

    setMinimumSize(1280,800);
    showMaximized();

    darkLISEM = false;

    op.nrRunsDone = 0;

    int ompt = omp_get_max_threads();
    nrUserCores->setMaximum(ompt);//omp_get_max_threads());

    helpbox = new QDialog();
    helpbox->resize(qApp->primaryScreen()->size().height()*2/3,qApp->primaryScreen()->size().height()*2/3);
    helpbox->setWindowTitle("option help");
    helpLayout = new QHBoxLayout(helpbox);
    helptxt = new QTextEdit();
    helpLayout->addWidget(helptxt);

    op.runfilename.clear();
    E_runFileList->clear();

    //TODO: check all options and default values
    resetAll();
    // all options and mapnames are reset to their default names and values
    // fill DEFmaps stringlist and make mapList with default names
    // mapList will be refilled with the runfile and user choices
    // so this contains the final list of maps

    SetToolBar();
    // slots and signals

    initMapTree();
    // initalize interface and make tree structure for map names (= DEFmaps stringlist)

    defaultRunFile();
    //fill namelist with default runfile names
    //use all actual mapnames from the mapList structure

    SetConnections();

    setupPlot();
    // set up the discharge graphs

    setupMapPlot();
    // set up the raster map drawing

    Ui_lisemqtClass::statusBar->addWidget(progressBar, 1);
    // put the progress bar into the statusbar
   // this->statusBar->addPermanentWidget(progressBar);

    SetStyleUI();
    // do some style things

    tabWidgetOptions->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

    lisMpeg = new lismpeg(this);

    tabWidgetOptions->setCurrentIndex(0);
    tabWidget_OutputMaps->setCurrentIndex(0);

    doBatchmode = doBatch; // save as global var in iface
    //batchRunname = runname;
    //doCheckRainfall(true); // ???????? why here
    op.doBatchmode = doBatch;  //copy batchmode for inside run

    // make the model world once, this structure is always needed regardless of the area
    W = new TWorld();
    connect(W, SIGNAL(show(bool)),this, SLOT(worldShow(bool)),Qt::BlockingQueuedConnection);
    connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
    connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
    connect(W, SIGNAL(timedb(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
    // connect emitted signals from the model thread to the interface routines that handle them
    //startplot = false; // start plotting
    stoprun = false;
    W->waitRequested = false;
    // run is not started so we don't accidentally do wrong things while W exists

    if(doBatch)
    {
        runfilelist.clear();
        runfilelist << runname;

        op.runfilename = runname;
        GetRunfile();     // get the nrunfile and fill namelist
        ParseInputData(); // fill interface with namelist data and fill mapList
                          // also update DEFmaps for map tree view in interface
        initMapTree();    // fill the tree strcuture on page 2 with DEFmaps
        RunAllChecks();   // activate the maps in the tree parts in response to checks
        E_runFileList->insertItem(0, runname);

        stopAct->setChecked(false);
        runAct->setChecked(true);
        pauseAct->setChecked(false);
        runmodel();
    }
    else
    {
        GetStorePath();
        // this gets the lisem.ini file and reads the run files, the first runfile
        // triggers also loading the interface with all variable values
        // if ini is empty nothing happens
        //triggered in on_E_runFileList_currentIndexChanged(int)
    }
}
//--------------------------------------------------------------------
lisemqt::~lisemqt()
{
    if (!doBatchmode)
        StorePath();
    delete W;
}
//--------------------------------------------------------------------
// NAMING convention void on_<widget name="">_<signal name="">(<signal parameters="">)
// works automatically. if included here may be executed twice!!! not sure...
void lisemqt::SetConnections()
{
    //connect(checkPesticides, SIGNAL(toggled(bool)), this, SLOT(doCheckPesticides(bool)));

    connect(toolButton_fileOpen, SIGNAL(clicked()), this, SLOT(openRunFile()));
    connect(toolButton_deleteRun, SIGNAL(clicked()), this, SLOT(deleteRunFileList()));
    connect(toolButton_MapDir, SIGNAL(clicked()), this, SLOT(setMapDir()));

    connect(treeView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(openMapname(QModelIndex)));
    connect(MapNameModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(editMapname(QModelIndex, QModelIndex)));
    connect(toolButton_ResultDir, SIGNAL(clicked()), this, SLOT(setResultDir()));

   // obsolete
   // connect(checkIncludeChannel, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));
   // connect(checkOverlandFlow1D, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));
   // connect(checkOverlandFlow2Dkindyn, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));
   // connect(checkOverlandFlow2Ddyn, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));

   // connect(checkDoErosion, SIGNAL(toggled(bool)), this, SLOT(setErosionTab(bool)));

    connect(spinBoxPointtoShow,SIGNAL(valueChanged(int)),this,SLOT(onOutletChanged(int)));

    connect(checkFormatGtiff, SIGNAL(toggled(bool)), this, SLOT(setFormatMaps(bool)));

}
//--------------------------------------------------------------------
void lisemqt::setFormatMaps(bool check)
{
    QString ext = ".map";
    if (check)
    {
        op.format = "GTiff";
        ext = ".tif";
    }
    else
    {
        op.format = "PCRaster";
        ext = ".map";
    }

    E_RainfallMap->setText(QFileInfo(E_RainfallMap->text()).baseName() + ext);
    E_InterceptionMap->setText(QFileInfo(E_InterceptionMap->text()).baseName() + ext);
    E_InfiltrationMap->setText(QFileInfo(E_InfiltrationMap->text()).baseName() + ext);
    E_RunoffMap->setText(QFileInfo(E_RunoffMap->text()).baseName() + ext);
    E_DetachmentMap->setText(QFileInfo(E_DetachmentMap->text()).baseName() + ext);
    E_DepositionMap->setText(QFileInfo(E_DepositionMap->text()).baseName() + ext);
    E_SoillossMap->setText(QFileInfo(E_SoillossMap->text()).baseName() + ext);
    E_ChanDetachmentMap->setText(QFileInfo(E_ChanDetachmentMap->text()).baseName() + ext);
    E_ChanDepositionMap->setText(QFileInfo(E_ChanDepositionMap->text()).baseName() + ext);
    //E_FloodlevelMap->setText(QFileInfo(E_FloodlevelMap->text()).baseName() + ext);
    E_FloodTimeMap->setText(QFileInfo(E_FloodTimeMap->text()).baseName() + ext);
    E_FloodFEW->setText(QFileInfo(E_FloodFEW->text()).baseName() + ext);
    E_FloodmaxVMap->setText(QFileInfo(E_FloodmaxVMap->text()).baseName() + ext);
    E_FloodmaxVHMap->setText(QFileInfo(E_FloodmaxVHMap->text()).baseName() + ext);
    E_ChannelMaxQ->setText(QFileInfo(E_ChannelMaxQ->text()).baseName() + ext);
    E_ChannelMaxWH->setText(QFileInfo(E_ChannelMaxWH->text()).baseName() + ext);
    E_ChannelQtotm3Map->setText(QFileInfo(E_ChannelQtotm3Map->text()).baseName() + ext);

    E_stormDrainMap->setText(QFileInfo(E_stormDrainMap->text()).baseName() + ext);
}

//--------------------------------------------------------------------
void lisemqt::on_tabWidget_out_currentChanged(int index)
{
    groupBox_drawMap->setEnabled(index == 1);
    groupBox_drawMap->setVisible(true);
    groupBox_info->setVisible(true);
/*
    if (this->height() > 800)
    {
        groupBox_drawMap->setVisible(true);
        groupBox_info->setVisible(true);
    }
    else
    {
        if (index == 0)//tabWidget_out->currentIndex() == 0)
        {
            groupBox_drawMap->setVisible(false);
            groupBox_info->setVisible(true);
        }

        else
        {
            groupBox_drawMap->setVisible(true);
            groupBox_info->setVisible(false);
        }
    }
    */
}
//--------------------------------------------------------------------
void lisemqt::setErosionMapOutput(bool doit)
{
    ComboMinSpinBox2->setEnabled(doit);
    ComboMaxSpinBox2->setEnabled(doit);
    DisplayComboBox2->setEnabled(doit);
    checkBoxComboMaps2->setEnabled(doit);
}

//--------------------------------------------------------------------
//gives values 0,1,2,4,6,8
void lisemqt::on_nrUserCores_valueChanged(int d)
{
    int cores = cpucores;
    cores = d;
    if (d > 1) {
        if (d % 2 == 1)  {
            if (d < cpucores) cores--;
            else
                if (d > cpucores) cores++;
        } else
            cores = d;
    }
    nrUserCores->setValue(cores);
    cpucores = cores;
}
//--------------------------------------------------------------------
void lisemqt::on_ComboMinSpinBox_valueChanged(double d)
{
    int i = DisplayComboBox->currentIndex();

    if( i > -1 && i < this->SymList.length())
    {
        if (op.userMaxV.at(i) == 0)
            op.userMinV.replace(i, d);

        if (op.userMaxV.at(i) > 0 && d < op.userMaxV.at(i))
            op.userMinV.replace(i, d);

        if(op.userMaxV.at(i) > 0 && d >= op.userMaxV.at(i))
            ComboMinSpinBox->setValue(op.userMinV.at(i));
    }
    this->showMap();

}
//--------------------------------------------------------------------
void lisemqt::on_ComboMaxSpinBox_valueChanged(double d)
{
    int i = DisplayComboBox->currentIndex();

    if( i > -1 && i < this->SymList.length())
    {
        op.userMaxV.replace(i, d);

        if (op.userMaxV.at(i) > 0)
            if(op.userMinV.at(i) > op.userMaxV.at(i))
            {
                op.userMinV.replace(i, 0);
                ComboMinSpinBox->setValue(0);
            }
    }
    this->showMap();

}
//--------------------------------------------------------------------
void lisemqt::on_ComboMinSpinBox2_valueChanged(double d)
{
    if (!DisplayComboBox2->isEnabled())
        return;

    int i = IndexList1.at(DisplayComboBox2->currentIndex());

    if( i > -1 && i < this->SymList.length())
    {
        ComboMinSpinBox2->setEnabled(true);
        if (op.userMaxV.at(i) == 0)
            op.userMinV.replace(i, d);

        if (op.userMaxV.at(i) > 0 && d < op.userMaxV.at(i))
            op.userMinV.replace(i, d);

        if(op.userMaxV.at(i) > 0 && d >= op.userMaxV.at(i))
            ComboMinSpinBox2->setValue(op.userMinV.at(i));
    }
    this->showMap();

}
//--------------------------------------------------------------------
void lisemqt::on_ComboMaxSpinBox2_valueChanged(double d)
{
    if (!DisplayComboBox2->isEnabled())
        return;

    int i = IndexList1.at(DisplayComboBox2->currentIndex());

    if( i > -1 && i < this->SymList.length())   //needed?
    {
        op.userMaxV.replace(i, d);
//        if (op.ComboSymColor.at(i))
//        {
//            op.userMinV.replace(i, -d);
//            ComboMinSpinBox2->setEnabled(false);
//        }
//        else
            ComboMinSpinBox2->setEnabled(true);


        if (op.userMaxV.at(i) > 0)
            if(op.userMinV.at(i) > op.userMaxV.at(i))
            {
                op.userMinV.replace(i, 0);
                ComboMinSpinBox2->setValue(0);
            }

    }
    this->showMap();

}
//--------------------------------------------------------------------
void lisemqt::on_checkBoxComboMaps_stateChanged(int)
{
    ActiveList = -1;
    if (checkBoxComboMaps->isChecked())
    {
        if (IndexList.count() > 0)
            ActiveList = 0;
        setDisplayComboBoxes();
    }
}
//--------------------------------------------------------------------
void lisemqt::on_checkBoxComboMaps2_stateChanged(int)
{
    ActiveList = -1;
    if (checkBoxComboMaps2->isChecked())
    {
        if (IndexList1.count() > 0)
            ActiveList = 1;
        setDisplayComboBoxes();
    }
}
//--------------------------------------------------------------------
void lisemqt::setDisplayComboBoxes()
{
    if (ActiveList == 0)
    {
    }
    if (ActiveList == 1)
    {
        if (IndexList1.count() > 0)
        {
            int j = IndexList1.at(DisplayComboBox2->currentIndex());
          //  ComboMinSpinBox2->setEnabled(!op.ComboSymColor.at(j));
        }
    }

    if (ActiveList > -1)
        this->showMap();
}
//--------------------------------------------------------------------
void lisemqt::on_DisplayComboBox_currentIndexChanged(int j)
{

    if (j < 0)
        return;

    int i = IndexList.at(j);

    ComboMaxSpinBox->setValue(op.userMaxV.at(i));
    ComboMinSpinBox->setValue(op.userMinV.at(i));
    this->showMap();
}
//--------------------------------------------------------------------
void lisemqt::on_DisplayComboBox2_currentIndexChanged(int j)
{

    if (j < 0)
        return;

    int i = IndexList1.at(j);

    ComboMaxSpinBox2->setValue(op.userMaxV.at(i));
    ComboMinSpinBox2->setValue(op.userMinV.at(i));
    this->showMap();
}
//--------------------------------------------------------------------
void lisemqt::setFloodTab(bool yes)
{
    yes = true;
    if (/*checkOverlandFlow2Dkindyn->isChecked()*/ E_OFWaveType->currentIndex() == 1 && !checkIncludeChannel->isChecked()) {
        yes = false;
        QMessageBox::warning(this,"openLISEM",QString("The combination of 1D overland flow and 2D flood can only be used with a channel activated."));
        //checkOverlandFlow1D->setChecked(true);
    }
    if (E_OFWaveType->currentIndex() == 0 /*checkOverlandFlow1D->isChecked()*/) {
        yes = false;
    }

    checkDiffusion->setEnabled(yes);

    groupFloodParams->setEnabled(yes);

    outputMapsFlood->setEnabled(yes);
    label_floodVolmm->setEnabled(yes);
    label_107->setEnabled(yes);

    //if (checkOverlandFlow2Ddyn->isChecked() || checkOverlandFlow2Dkindyn->isChecked()) {
    if (E_OFWaveType->currentIndex() > 0) {
        label_107->setText(QString("Flood(h>%1mm)").arg(E_floodMinHeight->value()*1000));
        label_40->setText(QString("Runoff(h<%1mm)").arg(E_floodMinHeight->value()*1000));
    }
    else
    {
        label_107->setText("Flood");
        label_40->setText("Runoff");
    }

}
//--------------------------------------------------------------------
void lisemqt::setErosionTab(bool yes)
{
    //  yes = checkDoErosion->isChecked();

  //  tab_erosion->setEnabled(checkDoErosion->isChecked());

   // outputMapsSediment->setEnabled(checkDoErosion->isChecked());

    // checkBox_OutConc->setEnabled(checkDoErosion->isChecked());
    // checkBox_OutDet->setEnabled(checkDoErosion->isChecked());
    // checkBox_OutDep->setEnabled(checkDoErosion->isChecked());
    // checkBox_OutSL->setEnabled(checkDoErosion->isChecked());
    // checkBox_OutSed->setEnabled(checkDoErosion->isChecked());
    // checkBox_OutTC->setEnabled(checkDoErosion->isChecked());
    // checkBox_OutSedSS->setEnabled(checkDoErosion->isChecked() && checkSed2Phase->isChecked());
    // checkBox_OutSedBL->setEnabled(checkDoErosion->isChecked() && checkSed2Phase->isChecked());

    // checkBoxComboMaps2->setEnabled(checkDoErosion->isChecked());
    // ComboMinSpinBox2->setEnabled(checkDoErosion->isChecked());
    // ComboMaxSpinBox2->setEnabled(checkDoErosion->isChecked());
    // DisplayComboBox2->setEnabled(checkDoErosion->isChecked());

    // checkDiffusion->setEnabled(E_OFWaveType->currentIndex() > 0);//!checkOverlandFlow1D->isChecked());

    // // reset output to 0
    // if (!checkDoErosion->isChecked())
    // {
    //     sedgrouptotals->setEnabled(true);
    //     int dig = E_DigitsOut->value();
    //     label_MBs->setText(QString::number(0,'e',dig));
    //     label_splashdet->setText(QString::number(0,'f',dig));
    //     label_flowdet->setText(QString::number(0,'f',dig));
    //     label_sedvol->setText(QString::number(0,'f',dig));
    //     label_dep->setText(QString::number(0,'f',dig));
    //     label_detch->setText(QString::number(0,'f',dig));
    //     label_depch->setText(QString::number(0,'f',dig));
    //     label_sedvolch->setText(QString::number(0,'f',dig));
    //     //        label_flooddet->setText(QString::number(0,'f',dig));
    //     //        label_flooddep->setText(QString::number(0,'f',dig));
    //     //        label_floodsed->setText(QString::number(0,'f',dig));
    //     label_soilloss->setText(QString::number(0,'f',dig));
    //     label_soillosskgha->setText(QString::number(0,'f',dig));
    //     label_SDR->setText(QString::number(0,'f',dig));
    //     label_soillosssub->setText(QString::number(0,'f',dig));
    //     label_Qssub->setText(QString::number(0,'f',dig));
    // }
    // sedgrouptotals->setEnabled(checkDoErosion->isChecked());

    // label_soillosskgha->setEnabled(checkDoErosion->isChecked());
    // label_soilloss->setEnabled(checkDoErosion->isChecked());
    // label_SDR->setEnabled(checkDoErosion->isChecked());

}

//--------------------------------------------------------------------
//OBSOLETE
void lisemqt::setWriteOutputSOBEK(bool doit)
{
    //   checkWriteSOBEK->setChecked(!doit);
    //checkWriteCommaDelimited->setChecked(!doit);
    //checkWritePCRaster->setChecked(!doit);
}
//--------------------------------------------------------------------
//OBSOLETE
void lisemqt::setWriteOutputCSV(bool doit)
{
//    checkWriteSOBEK->setChecked(!doit);
//    checkWriteCommaDelimited->setChecked(!doit);
//    checkWritePCRaster->setChecked(!doit);
}

void lisemqt::setOutputScreen()
{
  if (W) {
    W->noInterface = !W->noInterface;
    showAllAct->setChecked(!W->noInterface);
  }
}
//--------------------------------------------------------------------
void lisemqt::setOutputInfo(bool check)
{
    if (W) {
      W->showInfo = check;
      picker->setEnabled(check);
    }
}

//--------------------------------------------------------------------
void lisemqt::SetToolBar()
{
    toolBar->setIconSize(QSize(32,32));

    openAct = new QAction(QIcon(":/2X/Folder-Open-icon.png"), "&Open a run file...", this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip("Open a run file");
    connect(openAct, SIGNAL(triggered()), this, SLOT(openRunFile()));
    toolBar->addAction(openAct);

    saveAct = new QAction(QIcon(":/2X/filesave2X.png"), "&Save the run file...", this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip("Save a run file");
    connect(saveAct, SIGNAL(triggered()), this, SLOT(saveRunFile()));
    toolBar->addAction(saveAct);

    saveasAct = new QAction(QIcon(":/2X/filesaveas2X.png"), "Save &As...", this);
    saveasAct->setShortcuts(QKeySequence::SaveAs);
    saveasAct->setStatusTip("Save a run file as ...");
    connect(saveasAct, SIGNAL(triggered()), this, SLOT(savefileas()));
    toolBar->addAction(saveasAct);
    toolBar->addSeparator();

    shootscreenAct = new QAction(QIcon(":/2X/screenshots2X.png"), "make a screenshow of the current page", this);
    connect(shootscreenAct, SIGNAL(triggered()), this, SLOT(shootScreen()));
    toolBar->addAction(shootscreenAct);

    shootMscreenAct = new QAction(QIcon(":/2X/Mscreenshots2X.png"), "Save the run in multiple screenshots", this);
    shootMscreenAct->setCheckable(true);
    connect(shootMscreenAct, SIGNAL(triggered()), this, SLOT(shootMScreen()));
    toolBar->addAction(shootMscreenAct);

    makeMovieAct = new QAction(QIcon(":/2X/film.png"), "Save the run in multiple screenshots", this);
    makeMovieAct->setCheckable(true);
    connect(makeMovieAct, SIGNAL(triggered()), this, SLOT(convertScreenshotsToVideo()));
    toolBar->addAction(makeMovieAct);

    fontIncreaseAct = new QAction(QIcon(":/2X/fontbigger2X.png"), "&Increase font size", this);
    connect(fontIncreaseAct, SIGNAL(triggered()), this, SLOT(fontIncrease()));
    toolBar->addAction(fontIncreaseAct);
    fontDecreaseAct = new QAction(QIcon(":/2X/fontsmaller2X.png"), "&Decrease font size", this);
    connect(fontDecreaseAct, SIGNAL(triggered()), this, SLOT(fontDecrease()));
    toolBar->addAction(fontDecreaseAct);

    setBWAct = new QAction(QIcon(":/black-and-white.png"), "Save the run in multiple screenshots", this);
    setBWAct->setCheckable(true);
    connect(setBWAct, SIGNAL(triggered()), this, SLOT(setBWUI()));
    toolBar->addAction(setBWAct);

    toolBar->addSeparator();
    resizeAct = new QAction(QIcon(":/2X/resetmap.png"), "&Fit map to display", this);
    connect(resizeAct, SIGNAL(triggered()), this, SLOT(resizeMap()));
    toolBar->addAction(resizeAct);

    showAllAct = new QAction(QIcon(":/2X/noscreen.png"), "&no output to screen", this);
    showAllAct->setCheckable(true);
    connect(showAllAct, SIGNAL(triggered()), this, SLOT(setOutputScreen()));
    toolBar->addAction(showAllAct);

    showInfoAct = new QAction(QIcon(":/2X/noinfo.png"), "&no info under cursor", this);
    showInfoAct->setCheckable(true);
    connect(showInfoAct, SIGNAL(triggered(bool)), this, SLOT(setOutputInfo(bool)));
    toolBar->addAction(showInfoAct);

    toolBar->addSeparator();

    runAct = new QAction(QIcon(":/2X/play-icon.png"), "Run model...", this);
    runAct->setStatusTip("run the model ...");
    runAct->setCheckable(true);
    connect(runAct, SIGNAL(triggered()), this, SLOT(runmodel()));
    toolBar->addAction(runAct);

    pauseAct = new QAction(QIcon(":/2X/pause-icon.png"), "Pause the model...", this);
    pauseAct->setStatusTip("pause the model run ...");
    connect(pauseAct, SIGNAL(triggered()), this, SLOT(pausemodel()));
    pauseAct->setCheckable(true);
    toolBar->addAction(pauseAct);

    stopAct = new QAction(QIcon(":/2X/Stop-icon.png"), "Stop the model...", this);
    stopAct->setStatusTip("stop the model run ...");
    connect(stopAct, SIGNAL(triggered()), this, SLOT(stopmodel()));
    stopAct->setCheckable(true);
    toolBar->addAction(stopAct);

    QActionGroup *runGroup = new QActionGroup(this);
    runGroup->addAction(runAct);
    runGroup->addAction(pauseAct);
    runGroup->addAction(stopAct);
    stopAct->setChecked(true);

    toolBar->addSeparator();

    aboutActI = new QAction(QIcon(":/2X/question-mark-button2x.png"), "", this);
    connect(aboutActI, SIGNAL(triggered()), this, SLOT(aboutInfo()));
    toolBar_2->addAction(aboutActI);
    restartAct = new QAction(QIcon(":/2X/reset.png"), "&Reset interface and all options...", this);
    connect(restartAct, SIGNAL(triggered()), this, SLOT(resetAll()));
    toolBar_2->addAction(restartAct);

    //toolBar->addSeparator();

    connect(checkMapBuildings, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
    connect(checkMapRoads, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
    connect(checkMapChannels, SIGNAL(clicked(bool)), this, SLOT(showChannelVector(bool)));
    connect(checkMapImage, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
    connect(checkMapHardSurface, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));

    connect(ComboMaxSpinBox,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMinSpinBox,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMaxSpinBox2,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMinSpinBox2,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));

    connect(transparency, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlpha(int)));
    connect(transparencyHardSurface, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlphaHardSurface(int)));
    connect(transparencyMap, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlphaMap(int)));
    connect(spinChannelSize, SIGNAL(valueChanged(int)),this,SLOT(ssetAlphaChannel(int)));
    connect(spinCulvertSize, SIGNAL(valueChanged(int)),this,SLOT(ssetAlphaChannelOutlet(int)));
    connect(transparencyRoad, SIGNAL(valueChanged(int)),this,SLOT(ssetAlphaHardSurfaceW(int)));
}

//--------------------------------------------------------------------
void lisemqt::setMapDir()
{
    QString path;
    QString pathin;

    pathin = findValidDir(E_MapDir->text(), false);

    path = QFileDialog::getExistingDirectory(this, QString("Select maps directory"),pathin,QFileDialog::ShowDirsOnly);
                                             //QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks |
//    path = QFileDialog::getOpenFileName(this,QString("Select maps directory"),
//                                        pathin,"maps (*.map)");
    if(!path.isEmpty()) {
        if (path.count("/") > 0 && path.lastIndexOf("/") != path.size())
            path = path + "/";
        if (path.count("\\") > 0 && path.lastIndexOf("\\") != path.size())
            path = path + "\\";
        E_MapDir->setText( path );
    }
}
//--------------------------------------------------------------------
void lisemqt::setWorkDir()
{
    //    QString path;
    //    QString pathin;

    //    pathin = E_WorkDir;//findValidDir(E_WorkDir->text(), false);

    //    path = QFileDialog::getExistingDirectory(this, QString("Select work directory"),
    //                                             pathin,
    //                                             QFileDialog::ShowDirsOnly
    //                                             | QFileDialog::DontResolveSymlinks);
    //    if(!path.isEmpty())
    //        E_WorkDir = path;//->setText( path );
}
//--------------------------------------------------------------------
void lisemqt::setResultDir()
{
    QString path;
    QString pathin;

    pathin = findValidDir(E_ResultDir->text(), true);
    if (pathin.isEmpty())
        pathin = E_WorkDir;//findValidDir( E_WorkDir->text(), false);

    path = QFileDialog::getExistingDirectory(this, QString("Select a directory to write results"),
                                             pathin,
                                             QFileDialog::ShowDirsOnly
                                             | QFileDialog::DontResolveSymlinks);
    if(!path.isEmpty()) {
        if (path.count("/") > 0 && path.lastIndexOf("/") != path.size())
            path = path + "/";
        if (path.count("\\") > 0 && path.lastIndexOf("\\") != path.size())
            path = path + "\\";
        E_ResultDir->setText( path );
    }
}


//--------------------------------------------------------------------
void lisemqt::savefileas()
{
    if (W)
    {
        QMessageBox::warning(this, "openLISEM","Cannot save a file while model is running.");
        return;
    }

    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this, "openLISEM","This runfile will habe no pathnames.");
        //return;
    }

    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    QString("Give a new runfile name"),
                                                    op.runfilename,
                                                    QString("Text Files (*.run);;All Files (*)"),
                                                    &selectedFilter);

    if (!fileName.isEmpty()) {
        updateModelData();
        savefile(fileName);
        E_runFileList->insertItem(0,fileName);
        E_runFileList->setCurrentIndex(0);
    }

}
//--------------------------------------------------------------------
void lisemqt::saveRunFile()
{
//    if (W)
//    {
//        QMessageBox::warning(this, "openLISEM","Cannot save a file while model is running.");
//        return;
//    }

    updateModelData();
    // change runfile strings with current interface options
    savefile(op.runfilename);
}
//--------------------------------------------------------------------
void lisemqt::savefile(QString name)
{
//    if (W)
//    {
//        QMessageBox::warning(this, "openLISEM","Cannot save a file while model is running.");
//        return;
//    }

    QFile fp(name);
    if (!fp.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::warning(this, QString("openLISEM"),
                             QString("Cannot write file %1:\n%2.").arg(name).arg(fp.errorString()));
        return;
    }

    QTextStream out(&fp);
    out << QString("[openLISEM runfile version 6.0]\n");

    for (int i = 1; i < nrnamelist; i++)
    {
        if (namelist[i].name.contains("[") || namelist[i].name.isEmpty())
            out << namelist[i].name << "\n"; // already contains \n
        else
            out << namelist[i].name << "=" << namelist[i].value << "\n";
        //qDebug() << namelist[i].name << "=" << namelist[i].value;
    }
    fp.close();
}
//--------------------------------------------------------------------
void lisemqt::deleteRunFileList()
{
    E_runFileList->removeItem(E_runFileList->currentIndex());
}
//--------------------------------------------------------------------
void lisemqt::openRunFile()
{
    QString path;
    QString openDir;
    if (E_runFileList->count() > 0)
        openDir = QFileInfo(E_runFileList->currentText()).dir().absolutePath();
    else
        openDir = currentDir;
    path = QFileDialog::getOpenFileName(this,
                                        QString("Select run file(s)"),
                                        openDir,
                                        QString("*.run"));

    if (path.isEmpty())
        return;

    E_runFileList->setDuplicatesEnabled(false);
    E_runFileList->setInsertPolicy(QComboBox::InsertAtTop);

    // check if it exists, if not add
    bool exst = false;
    int nr = 0;
    for (int i = 0; i < E_runFileList->count(); i++)
        if (E_runFileList->itemText(i) == path)
        {
            exst = true;
            nr = i;
        }

    if (!exst)
        E_runFileList->insertItem(0,path);

    // renew runfilenames
//    RunFileNames.clear();
//    for (int i = 0; i <= E_runFileList->count(); i++)
//        RunFileNames << E_runFileList->itemText(i);

    op.runfilename = E_runFileList->itemText(nr);
    E_runFileList->setCurrentIndex(nr);
    /* !!! this triggers runfile loading in on_E_runFileList_currentIndexChanged:
    GetRunfile();
    ParseInputData();
    FillMapList();
    RunAllChecks();
    */

}
//---------------------------------------------------------------------------
void lisemqt::GetStorePath()
{
    runfilelist.clear();
    QFile fff(op.userAppDir + "openlisem.ini");

    if (!fff.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    //   QString S = fff.readLine();
    while (!fff.atEnd())
    {
        QString  line = fff.readLine();
        if (line.contains('\n'))
            line.remove(line.size()-1,1);
        //remove '/n'
        if (line.isEmpty())
            continue;
        QFile file(line);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            continue;

        runfilelist << QString(line);
        //E_runFileList->addItem(QString(line));
    }
    fff.close();

    if (runfilelist.count() == 0)
        return;

    if (!runfilelist[0].isEmpty())
    {
        QString S = runfilelist[0];//E_runFileList->currentText();

        QFileInfo fi(S);
        QDir dir = fi.absoluteDir();
        if (dir.exists())
        {
            currentDir = dir.absolutePath();
            dir.setPath(S);
        }
    }
    E_runFileList->addItems(runfilelist);

}
//---------------------------------------------------------------------------
void lisemqt::StorePath()
{
    if (op.runfilename.isEmpty())
        return;

    QFile fff(op.userAppDir + "openlisem.ini");
    if (!fff.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream ts( &fff );
    //  ts << op.runfilename << "\n";// << endl;

    for (int i = 0; i < E_runFileList->count(); i++)
        ts << E_runFileList->itemText(i) << "\n";

    fff.close();
}

//---------------------------------------------------------------------------
void lisemqt::on_E_runFileList_currentIndexChanged(int)
{
    if (E_runFileList->count() == 0)
        return;
    if (E_runFileList->currentText() == "")
        return;
    CurrentRunFile = E_runFileList->currentIndex();
    op.runfilename = E_runFileList->currentText();
    //RunFileNames.at(CurrentRunFile);

    GetRunfile();   // get the nrunfile and fill namelist

    ParseInputData(); // fill interface with namelist data and fill mapList
    // also update DEFmaps for map tree view in interface

    initMapTree();  // fill the tree strcuture on page 2 with DEFmaps
    RunAllChecks(); // activate the maps in the tree parts in response to checks
}
//--------------------------------------------------------------------
void lisemqt::on_E_MapDir_returnPressed()
{
    QFileInfo fin(E_MapDir->text());
    if(!fin.exists())
    {
        E_MapDir->setText("");
        QMessageBox::warning(this,"openLISEM",
                             QString("Map directory does not exist"));
    }
}
//--------------------------------------------------------------------
void lisemqt::on_E_ResultDir_returnPressed()
{
    if (E_ResultDir->text().isEmpty())
        return;
    QFileInfo fin(E_ResultDir->text());
    if(!fin.exists())
    {
        int ret =
                QMessageBox::question(this, QString("openLISEM"),
                                      QString("The directory \"%1\"does not exist.\n"
                                              "Do you want to create it (apply)?")
                                      .arg(fin.absoluteFilePath()),
                                      QMessageBox::Apply |QMessageBox::Cancel,QMessageBox::Cancel);
        if (ret == QMessageBox::Apply)
            QDir(E_ResultDir->text()).mkpath(E_ResultDir->text());

    }
}

//--------------------------------------------------------------------
void lisemqt::aboutQT()
{
    QMessageBox::aboutQt ( this, "openLISEM" );
}
//--------------------------------------------------------------------
void lisemqt::aboutInfo()
{
    QMessageBox::information ( this, "openLISEM",
                               QString("openLISEM verion %9 (%10) is created by Victor Jetten and Bastian van de Bout with:\n\n%1\n%2\n%3\n%4\n%5\n%6\n%7\n%8")
                               .arg("- MSYS2 with MingW64, Qt and CMake (https://www.msys2.org/,http://qt.nokia.com/,https://cmake.org/)")
                               .arg("- Qwt technical application widgets for Qt (http://qwt.sf.net)")
                               .arg("- Flood source code derived from fullSWOF2D (http://www.univ-orleans.fr/mapmo/soft/FullSWOF/)")
                               .arg("- Using openMP for parallel processing (https://www.openmp.org/)")
                               .arg("- Using GDAL for map handling (https://gdal.org/)")
                               .arg("- PCRaster lib map functions: http://pcraster.geo.uu.nl/")
                               .arg("Details can be found here: https://github.com/vjetten/openlisem")
                               .arg("This software is made available under GNU CPLv3.0")
                               .arg(VERSIONNR)
                               .arg(DATE)
                               );
}

//--------------------------------------------------------------------
void lisemqt::resetTabOptions()
{
    //checkOverlandFlow1D->setChecked(false);
    //checkOverlandFlow2Ddyn->setChecked(true);
    //checkOverlandFlow2Dkindyn->setChecked(false);
    E_OFWaveType->setCurrentIndex(2);
    checkDoErosion->setChecked(false);

    checkIncludeChannel->setChecked(true);
    checkChannelInfil->setChecked(false);
    //checkChannelBaseflow->setChecked(false);
    groupBaseflowParams->setEnabled(true);

    checkDischargeUser->setChecked(false);
    //checkChannelAdjustCHW->setChecked(true);

    checkInfrastructure->setChecked(false);
    checkRoadsystem->setChecked(false);
    checkHouses->setChecked(false);
    checkAddBuildingDEM->setChecked(false);
    checkHardsurface->setChecked(false);
    checkRaindrum->setChecked(false);
    checkStormDrains->setChecked(false);
}
//--------------------------------------------------------------------
void lisemqt::resetTabCalibration()
{
    //calibration
    E_CalibrateSmax->setValue(1.0);
    E_CalibrateRR->setValue(1.0);
    E_CalibrateKsat->setValue(1.0);
    E_CalibrateKsat2->setValue(1.0);
    E_CalibrateN->setValue(1.0);
    E_CalibrateTheta->setValue(1.0);
    E_CalibratePsi->setValue(1.0);
    E_CalibrateSD1->setValue(1.0);
    E_CalibrateSD2->setValue(1.0);
    E_CalibrateChKsat->setValue(1.0);
    E_CalibrateChN->setValue(1.0);
    E_CalibrateWave->setValue(0.0);
    E_CalibrateChTor->setValue(1.0);
    E_CalibrateAS->setValue(1.0);
    E_CalibrateCOH->setValue(1.0);
    E_CalibrateD50->setValue(1.0);
    E_CalibrateD90->setValue(1.0);
    E_CalibrateCHCOH->setValue(1.0);
    E_CalibrateCHUcr->setValue(1.0);
    E_CalibrateCHSV->setValue(1.0);
}


void lisemqt::resetTabInterception()
{
    checkInterception->setChecked(true);
    radioButton_1->setChecked(true); //<= crops interception
    E_CanopyOpeness->setValue(0.45);
    //    E_StemflowFraction->setValue(0.054);
    checkIncludeLitter->setChecked(false);
    E_LitterSmax->setValue(1.0);
}

void lisemqt::resetTabInfiltration()
{
    checkInfiltration->setChecked(true);

    E_InfiltrationMethod->clear();
   // E_InfiltrationMethod->addItem("no Infiltration");
    E_InfiltrationMethod->addItem("SWATRE");
    E_InfiltrationMethod->addItem("Green and Ampt");
    E_InfiltrationMethod->addItem("Smith and Parlange");
    E_InfiltrationMethod->addItem("Richards equation (experimental)");
    E_InfiltrationMethod->setCurrentIndex(1);

    checkInfilCompact->setChecked(false);
    checkInfilCrust->setChecked(false);
    //checkInfil2layer->setChecked(false);
    checkInfilImpermeable->setChecked(false);
    checkIncludeTiledrains->setChecked(false);
    checkGeometric->setChecked(true);
    E_SWATREDtsecFraction->setValue(0.2);
    E_SwatreTableName->setText("profile.inp");
    E_SwatreTableDir->setText("");
}
//--------------------------------------------------------------------
void lisemqt::resetTabFlow()
{
    E_FlowBoundary->setValue(1);
    E_floodMinHeight->setValue(0.05);
    checkFloodInitial->setChecked(false);
    checkFlowBarriers->setChecked(false);
    line_FlowBarriers->setText("flowbarriers.txt");
    checkBuffers->setChecked(false);
    check2DDiagonalFlow->setChecked(true);
    //check2DDiagonalFlowNew->setChecked(false);
    checkCorrectDem->setChecked(false);
    E_pitValue->setValue(1.0);
    E_TimestepMinFlood->setValue(0.2);
    E_courantFactor->setValue(0.2);

    GW_recharge->setValue(1.0);
    GW_flow->setValue(1.0);
    GW_slope->setValue(1.0);
    GW_deep->setValue(0.0);
    GW_threshold->setValue(0.2);
}
//--------------------------------------------------------------------
void lisemqt::resetTabErosion()
{
    E_splashEquation->setValue(1);
    radioButtonKE1->setChecked(true);
    spinKEparameterA1->setValue(28.3);
    spinKEparameterB1->setValue(0.52);
    spinKEparameterC1->setValue(0.042);
    radioButtonKE2->setChecked(false);
    spinKEparameterA2->setValue(8.9);
    spinKEparameterB2->setValue(8.44);
    radioButtonKE3->setChecked(false);
    spinKEparameterA3->setValue(7.6);
    spinKEparameterB3->setValue(0.22);

    checkKETimebased->setChecked(false);

    E_EfficiencyDET->setCurrentIndex(1);
    E_EfficiencyDETCH->setCurrentIndex(1);

    E_RBLMethod->setCurrentIndex(1);
    E_RSSMethod->setCurrentIndex(1);
    E_SSMethod->setCurrentIndex(0);
    E_BLMethod->setCurrentIndex(1);

    E_SigmaDiffusion->setValue(0.5);

    checkDiffusion->setChecked(false);
    checkDiffusionCH->setChecked(false);

    checkSed2Phase->setChecked(false);

    checkSedtrap->setChecked(false);

    checkMaterialDepth->setChecked(false);
    E_DepositedCohesion->setValue(0.5);
    E_BulkDens->setValue(1500);
    //E_BulkDens2->setText("1500.00");

    E_SplashDelibery->setValue(0.1);
    checkInfilGrass->setChecked(false);
    E_GrassStripN->setValue(0.2);
    E_SedTrapN->setValue(0.8);

}

void lisemqt::resetTabAdvanced()
{
    E_FloodMaxIter->setValue(200);
    E_minWHflow->setText("0.0001");
    E_FloodReconstruction->setValue(4);  //HLL2 etc
    //E_Z2Dcorrection->setValue(1.0);  //HLL2 etc
    E_FloodFluxLimiter->setValue(1);     //minmod etc
    E_courantFactorSed->setValue(0.2);
    //checkVariableTimestep->setChecked(false);
    checkTimeavgV->setChecked(true);
    checkMB_WH->setChecked(false);
    checkLinkedList->setChecked(false);
    //checkErosionInsideLoop->setChecked(true);
    checkKinWaveChannel->setChecked(false);
    E_ChannelKinWaveDt->setValue(60.0);
    nrUserCores->setValue(0);
    checkChanMaxVelocity->setChecked(true);
    checkChannel2DflowConnect->setChecked(false);
}

void lisemqt::resetAll()
{
    MapNameModel = nullptr;

    nrUserCores->setValue(0);
    doShootScreens = false;
    startShootScreens = true;

    tabWidgetOptions->setCurrentIndex(0);

    DefaultMapnames();
    // Make the default input map list, stringlist

    fillMapnames();
    // make mapList structure according to
    // DEFmaps stringlist that is used to build the map tree interface

    E_MapDir->setText("");
    E_RainfallName->setText("");
    //E_SnowmeltName->setText("");
    E_ResultDir->setText("");
    E_satImageName->setText("");
    checksatImage->setChecked(false);
    checkAdvancedOptions->setChecked(false);

    //   checkEventBased->setChecked(true);

    checkSeparateOutput->setChecked(false);
    E_DigitsOut->setValue(3);
   // checkWritePCRnames->setChecked(true);   //map series format
    checkWritePCRaster->setChecked(false); //timeplot format
    //checkWriteCommaDelimited->setChecked(true);

    E_RainfallMap->setText("rainfall.map");
    E_InterceptionMap->setText("interception.map");
    E_InfiltrationMap->setText("infiltration.map");
    E_RunoffMap->setText("runoff.map");
    E_ChannelQtotm3Map->setText("chandism3.map");
    E_DetachmentMap->setText("detachment.map");
    E_DepositionMap->setText("deposition.map");
    E_SoillossMap->setText("soilloss.map");
    E_ChanDetachmentMap->setText("chandet.map");
    E_ChanDepositionMap->setText("chandep.map");
    E_WHmaxMap->setText("whmax.map");
    E_FloodTimeMap->setText("floodtime.map");
    E_FloodStats->setText("floodstats.csv");
    E_ChannelMaxQ->setText("chanmaxq.map");
    E_ChannelMaxWH->setText("chanmaxwh.map");
    E_FloodFEW->setText("floodstart.map");
    E_FloodmaxVMap->setText("Vmax.map");
    E_FloodmaxVHMap->setText("VHmax.map");
    E_stormDrainMap->setText("Stormdrain.map");
    E_MainTotals->setText("totals.csv");
    E_SeriesTotals->setText("totalSeries.csv");
    E_PointResults->setText("hydrographs.csv");

    E_BeginTimeDay->setText("1");
  //  E_BeginTimeMin->setText("0");
    E_EndTimeDay->setText("1");
  //  E_EndTimeMin->setText("120");
    E_Timestep->setText("20");

    checkWritePCRaster->setChecked(true);

    checkBox_OutRunoff->setChecked(false);
    checkBox_OutConc->setChecked(false);
    checkBox_OutWH->setChecked(false);
    checkBox_OutTC->setChecked(false);
    checkBox_OutDet->setChecked(false);
    checkBox_OutDep->setChecked(false);
    checkBox_OutV->setChecked(false);
    checkBox_OutInf->setChecked(false);
    checkBox_OutSurfStor->setChecked(false);
    //checkBox_OutChanVol->setChecked(false);
    checkBox_OutTiledrain->setChecked(false);
   // checkBox_OutHmx->setChecked(false);
   // checkBox_OutQf->setChecked(false);
   // checkBox_OutVf->setChecked(false);
    //checkBox_OutHmxWH->setChecked(false);

    printinterval->setValue(1);

    initOP();

    progressBar->setValue(0);

    resetTabOptions();

    checkRainfall->setChecked(true);
    checkET->setChecked(false);

    resetTabInterception();

    resetTabInfiltration();

    //flow
    resetTabFlow();

    //erosion
    resetTabErosion();

    //calibration
    resetTabCalibration();

    // advanced hidden options
    resetTabAdvanced();

    tabWidget->setCurrentIndex(0);
 //   tabWidget_out->setCurrentIndex(1);
    tabWidget_out->setCurrentIndex(0);

    checkBoxComboMaps->setEnabled(true);
    checkBoxComboMaps->setChecked(true);
    checkBoxComboMaps2->setEnabled(true);
    checkBoxComboMaps2->setChecked(false);
    checkBoxComboMaps2->setEnabled(false);

    showOutputData();
}
//--------------------------------------------------------------------
QString lisemqt::findValidDir(QString path, bool up)
{
    if (!QFileInfo(path).exists() || path.isEmpty())
        path = E_MapDir->text();

    if (!QFileInfo(path).exists() || path.isEmpty())
    {
        QDir ddir(op.runfilename);
        ddir.cdUp();
        path = ddir.absolutePath();
    }
    if (up)
    {
        QDir dir = QDir(path + "/..");
        path = dir.absolutePath();
    }
    if (!QFileInfo(path).exists() || path.isEmpty())
        path = currentDir;

    if (path.indexOf("/",1) > 0)
        path.replace("\\","/");
    else
        if (path.indexOf("\\",1) > 0)
            path.replace("/","\\");

    return (path);
}
//---------------------------------------------------------------
void lisemqt::resizeMap()
{
      if (W && tabWidget_out->currentIndex() == 1)
            changeSize();

}



