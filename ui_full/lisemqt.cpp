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

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "global.h"

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
    resize(1280, 800);

    trayIcon = new QSystemTrayIcon(this);
    trayIcon->setIcon(QIcon(":/openLisem.ico"));
    trayIcon->show();

    //    QList<int> list;
    //    list << 300 << 600;
    //    splitter->setSizes(list);

    showMaximized();

    genfontsize = 8;

    tabWidgetOptions->setCurrentIndex(0);

    MapNameModel = NULL;
    HPlot = NULL;
    MPlot = NULL;
    // smallPlot = NULL;

    resetAll();
    // all options and mapnames are reset to their default names and values
    // fill DEFmaps stringlist and make mapList with default names
    // mapList will be refilled with the runfile and user choices
    // so this contains the final list of maps

    SetToolBar();
    initMapTree();
    // initalize interface and make tree structure for map names (= DEFmaps stringlist)

    defaultRunFile();
    //fill namelist with default runfile names
    //use all actual mapnames from the mapList structure

    SetConnections();

    W = NULL;
    // initalize pointer to the world, created when run button is pushed


    SetStyleUI();
    // do some style things

    Setup3DPlot();

    setupPlot();
    // setupSmallPlot();
    // set up the discharge graphs

    setupMapPlot();
    // set up the raster map drawing

    Ui_lisemqtClass::statusBar->addWidget(progressBar, 1);
    // put the progress bar into the statusbar

    E_runFileList->clear();


    doBatchmode = doBatch;
    batchRunname = runname;

    if(doBatchmode)
    {
        runfilelist.clear();
        runfilelist << batchRunname;
        op.runfilename = runname;
        GetRunfile();   // get the nrunfile and fill namelist
        ParseInputData(); // fill interface with namelist data and fill mapList
        // also update DEFmaps for map tree view in interface
        initMapTree();  // fill the tree strcuture on page 2 with DEFmaps
        RunAllChecks(); // activate the maps in the tree parts in response to checks
        E_runFileList->insertItem(0, batchRunname);


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
    }
}
//--------------------------------------------------------------------
lisemqt::~lisemqt()
{
    if (!doBatchmode)
        StorePath();
}
//--------------------------------------------------------------------
void lisemqt::SetConnections()
{
    connect(checkRainfall, SIGNAL(toggled(bool)), this, SLOT(doCheckRainfall(bool)));
    connect(checkSnowmelt, SIGNAL(toggled(bool)), this, SLOT(doCheckSnowmelt(bool)));
    connect(checkPesticides, SIGNAL(toggled(bool)), this, SLOT(doCheckPesticides(bool)));

    // NAMING convention void on_<widget name="">_<signal name="">(<signal parameters="">)
    // works automatically. if included here may be executed twice!!! not sure...

    // not needed here, is done automagically
    //   connect(DisplayComboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(on_DisplayComboBox_currentIndexChanged(int)));
    //   connect(DisplayComboBox2,SIGNAL(currentIndexChanged(int)),this,SLOT(on_DisplayComboBox2_currentIndexChanged(int)));
    //   connect(ComboMinSpinBox, SIGNAL(valueChanged(double)), this, SLOT(on_ComboMinSpinBox_valueChanged(double)));
    //   connect(ComboMaxSpinBox, SIGNAL(valueChanged(double)), this, SLOT(on_ComboMaxSpinBox_valueChanged(double)));

    connect(checkBoxComboMaps, SIGNAL(stateChanged(int)), this, SLOT(setDisplayComboBox(int)));
    connect(checkBoxComboMaps2, SIGNAL(stateChanged(int)), this, SLOT(setDisplayComboBox2(int)));

    connect(toolButton_fileOpen, SIGNAL(clicked()), this, SLOT(openRunFile()));
    connect(toolButton_deleteRun, SIGNAL(clicked()), this, SLOT(deleteRunFileList()));
    connect(toolButton_MapDir, SIGNAL(clicked()), this, SLOT(setMapDir()));
    //OBSOLETE connect(toolButton_WorkDir, SIGNAL(clicked()), this, SLOT(setWorkDir()));

    connect(treeView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(openMapname(QModelIndex)));
    // double click on mapname opens fileopen
    connect(MapNameModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(editMapname(QModelIndex, QModelIndex)));
    // doubleclick on mapname edits mapname

    connect(toolButton_ResultDir, SIGNAL(clicked()), this, SLOT(setResultDir()));

    connect(checkWritePCRaster,SIGNAL(toggled(bool)), this, SLOT(setWriteOutputPCR(bool)));
    connect(checkWriteCommaDelimited,SIGNAL(toggled(bool)), this, SLOT(setWriteOutputPCR(bool)));
    connect(checkWriteSOBEK,SIGNAL(toggled(bool)), this, SLOT(setWriteOutputPCR(bool)));

    //    connect(checkNoErosion, SIGNAL(toggled(bool)), this, SLOT(setErosionTab(bool)));
    connect(checkDoErosion, SIGNAL(toggled(bool)), this, SLOT(setErosionTab(bool)));
    connect(checkAdvancedSediment, SIGNAL(toggled(bool)), this, SLOT(setErosionTab(bool)));

    connect(spinBoxPointtoShow,SIGNAL(valueChanged(int)),this,SLOT(onOutletChanged(int)));
    // not needed here, is done automagically
    //   connect(checkEstimateGrainSizeDistribution, SIGNAL(toggled(bool)),this,SLOT(on_checkEstimateGrainSizeDistribution_toggled(bool)));
    //   connect(checkReadGrainSizeDistribution, SIGNAL(toggled(bool)),this,SLOT(on_checkReadGrainSizeDistribution_toggled(bool)));
    //connect(checkBox_SedSingleSingle,SIGNAL(toggled(bool)),this,SLOT(on_checkBox_SedSingleSingle_toggled(bool)));
    //connect(checkBox_SedMultiSingle,SIGNAL(toggled(bool)),this,SLOT(on_checkBox_SedMultiSingle_toggled(bool)));
    //connect(checkBox_SedMultiMulti,SIGNAL(toggled(bool)),this,SLOT(on_checkBox_SedMultiMulti_toggled(bool)));

}
//--------------------------------------------------------------------
void lisemqt::resizeEvent(QResizeEvent* event)
{
   QMainWindow::resizeEvent(event);
   groupBox_drawMap->setEnabled(tabWidget_out->currentIndex() == 1);
   if (this->height() >= 850)
   {
       groupBox_drawMap->setVisible(true);
       groupBox_info->setVisible(true);

   }
   else
   {
       if (tabWidget_out->currentIndex() == 0)
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

}
//--------------------------------------------------------------------
void lisemqt::on_tabWidget_out_currentChanged(int index)
{
    groupBox_drawMap->setEnabled(index == 1);
    if (this->height() >= 850)
    {
        groupBox_drawMap->setVisible(true);
        groupBox_info->setVisible(true);
        //tabWidget_out->currentIndex() == 1);
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
}

//--------------------------------------------------------------------

// disable and enable UI stuff when subcatch selected
void lisemqt::on_spinBoxPointtoShow_valueChanged(int i)
{
    //    bool subcatch = i > 1;

    //   // subcatchgroup->setEnabled(subcatch);
    //    label_soillosssub->setEnabled(/*subcatch &&*/ checkDoErosion->isChecked());
    //    label_94->setEnabled(/*subcatch && */checkDoErosion->isChecked());

    //    label_discharge->setEnabled(!subcatch);
    //    label_qtotm3->setEnabled(!subcatch);
    //    label_54->setEnabled(!subcatch);
    //    label_26->setEnabled(!subcatch);

    //    if(checkDoErosion->isChecked())
    //    {
    //        label_soillosskgha->setEnabled(!subcatch);
    //        label_soilloss->setEnabled(!subcatch);
    //        label_SDR->setEnabled(!subcatch);
    //        label_31->setEnabled(!subcatch);
    //        label_28->setEnabled(!subcatch);
    //        label_60->setEnabled(!subcatch);
    //    }
}

//--------------------------------------------------------------------
void lisemqt::on_checkImpermeable_stateChanged(int)
{
    checkPercolation->setChecked(!checkImpermeable->isChecked());
}
void lisemqt::on_checkPercolation_stateChanged(int)
{
    checkImpermeable->setChecked(!checkPercolation->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_ComboMinSpinBox_valueChanged(double d)
{
    int i = DisplayComboBox->currentIndex();

    //    if (i < 0)
    //        return;

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

    //    if (i < 0)
    //        return;

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
        if (!op.ComboSymColor.at(i))
        {
            ComboMinSpinBox2->setEnabled(true);
            if (op.userMaxV.at(i) == 0)
                op.userMinV.replace(i, d);

            if (op.userMaxV.at(i) > 0 && d < op.userMaxV.at(i))
                op.userMinV.replace(i, d);

            if(op.userMaxV.at(i) > 0 && d >= op.userMaxV.at(i))
                ComboMinSpinBox2->setValue(op.userMinV.at(i));
        }
        else
            ComboMinSpinBox2->setEnabled(false);

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
        if (op.ComboSymColor.at(i))
        {
            op.userMinV.replace(i, -d);
            ComboMinSpinBox2->setEnabled(false);
        }
        else
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
void lisemqt::setDisplayComboBox(int i)
{
    if (i == 2)
    {
        checkBoxComboMaps->setChecked(true);
        //        ComboMaxSpinBox->setEnabled(true);
        //        ComboMinSpinBox->setEnabled(true);
        //        DisplayComboBox->setEnabled(true);

        checkBoxComboMaps2->setChecked(false);
        //        DisplayComboBox2->setEnabled(false);
        //        ComboMaxSpinBox2->setEnabled(false);
        //        ComboMinSpinBox2->setEnabled(false);

        ActiveList = 0;

    }
    this->showMap();
}
//--------------------------------------------------------------------
void lisemqt::setDisplayComboBox2(int i)
{
    if (i == 2)
    {
        checkBoxComboMaps->setChecked(false);
        //        ComboMaxSpinBox->setEnabled(false);
        //        ComboMinSpinBox->setEnabled(false);
        //        DisplayComboBox->setEnabled(false);

        checkBoxComboMaps2->setChecked(true);
        //        ComboMaxSpinBox2->setEnabled(true);
        //        ComboMinSpinBox2->setEnabled(true);
        //        DisplayComboBox2->setEnabled(true);

        int j = IndexList1.at(DisplayComboBox2->currentIndex());

        ComboMinSpinBox2->setEnabled(!op.ComboSymColor.at(j));
        ActiveList = 1;

    }
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
    checkBoxComboMaps2->setChecked(false);
    checkBoxComboMaps->setChecked(true);
    this->showMap();
}
//--------------------------------------------------------------------
void lisemqt::on_DisplayComboBox2_currentIndexChanged(int j)
{

    if (j < 0)
        return;

    int i = IndexList1.at(j);

    ComboMinSpinBox2->setEnabled(!op.ComboSymColor.at(i));
    ComboMaxSpinBox2->setValue(op.userMaxV.at(i));
    ComboMinSpinBox2->setValue(op.userMinV.at(i));
    checkBoxComboMaps->setChecked(false);
    checkBoxComboMaps2->setChecked(true);
    this->showMap();
}
//--------------------------------------------------------------------
void lisemqt::on_checkBox_SedSingleSingle_toggled(bool v)
{
    if(v)
    {
        checkBox_SedMultiSingle->setChecked(false);
        checkBox_SedMultiMulti->setChecked(false);
        sedbox1->setEnabled(false);
        sedbox2->setEnabled(false);

        E_BLMethod->setValue(0);
        E_SSMethod->setValue(1);

        E_BLMethod->setMaximum(0);
        E_BLMethod->setMinimum(0);
        E_SSMethod->setMaximum(1);
        E_SSMethod->setMinimum(1);

    }else
    {
        if(!checkBox_SedMultiSingle->isChecked() && !checkBox_SedMultiMulti->isChecked())
        {
            checkBox_SedSingleSingle->setChecked(true);
        }
    }
}
//--------------------------------------------------------------------

void lisemqt::on_checkBox_SedMultiSingle_toggled(bool v)
{
    if(v)
    {
        checkBox_SedMultiMulti->setChecked(false);
        checkBox_SedSingleSingle->setChecked(false);
        sedbox1->setEnabled(false);
        sedbox2->setEnabled(true);

        E_BLMethod->setValue(1);
        E_SSMethod->setValue(1);

        E_BLMethod->setMaximum(2);
        E_BLMethod->setMinimum(0);
        E_SSMethod->setMaximum(2);
        E_SSMethod->setMinimum(1);

    }else
    {
        if(!checkBox_SedMultiMulti->isChecked() && !checkBox_SedSingleSingle->isChecked())
        {
            checkBox_SedMultiSingle->setChecked(true);
        }
    }
}
//--------------------------------------------------------------------

void lisemqt::on_checkBox_SedMultiMulti_toggled(bool v)
{
    if(v)
    {
        checkBox_SedMultiSingle->setChecked(false);
        checkBox_SedSingleSingle->setChecked(false);
        sedbox1->setEnabled(true);
        sedbox2->setEnabled(true);

        E_BLMethod->setValue(3);
        E_SSMethod->setValue(3);

        E_BLMethod->setMaximum(3);
        E_BLMethod->setMinimum(3);
        E_SSMethod->setMaximum(3);
        E_SSMethod->setMinimum(3);
    }else
    {
        if(!checkBox_SedMultiSingle->isChecked() && !checkBox_SedSingleSingle->isChecked())
        {
            checkBox_SedMultiMulti->setChecked(true);
        }
    }
}
//--------------------------------------------------------------------
void lisemqt::on_checkEstimateGrainSizeDistribution_toggled(bool v)
{
    checkReadGrainSizeDistribution->setChecked(!v);
}
//--------------------------------------------------------------------
void lisemqt::on_checkReadGrainSizeDistribution_toggled(bool v)
{
    checkEstimateGrainSizeDistribution->setChecked(!v);
}
//--------------------------------------------------------------------
//void lisemqt::on_E_NumberClasses_valueChanged(int v)
//{
//    if(checkEstimateGrainSizeDistribution->isChecked())
//    {
//        E_DisplaySedimentClass->setMaximum(E_NumberClasses->value());
//    }else
//    {
//        if(checkReadGrainSizeDistribution->isChecked())
//        {
//            E_DisplaySedimentClass->setMaximum(E_NumberClassesMaps->value());
//        }
//    }
//}
//--------------------------------------------------------------------
//void lisemqt::on_E_NumberClassesMaps_valueChanged(int v)
//{
//    if(checkEstimateGrainSizeDistribution->isChecked())
//    {
//        E_DisplaySedimentClass->setMaximum(E_NumberClasses->value());
//    }else
//    {
//        if(checkReadGrainSizeDistribution->isChecked())
//        {
//            E_DisplaySedimentClass->setMaximum(E_NumberClassesMaps->value());
//        }
//    }
//}
//--------------------------------------------------------------------
void lisemqt::setErosionTab(bool yes)
{
    //tabWidgetOptions->setTabEnabled(4, checkDoErosion->isChecked());
    //tabWidgetOptions->setTabEnabled(5, checkAdvancedSediment->isChecked() && checkDoErosion->isChecked());

    if (checkDoErosion->isChecked())
        checkBox_SedSingleSingle->setChecked(!checkAdvancedSediment->isChecked());
    // note checkBox_SedSingleSingle is not visible but still needed


    if (checkAdvancedSediment->isChecked())
    {
        if (!checkBox_SedMultiSingle->isChecked() && !checkBox_SedMultiMulti->isChecked())
        {
            checkBox_SedMultiMulti->setChecked(false);
            checkBox_SedMultiSingle->setChecked(true);
        }

    }
    yes = checkDoErosion->isChecked();
    outputMapsSediment->setEnabled(yes);
    checkBox_OutConc->setEnabled(yes);
    checkBox_OutDet->setEnabled(yes);
    checkBox_OutDep->setEnabled(yes);
    checkBox_OutSL->setEnabled(yes);
    checkBox_OutSed->setEnabled(yes);
    checkBox_OutTC->setEnabled(yes);
    groupKineticEnergy->setEnabled(yes);


    checkBoxComboMaps2->setEnabled(yes);
    ComboMinSpinBox2->setEnabled(yes);
    ComboMaxSpinBox2->setEnabled(yes);
    DisplayComboBox2->setEnabled(yes);

    sedgroup->setEnabled(yes);

    label_31->setEnabled(yes);
    label_28->setEnabled(yes);
    label_60->setEnabled(yes);

    label_soillosskgha->setEnabled(yes);
    label_soilloss->setEnabled(yes);
    label_SDR->setEnabled(yes);

}

//--------------------------------------------------------------------
void lisemqt::setWriteOutputSOBEK(bool doit)
{
    //   checkWriteSOBEK->setChecked(!doit);
    checkWriteCommaDelimited->setChecked(!doit);
    checkWritePCRaster->setChecked(!doit);
}
//--------------------------------------------------------------------
void lisemqt::setWriteOutputCSV(bool doit)
{
    checkWriteSOBEK->setChecked(!doit);
    //   checkWriteCommaDelimited->setChecked(!doit);
    checkWritePCRaster->setChecked(!doit);
}
//--------------------------------------------------------------------
void lisemqt::setWriteOutputPCR(bool /* doit */)
{
    if (checkWriteSOBEK->isChecked())
    {
        //checkWriteSOBEK->setChecked(false);
        checkWriteCommaDelimited->setChecked(false);
        checkWritePCRaster->setChecked(false);
        //checkSeparateOutput->setChecked(true);
    }
    else
        if (checkWritePCRaster->isChecked())
        {
            checkWriteSOBEK->setChecked(false);
            checkWriteCommaDelimited->setChecked(false);
            //checkWritePCRaster->setChecked(false);
        }
        else
            if (checkWriteCommaDelimited->isChecked())
            {
                checkWriteSOBEK->setChecked(false);
                //checkWriteCommaDelimited->setChecked(false);
                checkWritePCRaster->setChecked(false);
            }
}
//--------------------------------------------------------------------
void lisemqt::SetToolBar()
{
    restartAct = new QAction(QIcon(":/refresh_24.png"), "&Reset...", this);
    connect(restartAct, SIGNAL(triggered()), this, SLOT(resetAll()));
    toolBar->addAction(restartAct);
    toolBar->addSeparator();

    openAct = new QAction(QIcon(":/fileopen.png"), "&Open a run file...", this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip("Open a run file");
    connect(openAct, SIGNAL(triggered()), this, SLOT(openRunFile()));
    toolBar->addAction(openAct);

    saveAct = new QAction(QIcon(":/filesave.png"), "&Save the run file...", this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip("Save a run file");
    connect(saveAct, SIGNAL(triggered()), this, SLOT(saveRunFile()));
    toolBar->addAction(saveAct);

    saveasAct = new QAction(QIcon(":/filesaveas.png"), "Save &As...", this);
    saveasAct->setShortcuts(QKeySequence::SaveAs);
    saveasAct->setStatusTip("Save a run file as ...");
    connect(saveasAct, SIGNAL(triggered()), this, SLOT(savefileas()));
    toolBar->addAction(saveasAct);
    toolBar->addSeparator();

    shootscreenAct = new QAction(QIcon(":/screenshots.png"), "make a screendump of the current page", this);
    connect(shootscreenAct, SIGNAL(triggered()), this, SLOT(shootScreen()));
    toolBar->addAction(shootscreenAct);

    shootMscreenAct = new QAction(QIcon(":/Mscreenshots.png"), "Save the run in multiple screendumps", this);
    shootMscreenAct->setCheckable(true);
    connect(shootMscreenAct, SIGNAL(triggered()), this, SLOT(shootMScreen()));
    toolBar->addAction(shootMscreenAct);


    fontAct = new QAction(QIcon(":/fontselect.png"), "&Select font", this);
    connect(fontAct, SIGNAL(triggered()), this, SLOT(fontSelect()));
    toolBar->addAction(fontAct);

    fontIncreaseAct = new QAction(QIcon(":/fontbigger.png"), "&Increase font size", this);
    connect(fontIncreaseAct, SIGNAL(triggered()), this, SLOT(fontIncrease()));
    toolBar->addAction(fontIncreaseAct);
    fontDecreaseAct = new QAction(QIcon(":/fontsmaller.png"), "&Decrease font size", this);
    connect(fontDecreaseAct, SIGNAL(triggered()), this, SLOT(fontDecrease()));
    toolBar->addAction(fontDecreaseAct);

    toolBar->addSeparator();

    runAct = new QAction(QIcon(":/start1.png"), "Run model...", this);
    runAct->setStatusTip("run the model ...");
    runAct->setCheckable(true);
    connect(runAct, SIGNAL(triggered()), this, SLOT(runmodel()));
    toolBar->addAction(runAct);

    pauseAct = new QAction(QIcon(":/pause2.png"), "Pause the model...", this);
    pauseAct->setStatusTip("pause the model run ...");
    connect(pauseAct, SIGNAL(triggered()), this, SLOT(pausemodel()));
    pauseAct->setCheckable(true);
    toolBar->addAction(pauseAct);

    stopAct = new QAction(QIcon(":/stop1.png"), "Stop the model...", this);
    stopAct->setStatusTip("stop the model run ...");
    connect(stopAct, SIGNAL(triggered()), this, SLOT(stopmodel()));
    stopAct->setCheckable(true);
    toolBar->addAction(stopAct);

    QActionGroup *runGroup = new QActionGroup(this);
    runGroup->addAction(runAct);
    runGroup->addAction(pauseAct);
    runGroup->addAction(stopAct);
    stopAct->setChecked(true);

    aboutActI = new QAction(QIcon(":/Info.png"), "", this);
    connect(aboutActI, SIGNAL(triggered()), this, SLOT(aboutInfo()));
    toolBar_2->addAction(aboutActI);

    //aboutAct = new QAction(QIcon(":/Info.png"), "", this);
    //connect(aboutAct, SIGNAL(triggered()), this, SLOT(aboutQT()));
    //toolBar_2->addAction(aboutAct);

    toolBar_2->setMovable( false);
    toolBar->setMovable( false);

    connect(checkMapBuildings, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
    connect(checkMapRoads, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
    connect(checkMapChannels, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));

    connect(ComboMaxSpinBox,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMinSpinBox,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMaxSpinBox2,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMinSpinBox2,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));

    connect(transparency, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlpha(int)));

    connect(transparency2, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlpha2(int)));
    connect(transparency3, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlpha3(int)));
    connect(transparency4, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlpha4(int)));
}
//---------------------------------------------------------------------------
/// make some labels yellow
void lisemqt::SetStyleUI()
{
    checkBox_SedSingleSingle->setVisible(false);
    //label_103->setVisible(false);
    E_CalibratePsi->setVisible(false);
    label_77->setVisible(false);
    label_79->setVisible(false);
    E_floodMinHeight->setVisible(false);
    label_99->setVisible(false);
    checkRainfallFlood->setVisible(false);
    E_RainFloodGradient->setVisible(false);
    checkFloodInitial->setVisible(false);
    // interface elements that are not visible for now

    //groupBoxTime->setMaximumWidth(128);


    int w = 80, h = 15;//2*genfontsize;
    label_dx->setMinimumSize(w,h);
    label_area->setMinimumSize(w,h);
    label_time->setMinimumSize(w,h);
    label_endtime->setMinimumSize(w,h);
    label_raintot->setMinimumSize(w,h);
    label_watervoltot->setMinimumSize(w,h);
    label_qtot->setMinimumSize(w,h);
    label_infiltot->setMinimumSize(w,h);
    label_surfstor->setMinimumSize(w,h);
    label_interctot->setMinimumSize(w,h);
    //label_qtotm3->setMinimumSize(w,h);
    label_qpeaktime->setMinimumSize(w,h);
    label_ppeaktime->setMinimumSize(w,h);
    label_QPfrac->setMinimumSize(w,h);
    label_watervolchannel->setMinimumSize(w,h);
    label_litterstore->setMinimumSize(w,h);
    label_baseflowtot->setMinimumSize(w,h);

    label_qtotm3sub->setMinimumSize(w,h);
    label_dischargesub->setMinimumSize(w,h);
    label_qpeaksub->setMinimumSize(w,h);
    label_soillosssub->setMinimumSize(w,h);

    label_splashdet->setMinimumSize(w,h);
    label_flowdet->setMinimumSize(w,h);
    label_sedvol->setMinimumSize(w,h);
    label_dep->setMinimumSize(w,h);
    label_detch->setMinimumSize(w,h);
    label_depch->setMinimumSize(w,h);
    label_sedvolch->setMinimumSize(w,h);
    label_soilloss->setMinimumSize(w,h);
    label_soillosskgha->setMinimumSize(w,h);
    label_SDR->setMinimumSize(w,h);

    label_buffervol->setMinimumSize(w,h);
    label_buffersed->setMinimumSize(w,h);
    label_MBs->setMinimumSize(w,h);
    label_MB->setMinimumSize(w,h);


    label_dx->setStyleSheet("* { background-color: #ffffff }");
    label_area->setStyleSheet("* { background-color: #ffffff }");
    label_time->setStyleSheet("* { background-color: #ffffff }");
    label_endtime->setStyleSheet("* { background-color: #ffffff }");
    label_raintot->setStyleSheet("* { background-color: #ffff77 }");
    label_watervoltot->setStyleSheet("* { background-color: #ffff77 }");
    label_qtot->setStyleSheet("* { background-color: #ffff77 }");
    label_infiltot->setStyleSheet("* { background-color: #ffff77 }");
    label_surfstor->setStyleSheet("* { background-color: #ffff77 }");
    label_interctot->setStyleSheet("* { background-color: #ffff77 }");
    //label_qtotm3->setStyleSheet("* { background-color: #ffff77 }");
    label_qpeaktime->setStyleSheet("* { background-color: #ffff77 }");
    label_ppeaktime->setStyleSheet("* { background-color: #ffff77 }");
    label_QPfrac->setStyleSheet("* { background-color: #ffff77 }");
    //label_discharge->setStyleSheet("* { background-color: #ffff77 }");
    label_watervolchannel->setStyleSheet("* { background-color: #ffff77 }");
    label_litterstore->setStyleSheet("* { background-color: #ffff77 }");
    label_baseflowtot->setStyleSheet("* { background-color: #ffff77 }");

    label_qtotm3sub->setStyleSheet("* { background-color: #ffff77 }");
    label_dischargesub->setStyleSheet("* { background-color: #ffff77 }");
    label_qpeaksub->setStyleSheet("* { background-color: #ffff77 }");
    label_soillosssub->setStyleSheet("* { background-color: #ffff77 }");

    label_splashdet->setStyleSheet("* { background-color: #ffff77 }");
    label_flowdet->setStyleSheet("* { background-color: #ffff77 }");
    label_sedvol->setStyleSheet("* { background-color: #ffff77 }");
    label_dep->setStyleSheet("* { background-color: #ffff77 }");
    label_detch->setStyleSheet("* { background-color: #ffff77 }");
    label_depch->setStyleSheet("* { background-color: #ffff77 }");
    label_sedvolch->setStyleSheet("* { background-color: #ffff77 }");
    label_soilloss->setStyleSheet("* { background-color: #ffff77 }");
    label_soillosskgha->setStyleSheet("* { background-color: #ffff77 }");
    label_SDR->setStyleSheet("* { background-color: #ffff77 }");
    label_buffervol->setStyleSheet("* { background-color: #ffff77 }");
    label_buffersed->setStyleSheet("* { background-color: #ffff77 }");
}
//--------------------------------------------------------------------
void lisemqt::setMapDir()
{
    QString path;
    QString pathin;

    pathin = findValidDir(E_MapDir->text(), false);

    path = QFileDialog::getExistingDirectory(this, QString("Select maps directory"),
                                             pathin,
                                             QFileDialog::ShowDirsOnly
                                             | QFileDialog::DontResolveSymlinks);
    if(!path.isEmpty())
        E_MapDir->setText( path );
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
    if(!path.isEmpty())
        E_ResultDir->setText( path );
}

//--------------------------------------------------------------------
// this is for the directory with the table files
void lisemqt::on_toolButton_SwatreTableDir_clicked()
{
    QString path;
    QString pathin;

    pathin = findValidDir(E_SwatreTableDir->text(), false);

    path = QFileDialog::getExistingDirectory(this, QString("Select the directory with the Swatre tables"),
                                             pathin,
                                             QFileDialog::ShowDirsOnly
                                             | QFileDialog::DontResolveSymlinks);
    if(!path.isEmpty())
    {
        E_SwatreTableDir->setText( path );
        SwatreTableDir = path;
    }
}
//--------------------------------------------------------------------
// this is for the file profile.inp
void lisemqt::on_toolButton_SwatreTableFile_clicked()
{
    QString path;
    path = QFileDialog::getOpenFileName(this,
                                        QString("Select thee SWATRE profile definition file"),
                                        SwatreTableName,"Profiles (*.inp);;All files (*.*)");
    if(!path.isEmpty())
    {
        SwatreTableName = path;
        E_SwatreTableName->setText(path);
    }
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SwatreTableShow_clicked()
{
    QFile file(SwatreTableName);
    if (!file.open(QFile::ReadOnly | QFile::Text))
    {
        QMessageBox::warning(this,"openLISEM",
                             QString("Cannot read file %1:\n%2.")
                             .arg(SwatreTableName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);

    QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
    view->setWindowTitle(SwatreTableName);
    view->setMinimumWidth(400);
    view->setMinimumHeight(500);
    view->setAttribute(Qt::WA_DeleteOnClose);
    view->show();

    file.close();
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallName_clicked()
{
    QString path;

    RainFileDir = findValidDir(RainFileDir, false);

    path = QFileDialog::getOpenFileName(this,
                                        QString("Select rainfall file"),
                                        RainFileDir);
    if(!path.isEmpty())
    {
        QFileInfo fi(path);
        RainFileName = fi.fileName();
        RainFileDir = CheckDir(fi.absolutePath());//Dir().path());
        E_RainfallName->setText( RainFileDir + RainFileName );
    }
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltName_clicked()
{
    QString path;

    SnowmeltFileDir = findValidDir(SnowmeltFileDir, false);

    path = QFileDialog::getOpenFileName(this,
                                        QString("Select snow melt file"),
                                        SnowmeltFileDir);
    if(!path.isEmpty())
    {
        QFileInfo fi(path);
        SnowmeltFileName = fi.fileName();
        SnowmeltFileDir = CheckDir(fi.absolutePath());//Dir().path());
        E_SnowmeltName->setText( SnowmeltFileDir + SnowmeltFileName );
    }
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltShow_clicked()
{
    QFile file(SnowmeltFileDir + SnowmeltFileName);
    if (!file.open(QFile::ReadOnly | QFile::Text))
    {
        QMessageBox::warning(this,"openLISEM",
                             QString("Cannot read file %1:\n%2.")
                             .arg(SnowmeltFileDir + SnowmeltFileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);

    QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
    view->setWindowTitle(SnowmeltFileName);
    view->setMinimumWidth(400);
    view->setMinimumHeight(500);
    view->setAttribute(Qt::WA_DeleteOnClose);
    view->show();

    file.close();
}
//--------------------------------------------------------------------

void lisemqt::on_toolButton_RainfallShow_clicked()
{

    QFile file(RainFileDir + RainFileName);
    if (!file.open(QFile::ReadWrite | QFile::Text))
    {
        QMessageBox::warning(this, QString("openLISEM"),
                             QString("Cannot read file %1:\n%2.")
                             .arg(RainFileDir + RainFileName)
                             .arg(file.errorString()));
        return;
    }

    QString modifiedContents;
    QTextStream in(&file);

    QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
    view->setWindowTitle(RainFileName);
    view->setMinimumWidth(400);
    view->setMinimumHeight(500);
    view->setAttribute(Qt::WA_DeleteOnClose);
    view->show();
    if (view->document()->isModified())
    {
        int ret =
                QMessageBox::question(this, QString("openLISEM"),
                                      QString("You have modified the contents of this file.\n"
                                              "Do you want to save it?"),
                                      QMessageBox::Ok |QMessageBox::Cancel,QMessageBox::Cancel);
        if (ret == QMessageBox::Ok)
        {
            // Don't take the address of a temporary!
            // in.setString(&view->toPlainText());
            modifiedContents = view->toPlainText();
            in.setString(&modifiedContents);
        }
    }

    file.close();
}
//--------------------------------------------------------------------
void lisemqt::savefileas()
{
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
    //options);
    if (!fileName.isEmpty())
        savefile(fileName);

}
//--------------------------------------------------------------------
void lisemqt::saveRunFile()
{
    savefile(op.runfilename);
}
//--------------------------------------------------------------------
void lisemqt::savefile(QString name)
{
    updateModelData();
    // change runfile strings with current interface options

    QFile fp(name);
    if (!fp.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::warning(this, QString("openLISEM"),
                             QString("Cannot write file %1:\n%2.").arg(name).arg(fp.errorString()));
        return;
    }

    QTextStream out(&fp);
    out << QString("[openLISEM runfile version 1.0]\n");

    for (int i = 1; i < nrnamelist; i++)
    {
        if (namelist[i].name.contains("[") || namelist[i].name.isEmpty())
            out << namelist[i].name << "\n"; // already contains \n
        else
            out << namelist[i].name << "=" << namelist[i].value << "\n";
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

    RunFileNames.clear();
    for (int i = 0; i <= E_runFileList->count(); i++)
        RunFileNames << E_runFileList->itemText(i);
    //    RunFileNames.removeDuplicates();

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
    QFile fff(op.LisemDir + "openlisem.ini");
    if (!fff.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    //   QString S = fff.readLine();
    while (!fff.atEnd())
    {
        QString  line = fff.readLine();
        if (line.contains('\n'))
            line.remove(line.count()-1,1);
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

    QFile fff(op.LisemDir + "openlisem.ini");
    if (!fff.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream ts( &fff );
    //  ts << op.runfilename << "\n";// << endl;
    for (int i = 0; i < E_runFileList->count(); i++)
        ts << E_runFileList->itemText(i) << "\n";

    fff.close();
}
//---------------------------------------------------------------------------
void lisemqt::on_toolButton_ShowRunfile_clicked()
{

    QFile file(op.runfilename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "openLISEM",
                             QString("Cannot read file %1:\n%2.")
                             .arg(op.runfilename)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);



    QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
    view->createStandardContextMenu ();
    view->setWindowTitle(op.runfilename);
    view->setMinimumWidth(400);
    view->setMinimumHeight(500);
    view->setAttribute(Qt::WA_DeleteOnClose);

    view->show();

    file.close();
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
void lisemqt::shootMScreen()
{
    doShootScreens = shootMscreenAct->isChecked();
}

//--------------------------------------------------------------------
void lisemqt::shootScreen()
{

    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this, "openLISEM",QString("Select a run file first"));
        return;
    }

    QPixmap originalPixmap; // clear image for low memory situations
    QString format = "png";
    QFileInfo fi(op.runfilename);
    QString fileName;
    QString outdir = CheckDir(E_ResultDir->text());
    //HBITMAP	toWinHBITMAP ( HBitmapFormat format = NoAlpha ) const

    if (doShootScreens)
    {
        if (op.runstep % printinterval->value() > 0)
            return;

        originalPixmap = QPixmap::grabWidget(tabWidget->widget(2));
        fileName = outdir + fi.baseName() + QString("_q%1.png").arg(op.runstep,5,'d',0,'0');
        originalPixmap.save(fileName, format.toLatin1());

        originalPixmap = QPixmap::grabWidget(tabWidget->widget(3));
        QString type;

        int index = DisplayComboBox->currentIndex();
        if( index > -1 && index < NameList.length())
        {
            QString name = NameList.at(index);
            type = (name + QString("%1.png")).arg(op.runstep,5,'d',0,'0');
        }


        fileName = outdir + fi.baseName() + type;

        originalPixmap.save(fileName, format.toLatin1());
    }
    else
    {
        originalPixmap = QPixmap::grabWidget(tabWidget->currentWidget());

        QString format = "png";
        QString type = ".png";
        QString DT = QDateTime().currentDateTime().toString("yyMMdd-hhmm");//"hh.mm-yy.MM.dd");

        int index = DisplayComboBox->currentIndex();
        if( index > -1 && index < NameList.length())
        {
            QString name = NameList.at(index);
            type = (name + QString("%1.png")).arg(op.runstep,5,'d',0,'0');
        }

        fileName = outdir + fi.baseName() + type;
        originalPixmap.save(fileName, format.toLatin1());


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
                               QString("openLISEM verion %7 (%8) is created wih:\n\n%1\n%2\n%3\n%4\n%5\n")
                               .arg("- Qt cross platform application and UI framework version 4.8.X based on MingW 64bit and CMake (http://qt.nokia.com/).")
                               .arg("- Qwt technical application widgets for Qt (http://qwt.sf.net)")
                               .arg("- Flood source code based on fullSWOF2D (http://www.univ-orleans.fr/mapmo/soft/FullSWOF/)")
                               .arg("- PCRaster map functions: http://pcraster.geo.uu.nl/")
                               .arg("Details can be found at: http://lisem.sourceforge.net")
                               .arg(VERSIONNR)
                               .arg(DATE)
                               );
}
//--------------------------------------------------------------------
void lisemqt::resetAll()
{
    //E_runFileList->clear();
    // no because then never a list to build

    doShootScreens = false;

    DefaultMapnames();
    // Make the default input map list, stringlist

    fillMapnames();
    // make mapList structure according to
    // DEFmaps stringlist that is used to build the map tree interface

    RunFileNames.clear();
    op.runfilename.clear();

    E_MapDir->setText("");
    E_RainfallName->setText("");
    E_SnowmeltName->setText("");
    E_ResultDir->setText("");

    E_RainfallMap->setText("rainfall.map");
    E_InterceptionMap->setText("interception.map");
    E_InfiltrationMap->setText("infiltration.map");
    E_RunoffMap->setText("runoff.map");
    E_RunoffFractionMap->setText("rofraction.map");
    E_ChannelQtotm3Map->setText("chandism3.map");
    E_DetachmentMap->setText("detachment.map");
    E_DepositionMap->setText("deposition.map");
    E_SoillossMap->setText("soilloss.map");
    E_ChanDetachmentMap->setText("chandet.map");
    E_ChanDepositionMap->setText("chandep.map");
    E_WHmaxMap->setText("whmax.map");
    E_FloodlevelMap->setText("floodmax.map");
    E_FloodTimeMap->setText("floodtime.map");
    E_FloodStats->setText("floodstats.txt");
    E_ChannelMaxQ->setText("channelmaxq.map");
    E_ChannelMaxWH->setText("channelmaxhw.map");
    E_FloodFEW->setText("floodstart.map");
    E_FloodmaxVMap->setText("floodmaxv.map");
    E_MainTotals->setText("");
    E_PointResults->setText("");

    E_BeginTime->setText("");
    E_EndTime->setText("");
    E_Timestep->setText("");

    checkBox_OutRunoff->setChecked(false);
    checkBox_OutConc->setChecked(false);
    checkBox_OutWH->setChecked(false);

    checkBox_OutWHC->setChecked(false);
    checkBox_OutWHC->setVisible(false); // <== obsolete

    checkBox_OutTC->setChecked(false);
    checkBox_OutDet->setChecked(false);
    checkBox_OutDep->setChecked(false);
    checkBox_OutV->setChecked(false);
    checkBox_OutInf->setChecked(false);
    checkBox_OutSurfStor->setChecked(false);
    checkBox_OutChanVol->setChecked(false);
    checkBox_OutTiledrain->setChecked(false);
    checkBox_OutHmx->setChecked(false);
    checkBox_OutQf->setChecked(false);
    checkBox_OutVf->setChecked(false);
    checkBox_OutHmxWH->setChecked(false);

    printinterval->setValue(1);

    E_InfiltrationMethod->clear();
    E_InfiltrationMethod->addItem("no Infiltration");
    E_InfiltrationMethod->addItem("SWATRE");
    E_InfiltrationMethod->addItem("Green and Ampt");
    E_InfiltrationMethod->addItem("Smith and Parlange");
    //   E_InfiltrationMethod->addItem("Subtract Ksat");
    E_InfiltrationMethod->setCurrentIndex(0);

    initOP();
    progressBar->setValue(0);

    bool check = false;

    checkDoErosion->setChecked(check);

    checkIncludeChannel->setChecked(check);
    checkChannelInfil->setChecked(check);
    checkChannelBaseflow->setChecked(check);

    checkRoadsystem->setChecked(check);
    checkSnowmelt->setChecked(check);
    checkRainfall->setChecked(true);
    checkHardsurface->setChecked(false);

    checkBox_SedMultiSingle->setChecked(false);
    checkBox_SedMultiMulti->setChecked(false);
    checkBox_SedSingleSingle->setChecked(true);

    E_BLMethod->setValue(0);
    E_SSMethod->setValue(1);

    E_BLMethod->setMaximum(0);
    E_BLMethod->setMinimum(0);
    E_SSMethod->setMaximum(1);
    E_SSMethod->setMinimum(1);

    checkEstimateGrainSizeDistribution->setChecked(false); // if multiclass, estimate from D50 and D90
    checkReadGrainSizeDistribution->setChecked(false); // if multiclass, calculate from user series

    E_BLMethod->setValue(1);
    E_SSMethod->setValue(1);

    E_NumberClasses->setValue(6);
    //E_NumberClassesMaps->setValue(0);
    E_GrainSizes->setText("2,20,50,125,150,500");

    //houses
    checkHouses->setChecked(check);
    checkRaindrum->setChecked(check);
    // flooded areas
    checkChannelFlood->setChecked(check);
    checkInfilCompact->setChecked(check);
    checkInfilGrass->setChecked(check);
    checkInfilCrust->setChecked(check);
    checkImpermeable->setChecked(check);
    //	checkDumphead->setChecked(check);
    checkGeometric->setChecked(true);
    //	checkRunoffPerM->setChecked(check);
    checkPercolation->setChecked(check);

    checkWritePCRnames->setChecked(true);
    checkWritePCRaster->setChecked(check);
    checkWriteCommaDelimited->setChecked(true);
    //checkOutputTimeStep->setChecked(true);
    //checkOutputTimeUser->setChecked(check);
    //checkNoErosionOutlet->setChecked(check);
    //	checkDrainage->setChecked(check);

    checkSeparateOutput->setChecked(check);
    checkWriteSOBEK->setChecked(check);
    SOBEKdatestring->setText("10/01/01");
    //checkInterceptionLAI->setChecked(true);
    tabWidget->setCurrentIndex(0);
    tabWidget_out->setCurrentIndex(1);
    tabWidget_out->setCurrentIndex(0);

    buffergroup->setVisible(false);
    //   buffergroup->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());

    //    sedgroup->setEnabled(checkDoErosion->isChecked());
    //    label_31->setEnabled(checkDoErosion->isChecked());
    //    label_soillosskgha->setEnabled(checkDoErosion->isChecked());

    radioButton_1->setChecked(true); //<= grass interception
    E_CanopyOpeness->setValue(0.45);
    E_StemflowFraction->setValue(0.054);

    //VJ 110706 KE equations
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

    E_floodMinHeight->setValue(0.05);

    E_SSMethod->setValue(1);
    E_SigmaDiffusion->setValue(1);
    E_BLMethod->setValue(1);


}
//--------------------------------------------------------------------
QString lisemqt::findValidDir(QString path, bool up)
{
    if (!QFileInfo(path).exists() || path.isEmpty())
        path = E_MapDir->text();
    //    if (!QFileInfo(path).exists() || path.isEmpty())
    //        path = E_WorkDir->text();
    if (!QFileInfo(path).exists() || path.isEmpty())
    {
        //    path = QFileInfo(op.runfilename).absolutePath();
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

    return (path);
}
//---------------------------------------------------------------
void lisemqt::fontSelect()
{
    // bool ok;
    QFont font = QFontDialog::getFont(0, qApp->font());
    //         &ok, QFont("MS Shell Dlg 2", genfontsize), this);
    //  if (ok) {
    // the user clicked OK and font is set to the font the user selected
    qApp->setFont(font);
    this->setStyleSheet(QString("\
                                QLabel {font: %1pt;} \
                                QGroupBox {font: %1pt;} \
                                QLineEdit {font: %1pt;} \
                                QCheckBox {font: %1pt;} \
                                QRadioButton {font: %1pt;} \
                                QSpeedButton {font: %1pt;} \
                                QDoubleSpinBox {font: %1pt;} \
                                QSpinBox {font: %1pt;} \
                                QComboBox {font: %1pt;} \
                                QTabWidget {font: %1pt;} \
                                QTreeView {font: %1pt;} \
                                QPlainTextEdit {font: %1pt;} \
                                ").arg(genfontsize));
                                //   } else {

                                //   }
}
//---------------------------------------------------------------
void lisemqt::fontDecrease()
{
    genfontsize--;
    genfontsize = std::max(6, genfontsize);

    this->setStyleSheet(QString("\
                                QLabel {font: %1pt;} \
                                QGroupBox {font: %1pt;} \
                                QLineEdit {font: %1pt;} \
                                QCheckBox {font: %1pt;} \
                                QRadioButton {font: %1pt;} \
                                QSpeedButton {font: %1pt;} \
                                QDoubleSpinBox {font: %1pt;} \
                                QSpinBox {font: %1pt;} \
                                QComboBox {font: %1pt;} \
                                QTabWidget {font: %1pt;} \
                                QTreeView {font: %1pt;} \
                                QPlainTextEdit {font: %1pt;} \
                                ").arg(genfontsize));



}
//---------------------------------------------------------------
void lisemqt::fontIncrease()
{
    genfontsize++;
    genfontsize = std::min(18, genfontsize);

    this->setStyleSheet(QString("\
                                QLabel {font: %1pt;} \
                                QGroupBox {font: %1pt;} \
                                QLineEdit {font: %1pt;} \
                                QCheckBox {font: %1pt;} \
                                QRadioButton {font: %1pt;} \
                                QSpeedButton {font: %1pt;} \
                                QDoubleSpinBox {font: %1pt;} \
                                QSpinBox {font: %1pt;} \
                                QComboBox {font: %1pt;} \
                                QTabWidget {font: %1pt;} \
                                QTreeView {font: %1pt;} \
                                QPlainTextEdit {font: %1pt;} \
                                ").arg(genfontsize));
}
