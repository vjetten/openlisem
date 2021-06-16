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
    showMaximized();

    darkLISEM = false;//true;

    nrUserCores->setMaximum(omp_get_max_threads());

    helpbox = new QDialog();
    helpbox->resize(1080, 600);
    helpbox->setWindowTitle("option help");
    helpLayout = new QHBoxLayout(helpbox);
    helptxt = new QTextEdit();
    helpLayout->addWidget(helptxt);

    //checkAddDatetime->setVisible(false);

    //RunFileNames.clear();
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

    W = nullptr;
    // initalize pointer to the world, created when run button is pushed

    SetStyleUI();
    // do some style things

    doNewPlot = false;
    if(doNewPlot)
        setupNewPlot();
    else
        setupPlot();
    // set up the discharge graphs

    setupMapPlot();
    // set up the raster map drawing

//    QSplitter *splitter = new QSplitter(tabWidget->widget(2));
//    splitter->addWidget(tabWidget_out);
//    splitter->addWidget(widgetMB);

    Ui_lisemqtClass::statusBar->addWidget(progressBar, 1);
    // put the progress bar into the statusbar

    tabWidgetOptions->setCurrentIndex(0);
    tabWidget_OutputMaps->setCurrentIndex(0);

    doBatchmode = doBatch; // save in iface for later use
    batchRunname = runname;
    doCheckRainfall(true); // ???????? why here

    if(doBatch)
    {
        runfilelist.clear();
        runfilelist << runname;//.split(';');

//        for (int n = 0; n < runfilelist.size(); n++) {
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
        //triggered in on_E_runFileList_currentIndexChanged(int)
    }
}
//--------------------------------------------------------------------
lisemqt::~lisemqt()
{
    if (!doBatchmode)
        StorePath();
}
//--------------------------------------------------------------------
// NAMING convention void on_<widget name="">_<signal name="">(<signal parameters="">)
// works automatically. if included here may be executed twice!!! not sure...
void lisemqt::SetConnections()
{
    connect(checkPesticides, SIGNAL(toggled(bool)), this, SLOT(doCheckPesticides(bool)));

    connect(toolButton_fileOpen, SIGNAL(clicked()), this, SLOT(openRunFile()));
    connect(toolButton_deleteRun, SIGNAL(clicked()), this, SLOT(deleteRunFileList()));
    connect(toolButton_MapDir, SIGNAL(clicked()), this, SLOT(setMapDir()));

    connect(treeView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(openMapname(QModelIndex)));
    connect(MapNameModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(editMapname(QModelIndex, QModelIndex)));
    connect(toolButton_ResultDir, SIGNAL(clicked()), this, SLOT(setResultDir()));

    connect(checkIncludeChannel, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));
    connect(checkOverlandFlow1D, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));
    connect(checkOverlandFlow2Dkindyn, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));
    connect(checkOverlandFlow2Ddyn, SIGNAL(toggled(bool)), this, SLOT(setFloodTab(bool)));

    connect(spinBoxPointtoShow,SIGNAL(valueChanged(int)),this,SLOT(onOutletChanged(int)));

    connect(checkFormatGtiff, SIGNAL(toggled(bool)), this, SLOT(setFormatMaps(bool)));

    connect(E_BulkDens2,SIGNAL(editingFinished()),this, SLOT(updateBulkDens()));
    connect(E_BulkDens,SIGNAL(editingFinished()),this, SLOT(updateBulkDens2()));

}
void lisemqt::updateBulkDens()
{
    QString txt = E_BulkDens2->text();
    E_BulkDens->setText(txt);
}
void lisemqt::updateBulkDens2()
{
    QString txt = E_BulkDens->text();
    E_BulkDens2->setText(txt);
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
void lisemqt::resizeEvent(QResizeEvent* event)
{
    QMainWindow::resizeEvent(event);
    groupBox_drawMap->setEnabled(tabWidget_out->currentIndex() == 1);

    int h = QApplication::desktop()->height(); //event->size().height()
    if (h > 800)
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
}
//--------------------------------------------------------------------
// bad programming, checkboxes as radiobuttons, but needed to be square buttons!
void lisemqt::on_checkOverlandFlow1D_clicked()
{
    tabWidgetOptions->setTabEnabled(3, true);
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
     //   if (!op.ComboSymColor.at(i))
        {
            ComboMinSpinBox2->setEnabled(true);
            if (op.userMaxV.at(i) == 0)
                op.userMinV.replace(i, d);

            if (op.userMaxV.at(i) > 0 && d < op.userMaxV.at(i))
                op.userMinV.replace(i, d);

            if(op.userMaxV.at(i) > 0 && d >= op.userMaxV.at(i))
                ComboMinSpinBox2->setValue(op.userMinV.at(i));
        }
     //   else
        //    ComboMinSpinBox2->setEnabled(false);

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

  //  ComboMinSpinBox2->setEnabled(!op.ComboSymColor.at(i));
    ComboMaxSpinBox2->setValue(op.userMaxV.at(i));
    ComboMinSpinBox2->setValue(op.userMinV.at(i));
    this->showMap();
}
/*
void lisemqt::setSedimentText(int i, int j, int k)
{
    // i = TC, j= river or surface, k = BL or SS
    if (j == 0) {
        if ( k == 0) {
            if (i==1) label_RTBL->setText("Channel Bedload TC Van Rijn (simplified)");
            if (i==2) label_RTBL->setText("Channel Bedload TC Van Rijn (full)");
            if (i==3) label_RTBL->setText("Channel Bedload TC Wu, Wang & Jia (multiclass)");
        } else {
            if (i==0) label_RTSS->setText("Channel Suspended TC Govers");
            if (i==1) label_RTSS->setText("Channel Suspended TC Van Rijn (simplified)");
            if (i==2) label_RTSS->setText("Channel Suspended TC Van Rijn (full)");
            if (i==3) label_RTSS->setText("Channel Suspended TC Wu, Wang & Jia (multiclass)");
        }
    } else {
        if ( k == 0) {
            if (i==1) label_STBL->setText("Overland flow Bedload TC Van Rijn (simplified)");
            if (i==2) label_STBL->setText("Overland flow Bedload TC Van Rijn (full)");
            if (i==3) label_STBL->setText("Overland flow Bedload TC Wu, Wang & Jia (multiclass)");
        } else {
            if (i==0) label_STSS->setText("Overland flow Suspended TC Govers");
            if (i==1) label_STSS->setText("Overland flow Suspended TC Van Rijn (simplified)");
            if (i==2) label_STSS->setText("Overland flow Suspended TC Van Rijn (full)");
            if (i==3) label_STSS->setText("Overland flow Suspended TC Wu, Wang & Jia (multiclass)");
        }
    }
}
//--------------------------------------------------------------------
void lisemqt::on_E_RBLMethod_valueChanged(int i)
{
  //  setSedimentText(i, 0, 0);
}
//--------------------------------------------------------------------
void lisemqt::on_E_RSSMethod_valueChanged(int i)
{
    //setSedimentText(i, 0, 1);
}
//--------------------------------------------------------------------
void lisemqt::on_E_BLMethod_valueChanged(int i)
{
//    setSedimentText(i, 1, 0);
}
//--------------------------------------------------------------------
void lisemqt::on_E_SSMethod_valueChanged(int i)
{
  //  setSedimentText(i, 1, 1);
}
*/
//--------------------------------------------------------------------
void lisemqt::on_checkSed2Phase_toggled(bool v)
{
      sedbox3->setEnabled(v);
      E_RBLMethod->setEnabled(v);
      E_BLMethod->setEnabled(v);
}
//--------------------------------------------------------------------

//void lisemqt::on_checkSedMultiGrain_toggled(bool v)
//{
//    if(v) {
//        checkSed2Phase->setChecked(false);
//     //   sedbox1->setEnabled(true);
//        sedbox2->setEnabled(true);
//        sedbox3->setEnabled(true);
//        E_RBLMethod->setValue(3);
//        E_RSSMethod->setValue(4);
//        E_BLMethod->setValue(3);
//        E_SSMethod->setValue(4);
//        E_RBLMethod->setEnabled(false);
//        E_RSSMethod->setEnabled(false);
//        E_BLMethod->setEnabled(false);
//        E_SSMethod->setEnabled(false);

//    }
//}
//--------------------------------------------------------------------
void lisemqt::setFloodTab(bool yes)
{
    yes = true;
    if (checkOverlandFlow2Dkindyn->isChecked() && !checkIncludeChannel->isChecked()) {
        yes = false;
        QMessageBox::warning(this,"openLISEM",QString("The combination of 1D overland flow and 2D flood can only be used with a channel activated."));
        checkOverlandFlow1D->setChecked(true);
    }
    if (checkOverlandFlow1D->isChecked()) {
        yes = false;
    }

    tabWidgetOptions->setTabEnabled(3, yes);
    //frame_diffwave->setEnabled(checkOverlandFlow2D->isChecked());
    FloodParams->setEnabled(yes);

    outputMapsFlood->setEnabled(yes);
    label_floodVolmm->setEnabled(yes);
    label_107->setEnabled(yes);

    if (checkOverlandFlow2Ddyn->isChecked() || checkOverlandFlow2Dkindyn->isChecked()) {
        label_107->setText(QString("Flood (mm),h>%1)").arg(E_floodMinHeight->value()*1000));
        label_40->setText(QString("Runoff (mm),h<%1)").arg(E_floodMinHeight->value()*1000));
    }
    else
    {
        label_107->setText("Flood mm");
        label_40->setText("Runoff mm");
    }

}
//--------------------------------------------------------------------
void lisemqt::setErosionTab()
{
    tabWidgetOptions->setTabEnabled(4, checkDoErosion->isChecked());

    int i1 = E_RBLMethod->value();
    int i2 = E_RSSMethod->value();
    int i3 = E_BLMethod->value();
    int i4 = E_SSMethod->value();

    E_RBLMethod->setValue(0);
    E_RSSMethod->setValue(0);
    E_BLMethod->setValue(0);
    E_SSMethod->setValue(0);
    E_RBLMethod->setValue(i1);
    E_RSSMethod->setValue(i2);
    E_BLMethod->setValue(i3);
    E_SSMethod->setValue(i4);

    //  yes = checkDoErosion->isChecked();
    outputMapsSediment->setEnabled(checkDoErosion->isChecked());
    checkBox_OutConc->setEnabled(checkDoErosion->isChecked());
    checkBox_OutDet->setEnabled(checkDoErosion->isChecked());
    checkBox_OutDep->setEnabled(checkDoErosion->isChecked());
    checkBox_OutSL->setEnabled(checkDoErosion->isChecked());
    checkBox_OutSed->setEnabled(checkDoErosion->isChecked());
    checkBox_OutTC->setEnabled(checkDoErosion->isChecked());
    checkBox_OutSedSS->setEnabled(checkDoErosion->isChecked() && checkSed2Phase->isChecked());
    checkBox_OutSedBL->setEnabled(checkDoErosion->isChecked() && checkSed2Phase->isChecked());

    checkBoxComboMaps2->setEnabled(checkDoErosion->isChecked());
    ComboMinSpinBox2->setEnabled(checkDoErosion->isChecked());
    ComboMaxSpinBox2->setEnabled(checkDoErosion->isChecked());
    DisplayComboBox2->setEnabled(checkDoErosion->isChecked());

    if (!checkDoErosion->isChecked())
    {
        int dig = E_DigitsOut->value();
        label_MBs->setText(QString::number(0,'e',dig));
        label_splashdet->setText(QString::number(0,'f',dig));
        label_flowdet->setText(QString::number(0,'f',dig));
        label_sedvol->setText(QString::number(0,'f',dig));
        label_dep->setText(QString::number(0,'f',dig));
        label_detch->setText(QString::number(0,'f',dig));
        label_depch->setText(QString::number(0,'f',dig));
        label_sedvolch->setText(QString::number(0,'f',dig));
        //        label_flooddet->setText(QString::number(0,'f',dig));
        //        label_flooddep->setText(QString::number(0,'f',dig));
        //        label_floodsed->setText(QString::number(0,'f',dig));
        label_soilloss->setText(QString::number(0,'f',dig));
        label_soillosskgha->setText(QString::number(0,'f',dig));
        label_SDR->setText(QString::number(0,'f',dig));
        label_soillosssub->setText(QString::number(0,'f',dig));
    }
    sedgroup->setEnabled(checkDoErosion->isChecked());

    label_soillosskgha->setEnabled(checkDoErosion->isChecked());
    label_soilloss->setEnabled(checkDoErosion->isChecked());
    label_SDR->setEnabled(checkDoErosion->isChecked());

}
//--------------------------------------------------------------------
void lisemqt::setRunoffTab(bool)
{
    tabWidgetOptions->setTabEnabled(3, true);//checkChannelFlood->isChecked() || checkOverlandFlow2D->isChecked()  || checkOverlandFlow2Ddyn->isChecked());
}

//--------------------------------------------------------------------
void lisemqt::setWriteOutputSOBEK(bool doit)
{
    //   checkWriteSOBEK->setChecked(!doit);
    //checkWriteCommaDelimited->setChecked(!doit);
    //checkWritePCRaster->setChecked(!doit);
}
//--------------------------------------------------------------------
void lisemqt::setWriteOutputCSV(bool doit)
{
   // checkWriteSOBEK->setChecked(!doit);
    //   checkWriteCommaDelimited->setChecked(!doit);
   // checkWritePCRaster->setChecked(!doit);
}

void lisemqt::setOutputScreen()
{
  if (W) {
    W->noInterface = !W->noInterface;
    showAllAct->setChecked(!W->noInterface);
  }
}

//--------------------------------------------------------------------
void lisemqt::setWriteOutputPCR(bool /* doit */)
{
    //    if (checkWriteSOBEK->isChecked())
    //    {
    //        //checkWriteSOBEK->setChecked(false);
    //        checkWriteCommaDelimited->setChecked(false);
    //        checkWritePCRaster->setChecked(false);
    //        //checkSeparateOutput->setChecked(true);
    //    }
    //    else
//    if (checkWritePCRaster->isChecked())
//    {
//        //checkWriteSOBEK->setChecked(false);
//        //checkWriteCommaDelimited->setChecked(false);
//    }
//    else
//        if (checkWriteCommaDelimited->isChecked())
//        {
//            //checkWriteSOBEK->setChecked(false);
//            checkWritePCRaster->setChecked(false);
//        }
}
//--------------------------------------------------------------------
void lisemqt::SetToolBar()
{
    toolBar->setIconSize(QSize(32,32));

    //r
    restartAct = new QAction(QIcon(":/2X/reset.png"), "&Reset...", this);
    connect(restartAct, SIGNAL(triggered()), this, SLOT(resetAll()));
    toolBar->addAction(restartAct);
    toolBar->addSeparator();

    openAct = new QAction(QIcon(":/2X/Folder-Open-icon.png"), "&Open a run file...", this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip("Open a run file");
    connect(openAct, SIGNAL(triggered()), this, SLOT(openRunFile()));
    toolBar->addAction(openAct);

    saveAct = new QAction(QIcon(":/2X/filesave2X.png"), "&Save the run file...", this);
 //   saveAct = new QAction(QIcon(":/2X/savefile.png"), "&Save the run file...", this);
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

    fontIncreaseAct = new QAction(QIcon(":/2X/fontbigger2X.png"), "&Increase font size", this);
    connect(fontIncreaseAct, SIGNAL(triggered()), this, SLOT(fontIncrease()));
    toolBar->addAction(fontIncreaseAct);
    fontDecreaseAct = new QAction(QIcon(":/2X/fontsmaller2X.png"), "&Decrease font size", this);
    connect(fontDecreaseAct, SIGNAL(triggered()), this, SLOT(fontDecrease()));
    toolBar->addAction(fontDecreaseAct);

    toolBar->addSeparator();
    showAllAct = new QAction(QIcon(":/2X/noscreen.png"), "&no output to screen", this);
    showAllAct->setCheckable(true);
    connect(showAllAct, SIGNAL(triggered()), this, SLOT(setOutputScreen()));
    toolBar->addAction(showAllAct);

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
        toolBar->addSeparator();

    connect(checkMapBuildings, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
    connect(checkMapRoads, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
  //  connect(checkMapChannels, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));
    connect(checkMapChannels, SIGNAL(clicked(bool)), this, SLOT(hideChannelVector(bool)));
    connect(checkMapImage, SIGNAL(clicked(bool)), this, SLOT(showMapb(bool)));

    connect(ComboMaxSpinBox,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMinSpinBox,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMaxSpinBox2,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));
    connect(ComboMinSpinBox2,SIGNAL(valueChanged(double)),this,SLOT(showMapd(double)));

    connect(transparency, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlpha(int)));
  //  connect(transparencyChannel, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlphaChannel(int)));
    connect(transparencyRoad, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlphaRoad(int)));
    connect(transparencyHouse, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlphaHouse(int)));
    connect(transparencyImage, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlphaBarrier(int))); // barrier is use for image
    connect(transparencyMap, SIGNAL(sliderMoved(int)), this, SLOT(ssetAlphaMap(int)));
    //connect(showRiverSize, SIGNAL(valueChanged(int)),this,SLOT(ssetAlphaChannel(int)));
    connect(spinChannelSize, SIGNAL(valueChanged(int)),this,SLOT(ssetAlphaChannel(int)));
    connect(spinCulvertSize, SIGNAL(valueChanged(int)),this,SLOT(ssetAlphaChannelOutlet(int)));
    connect(toolShowMapDisplay, SIGNAL(pressed()),this,SLOT(showMapSettings()));
}
//---------------------------------------------------------------------------
void lisemqt::showMapSettings()
{
    if(groupBox_drawMap->isVisible())
        groupBox_drawMap->setVisible(false);
    else
        groupBox_drawMap->setVisible(true);
}
//---------------------------------------------------------------------------
int lisemqt::SetStyleUISize()
{
    QRect rect = QGuiApplication::primaryScreen()->availableGeometry();
    //qDebug() << rect << QGuiApplication::primaryScreen()->availableVirtualGeometry();
    int _H = rect.height();
    int disp = 1;
    if (abs(_H-800) > abs(_H-1080)) disp = 1;
    if (abs(_H-1080) > abs(_H-1200)) disp = 2;
    if (abs(_H-1200) > abs(_H-1440)) disp = 3;
      if (abs(_H-1440) > abs(_H-1600)) disp = 4;
    //qDebug() << disp << _H << abs(_H-800) << abs(_H-1080) << abs(_H-1200) << abs(_H-1440) << genfontsize;

    tabWidgetOptions->tabBar()->setExpanding(true);
  //  tabWidget_out->setIconSize(QSize(24, 24));

    // do a bit of size tweaking for large displays
    QSize iSize = QSize(16,16);
    if (disp == 1) {
        tabWidget_out->setIconSize(QSize(24, 24));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 48px; width: 32px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 16px 16px}"));
        iSize = QSize(16,16);

    }
    if (disp == 2 || disp == 3) {
        tabWidget_out->setIconSize(QSize(32, 32));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 64px; width: 48px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 24px 24px}"));
        iSize = QSize(24,24);
    }
    if (disp > 3) {
        tabWidget_out->setIconSize(QSize(48, 48));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 96px; width: 64px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 32px 32px}"));
        iSize = QSize(32,32);
    }

   // this->setStyleSheet(QString("QLabel::pixmap {height: 16px; width: 16px}"));

    toolBar->setIconSize(iSize);
    toolBar_2->setIconSize(iSize);

    return disp;
}
//---------------------------------------------------------------------------
/// make some labels yellow
void lisemqt::SetStyleUI()
{
    trayIcon = new QSystemTrayIcon(this);
    trayIcon->setIcon(QIcon(":/openLisem.ico"));
    trayIcon->show();
    tabWidgetOptions->tabBar()->setExpanding(true);

    QFont font = qApp->font();
    genfontsize=std::max(10,font.pointSize());
    int x = SetStyleUISize();
    genfontsize += x;
    font.setPointSize(genfontsize);
//    qApp->setFont(font);

    toolBar_2->setMovable( false);
    toolBar->setMovable( false);

    QString flat("QToolButton { background-color: white; border: none; }");

    groupBoxInput->setStyleSheet("QGroupBox::title{color: #4477aa;}");
    groupBoxOutput->setStyleSheet("QGroupBox::title{color: #4477aa;}");
    checkDoErosion->setStyleSheet("QCheckBox {color: #4477aa;}");
    checkSed2Phase->setStyleSheet("QCheckBox {color: #4477aa;}");
    label_55->setStyleSheet("QLabel {color: #4477aa;}");
    label_9->setStyleSheet("QLabel {color: #4477aa;}");
    label_10->setStyleSheet("QLabel {color: #4477aa;}");
    label_11->setStyleSheet("QLabel {color: #4477aa;}");

    // interface elements that are not visible for now
    frameSpare->setVisible(false);
    tabWidgetOptions->removeTab(8);
    tabWidgetOptions->removeTab(7);
    //frameNumerical->setVisible(false);

    //tabWidget_erosion->setCurrentIndex(0);
    if (darkLISEM) {
        QPalette darkPalette;
        darkPalette.setColor(QPalette::Window,QColor(53,53,53));
        darkPalette.setColor(QPalette::WindowText,Qt::white);
        darkPalette.setColor(QPalette::ButtonText,Qt::white);
        darkPalette.setColor(QPalette::HighlightedText,Qt::white);
        darkPalette.setColor(QPalette::Text,QColor("#ffffaa"));
        darkPalette.setColor(QPalette::Disabled,QPalette::WindowText,QColor(164,164,164));
        darkPalette.setColor(QPalette::Base,QColor(96,96,96));
        darkPalette.setColor(QPalette::AlternateBase,QColor(66,66,66));
        darkPalette.setColor(QPalette::ToolTipBase,Qt::white);
        darkPalette.setColor(QPalette::ToolTipText,Qt::white);

        darkPalette.setColor(QPalette::Disabled,QPalette::Text,QColor(164,164,164));
        darkPalette.setColor(QPalette::Dark,QColor(35,35,35));
        darkPalette.setColor(QPalette::Shadow,QColor(20,20,20));
        darkPalette.setColor(QPalette::Button,QColor(53,53,53));
        darkPalette.setColor(QPalette::Disabled,QPalette::ButtonText,QColor(164,164,164));
        darkPalette.setColor(QPalette::BrightText,Qt::red);
        darkPalette.setColor(QPalette::Link,QColor(42,130,218));
        darkPalette.setColor(QPalette::Highlight,QColor(42,130,218));
        darkPalette.setColor(QPalette::Disabled,QPalette::Highlight,QColor(80,80,80));
        darkPalette.setColor(QPalette::Disabled,QPalette::HighlightedText,QColor(127,127,127));

        qApp->setPalette(darkPalette);
    }

    int w = 80, h = 15;
    label_dx->setMinimumSize(w,h);
    label_area->setMinimumSize(w,h);
    label_time->setMinimumSize(w,h);
    label_endtime->setMinimumSize(w,h);
    label_raintot->setMinimumSize(w,h);
    label_watervoltot->setMinimumSize(w,h);
    label_stormdraintot->setMinimumSize(w,h);
    label_qtot->setMinimumSize(w,h);
    label_infiltot->setMinimumSize(w,h);
    label_surfstor->setMinimumSize(w,h);
    label_interctot->setMinimumSize(w,h);
    //label_qtotm3->setMinimumSize(w,h);
    label_qpeaktime->setMinimumSize(w,h);
    label_ppeaktime->setMinimumSize(w,h);
    label_QPfrac->setMinimumSize(w,h);
    //label_discharge->setMinimumSize(w,h);
    label_floodVolmm->setMinimumSize(w,h);
    label_watervolchannel->setMinimumSize(w,h);
    //   label_litterstore->setMinimumSize(w,h);
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

    label_MBs->setMinimumSize(w,h);
    label_MB->setMinimumSize(w,h);

    QString ly = "#ffff99";
    if (!darkLISEM) {
        label_dx->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_area->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_time->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_endtime->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_raintot->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_watervoltot->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_stormdraintot->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_qtot->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_infiltot->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_surfstor->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_interctot->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        //label_qtotm3->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_qpeaktime->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_ppeaktime->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_QPfrac->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        //label_discharge->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_floodVolmm->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_watervolchannel->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        //   label_litterstore->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_baseflowtot->setStyleSheet(QString("* { background-color: %1 }").arg(ly));

        label_qtotm3sub->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_dischargesub->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_qpeaksub->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_soillosssub->setStyleSheet(QString("* { background-color: %1 }").arg(ly));

        label_splashdet->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_flowdet->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_sedvol->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_dep->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_detch->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_depch->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_sedvolch->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_soilloss->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_soillosskgha->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
        label_SDR->setStyleSheet(QString("* { background-color: %1 }").arg(ly));
    }
    //Grouped Buttons become mututally exclusive
    GroupMapDisplay.addButton(checkBoxComboMaps, 1);
    GroupMapDisplay.addButton(checkBoxComboMaps2, 2);

 //   GroupImpermeable.addButton(checkImpermeable,1);
 //   GroupImpermeable.addButton(checkPercolation,2);

    GroupRunoff.addButton(checkOverlandFlow1D,1);
  //  GroupRunoff.addButton(checkOverlandFlow2D,2);
    GroupRunoff.addButton(checkOverlandFlow2Ddyn,3);
    GroupRunoff.addButton(checkOverlandFlow2Dkindyn,4);


    if (checkOverlandFlow2Ddyn->isChecked()) {
        label_107->setText(QString("Flood (mm),h>%1)").arg(E_floodMinHeight->value()*1000));
        label_40->setText(QString("Runoff (mm),h<%1)").arg(E_floodMinHeight->value()*1000));

    } else {
        label_107->setText("Flood mm");
        label_40->setText("Runoff mm");
    }
    bool yes = !checkOverlandFlow1D->isChecked();
    label_floodVolmm->setEnabled(yes);
    label_107->setEnabled(yes);

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
void lisemqt::on_E_floodMinHeight_valueChanged(double)
{
    label_107->setText(QString("Flood (mm),h>%1)").arg(E_floodMinHeight->value()*1000));
    label_40->setText(QString("Runoff (mm),h<%1)").arg(E_floodMinHeight->value()*1000));
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
                                        QString("Select the SWATRE profile definition file"),
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
void lisemqt::on_toolButton_satImageName_clicked()
{
    QString path;

    satImageFileDir = findValidDir(satImageFileDir, false);

    path = QFileDialog::getOpenFileName(this,
                                        QString("Select background satellite image file"),
                                        satImageFileDir,"GeoTiff (*.tif)");
    if(!path.isEmpty())
    {
        QFileInfo fi(path);
        satImageFileName = fi.fileName();
        satImageFileDir = CheckDir(fi.absolutePath());//Dir().path());
        E_satImageName->setText( satImageFileDir + satImageFileName );
    }
}
//--------------------------------------------------------------------
//void lisemqt::on_toolButton_SnowmeltName_clicked()
//{
//    QString path;

//    SnowmeltFileDir = findValidDir(SnowmeltFileDir, false);

//    path = QFileDialog::getOpenFileName(this,
//                                        QString("Select snow melt file"),
//                                        SnowmeltFileDir);
//    if(!path.isEmpty())
//    {
//        QFileInfo fi(path);
//        SnowmeltFileName = fi.fileName();
//        SnowmeltFileDir = CheckDir(fi.absolutePath());//Dir().path());
//        E_SnowmeltName->setText( SnowmeltFileDir + SnowmeltFileName );
//    }
//}
//--------------------------------------------------------------------
//void lisemqt::on_toolButton_SnowmeltShow_clicked()
//{
//    QFile file(SnowmeltFileDir + SnowmeltFileName);
//    if (!file.open(QFile::ReadOnly | QFile::Text))
//    {
//        QMessageBox::warning(this,"openLISEM",
//                             QString("Cannot read file %1:\n%2.")
//                             .arg(SnowmeltFileDir + SnowmeltFileName)
//                             .arg(file.errorString()));
//        return;
//    }

//    QTextStream in(&file);

//    QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
//    view->setWindowTitle(SnowmeltFileName);
//    view->setMinimumWidth(400);
//    view->setMinimumHeight(500);
//    view->setAttribute(Qt::WA_DeleteOnClose);
//    view->show();

//    file.close();
//}
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
    if (!fileName.isEmpty()) {
        updateModelData();
        savefile(fileName);
    }

}
//--------------------------------------------------------------------
void lisemqt::saveRunFile()
{
    updateModelData();
    // change runfile strings with current interface options
    savefile(op.runfilename);
}
//--------------------------------------------------------------------
void lisemqt::savefile(QString name)
{
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

    QString fileName = screenShotDir + fi.baseName();
    QString number = QString("-%1").arg(op.runstep,5,'d',0,'0');
    QString name;

    if (doShootScreens)
    {
        if (op.runstep % printinterval->value() > 0)
            return;

        tabWidget_out->setCurrentIndex(0);
        originalPixmap = tabWidget->widget(2)->grab(); //QPixmap::grabWidget(tabWidget->widget(2));
        // originalPixmap = QPixmap::grabWidget(tabWidget_out->widget(0));
        fileName = screenShotDir + fi.baseName()+ "_Q" + number + ".png";

        originalPixmap.save(fileName, format.toLatin1());

        tabWidget_out->setCurrentIndex(1);
        originalPixmap = tabWidget->widget(2)->grab(); //QPixmap::grabWidget(tabWidget->widget(2));
        QString name = "";
        if (checkBoxComboMaps->isChecked()) {
            int index = DisplayComboBox->currentIndex();
            if( index > -1 && index < NameList.length())
                name = NameList.at(index);
        } else if (checkBoxComboMaps2->isChecked()) {
            int index = DisplayComboBox2->currentIndex()+DisplayComboBox->count();
         //   qDebug() << index;
            if( index > -1 && index < NameList.length())
                name = NameList.at(index);
        }

        fileName = screenShotDir + fi.baseName()+ name + number  + ".png";
        originalPixmap.save(fileName, format.toLatin1());
    }
    else
    {
        originalPixmap = tabWidget->currentWidget()->grab(); //QPixmap::grabWidget(tabWidget->currentWidget());

        if (tabWidget->currentIndex() == 1)
            name = "inputmaps";
        if (tabWidget->currentIndex() == 2) // output
        {
            if (tabWidget_out->currentIndex() == 0) {
                name = "_Q";
            }
            if (tabWidget_out->currentIndex() == 1) {
                if (checkBoxComboMaps->isChecked()) {
                    int index = DisplayComboBox->currentIndex();
                    if( index > -1 && index < NameList.length())
                        name = NameList.at(index);
                } else if (checkBoxComboMaps2->isChecked()) {
                    int index = DisplayComboBox2->currentIndex()+DisplayComboBox->count();
                    if( index > -1 && index < NameList.length())
                        name = NameList.at(index);
                }
            }
        }
        fileName = screenShotDir + fi.baseName()+ name + number  + ".png";
                 //qDebug() << fileName;
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
                               QString("openLISEM verion %8 (%9) is created wih:\n\n%1\n%2\n%3\n%4\n%5\n%6\n%7\n%8")
                               .arg("- MSYS2 with MingW64, Qt and CMake (http://qt.nokia.com/).")
                               .arg("- Qwt technical application widgets for Qt (http://qwt.sf.net)")
                               .arg("- Flood source code partly based on fullSWOF2D (http://www.univ-orleans.fr/mapmo/soft/FullSWOF/)")
                               .arg("- Using openMP for parallel processing (https://www.openmp.org/)")
                               .arg("- Using GDAL for map handling (https://gdal.org/)")
                               .arg("- PCRaster lib map functions: http://pcraster.geo.uu.nl/")
                               .arg("Details can be found at: http://lisem.sourceforge.net")
                               .arg("This software is made available under GNU CPLv3.0")
                               .arg(VERSIONNR)
                               .arg(DATE)
                               );
}
//--------------------------------------------------------------------
void lisemqt::resetTabCalibration()
{
    //calibration
    E_CalibrateKsat->setValue(1.0);
    E_CalibrateN->setValue(1.0);
    E_CalibrateTheta->setValue(1.0);
    E_CalibratePsi->setValue(1.0);
    E_CalibrateChKsat->setValue(1.0);
    E_CalibrateChN->setValue(1.0);
    E_CalibrateAS->setValue(1.0);
    E_CalibrateCOH->setValue(1.0);
    E_CalibrateGS->setValue(1.0);
    E_CalibrateCHCOH->setValue(1.0);
}


void lisemqt::resetTabInterception()
{
    radioButton_1->setChecked(true); //<= crops interception
    E_CanopyOpeness->setValue(0.45);
    //    E_StemflowFraction->setValue(0.054);
    checkIncludeLitter->setChecked(false);
    E_LitterSmax->setValue(1.0);
}

void lisemqt::resetTabInfiltration()
{
    //infiltration
    checkInfilCompact->setChecked(false);
    checkInfilCrust->setChecked(false);
    checkInfil2layer->setChecked(false);
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
    E_pitValue->setValue(1.0);
    checkSWOFWatersheds->setChecked(false);
    rb_swof2->setChecked(true);
    E_TimestepMinFlood->setValue(0.2);
    E_courantFactor->setValue(0.2);
}
//--------------------------------------------------------------------
void lisemqt::resetTabSediment()
{
    //sediment
    checkSed2Phase->setChecked(false);
//    checkSedMultiGrain->setChecked(false);

    E_RBLMethod->setValue(1);
    E_RSSMethod->setValue(0);
    E_BLMethod->setValue(1);
    E_SSMethod->setValue(0);

    E_SigmaDiffusion->setValue(0.5);
//    E_RSigmaDiffusion->setValue(0.5);

//    checkEstimateGrainSizeDistribution->setChecked(false); // if multiclass, estimate from D50 and D90
//    checkReadGrainSizeDistribution->setChecked(false); // if multiclass, calculate from user series

//    E_NumberClasses->setValue(6);
//    E_GrainSizes->setText("2;20;50;125;150;500");

    checkDiffusion->setChecked(true);
    //checkRDiffusion->setChecked(true);

    checkSed2Phase->setChecked(false);

}
//--------------------------------------------------------------------
void lisemqt::resetTabErosion()
{
    //checkLimitTC->setChecked(false);
    checkSedtrap->setChecked(false);
    E_EfficiencyDET->setValue(1);
    checkMaterialDepth->setChecked(false);
    E_DepositedCohesion->setValue(0.5);
    E_BulkDens->setText("1500.00");
    E_BulkDens2->setText("1500.00");

    E_SplashDelibery->setValue(0.1);
    checkInfilGrass->setChecked(false);
    E_GrassStripN->setText("0.2");


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
    //checkFloodSedimentInterpolation->setChecked(true);

    E_SedTrapN->setText("0.8");
    E_GrassStripN->setText("0.2");

    resetTabSediment();
}

void lisemqt::resetTabAdvanced()
{
    E_mixingFactor->setValue(2.0);
    E_runoffPartitioning->setValue(1.0);
    E_FloodMaxIter->setValue(200);
    E_FloodReconstruction->setValue(4);  //HLL2 etc
    E_FloodFluxLimiter->setValue(1);     //minmod etc
    E_courantFactorSed->setValue(0.2);
    checkVariableTimestep->setChecked(false);
    checkHeun->setChecked(false);
    checkTimeavgV->setChecked(true);
    checkLinkedList->setChecked(false);
    checkErosionInsideLoop->setChecked(true);
    checkGravityToChannel->setChecked(false);
    checkKinWaveChannel->setChecked(true);
    E_ChannelKinWaveDt->setValue(60.0);
    checkKinWaveChannelAvg->setChecked(false);
    nrUserCores->setValue(0);
}

void lisemqt::resetAll()
{
    W = nullptr;
    MapNameModel = nullptr;
   // HPlot = nullptr;
   // MPlot = nullptr;

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
    E_FloodStats->setText("floodstats.txt");
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

    E_InfiltrationMethod->clear();
    E_InfiltrationMethod->addItem("no Infiltration");
    E_InfiltrationMethod->addItem("SWATRE");
    E_InfiltrationMethod->addItem("Green and Ampt");
    E_InfiltrationMethod->addItem("Smith and Parlange");
    E_InfiltrationMethod->setCurrentIndex(2);

    initOP();

    progressBar->setValue(0);

    //checkSnowmelt->setChecked(false);
    //checkRainfall->setChecked(true);

    //main
    checkOverlandFlow1D->setChecked(false);
  //  checkOverlandFlow2D->setChecked(false);
    checkOverlandFlow2Ddyn->setChecked(true);
    checkOverlandFlow2Dkindyn->setChecked(false);
    //frame_diffwave->setEnabled(checkOverlandFlow2D->isChecked());
    //frame_dynwave->setEnabled(checkOverlandFlow2Ddyn->isChecked());
    //groupBox_coupling->setEnabled(!checkOverlandFlow2Ddyn->isChecked());

    checkDoErosion->setChecked(false);
    //checkAdvancedSediment->setChecked(false);

    checkIncludeChannel->setChecked(true);
    //checkChannelFlood->setChecked(true);
    checkChannelInfil->setChecked(false);
    checkChannelBaseflow->setChecked(false);
    checkChannelInflow->setChecked(false);
    checkChannelAdjustCHW->setChecked(false);

    checkRoadsystem->setChecked(false);
    checkHouses->setChecked(false);
    checkHardsurface->setChecked(false);
    checkRaindrum->setChecked(false);
    checkStormDrains->setChecked(false);


    resetTabInterception();

    resetTabInfiltration();

    //flow
    resetTabFlow();

    //erosion
    resetTabErosion();

    //calibration
    resetTabCalibration();
    resetTabAdvanced();


    checkVariableTimestep->setChecked(false);
    checkHeun->setChecked(false);
    //checkSWOFomp->setChecked(true);
    //checkMuscl->setChecked(true);
    check2DDiagonalFlow->setChecked(true);
    E_courantFactor->setValue(0.2);
    //   if (p1.compare("Flooding courant factor diffusive")==0)        E_courantFactorSed->setValue(valc);

    E_floodMinHeight->setValue(0.05);
    E_mixingFactor->setValue(2.0);
    E_runoffPartitioning->setValue(1.0);
    checkFloodInitial->setChecked(false);
    E_FloodMaxIter->setValue(200);

    //E_gravityToChannel->setValue(0);
    //E_angleToChannel->setValue(0.02);

    tabWidget->setCurrentIndex(0);
    tabWidget_out->setCurrentIndex(1);
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
    setfontSize();
}
//---------------------------------------------------------------
void lisemqt::fontDecrease()
{
    genfontsize--;
    genfontsize = std::max(5, genfontsize);
    setfontSize();
}
//---------------------------------------------------------------
void lisemqt::fontIncrease()
{
    genfontsize++;
    genfontsize = std::min(20, genfontsize);
    setfontSize();

}
//---------------------------------------------------------------
void lisemqt::setfontSize()
{
   // int x = SetStyleUISize();
    int fs = genfontsize;

    qApp->setStyleSheet(QString("QLabel {font-size: %1px}\
                                QCheckBox::indicator {width: %1px; height: %1px}\
                                QRadioButton::indicator {width: %1px; height: %1px}\
                                QComboBox {font-size: %1px; padding: 1px 0px 1px 3px}\
                                QLineEdit {font-size: %1px; padding: 1px 1px 1px 1px}\
                                QToolButton {font-size: %1px}\
                                QCheckBox {font-size: %1px; padding:  1px 1px 1px 1px}\
                                QRadioButton {font-size: %1px; 1px 1px 1px 3px}\
                                QSpinBox {width: %1px; height: %1px; font-size: %1px; padding: 0px 0px 0px 0px}\
                                QDoubleSpinBox {width: %2px; height: %1px;font-size: %1px; padding: 0px 0px 0px 0px}\
                                ").arg(fs).arg(fs*2.4));

    tabWidgetOptions->setStyleSheet( QString("font-size: %1px; ").arg(fs) );
    tabWidget->setStyleSheet( QString("font-size: %1px; ").arg(fs) );

    QString S = QString("QGroupBox {font-size: %1px;font-weight: bold;color: black;}").arg(fs);
    //groupBox1->setStyleSheet(S);
    //groupBox2->setStyleSheet(S);
    //groupBox3->setStyleSheet(S);
    //groupBox4->setStyleSheet(S);
    //groupBox5->setStyleSheet(S);
    S = QString("QGroupBox {font-size: %1px;font-weight: bold;color: #1b6fb5;}").arg(fs);
    watergroup->setStyleSheet(S);
    sedgroup->setStyleSheet(S);
    outletgroup->setStyleSheet(S);
    groupTime->setStyleSheet(S);

    S = QString("QToolBox::tab {background-color: #1b6fb5}");
    //tabWidget_erosion->setStyleSheet(S);

}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetCalibration_clicked()
{
    resetTabCalibration();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetFlow_clicked()
{
    resetTabFlow();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetErosion_clicked()
{
    resetTabErosion();
}
//---------------------------------------------------------------

void lisemqt::on_toolButton_help1_clicked()
{
    on_toolButton_help(1);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help2_clicked()
{
    on_toolButton_help(2);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help3_clicked()
{
    on_toolButton_help(3);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help4_clicked()
{
    on_toolButton_help(4);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help5_clicked()
{
    on_toolButton_help(5);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help6_clicked()
{
    on_toolButton_help(6);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help7_clicked()
{
    on_toolButton_help(7);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help8_clicked()
{
    on_toolButton_help(8);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help1a_clicked()
{
    on_toolButton_help(6);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help(int page)
{
    QString filename;
    if (page == 1) filename=":/help1.html";
    if (page == 2) filename=":/help2.html";
    if (page == 3) filename=":/help3.html";
    if (page == 4) filename=":/help4.html";
    if (page == 5) filename=":/help5.html";
    if (page == 6) filename=":/help6.html";
    if (page == 7) filename=":/help7.html";
    if (page == 8) filename=":/help8.html";
    QFile file(filename);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream stream(&file);
    helptxt->setHtml(stream.readAll());
  //  helpbox->show();


    QTextEdit *view = new QTextEdit(helptxt->toHtml());
    view->createStandardContextMenu();
    view->setWindowTitle("Option help");
    view->setMinimumWidth(640);
    view->setMinimumHeight(480);
    view->setAttribute(Qt::WA_DeleteOnClose);

    view->show();
}

void lisemqt::on_check2DDiagonalFlow_toggled(bool checked)
{
    E_pitValue->setEnabled(checked);
    label_135->setEnabled(checked);
}

void lisemqt::on_checkDiffusion_toggled(bool checked)
{
    E_SigmaDiffusion->setEnabled(checked);
    label_101->setEnabled(checked);
    label_139->setEnabled(checked);
}

void lisemqt::on_checkHouses_toggled(bool checked)
{
    checkRaindrum->setEnabled(checked);
    label_157->setEnabled(checked);
}


void lisemqt::on_toolButton_DischargeInName_clicked()
{
    QString path;

    DischargeinDir = findValidDir(DischargeinDir, false);

    path = QFileDialog::getOpenFileName(this,
                                        QString("Select discarge input file"),
                                        DischargeinDir);
    if(!path.isEmpty())
    {
        QFileInfo fi(path);
        DischargeinFileName = fi.fileName();
        DischargeinDir = CheckDir(fi.absolutePath());//Dir().path());
        E_DischargeInName->setText( RainFileDir + DischargeinFileName  );
    }
}

// select a file or directory
// doFile = 0: select a directory;
// dofile = 1 select a file and return file name only;
// dofile = 2 return filename wioth full path
QString lisemqt::getFileorDir(QString inputdir,QString title, QStringList filters, int doFile)
{
    QFileDialog dialog;
    QString dirout = inputdir;
    if (doFile > 0) {
        dialog.setNameFilters(filters);
        dialog.setDirectory(QFileInfo(inputdir).absoluteDir());
        dialog.setFileMode(QFileDialog::ExistingFile);
    } else {
        filters.clear();
        dialog.setNameFilters(filters);
        dialog.setDirectory(QDir(inputdir));
        dialog.setFileMode(QFileDialog::DirectoryOnly);
    }

    dialog.setLabelText(QFileDialog::LookIn,title);
    dialog.exec();

    if (doFile > 0) {
        dirout = "";
        if (dialog.selectedFiles().count() > 0)
            dirout = dialog.selectedFiles().at(0);
        if (doFile == 1)
            dirout = QFileInfo(dirout).fileName();
        if (doFile == 2)
            dirout = QFileInfo(dirout).absoluteFilePath();
    } else {
        dirout = dialog.selectedUrls().at(0).path();
        dirout.remove(0,1);
        if (dirout.lastIndexOf('/') != dirout.length())
            dirout = dirout + "/";
    }

    return dirout;
}


void lisemqt::on_checkEventBased_clicked(bool checked)
{
    /*
    E_BeginTimeDay->setDisabled(checked);
    E_EndTimeDay->setDisabled(checked);
    if (checked) {
        label_6->setText("Begin time (min)");
        label_7->setText("End time (min)");
        E_BeginTimeDay->setInputMask("9999;0");
        E_EndTimeDay->setInputMask("9999;0");
    } else {
        label_6->setText("Begin time (daynr:minnr)");
        label_7->setText("End time (daynr:minnr)");
        E_BeginTimeDay->setInputMask("999:9999;0");
        E_EndTimeDay->setInputMask("999:9999;0");
    }
*/
}

void lisemqt::on_toolButton_rainsatName_clicked()
{
    RainSatFileDir = RainFileDir;
    if (!QFileInfo(RainSatFileDir).exists() || RainSatFileDir.isEmpty())
        RainSatFileDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});
    QString sss = getFileorDir(RainSatFileDir,"Select rainfall map list table", filters, 2);

    RainSatFileDir = QFileInfo(sss).absolutePath()+"/";
    RainSatFileName = QFileInfo(sss).fileName(); //baseName();

    E_rainsatName->setText(RainSatFileDir + RainSatFileName);
    RainFileDir = RainSatFileDir;
}


void lisemqt::on_toolButton_ETName_clicked()
{
    if (!QFileInfo(ETFileDir).exists() || ETFileDir.isEmpty())
        ETFileDir = RainFileDir;
    if (!QFileInfo(ETFileDir).exists() || ETFileDir.isEmpty())
        ETFileDir = currentDir;
    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    ETFileName = getFileorDir(ETFileDir,"Select ET stations file", filters, 2);

    E_ETName->setText(ETFileName);
}


void lisemqt::on_toolButton_ETsatName_clicked()
{
    ETSatFileDir = ETFileDir;

    if (!QFileInfo(ETSatFileDir).exists() || ETSatFileDir.isEmpty())
        ETSatFileDir = RainSatFileDir;
    if (!QFileInfo(ETSatFileDir).exists() || ETSatFileDir.isEmpty())
        ETSatFileDir = currentDir;
    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(ETSatFileDir,"Select ET map list table", filters, 2);

    ETSatFileDir = QFileInfo(sss).absolutePath()+"/";
    ETSatFileName = QFileInfo(sss).fileName(); //baseName();
    E_ETsatName->setText(ETSatFileDir + ETSatFileName);

    ETFileDir = ETSatFileDir;

}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeET_toggled(bool checked)
{
    radioGroupET->setEnabled(checked);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETShow_clicked()
{
    showTextfile(ETFileDir + ETFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallShow_clicked()
{
    showTextfile(RainFileDir + RainFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainmapShow_clicked()
{
    showTextfile(RainSatFileDir + RainSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETmapShow_clicked()
{
    showTextfile(ETSatFileDir + ETSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::showTextfile(QString name)
{

    QFile file(name);
    if (!file.open(QFile::ReadWrite | QFile::Text))
    {
        QMessageBox::warning(this, QString("openLISEM"),
                             QString("Cannot read file %1:\n%2.")
                             .arg(name)
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

void lisemqt::on_E_EndTimeDay_returnPressed()
{
    int daye = E_EndTimeDay->text().split(":")[0].toInt();
    int mine = E_EndTimeDay->text().split(":")[1].toInt();
    daye = std::max(1,std::min(daye, 366));
    if (mine > 1440) {
        daye = mine/1440 + 1;
        mine = mine % 1440;
    }
    E_EndTimeDay->setText(QString("%1:%2").arg(daye,3,10,QLatin1Char('0')).arg(mine,4,10,QLatin1Char('0')));
}


void lisemqt::on_E_BeginTimeDay_returnPressed()
{
       int daye = E_BeginTimeDay->text().split(":")[0].toInt();
       int mine = E_BeginTimeDay->text().split(":")[1].toInt();
       daye = std::max(1,std::min(daye, 366));
       if (mine > 1440) {
           daye = mine/1440 + 1;
           mine = mine % 1440;
       }
       E_BeginTimeDay->setText(QString("%1:%2").arg(daye,3,10,QLatin1Char('0')).arg(mine,4,10,QLatin1Char('0')));
}



void lisemqt::on_toolButton_resetInfiltration_clicked()
{
    resetTabInfiltration();
}

void lisemqt::on_toolButton_resetInterception_clicked()
{
    resetTabInterception();
}

void lisemqt::on_toolButton_resetAdvanced_clicked()
{
    resetTabAdvanced();
}
