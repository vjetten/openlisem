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

#include "lisemqt.h"
#include "global.h"

//---------------------------------------------------------------------------
int lisemqt::SetStyleUISize()
{
    //QRect rect = QGuiApplication::primaryScreen()->availableGeometry();
    int _H = QApplication::desktop()->height();//rect.height();

    int disp = 3;

    if(_H < 1400) disp = 2;
    if(_H < 1200) disp = 1;
    if(_H < 1080) disp = 0;
    if(_H < 800) disp = -1;
   // qDebug() << _H << disp;

    tabWidgetOptions->setMinimumSize(QSize(600, 500));
    scrollArea->setWidgetResizable(true);
    //scrollArea->setWidget(tabWidgetOptions);

    // do a bit of size tweaking for large displays
    QSize iSize = QSize(16,16);
    if (disp == -1) {
        iSize = QSize(16,16);
        tabWidget_out->setIconSize(iSize);
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 40px; width: 28px}");
        tabWidgetOptions->setIconSize(iSize);
        tabWidgetOptions->setStyleSheet("QTabBar::tab { height: 40px; width: 28px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 16px 16px}"));
    }
    if (disp == 0) {
        tabWidget_out->setIconSize(QSize(20, 20));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 48px; width: 32px}");
        tabWidgetOptions->setIconSize(QSize(20, 20));
        tabWidgetOptions->setStyleSheet("QTabBar::tab { height: 48px; width: 32px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 16px 16px}"));
        iSize = QSize(16,16);
    }
    if (disp == 1) {
        tabWidget_out->setIconSize(QSize(24, 24));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 48px; width: 40px}");
        tabWidgetOptions->setIconSize(QSize(24, 24));
        tabWidgetOptions->setStyleSheet("QTabBar::tab { height: 48px; width: 40px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 16px 16px}"));
        iSize = QSize(24,24);
    }
    if (disp == 2) {
        tabWidget_out->setIconSize(QSize(32, 32));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 64px; width: 48px}");
        tabWidgetOptions->setIconSize(QSize(32, 32));
        tabWidgetOptions->setStyleSheet("QTabBar::tab { height: 64px; width: 48px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 24px 24px}"));
        iSize = QSize(32,32);
    }
    if (disp == 3) {
        tabWidget_out->setIconSize(QSize(32, 32));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 64px; width: 48px}");
        tabWidgetOptions->setIconSize(QSize(32, 32));
        tabWidgetOptions->setStyleSheet("QTabBar::tab { height: 64px; width: 48px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 24px 24px}"));
        iSize = QSize(32,32);
    }
    if (disp > 3) {
        tabWidget_out->setIconSize(QSize(48, 48));
        tabWidget_out->setStyleSheet("QTabBar::tab { height: 96px; width: 64px}");
        tabWidgetOptions->setIconSize(QSize(48, 48));
        tabWidgetOptions->setStyleSheet("QTabBar::tab { height: 96px; width: 64px}");
        this->setStyleSheet(QString("QToolButton * {icon-size: 32px 32px}"));
        iSize = QSize(32,32);
    }

    toolBar->setIconSize(iSize);
    toolBar_2->setIconSize(iSize);

    return disp; //-1;
}

// labels in output tab
void lisemqt::setOutputTabStyle(QString bc, QString fc)
{
        label_dx->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
      label_area->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
      label_time->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_endtime->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_raintot->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_ETatot->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_watervoltot->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_stormdraintot->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
      label_qtot->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_infiltot->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_surfstor->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_interctot->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_qpeaktime->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_ppeaktime->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_QPfrac->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_floodVolmm->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_watervolchannel->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_qtotm3sub->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_dischargesub->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_qpeaksub->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_soillosssub->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
     label_Qssub->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_splashdet->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_flowdet->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_sedvol->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
       label_dep->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
     label_detch->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
     label_depch->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_sedvolch->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_soilloss->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
    label_soillosskgha->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));
       label_SDR->setStyleSheet(QString("QLabel { background-color: %1; color: %2; }").arg(bc).arg(fc));

}
//---------------------------------------------------------------------------

void lisemqt::lightStyleUI()
{

    QString sc1 = "#2266aa";
    QString sc = "#4488cc";
    QString bc = "#ffff99";
    QString fc = "#000000";

    tabWidgetOptions->setTabIcon(0,QIcon(":/settings2.png"));
    tabWidgetOptions->setTabIcon(1,QIcon(":/rain.png"));
    tabWidgetOptions->setTabIcon(2,QIcon(":/Plant-icon.png"));
    tabWidgetOptions->setTabIcon(3,QIcon(":/soil5.png"));
    tabWidgetOptions->setTabIcon(4,QIcon(":/water.png"));
    tabWidgetOptions->setTabIcon(5,QIcon(":/river3.png"));
    tabWidgetOptions->setTabIcon(6,QIcon(":/eros1bw.png"));
    tabWidgetOptions->setTabIcon(7,QIcon(":/advanced.png"));
    tabWidgetOptions->setTabIcon(8,QIcon(":/settings1.png"));

    qApp->setStyleSheet(QString(//"* { background-color: #f0f0f0; color: #000000;}"
                                "*:disabled { color: #a4a4a4; }"
                                "QTreeView::item:alternate {background-color: #e0e0e0;}"
                                "QAbstractItemView, QTreeView::item {"
                                   " background: #f0f0f0;"
                                   " background-color: #f0f0f0;"
                                    "selection-background-color: #2266aa;"
                                    "selection-color: #000000;"
                                    "alternate-background-color: #e0e0e0;"
                                "}"
                                // "QAbstractItemView{"   //,QTreeView *
                                //    " background: #f0f0f0;"
                                //    " background-color: #f0f0f0;"
                                //     "selection-background-color: #2266aa;"
                                //     "selection-color: #000000;"
                                //     "alternate-background-color: #e9e9e9;"
                                // "}"
                                "QComboBox {selection-background-color: #4488cc; selection-color: #f0f0f0;}"
                                "QLineEdit {background-color: #ffffff; color: black;}"
                                "QCheckBox::indicator:unchecked {background-color: #ffffff;border: 1px solid #646464; }"
                                //"QRadioButton::indicator {background-color: #ffffff;}"
                                "QTabWidget, QTabWidget::tab-bar {border-color: #292b2d; background-color: #fcfcfc;}"
                                "QGroupBox::title{color: %1;}"
                                "QTabWidget { background-color: #fcfcfc; }"
                                "QGroupBox#groupBoxOutput::title{color: %2;}"
                                "QGroupBox#groupBoxInput::title{color: %2;}"
                                ).arg(sc).arg(sc1)
                        );

    HPlot->setStyleSheet("*{background-color: #fcfcfc; color: #000000;}");
    MPlot->setStyleSheet("*{background-color: #fcfcfc; color: #000000;}");

    label_55->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_59->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_88->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_156->setStyleSheet(QString("QLabel {color:%1;}").arg(sc1));
    label_9->setStyleSheet(QString("QLabel {color:  %1;}").arg(sc1));
    label_10->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_128->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    checkWaveInUser->setStyleSheet(QString("QCheckBox {color: %1;}").arg(sc));

    setOutputTabStyle(bc, fc);
}
void lisemqt::darkStyleUI()
{
    QString sc1 = "#D67537";
    QString sc2 = "#888888";
    QString sc = "#dddd55";
    QString bc = "#a28000";
    QString fc = "#fcfcfc";

    tabWidgetOptions->setTabIcon(0,QIcon(":/d_settings2.png"));
    tabWidgetOptions->setTabIcon(1,QIcon(":/d_rain.png"));
    tabWidgetOptions->setTabIcon(2,QIcon(":/d_Plant-icon.png"));
    tabWidgetOptions->setTabIcon(3,QIcon(":/d_soil5.png"));
    tabWidgetOptions->setTabIcon(4,QIcon(":/d_water.png"));
    tabWidgetOptions->setTabIcon(5,QIcon(":/d_river3.png"));
    tabWidgetOptions->setTabIcon(6,QIcon(":/d_eros1bw.png"));
    tabWidgetOptions->setTabIcon(7,QIcon(":/d_advanced.png"));
    tabWidgetOptions->setTabIcon(8,QIcon(":/d_settings1.png"));

    qApp->setStyleSheet(QString("* { background-color: #2E2F30; color: #D6CF9A}"
                                "*:disabled { color: #888888; }"
                                "QTreeView::item:alternate {background-color: #464646;}"
                                "QAbstractItemView,QTreeView::item {"
                                   " background: #2E2F30;"
                                   " background-color: #2E2F30;"
                                    "selection-background-color: #2f65ca;"
                                    "selection-color: #fefefe;"
                                    "alternate-background-color: #464646;"
                                "}"
                                "QComboBox {selection-background-color: #1D545C; selection-color: #fcfcfc;}" //1D545C 3f76db
                                "QLineEdit {background-color: #606060; color: #fdfdfd;}"
                                //"QLineEdit:disabled {background-color: #4f4f4f;}"
                                "QGroupBox::title{color: %1;}"
                                "QGroupBox::title:disabled{color: %3;}"
                                "QGroupBox#groupBoxOutput::title{color: %2;}"
                                "QGroupBox#groupBoxInput::title{color: %2;}"
                                "QCheckBox::indicator:unchecked{background-color: #606060;border: 1px solid #303030; }"
                                //"QRadioButton::indicator:unchecked{background-color: #606060;border: 1px solid #303030; }"
                                //"QTabWidget, QTabWidget::tab-bar {border-color: #292b2d; background-color: #606000;}"
                                "QTabWidget#tabWidgetOptions:selected {background-color: #1D545C;}"
                                // "QWidget#tab_general  {background-color: #404244;} "
                                // "QWidget#tab_general  {background-color: #404244;} "
                                // "QWidget#tab_meteo    {background-color: #404244;} "
                                // "QWidget#tab_interc   {background-color: #404244;} "
                                // "QWidget#tab_infil    {background-color: #404244;} "
                                // "QWidget#tab_flooding {background-color: #404244;} "
                                // "QWidget#tab_Channel  {background-color: #404244;} "
                                // "QWidget#tab_calib    {background-color: #404244;} "
                                // "QWidget#tab_advanced {background-color: #404244;} "
                                "#label_55  {color: %2;}"
                                "#label_59  {color: %2;}"
                                "#label_88  {color: %2;}"
                                "#label_156 {color: %2;}"
                                "#label_9   {color: %2;}"
                                "#label_10  {color: %2;}"
                                "#label_128 {color: %2;}"
                                "QCheckBox#checkWaveInUser {color: %1;}"

                                ).arg(sc).arg(sc1).arg(sc2)
                        );

    HPlot->setStyleSheet("*{background-color: #e0e0e0; color: #000000;}");
    MPlot->setStyleSheet("*{background-color: #e0e0e0; color: #000000;}");

  //   label_55->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));

   //   label_59->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
   //   label_88->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
   //  label_156->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
   //    label_9->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
   //   label_10->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
   //  label_128->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));

   //  label_55->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //  label_59->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //  label_88->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   // label_156->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //   label_9->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //  label_10->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   // label_128->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));

  // checkWaveInUser->setStyleSheet(QString("QCheckBox {color: %1;}").arg(sc));

   setOutputTabStyle(bc, fc);


}

//---------------------------------------------------------------------------
// general sgtyle things, done once
void lisemqt::SetStyleUI()
{
    trayIcon = new QSystemTrayIcon(this);
    trayIcon->setIcon(QIcon(":/openLisemN.ico"));
    trayIcon->show();
    tabWidgetOptions->tabBar()->setExpanding(true);

    genfontsize = 8+SetStyleUISize();
    setfontSize();

    toolBar_2->setMovable(false);
    toolBar->setMovable(false);

    // interface elements that are not visible for now
    tabWidgetOptions->removeTab(9);

     int w = 80, h = 15;
    label_dx->setMinimumSize(w,h);
    label_area->setMinimumSize(w,h);
    label_time->setMinimumSize(w,h);
    label_endtime->setMinimumSize(w,h);
    label_raintot->setMinimumSize(w,h);
    label_ETatot->setMinimumSize(w,h);
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
    //label_baseflowtot->setMinimumSize(w,h);

    label_qtotm3sub->setMinimumSize(w,h);
    label_dischargesub->setMinimumSize(w,h);
    label_qpeaksub->setMinimumSize(w,h);
    label_soillosssub->setMinimumSize(w,h);
    label_Qssub->setMinimumSize(w,h);

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

    //Grouped Buttons become mututally exclusive
    GroupMapDisplay.addButton(checkBoxComboMaps, 1);
    GroupMapDisplay.addButton(checkBoxComboMaps2, 2);

    //if (checkOverlandFlow2Ddyn->isChecked()) {
    if (E_OFWaveType->currentIndex() == 2) {
        label_107->setText(QString("Flood (mm),h>%1)").arg(E_floodMinHeight->value()*1000));
        label_40->setText(QString("Runoff (mm),h<%1)").arg(E_floodMinHeight->value()*1000));

    } else {
        label_107->setText("Flood mm");
        label_40->setText("Runoff mm");
    }
    bool yes = E_OFWaveType->currentIndex() > 0; // !checkOverlandFlow1D->isChecked();
    label_floodVolmm->setEnabled(yes);
    label_107->setEnabled(yes);

    if (darkLISEM)
        darkStyleUI();
    else
        lightStyleUI();
    //  QString fileName = ":/qtdarcula.css";
    //  QFile file(fileName);
    //  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    //      qDebug() << "Unable to open the file!";
    //  }
    //  QTextStream in(&file);
    //  QString fileContent = in.readAll();
    //  file.close();
    // qApp->setStyleSheet(fileContent);
}
