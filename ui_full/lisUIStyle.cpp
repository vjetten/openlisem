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
void lisemqt::SetStyleUISize()
{
    QScreen *screen = QGuiApplication::primaryScreen();
    QRect screenGeometry = screen->geometry();
    int _H = screenGeometry.height();// * screen->devicePixelRatio();
    int disp = 3;
    // for (int i = 0; i < screens.size(); ++i) {
    //     QScreen *screen = screens.at(i);
    //     qreal logicalDpi = screen->logicalDotsPerInch();
    //     qreal physicalDpi = screen->physicalDotsPerInch();
    //     qreal devicePixelRatio = screen->devicePixelRatio();
    //     qDebug() << "Screen" << i << ":";
    //     qDebug() << "  Logical DPI:" << logicalDpi;
    //     qDebug() << "  Physical DPI:" << physicalDpi;
    //     qDebug() << "  Device Pixel Ratio:" << devicePixelRatio;
    // }
    if(_H < 1440) disp = 2;
    if(_H < 1280) disp = 1;
    if(_H < 1080) disp = 0;
    if(_H < 800) disp = -1;
    qDebug() << _H << disp;

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
    toolBar_2->setMovable(false);
    toolBar->setMovable(false);

    genfontsize = screen->devicePixelRatio()*(disp+8);
    qDebug() << genfontsize << screen->devicePixelRatio();
    //qApp->setStyleSheet(QString("* { font-size: %1px; }").arg(genfontsize));
    setfontSize();
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
    QString ac = "#e9e9e9";
    QString hc = "#fcfcfc";
    QString sbc = "#5599dd";   // dropdown background selected

    tabWidgetOptions->setTabIcon(0,QIcon(":/settings2.png"));
    tabWidgetOptions->setTabIcon(1,QIcon(":/rain.png"));
    tabWidgetOptions->setTabIcon(2,QIcon(":/Plant-icon.png"));
    tabWidgetOptions->setTabIcon(3,QIcon(":/soil5.png"));
    tabWidgetOptions->setTabIcon(4,QIcon(":/water1.png"));
    tabWidgetOptions->setTabIcon(5,QIcon(":/river4.png"));
    tabWidgetOptions->setTabIcon(6,QIcon(":/house.png"));
    tabWidgetOptions->setTabIcon(7,QIcon(":/eros1bw.png"));
    tabWidgetOptions->setTabIcon(8,QIcon(":/advanced.png"));
    tabWidgetOptions->setTabIcon(9,QIcon(":/settings1.png"));

    setBWAct->setIcon(QIcon(":/black-and-white.png"));
    fontIncreaseAct->setIcon(QIcon(":/2X/fontbigger2X.png"));
    fontDecreaseAct->setIcon(QIcon(":/2X/fontsmaller2X.png"));

    //qApp->setStyleSheet(QString("* { font-size: %1px; }").arg(genfontsize));

    qApp->setStyleSheet(//"* { background-color: #f0f0f0; color: #000000;}"
                        "*:disabled { color: #a4a4a4; }"
                        //"QAbstractItemView, QTreeView *{"
                        "QAbstractItemView, QTreeView *{"
                            "background-color: #f0f0f0;"
                            "selection-background-color: #5599dd;"
                            "selection-color: #fcfcfc;"
                            "alternate-background-color: #e9e9e9;"
                        "}"
                        "QComboBox {selection-background-color: #5599dd; selection-color: #fcfcfc;}"
                        "QComboBox QAbstractItemView {border: 1px solid #4488cc;selection-background-color: #5599dd; selection-color: #fcfcfc;}"

                        "QComboBox:editable {selection-background-color: #e0e0e0; selection-color: #f0f0f0;}"
                        "QLineEdit {background-color: #fcfcfc; color: black;}"
                        "QCheckBox::indicator:unchecked {background-color: #ffffff;border: 1px solid #646464; }"
                        //"QRadioButton::indicator {background-color: #ffffff;}"
                        "QTabWidget, QTabWidget::tab-bar {border-color: #292b2d; background-color: #fcfcfc;}"
                        "QGroupBox::title{color: #4488cc;}"
                        "QTabWidget { background-color: #fcfcfc; }"
                        "QGroupBox#groupBoxOutput::title{color: #2266aa;}"
                        "QGroupBox#groupBoxInput::title{color: #2266aa;}"
                        "QGroupBox#groupRainfall::title{color: #2266aa;}"
                        "QGroupBox#groupInfiltration::title{color: #2266aa;}"
                        "QGroupBox#groupInterception::title{color: #2266aa;}"
                        );

    HPlot->setStyleSheet("*{background-color: #fcfcfc; color: #000000;}");
    MPlot->setStyleSheet("*{background-color: #fcfcfc; color: #000000;}");

    label_55->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_59->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_88->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_89->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_156->setStyleSheet(QString("QLabel {color:%1;}").arg(sc1));
    label_9->setStyleSheet(QString("QLabel {color:  %1;}").arg(sc1));
    label_10->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    label_128->setStyleSheet(QString("QLabel {color: %1;}").arg(sc1));
    //groupWaveUser->setStyleSheet(QString("QCheckBox {color: %1;}").arg(sc));

    setOutputTabStyle(bc, fc);

}
void lisemqt::darkStyleUI()
{
    QString sc1 = "#e57537";
    QString sc2 = "#888888";
    QString sc = "#dddd55";
    QString bc = "##a28000";
    QString fc = "#fcfcfc";

    tabWidgetOptions->setTabIcon(0,QIcon(":/d_settings2.png"));
    tabWidgetOptions->setTabIcon(1,QIcon(":/d_rain.png"));
    tabWidgetOptions->setTabIcon(2,QIcon(":/d_Plant-icon.png"));
    tabWidgetOptions->setTabIcon(3,QIcon(":/d_soil5.png"));
    tabWidgetOptions->setTabIcon(4,QIcon(":/d_water1.png"));
    tabWidgetOptions->setTabIcon(5,QIcon(":/d_river3.png"));
    tabWidgetOptions->setTabIcon(6,QIcon(":/d_eros1bw.png"));
    tabWidgetOptions->setTabIcon(7,QIcon(":/house.png"));
    tabWidgetOptions->setTabIcon(8,QIcon(":/d_advanced.png"));
    tabWidgetOptions->setTabIcon(9,QIcon(":/d_settings1.png"));

    setBWAct->setIcon(QIcon(":/d_black-and-white.png"));
    fontIncreaseAct->setIcon(QIcon(":/2X/d_fontbigger2X.png"));
    fontDecreaseAct->setIcon(QIcon(":/2X/d_fontsmaller2X.png"));

 //   qApp->setStyleSheet(QString("* { font-size: %1px; }").arg(genfontsize));

    qApp->setStyleSheet("* { background-color: #2E2F30; color: #D6CF9A;}"
                        "*:disabled { color: #888888; }"

                        "QTreeView, QAbstractItemView *{ background: #2E2F30; background-color: #2E2F30;"
                            "selection-background-color: #2266aa; selection-color: #f0f0f0;" //#2266aa
                            "alternate-background-color: #363636;}"

                        //"QComboBox {selection-background-color: #1D545C; selection-color: #f0f0f0;}"
                        "QComboBox:selected {background-color: #1D545C; color: #f0f0f0;}"
                        //"QComboBox, QAbstractItemView  *{background-color: #464646; color: #f0f0f0;}"
                        "QComboBox QAbstractItemView {border: 1px solid #1D545C;selection-background-color: #1D545C;"
                                "background-color: #464646; color: #f0f0f0;}"
                        "QComboBox, QLineEdit *{background-color: #464646; color: #f0f0f0}"
                        "QComboBox, QLineEdit:selected *{border: 1px solid #1D545C;}"

                        "QLineEdit {background-color: #464646; color: #f0f0f0;}"
                        //"QLineEdit:selected *{border: 1px solid #1D545C;}"

                        "QGroupBox::title{color: #dddd55;}"
                        "QGroupBox::title:disabled{color: #888888;}"
                        "QGroupBox#groupBoxOutput::title{color: #e57537;}"
                        "QGroupBox#groupBoxInput::title{color: #e57537;}"
                        "QGroupBox#groupRainfall::title{color: #e57537;}"
                        "QGroupBox#groupInfiltration::title{color: #e57537;}"
                        "QGroupBox#groupInterception::title{color: #e57537;}"


                        "QCheckBox::indicator:unchecked{background-color: #606060;border: 1px solid #1D545C; }"

                        "QToolBar {background: #464646}"
                        "QToolBar#toolBar QToolButton {background : #464646;}"
                        "QToolBar#toolBar_2 QToolButton {background : #464646;}"
                        "QToolBar#toolBar QToolButton:checked {background : #565656;}"
                        "QToolBar#toolBar_2 QToolButton:checked {background : #565656;}"

                        "QTabWidget#tabWidgetOptions QTabBar::tab:selected {background-color: #1D545C; "
                                    "height:32px; width: 42px;margin: 0px; padding-top: -15px; padding-bottom: 15px}"
                        "QTabWidget#tabWidget QTabBar::tab:selected {background-color: #1D545C; color #e57537;}"
                        );

    HPlot->setStyleSheet("*{background-color: #e0e0e0; color: #000000;}");
    MPlot->setStyleSheet("*{background-color: #e0e0e0; color: #000000;}");

     label_55->setStyleSheet( "QLabel {color: #42db00;}");  //e57537
     label_59->setStyleSheet( "QLabel {color: #42db00;}");  //e57537
     label_88->setStyleSheet( "QLabel {color: #42db00;}");  //e57537
     label_89->setStyleSheet( "QLabel {color: #42db00;}");  //e57537
    label_156->setStyleSheet( "QLabel {color: #42db00;}");  //e57537
      label_9->setStyleSheet( "QLabel {color: #42db00;}");  //e57537
     label_10->setStyleSheet( "QLabel {color: #42db00;}");  //e57537
    label_128->setStyleSheet( "QLabel {color: #42db00;}");  //e57537

    setOutputTabStyle("#a28000", "#f0f0f0");


    // "QWidget#tab_general  {background-color: #363636;} "
    // "QWidget#tab_meteo    {background-color: #363636;} "
    // "QWidget#tab_interc   {background-color: #363636;} "
    // "QWidget#tab_infil    {background-color: #363636;} "
    // "QWidget#tab_flooding {background-color: #363636;} "
    // "QWidget#tab_Channel  {background-color: #363636;} "
    // "QWidget#tab_calib    {background-color: #363636;} "
    // "QWidget#tab_advanced {background-color: #363636;} "
    //  label_55->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //  label_59->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //  label_88->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   // label_156->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //   label_9->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   //  label_10->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));
   // label_128->setStyleSheet(QString("QLabel:disabled {color: %1;}").arg(sc2));

   //checkWaveInUser->setStyleSheet("QCheckBox {color: #e57537;}");

}

//---------------------------------------------------------------------------
// general sgtyle things, done once
void lisemqt::SetStyleUI()
{

   trayIcon = new QSystemTrayIcon(this);
   trayIcon->setIcon(QIcon(":/openLisemN.ico"));
   trayIcon->show();

   tabWidgetOptions->tabBar()->setExpanding(true);

   SetStyleUISize();

   // interface elements that are not visible for now
   tabWidgetOptions->removeTab(10);

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
}
//---------------------------------------------------------------
void lisemqt::setBWUI()
{
    if(darkLISEM)
        darkLISEM = false;
    else
        darkLISEM = true;

    if (darkLISEM)
        darkStyleUI();
    else
        lightStyleUI();
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
    if (darkLISEM)
        darkStyleUI();
    else
        lightStyleUI();
    //setfontSize();
}
//---------------------------------------------------------------
void lisemqt::fontDecrease()
{
    genfontsize--;
    genfontsize = std::max(6, genfontsize);
    //qApp->setStyleSheet(QString("* { font-size: %1px; }").arg(genfontsize));
    setfontSize();
}
//---------------------------------------------------------------
void lisemqt::fontIncrease()
{
    genfontsize++;
    genfontsize = std::min(32, genfontsize);
    //qApp->setStyleSheet(QString("* { font-size: %1px; }").arg(genfontsize));
    setfontSize();
}
//---------------------------------------------------------------
void lisemqt::setfontSize()
{
    // if (darkLISEM)
    //     darkStyleUI();
    // else
    //     lightStyleUI();
    QFont font = qApp->font();
    const QWidgetList allWidgets = QApplication::allWidgets();
    for (QWidget *widget : allWidgets) {
        QFont font = widget->font();
        font.setPointSize(genfontsize);
        widget->setFont(font);
        widget->update();
    }

    qDebug() <<"F"<<genfontsize;
}


void lisemqt::on_checkMB_WH_toggled(bool checked)
{
    op.SwitchCorrectMB_WH = checked;
}

