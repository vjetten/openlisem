/********************************************************************************
** Form generated from reading UI file 'lisemqt.ui'
**
** Created: Sun May 30 12:32:26 2010
**      by: Qt User Interface Compiler version 4.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_LISEMQT_H
#define UI_LISEMQT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QRadioButton>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QToolBar>
#include <QtGui/QToolBox>
#include <QtGui/QToolButton>
#include <QtGui/QTreeView>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_lisemqtClass
{
public:
    QAction *action_Open_runfile;
    QAction *action_Save_runfuile;
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    QTabWidget *tabWidget;
    QWidget *tab;
    QFrame *line;
    QToolBox *toolBox;
    QWidget *page_3;
    QFrame *frame;
    QCheckBox *checkNoErosion;
    QCheckBox *checkIncludeChannel;
    QCheckBox *checkSnowmelt;
    QCheckBox *checkNoErosionOutlet;
    QCheckBox *checkAltErosion;
    QCheckBox *checkAltDepression;
    QCheckBox *checkHardsurface;
    QCheckBox *checkChannelInfil;
    QCheckBox *checkChannelBaseflow;
    QWidget *page_4;
    QComboBox *E_InfiltrationMethod;
    QFrame *frame1;
    QGridLayout *gridLayout_8;
    QCheckBox *checkInfilCompact;
    QCheckBox *checkInfilCrust;
    QCheckBox *checkInfilClosebottom;
    QCheckBox *checkInfil2layer;
    QGroupBox *groupBox_SwatreOptions;
    QGridLayout *gridLayout_7;
    QCheckBox *checkGeometric;
    QCheckBox *checkDumphead;
    QLabel *label_13;
    QLineEdit *E_SWATRETableDir;
    QToolButton *toolButton_SwatreTable;
    QLabel *label_14;
    QDoubleSpinBox *E_SWATREDtsecFraction;
    QWidget *page_5;
    QFrame *frame2;
    QGridLayout *gridLayout_11;
    QCheckBox *checkBuffers;
    QCheckBox *checkSedtrap;
    QCheckBox *checkInfilGrass;
    QLabel *label_15;
    QLineEdit *E_GrassStripN;
    QWidget *page;
    QFrame *frame3;
    QGridLayout *gridLayout_6;
    QLabel *label_9;
    QDoubleSpinBox *E_CalibrateKsat;
    QLabel *label_10;
    QDoubleSpinBox *E_CailbrateN;
    QLabel *label_11;
    QDoubleSpinBox *E_calibrateChKsat;
    QLabel *label_12;
    QDoubleSpinBox *E_CalibrateChN;
    QWidget *page_2;
    QTabWidget *tabWidget_OutputMaps;
    QWidget *tab_4;
    QWidget *widget;
    QVBoxLayout *verticalLayout;
    QCheckBox *checkBox_OutRunoff;
    QCheckBox *checkBox_OutWH;
    QCheckBox *checkBox_OutWHC;
    QCheckBox *checkBox_OutInf;
    QCheckBox *checkBox_OutV;
    QCheckBox *checkBox_OutDep;
    QCheckBox *checkBox_OutDet;
    QCheckBox *checkBox_OutConc;
    QCheckBox *checkBox_OutTC;
    QCheckBox *checkBox_OutSurfStor;
    QCheckBox *checkBox_OutChanVol;
    QWidget *tab_6;
    QWidget *tab_5;
    QWidget *tab_7;
    QGroupBox *groupBox_4;
    QSpinBox *spinBox;
    QLabel *label_16;
    QRadioButton *radioButton;
    QGroupBox *groupBoxTime;
    QGridLayout *gridLayout_3;
    QLabel *label_6;
    QLineEdit *E_BeginTime;
    QLabel *label_7;
    QLineEdit *E_EndTIme;
    QLabel *label_8;
    QLineEdit *E_TimeStep;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_2;
    QLabel *label_2;
    QLineEdit *E_MapDir;
    QToolButton *toolButton_MapDir;
    QLabel *label_4;
    QLineEdit *E_RainfallName;
    QToolButton *toolButton_RainfallShow;
    QToolButton *toolButton_RainfallName;
    QLabel *label_5;
    QLineEdit *E_SnowmeltName;
    QToolButton *toolButton_SnowmeltShow;
    QToolButton *toolButton_SnowmeltNameS;
    QLabel *label;
    QToolButton *toolButton_ShowRunfile;
    QComboBox *E_runFileList;
    QGroupBox *groupOutputMain;
    QLabel *label_17;
    QLineEdit *E_DetachmentMap;
    QLabel *label_18;
    QLineEdit *E_DepositionMap;
    QLabel *label_19;
    QLineEdit *E_SoillossMap;
    QLabel *label_22;
    QLineEdit *E_DetachmentMap_2;
    QLineEdit *E_DepositionMap_2;
    QLabel *label_20;
    QLabel *label_3;
    QLineEdit *E_ResultDir;
    QToolButton *toolButton_ResultDir;
    QFrame *line_2;
    QCheckBox *checkBox_4;
    QGroupBox *groupBox_2;
    QRadioButton *checkUnits_kgcell;
    QRadioButton *checkUnits_kgm2;
    QRadioButton *checkUnits_tonha;
    QGroupBox *groupOutputFormat;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *checkBox;
    QCheckBox *checkBox_2;
    QCheckBox *checkBox_3;
    QWidget *tab_12;
    QGridLayout *gridLayout_12;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_10;
    QTreeView *treeView;
    QCheckBox *checkExpandActive;
    QWidget *tab_3;
    QMenuBar *menubar;
    QMenu *menu_File;
    QToolBar *toolBar;
    QStatusBar *statusBar;
    QButtonGroup *buttonGroup;

    void setupUi(QMainWindow *lisemqtClass)
    {
        if (lisemqtClass->objectName().isEmpty())
            lisemqtClass->setObjectName(QString::fromUtf8("lisemqtClass"));
        lisemqtClass->resize(920, 700);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(lisemqtClass->sizePolicy().hasHeightForWidth());
        lisemqtClass->setSizePolicy(sizePolicy);
        QPalette palette;
        QBrush brush(QColor(0, 0, 0, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::ButtonText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::ButtonText, brush);
        QBrush brush1(QColor(122, 122, 122, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Disabled, QPalette::ButtonText, brush1);
        lisemqtClass->setPalette(palette);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/openLisem25.png"), QSize(), QIcon::Normal, QIcon::Off);
        lisemqtClass->setWindowIcon(icon);
        action_Open_runfile = new QAction(lisemqtClass);
        action_Open_runfile->setObjectName(QString::fromUtf8("action_Open_runfile"));
        action_Save_runfuile = new QAction(lisemqtClass);
        action_Save_runfuile->setObjectName(QString::fromUtf8("action_Save_runfuile"));
        centralwidget = new QWidget(lisemqtClass);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout = new QGridLayout(centralwidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        tabWidget = new QTabWidget(centralwidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setEnabled(true);
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        line = new QFrame(tab);
        line->setObjectName(QString::fromUtf8("line"));
        line->setGeometry(QRect(380, 0, 20, 591));
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);
        toolBox = new QToolBox(tab);
        toolBox->setObjectName(QString::fromUtf8("toolBox"));
        toolBox->setGeometry(QRect(10, 240, 371, 351));
        QFont font;
        font.setBold(true);
        font.setUnderline(false);
        font.setWeight(75);
        font.setKerning(true);
        toolBox->setFont(font);
        toolBox->setAutoFillBackground(true);
        toolBox->setFrameShape(QFrame::Box);
        toolBox->setFrameShadow(QFrame::Sunken);
        toolBox->setLineWidth(1);
        toolBox->setMidLineWidth(1);
        page_3 = new QWidget();
        page_3->setObjectName(QString::fromUtf8("page_3"));
        page_3->setGeometry(QRect(0, 0, 369, 244));
        frame = new QFrame(page_3);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setGeometry(QRect(10, 10, 348, 201));
        checkNoErosion = new QCheckBox(frame);
        checkNoErosion->setObjectName(QString::fromUtf8("checkNoErosion"));
        checkNoErosion->setGeometry(QRect(1, 2, 140, 17));
        QFont font1;
        font1.setBold(false);
        font1.setWeight(50);
        checkNoErosion->setFont(font1);
        checkIncludeChannel = new QCheckBox(frame);
        checkIncludeChannel->setObjectName(QString::fromUtf8("checkIncludeChannel"));
        checkIncludeChannel->setGeometry(QRect(1, 24, 231, 17));
        checkIncludeChannel->setFont(font1);
        checkIncludeChannel->setChecked(false);
        checkSnowmelt = new QCheckBox(frame);
        checkSnowmelt->setObjectName(QString::fromUtf8("checkSnowmelt"));
        checkSnowmelt->setGeometry(QRect(1, 90, 303, 17));
        checkSnowmelt->setFont(font1);
        checkNoErosionOutlet = new QCheckBox(frame);
        checkNoErosionOutlet->setObjectName(QString::fromUtf8("checkNoErosionOutlet"));
        checkNoErosionOutlet->setGeometry(QRect(1, 112, 185, 17));
        checkNoErosionOutlet->setFont(font1);
        checkAltErosion = new QCheckBox(frame);
        checkAltErosion->setObjectName(QString::fromUtf8("checkAltErosion"));
        checkAltErosion->setGeometry(QRect(1, 134, 312, 17));
        checkAltErosion->setFont(font1);
        checkAltDepression = new QCheckBox(frame);
        checkAltDepression->setObjectName(QString::fromUtf8("checkAltDepression"));
        checkAltDepression->setGeometry(QRect(1, 156, 330, 17));
        checkAltDepression->setFont(font1);
        checkHardsurface = new QCheckBox(frame);
        checkHardsurface->setObjectName(QString::fromUtf8("checkHardsurface"));
        checkHardsurface->setGeometry(QRect(1, 178, 346, 17));
        checkHardsurface->setFont(font1);
        checkChannelInfil = new QCheckBox(frame);
        checkChannelInfil->setObjectName(QString::fromUtf8("checkChannelInfil"));
        checkChannelInfil->setEnabled(false);
        checkChannelInfil->setGeometry(QRect(24, 46, 241, 17));
        checkChannelInfil->setFont(font1);
        checkChannelBaseflow = new QCheckBox(frame);
        checkChannelBaseflow->setObjectName(QString::fromUtf8("checkChannelBaseflow"));
        checkChannelBaseflow->setEnabled(false);
        checkChannelBaseflow->setGeometry(QRect(24, 68, 225, 17));
        checkChannelBaseflow->setFont(font1);
        toolBox->addItem(page_3, QString::fromUtf8("Global Options"));
        page_4 = new QWidget();
        page_4->setObjectName(QString::fromUtf8("page_4"));
        page_4->setGeometry(QRect(0, 0, 369, 244));
        E_InfiltrationMethod = new QComboBox(page_4);
        E_InfiltrationMethod->setObjectName(QString::fromUtf8("E_InfiltrationMethod"));
        E_InfiltrationMethod->setGeometry(QRect(9, 4, 351, 22));
        E_InfiltrationMethod->setFont(font1);
        frame1 = new QFrame(page_4);
        frame1->setObjectName(QString::fromUtf8("frame1"));
        frame1->setGeometry(QRect(10, 30, 260, 91));
        gridLayout_8 = new QGridLayout(frame1);
        gridLayout_8->setSpacing(4);
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        checkInfilCompact = new QCheckBox(frame1);
        checkInfilCompact->setObjectName(QString::fromUtf8("checkInfilCompact"));
        checkInfilCompact->setFont(font1);

        gridLayout_8->addWidget(checkInfilCompact, 0, 0, 1, 1);

        checkInfilCrust = new QCheckBox(frame1);
        checkInfilCrust->setObjectName(QString::fromUtf8("checkInfilCrust"));
        checkInfilCrust->setFont(font1);

        gridLayout_8->addWidget(checkInfilCrust, 1, 0, 1, 1);

        checkInfilClosebottom = new QCheckBox(frame1);
        checkInfilClosebottom->setObjectName(QString::fromUtf8("checkInfilClosebottom"));
        checkInfilClosebottom->setFont(font1);

        gridLayout_8->addWidget(checkInfilClosebottom, 2, 0, 1, 1);

        checkInfil2layer = new QCheckBox(frame1);
        checkInfil2layer->setObjectName(QString::fromUtf8("checkInfil2layer"));
        checkInfil2layer->setFont(font1);

        gridLayout_8->addWidget(checkInfil2layer, 3, 0, 1, 1);

        groupBox_SwatreOptions = new QGroupBox(page_4);
        groupBox_SwatreOptions->setObjectName(QString::fromUtf8("groupBox_SwatreOptions"));
        groupBox_SwatreOptions->setGeometry(QRect(10, 120, 351, 111));
        gridLayout_7 = new QGridLayout(groupBox_SwatreOptions);
        gridLayout_7->setSpacing(4);
        gridLayout_7->setContentsMargins(4, 4, 4, 4);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        checkGeometric = new QCheckBox(groupBox_SwatreOptions);
        checkGeometric->setObjectName(QString::fromUtf8("checkGeometric"));
        checkGeometric->setEnabled(true);
        checkGeometric->setFont(font1);

        gridLayout_7->addWidget(checkGeometric, 0, 0, 1, 2);

        checkDumphead = new QCheckBox(groupBox_SwatreOptions);
        checkDumphead->setObjectName(QString::fromUtf8("checkDumphead"));
        checkDumphead->setEnabled(true);
        QFont font2;
        font2.setBold(false);
        font2.setWeight(50);
        font2.setKerning(false);
        checkDumphead->setFont(font2);
        checkDumphead->setChecked(false);

        gridLayout_7->addWidget(checkDumphead, 1, 0, 1, 2);

        label_13 = new QLabel(groupBox_SwatreOptions);
        label_13->setObjectName(QString::fromUtf8("label_13"));
        label_13->setEnabled(true);
        label_13->setFont(font1);
        label_13->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        gridLayout_7->addWidget(label_13, 2, 0, 1, 1);

        E_SWATRETableDir = new QLineEdit(groupBox_SwatreOptions);
        E_SWATRETableDir->setObjectName(QString::fromUtf8("E_SWATRETableDir"));
        E_SWATRETableDir->setEnabled(true);
        E_SWATRETableDir->setFont(font1);

        gridLayout_7->addWidget(E_SWATRETableDir, 2, 1, 1, 1);

        toolButton_SwatreTable = new QToolButton(groupBox_SwatreOptions);
        toolButton_SwatreTable->setObjectName(QString::fromUtf8("toolButton_SwatreTable"));
        toolButton_SwatreTable->setEnabled(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/fileopen.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButton_SwatreTable->setIcon(icon1);

        gridLayout_7->addWidget(toolButton_SwatreTable, 2, 2, 1, 2);

        label_14 = new QLabel(groupBox_SwatreOptions);
        label_14->setObjectName(QString::fromUtf8("label_14"));
        label_14->setEnabled(true);
        label_14->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        gridLayout_7->addWidget(label_14, 3, 0, 1, 2);

        E_SWATREDtsecFraction = new QDoubleSpinBox(groupBox_SwatreOptions);
        E_SWATREDtsecFraction->setObjectName(QString::fromUtf8("E_SWATREDtsecFraction"));
        E_SWATREDtsecFraction->setEnabled(true);
        E_SWATREDtsecFraction->setFont(font1);
        E_SWATREDtsecFraction->setMaximum(1);
        E_SWATREDtsecFraction->setSingleStep(0.01);
        E_SWATREDtsecFraction->setValue(0.2);

        gridLayout_7->addWidget(E_SWATREDtsecFraction, 3, 3, 1, 1);

        toolBox->addItem(page_4, QString::fromUtf8("Infiltration"));
        page_5 = new QWidget();
        page_5->setObjectName(QString::fromUtf8("page_5"));
        page_5->setGeometry(QRect(0, 0, 369, 244));
        frame2 = new QFrame(page_5);
        frame2->setObjectName(QString::fromUtf8("frame2"));
        frame2->setGeometry(QRect(19, 13, 283, 91));
        gridLayout_11 = new QGridLayout(frame2);
        gridLayout_11->setContentsMargins(4, 4, 4, 4);
        gridLayout_11->setObjectName(QString::fromUtf8("gridLayout_11"));
        checkBuffers = new QCheckBox(frame2);
        checkBuffers->setObjectName(QString::fromUtf8("checkBuffers"));
        sizePolicy.setHeightForWidth(checkBuffers->sizePolicy().hasHeightForWidth());
        checkBuffers->setSizePolicy(sizePolicy);
        checkBuffers->setFont(font1);

        gridLayout_11->addWidget(checkBuffers, 0, 0, 1, 2);

        checkSedtrap = new QCheckBox(frame2);
        checkSedtrap->setObjectName(QString::fromUtf8("checkSedtrap"));
        sizePolicy.setHeightForWidth(checkSedtrap->sizePolicy().hasHeightForWidth());
        checkSedtrap->setSizePolicy(sizePolicy);
        checkSedtrap->setFont(font1);

        gridLayout_11->addWidget(checkSedtrap, 1, 0, 1, 2);

        checkInfilGrass = new QCheckBox(frame2);
        checkInfilGrass->setObjectName(QString::fromUtf8("checkInfilGrass"));
        sizePolicy.setHeightForWidth(checkInfilGrass->sizePolicy().hasHeightForWidth());
        checkInfilGrass->setSizePolicy(sizePolicy);
        checkInfilGrass->setFont(font1);

        gridLayout_11->addWidget(checkInfilGrass, 2, 0, 1, 1);

        label_15 = new QLabel(frame2);
        label_15->setObjectName(QString::fromUtf8("label_15"));
        label_15->setEnabled(true);
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label_15->sizePolicy().hasHeightForWidth());
        label_15->setSizePolicy(sizePolicy1);
        label_15->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        gridLayout_11->addWidget(label_15, 3, 0, 1, 1);

        E_GrassStripN = new QLineEdit(frame2);
        E_GrassStripN->setObjectName(QString::fromUtf8("E_GrassStripN"));
        E_GrassStripN->setEnabled(false);
        sizePolicy.setHeightForWidth(E_GrassStripN->sizePolicy().hasHeightForWidth());
        E_GrassStripN->setSizePolicy(sizePolicy);

        gridLayout_11->addWidget(E_GrassStripN, 3, 1, 1, 1);

        toolBox->addItem(page_5, QString::fromUtf8("Conservation measures"));
        page = new QWidget();
        page->setObjectName(QString::fromUtf8("page"));
        page->setGeometry(QRect(0, 0, 369, 244));
        frame3 = new QFrame(page);
        frame3->setObjectName(QString::fromUtf8("frame3"));
        frame3->setGeometry(QRect(10, 10, 281, 108));
        gridLayout_6 = new QGridLayout(frame3);
        gridLayout_6->setContentsMargins(4, 4, 4, 4);
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        label_9 = new QLabel(frame3);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setFont(font1);
        label_9->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_6->addWidget(label_9, 0, 0, 1, 1);

        E_CalibrateKsat = new QDoubleSpinBox(frame3);
        E_CalibrateKsat->setObjectName(QString::fromUtf8("E_CalibrateKsat"));
        sizePolicy.setHeightForWidth(E_CalibrateKsat->sizePolicy().hasHeightForWidth());
        E_CalibrateKsat->setSizePolicy(sizePolicy);
        E_CalibrateKsat->setFont(font1);
        E_CalibrateKsat->setSingleStep(0.01);
        E_CalibrateKsat->setValue(1);

        gridLayout_6->addWidget(E_CalibrateKsat, 0, 1, 1, 1);

        label_10 = new QLabel(frame3);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setFont(font1);
        label_10->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_6->addWidget(label_10, 1, 0, 1, 1);

        E_CailbrateN = new QDoubleSpinBox(frame3);
        E_CailbrateN->setObjectName(QString::fromUtf8("E_CailbrateN"));
        E_CailbrateN->setFont(font1);
        E_CailbrateN->setSingleStep(0.01);
        E_CailbrateN->setValue(1);

        gridLayout_6->addWidget(E_CailbrateN, 1, 1, 1, 1);

        label_11 = new QLabel(frame3);
        label_11->setObjectName(QString::fromUtf8("label_11"));
        label_11->setFont(font1);
        label_11->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_6->addWidget(label_11, 2, 0, 1, 1);

        E_calibrateChKsat = new QDoubleSpinBox(frame3);
        E_calibrateChKsat->setObjectName(QString::fromUtf8("E_calibrateChKsat"));
        E_calibrateChKsat->setFont(font1);
        E_calibrateChKsat->setSingleStep(0.01);
        E_calibrateChKsat->setValue(1);

        gridLayout_6->addWidget(E_calibrateChKsat, 2, 1, 1, 1);

        label_12 = new QLabel(frame3);
        label_12->setObjectName(QString::fromUtf8("label_12"));
        label_12->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_6->addWidget(label_12, 3, 0, 1, 1);

        E_CalibrateChN = new QDoubleSpinBox(frame3);
        E_CalibrateChN->setObjectName(QString::fromUtf8("E_CalibrateChN"));
        E_CalibrateChN->setFont(font1);
        E_CalibrateChN->setSingleStep(0.01);
        E_CalibrateChN->setValue(1);

        gridLayout_6->addWidget(E_CalibrateChN, 3, 1, 1, 1);

        toolBox->addItem(page, QString::fromUtf8("Calibration"));
        page_2 = new QWidget();
        page_2->setObjectName(QString::fromUtf8("page_2"));
        page_2->setGeometry(QRect(0, 0, 369, 244));
        toolBox->addItem(page_2, QString::fromUtf8("Interception"));
        tabWidget_OutputMaps = new QTabWidget(tab);
        tabWidget_OutputMaps->setObjectName(QString::fromUtf8("tabWidget_OutputMaps"));
        tabWidget_OutputMaps->setEnabled(true);
        tabWidget_OutputMaps->setGeometry(QRect(400, 290, 291, 301));
        tab_4 = new QWidget();
        tab_4->setObjectName(QString::fromUtf8("tab_4"));
        widget = new QWidget(tab_4);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(10, 10, 271, 246));
        verticalLayout = new QVBoxLayout(widget);
        verticalLayout->setSpacing(4);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        checkBox_OutRunoff = new QCheckBox(widget);
        checkBox_OutRunoff->setObjectName(QString::fromUtf8("checkBox_OutRunoff"));

        verticalLayout->addWidget(checkBox_OutRunoff);

        checkBox_OutWH = new QCheckBox(widget);
        checkBox_OutWH->setObjectName(QString::fromUtf8("checkBox_OutWH"));

        verticalLayout->addWidget(checkBox_OutWH);

        checkBox_OutWHC = new QCheckBox(widget);
        checkBox_OutWHC->setObjectName(QString::fromUtf8("checkBox_OutWHC"));

        verticalLayout->addWidget(checkBox_OutWHC);

        checkBox_OutInf = new QCheckBox(widget);
        checkBox_OutInf->setObjectName(QString::fromUtf8("checkBox_OutInf"));

        verticalLayout->addWidget(checkBox_OutInf);

        checkBox_OutV = new QCheckBox(widget);
        checkBox_OutV->setObjectName(QString::fromUtf8("checkBox_OutV"));

        verticalLayout->addWidget(checkBox_OutV);

        checkBox_OutDep = new QCheckBox(widget);
        checkBox_OutDep->setObjectName(QString::fromUtf8("checkBox_OutDep"));

        verticalLayout->addWidget(checkBox_OutDep);

        checkBox_OutDet = new QCheckBox(widget);
        checkBox_OutDet->setObjectName(QString::fromUtf8("checkBox_OutDet"));

        verticalLayout->addWidget(checkBox_OutDet);

        checkBox_OutConc = new QCheckBox(widget);
        checkBox_OutConc->setObjectName(QString::fromUtf8("checkBox_OutConc"));

        verticalLayout->addWidget(checkBox_OutConc);

        checkBox_OutTC = new QCheckBox(widget);
        checkBox_OutTC->setObjectName(QString::fromUtf8("checkBox_OutTC"));

        verticalLayout->addWidget(checkBox_OutTC);

        checkBox_OutSurfStor = new QCheckBox(widget);
        checkBox_OutSurfStor->setObjectName(QString::fromUtf8("checkBox_OutSurfStor"));

        verticalLayout->addWidget(checkBox_OutSurfStor);

        checkBox_OutChanVol = new QCheckBox(widget);
        checkBox_OutChanVol->setObjectName(QString::fromUtf8("checkBox_OutChanVol"));

        verticalLayout->addWidget(checkBox_OutChanVol);

        tabWidget_OutputMaps->addTab(tab_4, QString());
        tab_6 = new QWidget();
        tab_6->setObjectName(QString::fromUtf8("tab_6"));
        tabWidget_OutputMaps->addTab(tab_6, QString());
        tab_5 = new QWidget();
        tab_5->setObjectName(QString::fromUtf8("tab_5"));
        tabWidget_OutputMaps->addTab(tab_5, QString());
        tab_7 = new QWidget();
        tab_7->setObjectName(QString::fromUtf8("tab_7"));
        tabWidget_OutputMaps->addTab(tab_7, QString());
        groupBox_4 = new QGroupBox(tab);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        groupBox_4->setGeometry(QRect(700, 290, 101, 311));
        groupBox_4->setFlat(false);
        spinBox = new QSpinBox(groupBox_4);
        spinBox->setObjectName(QString::fromUtf8("spinBox"));
        spinBox->setGeometry(QRect(60, 20, 42, 22));
        label_16 = new QLabel(groupBox_4);
        label_16->setObjectName(QString::fromUtf8("label_16"));
        label_16->setGeometry(QRect(10, 11, 51, 41));
        label_16->setWordWrap(true);
        radioButton = new QRadioButton(groupBox_4);
        radioButton->setObjectName(QString::fromUtf8("radioButton"));
        radioButton->setGeometry(QRect(10, 50, 91, 31));
        groupBoxTime = new QGroupBox(tab);
        groupBoxTime->setObjectName(QString::fromUtf8("groupBoxTime"));
        groupBoxTime->setGeometry(QRect(10, 140, 181, 91));
        QFont font3;
        font3.setBold(true);
        font3.setWeight(75);
        groupBoxTime->setFont(font3);
        gridLayout_3 = new QGridLayout(groupBoxTime);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setSizeConstraint(QLayout::SetNoConstraint);
        gridLayout_3->setHorizontalSpacing(4);
        gridLayout_3->setVerticalSpacing(2);
        gridLayout_3->setContentsMargins(8, 4, 8, 4);
        label_6 = new QLabel(groupBoxTime);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setFont(font1);
        label_6->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_6, 0, 0, 1, 1);

        E_BeginTime = new QLineEdit(groupBoxTime);
        E_BeginTime->setObjectName(QString::fromUtf8("E_BeginTime"));
        E_BeginTime->setFont(font1);

        gridLayout_3->addWidget(E_BeginTime, 0, 1, 1, 1);

        label_7 = new QLabel(groupBoxTime);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setFont(font1);
        label_7->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_7, 1, 0, 1, 1);

        E_EndTIme = new QLineEdit(groupBoxTime);
        E_EndTIme->setObjectName(QString::fromUtf8("E_EndTIme"));
        E_EndTIme->setFont(font1);

        gridLayout_3->addWidget(E_EndTIme, 1, 1, 1, 1);

        label_8 = new QLabel(groupBoxTime);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setFont(font1);
        label_8->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_8, 2, 0, 1, 1);

        E_TimeStep = new QLineEdit(groupBoxTime);
        E_TimeStep->setObjectName(QString::fromUtf8("E_TimeStep"));
        E_TimeStep->setFont(font1);

        gridLayout_3->addWidget(E_TimeStep, 2, 1, 1, 1);

        groupBox = new QGroupBox(tab);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(10, 10, 371, 121));
        groupBox->setFont(font3);
        gridLayout_2 = new QGridLayout(groupBox);
        gridLayout_2->setSpacing(2);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setContentsMargins(4, 2, 4, 2);
        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setFont(font1);
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_2, 1, 0, 1, 1);

        E_MapDir = new QLineEdit(groupBox);
        E_MapDir->setObjectName(QString::fromUtf8("E_MapDir"));
        E_MapDir->setFont(font1);

        gridLayout_2->addWidget(E_MapDir, 1, 1, 1, 2);

        toolButton_MapDir = new QToolButton(groupBox);
        toolButton_MapDir->setObjectName(QString::fromUtf8("toolButton_MapDir"));
        toolButton_MapDir->setIcon(icon1);

        gridLayout_2->addWidget(toolButton_MapDir, 1, 3, 1, 1);

        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setFont(font1);
        label_4->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_4, 2, 0, 1, 1);

        E_RainfallName = new QLineEdit(groupBox);
        E_RainfallName->setObjectName(QString::fromUtf8("E_RainfallName"));
        E_RainfallName->setFont(font1);

        gridLayout_2->addWidget(E_RainfallName, 2, 1, 1, 1);

        toolButton_RainfallShow = new QToolButton(groupBox);
        toolButton_RainfallShow->setObjectName(QString::fromUtf8("toolButton_RainfallShow"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/onewfile.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButton_RainfallShow->setIcon(icon2);

        gridLayout_2->addWidget(toolButton_RainfallShow, 2, 2, 1, 1);

        toolButton_RainfallName = new QToolButton(groupBox);
        toolButton_RainfallName->setObjectName(QString::fromUtf8("toolButton_RainfallName"));
        toolButton_RainfallName->setEnabled(true);
        toolButton_RainfallName->setIcon(icon1);

        gridLayout_2->addWidget(toolButton_RainfallName, 2, 3, 1, 1);

        label_5 = new QLabel(groupBox);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setEnabled(false);
        label_5->setFont(font1);
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_5, 3, 0, 1, 1);

        E_SnowmeltName = new QLineEdit(groupBox);
        E_SnowmeltName->setObjectName(QString::fromUtf8("E_SnowmeltName"));
        E_SnowmeltName->setEnabled(false);
        E_SnowmeltName->setFont(font1);

        gridLayout_2->addWidget(E_SnowmeltName, 3, 1, 1, 1);

        toolButton_SnowmeltShow = new QToolButton(groupBox);
        toolButton_SnowmeltShow->setObjectName(QString::fromUtf8("toolButton_SnowmeltShow"));
        toolButton_SnowmeltShow->setEnabled(false);
        toolButton_SnowmeltShow->setIcon(icon2);

        gridLayout_2->addWidget(toolButton_SnowmeltShow, 3, 2, 1, 1);

        toolButton_SnowmeltNameS = new QToolButton(groupBox);
        toolButton_SnowmeltNameS->setObjectName(QString::fromUtf8("toolButton_SnowmeltNameS"));
        toolButton_SnowmeltNameS->setEnabled(false);
        toolButton_SnowmeltNameS->setIcon(icon1);

        gridLayout_2->addWidget(toolButton_SnowmeltNameS, 3, 3, 1, 1);

        label = new QLabel(groupBox);
        label->setObjectName(QString::fromUtf8("label"));
        QSizePolicy sizePolicy2(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy2);
        label->setFont(font1);
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label, 0, 0, 1, 1);

        toolButton_ShowRunfile = new QToolButton(groupBox);
        toolButton_ShowRunfile->setObjectName(QString::fromUtf8("toolButton_ShowRunfile"));
        toolButton_ShowRunfile->setIcon(icon2);

        gridLayout_2->addWidget(toolButton_ShowRunfile, 0, 3, 1, 1);

        E_runFileList = new QComboBox(groupBox);
        E_runFileList->setObjectName(QString::fromUtf8("E_runFileList"));
        E_runFileList->setFont(font1);

        gridLayout_2->addWidget(E_runFileList, 0, 1, 1, 2);

        groupOutputMain = new QGroupBox(tab);
        groupOutputMain->setObjectName(QString::fromUtf8("groupOutputMain"));
        groupOutputMain->setGeometry(QRect(400, 10, 401, 181));
        groupOutputMain->setFont(font3);
        label_17 = new QLabel(groupOutputMain);
        label_17->setObjectName(QString::fromUtf8("label_17"));
        label_17->setGeometry(QRect(5, 54, 81, 16));
        label_17->setFont(font1);
        label_17->setLayoutDirection(Qt::RightToLeft);
        E_DetachmentMap = new QLineEdit(groupOutputMain);
        E_DetachmentMap->setObjectName(QString::fromUtf8("E_DetachmentMap"));
        E_DetachmentMap->setGeometry(QRect(91, 54, 271, 20));
        label_18 = new QLabel(groupOutputMain);
        label_18->setObjectName(QString::fromUtf8("label_18"));
        label_18->setGeometry(QRect(5, 79, 73, 16));
        label_18->setFont(font1);
        label_18->setLayoutDirection(Qt::RightToLeft);
        E_DepositionMap = new QLineEdit(groupOutputMain);
        E_DepositionMap->setObjectName(QString::fromUtf8("E_DepositionMap"));
        E_DepositionMap->setGeometry(QRect(91, 79, 271, 20));
        label_19 = new QLabel(groupOutputMain);
        label_19->setObjectName(QString::fromUtf8("label_19"));
        label_19->setGeometry(QRect(5, 104, 57, 16));
        label_19->setFont(font1);
        label_19->setLayoutDirection(Qt::RightToLeft);
        E_SoillossMap = new QLineEdit(groupOutputMain);
        E_SoillossMap->setObjectName(QString::fromUtf8("E_SoillossMap"));
        E_SoillossMap->setGeometry(QRect(91, 104, 271, 20));
        label_22 = new QLabel(groupOutputMain);
        label_22->setObjectName(QString::fromUtf8("label_22"));
        label_22->setGeometry(QRect(5, 129, 69, 16));
        label_22->setFont(font1);
        label_22->setLayoutDirection(Qt::RightToLeft);
        E_DetachmentMap_2 = new QLineEdit(groupOutputMain);
        E_DetachmentMap_2->setObjectName(QString::fromUtf8("E_DetachmentMap_2"));
        E_DetachmentMap_2->setGeometry(QRect(91, 129, 271, 20));
        E_DepositionMap_2 = new QLineEdit(groupOutputMain);
        E_DepositionMap_2->setObjectName(QString::fromUtf8("E_DepositionMap_2"));
        E_DepositionMap_2->setGeometry(QRect(91, 154, 211, 20));
        label_20 = new QLabel(groupOutputMain);
        label_20->setObjectName(QString::fromUtf8("label_20"));
        label_20->setGeometry(QRect(5, 154, 59, 16));
        label_20->setFont(font1);
        label_20->setLayoutDirection(Qt::RightToLeft);
        label_3 = new QLabel(groupOutputMain);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(5, 19, 82, 16));
        label_3->setFont(font1);
        label_3->setLayoutDirection(Qt::RightToLeft);
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        E_ResultDir = new QLineEdit(groupOutputMain);
        E_ResultDir->setObjectName(QString::fromUtf8("E_ResultDir"));
        E_ResultDir->setGeometry(QRect(91, 20, 271, 20));
        QSizePolicy sizePolicy3(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(E_ResultDir->sizePolicy().hasHeightForWidth());
        E_ResultDir->setSizePolicy(sizePolicy3);
        toolButton_ResultDir = new QToolButton(groupOutputMain);
        toolButton_ResultDir->setObjectName(QString::fromUtf8("toolButton_ResultDir"));
        toolButton_ResultDir->setGeometry(QRect(370, 20, 23, 22));
        QSizePolicy sizePolicy4(QSizePolicy::Fixed, QSizePolicy::Minimum);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(toolButton_ResultDir->sizePolicy().hasHeightForWidth());
        toolButton_ResultDir->setSizePolicy(sizePolicy4);
        toolButton_ResultDir->setIcon(icon1);
        line_2 = new QFrame(groupOutputMain);
        line_2->setObjectName(QString::fromUtf8("line_2"));
        line_2->setGeometry(QRect(5, 40, 391, 16));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);
        checkBox_4 = new QCheckBox(groupOutputMain);
        checkBox_4->setObjectName(QString::fromUtf8("checkBox_4"));
        checkBox_4->setGeometry(QRect(307, 155, 89, 17));
        checkBox_4->setFont(font1);
        groupBox_2 = new QGroupBox(tab);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setGeometry(QRect(720, 200, 81, 80));
        checkUnits_kgcell = new QRadioButton(groupBox_2);
        buttonGroup = new QButtonGroup(lisemqtClass);
        buttonGroup->setObjectName(QString::fromUtf8("buttonGroup"));
        buttonGroup->addButton(checkUnits_kgcell);
        checkUnits_kgcell->setObjectName(QString::fromUtf8("checkUnits_kgcell"));
        checkUnits_kgcell->setGeometry(QRect(10, 55, 61, 17));
        checkUnits_kgm2 = new QRadioButton(groupBox_2);
        buttonGroup->addButton(checkUnits_kgm2);
        checkUnits_kgm2->setObjectName(QString::fromUtf8("checkUnits_kgm2"));
        checkUnits_kgm2->setGeometry(QRect(10, 35, 51, 17));
        checkUnits_tonha = new QRadioButton(groupBox_2);
        buttonGroup->addButton(checkUnits_tonha);
        checkUnits_tonha->setObjectName(QString::fromUtf8("checkUnits_tonha"));
        checkUnits_tonha->setGeometry(QRect(10, 15, 61, 17));
        checkUnits_tonha->setChecked(true);
        groupOutputFormat = new QGroupBox(tab);
        groupOutputFormat->setObjectName(QString::fromUtf8("groupOutputFormat"));
        groupOutputFormat->setGeometry(QRect(400, 200, 311, 81));
        verticalLayout_2 = new QVBoxLayout(groupOutputFormat);
        verticalLayout_2->setContentsMargins(4, 4, 4, 4);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        checkBox = new QCheckBox(groupOutputFormat);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));

        verticalLayout_2->addWidget(checkBox);

        checkBox_2 = new QCheckBox(groupOutputFormat);
        checkBox_2->setObjectName(QString::fromUtf8("checkBox_2"));

        verticalLayout_2->addWidget(checkBox_2);

        checkBox_3 = new QCheckBox(groupOutputFormat);
        checkBox_3->setObjectName(QString::fromUtf8("checkBox_3"));

        verticalLayout_2->addWidget(checkBox_3);

        tabWidget->addTab(tab, QString());
        tab_12 = new QWidget();
        tab_12->setObjectName(QString::fromUtf8("tab_12"));
        gridLayout_12 = new QGridLayout(tab_12);
        gridLayout_12->setObjectName(QString::fromUtf8("gridLayout_12"));
        groupBox_3 = new QGroupBox(tab_12);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        QSizePolicy sizePolicy5(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(groupBox_3->sizePolicy().hasHeightForWidth());
        groupBox_3->setSizePolicy(sizePolicy5);
        gridLayout_10 = new QGridLayout(groupBox_3);
        gridLayout_10->setObjectName(QString::fromUtf8("gridLayout_10"));
        treeView = new QTreeView(groupBox_3);
        treeView->setObjectName(QString::fromUtf8("treeView"));
        treeView->setFrameShape(QFrame::StyledPanel);
        treeView->setFrameShadow(QFrame::Sunken);
        treeView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        treeView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);

        gridLayout_10->addWidget(treeView, 0, 1, 1, 1);

        checkExpandActive = new QCheckBox(groupBox_3);
        checkExpandActive->setObjectName(QString::fromUtf8("checkExpandActive"));

        gridLayout_10->addWidget(checkExpandActive, 1, 1, 1, 1);


        gridLayout_12->addWidget(groupBox_3, 0, 0, 1, 1);

        tabWidget->addTab(tab_12, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        tabWidget->addTab(tab_3, QString());

        gridLayout->addWidget(tabWidget, 0, 1, 1, 1);

        lisemqtClass->setCentralWidget(centralwidget);
        menubar = new QMenuBar(lisemqtClass);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 920, 21));
        menu_File = new QMenu(menubar);
        menu_File->setObjectName(QString::fromUtf8("menu_File"));
        lisemqtClass->setMenuBar(menubar);
        toolBar = new QToolBar(lisemqtClass);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        QSizePolicy sizePolicy6(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(toolBar->sizePolicy().hasHeightForWidth());
        toolBar->setSizePolicy(sizePolicy6);
        toolBar->setIconSize(QSize(18, 18));
        toolBar->setFloatable(false);
        lisemqtClass->addToolBar(Qt::TopToolBarArea, toolBar);
        statusBar = new QStatusBar(lisemqtClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        lisemqtClass->setStatusBar(statusBar);
        QWidget::setTabOrder(toolButton_MapDir, toolButton_RainfallShow);
        QWidget::setTabOrder(toolButton_RainfallShow, toolButton_SnowmeltShow);
        QWidget::setTabOrder(toolButton_SnowmeltShow, E_MapDir);
        QWidget::setTabOrder(E_MapDir, E_RainfallName);
        QWidget::setTabOrder(E_RainfallName, E_SnowmeltName);
        QWidget::setTabOrder(E_SnowmeltName, E_BeginTime);
        QWidget::setTabOrder(E_BeginTime, E_TimeStep);
        QWidget::setTabOrder(E_TimeStep, E_InfiltrationMethod);
        QWidget::setTabOrder(E_InfiltrationMethod, checkInfilCrust);
        QWidget::setTabOrder(checkInfilCrust, checkInfilCompact);
        QWidget::setTabOrder(checkInfilCompact, checkInfilClosebottom);
        QWidget::setTabOrder(checkInfilClosebottom, E_EndTIme);

        menubar->addAction(menu_File->menuAction());

        retranslateUi(lisemqtClass);

        tabWidget->setCurrentIndex(0);
        toolBox->setCurrentIndex(0);
        toolBox->layout()->setSpacing(0);
        tabWidget_OutputMaps->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(lisemqtClass);
    } // setupUi

    void retranslateUi(QMainWindow *lisemqtClass)
    {
        lisemqtClass->setWindowTitle(QApplication::translate("lisemqtClass", "MainWindow", 0, QApplication::UnicodeUTF8));
        action_Open_runfile->setText(QApplication::translate("lisemqtClass", "&Open runfile", 0, QApplication::UnicodeUTF8));
        action_Save_runfuile->setText(QApplication::translate("lisemqtClass", "&Save runfuile", 0, QApplication::UnicodeUTF8));
        checkNoErosion->setText(QApplication::translate("lisemqtClass", "Runoff only (no erosion)", 0, QApplication::UnicodeUTF8));
        checkIncludeChannel->setText(QApplication::translate("lisemqtClass", "Include main channel (needs channel maps)", 0, QApplication::UnicodeUTF8));
        checkSnowmelt->setText(QApplication::translate("lisemqtClass", "Include snowmelt (needs sowmelt file and additional maps)", 0, QApplication::UnicodeUTF8));
        checkNoErosionOutlet->setText(QApplication::translate("lisemqtClass", "Exclude outlet erosion/deposition ", 0, QApplication::UnicodeUTF8));
        checkAltErosion->setText(QApplication::translate("lisemqtClass", "Alternative flow detachment (no settling vel. in detachment)", 0, QApplication::UnicodeUTF8));
        checkAltDepression->setText(QApplication::translate("lisemqtClass", "Simple surface storage (no runoff below max storage threshold)", 0, QApplication::UnicodeUTF8));
        checkHardsurface->setText(QApplication::translate("lisemqtClass", "Include hard surfaces (no interc, infil, detach, needs hardsurf.map)", 0, QApplication::UnicodeUTF8));
        checkChannelInfil->setText(QApplication::translate("lisemqtClass", "Channel infiltration (needs aditional map)", 0, QApplication::UnicodeUTF8));
        checkChannelBaseflow->setText(QApplication::translate("lisemqtClass", "Channel baseflow (needs additional maps)", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page_3), QApplication::translate("lisemqtClass", "Global Options", 0, QApplication::UnicodeUTF8));
        checkInfilCompact->setText(QApplication::translate("lisemqtClass", "Include compacted areas", 0, QApplication::UnicodeUTF8));
        checkInfilCrust->setText(QApplication::translate("lisemqtClass", "Include crusts", 0, QApplication::UnicodeUTF8));
        checkInfilClosebottom->setText(QApplication::translate("lisemqtClass", "Impermeable lower soil boundary", 0, QApplication::UnicodeUTF8));
        checkInfil2layer->setText(QApplication::translate("lisemqtClass", "Include 2nd soil layer", 0, QApplication::UnicodeUTF8));
        groupBox_SwatreOptions->setTitle(QApplication::translate("lisemqtClass", "Swatre", 0, QApplication::UnicodeUTF8));
        checkGeometric->setText(QApplication::translate("lisemqtClass", "geometric average Ksat", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkDumphead->setToolTip(QApplication::translate("lisemqtClass", "Specify pixels 1,2,3,4 in HEADOUT.MAP", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        checkDumphead->setText(QApplication::translate("lisemqtClass", "Save matrix head at 4 locations", 0, QApplication::UnicodeUTF8));
        label_13->setText(QApplication::translate("lisemqtClass", "Table directory", 0, QApplication::UnicodeUTF8));
        toolButton_SwatreTable->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("lisemqtClass", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">Fraction internal timestep versus LISEMtimestep</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page_4), QApplication::translate("lisemqtClass", "Infiltration", 0, QApplication::UnicodeUTF8));
        checkBuffers->setText(QApplication::translate("lisemqtClass", "Include water buffers", 0, QApplication::UnicodeUTF8));
        checkSedtrap->setText(QApplication::translate("lisemqtClass", "Include sediment traps", 0, QApplication::UnicodeUTF8));
        checkInfilGrass->setText(QApplication::translate("lisemqtClass", "Include grasstrips", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("lisemqtClass", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">Mannings n grass strip</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page_5), QApplication::translate("lisemqtClass", "Conservation measures", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("lisemqtClass", "Multiplication factor Ksat slopes", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("lisemqtClass", "Multiplication factor Manning N slopes", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("lisemqtClass", "Multiplication factor Ksat Channel", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("lisemqtClass", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">Multiplication factor Manning N Channel</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page), QApplication::translate("lisemqtClass", "Calibration", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page_2), QApplication::translate("lisemqtClass", "Interception", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkBox_OutRunoff->setToolTip(QApplication::translate("lisemqtClass", "mapseries saved as is \"ro\"", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        checkBox_OutRunoff->setText(QApplication::translate("lisemqtClass", "ro - runoff slope & channel (l/s)", 0, QApplication::UnicodeUTF8));
        checkBox_OutWH->setText(QApplication::translate("lisemqtClass", "wh - water depth on surface (mm)", 0, QApplication::UnicodeUTF8));
        checkBox_OutWHC->setText(QApplication::translate("lisemqtClass", "whc - cumulative runof depth (mm", 0, QApplication::UnicodeUTF8));
        checkBox_OutInf->setText(QApplication::translate("lisemqtClass", "inf - cumulative infilktration (mm)", 0, QApplication::UnicodeUTF8));
        checkBox_OutV->setText(QApplication::translate("lisemqtClass", "velo - velocity (m/s)", 0, QApplication::UnicodeUTF8));
        checkBox_OutDep->setText(QApplication::translate("lisemqtClass", "dep - deposition (select units)", 0, QApplication::UnicodeUTF8));
        checkBox_OutDet->setText(QApplication::translate("lisemqtClass", "det - detachment (selected units)", 0, QApplication::UnicodeUTF8));
        checkBox_OutConc->setText(QApplication::translate("lisemqtClass", "conc - sediment concentration (g/l)", 0, QApplication::UnicodeUTF8));
        checkBox_OutTC->setText(QApplication::translate("lisemqtClass", "tc - transport capacity (kg/m3)", 0, QApplication::UnicodeUTF8));
        checkBox_OutSurfStor->setText(QApplication::translate("lisemqtClass", "sstor - surface storage (mm)", 0, QApplication::UnicodeUTF8));
        checkBox_OutChanVol->setText(QApplication::translate("lisemqtClass", "chanvol - channel volume (m3)", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_4), QApplication::translate("lisemqtClass", "Basic Output", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_6), QApplication::translate("lisemqtClass", "Multiclass Output", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_5), QApplication::translate("lisemqtClass", "Nutrient Output", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_7), QApplication::translate("lisemqtClass", "Gully Output", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("lisemqtClass", "Output times", 0, QApplication::UnicodeUTF8));
        label_16->setText(QApplication::translate("lisemqtClass", "Skip timesteps", 0, QApplication::UnicodeUTF8));
        radioButton->setText(QApplication::translate("lisemqtClass", "Skip timesteps", 0, QApplication::UnicodeUTF8));
        groupBoxTime->setTitle(QApplication::translate("lisemqtClass", "Simulation times", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("lisemqtClass", "Begin time (min)", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("lisemqtClass", "End time (min)", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("lisemqtClass", "Timestep (sec)", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("lisemqtClass", "Input", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("lisemqtClass", "Map directory  ", 0, QApplication::UnicodeUTF8));
        toolButton_MapDir->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("lisemqtClass", "Rainfall file  ", 0, QApplication::UnicodeUTF8));
        toolButton_RainfallShow->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        toolButton_RainfallName->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("lisemqtClass", "Snowmelt file  ", 0, QApplication::UnicodeUTF8));
        toolButton_SnowmeltShow->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        toolButton_SnowmeltNameS->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("lisemqtClass", "Run file(s)  ", 0, QApplication::UnicodeUTF8));
        toolButton_ShowRunfile->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        groupOutputMain->setTitle(QApplication::translate("lisemqtClass", "Output", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("lisemqtClass", "Detachment Map", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("lisemqtClass", "Deposition Map", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("lisemqtClass", "Soilloss Map", 0, QApplication::UnicodeUTF8));
        label_22->setText(QApplication::translate("lisemqtClass", "Text file Totals", 0, QApplication::UnicodeUTF8));
        label_20->setText(QApplication::translate("lisemqtClass", "Point output", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("lisemqtClass", "Result directory  ", 0, QApplication::UnicodeUTF8));
        toolButton_ResultDir->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        checkBox_4->setText(QApplication::translate("lisemqtClass", "Separate files", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("lisemqtClass", "Map units", 0, QApplication::UnicodeUTF8));
        checkUnits_kgcell->setText(QApplication::translate("lisemqtClass", "kg/cell", 0, QApplication::UnicodeUTF8));
        checkUnits_kgm2->setText(QApplication::translate("lisemqtClass", "kg/m2", 0, QApplication::UnicodeUTF8));
        checkUnits_tonha->setText(QApplication::translate("lisemqtClass", "ton/ha", 0, QApplication::UnicodeUTF8));
        groupOutputFormat->setTitle(QApplication::translate("lisemqtClass", "Output format", 0, QApplication::UnicodeUTF8));
        checkBox->setText(QApplication::translate("lisemqtClass", "hydrographs PCRaster format (else comma delimited)", 0, QApplication::UnicodeUTF8));
        checkBox_2->setText(QApplication::translate("lisemqtClass", "Hydrographs SOBEK format, date string: ", 0, QApplication::UnicodeUTF8));
        checkBox_3->setText(QApplication::translate("lisemqtClass", "Mapseries as PCRaster filenames (.001, .002)", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("lisemqtClass", "Input/Output", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("lisemqtClass", "Input Maps", 0, QApplication::UnicodeUTF8));
        checkExpandActive->setText(QApplication::translate("lisemqtClass", "Expand activated maps categories", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_12), QApplication::translate("lisemqtClass", "Maps", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("lisemqtClass", "Simulation", 0, QApplication::UnicodeUTF8));
        menu_File->setTitle(QApplication::translate("lisemqtClass", "&File", 0, QApplication::UnicodeUTF8));
        toolBar->setWindowTitle(QApplication::translate("lisemqtClass", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class lisemqtClass: public Ui_lisemqtClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_LISEMQT_H
