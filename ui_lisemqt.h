/********************************************************************************
** Form generated from reading UI file 'lisemqt.ui'
**
** Created: Sun Jun 6 11:03:17 2010
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
#include <QtGui/QPlainTextEdit>
#include <QtGui/QProgressBar>
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
    QCheckBox *checkSimpleDepression;
    QCheckBox *checkHardsurface;
    QCheckBox *checkChannelInfil;
    QCheckBox *checkChannelBaseflow;
    QWidget *page_4;
    QComboBox *E_InfiltrationMethod;
    QFrame *frame1;
    QGridLayout *gridLayout_8;
    QCheckBox *checkInfilCompact;
    QCheckBox *checkInfilCrust;
    QCheckBox *checkImpermeable;
    QCheckBox *checkInfil2layer;
    QGroupBox *groupBox_SwatreOptions;
    QGridLayout *gridLayout_7;
    QCheckBox *checkGeometricMean;
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
    QWidget *page_2;
    QCheckBox *checkInterceptionLAI;
    QGroupBox *groupInterception;
    QGridLayout *gridLayout_5;
    QRadioButton *radioButton_1;
    QRadioButton *radioButton_2;
    QRadioButton *radioButton_3;
    QRadioButton *radioButton_4;
    QRadioButton *radioButton_5;
    QRadioButton *radioButton_6;
    QRadioButton *radioButton_7;
    QRadioButton *radioButton_8;
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
    QTabWidget *tabWidget_OutputMaps;
    QWidget *tab_4;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout;
    QCheckBox *checkBox_OutRunoff;
    QCheckBox *checkBox_OutWH;
    QCheckBox *checkBox_OutWHC;
    QCheckBox *checkBox_OutInf;
    QCheckBox *checkBox_OutV;
    QCheckBox *checkBox_OutDet;
    QCheckBox *checkBox_OutDep;
    QCheckBox *checkBox_OutConc;
    QCheckBox *checkBox_OutTC;
    QCheckBox *checkBox_OutSurfStor;
    QCheckBox *checkBox_OutChanVol;
    QWidget *tab_6;
    QWidget *tab_5;
    QWidget *tab_7;
    QGroupBox *groupBox_4;
    QSpinBox *printinterval;
    QRadioButton *checkOutputTimeStep;
    QRadioButton *checkOutputTimeUser;
    QLabel *label_16;
    QPlainTextEdit *plainTextEdit;
    QGroupBox *groupBoxTime;
    QGridLayout *gridLayout_3;
    QLabel *label_6;
    QLineEdit *E_BeginTime;
    QLabel *label_7;
    QLineEdit *E_EndTime;
    QLabel *label_8;
    QLineEdit *E_Timestep;
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
    QComboBox *E_runFileList;
    QToolButton *toolButton_ShowRunfile;
    QToolButton *toolButton_fileOpen;
    QGroupBox *groupOutputMain;
    QLabel *label_17;
    QLineEdit *E_DetachmentMap;
    QLabel *label_18;
    QLineEdit *E_DepositionMap;
    QLabel *label_19;
    QLineEdit *E_SoillossMap;
    QLabel *label_22;
    QLineEdit *E_MainTotals;
    QLineEdit *E_PointResults;
    QLabel *label_20;
    QLabel *label_3;
    QLineEdit *E_ResultDir;
    QToolButton *toolButton_ResultDir;
    QFrame *line_2;
    QCheckBox *checkSeparateOutput;
    QGroupBox *groupBox_2;
    QRadioButton *checkUnits_kgcell;
    QRadioButton *checkUnits_kgm2;
    QRadioButton *checkUnits_tonha;
    QGroupBox *groupOutputFormat;
    QVBoxLayout *verticalLayout_2;
    QFrame *frame4;
    QGridLayout *gridLayout_4;
    QCheckBox *checkSOBEKOutput;
    QLineEdit *SOBEKdatestring;
    QCheckBox *checkWritePCRtimeplot;
    QCheckBox *checkWritePCRnames;
    QLabel *label_49;
    QWidget *tab_12;
    QGridLayout *gridLayout_12;
    QGroupBox *groupBox_InputMaps;
    QGridLayout *gridLayout_10;
    QTreeView *treeView;
    QCheckBox *checkExpandActive;
    QLabel *label_47;
    QWidget *tab_3;
    QGroupBox *sedgroup;
    QGridLayout *gridLayout_9;
    QLabel *label_50;
    QLabel *label_21;
    QLabel *label_splashdet;
    QLabel *label_41;
    QLabel *label_detch;
    QLabel *label_23;
    QLabel *label_flowdet;
    QLabel *label_43;
    QLabel *label_depch;
    QLabel *label_27;
    QLabel *label_dep;
    QLabel *label_39;
    QLabel *label_sedvolch;
    QLabel *label_25;
    QLabel *label_sedvol;
    QLabel *label_45;
    QLabel *label_MBs;
    QLabel *label_51;
    QGroupBox *groupTime;
    QGridLayout *gridLayout_13;
    QLabel *label_endtime;
    QLabel *label_area;
    QLabel *label_time;
    QLabel *label_29;
    QLabel *label_dx;
    QLabel *label_46;
    QLabel *label_38;
    QLabel *label_48;
    QLabel *label_30;
    QLabel *label_endruntime;
    QLabel *label_runtime;
    QGroupBox *watergroup;
    QGridLayout *gridLayout_15;
    QLabel *label_32;
    QLabel *label_raintot;
    QLabel *label_34;
    QLabel *label_interctot;
    QLabel *label_36;
    QLabel *label_infiltot;
    QLabel *label_37;
    QLabel *label_surfstor;
    QLabel *label_40;
    QLabel *label_watervoltot;
    QLabel *label_42;
    QLabel *label_qtot;
    QLabel *label_MB;
    QLabel *label_44;
    QLabel *label_24;
    QLabel *label_qpeaktime;
    QLabel *label_52;
    QLabel *label_QPfrac;
    QLabel *label_35;
    QLabel *label_qpeak;
    QLabel *label_53;
    QLabel *label_ppeaktime;
    QProgressBar *progressBar;
    QLabel *label_debug;
    QWidget *widgetGraph;
    QPlainTextEdit *textGraph;
    QGroupBox *outletgroup_2;
    QGridLayout *gridLayout_16;
    QLabel *label_buffervol;
    QLabel *label_56;
    QLabel *label_buffersed;
    QLabel *label_57;
    QGroupBox *outletgroup;
    QGridLayout *gridLayout_14;
    QLabel *label_26;
    QLabel *label_qtotm3;
    QLabel *label_28;
    QLabel *label_soilloss;
    QLabel *label_discharge;
    QLabel *label_54;
    QLabel *label_soillosskgha;
    QLabel *label_31;
    QToolBar *toolBar;
    QStatusBar *statusBar;
    QButtonGroup *buttonGroup_2;
    QButtonGroup *buttonGroup_3;
    QButtonGroup *buttonGroup;

    void setupUi(QMainWindow *lisemqtClass)
    {
        if (lisemqtClass->objectName().isEmpty())
            lisemqtClass->setObjectName(QString::fromUtf8("lisemqtClass"));
        lisemqtClass->resize(852, 676);
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
        line->setGeometry(QRect(400, 0, 20, 581));
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);
        toolBox = new QToolBox(tab);
        toolBox->setObjectName(QString::fromUtf8("toolBox"));
        toolBox->setGeometry(QRect(10, 240, 391, 341));
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
        page_3->setGeometry(QRect(0, 0, 389, 229));
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
        checkNoErosion->setCheckable(true);
        checkNoErosion->setChecked(false);
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
        checkSimpleDepression = new QCheckBox(frame);
        checkSimpleDepression->setObjectName(QString::fromUtf8("checkSimpleDepression"));
        checkSimpleDepression->setGeometry(QRect(1, 156, 330, 17));
        checkSimpleDepression->setFont(font1);
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
        page_4->setGeometry(QRect(0, 0, 389, 229));
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

        checkImpermeable = new QCheckBox(frame1);
        checkImpermeable->setObjectName(QString::fromUtf8("checkImpermeable"));
        checkImpermeable->setFont(font1);

        gridLayout_8->addWidget(checkImpermeable, 2, 0, 1, 1);

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
        checkGeometricMean = new QCheckBox(groupBox_SwatreOptions);
        checkGeometricMean->setObjectName(QString::fromUtf8("checkGeometricMean"));
        checkGeometricMean->setEnabled(true);
        checkGeometricMean->setFont(font1);

        gridLayout_7->addWidget(checkGeometricMean, 0, 0, 1, 2);

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
        page_5->setGeometry(QRect(0, 0, 389, 229));
        frame2 = new QFrame(page_5);
        frame2->setObjectName(QString::fromUtf8("frame2"));
        frame2->setGeometry(QRect(19, 13, 231, 91));
        gridLayout_11 = new QGridLayout(frame2);
        gridLayout_11->setContentsMargins(4, 4, 4, 4);
        gridLayout_11->setObjectName(QString::fromUtf8("gridLayout_11"));
        checkBuffers = new QCheckBox(frame2);
        checkBuffers->setObjectName(QString::fromUtf8("checkBuffers"));
        sizePolicy.setHeightForWidth(checkBuffers->sizePolicy().hasHeightForWidth());
        checkBuffers->setSizePolicy(sizePolicy);
        checkBuffers->setFont(font1);

        gridLayout_11->addWidget(checkBuffers, 0, 0, 1, 3);

        checkSedtrap = new QCheckBox(frame2);
        checkSedtrap->setObjectName(QString::fromUtf8("checkSedtrap"));
        sizePolicy.setHeightForWidth(checkSedtrap->sizePolicy().hasHeightForWidth());
        checkSedtrap->setSizePolicy(sizePolicy);
        checkSedtrap->setFont(font1);

        gridLayout_11->addWidget(checkSedtrap, 1, 0, 1, 3);

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
        label_15->setFont(font1);
        label_15->setTextFormat(Qt::PlainText);
        label_15->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignVCenter);

        gridLayout_11->addWidget(label_15, 3, 0, 1, 1);

        E_GrassStripN = new QLineEdit(frame2);
        E_GrassStripN->setObjectName(QString::fromUtf8("E_GrassStripN"));
        E_GrassStripN->setEnabled(false);
        sizePolicy.setHeightForWidth(E_GrassStripN->sizePolicy().hasHeightForWidth());
        E_GrassStripN->setSizePolicy(sizePolicy);
        E_GrassStripN->setFont(font1);

        gridLayout_11->addWidget(E_GrassStripN, 3, 1, 1, 1);

        toolBox->addItem(page_5, QString::fromUtf8("Conservation measures"));
        page_2 = new QWidget();
        page_2->setObjectName(QString::fromUtf8("page_2"));
        page_2->setGeometry(QRect(0, 0, 389, 229));
        checkInterceptionLAI = new QCheckBox(page_2);
        checkInterceptionLAI->setObjectName(QString::fromUtf8("checkInterceptionLAI"));
        checkInterceptionLAI->setGeometry(QRect(10, 10, 341, 17));
        checkInterceptionLAI->setFont(font1);
        checkInterceptionLAI->setChecked(true);
        groupInterception = new QGroupBox(page_2);
        groupInterception->setObjectName(QString::fromUtf8("groupInterception"));
        groupInterception->setGeometry(QRect(10, 40, 331, 181));
        gridLayout_5 = new QGridLayout(groupInterception);
        gridLayout_5->setContentsMargins(4, 4, 4, 4);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        radioButton_1 = new QRadioButton(groupInterception);
        buttonGroup_3 = new QButtonGroup(lisemqtClass);
        buttonGroup_3->setObjectName(QString::fromUtf8("buttonGroup_3"));
        buttonGroup_3->addButton(radioButton_1);
        radioButton_1->setObjectName(QString::fromUtf8("radioButton_1"));
        radioButton_1->setFont(font1);

        gridLayout_5->addWidget(radioButton_1, 0, 0, 1, 1);

        radioButton_2 = new QRadioButton(groupInterception);
        buttonGroup_3->addButton(radioButton_2);
        radioButton_2->setObjectName(QString::fromUtf8("radioButton_2"));
        radioButton_2->setFont(font1);

        gridLayout_5->addWidget(radioButton_2, 1, 0, 1, 1);

        radioButton_3 = new QRadioButton(groupInterception);
        buttonGroup_3->addButton(radioButton_3);
        radioButton_3->setObjectName(QString::fromUtf8("radioButton_3"));
        radioButton_3->setFont(font1);

        gridLayout_5->addWidget(radioButton_3, 2, 0, 1, 1);

        radioButton_4 = new QRadioButton(groupInterception);
        buttonGroup_3->addButton(radioButton_4);
        radioButton_4->setObjectName(QString::fromUtf8("radioButton_4"));
        radioButton_4->setFont(font1);

        gridLayout_5->addWidget(radioButton_4, 3, 0, 1, 1);

        radioButton_5 = new QRadioButton(groupInterception);
        buttonGroup_3->addButton(radioButton_5);
        radioButton_5->setObjectName(QString::fromUtf8("radioButton_5"));
        radioButton_5->setFont(font1);

        gridLayout_5->addWidget(radioButton_5, 4, 0, 1, 1);

        radioButton_6 = new QRadioButton(groupInterception);
        buttonGroup_3->addButton(radioButton_6);
        radioButton_6->setObjectName(QString::fromUtf8("radioButton_6"));
        radioButton_6->setFont(font1);

        gridLayout_5->addWidget(radioButton_6, 5, 0, 1, 1);

        radioButton_7 = new QRadioButton(groupInterception);
        buttonGroup_3->addButton(radioButton_7);
        radioButton_7->setObjectName(QString::fromUtf8("radioButton_7"));
        radioButton_7->setFont(font1);

        gridLayout_5->addWidget(radioButton_7, 6, 0, 1, 1);

        radioButton_8 = new QRadioButton(groupInterception);
        buttonGroup_3->addButton(radioButton_8);
        radioButton_8->setObjectName(QString::fromUtf8("radioButton_8"));
        radioButton_8->setFont(font1);

        gridLayout_5->addWidget(radioButton_8, 7, 0, 1, 1);

        toolBox->addItem(page_2, QString::fromUtf8("Interception"));
        page = new QWidget();
        page->setObjectName(QString::fromUtf8("page"));
        page->setGeometry(QRect(0, 0, 389, 229));
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
        tabWidget_OutputMaps = new QTabWidget(tab);
        tabWidget_OutputMaps->setObjectName(QString::fromUtf8("tabWidget_OutputMaps"));
        tabWidget_OutputMaps->setEnabled(true);
        tabWidget_OutputMaps->setGeometry(QRect(420, 300, 291, 281));
        tab_4 = new QWidget();
        tab_4->setObjectName(QString::fromUtf8("tab_4"));
        layoutWidget = new QWidget(tab_4);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(10, 10, 281, 246));
        verticalLayout = new QVBoxLayout(layoutWidget);
        verticalLayout->setSpacing(4);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        checkBox_OutRunoff = new QCheckBox(layoutWidget);
        checkBox_OutRunoff->setObjectName(QString::fromUtf8("checkBox_OutRunoff"));

        verticalLayout->addWidget(checkBox_OutRunoff);

        checkBox_OutWH = new QCheckBox(layoutWidget);
        checkBox_OutWH->setObjectName(QString::fromUtf8("checkBox_OutWH"));

        verticalLayout->addWidget(checkBox_OutWH);

        checkBox_OutWHC = new QCheckBox(layoutWidget);
        checkBox_OutWHC->setObjectName(QString::fromUtf8("checkBox_OutWHC"));

        verticalLayout->addWidget(checkBox_OutWHC);

        checkBox_OutInf = new QCheckBox(layoutWidget);
        checkBox_OutInf->setObjectName(QString::fromUtf8("checkBox_OutInf"));

        verticalLayout->addWidget(checkBox_OutInf);

        checkBox_OutV = new QCheckBox(layoutWidget);
        checkBox_OutV->setObjectName(QString::fromUtf8("checkBox_OutV"));

        verticalLayout->addWidget(checkBox_OutV);

        checkBox_OutDet = new QCheckBox(layoutWidget);
        checkBox_OutDet->setObjectName(QString::fromUtf8("checkBox_OutDet"));

        verticalLayout->addWidget(checkBox_OutDet);

        checkBox_OutDep = new QCheckBox(layoutWidget);
        checkBox_OutDep->setObjectName(QString::fromUtf8("checkBox_OutDep"));

        verticalLayout->addWidget(checkBox_OutDep);

        checkBox_OutConc = new QCheckBox(layoutWidget);
        checkBox_OutConc->setObjectName(QString::fromUtf8("checkBox_OutConc"));

        verticalLayout->addWidget(checkBox_OutConc);

        checkBox_OutTC = new QCheckBox(layoutWidget);
        checkBox_OutTC->setObjectName(QString::fromUtf8("checkBox_OutTC"));

        verticalLayout->addWidget(checkBox_OutTC);

        checkBox_OutSurfStor = new QCheckBox(layoutWidget);
        checkBox_OutSurfStor->setObjectName(QString::fromUtf8("checkBox_OutSurfStor"));

        verticalLayout->addWidget(checkBox_OutSurfStor);

        checkBox_OutChanVol = new QCheckBox(layoutWidget);
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
        groupBox_4->setGeometry(QRect(720, 300, 101, 281));
        groupBox_4->setFlat(false);
        printinterval = new QSpinBox(groupBox_4);
        printinterval->setObjectName(QString::fromUtf8("printinterval"));
        printinterval->setGeometry(QRect(8, 44, 37, 22));
        printinterval->setMinimum(1);
        checkOutputTimeStep = new QRadioButton(groupBox_4);
        buttonGroup_2 = new QButtonGroup(lisemqtClass);
        buttonGroup_2->setObjectName(QString::fromUtf8("buttonGroup_2"));
        buttonGroup_2->addButton(checkOutputTimeStep);
        checkOutputTimeStep->setObjectName(QString::fromUtf8("checkOutputTimeStep"));
        checkOutputTimeStep->setGeometry(QRect(8, 16, 91, 31));
        checkOutputTimeStep->setChecked(true);
        checkOutputTimeUser = new QRadioButton(groupBox_4);
        buttonGroup_2->addButton(checkOutputTimeUser);
        checkOutputTimeUser->setObjectName(QString::fromUtf8("checkOutputTimeUser"));
        checkOutputTimeUser->setEnabled(true);
        checkOutputTimeUser->setGeometry(QRect(9, 72, 91, 31));
        label_16 = new QLabel(groupBox_4);
        label_16->setObjectName(QString::fromUtf8("label_16"));
        label_16->setGeometry(QRect(48, 48, 46, 13));
        plainTextEdit = new QPlainTextEdit(groupBox_4);
        plainTextEdit->setObjectName(QString::fromUtf8("plainTextEdit"));
        plainTextEdit->setEnabled(false);
        plainTextEdit->setGeometry(QRect(10, 100, 81, 171));
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

        E_EndTime = new QLineEdit(groupBoxTime);
        E_EndTime->setObjectName(QString::fromUtf8("E_EndTime"));
        E_EndTime->setFont(font1);

        gridLayout_3->addWidget(E_EndTime, 1, 1, 1, 1);

        label_8 = new QLabel(groupBoxTime);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setFont(font1);
        label_8->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_3->addWidget(label_8, 2, 0, 1, 1);

        E_Timestep = new QLineEdit(groupBoxTime);
        E_Timestep->setObjectName(QString::fromUtf8("E_Timestep"));
        E_Timestep->setFont(font1);

        gridLayout_3->addWidget(E_Timestep, 2, 1, 1, 1);

        groupBox = new QGroupBox(tab);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(10, 10, 391, 121));
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

        E_runFileList = new QComboBox(groupBox);
        E_runFileList->setObjectName(QString::fromUtf8("E_runFileList"));
        E_runFileList->setFont(font1);
        E_runFileList->setEditable(true);
        E_runFileList->setDuplicatesEnabled(false);

        gridLayout_2->addWidget(E_runFileList, 0, 1, 1, 1);

        toolButton_ShowRunfile = new QToolButton(groupBox);
        toolButton_ShowRunfile->setObjectName(QString::fromUtf8("toolButton_ShowRunfile"));
        toolButton_ShowRunfile->setIcon(icon2);

        gridLayout_2->addWidget(toolButton_ShowRunfile, 0, 2, 1, 1);

        toolButton_fileOpen = new QToolButton(groupBox);
        toolButton_fileOpen->setObjectName(QString::fromUtf8("toolButton_fileOpen"));
        toolButton_fileOpen->setIcon(icon1);

        gridLayout_2->addWidget(toolButton_fileOpen, 0, 3, 1, 1);

        groupOutputMain = new QGroupBox(tab);
        groupOutputMain->setObjectName(QString::fromUtf8("groupOutputMain"));
        groupOutputMain->setGeometry(QRect(420, 10, 401, 281));
        groupOutputMain->setFont(font3);
        label_17 = new QLabel(groupOutputMain);
        label_17->setObjectName(QString::fromUtf8("label_17"));
        label_17->setGeometry(QRect(5, 56, 81, 16));
        label_17->setFont(font1);
        label_17->setLayoutDirection(Qt::LeftToRight);
        label_17->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        E_DetachmentMap = new QLineEdit(groupOutputMain);
        E_DetachmentMap->setObjectName(QString::fromUtf8("E_DetachmentMap"));
        E_DetachmentMap->setGeometry(QRect(91, 54, 271, 20));
        E_DetachmentMap->setFont(font1);
        label_18 = new QLabel(groupOutputMain);
        label_18->setObjectName(QString::fromUtf8("label_18"));
        label_18->setGeometry(QRect(5, 81, 81, 16));
        label_18->setFont(font1);
        label_18->setLayoutDirection(Qt::LeftToRight);
        label_18->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        E_DepositionMap = new QLineEdit(groupOutputMain);
        E_DepositionMap->setObjectName(QString::fromUtf8("E_DepositionMap"));
        E_DepositionMap->setGeometry(QRect(91, 79, 271, 20));
        E_DepositionMap->setFont(font1);
        label_19 = new QLabel(groupOutputMain);
        label_19->setObjectName(QString::fromUtf8("label_19"));
        label_19->setGeometry(QRect(5, 106, 81, 16));
        label_19->setFont(font1);
        label_19->setLayoutDirection(Qt::LeftToRight);
        label_19->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        E_SoillossMap = new QLineEdit(groupOutputMain);
        E_SoillossMap->setObjectName(QString::fromUtf8("E_SoillossMap"));
        E_SoillossMap->setGeometry(QRect(91, 104, 271, 20));
        E_SoillossMap->setFont(font1);
        label_22 = new QLabel(groupOutputMain);
        label_22->setObjectName(QString::fromUtf8("label_22"));
        label_22->setGeometry(QRect(5, 131, 81, 16));
        label_22->setFont(font1);
        label_22->setLayoutDirection(Qt::LeftToRight);
        label_22->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        E_MainTotals = new QLineEdit(groupOutputMain);
        E_MainTotals->setObjectName(QString::fromUtf8("E_MainTotals"));
        E_MainTotals->setGeometry(QRect(91, 129, 271, 20));
        E_MainTotals->setFont(font1);
        E_PointResults = new QLineEdit(groupOutputMain);
        E_PointResults->setObjectName(QString::fromUtf8("E_PointResults"));
        E_PointResults->setGeometry(QRect(91, 154, 211, 20));
        E_PointResults->setFont(font1);
        label_20 = new QLabel(groupOutputMain);
        label_20->setObjectName(QString::fromUtf8("label_20"));
        label_20->setGeometry(QRect(5, 156, 81, 16));
        label_20->setFont(font1);
        label_20->setLayoutDirection(Qt::LeftToRight);
        label_20->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_3 = new QLabel(groupOutputMain);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(5, 22, 82, 16));
        label_3->setFont(font1);
        label_3->setLayoutDirection(Qt::RightToLeft);
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        E_ResultDir = new QLineEdit(groupOutputMain);
        E_ResultDir->setObjectName(QString::fromUtf8("E_ResultDir"));
        E_ResultDir->setGeometry(QRect(91, 21, 271, 20));
        QSizePolicy sizePolicy3(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(E_ResultDir->sizePolicy().hasHeightForWidth());
        E_ResultDir->setSizePolicy(sizePolicy3);
        E_ResultDir->setFont(font1);
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
        line_2->setGeometry(QRect(9, 40, 381, 16));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);
        checkSeparateOutput = new QCheckBox(groupOutputMain);
        checkSeparateOutput->setObjectName(QString::fromUtf8("checkSeparateOutput"));
        checkSeparateOutput->setGeometry(QRect(310, 155, 89, 17));
        checkSeparateOutput->setFont(font1);
        groupBox_2 = new QGroupBox(groupOutputMain);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setGeometry(QRect(320, 180, 71, 80));
        groupBox_2->setFont(font1);
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
        groupOutputFormat = new QGroupBox(groupOutputMain);
        groupOutputFormat->setObjectName(QString::fromUtf8("groupOutputFormat"));
        groupOutputFormat->setGeometry(QRect(10, 180, 301, 91));
        groupOutputFormat->setFont(font1);
        verticalLayout_2 = new QVBoxLayout(groupOutputFormat);
        verticalLayout_2->setContentsMargins(4, 4, 4, 4);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        frame4 = new QFrame(groupOutputFormat);
        frame4->setObjectName(QString::fromUtf8("frame4"));
        gridLayout_4 = new QGridLayout(frame4);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        checkSOBEKOutput = new QCheckBox(frame4);
        checkSOBEKOutput->setObjectName(QString::fromUtf8("checkSOBEKOutput"));
        sizePolicy.setHeightForWidth(checkSOBEKOutput->sizePolicy().hasHeightForWidth());
        checkSOBEKOutput->setSizePolicy(sizePolicy);

        gridLayout_4->addWidget(checkSOBEKOutput, 3, 0, 1, 1);

        SOBEKdatestring = new QLineEdit(frame4);
        SOBEKdatestring->setObjectName(QString::fromUtf8("SOBEKdatestring"));
        SOBEKdatestring->setEnabled(false);
        SOBEKdatestring->setMaxLength(10);
        SOBEKdatestring->setFrame(true);
        SOBEKdatestring->setCursorPosition(10);

        gridLayout_4->addWidget(SOBEKdatestring, 3, 1, 1, 1);

        checkWritePCRtimeplot = new QCheckBox(frame4);
        checkWritePCRtimeplot->setObjectName(QString::fromUtf8("checkWritePCRtimeplot"));

        gridLayout_4->addWidget(checkWritePCRtimeplot, 0, 0, 1, 2);

        checkWritePCRnames = new QCheckBox(frame4);
        checkWritePCRnames->setObjectName(QString::fromUtf8("checkWritePCRnames"));

        gridLayout_4->addWidget(checkWritePCRnames, 1, 0, 1, 2);


        verticalLayout_2->addWidget(frame4);

        label_49 = new QLabel(tab);
        label_49->setObjectName(QString::fromUtf8("label_49"));
        label_49->setGeometry(QRect(195, 150, 211, 20));
        tabWidget->addTab(tab, QString());
        tab_12 = new QWidget();
        tab_12->setObjectName(QString::fromUtf8("tab_12"));
        gridLayout_12 = new QGridLayout(tab_12);
        gridLayout_12->setObjectName(QString::fromUtf8("gridLayout_12"));
        groupBox_InputMaps = new QGroupBox(tab_12);
        groupBox_InputMaps->setObjectName(QString::fromUtf8("groupBox_InputMaps"));
        QSizePolicy sizePolicy5(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(groupBox_InputMaps->sizePolicy().hasHeightForWidth());
        groupBox_InputMaps->setSizePolicy(sizePolicy5);
        gridLayout_10 = new QGridLayout(groupBox_InputMaps);
        gridLayout_10->setObjectName(QString::fromUtf8("gridLayout_10"));
        treeView = new QTreeView(groupBox_InputMaps);
        treeView->setObjectName(QString::fromUtf8("treeView"));
        treeView->setFrameShape(QFrame::StyledPanel);
        treeView->setFrameShadow(QFrame::Sunken);
        treeView->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        treeView->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);

        gridLayout_10->addWidget(treeView, 0, 1, 1, 1);

        checkExpandActive = new QCheckBox(groupBox_InputMaps);
        checkExpandActive->setObjectName(QString::fromUtf8("checkExpandActive"));

        gridLayout_10->addWidget(checkExpandActive, 1, 1, 1, 1);


        gridLayout_12->addWidget(groupBox_InputMaps, 1, 0, 1, 1);

        label_47 = new QLabel(tab_12);
        label_47->setObjectName(QString::fromUtf8("label_47"));

        gridLayout_12->addWidget(label_47, 0, 0, 1, 1);

        tabWidget->addTab(tab_12, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        sedgroup = new QGroupBox(tab_3);
        sedgroup->setObjectName(QString::fromUtf8("sedgroup"));
        sedgroup->setGeometry(QRect(10, 220, 351, 111));
        sedgroup->setFont(font3);
        gridLayout_9 = new QGridLayout(sedgroup);
        gridLayout_9->setObjectName(QString::fromUtf8("gridLayout_9"));
        gridLayout_9->setHorizontalSpacing(6);
        gridLayout_9->setVerticalSpacing(2);
        gridLayout_9->setContentsMargins(6, 0, 6, 6);
        label_50 = new QLabel(sedgroup);
        label_50->setObjectName(QString::fromUtf8("label_50"));
        label_50->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_50, 0, 2, 1, 1);

        label_21 = new QLabel(sedgroup);
        label_21->setObjectName(QString::fromUtf8("label_21"));
        label_21->setFont(font1);
        label_21->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_21, 1, 0, 1, 1);

        label_splashdet = new QLabel(sedgroup);
        label_splashdet->setObjectName(QString::fromUtf8("label_splashdet"));
        label_splashdet->setFont(font1);
        label_splashdet->setFrameShape(QFrame::StyledPanel);
        label_splashdet->setFrameShadow(QFrame::Sunken);
        label_splashdet->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_splashdet, 1, 1, 1, 1);

        label_41 = new QLabel(sedgroup);
        label_41->setObjectName(QString::fromUtf8("label_41"));
        label_41->setFont(font1);
        label_41->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_41, 1, 2, 1, 1);

        label_detch = new QLabel(sedgroup);
        label_detch->setObjectName(QString::fromUtf8("label_detch"));
        label_detch->setFont(font1);
        label_detch->setFrameShape(QFrame::StyledPanel);
        label_detch->setFrameShadow(QFrame::Sunken);
        label_detch->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_detch, 1, 3, 1, 1);

        label_23 = new QLabel(sedgroup);
        label_23->setObjectName(QString::fromUtf8("label_23"));
        label_23->setFont(font1);
        label_23->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_23, 2, 0, 1, 1);

        label_flowdet = new QLabel(sedgroup);
        label_flowdet->setObjectName(QString::fromUtf8("label_flowdet"));
        label_flowdet->setFont(font1);
        label_flowdet->setFrameShape(QFrame::StyledPanel);
        label_flowdet->setFrameShadow(QFrame::Sunken);
        label_flowdet->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_flowdet, 2, 1, 1, 1);

        label_43 = new QLabel(sedgroup);
        label_43->setObjectName(QString::fromUtf8("label_43"));
        label_43->setFont(font1);
        label_43->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_43, 2, 2, 1, 1);

        label_depch = new QLabel(sedgroup);
        label_depch->setObjectName(QString::fromUtf8("label_depch"));
        label_depch->setFont(font1);
        label_depch->setFrameShape(QFrame::StyledPanel);
        label_depch->setFrameShadow(QFrame::Sunken);
        label_depch->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_depch, 2, 3, 1, 1);

        label_27 = new QLabel(sedgroup);
        label_27->setObjectName(QString::fromUtf8("label_27"));
        label_27->setFont(font1);
        label_27->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_27, 3, 0, 1, 1);

        label_dep = new QLabel(sedgroup);
        label_dep->setObjectName(QString::fromUtf8("label_dep"));
        label_dep->setFont(font1);
        label_dep->setFrameShape(QFrame::StyledPanel);
        label_dep->setFrameShadow(QFrame::Sunken);
        label_dep->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_dep, 3, 1, 1, 1);

        label_39 = new QLabel(sedgroup);
        label_39->setObjectName(QString::fromUtf8("label_39"));
        label_39->setFont(font1);
        label_39->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_39, 3, 2, 1, 1);

        label_sedvolch = new QLabel(sedgroup);
        label_sedvolch->setObjectName(QString::fromUtf8("label_sedvolch"));
        label_sedvolch->setFont(font1);
        label_sedvolch->setFrameShape(QFrame::StyledPanel);
        label_sedvolch->setFrameShadow(QFrame::Sunken);
        label_sedvolch->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_sedvolch, 3, 3, 1, 1);

        label_25 = new QLabel(sedgroup);
        label_25->setObjectName(QString::fromUtf8("label_25"));
        label_25->setFont(font1);
        label_25->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_25, 4, 0, 1, 1);

        label_sedvol = new QLabel(sedgroup);
        label_sedvol->setObjectName(QString::fromUtf8("label_sedvol"));
        label_sedvol->setFont(font1);
        label_sedvol->setFrameShape(QFrame::StyledPanel);
        label_sedvol->setFrameShadow(QFrame::Sunken);
        label_sedvol->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_sedvol, 4, 1, 1, 1);

        label_45 = new QLabel(sedgroup);
        label_45->setObjectName(QString::fromUtf8("label_45"));
        QPalette palette1;
        QBrush brush2(QColor(128, 128, 128, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette1.setBrush(QPalette::Active, QPalette::WindowText, brush2);
        palette1.setBrush(QPalette::Inactive, QPalette::WindowText, brush2);
        palette1.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_45->setPalette(palette1);
        label_45->setFont(font1);
        label_45->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_45, 4, 2, 1, 1);

        label_MBs = new QLabel(sedgroup);
        label_MBs->setObjectName(QString::fromUtf8("label_MBs"));
        QPalette palette2;
        palette2.setBrush(QPalette::Active, QPalette::WindowText, brush2);
        palette2.setBrush(QPalette::Inactive, QPalette::WindowText, brush2);
        palette2.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_MBs->setPalette(palette2);
        label_MBs->setFont(font1);
        label_MBs->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_MBs, 4, 3, 1, 1);

        label_51 = new QLabel(sedgroup);
        label_51->setObjectName(QString::fromUtf8("label_51"));
        label_51->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_9->addWidget(label_51, 0, 0, 1, 1);

        groupTime = new QGroupBox(tab_3);
        groupTime->setObjectName(QString::fromUtf8("groupTime"));
        groupTime->setGeometry(QRect(10, 0, 351, 76));
        groupTime->setFont(font3);
        gridLayout_13 = new QGridLayout(groupTime);
        gridLayout_13->setObjectName(QString::fromUtf8("gridLayout_13"));
        gridLayout_13->setHorizontalSpacing(6);
        gridLayout_13->setVerticalSpacing(2);
        gridLayout_13->setContentsMargins(6, 0, 6, 6);
        label_endtime = new QLabel(groupTime);
        label_endtime->setObjectName(QString::fromUtf8("label_endtime"));
        label_endtime->setFont(font1);
        label_endtime->setFrameShape(QFrame::StyledPanel);
        label_endtime->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_endtime, 0, 3, 1, 1);

        label_area = new QLabel(groupTime);
        label_area->setObjectName(QString::fromUtf8("label_area"));
        label_area->setFont(font1);
        label_area->setFrameShape(QFrame::StyledPanel);
        label_area->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_area, 2, 3, 1, 1);

        label_time = new QLabel(groupTime);
        label_time->setObjectName(QString::fromUtf8("label_time"));
        label_time->setFont(font1);
        label_time->setFrameShape(QFrame::StyledPanel);
        label_time->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_time, 0, 1, 1, 1);

        label_29 = new QLabel(groupTime);
        label_29->setObjectName(QString::fromUtf8("label_29"));
        label_29->setFont(font1);
        label_29->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_29, 0, 0, 1, 1);

        label_dx = new QLabel(groupTime);
        label_dx->setObjectName(QString::fromUtf8("label_dx"));
        label_dx->setFont(font1);
        label_dx->setFrameShape(QFrame::StyledPanel);
        label_dx->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_dx, 2, 1, 1, 1);

        label_46 = new QLabel(groupTime);
        label_46->setObjectName(QString::fromUtf8("label_46"));
        label_46->setFont(font1);
        label_46->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_46, 2, 0, 1, 1);

        label_38 = new QLabel(groupTime);
        label_38->setObjectName(QString::fromUtf8("label_38"));
        label_38->setFont(font1);
        label_38->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_38, 0, 2, 1, 1);

        label_48 = new QLabel(groupTime);
        label_48->setObjectName(QString::fromUtf8("label_48"));
        label_48->setFont(font1);
        label_48->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_48, 2, 2, 1, 1);

        label_30 = new QLabel(groupTime);
        label_30->setObjectName(QString::fromUtf8("label_30"));
        QPalette palette3;
        QBrush brush3(QColor(129, 129, 129, 255));
        brush3.setStyle(Qt::SolidPattern);
        palette3.setBrush(QPalette::Active, QPalette::WindowText, brush3);
        palette3.setBrush(QPalette::Inactive, QPalette::WindowText, brush3);
        palette3.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_30->setPalette(palette3);
        label_30->setFont(font1);
        label_30->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_30, 3, 0, 1, 2);

        label_endruntime = new QLabel(groupTime);
        label_endruntime->setObjectName(QString::fromUtf8("label_endruntime"));
        QPalette palette4;
        palette4.setBrush(QPalette::Active, QPalette::WindowText, brush3);
        palette4.setBrush(QPalette::Inactive, QPalette::WindowText, brush3);
        palette4.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_endruntime->setPalette(palette4);
        label_endruntime->setFont(font1);
        label_endruntime->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_endruntime, 3, 3, 1, 1);

        label_runtime = new QLabel(groupTime);
        label_runtime->setObjectName(QString::fromUtf8("label_runtime"));
        QPalette palette5;
        palette5.setBrush(QPalette::Active, QPalette::WindowText, brush3);
        palette5.setBrush(QPalette::Inactive, QPalette::WindowText, brush3);
        palette5.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_runtime->setPalette(palette5);
        label_runtime->setFont(font1);
        label_runtime->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_13->addWidget(label_runtime, 3, 2, 1, 1);

        watergroup = new QGroupBox(tab_3);
        watergroup->setObjectName(QString::fromUtf8("watergroup"));
        watergroup->setGeometry(QRect(10, 80, 351, 131));
        QSizePolicy sizePolicy6(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(watergroup->sizePolicy().hasHeightForWidth());
        watergroup->setSizePolicy(sizePolicy6);
        watergroup->setFont(font3);
        watergroup->setFlat(false);
        gridLayout_15 = new QGridLayout(watergroup);
        gridLayout_15->setSpacing(2);
        gridLayout_15->setObjectName(QString::fromUtf8("gridLayout_15"));
        gridLayout_15->setContentsMargins(6, 0, 6, 6);
        label_32 = new QLabel(watergroup);
        label_32->setObjectName(QString::fromUtf8("label_32"));
        label_32->setFont(font1);
        label_32->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_32, 0, 0, 1, 1);

        label_raintot = new QLabel(watergroup);
        label_raintot->setObjectName(QString::fromUtf8("label_raintot"));
        label_raintot->setFont(font1);
        label_raintot->setAutoFillBackground(false);
        label_raintot->setFrameShape(QFrame::StyledPanel);
        label_raintot->setFrameShadow(QFrame::Sunken);
        label_raintot->setLineWidth(1);
        label_raintot->setMidLineWidth(0);
        label_raintot->setScaledContents(false);
        label_raintot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_raintot->setMargin(0);

        gridLayout_15->addWidget(label_raintot, 0, 1, 1, 1);

        label_34 = new QLabel(watergroup);
        label_34->setObjectName(QString::fromUtf8("label_34"));
        label_34->setFont(font1);
        label_34->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_34, 1, 0, 1, 1);

        label_interctot = new QLabel(watergroup);
        label_interctot->setObjectName(QString::fromUtf8("label_interctot"));
        label_interctot->setFont(font1);
        label_interctot->setFrameShape(QFrame::StyledPanel);
        label_interctot->setFrameShadow(QFrame::Sunken);
        label_interctot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_interctot->setMargin(0);

        gridLayout_15->addWidget(label_interctot, 1, 1, 1, 1);

        label_36 = new QLabel(watergroup);
        label_36->setObjectName(QString::fromUtf8("label_36"));
        label_36->setFont(font1);
        label_36->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_36, 2, 0, 1, 1);

        label_infiltot = new QLabel(watergroup);
        label_infiltot->setObjectName(QString::fromUtf8("label_infiltot"));
        label_infiltot->setFont(font1);
        label_infiltot->setFrameShape(QFrame::StyledPanel);
        label_infiltot->setFrameShadow(QFrame::Sunken);
        label_infiltot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_infiltot->setMargin(0);

        gridLayout_15->addWidget(label_infiltot, 2, 1, 1, 1);

        label_37 = new QLabel(watergroup);
        label_37->setObjectName(QString::fromUtf8("label_37"));
        label_37->setFont(font1);
        label_37->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_37, 3, 0, 1, 1);

        label_surfstor = new QLabel(watergroup);
        label_surfstor->setObjectName(QString::fromUtf8("label_surfstor"));
        label_surfstor->setFont(font1);
        label_surfstor->setFrameShape(QFrame::StyledPanel);
        label_surfstor->setFrameShadow(QFrame::Sunken);
        label_surfstor->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_surfstor->setMargin(0);

        gridLayout_15->addWidget(label_surfstor, 3, 1, 1, 1);

        label_40 = new QLabel(watergroup);
        label_40->setObjectName(QString::fromUtf8("label_40"));
        label_40->setFont(font1);
        label_40->setFrameShape(QFrame::NoFrame);
        label_40->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_40, 4, 0, 1, 1);

        label_watervoltot = new QLabel(watergroup);
        label_watervoltot->setObjectName(QString::fromUtf8("label_watervoltot"));
        label_watervoltot->setFont(font1);
        label_watervoltot->setFrameShape(QFrame::StyledPanel);
        label_watervoltot->setFrameShadow(QFrame::Sunken);
        label_watervoltot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_watervoltot->setMargin(0);

        gridLayout_15->addWidget(label_watervoltot, 4, 1, 1, 1);

        label_42 = new QLabel(watergroup);
        label_42->setObjectName(QString::fromUtf8("label_42"));
        label_42->setFont(font1);
        label_42->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_42, 5, 0, 1, 1);

        label_qtot = new QLabel(watergroup);
        label_qtot->setObjectName(QString::fromUtf8("label_qtot"));
        label_qtot->setFont(font1);
        label_qtot->setFrameShape(QFrame::StyledPanel);
        label_qtot->setFrameShadow(QFrame::Sunken);
        label_qtot->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_qtot->setMargin(0);

        gridLayout_15->addWidget(label_qtot, 5, 1, 1, 1);

        label_MB = new QLabel(watergroup);
        label_MB->setObjectName(QString::fromUtf8("label_MB"));
        QPalette palette6;
        palette6.setBrush(QPalette::Active, QPalette::WindowText, brush3);
        palette6.setBrush(QPalette::Inactive, QPalette::WindowText, brush3);
        palette6.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_MB->setPalette(palette6);
        label_MB->setFont(font1);
        label_MB->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_MB, 5, 3, 1, 1);

        label_44 = new QLabel(watergroup);
        label_44->setObjectName(QString::fromUtf8("label_44"));
        QPalette palette7;
        palette7.setBrush(QPalette::Active, QPalette::WindowText, brush3);
        palette7.setBrush(QPalette::Inactive, QPalette::WindowText, brush3);
        palette7.setBrush(QPalette::Disabled, QPalette::WindowText, brush1);
        label_44->setPalette(palette7);
        label_44->setFont(font1);
        label_44->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_44, 5, 2, 1, 1);

        label_24 = new QLabel(watergroup);
        label_24->setObjectName(QString::fromUtf8("label_24"));
        sizePolicy6.setHeightForWidth(label_24->sizePolicy().hasHeightForWidth());
        label_24->setSizePolicy(sizePolicy6);
        label_24->setFont(font1);
        label_24->setLayoutDirection(Qt::LeftToRight);
        label_24->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_24, 2, 2, 1, 1);

        label_qpeaktime = new QLabel(watergroup);
        label_qpeaktime->setObjectName(QString::fromUtf8("label_qpeaktime"));
        sizePolicy6.setHeightForWidth(label_qpeaktime->sizePolicy().hasHeightForWidth());
        label_qpeaktime->setSizePolicy(sizePolicy6);
        label_qpeaktime->setFont(font1);
        label_qpeaktime->setLayoutDirection(Qt::RightToLeft);
        label_qpeaktime->setFrameShape(QFrame::StyledPanel);
        label_qpeaktime->setFrameShadow(QFrame::Sunken);

        gridLayout_15->addWidget(label_qpeaktime, 2, 3, 1, 1);

        label_52 = new QLabel(watergroup);
        label_52->setObjectName(QString::fromUtf8("label_52"));
        label_52->setFont(font1);
        label_52->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_52, 0, 2, 1, 1);

        label_QPfrac = new QLabel(watergroup);
        label_QPfrac->setObjectName(QString::fromUtf8("label_QPfrac"));
        label_QPfrac->setFont(font1);
        label_QPfrac->setFrameShape(QFrame::StyledPanel);
        label_QPfrac->setFrameShadow(QFrame::Plain);
        label_QPfrac->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_QPfrac, 0, 3, 1, 1);

        label_35 = new QLabel(watergroup);
        label_35->setObjectName(QString::fromUtf8("label_35"));
        sizePolicy6.setHeightForWidth(label_35->sizePolicy().hasHeightForWidth());
        label_35->setSizePolicy(sizePolicy6);
        label_35->setFont(font1);
        label_35->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_35, 3, 2, 1, 1);

        label_qpeak = new QLabel(watergroup);
        label_qpeak->setObjectName(QString::fromUtf8("label_qpeak"));
        sizePolicy6.setHeightForWidth(label_qpeak->sizePolicy().hasHeightForWidth());
        label_qpeak->setSizePolicy(sizePolicy6);
        label_qpeak->setFont(font1);
        label_qpeak->setFrameShape(QFrame::StyledPanel);
        label_qpeak->setFrameShadow(QFrame::Sunken);
        label_qpeak->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_qpeak, 3, 3, 1, 1);

        label_53 = new QLabel(watergroup);
        label_53->setObjectName(QString::fromUtf8("label_53"));
        label_53->setFont(font1);
        label_53->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_53, 1, 2, 1, 1);

        label_ppeaktime = new QLabel(watergroup);
        label_ppeaktime->setObjectName(QString::fromUtf8("label_ppeaktime"));
        label_ppeaktime->setFont(font1);
        label_ppeaktime->setFrameShape(QFrame::StyledPanel);
        label_ppeaktime->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_15->addWidget(label_ppeaktime, 1, 3, 1, 1);

        progressBar = new QProgressBar(tab_3);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(10, 560, 811, 21));
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignCenter);
        label_debug = new QLabel(tab_3);
        label_debug->setObjectName(QString::fromUtf8("label_debug"));
        label_debug->setGeometry(QRect(10, 470, 351, 81));
        label_debug->setFrameShape(QFrame::Panel);
        label_debug->setFrameShadow(QFrame::Sunken);
        label_debug->setWordWrap(true);
        widgetGraph = new QWidget(tab_3);
        widgetGraph->setObjectName(QString::fromUtf8("widgetGraph"));
        widgetGraph->setGeometry(QRect(370, 10, 451, 401));
        QFont font4;
        font4.setFamily(QString::fromUtf8("Calibri"));
        font4.setPointSize(10);
        widgetGraph->setFont(font4);
        textGraph = new QPlainTextEdit(tab_3);
        textGraph->setObjectName(QString::fromUtf8("textGraph"));
        textGraph->setGeometry(QRect(400, 440, 391, 111));
        outletgroup_2 = new QGroupBox(tab_3);
        outletgroup_2->setObjectName(QString::fromUtf8("outletgroup_2"));
        outletgroup_2->setGeometry(QRect(10, 340, 351, 38));
        outletgroup_2->setFont(font3);
        gridLayout_16 = new QGridLayout(outletgroup_2);
        gridLayout_16->setSpacing(6);
        gridLayout_16->setObjectName(QString::fromUtf8("gridLayout_16"));
        gridLayout_16->setContentsMargins(6, 0, 6, 6);
        label_buffervol = new QLabel(outletgroup_2);
        label_buffervol->setObjectName(QString::fromUtf8("label_buffervol"));
        label_buffervol->setFont(font1);
        label_buffervol->setFrameShape(QFrame::StyledPanel);
        label_buffervol->setFrameShadow(QFrame::Plain);
        label_buffervol->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_16->addWidget(label_buffervol, 0, 1, 1, 1);

        label_56 = new QLabel(outletgroup_2);
        label_56->setObjectName(QString::fromUtf8("label_56"));
        label_56->setFont(font1);
        label_56->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_16->addWidget(label_56, 0, 0, 1, 1);

        label_buffersed = new QLabel(outletgroup_2);
        label_buffersed->setObjectName(QString::fromUtf8("label_buffersed"));
        label_buffersed->setFont(font1);
        label_buffersed->setFrameShape(QFrame::StyledPanel);
        label_buffersed->setFrameShadow(QFrame::Sunken);
        label_buffersed->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_16->addWidget(label_buffersed, 0, 3, 1, 1);

        label_57 = new QLabel(outletgroup_2);
        label_57->setObjectName(QString::fromUtf8("label_57"));
        label_57->setFont(font1);
        label_57->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_16->addWidget(label_57, 0, 2, 1, 1);

        outletgroup = new QGroupBox(tab_3);
        outletgroup->setObjectName(QString::fromUtf8("outletgroup"));
        outletgroup->setGeometry(QRect(10, 380, 351, 56));
        outletgroup->setFont(font3);
        gridLayout_14 = new QGridLayout(outletgroup);
        gridLayout_14->setObjectName(QString::fromUtf8("gridLayout_14"));
        gridLayout_14->setHorizontalSpacing(6);
        gridLayout_14->setVerticalSpacing(2);
        gridLayout_14->setContentsMargins(6, 0, 6, 6);
        label_26 = new QLabel(outletgroup);
        label_26->setObjectName(QString::fromUtf8("label_26"));
        sizePolicy6.setHeightForWidth(label_26->sizePolicy().hasHeightForWidth());
        label_26->setSizePolicy(sizePolicy6);
        label_26->setFont(font1);
        label_26->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_26, 1, 0, 1, 1);

        label_qtotm3 = new QLabel(outletgroup);
        label_qtotm3->setObjectName(QString::fromUtf8("label_qtotm3"));
        sizePolicy6.setHeightForWidth(label_qtotm3->sizePolicy().hasHeightForWidth());
        label_qtotm3->setSizePolicy(sizePolicy6);
        label_qtotm3->setMinimumSize(QSize(0, 15));
        label_qtotm3->setFont(font1);
        label_qtotm3->setFrameShape(QFrame::StyledPanel);
        label_qtotm3->setFrameShadow(QFrame::Sunken);
        label_qtotm3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_qtotm3, 1, 1, 1, 1);

        label_28 = new QLabel(outletgroup);
        label_28->setObjectName(QString::fromUtf8("label_28"));
        label_28->setFont(font1);
        label_28->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_28, 1, 2, 1, 1);

        label_soilloss = new QLabel(outletgroup);
        label_soilloss->setObjectName(QString::fromUtf8("label_soilloss"));
        label_soilloss->setFont(font1);
        label_soilloss->setFrameShape(QFrame::StyledPanel);
        label_soilloss->setFrameShadow(QFrame::Sunken);
        label_soilloss->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_soilloss, 1, 3, 1, 1);

        label_discharge = new QLabel(outletgroup);
        label_discharge->setObjectName(QString::fromUtf8("label_discharge"));
        label_discharge->setFont(font1);
        label_discharge->setFrameShape(QFrame::StyledPanel);
        label_discharge->setFrameShadow(QFrame::Plain);
        label_discharge->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_discharge, 0, 1, 1, 1);

        label_54 = new QLabel(outletgroup);
        label_54->setObjectName(QString::fromUtf8("label_54"));
        label_54->setFont(font1);
        label_54->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_54, 0, 0, 1, 1);

        label_soillosskgha = new QLabel(outletgroup);
        label_soillosskgha->setObjectName(QString::fromUtf8("label_soillosskgha"));
        label_soillosskgha->setFont(font1);
        label_soillosskgha->setFrameShape(QFrame::StyledPanel);
        label_soillosskgha->setFrameShadow(QFrame::Sunken);
        label_soillosskgha->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_soillosskgha, 0, 3, 1, 1);

        label_31 = new QLabel(outletgroup);
        label_31->setObjectName(QString::fromUtf8("label_31"));
        label_31->setFont(font1);
        label_31->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_14->addWidget(label_31, 0, 2, 1, 1);

        tabWidget->addTab(tab_3, QString());

        gridLayout->addWidget(tabWidget, 0, 0, 1, 1);

        lisemqtClass->setCentralWidget(centralwidget);
        toolBar = new QToolBar(lisemqtClass);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
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
        QWidget::setTabOrder(E_BeginTime, E_Timestep);
        QWidget::setTabOrder(E_Timestep, E_InfiltrationMethod);
        QWidget::setTabOrder(E_InfiltrationMethod, checkInfilCrust);
        QWidget::setTabOrder(checkInfilCrust, checkInfilCompact);
        QWidget::setTabOrder(checkInfilCompact, checkImpermeable);
        QWidget::setTabOrder(checkImpermeable, E_EndTime);

        retranslateUi(lisemqtClass);
        QObject::connect(checkIncludeChannel, SIGNAL(toggled(bool)), checkChannelInfil, SLOT(setEnabled(bool)));
        QObject::connect(checkIncludeChannel, SIGNAL(toggled(bool)), checkChannelBaseflow, SLOT(setEnabled(bool)));
        QObject::connect(checkChannelInfil, SIGNAL(toggled(bool)), checkChannelBaseflow, SLOT(setDisabled(bool)));
        QObject::connect(checkChannelBaseflow, SIGNAL(toggled(bool)), checkChannelInfil, SLOT(setDisabled(bool)));
        QObject::connect(checkSnowmelt, SIGNAL(toggled(bool)), E_SnowmeltName, SLOT(setEnabled(bool)));
        QObject::connect(checkSnowmelt, SIGNAL(toggled(bool)), label_5, SLOT(setEnabled(bool)));
        QObject::connect(checkSnowmelt, SIGNAL(toggled(bool)), toolButton_SnowmeltShow, SLOT(setEnabled(bool)));
        QObject::connect(groupBox, SIGNAL(toggled(bool)), toolButton_SnowmeltNameS, SLOT(setEnabled(bool)));
        QObject::connect(checkSnowmelt, SIGNAL(toggled(bool)), toolButton_SnowmeltNameS, SLOT(setEnabled(bool)));
        QObject::connect(checkNoErosion, SIGNAL(toggled(bool)), checkAltErosion, SLOT(setDisabled(bool)));
        QObject::connect(checkOutputTimeStep, SIGNAL(toggled(bool)), printinterval, SLOT(setEnabled(bool)));
        QObject::connect(checkOutputTimeStep, SIGNAL(toggled(bool)), label_16, SLOT(setEnabled(bool)));
        QObject::connect(checkInfilGrass, SIGNAL(toggled(bool)), label_15, SLOT(setEnabled(bool)));
        QObject::connect(checkInfilGrass, SIGNAL(toggled(bool)), E_GrassStripN, SLOT(setEnabled(bool)));
        QObject::connect(checkOutputTimeUser, SIGNAL(toggled(bool)), plainTextEdit, SLOT(setEnabled(bool)));
        QObject::connect(checkInterceptionLAI, SIGNAL(toggled(bool)), groupInterception, SLOT(setEnabled(bool)));

        tabWidget->setCurrentIndex(2);
        toolBox->setCurrentIndex(2);
        toolBox->layout()->setSpacing(1);
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
        checkSimpleDepression->setText(QApplication::translate("lisemqtClass", "Simple surface storage (no runoff below max storage threshold)", 0, QApplication::UnicodeUTF8));
        checkHardsurface->setText(QApplication::translate("lisemqtClass", "Include hard surfaces (no interc, infil, detach, needs hardsurf.map)", 0, QApplication::UnicodeUTF8));
        checkChannelInfil->setText(QApplication::translate("lisemqtClass", "Channel infiltration (needs aditional map)", 0, QApplication::UnicodeUTF8));
        checkChannelBaseflow->setText(QApplication::translate("lisemqtClass", "Channel baseflow (needs additional maps)", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page_3), QApplication::translate("lisemqtClass", "Global Options", 0, QApplication::UnicodeUTF8));
        checkInfilCompact->setText(QApplication::translate("lisemqtClass", "Include compacted areas", 0, QApplication::UnicodeUTF8));
        checkInfilCrust->setText(QApplication::translate("lisemqtClass", "Include crusts", 0, QApplication::UnicodeUTF8));
        checkImpermeable->setText(QApplication::translate("lisemqtClass", "Impermeable lower soil boundary", 0, QApplication::UnicodeUTF8));
        checkInfil2layer->setText(QApplication::translate("lisemqtClass", "Include 2nd soil layer", 0, QApplication::UnicodeUTF8));
        groupBox_SwatreOptions->setTitle(QApplication::translate("lisemqtClass", "Swatre", 0, QApplication::UnicodeUTF8));
        checkGeometricMean->setText(QApplication::translate("lisemqtClass", "geometric average Ksat", 0, QApplication::UnicodeUTF8));
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
        label_15->setText(QApplication::translate("lisemqtClass", "Mannings n grass strip", 0, QApplication::UnicodeUTF8));
        E_GrassStripN->setInputMask(QApplication::translate("lisemqtClass", "9.99; ", 0, QApplication::UnicodeUTF8));
        E_GrassStripN->setText(QApplication::translate("lisemqtClass", "0.30", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page_5), QApplication::translate("lisemqtClass", "Conservation measures", 0, QApplication::UnicodeUTF8));
        checkInterceptionLAI->setText(QApplication::translate("lisemqtClass", "Use storage equation below (else use Smax map directly)", 0, QApplication::UnicodeUTF8));
        groupInterception->setTitle(QApplication::translate("lisemqtClass", "Storage equation", 0, QApplication::UnicodeUTF8));
        radioButton_1->setText(QApplication::translate("lisemqtClass", "Orig. LISEM (crops): S = 0.935+0.498*LAI-0.00575*LAI^2", 0, QApplication::UnicodeUTF8));
        radioButton_2->setText(QApplication::translate("lisemqtClass", "Pinus:         S = 0.2331*LAI          (n=12,R2=0.88)", 0, QApplication::UnicodeUTF8));
        radioButton_3->setText(QApplication::translate("lisemqtClass", "Douglas Fir:   S = 0.3165*LAI          (n=4, R2=0.83)", 0, QApplication::UnicodeUTF8));
        radioButton_4->setText(QApplication::translate("lisemqtClass", "Olive:         S = 1.46 * LAI^0.56     (n=5, R2=0.87)", 0, QApplication::UnicodeUTF8));
        radioButton_5->setText(QApplication::translate("lisemqtClass", "Eucalypt:      S = 0.0918*LAI^1.04     (n=8, R2=0.51)", 0, QApplication::UnicodeUTF8));
        radioButton_6->setText(QApplication::translate("lisemqtClass", "Rainforest:    S = 0.2856*LAI          (n=5, R2=0.60)", 0, QApplication::UnicodeUTF8));
        radioButton_7->setText(QApplication::translate("lisemqtClass", "Bracken:       S = 0.1713*LAI          (n=8, R2=0.98)", 0, QApplication::UnicodeUTF8));
        radioButton_8->setText(QApplication::translate("lisemqtClass", "Clumped grass: S = 0.59 * LAI^0.88     (n=6, R2=0.82)", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page_2), QApplication::translate("lisemqtClass", "Interception", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("lisemqtClass", "Multiplication factor Ksat slopes", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("lisemqtClass", "Multiplication factor Manning N slopes", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("lisemqtClass", "Multiplication factor Ksat Channel", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("lisemqtClass", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">Multiplication factor Manning N Channel</span></p></body></html>", 0, QApplication::UnicodeUTF8));
        toolBox->setItemText(toolBox->indexOf(page), QApplication::translate("lisemqtClass", "Calibration", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkBox_OutRunoff->setToolTip(QApplication::translate("lisemqtClass", "mapseries saved as is \"ro\"", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        checkBox_OutRunoff->setText(QApplication::translate("lisemqtClass", "ro - runoff slope & channel (l/s)", 0, QApplication::UnicodeUTF8));
        checkBox_OutWH->setText(QApplication::translate("lisemqtClass", "wh - water depth on surface (mm)", 0, QApplication::UnicodeUTF8));
        checkBox_OutWHC->setText(QApplication::translate("lisemqtClass", "whc - cumulative runof depth (mm", 0, QApplication::UnicodeUTF8));
        checkBox_OutInf->setText(QApplication::translate("lisemqtClass", "inf - cumulative infilktration (mm)", 0, QApplication::UnicodeUTF8));
        checkBox_OutV->setText(QApplication::translate("lisemqtClass", "velo - velocity (m/s)", 0, QApplication::UnicodeUTF8));
        checkBox_OutDet->setText(QApplication::translate("lisemqtClass", "det - detachment (selected units)", 0, QApplication::UnicodeUTF8));
        checkBox_OutDep->setText(QApplication::translate("lisemqtClass", "dep - deposition (select units)", 0, QApplication::UnicodeUTF8));
        checkBox_OutConc->setText(QApplication::translate("lisemqtClass", "conc - sediment concentration (g/l)", 0, QApplication::UnicodeUTF8));
        checkBox_OutTC->setText(QApplication::translate("lisemqtClass", "tc - transport capacity (kg/m3)", 0, QApplication::UnicodeUTF8));
        checkBox_OutSurfStor->setText(QApplication::translate("lisemqtClass", "sstor - surface storage (mm)", 0, QApplication::UnicodeUTF8));
        checkBox_OutChanVol->setText(QApplication::translate("lisemqtClass", "chanvol - channel volume (m3)", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_4), QApplication::translate("lisemqtClass", "Basic Output", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_6), QApplication::translate("lisemqtClass", "Multiclass Output", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_5), QApplication::translate("lisemqtClass", "Nutrient Output", 0, QApplication::UnicodeUTF8));
        tabWidget_OutputMaps->setTabText(tabWidget_OutputMaps->indexOf(tab_7), QApplication::translate("lisemqtClass", "Gully Output", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("lisemqtClass", "Output times", 0, QApplication::UnicodeUTF8));
        checkOutputTimeStep->setText(QApplication::translate("lisemqtClass", "Report every:", 0, QApplication::UnicodeUTF8));
        checkOutputTimeUser->setText(QApplication::translate("lisemqtClass", "User defined:", 0, QApplication::UnicodeUTF8));
        label_16->setText(QApplication::translate("lisemqtClass", "timesteps", 0, QApplication::UnicodeUTF8));
        plainTextEdit->setPlainText(QString());
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
        toolButton_fileOpen->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        groupOutputMain->setTitle(QApplication::translate("lisemqtClass", "Output", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("lisemqtClass", "Detachment Map", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("lisemqtClass", "Deposition Map", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("lisemqtClass", "Soilloss Map", 0, QApplication::UnicodeUTF8));
        label_22->setText(QApplication::translate("lisemqtClass", "Text file Totals", 0, QApplication::UnicodeUTF8));
        label_20->setText(QApplication::translate("lisemqtClass", "Point output", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("lisemqtClass", "Result directory  ", 0, QApplication::UnicodeUTF8));
        toolButton_ResultDir->setText(QApplication::translate("lisemqtClass", "...", 0, QApplication::UnicodeUTF8));
        checkSeparateOutput->setText(QApplication::translate("lisemqtClass", "Separate files", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("lisemqtClass", "Map units", 0, QApplication::UnicodeUTF8));
        checkUnits_kgcell->setText(QApplication::translate("lisemqtClass", "kg/cell", 0, QApplication::UnicodeUTF8));
        checkUnits_kgm2->setText(QApplication::translate("lisemqtClass", "kg/m2", 0, QApplication::UnicodeUTF8));
        checkUnits_tonha->setText(QApplication::translate("lisemqtClass", "ton/ha", 0, QApplication::UnicodeUTF8));
        groupOutputFormat->setTitle(QApplication::translate("lisemqtClass", "Output format", 0, QApplication::UnicodeUTF8));
        checkSOBEKOutput->setText(QApplication::translate("lisemqtClass", "Hydrographs SOBEK format, date string:", 0, QApplication::UnicodeUTF8));
        SOBEKdatestring->setInputMask(QApplication::translate("lisemqtClass", "99\\\\99\\\\9999; ", 0, QApplication::UnicodeUTF8));
        SOBEKdatestring->setText(QApplication::translate("lisemqtClass", "01\\01\\2000", 0, QApplication::UnicodeUTF8));
        checkWritePCRtimeplot->setText(QApplication::translate("lisemqtClass", "hydrographs PCRaster format (else comma delimited)", 0, QApplication::UnicodeUTF8));
        checkWritePCRnames->setText(QApplication::translate("lisemqtClass", "Mapseries as PCRaster filenames (.001, .002)", 0, QApplication::UnicodeUTF8));
        label_49->setText(QApplication::translate("lisemqtClass", "TextLabel", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("lisemqtClass", "Input/Output", 0, QApplication::UnicodeUTF8));
        groupBox_InputMaps->setTitle(QApplication::translate("lisemqtClass", "Input Maps", 0, QApplication::UnicodeUTF8));
        checkExpandActive->setText(QApplication::translate("lisemqtClass", "Expand activated maps categories", 0, QApplication::UnicodeUTF8));
        label_47->setText(QApplication::translate("lisemqtClass", "TextLabel", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_12), QApplication::translate("lisemqtClass", "Maps", 0, QApplication::UnicodeUTF8));
        sedgroup->setTitle(QApplication::translate("lisemqtClass", "Sediment", 0, QApplication::UnicodeUTF8));
        label_50->setText(QApplication::translate("lisemqtClass", "Channel", 0, QApplication::UnicodeUTF8));
        label_21->setText(QApplication::translate("lisemqtClass", "Splash (ton)", 0, QApplication::UnicodeUTF8));
        label_splashdet->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_41->setText(QApplication::translate("lisemqtClass", "Detachment (ton)", 0, QApplication::UnicodeUTF8));
        label_detch->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_23->setText(QApplication::translate("lisemqtClass", "Flow detachment (ton)", 0, QApplication::UnicodeUTF8));
        label_flowdet->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_43->setText(QApplication::translate("lisemqtClass", "Deposition (ton)", 0, QApplication::UnicodeUTF8));
        label_depch->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_27->setText(QApplication::translate("lisemqtClass", "Deposition (ton)", 0, QApplication::UnicodeUTF8));
        label_dep->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_39->setText(QApplication::translate("lisemqtClass", "Sed in flow (ton)", 0, QApplication::UnicodeUTF8));
        label_sedvolch->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_25->setText(QApplication::translate("lisemqtClass", "Sediment in flow (ton)", 0, QApplication::UnicodeUTF8));
        label_sedvol->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_45->setText(QApplication::translate("lisemqtClass", "MBs (%)", 0, QApplication::UnicodeUTF8));
        label_MBs->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_51->setText(QApplication::translate("lisemqtClass", "Slope", 0, QApplication::UnicodeUTF8));
        groupTime->setTitle(QApplication::translate("lisemqtClass", "Time", 0, QApplication::UnicodeUTF8));
        label_endtime->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_area->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_time->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_29->setText(QApplication::translate("lisemqtClass", "                     time (min)", 0, QApplication::UnicodeUTF8));
        label_dx->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_46->setText(QApplication::translate("lisemqtClass", "Cellsize (m)", 0, QApplication::UnicodeUTF8));
        label_38->setText(QApplication::translate("lisemqtClass", "   End time (min)", 0, QApplication::UnicodeUTF8));
        label_48->setText(QApplication::translate("lisemqtClass", "ha", 0, QApplication::UnicodeUTF8));
        label_30->setText(QApplication::translate("lisemqtClass", "model runtime & est.end (min)", 0, QApplication::UnicodeUTF8));
        label_endruntime->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_runtime->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        watergroup->setTitle(QApplication::translate("lisemqtClass", "Water", 0, QApplication::UnicodeUTF8));
        label_32->setText(QApplication::translate("lisemqtClass", "Rain tot (mm)", 0, QApplication::UnicodeUTF8));
        label_raintot->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_34->setText(QApplication::translate("lisemqtClass", "Interception (mm)", 0, QApplication::UnicodeUTF8));
        label_interctot->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_36->setText(QApplication::translate("lisemqtClass", "Infil (mm)", 0, QApplication::UnicodeUTF8));
        label_infiltot->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_37->setText(QApplication::translate("lisemqtClass", "Surface Store (mm)", 0, QApplication::UnicodeUTF8));
        label_surfstor->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_40->setText(QApplication::translate("lisemqtClass", "Runoff (incl chan) (mm)", 0, QApplication::UnicodeUTF8));
        label_watervoltot->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_42->setText(QApplication::translate("lisemqtClass", "Qtot (mm)", 0, QApplication::UnicodeUTF8));
        label_qtot->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_MB->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_44->setText(QApplication::translate("lisemqtClass", "MB (%)", 0, QApplication::UnicodeUTF8));
        label_24->setText(QApplication::translate("lisemqtClass", "peak time Q (min)", 0, QApplication::UnicodeUTF8));
        label_qpeaktime->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_52->setText(QApplication::translate("lisemqtClass", "Q/P (%)", 0, QApplication::UnicodeUTF8));
        label_QPfrac->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_35->setText(QApplication::translate("lisemqtClass", "Qpeak (l/s)", 0, QApplication::UnicodeUTF8));
        label_qpeak->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_53->setText(QApplication::translate("lisemqtClass", "peak time P (min)", 0, QApplication::UnicodeUTF8));
        label_ppeaktime->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        progressBar->setFormat(QApplication::translate("lisemqtClass", "%p% (%v-%m)", 0, QApplication::UnicodeUTF8));
        label_debug->setText(QApplication::translate("lisemqtClass", "!!!", 0, QApplication::UnicodeUTF8));
        outletgroup_2->setTitle(QApplication::translate("lisemqtClass", "Buffers", 0, QApplication::UnicodeUTF8));
        label_buffervol->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_56->setText(QApplication::translate("lisemqtClass", "                   Water (m3)", 0, QApplication::UnicodeUTF8));
        label_buffersed->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_57->setText(QApplication::translate("lisemqtClass", "     Sediment (ton)", 0, QApplication::UnicodeUTF8));
        outletgroup->setTitle(QApplication::translate("lisemqtClass", "Outlet", 0, QApplication::UnicodeUTF8));
        label_26->setText(QApplication::translate("lisemqtClass", "Qtot (m3)", 0, QApplication::UnicodeUTF8));
        label_qtotm3->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_28->setText(QApplication::translate("lisemqtClass", "Soil loss (ton)", 0, QApplication::UnicodeUTF8));
        label_soilloss->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_discharge->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_54->setText(QApplication::translate("lisemqtClass", "                           Q (l/s)", 0, QApplication::UnicodeUTF8));
        label_soillosskgha->setText(QApplication::translate("lisemqtClass", "0", 0, QApplication::UnicodeUTF8));
        label_31->setText(QApplication::translate("lisemqtClass", "    Soil loss (kg/ha)", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("lisemqtClass", "Simulation", 0, QApplication::UnicodeUTF8));
        toolBar->setWindowTitle(QApplication::translate("lisemqtClass", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class lisemqtClass: public Ui_lisemqtClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_LISEMQT_H
