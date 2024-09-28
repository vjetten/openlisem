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
 \file lisemqt.h
 \brief main interface header file containing UI class, linked to lisemqt.ui
 */

#ifndef LISEMQT_H
#define LISEMQT_H

#include <QtGui>
#include <QApplication>
#include <QtWidgets>
#include <QSystemTrayIcon>
#include <QTranslator>


#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_marker.h>
#include <qwt_legend.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_magnifier.h>
#include <qwt_color_map.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_plot_layout.h>
#include <qwt_matrix_raster_data.h>
#include <qwt_scale_widget.h>
#include <qwt_plot_renderer.h>
#include <qwt_plot_rescaler.h>
#include <qwt_scale_engine.h>
#include <qwt_plot_zoomer.h>
#include <qwt_picker.h>
#include <qwt_symbol.h>
#include <qwt_axis_id.h>

#include "version.h"
#include "ui_lisemqt.h"
#include "model.h"
#include "LisUIoutput.h"
#include "LisUItreemodel.h"
#include "LisUImapplot.h"
#include "lismpeg.h"



// constants to define the place of the main parts in the map tree structure
#define RAINFALLMAPS 0
#define CATCHMENTMAPS 1
#define LANDUSEMAPS 2
#define SURFACEMAPS 3
#define INFILTRATIONMAPS 4
#define CHANNELMAPS 5
#define HOUSESMAPS 6
#define EROSIONMAPS 7
#define CONSERVATIONMAPS 7
#define TILEDRAINMAPS 8
#define TILEMAPS 10
#define PESTMAPS 11


//---------------------------------------------------------------------------
/// map name list structure for interaction with interface
typedef struct MAP_LIST {
    QString name;
    QString value;
    QString dir;
    int groupnr;
    int varnr;
} MAP_LIST;


/// Exteneded interface class
class lisemqt : public QMainWindow, private Ui::lisemqtClass
{
    Q_OBJECT

public:
    lisemqt(QWidget *parent = 0, bool doBatch = false, QString runName = "");
    ~lisemqt();

    int genfontsize;
    double dpiscale;
    long nrValidCells;

    QDialog *helpbox;
    QTextEdit *helptxt;
    QHBoxLayout *helpLayout;

    QProgressBar *pb;

    lismpeg *lisMpeg;
    QString mencoderDir;

    bool darkLISEM;
    bool doBatchmode;
    QString batchRunname;

   // bool WhasStopped;

    void initMapTree();
    void DefaultMapnames();
    void fillMapnames();
    void updateDEFmaps();
    void fillNamelistMapnames(bool to);
    void checkMapNameModel(int parentrow, int selrow, bool setit);
    void SetToolBar();
    void GetStorePath();
    void StorePath();
    void SetStyleUI();
    void lightStyleUI();
    void darkStyleUI();
    void setOutputTabStyle(QString bc, QString fc);
    void SetStyleUISize();
    void GetRunfile();
    void ParseInputData();
    void updateModelData();
    void defaultRunFile();
    QString CheckDir(QString p, bool makeit = false);
    void RunAllChecks();
    void savefile(QString name);
    void SetConnections();
    QStringList runfilelist;


//    bool doNewPlot;
//    void newPlot(bool refresh);
//    void setupNewPlot();
//    void initNewPlot();

    QList <QPointF> dataRain;
    QList <QPointF> dataQ;
    QList <QPointF> dataQs;
    QList <QPointF> dataC;

    // structure output
    void initOP();

    // graph functions
    void setupPlot();
    void startPlots();
    void showPlot();
    void initPlot();

    void showOutputData();
    void showOutputDataZero();
    // Map drawing variable
    void setupMapPlot();
    void initMapPlot();

    void showMap();

    void showComboMap(int i);

    void showBaseMap();
    void getOutletMap();
    void showChannelVectorNew();
    void showRoadMap();
    void showHouseMap();
    void showHardSurfaceMap();
    void showImageMap();
    void changeSize();
    double Masp;
    double fillDrawMapData(cTMap *_M, double scale, QwtMatrixRasterData *_RD, double *minv, double *maxv);
    double fillDrawMapDataRGB(cTMap * base,  cTRGBMap *_M, QwtMatrixRasterData *_RD);

    QwtPlot *MPlot;               // plot in which the raster map is drawn
    QwtText title;
    //OSMPlot *OSMplot;
    QwtPlotSpectrogram *drawMap;  // raster map drawing
    QwtPlotSpectrogram *baseMap;  // raster map drawing
    QwtPlotSpectrogram *baseMapDEM;  // raster map drawing
    QwtPlotSpectrogram *baseMapImage;  // raster map drawing
    QwtPlotSpectrogram *contourDEM;  // raster map drawing
    QwtPlotSpectrogram *hardsurfMap;  // raster map drawing
    QwtPlotSpectrogram *roadMap;  // raster map drawing
    QwtPlotSpectrogram *houseMap;  // raster map drawing
    QwtPlotSpectrogram *imageMap;
    QwtPlotSpectrogram *outletMap;
    QwtMatrixRasterData *RD;      // data for thematic raster maps
    QwtMatrixRasterData *RDb;
    QwtMatrixRasterData *RDbb;
    QwtMatrixRasterData *RDc;
    QwtMatrixRasterData *RDd;
    QwtMatrixRasterData *RDe;
    QwtMatrixRasterData *RDf;
    QwtMatrixRasterData *RImage;
    QList<double> contourLevels;
    // QwtAxisId *axisYL1;
    // QwtAxisId *axisYL2;
    // QwtAxisId *axisYR1;
    // QwtAxisId *axisYR2;
    // QwtAxisId *axisX;
    QList <QVector <double>> Xa;
    QList <QVector <double>> Ya;
    QList <QVector <double>> Xc;
    QList <QVector <double>> Yc;
    QList<QwtPlotCurve*> rivers; // black river
    QList<QwtPlotCurve*> culverts;  //white culvert part in river
    QwtPlotCurve obspoints;
    QwtPlotCurve outlets;


    double contourmin, contourmax;
    // vars for store map display in runfile
    int MapDisplayMapSelection, MapDisplayBuilding, MapDisplayRoads, MapDisplayChannels, MapDisplayHydrology;
    double MapDisplayRunoffMax, MapDisplayInfiltrationMax, MapDisplaySoillossMax, MapDisplayFlooddepthMax;
    int MapDisplayIncludeRunoff, MapDisplayMinimumDepth, MapDisplayScreenDumps;

    QString E_WorkDir;

    QVector<double> mapData;
    QVector<double> RGBData;
    QwtInterval legend;
    QwtScaleWidget *rightAxis;
    QwtPlotRescaler *mapRescaler;
    double maxAxis1, maxAxis2, maxAxis3, maxAxis4, maxAxis5, maxAxis6, maxAxis7, maxAxis8;
    int pstep;
    QwtPlotMagnifier *magnifier;
    QwtPlotPanner *panner;
    QwtPlotZoomer* zoomer;
    MyPicker *picker;


    //hydrograph options
    QList<int> OutletIndices;
    QList<int> OutletLocationX;
    QList<int> OutletLocationY;
    double mult;

    int outletpoint = 1;

    //max va;ues for each outlet point
    double pmax;
    QVector<double> qmax;
    QVector<double> qsmax;
    QVector<double> cmax;

    //Map display options

    QList <QwtComboColorMap*> cmMap;
    QList <QwtComboColorMap*> cmLeg;
    QList <QwtComboColorMap*> ColorMapList;
    QList<QString> NameList;
    QList<QString> UnitList;
    QList<bool> SymList;
    QList<bool> LogList;
    QList<int> IndexList;
    QList<int> IndexList1;
    int ActiveList = 0;

    // graph variables
    QwtPlot *HPlot;
    QwtPlotCurve *QGraph;
    QwtPlotCurve *QbGraph;
    QwtPlotCurve *QsGraph;
    QwtPlotCurve *CGraph;
    QwtPlotCurve *PGraph;
    QwtPlotCurve *QtileGraph;

    bool startplot;
    bool stoprun;
    QVector <double> times;
    int lastOptionSceen;

    bool oldRunfile; // check is old runfile for ksat calibration
    bool saveRunFileOnce; // check is old runfile for ksat calibration

    //interface names
    TreeModel *MapNameModel;
    QString currentDir;
    QString RainFileName;
    QString RainFileDir;
    QString RainSatFileName;
    QString RainSatFileDir;
    QString ETFileName;
    QString ETFileDir;
    QString ETSatFileName;
    QString ETSatFileDir;
    QString DischargeinDir;
    QString DischargeinFileName;
    QString WaveinDir;
    QString WaveinFileName;
    QString SnowmeltFileName;
    QString SnowmeltFileDir;
    QString SwatreTableName;
    QString SwatreTableDir;
    QStringList DEFmaps;
    //QStringList RunFileNames;
    QString satImageFileDir;
    QString satImageFileName;
    int CurrentRunFile;
    int uiInfilMethod;
    double swatreDT;
    QString screenShotDir;
    bool doChannelBaseflow;

 //   QList<cTMap *> ComboMapsSafe;

    MAP_LIST mapList[NUMMAPS]; /// structure for current map names, can be edited by user
    int nrmaplist;
    NAME_LIST namelist[NUMNAMES]; /// structure to read all runfile variables and names
    int nrnamelist;
    QStringList outputcheck; /// list of '0' and '1' to see which output mapseries are checled by the nuser
    int InterceptionEqNr;
    int mapstartnr;
    bool doShootScreens, startShootScreens;
    QString mapFormat;
    int cpucores;

    // buttongroups to make checkboxes mutually exclusive
    QButtonGroup GroupMapDisplay;
    QButtonGroup GroupImpermeable;
    QButtonGroup GroupBaseflow;
    QButtonGroup GroupRunoff;
    QButtonGroup GroupChannel;

    void setDisplayComboBoxes();
    void on_toolButton_help(int page);
    void resetTabRainfall();
    void resetTabOptions();
    void resetTabErosion();
    void resetTabFlow();
    void resetTabCalibration();
    void resetTabAdvanced();
    void resetTabInterception();
    void resetTabInfiltration();
    void doCheckRainfall(bool);
    void showTextfile(QString name);
    void showTextfileOld(QString name);


public slots:
    // functions linked to actions
    void saveRunFile();
    void savefileas();
    void openRunFile();
    void deleteRunFileList();
    void runmodel();
    void ClearOP();
    void stopmodel();
    void pausemodel();
    void shootScreen();
    void shootSingleScreen(int options);
    void shootMultipleScreens();
    void shootMScreen();
    void convertScreenshotsToVideo();
    void aboutQT();
    void aboutInfo();
    void resetAll();
    void setOutputScreen();
    void setOutputInfo(bool check);

    void onOutletChanged(int);
    void editMapname(QModelIndex topLeft, QModelIndex bottomRight );
    void openMapname(QModelIndex topLeft);

    QString findValidDir(QString path, bool up);
    //void on_toolButton_MapDir_clicked();
    void setMapDir();
    void setWorkDir();
    //void on_toolButton_ResultDir_clicked();
    void setResultDir();


    void doResetAll();
    void on_toolButton_resetOptions_clicked();
    void on_toolButton_resetRainfall_clicked();
    void on_toolButton_resetInterception_clicked();
    void on_toolButton_resetInfiltration_clicked();
    void on_toolButton_resetFlow_clicked();
    void on_toolButton_resetChannel_clicked();
    void on_toolButton_resetInfra_clicked();
    void on_toolButton_resetErosion_clicked();
    void on_toolButton_resetCalibration_clicked();
    void on_toolButton_resetAdvanced_clicked();

    void on_toolButton_helpOptions_clicked();
    void on_toolButton_helpRainfall_clicked();
    void on_toolButton_helpInterception_clicked();
    void on_toolButton_helpInfiltration_clicked();
    void on_toolButton_helpFlow_clicked();
    void on_toolButton_helpChannel_clicked();
    void on_toolButton_helpInfra_clicked();
    void on_toolButton_helpErosion_clicked();
    void on_toolButton_helpCalibration_clicked();
    void on_toolButton_helpAdvanced_clicked();

    void on_toolButton_RainfallName_clicked();
   // void on_toolButton_DichargeInName_clicked();

    //void on_toolButton_SnowmeltName_clicked();
    void on_toolButton_RainfallShow_clicked();
//    void on_toolButton_SnowmeltShow_clicked();
//    void on_toolButton_ShowRunfile_clicked();
    void on_toolButton_satImageName_clicked();
    //void on_toolButton_fileOpen_clicked();
    void on_toolButton_SwatreTableDir_clicked();
    void on_toolButton_SwatreTableFile_clicked();
    void on_toolButton_SwatreTableShow_clicked();
    void on_E_floodMinHeight_valueChanged(double);

    void on_DisplayComboBox_currentIndexChanged(int);
    void on_DisplayComboBox2_currentIndexChanged(int);

   // void doCheckSnowmelt(bool);

    void doCheckPesticides(bool check);

    void on_E_InfiltrationMethod_currentIndexChanged(int inr);
    void on_E_runFileList_currentIndexChanged(int);

    void on_checkFlowBarriers_clicked();
    void on_checkChannelInfil_clicked();
 //  void on_checkChannelBaseflow_clicked();
    void on_checkDoErosion_clicked();
    void on_checkIncludeChannel_clicked();
    void on_checkIncludeTiledrains_clicked();
    void on_checkBoxComboMaps_stateChanged(int);
    void on_checkBoxComboMaps2_stateChanged(int);
    void on_nrUserCores_valueChanged(int d);
    void on_ComboMinSpinBox_valueChanged(double);
    void on_ComboMaxSpinBox_valueChanged(double);
    void on_ComboMinSpinBox2_valueChanged(double);
    void on_ComboMaxSpinBox2_valueChanged(double);
    void onImageToggled(bool b);

    void setErosionMapOutput(bool doit);
    //void on_spinBoxPointtoShow_valueChanged(int);
    void on_tabWidget_out_currentChanged(int);

    //houses
    void on_checkHouses_clicked();
    void on_checkInfilCompact_clicked();
    void on_checkInfilCrust_clicked();
    void on_checkInfilGrass_clicked();
    //void on_checkInfil2layer_clicked();
    void on_checkSedtrap_clicked();
 //   void on_checkMaterialDepth_clicked();
//    void on_E_BulkDens2_editingFinished();
//    void on_E_BulkDens_editingFinished();
  //  void on_checkSnowmelt_clicked();
    void on_checkExpandActive_clicked();
    void on_E_MapDir_returnPressed();
    void on_E_ResultDir_returnPressed();

 //   void on_checkEstimateGrainSizeDistribution_toggled(bool v);
  //  void on_checkReadGrainSizeDistribution_toggled(bool v);

    void ssetAlpha(int v);
    void ssetAlphaChannel(int v);
    void ssetAlphaChannelOutlet(int v);
    void ssetAlphaRoad(int v);
    void ssetAlphaHouse(int v);
    void ssetAlphaHardSurface(int v);
    void ssetAlphaHardSurfaceW(int v);
    void ssetAlphaMap(int v);

    void setWriteOutputSOBEK(bool);
    void setWriteOutputCSV(bool);
    //void setWriteOutputPCR(bool);

    void setFloodTab(bool);
    void setErosionTab(bool);

    void setBWUI();
    void resizeMap();
    void fontSelect();
    void fontDecrease();
    void fontIncrease();
    void setfontSize();

    void setFormatMaps(bool);

    QString getFileorDir(QString inputdir,QString title, QStringList filters, int doFile);

private slots:

    void showMapb(bool);
    void showMapd(double);
    void showChannelVector(bool);

    // functions that interact with the world thread signals
    void worldShow(bool showall);
    void worldDone(const QString &results);
    void worldDebug(const QString &results);
    void worldTimedb(const QString &results);


    void on_check2DDiagonalFlow_toggled(bool checked);

    //void on_checkDiffusion_toggled(bool checked);

    void on_checkHouses_toggled(bool checked);

    void on_toolButton_rainsatName_clicked();

    void on_toolButton_RainmapShow_clicked();

    void on_toolButton_ETName_clicked();

    void on_toolButton_ETsatName_clicked();

   // void on_checkIncludeET_toggled(bool checked);

    void on_toolButton_ETShow_clicked();

   // void on_checkWaveInUser_toggled(bool checked);

    void on_toolButton_DischargeShow_clicked();

    void on_toolButton_DischargeName_clicked();

    void on_toolButton_WaveInName_clicked();

    void on_toolButton_WaveShow_clicked();

    void on_toolButton_ETmapShow_clicked();

    void on_E_EndTimeDay_returnPressed();

    void on_E_BeginTimeDay_returnPressed();

    void on_checkStationaryBaseflow_toggled(bool checked);

    void on_checkChannelInfil_toggled(bool checked);

    void on_E_EfficiencyDETCH_currentIndexChanged(int index);

    void on_checkGWflow_toggled(bool checked);


    void on_checkMB_WH_toggled(bool checked);

    void on_checkRainfall_toggled(bool checked);

    void on_checkET_toggled(bool checked);

    void on_checkInterception_toggled(bool checked);

    void on_E_OFWaveType_currentIndexChanged(int index);

    void on_checkIncludeChannel_toggled(bool checked);

    void on_checkDoErosion_toggled(bool checked);

    void on_checkConservation_toggled(bool checked);

    void on_checkInfiltration_toggled(bool checked);

    void on_toolButton_ShowRunfile_clicked();

    void on_checkInfrastructure_toggled(bool checked);

    void on_E_ETName_returnPressed();

    void on_E_RainsatName_returnPressed();

    void on_E_RainfallName_returnPressed();


private:

    QSystemTrayIcon *trayIcon;
    //toolbar actions
    QAction *openAct;
    QAction *saveAct;
    QAction *saveasAct;
    QAction *runAct;
    QAction *pauseAct;
    QAction *stopAct;
    QAction *shootscreenAct;
    QAction *shootMscreenAct;
    QAction *makeMovieAct;
    QAction *aboutAct;
    QAction *aboutActI;
    QAction *resetAllAct;
    QAction *showAllAct;
    QAction *showInfoAct;
    QAction *resizeAct;
    QAction *setBWAct;

    QAction *fontAct;
    QAction *fontIncreaseAct;
    QAction *fontDecreaseAct;
    // the model world
    TWorld *W;

};


#endif // LISEMQT_H
