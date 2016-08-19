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
 \file lisemqt.h
 \brief main interface header file containing UI class, linked to lisemqt.ui
 */

#ifndef LISEMQT_H
#define LISEMQT_H

#include <QtGui>
#include <QApplication>

//QWT library files
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
#include <qwt_plot_grid.h>
#include <qwt_plot_rescaler.h>
#include <qwt_scale_engine.h>
#include <qwt_plot_zoomer.h>
#include <qwt_picker.h>


#include "version.h"
#include "ui_lisemqt.h"
#include "model.h"
#include "LisUIoutput.h"
#include "LisUItreemodel.h"
#include "LisUImapplot.h"


// constants to define the place of the main parts in the map tree structure
#define RAINFALLMAPS 0
#define CATCHMENTMAPS 1
#define LANDUSEMAPS 2
#define SURFACEMAPS 3
#define EROSIONMAPS 4
#define INFILTRATIONMAPS 5
#define CHANNELMAPS 6
#define CHANNELFLOODMAPS 7
#define BUFFERSMAPS 8
#define SNOWMELTMAPS 9
#define TILEDRAINMAPS 10   //VJ 110111
#define HOUSESMAPS 11  //VJ 120314
#define NUTRIENTSMAPS 12
//#define TEXTUREMAPS 13
//#define WHEELTRACKSMAPS 14
//#define GULLIESMAPS 15


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

    QProgressBar *pb;
    void resizeEvent(QResizeEvent* event);

    bool doBatchmode;
    QString batchRunname;

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
    void GetRunfile();
    void ParseInputData();
    void updateModelData();
    void defaultRunFile();
    QString CheckDir(QString p, bool makeit = false);
    void RunAllChecks();
    void savefile(QString name);
    void SetConnections();
    QStringList runfilelist;

    // structure output
    void initOP();
    // graph functions
    void setupPlot();
  //  void setupSmallPlot();
    void startPlots();
    void SetPlotsData();
    void showPlot();
    void GetPlotData();
    void initPlot();
    void killPlot();
 //   void showSmallPlot();
    void initOutputData();
    void showOutputData();
    void SetTextHydrographs();
    // Map drawing variable
    void setupMapPlot();
    void initMapPlot();

    void showMap();

    void showComboMap(int i);

    void showBaseMap();
    void showChannelMap();
    void showRoadMap();
    void showHouseMap();
    double fillDrawMapData(cTMap *_M, QwtMatrixRasterData *_RD, double type);

    QwtText title;
    QwtPlotSpectrogram *drawMap;  // raster map drawing
    QwtPlotSpectrogram *baseMap;  // raster map drawing
    QwtPlotSpectrogram *baseMapDEM;  // raster map drawing
    QwtPlotSpectrogram *channelMap;  // raster map drawing
    QwtPlotSpectrogram *roadMap;  // raster map drawing
    QwtPlotSpectrogram *houseMap;  // raster map drawing
    QwtPlot *MPlot;               // plot in which the raster map is drawn
    QwtMatrixRasterData *RD;
    QwtMatrixRasterData *RDb;
    QwtMatrixRasterData *RDbb;
    QwtMatrixRasterData *RDc;
    QwtMatrixRasterData *RDd;
    QwtMatrixRasterData *RDe;
    //   double drawNrCols;
    //   double drawNrRows;
    // vars for store map display in runfile
    int MapDisplayMapSelection, MapDisplayBuilding, MapDisplayRoads, MapDisplayChannels, MapDisplayHydrology;
    double MapDisplayRunoffMax, MapDisplayInfiltrationMax, MapDisplaySoillossMax, MapDisplayFlooddepthMax;
    int MapDisplayIncludeRunoff, MapDisplayMinimumDepth, MapDisplayScreenDumps;

    QString E_WorkDir;

    QVector<double> mapData;
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
    QList<QList<double>*> OutletQ;
    QList<QList<double>*> OutletQs;
    QList<QList<double>*> OutletC;
    QList<QList<double>*> OutletChannelWH;
    QList<double> OutletQpeak;
    QList<double> OutletQpeaktime;
    QList<double> OutletQtot;
    QList<double> OutletQstot;
    QList<double> Rainfall;
    double timestep;

    int outletpoint = 1;

    QList<double> qmax;
    QList<double> qsmax;
    QList<double> cmax;

    //Map display options
    QList<QwtComboColorMap *> ColorMapList;
    QList<QString> NameList;
    QList<QString> UnitList;
    QList<bool> SymList;
    QList<bool> LogList;
    //QList<int> ListList;
    QList<int> IndexList;
    QList<int> IndexList1;
    int ActiveList = 0;

    // graph variables
    QwtPlot *HPlot;
   // QwtPlot *smallPlot;
    QwtPlotCurve *QGraph;
    QwtPlotCurve *QsGraph;
    QwtPlotCurve *CGraph;
    QwtPlotCurve *PGraph;
    QwtPlotCurve *QtileGraph;
//    QwtPlotCurve *sPGraph;
//    QwtPlotCurve *sQGraph;
//    QwtPlotCurve *sQsGraph;
    bool startplot = true;
    bool stopplot;
    double yas, yasP, y2as;
    QVector <double> QData;
    QVector <double> QData1;
    QVector <double> QData2;
    QVector <double> QtileData;
    QVector <double> QsData;
    QVector <double> CData;
    QVector <double> PData;
    QVector <double> TData;
    long stepP;

    bool oldRunfile; // check is old runfile for ksat calibration

    //interface names
    TreeModel *MapNameModel;
    QString currentDir;
    QString RainFileName;
    QString RainFileDir;
    QString SnowmeltFileName;
    QString SnowmeltFileDir;
    QString SwatreTableName;
    QString SwatreTableDir;
    QStringList DEFmaps;
    QStringList RunFileNames;
    int CurrentRunFile;
    int uiInfilMethod;
    double swatreDT;

    QList<cTMap *> ComboMapsSafe;

    MAP_LIST mapList[NUMMAPS]; /// structure for current map names, can be edited by user
    int nrmaplist;
    NAME_LIST namelist[NUMNAMES]; /// structure to read all runfile variables and names
    int nrnamelist;
    QStringList outputcheck; /// list of '0' and '1' to see which output mapseries are checled by the nuser
    int InterceptionEqNr;
    int mapstartnr;
    bool doShootScreens;


public slots:
    // functions linked to actions
    void saveRunFile();
    void savefileas();
    void openRunFile();
    void deleteRunFileList();
    void runmodel();
    void stopmodel();
    void pausemodel();
    void shootScreen();
    void shootMScreen();
    void aboutQT();
    void aboutInfo();
    void resetAll();

    void onOutletChanged(int);
    void editMapname(QModelIndex topLeft, QModelIndex bottomRight );
    void openMapname(QModelIndex topLeft);

    QString findValidDir(QString path, bool up);
    //void on_toolButton_MapDir_clicked();
    void setMapDir();
    void setWorkDir();
    //void on_toolButton_ResultDir_clicked();
    void setResultDir();
    void on_toolButton_RainfallName_clicked();
    void on_toolButton_SnowmeltName_clicked();
    void on_toolButton_RainfallShow_clicked();
    void on_toolButton_SnowmeltShow_clicked();
    void on_toolButton_ShowRunfile_clicked();
    //void on_toolButton_fileOpen_clicked();
    void on_toolButton_SwatreTableDir_clicked();
    void on_toolButton_SwatreTableFile_clicked();
    void on_toolButton_SwatreTableShow_clicked();

    void on_checkBox_SedSingleSingle_toggled(bool v);
    void on_checkBox_SedMultiSingle_toggled(bool v);
    void on_checkBox_SedMultiMulti_toggled(bool v);

    void on_DisplayComboBox_currentIndexChanged(int);
    void on_DisplayComboBox2_currentIndexChanged(int);

    void doCheckSnowmelt(bool);
    void doCheckRainfall(bool);
    void doCheckPesticides(bool check);

    void on_E_InfiltrationMethod_currentIndexChanged(int inr);
    void on_E_runFileList_currentIndexChanged(int);

    void on_checkChannelInfil_clicked();
    void on_checkChannelBaseflow_clicked();
    void on_checkChannelFlood_clicked();
    void on_checkDoErosion_clicked();
    void on_checkOverlandFlow1D_clicked();
    void on_checkOverlandFlow2D_clicked();
    void on_checkIncludeChannel_clicked();
    void on_checkIncludeTiledrains_clicked();
    void on_checkImpermeable_stateChanged(int);
    void on_checkPercolation_stateChanged(int);
    void on_ComboMinSpinBox_valueChanged(double);
    void on_ComboMaxSpinBox_valueChanged(double);
    void on_ComboMinSpinBox2_valueChanged(double);
    void on_ComboMaxSpinBox2_valueChanged(double);
    void setDisplayComboBox(int);
    void setDisplayComboBox2(int);
    void on_spinBoxPointtoShow_valueChanged(int);
    void on_tabWidget_out_currentChanged(int);

    //houses
    void on_checkHouses_clicked();
    void on_checkInfilCompact_clicked();
    void on_checkInfilCrust_clicked();
    void on_checkInfilGrass_clicked();
    void on_checkInfil2layer_clicked();
    void on_checkBuffers_clicked();
    void on_checkSedtrap_clicked();
    void on_checkSnowmelt_clicked();
    void on_checkExpandActive_clicked();
    void on_E_MapDir_returnPressed();
    void on_E_ResultDir_returnPressed();

    void on_checkEstimateGrainSizeDistribution_toggled(bool v);
    void on_checkReadGrainSizeDistribution_toggled(bool v);

    void ssetAlpha(int v);
    void ssetAlpha2(int v);
    void ssetAlpha3(int v);
    void ssetAlpha4(int v);

    void setWriteOutputSOBEK(bool);
    void setWriteOutputCSV(bool);
    void setWriteOutputPCR(bool);

    void setErosionTab(bool);

    void fontSelect();
    void fontDecrease();
    void fontIncrease();


private slots:

    void showMapb(bool);
    void showMapd(double);


    // functions that interact with the world thread signals
    void worldShow();
    void worldDone(const QString &results);
    void worldDebug(const QString &results);


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
    QAction *aboutAct;
    QAction *aboutActI;
    QAction *restartAct;

    QAction *fontAct;
    QAction *fontIncreaseAct;
    QAction *fontDecreaseAct;
    // the model world
    TWorld *W;

};


#endif // LISEMQT_H
