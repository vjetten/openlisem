#ifndef LISEMQT_H
#define LISEMQT_H

#include <QtGui>
#include <QApplication>

#include "ui_lisemqt.h"
#include "LisUItreemodel.h"
#include "model.h"

struct output{
	int runstep;
	int printstep;
	int maxstep;
	double CatchmentArea, dx, t,
	time, maxtime, EndTime;

	double MB, Qtot, Qtotmm, Qpeak, IntercTotmm, WaterVolTotmm, InfilTotmm,
	RainTotmm, SurfStormm, InfilKWTotmm,
	MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedTot,
	ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot,
	RunoffFraction, RainpeakTime, QpeakTime;

	bool SwitchErosion;
	bool SwitchIncludeChannel;
	QString runfilename;
	QString LisemDir;
};


class lisemqt : public QMainWindow, private Ui::lisemqtClass
{
	Q_OBJECT

public:
	lisemqt(QWidget *parent = 0);
	~lisemqt();

	void FillMapList();
	void LisemReset();
	void DefaultMapnames();
	void change_MapNameModel(int parentrow, int selrow, bool setit);
	void SetMenuandToolBar();

	int prevsel, prevselinf;
	QString RainFileName;
	QString SnowmeltFileName;
	TreeModel *MapNameModel;
	QStringList DEFmaps;
	QStringList OUTmaps;
	QStringList OUTMCmaps;
	QStringList OUTNUTmaps;
	QStringList OUTGULmaps;

public slots:

void GetStorePath();
void StorePath();

void on_toolButton_MapDir_clicked();
void on_toolButton_ResultDir_clicked();
void on_toolButton_RainfallName_clicked();
void on_toolButton_SnowmeltName_clicked();
void on_toolButton_RainfallShow_clicked();
void on_toolButton_SnowmeltShow_clicked();
void on_toolButton_ShowRunfile_clicked();

void on_E_LisemType_currentIndexChanged(int);
void on_E_InfiltrationMethod_currentIndexChanged(int);

void on_checkChannelInfil_clicked();
void on_checkChannelBaseflow_clicked();
void on_checkNoErosion_clicked();
void on_checkIncludeChannel_clicked();
void on_checkInfilCompact_clicked();
void on_checkInfilCrust_clicked();
void on_checkInfilGrass_clicked();
void on_checkBuffers_clicked();
void on_checkSedtrap_clicked();
void on_checkSnowmelt_clicked();
void on_checkExpandActive_clicked();
void on_toolButton_SwatreTable_clicked();

void savefile();
void openRunFile();
void runmodel();

private:
QAction *openAct;
QAction *saveAct;
QAction *saveasAct;
QAction *runAct;
QAction *pauseAct;
QAction *stopAct;

TWorld *W;

};


#endif // LISEMQT_H
