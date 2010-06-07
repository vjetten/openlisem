/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * This file maintains the link between the interface and the LISEM model
 * with data and graph displayed on the screen
 * it creates the model "world" W and activates and kills the model thread
 */


#include "lisemqt.h"
#include "model.h"
#include "global.h"


//---------------------------------------------------------------------------
void lisemqt::stopmodel()
{
	if(W)
		W->stopRequested = true;
}
//---------------------------------------------------------------------------
void lisemqt::runmodel()
{
	tabWidget->setCurrentIndex(2);

   startplot = true;
   QData = NULL;
	QsData = NULL;
	CData = NULL;
	PData = NULL;
	timeData = NULL;
	//intialize plot stuff for this run

	//TODO replace with temp run file
	QFile file(op.runfilename);
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this, "openLISEM",
				QString("Cannot read file \"%1\":\n%2.")
				.arg(op.runfilename)
				.arg(file.errorString()));
		return;
	}
	// check run file

	W = new TWorld();
	// make the world

	connect(W, SIGNAL(show(void)),this, SLOT(Showit(void)),Qt::QueuedConnection);
	connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
	connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
	// connect emitted signals from the model thread to the interface routines that handle them

	W->stopRequested = false;
	// stoprequested is used to stop the thread with the interface
	W->start();
	// start the model thread, executes W->run()
}
//---------------------------------------------------------------------------
void lisemqt::Showit()
{
	// copy the run results from the "output structure op" to the ui labels
	// "op" is filled in the model run each timestep
	// "op" struct is declared in lisUIoutput.h

	label_dx->setText(QString::number(op.dx,'f',3));
	label_area->setText(QString::number(op.CatchmentArea/10000,'f',3));
	label_time->setText(QString::number(op.time,'f',3));
	label_endtime->setText(QString::number(op.EndTime,'f',3));
	label_runtime->setText(QString::number(op.t,'f',3));
	label_endruntime->setText(QString::number(op.maxtime,'f',3));

	label_MB->setText(QString::number(op.MB,'e',3));
	label_raintot->setText(QString::number(op.RainTotmm,'f',3));
	label_watervoltot->setText(QString::number(op.WaterVolTotmm,'f',3));
	label_qtot->setText(QString::number(op.Qtotmm,'f',3));
	label_infiltot->setText(QString::number(op.InfilTotmm,'f',3));
	label_surfstor->setText(QString::number(op.SurfStormm,'f',3));
	label_interctot->setText(QString::number(op.IntercTotmm,'f',3));
	label_qtotm3->setText(QString::number(op.Qtot,'f',3));
	label_qpeak->setText(QString::number(op.Qpeak,'f',3));
	label_qpeaktime->setText(QString::number(op.QpeakTime,'f',3));
	label_ppeaktime->setText(QString::number(op.RainpeakTime,'f',3));
	label_discharge->setText(QString::number(op.Q,'f',3));
	label_QPfrac->setText(QString::number((op.RainTotmm > 0 ? op.Qtotmm/op.RainTotmm*100 : 0),'f',3));
	if (op.SwitchErosion)
	{
		label_MBs->setText(QString::number(op.MBs,'e',3));
		label_splashdet->setText(QString::number(op.DetTotSplash,'f',3));
		label_flowdet->setText(QString::number(op.DetTotFlow,'f',3));
		label_sedvol->setText(QString::number(op.SedTot,'f',3));
		label_dep->setText(QString::number(op.DepTot,'f',3));

		label_detch->setText(QString::number(op.ChannelDetTot,'f',3));
		label_depch->setText(QString::number(op.ChannelDepTot,'f',3));
		label_sedvolch->setText(QString::number(op.ChannelSedTot,'f',3));

		label_soilloss->setText(QString::number(op.SoilLossTot,'f',3));
		label_soillosskgha->setText(QString::number(op.SoilLossTot/(op.CatchmentArea/10000)*1000,'f',3));
	}

	if (checkBuffers->isChecked())
	{
		label_buffervol->setText(QString::number(op.BufferVolTot,'f',3));
		label_buffersed->setText(QString::number(op.BufferSedTot,'f',3));
	}

	progressBar->setMaximum(op.maxstep);
	progressBar->setValue(op.runstep);

	ShowGraph();

	QString SSS;
	SSS = QString("%1 %2 %3 %4 %5").arg(op.time,15,'f',3,' ').arg(op.P,15,'f',3,' ').arg(op.Q,15,'f',3,' ').arg(op.Qs,12,'f',3).arg(op.C,15,'f',3,' ');
	textGraph->setMaximumBlockCount(6);
	textGraph->appendPlainText(SSS);

}
//---------------------------------------------------------------------------
void lisemqt::ShowGraph()
{
	if (startplot)
	{

		startplot = false;

		yas = 0.1;
		y2as = 0.1;

		timeData = new double[op.maxstep+2];
		QData = new double[op.maxstep+2];
		QsData = new double[op.maxstep+2];
		CData = new double[op.maxstep+2];
		PData = new double[op.maxstep+2];
		for (int i = 0; i < op.maxstep; i++)
		{
			timeData[i] = 0;
			QData[i] = 0;
			QsData[i] = 0;
			CData[i] = 0;
			PData[i] = 0;
		}

		HPlot->setAxisScale(HPlot->xBottom, op.BeginTime, op.EndTime);
	}
	QGraph->setRawData(timeData,QData,op.runstep);
	QsGraph->setRawData(timeData,QsData,op.runstep);
	CGraph->setRawData(timeData,CData,op.runstep);
	PGraph->setRawData(timeData,PData,op.runstep);

	timeData[op.runstep] = op.time;
	QData[op.runstep] = op.Q;
	QsData[op.runstep] = op.Qs;
	CData[op.runstep] = op.C;
	PData[op.runstep] = op.P;


	y2as = max(y2as, op.Qs);
	y2as = max(y2as, op.C);
	HPlot->setAxisScale(HPlot->yRight, 0, y2as*1.05);
	yas = max(yas, op.Q);
	yas = max(yas, op.P);
	HPlot->setAxisScale(HPlot->yLeft, 0, yas*1.05);

	HPlot->replot();

}
//---------------------------------------------------------------------------
void lisemqt::worldDone(const QString &results)
{
	label_28->setText(results);
	// arrive here after model emits done signal
	if (W)
	{
		delete W;
		W=NULL;
	}
	//free the world instance

	delete QData;
	delete QsData;
	delete CData;
	delete PData;
	delete timeData;
   QData = NULL;
	QsData = NULL;
	CData = NULL;
	PData = NULL;
	timeData = NULL;
	// free data structures graph

}
//---------------------------------------------------------------------------
void lisemqt::worldDebug(const QString &results)
{

	label_debug->setText(results);
	// show messages from the World model on the screen
}
