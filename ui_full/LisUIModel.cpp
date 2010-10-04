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
void lisemqt::pausemodel()
{
	if(W)
	{
		W->waitRequested = !W->waitRequested;
		if (!W->waitRequested)
		{
			label_debug->setText("User continue...");
			W->condition.wakeAll();
		}
	}
}
//---------------------------------------------------------------------------
void lisemqt::stopmodel()
{
	if(W)
		W->stopRequested = true;
}
//---------------------------------------------------------------------------
void lisemqt::runmodel()
{
	if (op.runfilename.isEmpty())
	{
		QMessageBox::warning(this,"openLISEM",
				QString("Load a runfile first!"));
		return;
	}
	label_runfilename->setText(op.runfilename);
	//TODO if run from commandline this name must exist
	// this assumes runfilename is correct
	savefile(QString(op.LisemDir+"openlisemtmp.run"));

	tabWidget->setCurrentIndex(2);
	//switch to output screen

	startplot = true;
	QData = NULL;
	QsData = NULL;
	CData = NULL;
	PData = NULL;
	timeData = NULL;
	//intialize plot stuff for this run
	InitOP();
	// wipe the result screen

	W = new TWorld();
	// make the world

	connect(W, SIGNAL(show(void)),this, SLOT(Showit(void)),Qt::QueuedConnection);
	connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
	connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
	// connect emitted signals from the model thread to the interface routines that handle them

	W->stopRequested = false;
	// stoprequested is used to stop the thread with the interface
	W->waitRequested = false;
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
	if (checkBuffers->isChecked())
		label_buffervol->setText(QString::number(op.BufferVolTot,'f',3));

	if (!checkNoErosion->isChecked())
	{
		int dig = 2;
		label_MBs->setText(QString::number(op.MBs,'e',dig));
		label_splashdet->setText(QString::number(op.DetTotSplash,'f',dig));
		label_flowdet->setText(QString::number(op.DetTotFlow,'f',dig));
		label_sedvol->setText(QString::number(op.SedTot,'f',dig));
		label_dep->setText(QString::number(op.DepTot,'f',dig));

		label_detch->setText(QString::number(op.ChannelDetTot,'f',dig));
		label_depch->setText(QString::number(op.ChannelDepTot,'f',dig));
		label_sedvolch->setText(QString::number(op.ChannelSedTot,'f',dig));

		label_soilloss->setText(QString::number(op.SoilLossTot,'f',dig));
		label_soillosskgha->setText(QString::number(op.SoilLossTot/(op.CatchmentArea/10000)*1000,'f',dig));

		double SDR = op.DetTotSplash + op.ChannelDetTot + op.DetTotFlow;
		SDR = (SDR > 0? 100*op.SoilLossTot/(SDR) : 0);
		SDR = min(SDR ,100);
		label_SDR->setText(QString::number(SDR,'f',dig));
		if (checkBuffers->isChecked() || checkSedtrap->isChecked())
			label_buffersed->setText(QString::number(op.BufferSedTot,'f',dig));
	}

	progressBar->setMaximum(op.maxstep);
	progressBar->setValue(op.runstep);

	ShowGraph();

	ShowMap();

	QString SSS;
	SSS = QString("%1 %2 %3 %4 %5").arg(op.time,15,'f',3,' ').arg(op.P,15,'f',3,' ').arg(op.Q,15,'f',3,' ').arg(op.Qs,12,'f',3).arg(op.C,15,'f',3,' ');
	textGraph->setMaximumBlockCount(6);
	textGraph->appendPlainText(SSS);

}
//---------------------------------------------------------------------------
void lisemqt::ShowMap()
{
	//MapPlot->replot();
	//if all calcs are done in the model only replot here?
	//op needs a poiinter to TMMap
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
		HPlot->setAxisScale(HPlot->xBottom, op.BeginTime, op.EndTime);
		qDebug() << op.BeginTime;
	}

	timeData[op.runstep] = op.time;
	PData[op.runstep] = op.P;
	QData[op.runstep] = op.Q;
	QsData[op.runstep] = op.Qs;
	CData[op.runstep] = op.C;
	// to avoid strange graphs
	timeData[op.runstep+1] = op.time;
	PData[op.runstep+1] = op.P;
	QData[op.runstep+1] = op.Q;
	QsData[op.runstep+1] = op.Qs;
	CData[op.runstep+1] = op.C;

	QGraph->setRawData(timeData,QData,op.runstep);
	PGraph->setRawData(timeData,PData,op.runstep);
	if(!checkNoErosion->isChecked())
	{
		QsGraph->setRawData(timeData,QsData,op.runstep);
		CGraph->setRawData(timeData,CData,op.runstep);
	}

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
	label_debug->setText(results);
	if (results.contains("ERROR"))
		QMessageBox::critical(this,QString("openLISEM"), results, QMessageBox::Ok );

	shootScreen();

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

	QFile(QString(op.LisemDir+"openlisemtmp.run")).remove();
}
//---------------------------------------------------------------------------
void lisemqt::worldDebug(const QString &results)
{

	label_debug->setText(results);
	// show messages from the World model on the screen
}
//---------------------------------------------------------------------------
void lisemqt::InitOP()
{
	op.runstep = 0;
	op.printstep = 0;
	op.maxstep = 0;
	op.CatchmentArea = 0;
	op.dx = 0;
	op.t = 0;
	op.time = 0;
	op.maxtime = 0;
	op.EndTime = 0;
	op.BeginTime = 0;
	op.MB = 0;
	op.Qtot = 0;
	op.Qtotmm = 0;
	op.Qpeak = 0;
	op.IntercTotmm = 0;
	op.WaterVolTotmm = 0;
	op.InfilTotmm = 0;
	op.RainTotmm = 0;
	op.SurfStormm = 0;
	op.InfilKWTotmm = 0;
	op.MBs = 0;
	op.DetTot = 0;
	op.DetTotSplash = 0;
	op.DetTotFlow = 0;
	op.DepTot = 0;
	op.SoilLossTot = 0;
	op.SedTot = 0;
	op.ChannelVolTot = 0;
	op.ChannelSedTot = 0;
	op.ChannelDepTot = 0;
	op.ChannelDetTot = 0;
	op.RunoffFraction = 0;
	op.RainpeakTime = 0;
	op.QpeakTime = 0;
	op.Q = 0;
	op.Qs = 0;
	op.C = 0;
	op.P = 0;
	op.BufferVolTot = 0;
	op.BufferSedTot = 0;
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
	if (checkBuffers->isChecked())
		label_buffervol->setText(QString::number(op.BufferVolTot,'f',3));

	if (!checkNoErosion->isChecked())
	{
		int dig = 2;
		label_MBs->setText(QString::number(op.MBs,'e',dig));
		label_splashdet->setText(QString::number(op.DetTotSplash,'f',dig));
		label_flowdet->setText(QString::number(op.DetTotFlow,'f',dig));
		label_sedvol->setText(QString::number(op.SedTot,'f',dig));
		label_dep->setText(QString::number(op.DepTot,'f',dig));

		label_detch->setText(QString::number(op.ChannelDetTot,'f',dig));
		label_depch->setText(QString::number(op.ChannelDepTot,'f',dig));
		label_sedvolch->setText(QString::number(op.ChannelSedTot,'f',dig));

		label_soilloss->setText(QString::number(op.SoilLossTot,'f',dig));
		label_soillosskgha->setText(QString::number(0,'f',dig));
		label_SDR->setText(QString::number(0,'f',dig));
		if (checkBuffers->isChecked() || checkSedtrap->isChecked())
			label_buffersed->setText(QString::number(op.BufferSedTot,'f',dig));
	}

	textGraph->clear();

}
//---------------------------------------------------------------------------
