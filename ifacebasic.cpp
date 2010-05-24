/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

#include "ifacebasic.h"
#include "model.h"
#include "global.h"

//---------------------------------------------------------------------------
output op;
// in and output structure to link model output to interface

//---------------------------------------------------------------------------
ifacebasic::ifacebasic(QWidget *parent)
: QWidget(parent)
{
	setupUi(this);
	W = NULL;
	// set up interface
	GetStorePath();

	E_runfilename->setText(op.runfilename);

	SetStyleUI();
}
//---------------------------------------------------------------------------
ifacebasic::~ifacebasic()
{
	StorePath();
}
//---------------------------------------------------------------------------
void ifacebasic::SetStyleUI()
{
	label_dx->setStyleSheet("* { background-color: #ffffff }");
	label_area->setStyleSheet("* { background-color: #ffffff }");
	label_time->setStyleSheet("* { background-color: #ffffff }");
	label_endtime->setStyleSheet("* { background-color: #ffffff }");
	label_runtime->setStyleSheet("* { background-color: #ffffff }");
	label_endruntime->setStyleSheet("* { background-color: #ffffff }");
	label_raintot->setStyleSheet("* { background-color: #ffff77 }");
	label_watervoltot->setStyleSheet("* { background-color: #ffff77 }");
	label_qtot->setStyleSheet("* { background-color: #ffff77 }");
	label_infiltot->setStyleSheet("* { background-color: #ffff77 }");
	label_interctot->setStyleSheet("* { background-color: #ffff77 }");
	label_qtotm3->setStyleSheet("* { background-color: #ffff77 }");
	label_qpeak->setStyleSheet("* { background-color: #ffff77 }");
	label_qpeaktime->setStyleSheet("* { background-color: #ffff77 }");

	label_splashdet->setStyleSheet("* { background-color: #ffff77 }");
	label_flowdet->setStyleSheet("* { background-color: #ffff77 }");
	label_sedvol->setStyleSheet("* { background-color: #ffff77 }");
	label_dep->setStyleSheet("* { background-color: #ffff77 }");
	label_detch->setStyleSheet("* { background-color: #ffff77 }");
	label_depch->setStyleSheet("* { background-color: #ffff77 }");
	label_sedvolch->setStyleSheet("* { background-color: #ffff77 }");
	label_soilloss->setStyleSheet("* { background-color: #ffff77 }");
	label_soillosskgha->setStyleSheet("* { background-color: #ffff77 }");
}
//---------------------------------------------------------------------------
/*
void ifacebasic::on_checkChannel_clicked()
{
	op.SwitchIncludeChannel = checkChannel->isChecked();
}
//---------------------------------------------------------------------------
void ifacebasic::on_checkErosion_clicked()
{
	op.SwitchErosion = checkErosion->isChecked();
}*/
//---------------------------------------------------------------------------
void ifacebasic::on_runButton_clicked()
{
	if(W)
		W->stopRequested = true;
	else
	{
		QFile file(op.runfilename);
		if (!file.open(QFile::ReadOnly | QFile::Text)) {
			QMessageBox::warning(this, "openLISEM",
					QString("Cannot read file \"%1\":\n%2.")
					.arg(op.runfilename)
					.arg(file.errorString()));
			return;
		}

		W = new TWorld();
		label_28->setText("debug");

		connect(W, SIGNAL(show(void)),this, SLOT(Showit(void)),Qt::QueuedConnection);
		connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
		connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
		// connect emitted signals from the model to the interface routines that handle them
		W->stopRequested = false;
		// stoprequested is used to stop the thread with the interface
		W->start();
		// start the model thread, executes W->run()
		runButton->setText("stop");
	}
}
//---------------------------------------------------------------------------
void ifacebasic::Showit()
{
	// copy the run results from the output structure op to the ui labels
	// op is filled in the model run each timestep
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
	label_interctot->setText(QString::number(op.IntercTotmm,'f',3));
	label_qtotm3->setText(QString::number(op.Qtot,'f',3));
	label_qpeak->setText(QString::number(op.Qpeak,'f',3));
	label_qpeaktime->setText(QString::number(op.QpeakTime,'f',3));

	if (op.SwitchErosion)
	{
		label_MBs->setText(QString::number(op.MBs,'e',3));
		label_splashdet->setText(QString::number(op.DetTotSplash,'f',3));
		label_flowdet->setText(QString::number(op.DetTotFlow,'f',3));
		label_sedvol->setText(QString::number(op.SedVolTot,'f',3));
		label_dep->setText(QString::number(op.DepTot,'f',3));

		label_detch->setText(QString::number(op.ChannelDetTot,'f',3));
		label_depch->setText(QString::number(op.ChannelDepTot,'f',3));
		label_sedvolch->setText(QString::number(op.ChannelSedTot,'f',3));

		label_soilloss->setText(QString::number(op.SoilLossTot,'f',3));
		label_soillosskgha->setText(QString::number(op.SoilLossTot/(op.CatchmentArea/10000)*1000,'f',3));
	}

	progressBar->setMaximum(op.maxstep);
	progressBar->setValue(op.runstep);

}
//---------------------------------------------------------------------------
void ifacebasic::worldDone(const QString &results)
{
	label_28->setText(results);
	// arrive here after model emits done signal
	if (W)
	{
		delete W;
		W=NULL;
	}
	runButton->setText("run");
	//free the model instance
}
//---------------------------------------------------------------------------
void ifacebasic::worldDebug(const QString &results)
{

	label_28->setText(results);
	// arrive here after model emits done signal
}
//---------------------------------------------------------------------------
void ifacebasic::on_toolButton_runfilename_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,"Select runfile","*.run");
	E_runfilename->setText( path );
	op.runfilename = path;
}

//---------------------------------------------------------------------------
void ifacebasic::StorePath()
{
	QFile fff(op.LisemDir + "openlisem.ini");
	if (!fff.open(QIODevice::WriteOnly | QIODevice::Text))
		return;

	QTextStream ts( &fff );

	ts << op.runfilename;// << endl;

	fff.close();
}
//---------------------------------------------------------------------------
void ifacebasic::GetStorePath()
{
	QFile fff(op.LisemDir + "openlisem.ini");
	if (!fff.open(QIODevice::ReadOnly | QIODevice::Text))
		return;

	op.runfilename = fff.readLine();

	fff.close();
}
//---------------------------------------------------------------------------
void ifacebasic::on_toolButton_ShowRunfile_clicked()
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
	view->setWindowTitle(op.runfilename);
	view->setMinimumWidth(400);
	view->setMinimumHeight(500);
	view->setAttribute(Qt::WA_DeleteOnClose);
	view->show();

	file.close();
}
