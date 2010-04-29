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

      //E_runfilename->setText("D:\\data\\yen bai\\data\\run\\storm_15.run");
      E_runfilename->setText(op.runfilename);
//      label_28->setText(op.LisemDir);
}
//---------------------------------------------------------------------------
ifacebasic::~ifacebasic()
{
  StorePath();
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
    W = new TWorld();
    label_28->setText("debug");

     connect(W, SIGNAL(show(int)),this, SLOT(Showit(int)),Qt::QueuedConnection);
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
void ifacebasic::Showit(const int step)
{
// copy the run results from the output structure op to the ui labels
// op is filled in the model run each timestep

	label->setNum(op.MB);
	label_2->setNum(op.RainTot);
	label_3->setNum(op.WaterVolTot);
	label_4->setNum(op.Qtotmm);
	label_5->setNum(op.InfilTot);
	label_6->setNum(op.IntercTot);
	//label_7->setNum(op.InfilKWTot);
	//label_8->setNum(step);
    label_34->setNum(op.Qtot);
    label_36->setNum(op.Qpeak);

	if (op.SwitchErosion)
	{
          label_13->setNum(op.MBs);
          label_18->setNum(op.DetTotSplash);
          label_20->setNum(op.DetTotFlow);
          label_22->setNum(op.SoilLossTot);
          label_24->setNum(op.SedVolTot);
          label_26->setNum(op.DepTot);
	}
	label_9->setNum(op.t);
	label_8->setNum(op.time);
        label_32->setNum(op.maxtime);
	progressBar->setMaximum(op.maxstep);
	progressBar->setValue(step);
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
                QMessageBox::warning(this, tr("Application"),
                                                         tr("Cannot read file %1:\n%2.")
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
