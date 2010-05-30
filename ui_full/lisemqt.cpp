/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * This is the extensive interface
 */

#include "lisemqt.h"
#include "model.h"
#include "global.h"


//--------------------------------------------------------------------
lisemqt::lisemqt(QWidget *parent)
: QMainWindow(parent)
{
	setupUi(this);
	// set up interface
   resize(920, 700);

	MapNameModel = NULL;
	LisemReset();
	DefaultMapnames();
	FillMapList();
	// initalize interface

	W = NULL;
	// initalize pointer to the world, created when run button is pushed

	GetStorePath();
//	E_runfilename->setText(op.runfilename);
	// open openlisem.ini and read last runfile name
	E_runFileList->addItem(op.runfilename);

//	SetStyleUI();
	// do some style things
}
//--------------------------------------------------------------------
lisemqt::~lisemqt()
{
	StorePath();
}
//--------------------------------------------------------------------
void lisemqt::LisemReset()
{
	/*
	E_LisemType->addItem("LISEM Basic");
	E_LisemType->addItem("LISEM Wheeltracks");
	E_LisemType->addItem("LISEM Multiclass Sediment");
	E_LisemType->addItem("LISEM Nutrients");
	E_LisemType->addItem("LISEM Gullies");
	 */
	E_InfiltrationMethod->addItem("no Infiltration");
	E_InfiltrationMethod->addItem("SWATRE");
	E_InfiltrationMethod->addItem("Green and Ampt");
	E_InfiltrationMethod->addItem("Smith and Parlange");
	E_InfiltrationMethod->addItem("Subtract Ksat");

	prevsel = 0;
	prevselinf = 0;
	RainFileName = "";

	SetMenuandToolBar();


}
//--------------------------------------------------------------------
void lisemqt::on_E_LisemType_currentIndexChanged(int)
{
	if (MapNameModel)
	{
		/*
		int selrow = E_LisemType->currentIndex()+7;
		if (E_LisemType->currentIndex() > 0)
			change_MapNameModel(selrow, 0, true);
		if (prevsel >0)
			change_MapNameModel(prevsel, 0, false);
		prevsel = selrow;
		 */
	}
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_MapDir_clicked()
{
	QString path;
	path = QFileDialog::getExistingDirectory(this,
			tr("Select maps directory"),
			QString::null,
			QFileDialog::ShowDirsOnly);

	E_MapDir->setText( path );
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ResultDir_clicked()
{
	QString path;
	path = QFileDialog::getExistingDirectory(this,
			tr("Select or create a directory to write results"),
			QString::null,
			QFileDialog::ShowDirsOnly);

	E_ResultDir->setText( path );
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallName_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
			tr("Select rainfall file"),
			QString::null);
	E_RainfallName->setText( path );
	RainFileName = path;
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltName_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
			tr("Select snow melt file"),
			QString::null);

	E_SnowmeltName->setText( path );
	SnowmeltFileName = path;
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltShow_clicked()
{
	QFile file(SnowmeltFileName);
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this,
			QString("Snowmelt fluxes"),
			QString("Cannot read file %1:\n%2.")
			.arg(SnowmeltFileName)
			.arg(file.errorString()));
		return;
	}

	QTextStream in(&file);

	QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
	view->setWindowTitle(SnowmeltFileName);
	view->setMinimumWidth(400);
	view->setMinimumHeight(500);
	view->setAttribute(Qt::WA_DeleteOnClose);
	view->show();

	file.close();
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallShow_clicked()
{

	QFile file(RainFileName);
	if (!file.open(QFile::ReadOnly | QFile::Text)) {
		QMessageBox::warning(this, tr("Application"),
				tr("Cannot read file %1:\n%2.")
				.arg(RainFileName)
				.arg(file.errorString()));
		return;
	}

	QTextStream in(&file);

	QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
	view->setWindowTitle(RainFileName);
	view->setMinimumWidth(400);
	view->setMinimumHeight(500);
	view->setAttribute(Qt::WA_DeleteOnClose);
	view->show();

	file.close();
}
//--------------------------------------------------------------------
void lisemqt::on_checkNoErosion_clicked()
{
	change_MapNameModel(3, 0, !checkNoErosion->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_E_InfiltrationMethod_currentIndexChanged(int)
{
	if (MapNameModel)
	{
		int selrow = E_InfiltrationMethod->currentIndex() + 9;
		if (E_InfiltrationMethod->currentIndex() == 0)
		{
			selrow = -1;
			treeView->collapse(MapNameModel->index(4,0));
		}

		if (prevselinf > 0)
			change_MapNameModel(4, prevselinf, false);
		prevselinf = selrow;

		change_MapNameModel(4,selrow, true);

	}
	groupBox_SwatreOptions->setEnabled(E_InfiltrationMethod->currentIndex() == 1);

}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SwatreTable_clicked()
{
	QString path;
	path = QFileDialog::getExistingDirectory(
			this,
			tr("Select directory with SWATRE tables"),
			QString::null,
			QFileDialog::ShowDirsOnly);

	E_SWATRETableDir->setText( path );
	//SWATRETableDir = path;
}   // E_SwatreDTSEC->Text = value;

//--------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_clicked()
{
	if (checkIncludeChannel->isChecked())
	{
		change_MapNameModel(5, 11, checkChannelInfil->isChecked());
		change_MapNameModel(5, 12, checkChannelBaseflow->isChecked());
	}
	change_MapNameModel(5, 10, checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelInfil_clicked()
{
	if (checkChannelBaseflow->isChecked())
		checkChannelBaseflow->setChecked(false);
	change_MapNameModel(5, 12, checkChannelBaseflow->isChecked());
	change_MapNameModel(5, 11, checkChannelInfil->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelBaseflow_clicked()
{
	if (checkChannelInfil->isChecked())
		checkChannelInfil->setChecked(false);
	change_MapNameModel(5, 11, checkChannelInfil->isChecked());
	change_MapNameModel(5, 12, checkChannelBaseflow->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilCompact_clicked()
{
	if (E_InfiltrationMethod->currentIndex() > 0)
		change_MapNameModel(4, 15, checkInfilCompact->isChecked());
	else
	{
		checkInfilCompact->setChecked(false);
		QMessageBox msgBox;
		msgBox.setText("Select an infiltration method first.");
		msgBox.exec();
	}
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilCrust_clicked()
{
	if (E_InfiltrationMethod->currentIndex()> 0)
		change_MapNameModel(4, 15, checkInfilCrust->isChecked());
	else
	{
		checkInfilCrust->setChecked(false);
		QMessageBox msgBox;
		msgBox.setText("Select an infiltration method first.");
		msgBox.exec();
	}

}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilGrass_clicked()
{
	if (E_InfiltrationMethod->currentIndex()> 0)
		change_MapNameModel(4, 15, checkInfilGrass->isChecked());
	else
	{
		checkInfilGrass->setChecked(false);
		QMessageBox msgBox;
		msgBox.setText("Select an infiltration method first.");
		msgBox.exec();
	}
}
//--------------------------------------------------------------------
void lisemqt::on_checkBuffers_clicked()
{
	change_MapNameModel(6, 0, checkBuffers->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkSedtrap_clicked()
{
	change_MapNameModel(6, 0, checkSedtrap->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkSnowmelt_clicked()
{
	change_MapNameModel(7, 0, checkSnowmelt->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkExpandActive_clicked()
{
	if (!checkExpandActive->isChecked())
		treeView->collapseAll();
	else
		for (int i = 0; i < MapNameModel->rowCount(); i++)
		{
			if (MapNameModel->getflag(i))
				treeView->expand(MapNameModel->index(i,0));
			/*
                QModelIndex indexParent = MapNameModel->index(i, 0);
                QModelIndex indexChild = MapNameModel->index(i, 0, indexParent);

                for (int k = 0; k < MapNameModel->rowCount(indexChild); k++)
                        if (MapNameModel->getflag(k, indexChild))
                                treeView->expand(MapNameModel->index(k,0));
			 */
		}

	//  treeView->resizeColumnToContents(0);
	//    treeView->resizeColumnToContents(1);
	//treeView->resizeColumnToContents(1);
}

//--------------------------------------------------------------------
void lisemqt::SetMenuandToolBar()
{
	//QSize is = QSize(21,21);
	//toolBar->setIconSize(is);

	openAct = new QAction(QIcon(":/fileopen.png"), "&Open...", this);
	openAct->setShortcuts(QKeySequence::Open);
	openAct->setStatusTip("Open a run file");
	connect(openAct, SIGNAL(triggered()), this, SLOT(openRunFile()));
	menu_File->addAction(openAct);
	toolBar->addAction(openAct);

	saveAct = new QAction(QIcon(":/filesave.png"), "&Save...", this);
	saveAct->setShortcuts(QKeySequence::Save);
	saveAct->setStatusTip("Save a run file");
	connect(saveAct, SIGNAL(triggered()), this, SLOT(savefile()));
	menu_File->addAction(saveAct);
	toolBar->addAction(saveAct);

	saveasAct = new QAction(QIcon(":/filesaveas.png"), "Save &As...", this);
	saveasAct->setShortcuts(QKeySequence::SaveAs);
	saveasAct->setStatusTip("Save a run file as ...");
	connect(saveasAct, SIGNAL(triggered()), this, SLOT(savefile()));
	menu_File->addAction(saveasAct);
	toolBar->addAction(saveasAct);

	runAct = new QAction(QIcon(":/start1.png"), "Run model...", this);
	//	runAct->setShortcuts(QKeySequence(Qt::CTRL + Qt::Key_R));
	runAct->setStatusTip("run the model ...");
	connect(runAct, SIGNAL(triggered()), this, SLOT(runmodel()));
	menu_File->addAction(runAct);
	toolBar->addAction(runAct);

}
//--------------------------------------------------------------------
void lisemqt::runmodel()
{
	StorePath();
}
//--------------------------------------------------------------------
void lisemqt::savefile()
{
	//label->setText("saving ...");
	QFile fout("hup.txt");
	fout.open(QIODevice::ReadWrite);
	for (int i = 0; i < MapNameModel->rowCount(); i++)
	{
		QVariant d;
		int _r = MapNameModel->index(i,0).row();
		QString S = "["+QString::number(_r)+"-" + QString::number(0) + "] ";
		for (int j = 0; j < 2; j++)
		{
			d = MapNameModel->data(MapNameModel->index(i,j),0); // parent
			S = S + d.toString()+";";
		}
		S = S + "\n";
		QModelIndex indexParent = MapNameModel->index(i, 0);
		for (int j = 0; j < MapNameModel->rowCount(indexParent); j++)
		{
			int _r = MapNameModel->index(i,j).row();
			int _c = MapNameModel->index(i,j).column();
			S = S + "["+QString::number(_r)+"-" + QString::number(_c) + "] ";

			for (int k = 0; k < MapNameModel->columnCount(indexParent); k++)
			{
				int _rr = MapNameModel->index(j, k, indexParent).row();
				int _cc = MapNameModel->index(j, k, indexParent).column();
				d = MapNameModel->data(MapNameModel->index(j, k, indexParent),0);
				S = S + "["+QString::number(_rr)+"-" + QString::number(_cc) + "] ";
				S = S + d.toString()+";";
			}
			S = S + "\n";
		}
		S = S + "\n";
		QByteArray line(S.toAscii());
		fout.write(line);
	}
	fout.close();

}
//--------------------------------------------------------------------
void lisemqt::openRunFile()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
			QString("Select run file(s)"),
			QString("c:/lisemcourse"),
			QString("*.run"));
	E_RainfallName->setText( path );
	RainFileName = path;
}
//---------------------------------------------------------------------------
void lisemqt::GetStorePath()
{
	QFile fff(op.LisemDir + "openlisem.ini");
	if (!fff.open(QIODevice::ReadOnly | QIODevice::Text))
		return;

	op.runfilename = fff.readLine();

	fff.close();
}
//---------------------------------------------------------------------------
void lisemqt::StorePath()
{
	QFile fff(op.LisemDir + "openlisem.ini");
	if (!fff.open(QIODevice::WriteOnly | QIODevice::Text))
		return;

	QTextStream ts( &fff );

	ts << op.runfilename;// << endl;

	fff.close();
}
//---------------------------------------------------------------------------
void lisemqt::on_toolButton_ShowRunfile_clicked()
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
//---------------------------------------------------------------------------
//--------------------------------------------------------------------



