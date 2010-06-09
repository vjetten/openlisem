/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * This is the extensive interface, with toolbar, initialization of graph
 * and options
 */

#include "lisemqt.h"
#include "model.h"
#include "global.h"

output op;

//--------------------------------------------------------------------
lisemqt::lisemqt(QWidget *parent)
: QMainWindow(parent)
{
	setupUi(this);
	// set up interface
	resize(856, 674);

	MapNameModel = NULL;
	HPlot = NULL;
	LisemReset();
	SetToolBar();
	FillMapList();
	// initalize interface

	W = NULL;
	// initalize pointer to the world, created when run button is pushed

	GetStorePath();
	// get last place visited by opneLISEM and go there

	SetStyleUI();
	// do some style things

	SetGraph();

}
//--------------------------------------------------------------------
lisemqt::~lisemqt()
{
	StorePath();

	if (HPlot)
		delete HPlot;
	//delete QGraph;
	//delete QsGraph;
	//delete CGraph;

}
//--------------------------------------------------------------------
void lisemqt::SetToolBar()
{
	openAct = new QAction(QIcon(":/fileopen.png"), "&Open...", this);
	openAct->setShortcuts(QKeySequence::Open);
	openAct->setStatusTip("Open a run file");
	connect(openAct, SIGNAL(triggered()), this, SLOT(openRunFile()));
	toolBar->addAction(openAct);

	saveAct = new QAction(QIcon(":/filesave.png"), "&Save...", this);
	saveAct->setShortcuts(QKeySequence::Save);
	saveAct->setStatusTip("Save a run file");
	connect(saveAct, SIGNAL(triggered()), this, SLOT(SaveRunFile()));
	toolBar->addAction(saveAct);

	saveasAct = new QAction(QIcon(":/filesaveas.png"), "Save &As...", this);
	saveasAct->setShortcuts(QKeySequence::SaveAs);
	saveasAct->setStatusTip("Save a run file as ...");
	connect(saveasAct, SIGNAL(triggered()), this, SLOT(savefileas()));
	toolBar->addAction(saveasAct);

	runAct = new QAction(QIcon(":/start1.png"), "Run model...", this);
	//	runAct->setShortcuts(QKeySequence(QString("Ctrl+R")));
	runAct->setStatusTip("run the model ...");
	connect(runAct, SIGNAL(triggered()), this, SLOT(runmodel()));
	toolBar->addAction(runAct);

	stopAct = new QAction(QIcon(":/stop16_2.png"), "Stop the model...", this);
	//	runAct->setShortcuts(QKeySequence(Qt::CTRL + Qt::Key_R));
	stopAct->setStatusTip("stop the model run ...");
	connect(stopAct, SIGNAL(triggered()), this, SLOT(stopmodel()));
	toolBar->addAction(stopAct);
	toolBar->addSeparator();

}
//--------------------------------------------------------------------
void lisemqt::LisemReset()
{
	DefaultMapnames();

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

	RunFileNames.clear();

	RainFileName.clear();
	rainFileDir.clear();
	SnowmeltFileName.clear();
	snowmeltFileDir.clear();

}
//---------------------------------------------------------------------------
void lisemqt::SetGraph()
{
	QwtText title;
	title.setText("Hydrograph/Sedigraph outlet");
	HPlot = new QwtPlot(title, widgetGraph);
	// make the plot window

	PGraph = new QwtPlotCurve("Rainfall");
	QGraph = new QwtPlotCurve("Discharge");
	QsGraph = new QwtPlotCurve("Sediment discharge");
	CGraph = new QwtPlotCurve("Concentration");
	PGraph->attach(HPlot);
	QGraph->attach(HPlot);
	if(!checkNoErosion->isChecked())
	{
		QsGraph->attach(HPlot);
		CGraph->attach(HPlot);
	}
	// order determines order of display in Legend
	PGraph->setAxis(HPlot->xBottom, HPlot->yLeft);
	QGraph->setAxis(HPlot->xBottom, HPlot->yLeft);
	QsGraph->setAxis(HPlot->xBottom, HPlot->yRight);
	CGraph->setAxis(HPlot->xBottom, HPlot->yRight);
	QColor col;
	col.setRgb( 200,0,0,255 );
	CGraph->setPen(QPen(col));
	QsGraph->setPen(QPen(Qt::red));
	col.setRgb( 60,100,160,255 );
	QGraph->setPen(QPen(col));
	PGraph->setPen(QPen("#000000"));
	PGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
	QGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
	QsGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
	CGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
	// make all graphs to be drawn and link them to HPlot
	// set colors

	QwtLegend *legend = new QwtLegend(widgetGraph);
	legend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
	HPlot->insertLegend(legend, QwtPlot::BottomLegend);

	//legend

	HPlot->resize(450,380);
	HPlot->setCanvasBackground("#FFFFFF");
	// size and white graph

	HPlot->enableAxis(HPlot->yRight,true);
	HPlot->enableAxis(HPlot->yLeft,true);
	HPlot->enableAxis(HPlot->xBottom,true);
	HPlot->setAxisTitle(HPlot->xBottom, "time (min)");
	HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s)/P (mm/h)");
	HPlot->setAxisTitle(HPlot->yRight, "Qs(kg/s)/C(g/l)");
	HPlot->setAxisScale(HPlot->yRight, 0, 1);
	HPlot->setAxisScale(HPlot->yLeft, 0, 100);
	HPlot->setAxisScale(HPlot->xBottom, 0, 100);

	// set axes

	QwtPlotGrid *grid = new QwtPlotGrid();
	grid->enableXMin(true);
	grid->enableYMin(true);
	col.setRgb( 180,180,180,180 );
	grid->setMajPen(QPen(col, 0, Qt::DashLine));
	col.setRgb( 210,210,210,180 );
	grid->setMinPen(QPen(col, 0 , Qt::DotLine));
	grid->attach(HPlot);
	// set gridlines


	HPlot->replot();
	// draw empty plot

	QData = NULL; //discharge
	QsData = NULL;  //sed discharge
	CData = NULL; //conc
	PData = NULL; //rainfall
	timeData = NULL;  //time
	// init data arrays for plot data

}
//---------------------------------------------------------------------------
void lisemqt::SetStyleUI()
{
	// make some labels yellow

	label_dx->setStyleSheet("* { background-color: #ffffff }");
	label_area->setStyleSheet("* { background-color: #ffffff }");
	label_time->setStyleSheet("* { background-color: #ffffff }");
	label_endtime->setStyleSheet("* { background-color: #ffffff }");
	label_raintot->setStyleSheet("* { background-color: #ffff77 }");
	label_watervoltot->setStyleSheet("* { background-color: #ffff77 }");
	label_qtot->setStyleSheet("* { background-color: #ffff77 }");
	label_infiltot->setStyleSheet("* { background-color: #ffff77 }");
	label_surfstor->setStyleSheet("* { background-color: #ffff77 }");
	label_interctot->setStyleSheet("* { background-color: #ffff77 }");
	label_qtotm3->setStyleSheet("* { background-color: #ffff77 }");
	label_qpeak->setStyleSheet("* { background-color: #ffff77 }");
	label_qpeaktime->setStyleSheet("* { background-color: #ffff77 }");
	label_ppeaktime->setStyleSheet("* { background-color: #ffff77 }");
	label_QPfrac->setStyleSheet("* { background-color: #ffff77 }");
	label_discharge->setStyleSheet("* { background-color: #ffff77 }");

	label_splashdet->setStyleSheet("* { background-color: #ffff77 }");
	label_flowdet->setStyleSheet("* { background-color: #ffff77 }");
	label_sedvol->setStyleSheet("* { background-color: #ffff77 }");
	label_dep->setStyleSheet("* { background-color: #ffff77 }");
	label_detch->setStyleSheet("* { background-color: #ffff77 }");
	label_depch->setStyleSheet("* { background-color: #ffff77 }");
	label_sedvolch->setStyleSheet("* { background-color: #ffff77 }");
	label_soilloss->setStyleSheet("* { background-color: #ffff77 }");
	label_soillosskgha->setStyleSheet("* { background-color: #ffff77 }");

	label_buffervol->setStyleSheet("* { background-color: #ffff77 }");
	label_buffersed->setStyleSheet("* { background-color: #ffff77 }");
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_fileOpen_clicked()
{
	openRunFile();
	GetRunfile();
	ParseInputData();
	FillMapList();
	RunAllChecks();
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_MapDir_clicked()
{
	QString path;
	path = QFileDialog::getExistingDirectory(this,
			tr("Select maps directory"),
			E_MapDir->text(),
			QFileDialog::ShowDirsOnly);
	if(!path.isEmpty())
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

	if(!path.isEmpty())
		E_ResultDir->setText( path );
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallName_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
			tr("Select rainfall file"),
			QString::null);
	if(!path.isEmpty())
	{
		E_RainfallName->setText( path );
		RainFileName = path;
	}
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltName_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
			tr("Select snow melt file"),
			QString::null);
	if(!path.isEmpty())
	{
		E_SnowmeltName->setText( path );
		SnowmeltFileName = path;
	}
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltShow_clicked()
{
	QFile file(SnowmeltFileName);
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this,"openLISEM",
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
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this, QString("openLISEM"),
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
void lisemqt::savefileas()
{
	if (op.runfilename.isEmpty())
	{
		QMessageBox::warning(this, "openLISEM","Select a runfile first.");
		return;
	}

	QString selectedFilter;
	QString fileName = QFileDialog::getSaveFileName(this,
			tr("Give a new runfile name"),
			op.runfilename,
			tr("Text Files (*.run);;All Files (*)"),
			&selectedFilter);
	//options);
	if (!fileName.isEmpty())
		savefile(fileName);

}
//--------------------------------------------------------------------
void lisemqt::SaveRunFile()
{
	savefile(op.runfilename);
}
//--------------------------------------------------------------------
void lisemqt::savefile(QString name)
{
	UpdateModelData();
	// change runfile strings with current interface options

	QFile fp(name);
	if (!fp.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		QMessageBox::warning(this, QString("openLISEM"),
				tr("Cannot write file %1:\n%2.").arg(name).arg(fp.errorString()));
		return;
	}

	QTextStream out(&fp);
	out << QString("[openLISEM runfile version 4]\n");
	for (int i = 1; i < nrnamelist; i++)
	{
		if (namelist[i].name.contains("["))
			out << "\n" << namelist[i].name; // already contains \n
		else
			out << namelist[i].name << "=" << namelist[i].value << "\n";
	}
	fp.close();
}
//--------------------------------------------------------------------
void lisemqt::openRunFile()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
			QString("Select run file(s)"),
			currentDir,
			QString("*.run"));

	if (path.isEmpty())
		return;
	E_runFileList->setInsertPolicy(QComboBox::InsertAtTop);

	bool exst = false;
	for (int i = 0; i < E_runFileList->count(); i++)
		if (E_runFileList->itemText(i) == path)
			exst = true;
	if (!exst)
		E_runFileList->insertItem(0,path);
	E_runFileList->setCurrentIndex(0);

	RunFileNames.clear();
	for (int i = 0; i <= E_runFileList->count(); i++)
		RunFileNames << E_runFileList->itemText(i);

	RunFileNames.removeDuplicates();
	op.runfilename = E_runFileList->itemText(0);
}
//---------------------------------------------------------------------------
void lisemqt::GetStorePath()
{
	QFile fff(op.LisemDir + "openlisem.ini");
	if (!fff.open(QIODevice::ReadOnly | QIODevice::Text))
		return;
	QFileInfo fi(fff.readLine());
	QDir dir = fi.absoluteDir();
	currentDir = dir.absolutePath();
	//label_debug->setText(currentDir);
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
void lisemqt::on_E_runFileList_currentIndexChanged(int)
{
	CurrentRunFile = E_runFileList->currentIndex();
	op.runfilename = E_runFileList->currentText();
	//RunFileNames.at(CurrentRunFile);
	GetRunfile();
	ParseInputData();
	FillMapList();
	RunAllChecks();
}
//--------------------------------------------------------------------
void lisemqt::on_E_MapDir_textEdited()
{
	QFileInfo fin(E_MapDir->text());
	if(!fin.exists())
	{
		E_MapDir->setText("");
		QMessageBox::warning(this,"openLISEM",
				QString("Map directory does not exist"));
	}
}
//--------------------------------------------------------------------
void lisemqt::on_E_ResultDir_textEdited()
{
	if (E_ResultDir->text().isEmpty())
		return;
	QFileInfo fin(E_ResultDir->text());
	if(!fin.exists())
	{
		int ret = QMessageBox::question(this, QString("openLISEM"),
				QString("The directory \"%1\"does not exist.\n"
						"Do you want to create it (apply)?").arg(fin.absoluteFilePath()),
						QMessageBox::Apply |QMessageBox::Cancel,
						QMessageBox::Cancel);
		if (ret == QMessageBox::Apply)
			QDir(E_ResultDir->text()).mkpath(E_ResultDir->text());

	}
}


