/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
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
//	MapPlot = NULL;
	resetAll();
	//LisemReset();

	SetToolBar();
	FillMapList();
	// initalize interface and make tree structure for map names

   DefaultRunFile();
   // fill defnamelist with default runfile names

   SetConnections();

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
//	if (MapPlot)
//		delete MapPlot;
	//delete QGraph;
	//delete QsGraph;
	//delete CGraph;

}
//--------------------------------------------------------------------
void lisemqt::SetConnections()
{
    connect(checkRainfall, SIGNAL(toggled(bool)), this, SLOT(doCheckRainfall(bool)));
    connect(checkSnowmelt, SIGNAL(toggled(bool)), this, SLOT(doCheckSnowmelt(bool)));
    connect(toolButton_fileOpen, SIGNAL(clicked()), this, SLOT(openRunFile()));

    connect(treeView, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(doOpenMapname(QModelIndex)));
    // double click on mapnake opens fileopen
    //connect(MapNameModel, SIGNAL(dataChanged(QModelIndex, QModelIndex)), this, SLOT(doChangeMapname(QModelIndex, QModelIndex)));
    // doubleclick on mapname edits mapname
}
//--------------------------------------------------------------------
void lisemqt::SetToolBar()
{
	restartAct = new QAction(QIcon(":/Undo-icon.png"), "&Reset...", this);
	connect(restartAct, SIGNAL(triggered()), this, SLOT(resetAll()));
	toolBar->addAction(restartAct);
	toolBar->addSeparator();

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
	toolBar->addSeparator();

	shootscreenAct = new QAction(QIcon(":/screenshots.png"), "Stop the model...", this);
	//	runAct->setShortcuts(QKeySequence(Qt::CTRL + Qt::Key_R));
	shootscreenAct->setStatusTip("make a screendump ...");
	connect(shootscreenAct, SIGNAL(triggered()), this, SLOT(shootScreen()));
	toolBar->addAction(shootscreenAct);
	toolBar->addSeparator();

	//runAct = new QAction(QIcon(":/Button-Play-icon.png"), "Run model...", this);
	runAct = new QAction(QIcon(":/start1.png"), "Run model...", this);
	//	runAct->setShortcuts(QKeySequence(QString("Ctrl+R")));
	runAct->setStatusTip("run the model ...");
	connect(runAct, SIGNAL(triggered()), this, SLOT(runmodel()));
	toolBar->addAction(runAct);

	//	pauseAct = new QAction(QIcon(":/Button-Pause-icon.png"), "Pause the model...", this);
	pauseAct = new QAction(QIcon(":/pause2.png"), "Pause the model...", this);
	pauseAct->setStatusTip("pause the model run ...");
	connect(pauseAct, SIGNAL(triggered()), this, SLOT(pausemodel()));
	toolBar->addAction(pauseAct);

	//stopAct = new QAction(QIcon(":/Button-Stop-icon.png"), "Stop the model...", this);
	stopAct = new QAction(QIcon(":/stop1.png"), "Stop the model...", this);
	//	runAct->setShortcuts(QKeySequence(Qt::CTRL + Qt::Key_R));
	stopAct->setStatusTip("stop the model run ...");
	connect(stopAct, SIGNAL(triggered()), this, SLOT(stopmodel()));
	toolBar->addAction(stopAct);

	aboutActI = new QAction(QIcon(":/Info.png"), "", this);
	connect(aboutActI, SIGNAL(triggered()), this, SLOT(aboutInfo()));
	toolBar_2->addAction(aboutActI);

	//aboutAct = new QAction(QIcon(":/Info.png"), "", this);
	//connect(aboutAct, SIGNAL(triggered()), this, SLOT(aboutQT()));
	//toolBar_2->addAction(aboutAct);

	toolBar_2->setMovable( false);
	toolBar->setMovable( false);

}
//--------------------------------------------------------------------
void lisemqt::LisemReset()
{
	DefaultMapnames();

	E_InfiltrationMethod->addItem("no Infiltration");
	E_InfiltrationMethod->addItem("SWATRE");
	E_InfiltrationMethod->addItem("Green and Ampt");
	E_InfiltrationMethod->addItem("Smith and Parlange");
	E_InfiltrationMethod->addItem("Subtract Ksat");

	RunFileNames.clear();

	RainFileName.clear();
	RainFileDir.clear();
	SnowmeltFileName.clear();
	SnowmeltFileDir.clear();
	SwatreTableName.clear();
	SwatreTableDir.clear();
}
//---------------------------------------------------------------------------
void lisemqt::SetMapPlot()
{
    /*
	QwtText title;
	title.setText("something");
	MapDrawing = new QwtPlotSpectrogram();
	MapPlot = new QwtPlot(title, widgetMap);
	// make the plot window and link it to the histogram

	MapDrawData = new SpectrogramData();//QwtRasterData();
	//mapDrawData->DMap->

	QwtLinearColorMap colorMap(Qt::darkCyan, Qt::red);
	colorMap.addColorStop(0.1, Qt::cyan);
	colorMap.addColorStop(0.6, Qt::green);
	colorMap.addColorStop(0.95, Qt::yellow);

	MapDrawing->setColorMap(colorMap);
	MapDrawing->setDisplayMode(QwtPlotSpectrogram::ImageMode);
	MapDrawing->attach(MapPlot);

   //MapDrawing->setData(MapDrawData);



   //	MapPlot->replot();
   */
}
//---------------------------------------------------------------------------
void lisemqt::SetGraph()
{
	textGraph->setMaximumBlockCount(6);
	textGraph->setWordWrapMode(QTextOption::NoWrap);
	textGraph->setMaximumHeight(96);

	QwtText title;
	title.setText("Hydrograph/Sedigraph outlet");
	HPlot = new QwtPlot(title, this);//, widgetGraph);
	// make the plot window
	verticalLayout_6->insertWidget(0, HPlot);

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

	QwtLegend *legend = new QwtLegend(HPlot);//this);//widgetGraph);
	legend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
	HPlot->insertLegend(legend, QwtPlot::BottomLegend);

	//legend

	//HPlot->resize(450,380);
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
	label_SDR->setStyleSheet("* { background-color: #ffff77 }");

	label_buffervol->setStyleSheet("* { background-color: #ffff77 }");
	label_buffersed->setStyleSheet("* { background-color: #ffff77 }");
}
//--------------------------------------------------------------------
//void lisemqt::on_toolButton_fileOpen_clicked()
//{
//	openRunFile();
//}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_MapDir_clicked()
{
	QString path;
	QString pathin;

	pathin = E_MapDir->text();
	if (pathin.isEmpty())
		pathin = currentDir;

	path = QFileDialog::getExistingDirectory(this,
	QString("Select maps directory"),
	pathin,
	QFileDialog::ShowDirsOnly);
	if(!path.isEmpty())
		E_MapDir->setText( path );
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ResultDir_clicked()
{
	QString path;
	QString pathin;

	pathin = E_ResultDir->text();
	if (pathin.isEmpty())
		pathin = currentDir;

	path = QFileDialog::getExistingDirectory(this,
	QString("Select or create a directory to write results"),
	pathin,
	QFileDialog::ShowDirsOnly);

	if(!path.isEmpty())
		E_ResultDir->setText( path );
}
//--------------------------------------------------------------------
// this is for the directory with the table files
void lisemqt::on_toolButton_SwatreTableDir_clicked()
{
	QString path;
	QString pathdir = E_SwatreTableDir->text() + "/..";
	path = QFileDialog::getExistingDirectory(this,
	QString("Select the directory with the Swatre profile tables"),
	pathdir,
	QFileDialog::HideNameFilterDetails);//ShowDirsOnly);

	if(!path.isEmpty())
	{
		E_SwatreTableDir->setText( path );
		SwatreTableDir = path;
	}
}
//--------------------------------------------------------------------
// this is for the file profile.inp
void lisemqt::on_toolButton_SwatreTableFile_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
	QString("Select SWATRE table"),
	SwatreTableName);
	if(!path.isEmpty())
	{

		QFileInfo fi(path);
		SwatreTableName = path;
		E_SwatreTableName->setText(path);
	}
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SwatreTableShow_clicked()
{
	QFile file(SwatreTableName);
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this,"openLISEM",
		QString("Cannot read file %1:\n%2.")
		.arg(SwatreTableName)
		.arg(file.errorString()));
		return;
	}

	QTextStream in(&file);

	QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
	view->setWindowTitle(SwatreTableName);
	view->setMinimumWidth(400);
	view->setMinimumHeight(500);
	view->setAttribute(Qt::WA_DeleteOnClose);
	view->show();

	file.close();
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallName_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
	QString("Select rainfall file"),
	RainFileDir);
	//QString::null);
	if(!path.isEmpty())
	{
		QFileInfo fi(path);
		RainFileName = fi.fileName();
		RainFileDir = CheckDir(fi.absoluteDir().path());
		E_RainfallName->setText( RainFileName );
	}
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltName_clicked()
{
	QString path;
	path = QFileDialog::getOpenFileName(this,
	QString("Select snow melt file"),
	SnowmeltFileDir);
	if(!path.isEmpty())
	{
		QFileInfo fi(path);
		SnowmeltFileName = fi.fileName();
		SnowmeltFileDir = CheckDir(fi.absoluteDir().path());
		E_SnowmeltName->setText( SnowmeltFileName );
	}
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SnowmeltShow_clicked()
{
	QFile file(SnowmeltFileDir + SnowmeltFileName);
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this,"openLISEM",
		QString("Cannot read file %1:\n%2.")
		.arg(SnowmeltFileDir + SnowmeltFileName)
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

	QFile file(RainFileDir + RainFileName);
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this, QString("openLISEM"),
		QString("Cannot read file %1:\n%2.")
		.arg(RainFileDir + RainFileName)
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
	QString("Give a new runfile name"),
	op.runfilename,
	QString("Text Files (*.run);;All Files (*)"),
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
		QString("Cannot write file %1:\n%2.").arg(name).arg(fp.errorString()));
		return;
	}

	QTextStream out(&fp);
	out << QString("[openLISEM runfile version 1.0]\n");
	for (int i = 1; i < nrdefnamelist; i++)
	{
		if (defnamelist[i].name.contains("[") || defnamelist[i].name.isEmpty())
			out << defnamelist[i].name << "\n"; // already contains \n
		else
			out << defnamelist[i].name << "=" << defnamelist[i].value << "\n";
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
	// this triggers a runfile in on_E_runFileList_currentIndexChanged

	RunFileNames.clear();
	for (int i = 0; i <= E_runFileList->count(); i++)
		RunFileNames << E_runFileList->itemText(i);

	RunFileNames.removeDuplicates();
	op.runfilename = E_runFileList->itemText(0);

	/* this is done in  E_runFileList change
 GetRunfile();
 ParseInputData();
 FillMapList();
 RunAllChecks();
 */
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
	if (E_runFileList->count() == 0)
		return;
	CurrentRunFile = E_runFileList->currentIndex();
	op.runfilename = E_runFileList->currentText();
	//RunFileNames.at(CurrentRunFile);

	GetRunfile();
	ParseInputData();
	FillMapList();  // fill the tree strcuture on page 2
	RunAllChecks(); // activate the maps in the tree parts in response to checks
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
		int ret =
		QMessageBox::question(this, QString("openLISEM"),
		QString("The directory \"%1\"does not exist.\n"
		"Do you want to create it (apply)?")
		.arg(fin.absoluteFilePath()),
		QMessageBox::Apply |QMessageBox::Cancel,QMessageBox::Cancel);
		if (ret == QMessageBox::Apply)
			QDir(E_ResultDir->text()).mkpath(E_ResultDir->text());

	}
}
//--------------------------------------------------------------------
void lisemqt::shootScreen()
{
	if (op.runfilename.isEmpty())
	{
		QMessageBox::warning(this, "openLISEM",QString("Select a run file first"));
		return;
	}
	QPixmap originalPixmap; // clear image for low memory situations
	// on embedded devices.
	originalPixmap = QPixmap::grabWidget(tabWidget->currentWidget());

	QString format = "png";
	QFileInfo fi(op.runfilename);

	QString fileName = CheckDir(E_ResultDir->text()) + fi.baseName() + "." + format;

	originalPixmap.save(fileName, format.toAscii());
}
//--------------------------------------------------------------------
void lisemqt::aboutQT()
{
	QMessageBox::aboutQt ( this, "openLISEM" );
}
//--------------------------------------------------------------------
void lisemqt::aboutInfo()
{
	QMessageBox::information ( this, "openLISEM",
	QString("openLISEM is created wih:\n\n%1\n%2\n%3\n%4\n%5\n\n%6")
	.arg("- MingW C/C++ compiler (http://www.mingw.org);")
	.arg("- Qt cross platform application and UI framework (http://qt.nokia.com/)")
	.arg("- QtCreator IDE (http://qt.nokia.com/)")
	.arg("- Qwt technical application widgets for Qt (http://qwt.sf.net)")
	.arg("- Tortoise SVN for version control: (http://tortoisesvn.net/)")
	.arg("Details can be found at: http://sourceforge.net/projects/lisem/")
	);
}
//--------------------------------------------------------------------
void lisemqt::resetAll()
{
	E_runFileList->clear();

	E_InfiltrationMethod->clear();
	E_InfiltrationMethod->addItem("no Infiltration");
	E_InfiltrationMethod->addItem("SWATRE");
	E_InfiltrationMethod->addItem("Green and Ampt");
	E_InfiltrationMethod->addItem("Smith and Parlange");
	E_InfiltrationMethod->addItem("Subtract Ksat");

	DefaultMapnames();
	RunFileNames.clear();
	op.runfilename.clear();

	E_MapDir->setText("");
	E_RainfallName->setText("");
	E_SnowmeltName->setText("");
	E_ResultDir->setText("");
	E_DetachmentMap->setText("");
	E_DepositionMap->setText("");
	E_SoillossMap->setText("");
	E_MainTotals->setText("");
	E_PointResults->setText("");

	E_BeginTime->setText("");
	E_EndTime->setText("");
	E_Timestep->setText("");

	checkBox_OutRunoff->setChecked(false);
	checkBox_OutConc->setChecked(false);
	checkBox_OutWH->setChecked(false);
	checkBox_OutWHC->setChecked(false);
	checkBox_OutTC->setChecked(false);
	checkBox_OutDet->setChecked(false);
	checkBox_OutDep->setChecked(false);
	checkBox_OutV->setChecked(false);
	checkBox_OutInf->setChecked(false);
	checkBox_OutSurfStor->setChecked(false);
	checkBox_OutChanVol->setChecked(false);

	printinterval->setValue(1);


	E_InfiltrationMethod->setCurrentIndex(0);

	InitOP();
	progressBar->setValue(0);

	bool check = false;
	checkNoErosion->setChecked(check);
	checkIncludeChannel->setChecked(check);
	checkChannelInfil->setChecked(check);
	checkChannelBaseflow->setChecked(check);
	//	checkAllinChannel->setChecked(check);
	checkSnowmelt->setChecked(check);
	checkRainfall->setChecked(true);
	checkAltErosion->setChecked(check);
	checkSimpleDepression->setChecked(check);
	checkHardsurface->setChecked(check);
	checkBuffers->setChecked(check);
	checkSedtrap->setChecked(check);
	checkInfilCompact->setChecked(check);
	checkInfilGrass->setChecked(check);
	checkInfilCrust->setChecked(check);
	checkImpermeable->setChecked(check);
	//	checkDumphead->setChecked(check);
	checkGeometric->setChecked(true);
	//	checkRunoffPerM->setChecked(check);
	checkWritePCRnames->setChecked(true);
	checkWritePCRtimeplot->setChecked(check);
	checkOutputTimeStep->setChecked(true);
	checkOutputTimeUser->setChecked(check);
	checkNoErosionOutlet->setChecked(check);
	//	checkDrainage->setChecked(check);
	//	checkGullyInfil->setChecked(check);
	//	checkGullyInit->setChecked(check);
	checkSeparateOutput->setChecked(check);
	checkSOBEKOutput->setChecked(check);
	SOBEKdatestring->setText("10/01/01");
	//checkInterceptionLAI->setChecked(true);
	E_BulkDens->setText("1200.00");

	tabWidget->setCurrentIndex(0);

}
//--------------------------------------------------------------------
