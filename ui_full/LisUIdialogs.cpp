/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2024  Victor Jetten
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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#include "lisemqt.h"
//#include "model.h"
#include "global.h"


//---------------------------------------------------------------
void lisemqt::on_toolButton_resetCalibration_clicked()
{
    resetTabCalibration();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetFlow_clicked()
{
    resetTabFlow();
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_resetErosion_clicked()
{
    resetTabErosion();
}
//---------------------------------------------------------------

void lisemqt::on_toolButton_help1_clicked()
{
    on_toolButton_help(1);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help2_clicked()
{
    on_toolButton_help(2);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help3_clicked()
{
    on_toolButton_help(3);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help4_clicked()
{
    on_toolButton_help(4);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help5_clicked()
{
    on_toolButton_help(5);
}
//---------------------------------------------------------------
//void lisemqt::on_toolButton_help6_clicked()
//{
//    on_toolButton_help(6);
//}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help7_clicked()
{
    on_toolButton_help(7);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help8_clicked()
{
    on_toolButton_help(8);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help1a_clicked()
{
    on_toolButton_help(6);
}
//---------------------------------------------------------------
void lisemqt::on_toolButton_help(int page)
{
    QString filename;
    if (page == 1) filename=":/help1.html";
    if (page == 2) filename=":/help2.html";
    if (page == 3) filename=":/help3.html";
    if (page == 4) filename=":/help4.html";
    if (page == 5) filename=":/help5.html";
    if (page == 6) filename=":/help6.html";
    if (page == 7) filename=":/help7.html";
    if (page == 8) filename=":/help8.html";
    QFile file(filename);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream stream(&file);
    helptxt->setHtml(stream.readAll());
  //  helpbox->show();

    QScreen *screen = QGuiApplication::primaryScreen();
    QRect screenGeometry = screen->geometry();
    int screenWidth = screenGeometry.width();
    int screenHeight = screenGeometry.height();

    QTextEdit *view = new QTextEdit(helptxt->toHtml());
    view->createStandardContextMenu();
    view->setWindowTitle("Option help");

    view->setMinimumWidth(screenWidth/3*2);
    view->setMinimumHeight(screenHeight/3*2);
    view->setAttribute(Qt::WA_DeleteOnClose);

    view->show();
}
//--------------------------------------------------------------------

void lisemqt::on_check2DDiagonalFlow_toggled(bool checked)
{
    E_pitValue->setEnabled(checked);
    label_135->setEnabled(checked);
}
//--------------------------------------------------------------------

//void lisemqt::on_checkDiffusion_toggled(bool checked)
//{
//    E_SigmaDiffusion->setEnabled(checked);
//    label_101->setEnabled(checked);
//    label_139->setEnabled(checked);
//}
//--------------------------------------------------------------------

void lisemqt::on_checkHouses_toggled(bool checked)
{
    checkRaindrum->setEnabled(checked);
    label_78->setEnabled(checked);
}
//--------------------------------------------------------------------

// select a file or directory
// doFile = 0: select a directory;
// dofile = 1 select a file and return file name only;
// dofile = 2 return filename with full path
QString lisemqt::getFileorDir(QString inputdir,QString title, QStringList filters, int doFile)
{
    QFileDialog dialog;

    QString dirout = inputdir;

    if (doFile > 0) {
        dialog.setNameFilters(filters);
        dialog.setDirectory(QFileInfo(inputdir).absoluteDir());
        dialog.setFileMode(QFileDialog::ExistingFile);
    } else {
        filters.clear();
        dialog.setNameFilters(filters);
        dialog.setDirectory(QDir(inputdir));
        dialog.setFileMode(QFileDialog::Directory);
    }

    dialog.setLabelText(QFileDialog::LookIn,title);
    dialog.exec();
    if (dialog.selectedFiles().isEmpty())
        dirout = inputdir;
    else
        dirout = dialog.selectedFiles().at(0);

    if (doFile > 0) {
        if (doFile == 1)
            dirout = QFileInfo(dirout).fileName();
        if (doFile == 2)
            dirout = QFileInfo(dirout).absoluteFilePath();
    } else {
        dirout = dialog.selectedUrls().at(0).path();
        dirout.remove(0,1);
        if (dirout.lastIndexOf('/') != dirout.length())
            dirout = dirout + "/";
    }

    return dirout;
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_rainsatName_clicked()
{
    //RainSatFileDir = RainFileDir;
    if (!QFileInfo(RainSatFileDir).exists() || RainSatFileDir.isEmpty())
        RainSatFileDir = RainFileDir;
    if (!QFileInfo(RainSatFileDir).exists() || RainSatFileDir.isEmpty())
        RainSatFileDir = currentDir;
  //  qDebug() << RainSatFileDir << RainSatFileName << currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});
    QString sss = getFileorDir(RainSatFileDir,"Select rainfall map list table", filters, 2);

    RainSatFileDir = QFileInfo(sss).absolutePath()+"/";
    RainSatFileName = QFileInfo(sss).fileName(); //baseName();

    E_rainsatName->setText(RainSatFileDir + RainSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETsatName_clicked()
{
  //  ETSatFileDir = ETFileDir;

    if (!QFileInfo(ETSatFileDir).exists() || ETSatFileDir.isEmpty())
        ETSatFileDir = RainSatFileDir;
    if (!QFileInfo(ETSatFileDir).exists() || ETSatFileDir.isEmpty())
        ETSatFileDir = currentDir;
    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(ETSatFileDir,"Select ET map list table", filters, 2);

    ETSatFileDir = QFileInfo(sss).absolutePath()+"/";
    ETSatFileName = QFileInfo(sss).fileName(); //baseName();

    E_ETsatName->setText(ETSatFileDir + ETSatFileName);
}

//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallName_clicked()
{
    if (!QFileInfo(RainFileDir).exists() || RainFileDir.isEmpty())
        RainFileDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});
    QString sss = getFileorDir(RainFileDir,"Select rainfall station table", filters, 2);

    RainFileDir = QFileInfo(sss).absolutePath()+"/";
    RainFileName = QFileInfo(sss).fileName(); //baseName();

    E_RainfallName->setText(RainFileDir + RainFileName);

}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETName_clicked()
{
    if (!QFileInfo(ETFileDir).exists() || ETFileDir.isEmpty())
        ETFileDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(ETFileDir,"Select ET stations file", filters, 2);

    ETFileDir = QFileInfo(sss).absolutePath()+"/";
    ETFileName = QFileInfo(sss).fileName(); //baseName();
    E_ETName->setText(ETFileDir + ETFileName);
}
//--------------------------------------------------------------------
// void lisemqt::on_checkWaveInUser_toggled(bool checked)
// {
//     groupWaveUser->setEnabled(checked);
// }
//--------------------------------------------------------------------
void lisemqt::on_toolButton_DischargeShow_clicked()
{
    showTextfile(DischargeinDir + DischargeinFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_WaveShow_clicked()
{
    qDebug() <<WaveinDir + WaveinFileName;
    showTextfile(WaveinDir + WaveinFileName);
}
//--------------------------------------------------------------------
// void lisemqt::on_checkIncludeET_toggled(bool checked)
// {
//     widgetEToptions->setEnabled(checked);
// }
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETShow_clicked()
{
    //qDebug() <<ETFileDir + ETFileName;
    showTextfile(ETFileDir + ETFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_DischargeName_clicked()
{
    if (!QFileInfo(DischargeinDir).exists() || DischargeinDir.isEmpty())
        DischargeinDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(DischargeinDir,"Select ET stations file", filters, 2);

    DischargeinDir = QFileInfo(sss).absolutePath()+"/";
    DischargeinFileName = QFileInfo(sss).fileName(); //baseName();
    E_DischargeInName->setText(DischargeinDir + DischargeinFileName);

}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_WaveInName_clicked()
{
    if (!QFileInfo(WaveinDir).exists() || WaveinDir.isEmpty())
        WaveinDir = currentDir;

    QStringList filters({"Text file (*.txt *.tbl *.tss)","Any files (*)"});

    QString sss = getFileorDir(WaveinDir,"Select file with water heights", filters, 2);

    WaveinDir = QFileInfo(sss).absolutePath()+"/";
    WaveinFileName = QFileInfo(sss).fileName();
    E_WaveInName->setText(WaveinDir + WaveinFileName);

}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainfallShow_clicked()
{
    showTextfile(RainFileDir + RainFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_RainmapShow_clicked()
{
    showTextfile(RainSatFileDir + RainSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_ETmapShow_clicked()
{
    showTextfile(ETSatFileDir + ETSatFileName);
}
//--------------------------------------------------------------------
void lisemqt::showTextfile(QString name)
{
   // Read text from file
    QFile file(name);
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);
    QString initialText = in.readAll();
    file.close();

    // Find the longest line in the text
    QStringList lines = initialText.split("\n");
    int maxWidth = 0;
    for (const QString& line : lines) {
        maxWidth = qMax(maxWidth, QFontMetrics(QFont()).horizontalAdvance(line));
    }
    int maxHeight = qMin(qApp->primaryScreen()->size().height() * 2 / 3, lines.size() * 40);
    // Create a QDialog instance
    QDialog dialog;
    dialog.setWindowTitle(name);
    dialog.resize(maxWidth, maxHeight);//qApp->desktop()->width()/3*2,qApp->desktop()->height()/3*2);

    // Create a layout for the dialog
    QVBoxLayout *layout = new QVBoxLayout(&dialog);

    // Create a QTextEdit widget for editing text and load the initial text
    QTextEdit *textEdit = new QTextEdit(&dialog);
    textEdit->setPlainText(initialText); // Load initial text
    layout->addWidget(textEdit);

    // Create buttons layout
    QHBoxLayout *buttonLayout = new QHBoxLayout;

    // Create a QPushButton to save changes
    QPushButton *saveButton = new QPushButton("Save", &dialog);
    buttonLayout->addWidget(saveButton);

    // Create a QPushButton to cancel changes
    QPushButton *cancelButton = new QPushButton("Cancel", &dialog);
    buttonLayout->addWidget(cancelButton);

    QPushButton *closeButton = new QPushButton("Close", &dialog);
    buttonLayout->addWidget(closeButton);

    layout->addLayout(buttonLayout);

    // Connect the saveButton's clicked signal to a lambda function
    QObject::connect(saveButton, &QPushButton::clicked, [&]() {
        // Retrieve the text from the textEdit widget
        QString editedText = textEdit->toPlainText();
        // Open the file in write mode to save changes
        QFile writeFile(name);
        if (writeFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&writeFile);
            out << editedText;
            writeFile.close();
        }
        // Close the dialog when saving is done
        //dialog.close();
    });

    // Connect the cancelButton's clicked signal to close the dialog without saving changes
    QObject::connect(cancelButton, &QPushButton::clicked, [&]() {
        // Close the dialog without saving changes
        dialog.close();
    });

    QObject::connect(closeButton, &QPushButton::clicked, [&]() {
        dialog.close();
    });

    // Show the dialog
    dialog.show();

    dialog.exec();
}

void lisemqt::showTextfileOld(QString name)
{

    QFile file(name);
    if (!file.open(QFile::ReadWrite | QFile::Text))
    {
        QMessageBox::warning(this, QString("openLISEM"),
                             QString("Cannot read file %1:\n%2.")
                             .arg(name)
                             .arg(file.errorString()));
        return;
    }

    QString modifiedContents;
    QTextStream in(&file);

    QPlainTextEdit *view = new QPlainTextEdit(in.readAll());
    view->setWindowTitle(RainFileName);
    view->setMinimumWidth(400);
    view->setMinimumHeight(500);
    view->setAttribute(Qt::WA_DeleteOnClose);
    view->show();
    if (view->document()->isModified())
    {
        int ret =
                QMessageBox::question(this, QString("openLISEM"),
                                      QString("You have modified the contents of this file.\n"
                                              "Do you want to save it?"),
                                      QMessageBox::Ok |QMessageBox::Cancel,QMessageBox::Cancel);
        if (ret == QMessageBox::Ok)
        {
            // Don't take the address of a temporary!
            // in.setString(&view->toPlainText());
            modifiedContents = view->toPlainText();
            in.setString(&modifiedContents);
        }
    }

    file.close();
}

void lisemqt::on_E_EndTimeDay_returnPressed()
{
    int daye = E_EndTimeDay->text().split(":")[0].toInt();
    int mine = E_EndTimeDay->text().split(":")[1].toInt();
    daye = std::max(1,std::min(daye, 366));
    if (mine > 1440) {
        daye = mine/1440 + 1;
        mine = mine % 1440;
    }
    E_EndTimeDay->setText(QString("%1:%2").arg(daye,3,10,QLatin1Char('0')).arg(mine,4,10,QLatin1Char('0')));
}


void lisemqt::on_E_BeginTimeDay_returnPressed()
{
       int daye = E_BeginTimeDay->text().split(":")[0].toInt();
       int mine = E_BeginTimeDay->text().split(":")[1].toInt();
       daye = std::max(1,std::min(daye, 366));
       if (mine > 1440) {
           daye = mine/1440 + 1;
           mine = mine % 1440;
       }
       E_BeginTimeDay->setText(QString("%1:%2").arg(daye,3,10,QLatin1Char('0')).arg(mine,4,10,QLatin1Char('0')));
}



void lisemqt::on_toolButton_resetInfiltration_clicked()
{
    resetTabInfiltration();
}

void lisemqt::on_toolButton_resetInterception_clicked()
{
    resetTabInterception();
}

void lisemqt::on_toolButton_resetAdvanced_clicked()
{
    resetTabAdvanced();
}

void lisemqt::on_toolButton_resetOptions_clicked()
{
    resetTabOptions();
}

void lisemqt::on_checkStationaryBaseflow_toggled(bool checked)
{
    if (checked) checkChannelInfil->setChecked(false);
   // doChannelBaseflow = checked;
}

void lisemqt::on_checkChannelInfil_toggled(bool checked)
{
    if (checked) checkStationaryBaseflow->setChecked(false);
}

void lisemqt::on_E_EfficiencyDETCH_currentIndexChanged(int index)
{
    E_EfficiencyDirect->setEnabled(index == 3);
}

void lisemqt::on_checkGWflow_toggled(bool checked)
{
    GW_widget->setEnabled(checked);
    widget_GWparams->setEnabled(checked);
    groupBaseflowParams->setEnabled(checked);
    //qDebug() << checked;
}

//--------------------------------------------------------------------
void lisemqt::on_E_floodMinHeight_valueChanged(double)
{
    label_107->setText(QString("Flood (mm),h>%1)").arg(E_floodMinHeight->value()*1000));
    label_40->setText(QString("Runoff (mm),h<%1)").arg(E_floodMinHeight->value()*1000));
}
//--------------------------------------------------------------------
// this is for the directory with the table files
void lisemqt::on_toolButton_SwatreTableDir_clicked()
{
    QString path;
    QString pathin;

    pathin = findValidDir(E_SwatreTableDir->text(), false);

    path = QFileDialog::getExistingDirectory(this, QString("Select the directory with the Swatre tables"),
                                             pathin,
                                             QFileDialog::ShowDirsOnly
                                             | QFileDialog::DontResolveSymlinks);
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
                                        QString("Select the SWATRE profile definition file"),
                                        SwatreTableName,"Profiles (*.inp);;All files (*.*)");
    if(!path.isEmpty())
    {
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
void lisemqt::on_toolButton_satImageName_clicked()
{
    QString path;

    satImageFileDir = findValidDir(satImageFileDir, false);

    path = QFileDialog::getOpenFileName(this,
                                        QString("Select background satellite image file"),
                                        satImageFileDir,"GeoTiff (*.tif)");
    if(!path.isEmpty())
    {
        QFileInfo fi(path);
        satImageFileName = fi.fileName();
        satImageFileDir = CheckDir(fi.absolutePath());//Dir().path());
        E_satImageName->setText( satImageFileDir + satImageFileName );
    }
}
//---------------------------------------------------------------------------
void lisemqt::on_checkRainfall_toggled(bool checked)
{
    groupRainfall->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkET_toggled(bool checked)
{
    groupET->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkInterception_toggled(bool checked)
{
    groupInterception->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_E_OFWaveType_currentIndexChanged(int index)
{
    groupFloodParams->setEnabled(index > 0);
    groupWaveUser->setEnabled(index > 0);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkInfiltration_toggled(bool checked)
{
    groupInfiltration->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_toggled(bool checked)
{
    groupChannelParams->setEnabled(checked);
    checkMapChannels->setEnabled(checked);

    checkMapNameModel(CHANNELMAPS, 0, checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkDoErosion_toggled(bool checked)
{
    widgetErosion->setEnabled(checked);

    groupConservationSed->setEnabled(checked);

    ComboMinSpinBox2->setEnabled(checked);
    ComboMaxSpinBox2->setEnabled(checked);
    DisplayComboBox2->setEnabled(checked);
    checkBoxComboMaps2->setEnabled(checked);

    outputMapsSediment->setEnabled(checked);
    groupSedMapseriesout->setEnabled(checked);
    tabWidget_totout->setTabEnabled(1,checked);

    groupCalErosion->setEnabled(checked);

    label_soillosskgha->setEnabled(checked);
    label_soilloss->setEnabled(checked);
    label_SDR->setEnabled(checked);
    label_94->setEnabled(checked);
    label_105->setEnabled(checked);
    label_45->setEnabled(checked);

    checkMapNameModel(EROSIONMAPS, 0, checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkInfrastructure_toggled(bool checked)
{
    widgetInfra->setEnabled(checked);
    transparencyHardSurface->setEnabled(checked);
    transparencyRoad->setEnabled(checked);
    checkMapHardSurface->setEnabled(checked);
    checkMapBuildings->setEnabled(checked);
    checkMapRoads->setEnabled(checked);
}
//---------------------------------------------------------------------------
void lisemqt::on_checkConservation_toggled(bool checked)
{
    groupMitigationWater->setEnabled(checked);
    groupConservationSed->setEnabled(checked);

}
//---------------------------------------------------------------------------
void lisemqt::on_toolButton_ShowRunfile_clicked()
{
    //qDebug() << E_runFileList->currentText();
    showTextfile(E_runFileList->currentText());
}
