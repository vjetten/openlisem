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
#include "model.h"
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


    QTextEdit *view = new QTextEdit(helptxt->toHtml());
    view->createStandardContextMenu();
    view->setWindowTitle("Option help");
    view->setMinimumWidth(qApp->desktop()->width()/3*2);
    view->setMinimumHeight(qApp->desktop()->height()/3*2);
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
    label_157->setEnabled(checked);
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
        dialog.setFileMode(QFileDialog::DirectoryOnly);
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
void lisemqt::on_checkDischargeUser_toggled(bool checked)
{
    groupDischargeUser->setEnabled(checked);
}
//--------------------------------------------------------------------
void lisemqt::on_checkWaveInUser_toggled(bool checked)
{
    groupWaveUser->setEnabled(checked);
}
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
void lisemqt::on_checkIncludeET_toggled(bool checked)
{
    radioGroupET->setEnabled(checked);
}
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
        maxWidth = qMax(maxWidth, QFontMetrics(QFont()).width(line));
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
 //   BaseflowParams->setEnabled(checked);
    if (checked) checkChannelInfil->setChecked(false);
}

void lisemqt::on_checkChannelInfil_toggled(bool checked)
{
   // BaseflowParams->setEnabled(!checked);
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
    BaseflowParams->setEnabled(checked);
    qDebug() << checked;
}
