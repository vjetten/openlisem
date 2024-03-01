#include "lisemqt.h"

//---------------------------------------------------------------------------
lismpeg::lismpeg(QWidget *parent) :
        QDialog(parent)
{
    setupUi(this);
    this->setWindowTitle("Create mpeg from screenshots");
   // this->setWindowFlags(Qt::WindowSystemMenuHint | Qt::WindowTitleHint);

  //  mencoderDir = QCoreApplication::applicationDirPath() + "/mencoder.exe";
  //  E_mencoderDir->setText(mencoderDir);
    resize(qApp->primaryScreen()->size().height()*2/3,qApp->primaryScreen()->size().height()*2/3);


    toolButton_createMP4->setStyleSheet("QToolButton { text-align: center; }");
    QStringList sss;
    sss << "2560x1440" << "1920x1080" << "1280x720";
    comboBox_mpegResolution->addItems(sss);
    comboBox_mpegResolution->setCurrentIndex(1);

}
//---------------------------------------------------------------------------
lismpeg::~lismpeg()
{
  //  delete ui;
}


// select a file or directory
// doFile = 0: select a directory;
// dofile = 1 select a file and return file name only;
// dofile = 2 return filename with full path
QString lismpeg::getFileorDir(QString inputdir,QString title, QStringList filters, int doFile)
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

void lismpeg::setMencoderDir(QString d)
{
    mencoderDir = d;
    E_mencoderDir->setText(mencoderDir);
}

// QString lismpeg::getMencoderDir()
// {
//     retrun(mencoderDir);
//     reurn E_resultsDir->setText(mencoderDir);
// }

void lismpeg::setWorkDir(QString d)
{
    resultsDir = d;
    E_resultsDir->setText(resultsDir);
}
//---------------------------------------------------------------------------

void lismpeg::on_toolButton_resultDir_clicked()
{
    QStringList filters({"Any files (*)"});
    QString sss = getFileorDir(resultsDir,"Select result dir", filters, 0);

    resultsDir = QFileInfo(sss).absolutePath()+"/";

    E_resultsDir->setText(resultsDir);

}


void lismpeg::on_toolButton_mencoderDir_clicked()
{
    QStringList filters({"Exe files (*.exe)"});
    QString sss = getFileorDir(mencoderDir,"Select mencoder.exe", filters, 2);
    mencoderDir = QFileInfo(sss).absoluteFilePath();
    qDebug() << sss << mencoderDir;
   if (!QFileInfo(mencoderDir).exists() || !mencoderDir.contains("mencoder.exe")) {
       QMessageBox::warning(this, QString("openLISEM"),
       QString("Download the latest mplayer zip at http://www.mplayerhq.hu and copy mencoder.exe from this package into the lisem.exe folder! "));
   }
   else
    E_mencoderDir->setText(mencoderDir);
   qDebug() << mencoderDir;
}

void lismpeg::on_toolButton_createMP4_clicked()
{
    E_mpegProcessOutput->clear();
    screenDir = resultsDir;

    if (!screenDir.contains("screens"))
        screenDir = screenDir+"screens";
    qDebug() << screenDir;
    QString filePattern = "*.png";
    QDir directory(screenDir);
    QStringList files = directory.entryList(QStringList(filePattern), QDir::Files);
    QString listName = screenDir+"/list.txt";

    if (files.count() == 0) {
        E_mpegProcessOutput->append("No screenshots found.");
        return;
    }


    qDebug() << listName;

    QFile outputFile(listName); // Replace this with the path for the output file
    if (!outputFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        E_mpegProcessOutput->append("No screenshots found.");
        return;
    }


    QString vidname;
    QTextStream outStream(&outputFile);
    int count = 0;
    for (const QString &file : files) {
        if(!file.contains("_")) {
            outStream << screenDir+"/"+file << "\n";
            vidname = QFileInfo(file).fileName();
            count++;
        }
    }
    if (count == 0) {
        E_mpegProcessOutput->append("No map screenshots found.");
        return;
    }
    E_mpegProcessOutput->append(QString("Start: %1 map screenshots found.").arg(count));

    vidname.remove(vidname.indexOf("-"),20).append(".mp4");

    outputFile.close();
    QString prog = mencoderDir;
    int wid = 1920;
    int hei = 1080;
    if (comboBox_mpegResolution->currentIndex() == 0) {
        wid = 2560;
        hei = 1440;
    } else
        if (comboBox_mpegResolution->currentIndex() == 2) {
            wid = 1280;
            hei = 720;
        }

    QString S = QString("mf://@%1 -mf w=%2:h=%3:fps=%4:type=png -ovc x264 -x264encopts crf=20 -oac copy -of lavf -lavfopts format=mp4 -o").arg(listName).arg(wid).arg(hei).arg(doubleSpin_mpegFps->value());

    QStringList args;
    args << S.split(" ");
    args << screenDir+"/" + vidname;
    qDebug() << prog << args;
    E_mpegProcessOutput->append("Creating: "+screenDir+"/"+vidname+"\n");

    mpegProcess = new QProcess(this);
    //mpegProcess->setReadChannel ( QProcess::StandardError );
    // pcrcalc outputs on the error channel

    connect(mpegProcess, SIGNAL(readyReadStandardOutput()),this, SLOT(readFromStderr()) );
    connect(mpegProcess, SIGNAL(finished(int)),this, SLOT(finishedModel(int)) );

    mpegProcess->start(prog, args);

  //  mpegProcess->waitForFinished(-1);
  //  qDebug() << mpegProcess->readAllStandardOutput();
  //  qDebug() << mpegProcess->readAllStandardError();
    //while (mpegProcess.readyReadStandardOutput());


}

void lismpeg::readFromStderr()
{
    // Read the output of the process
    QByteArray data = mpegProcess->readAllStandardOutput();

    // Convert the data to QString
    QString output = QString::fromUtf8(data);

    // Update the current line in the QTextEdit
//    E_mpegProcessOutput->setPlainText(E_mpegProcessOutput->toPlainText() + output);

   QStringList lines = output.split("\r\n");

   // Iterate through the lines
   for (const QString& line : lines) {
       // Check if the line contains "Pos:"
       if (line.contains("Pos:")) {
           // Update the current line in the QTextEdit
           QTextCursor cursor = E_mpegProcessOutput->textCursor();
           cursor.movePosition(QTextCursor::End);
           cursor.movePosition(QTextCursor::PreviousBlock, QTextCursor::KeepAnchor);
           cursor.removeSelectedText();
           cursor.insertText(line);
       } else {
           // Append the output to the QTextEdit
           if(!line.isEmpty())
           E_mpegProcessOutput->append(line);
       }
   }

}
//---------------------------------------------------------------

void lismpeg::finishedModel(int)
{
     E_mpegProcessOutput->append("\nDone.");

}

void lismpeg::on_toolButton_stpMP4_clicked()
{
    E_mpegProcessOutput->append("\nUser interrupt, invalid MP4.");
    mpegProcess->kill();
}

