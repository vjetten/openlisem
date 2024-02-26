#include "lisemqt.h"

//---------------------------------------------------------------------------
lismpeg::lismpeg(QWidget *parent) :
        QDialog(parent)
{
    setupUi(this);
    this->setWindowTitle("Create mpeg from screenshots");
   // this->setWindowFlags(Qt::WindowSystemMenuHint | Qt::WindowTitleHint);

    mencoderDir = QCoreApplication::applicationDirPath() + "/mencoder.exe";
    E_mencoderDir->setText(mencoderDir);
    resize(qApp->primaryScreen()->size().height()*2/3,qApp->primaryScreen()->size().height()*2/3);

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
    QString sss = getFileorDir(mencoderDir,"Select dir with mencoder.exe", filters, 2);

    mencoderDir = QFileInfo(sss).absoluteFilePath();
//    if (!QFileInfo(mencoderDir).exists() || !mencoderDir.contains("mencoder.exe")) {
//        qWarning(" download the latest mplayer zip at http://www.mplayerhq.hu and copy mencoder.exe from this package into the lisem.exe folder! ");
//    }
//    else
        E_mencoderDir->setText(mencoderDir);

}

void lismpeg::on_toolButton_createMP4_clicked()
{
    screenDir = resultsDir+"screens";
    QString filePattern = "*.png";
    QDir directory(screenDir);
    QStringList files = directory.entryList(QStringList(filePattern), QDir::Files);
    QString listName = screenDir+"/list.txt";

    qDebug() << listName;

    QFile outputFile(listName); // Replace this with the path for the output file
    if (!outputFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qWarning() << "Failed to open output file list.txt";
        return;
    }
    QString vidname;
    QTextStream outStream(&outputFile);
    for (const QString &file : files) {
        if(!file.contains("_")) {
            outStream << screenDir+"/"+file << "\n";
            vidname = QFileInfo(file).fileName();
        }
    }
    vidname.remove(vidname.indexOf("-"),20).append(".mp4");

    outputFile.close();
    QString prog = mencoderDir;
    QString S = QString("mf://@%1 -mf w=1920:h=1080:fps=12:type=png -ovc x264 -x264encopts crf=20:threads=6 -oac copy -o").arg(listName);
    QStringList args;
    args << S.split(" ");
    args << screenDir+"/" + vidname;
qDebug() << prog << args;
    QProcess mpegProcess;
    mpegProcess.start(prog, args);

    mpegProcess.waitForFinished(-1);
    qDebug() << mpegProcess.readAllStandardOutput();
    qDebug() << mpegProcess.readAllStandardError();
    //while (mpegProcess.readyReadStandardOutput());
E_mpegProcessOutput->insertPlainText(mpegProcess.readAllStandardError());
    qDebug() << "done";
}

