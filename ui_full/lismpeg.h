#ifndef LISMPEG_H
#define LISMPEG_H

#include <QDialog>
#include "ui_lismpeg.h"


class lismpeg : public QDialog, private Ui::lismpeg
{
    Q_OBJECT

public:
    explicit lismpeg(QWidget *parent = nullptr);
    ~lismpeg();

    QString mencoderDir;
    QString screenDir;
    QString resultsDir;
    QString baseDir;
    QString vidname;

    QProcess *mpegProcess;
    QProcess *mp4Process;

    void setMencoderDir(QString d);
    void setWorkDir(QString d);
    QString getFileorDir(QString inputdir,QString title, QStringList filters, int doFile);

private slots:
    void on_toolButton_resultDir_clicked();
    void on_toolButton_mencoderDir_clicked();

    void on_toolButton_createMP4_clicked();

    void readFromStderr();
    void finishedModel(int);

    void on_toolButton_stpMP4_clicked();
    void on_toolButton_showMP4_clicked();
};

#endif // LISMPEG_H

