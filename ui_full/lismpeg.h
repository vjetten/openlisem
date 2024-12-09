/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#ifndef LISMPEG_H
#define LISMPEG_H

#include <QDialog>
#include <QProcess>
#include "ui_lismpeg.h"
//#include "lisemqt.h"

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

