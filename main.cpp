
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
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

/*!
  \file main.cpp
  \brief main function, call the app based on 2 options. If in the command line '-ni'
    is found then no GUI is loaded, otherwise the interface is made and called

functions: \n
- int main(int argc, char *argv[]) \n
 */
#include <stdlib.h>
#include <QtGui>
#include <QApplication>

#include "fixture.h"
#include "lisemqt.h"
#include "global.h"

#include <iostream>

QStringList optionList;

int main(int argc, char *argv[])
{
    Fixture fixture; // <= necessary for GDAL
    QString runFileName;
    bool noInterface = false;
    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        QString arg = argv[i];
        if (arg == "-ni") {
            noInterface = true;
        } else if (arg == "-r" && i + 1 < argc) {
            runFileName = argv[++i];
        } else {
            printf("syntax:\nlisem [-ni] -r runfile \n-ni = no graphical user interface, uses runfile directly!\n");
            return 0;
        }
    }

    // Initialize application based on noInterface flag
    if (noInterface) {
        QCoreApplication app(argc, argv); // Use QCoreApplication for headless mode

        op.LisemDir = QCoreApplication::applicationDirPath() + "/";
        // exe path, used for ini file

        QString appDataLocalPath = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
        QFileInfo appDataLocalFileInfo(appDataLocalPath);
        QString localPath = appDataLocalFileInfo.absolutePath() + "/lisem";
        QDir dir;
        if (!dir.exists(localPath))
            dir.mkpath(localPath);
        op.userAppDir = localPath + "/";
        QLocale loc = QLocale::system(); // current locale
        loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
        QLocale::setDefault(loc);

        //Start the model based on the specified runfile directly
        if (!runFileName.isEmpty()) {
            op.runfilename = runFileName;
            op.doBatchmode = true;

            TWorld *W = new TWorld();

            W->stopRequested = false;
            W->waitRequested = false;
            W->noInterface = noInterface;
            W->start();
            qDebug() << "\nrunning OpenLISEM with:" << runFileName;
            return app.exec();
        } else {
            printf("syntax:\nLisem [-ni] -r runfile \n"
                   "-ni = no graphical user interface, uses runfile directly!\n");
            return 0;
        }
    } else {
        QApplication app(argc, argv); // Use QApplication for GUI mode
        app.setWindowIcon(QIcon(":/openlisemN.ico"));
        app.setStyle(QStyleFactory::create("Fusion"));

        op.LisemDir = QCoreApplication::applicationDirPath() + "/";
        // exe path, used for ini file

        QString appDataLocalPath = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
        QFileInfo appDataLocalFileInfo(appDataLocalPath);
        QString localPath = appDataLocalFileInfo.absolutePath() + "/lisem";
        QDir dir;
        if (!dir.exists(localPath))
            dir.mkpath(localPath);
        op.userAppDir = localPath + "/";
        QLocale loc = QLocale::system(); // current locale
        loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
        QLocale::setDefault(loc);

        // select between a standard run with GUI or a run with GUI based on a specified runfile from the command line
        if (argc <= 1) {
            lisemqt iface;
            iface.setWindowTitle(VERSION);
            iface.show();
            return app.exec();
        } else {
            if (!runFileName.isEmpty()) {
                lisemqt iface(0, true, runFileName);
                iface.setWindowTitle(VERSION);
                iface.show();
                return app.exec();
            } else {
                printf("syntax:\nlisem [-ni] -r runfile \n"
                       "-ni = no graphical user interface, uses runfile directly!\n");
                return 0;
            }
        }
    }
}
