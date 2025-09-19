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

#include "fixture.h"  // for gdal
#include "lisemqt.h"
#include "global.h"
#include <iostream>

#ifdef Q_OS_WIN
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#endif

QStringList optionList;

int main(int argc, char *argv[])
{
    Fixture fixture; // <= necessary for GDAL

    QString runFileName;
    bool noInterface = false;
    bool forceRes = false;
    bool doBatch = false;
    bool syntax = true;
    double ksat1cal = -999;
    double mancal = -999;
    double manchcal = -999;
    double hcal = -999;
    QString explanation;


    if (argc == 1)
        syntax = false; // run with GUI

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        QString arg = argv[i];
        if (arg == "-ni") {
            noInterface = true;
        }
        if (arg == "-f") {
            forceRes = true;
        }
        if (arg == "-r" && i+1 < argc) {
            runFileName = argv[++i];
            syntax = false;
            doBatch = true;
        }
        if (arg == "-cks") {
            QString S = argv[++i];
            ksat1cal = S.toDouble();
        }
        if (arg == "-cn") {
            QString S = argv[++i];
            mancal = S.toDouble();
        }
        if (arg == "-cchn") {
            QString S = argv[++i];
            manchcal = S.toDouble();
        }
        if (arg == "-ch") {
            QString S = argv[++i];
            hcal = S.toDouble();
        }
        if (arg == "-S") {
            explanation = argv[++i];
        }

    }

    // this path is needed for openlisemtemp.run and openlisem.ini
    QString appDataLocalPath = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);
    QFileInfo appDataLocalFileInfo(appDataLocalPath);
    QString localPath = appDataLocalFileInfo.absolutePath() + "/lisem";
    QDir dir;
    if (!dir.exists(localPath))
        dir.mkpath(localPath);
    op.userAppDir = localPath + "/";

    op.ksat1cal = ksat1cal;
    op.mancal = mancal;
    op.chmancal = manchcal;
    op.hcal = hcal;
    op.explanation = explanation;

    if (noInterface || syntax) {
    #ifdef Q_OS_WIN
      // open a console in windows for the headless output, this works for running from a batch file and from cmd.exe
        AllocConsole();
        FILE* fp;
        freopen_s(&fp, "CONOUT$", "w", stdout);
        freopen_s(&fp, "CONOUT$", "w", stderr);
        freopen_s(&fp, "CONIN$", "r", stdin);
    #endif
        QTextStream consoleout(stdout); // text to console

        if (syntax) {
            consoleout << "syntax:\nlisem [-ni] [-f] -r runfile \n"
                       <<   "-f = create the result directory it does not exist. \n"
                       <<   "-ni = no graphical user interface, uses runfile directly.\n"
                       <<   "-r runfile = Give the full path tot he runfile.\n\n" ;
            consoleout.flush();
            #ifdef Q_OS_WIN
            system("pause"); // waits for a key press
            #endif
            return 0;
        }

        if (!runFileName.isEmpty()) {

            if (!QFileInfo(runFileName).exists()) {
                consoleout << "\nCannot find the runfile:" << runFileName << "\n";
                consoleout.flush();
                #ifdef Q_OS_WIN
                system("pause"); // waits for a key press
                #endif
                return 0;
            }

            QCoreApplication app(argc, argv); // Use QCoreApplication for headless mode

            op.runfilename = runFileName;
            op.doBatchmode = true;
            op.forceResDir = forceRes;

            //TWorld *W = new TWorld(); // pointer is not deleted so mem leak, declare directly
            TWorld W;

            // deal with different digit symbols dot or comma
            W.loc = QLocale::system(); // current locale
            W.loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
            QLocale::setDefault(W.loc);

            W.stopRequested = false;
            W.waitRequested = false;
            W.noInterface = noInterface;

            // don't use the QThread worldThread because there is no GUI, call DoModel directly
            W.DoModel();
            return 0;
            // return app.exec(); // DoModel has quit(); but that is not called properly, prevents DoModel from quiting properly
        }
    } else {
        // Use QApplication for GUI mode
        QApplication app(argc, argv);
        app.setStyle(QStyleFactory::create("Fusion"));

        // select between a standard run with GUI or a run with GUI based on a specified runfile from the command line
        if (argc <= 1) {
            lisemqt iface;
            iface.setWindowTitle(VERSION);
            iface.show();
            return app.exec();
        } else {
            if (!runFileName.isEmpty()) {
                lisemqt iface(0, doBatch, forceRes, runFileName);
                iface.setWindowTitle(VERSION);
                iface.show();
                return app.exec();
            }
        }
    }
    return 0;
}
