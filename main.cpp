
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
  \brief main function, making and calling interface

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
    // QApplication* temp = new QApplication(argc, argv);
    // double width = QApplication::desktop()->width();
    // double height = QApplication::desktop()->height();

    // if (height < 1080) {
    //     // assumes that the default desktop resolution is 720p (scale of 1)
    //  //   int minWidth = 1280;


    //  //   double scale = width / minWidth;
    //  //   std::string scaleAsString = std::to_string(scale);
    //   //  QByteArray scaleAsQByteArray(scaleAsString.c_str(), scaleAsString.length());
    //    // qputenv("QT_SCALE_FACTOR", scaleAsQByteArray);
    //   //  qputenv("QT_AUTO_SCREEN_SCALE_FACTOR",scaleAsQByteArray);
    //    // qDebug() <<scaleAsQByteArray << height << width;
    // }
    //     delete temp;

    // open console but only if run from cmd.exe in win
// #ifdef _WIN32
//     if (AttachConsole(ATTACH_PARENT_PROCESS)) {
//         freopen("CONOUT$", "w", stdout);
//         freopen("CONOUT$", "w", stderr);
//     }
// #endif

    Fixture fixture; // <= necessary for GDAL
 //   QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
  //  argc += 2;

  //  argv[] = {(char*)"Appname", (char*)"--platform", (char*)"windows:dpiawareness=0"};


    QApplication app(argc, argv);

    app.setWindowIcon(QIcon(":/openlisemN.ico"));

    app.setStyle(QStyleFactory::create("Fusion"));
   //  app.setStyle(QStyleFactory::create("Windows"));

     // modify palette to dark

    //qputenv("QT_SCALE_FACTOR", "1.0");

  //  qputenv("QT_AUTO_SCREEN_SCALE_FACTOR","1");

    op.LisemDir = QCoreApplication::applicationDirPath()+"/";
    // exe path, used for ini file

    QStringList args=QCoreApplication::arguments();

    QLocale loc = QLocale::system(); // current locale
    loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
    QLocale::setDefault(loc);


    if (argc <= 1)
    {
        lisemqt iface;

        iface.setWindowTitle(VERSION);
        iface.show();

        return app.exec();

    } else {
        // 2 options:
        // noInterface = run without GUI in console
        // batchmode = run with GUI but start run automatically (default)
        bool noInterface = false;

        QString ag = args.join(" ");
        QString name;

        if (ag.contains("?")) {
            printf("syntax:\nlisem [-ni] -r runfile \n"
                   "-ni = no graphical user interface, uses runfile directly!\n");
            return 0;
        }

        // run from console with or without GUI
        if (ag.contains("-r")) {
            QStringList sl = ag.split("-r");
            name = sl[1].simplified();           


            if (ag.contains("-ni")) {
                noInterface = true;
                op.runfilename = name;
                op.doBatchmode = true;

                TWorld *W = new TWorld();
                // make the model world
				op.timeStartRun = QDateTime().currentDateTime().toString("yyMMdd-hhmm");
                W->stopRequested = false;
                W->waitRequested = false;
                W->noInterface = noInterface;
                W->start();
                qDebug() << "\nrunning OpenLISEM with:" << name;
                return app.exec();
            } else {

                //qDebug() << "running: " << name;

                lisemqt iface(0, true, name);
                iface.setWindowTitle(VERSION);
                iface.show();

            return app.exec();
            }

        } else {
            printf("syntax:\nlisem [-ni] -r runfile \n"
                   "-ni = no graphical user interface, uses runfile directly!\n");
            return 0;
        }
    }
}
