
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
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file main.cpp
  \brief main function, making and calling interface

functions: \n
- int main(int argc, char *argv[]) \n
 */

#include <QtGui>
#include <QApplication>

#include "fixture.h"
#include "lisemqt.h"
#include "global.h"

QStringList optionList;

int main(int argc, char *argv[])
{
    // open console but only if run from cmd.exe in win
#ifdef _WIN32
    if (AttachConsole(ATTACH_PARENT_PROCESS)) {
        freopen("CONOUT$", "w", stdout);
        freopen("CONOUT$", "w", stderr);
    }
#endif

    Fixture fixture; // <= necessary for GDAL
 //   QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);

    app.setWindowIcon(QIcon(":/openlisem.ico"));

    app.setStyle(QStyleFactory::create("Fusion"));
   //  app.setStyle(QStyleFactory::create("Windows"));

    //    QFont font = qApp->font();
    //    font.setPixelSize(11);
    //    qApp->setFont(font);
       // increase font size for better reading
//       QFont defaultFont = QApplication::font();
//       defaultFont.setPointSize(defaultFont.pointSize()+2);
//       qApp->setFont(defaultFont);
       // modify palette to dark

    //qputenv("QT_SCALE_FACTOR", "1.0");

    qputenv("QT_AUTO_SCREEN_SCALE_FACTOR","1");

    op.LisemDir = QCoreApplication::applicationDirPath()+"/";
    // exe path, used for ini file

    QStringList args=QCoreApplication::arguments();

    if (argc <= 1)
    {
        lisemqt iface;

        iface.setWindowTitle(VERSION);
        iface.show();

        return app.exec();

    } else {

        QString ag = args.join(" ");
        // run without interface nin console
        QString name;
        if (ag.contains("?")) {
            printf("syntax:\nlisem -r runfile \n");
            return 0;
        }

        if (ag.contains("-r")) {
            QStringList sl = ag.split("-r");
            name = sl[1].simplified();
            qDebug() << "running: " << name;

            lisemqt iface(0, true, name);
            iface.setWindowTitle(VERSION);
            iface.show();

            return app.exec();

        } else {
            printf("syntax:\nlisem -r runfile \n");
            return 0;
        }
    }
}
