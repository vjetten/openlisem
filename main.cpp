
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
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
**  Author: Victor Jetten
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
    // https://doc.qt.io/qt-5/highdpi.html
  //   qputenv("QT_ENABLE_HIGHDPI_SCALING","1");
  //   qputenv("QT_SCALE_FACTOR_ROUNDING_POLICY","1");

     Fixture fixture; // <= necessary for GDAL
     QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
     QApplication app(argc, argv);

     app.setWindowIcon(QIcon(":/openlisem.ico"));
     QLocale loc = QLocale::system(); // current locale
     loc.setNumberOptions(QLocale::c().numberOptions()); // borrow number options from the "C" locale
     QLocale::setDefault(loc);

    app.setStyle(QStyleFactory::create("Fusion"));
    //  app.setStyle(QStyleFactory::create("Windows"));

    op.LisemDir = QCoreApplication::applicationDirPath()+"/";
    // exe path, used for ini file

    if (argc <= 1)
    {
        lisemqt iface;

        iface.setWindowTitle(VERSION);
        iface.show();
        // make an instance of the interface and show it

        return app.exec();
    }
    else
    {
        // run without interface
        bool noInterface;
        bool noOutput;
        bool batchmode;
        bool options;

        QString name;
        QString optionStr;

        QStringList args;
        for (int i = 1; i < argc+1; i++)
            args << QString(argv[i]);

        bool runfound = false;

        //get run filename
        int n = 0;
        for(QString const& str: args)
        {
            if (str.contains("-r"))
            {
                runfound = true;
                name = str;
                name.remove("-r");
                if (name.isEmpty())
                    name= args[n+1];
            }
            n++;
        }

        // get user defined individual options
        optionList.clear();
        for(QString str: args)
        {
            if (runfound)
            {
                optionStr += str+ " ";;
                if (str.contains("]"))
                {
                    str.remove(']');
                    options = false;
                }
            }
            if (str.contains("["))
            {
                str.remove('[');
                options = true;
                optionStr = str+ " ";
            }

        }

        //qDebug() << "-c" << optionStr;

        int count1 = args.indexOf("-ni");
        int count2 = args.indexOf("-no");
        int count3 = args.indexOf("-b");
        int count4 = args.indexOf("-c");

        //   qDebug() << count0 << count1 << count2;

        noInterface = count1 > -1;
        noOutput = count2 > -1;
        batchmode = count3 > -1;
        options = count4 > -1;
        if (!batchmode)
        {
            if (noInterface)
                noOutput = false;
            if (noOutput)
                noInterface = true;
        }

        if (!runfound)
        {
            printf("syntax:\nopenlisem [-b,-ni,-no] -r runfilename -c options\n"
                   "-b batch mode with interface, run multiple files\n"
                   "-ni = no inteface, with counter and info\n"
                   "-no = only error output\n"
                   "-c = seperate run file override options"
                   "options = exact strings from the run file within [ ] and separated with ; \n"
                   "          for instance: -c [Subsoil drainage=0]\n");
            return 0;
        }

        if (runfound && options)
        {
            optionList = optionStr.split(';');
        }


        if (runfound && !noInterface)
        {
            lisemqt iface(0, batchmode, name);
            //            iface.doBatchmode = batchmode;
            //            iface.batchRunname = name;
            iface.setWindowTitle(VERSION);
            iface.show();

            return app.exec();
        }
        else
        {
            op.runfilename = name;

            TWorld *W = new TWorld();
            // make the model world
            W->stopRequested = false;
            W->waitRequested = false;
            W->noInterface = noInterface;
            W->noOutput = noOutput;
            W->start();

            return app.exec();
        }
    }

}
