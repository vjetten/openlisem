
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

#include "lisemqt.h"
#include "global.h"

QString ErrorString;    // declare here, referenced by error.h

void Error(QString s)
{
   ErrorString = s;
   throw 1;
}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setWindowIcon(QIcon("openlisem.ico"));

    op.LisemDir = QCoreApplication::applicationDirPath()+"/";
    // exe path, used for ini file

    lisemqt iface;

     iface.setWindowTitle(VERSION);

    if (argc <= 1)
    {
         iface.show();
         // make an instance of the interface and show it
    }
    else
    {
       // run without interface

       QStringList args;
       for (int i = 1; i < argc+1; i++)
            args << QString(argv[i]);

       TWorld *W = new TWorld();
       // make the model world

       if (args.contains("-ni"))
       {
          W->noInterface = true;
          W->noOutput = false;
       }
       else
       if (args.contains("-no"))
       {
          W->noInterface = true;
          W->noOutput = true;
       }
       else
       {
          printf("syntax:\nopenlisem [-ni,-no] runfilename\n-ni = no inteface but counter and info\n-no = only error output\n");
          delete W;
          return 0;
       }

       W->stopRequested = false;
       // stoprequested is used to stop the thread with the interface
       W->waitRequested = false;
       // waitrequested is used to pause the thread with the interface, only on windows machines!
       int l = args.count()-2;

       W->temprunname = args[l];

       W->start();
    }

    return app.exec();
}
