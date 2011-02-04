
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
//#include "version.h" // moved to lisemqt.h

QString ErrorString;    // declare here, referenced by error.h


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    op.LisemDir = QCoreApplication::applicationDirPath() + QDir::separator();
    // exe path, used for ini file

    lisemqt iface;

    iface.setWindowTitle(VERSION);
			 //QString("openLISEM ") + VERSION + DATE);

    iface.show();
    // make an instance of the interface and show it

    return a.exec();
}
