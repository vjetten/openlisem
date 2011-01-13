/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

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

  //  QApplication::setStyle("WindowsVista");
 //   QApplication::setStyle("Cleanlooks");
 //   QApplication::addLibraryPath(QString("INSTALL/libs/"));
    //QSystemTrayIcon icon;

    lisemqt iface;

    iface.setWindowTitle(VERSION);
			 //QString("openLISEM ") + VERSION + DATE);

    iface.show();
    // make an instance of the interface and show it

    return a.exec();
}
