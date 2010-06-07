/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

#include <QtGui>
#include <QApplication>

#ifdef __BASIC__
#include "ifacebasic.h"
#endif
//TODO do this better elseif?
#include "lisemqt.h"
#include "global.h"
#include "version.h"

QString ErrorString;    // declare here, referenced by error.h


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    op.LisemDir = QCoreApplication::applicationDirPath() + QDir::separator();
    // exe path, used for ini file

    QApplication::setStyle("WindowsVista");
 //   QApplication::setStyle("Cleanlooks");

#ifdef __BASIC__
    ifacebasic iface;
#endif
    lisemqt iface;

    iface.setWindowTitle(QString("openLISEM ") + VERSION + DATE);

    iface.show();
    // make an instance of the interface and show it

    return a.exec();
}
