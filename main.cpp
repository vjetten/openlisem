/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

#include <QtGui>
#include <QApplication>

#include "ifacebasic.h"
#include "global.h"

QString ErrorString;    // declare here, referenced by error.h

//ifacebasic *iface;
// glocal declaration of window form

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
//    iface = new ifacebasic();

    op.LisemDir = QCoreApplication::applicationDirPath() + QDir::separator();
    // exe path, used for ini file

    QApplication::setStyle("WindowsVista");
    //QApplication::setStyle("Cleanlooks");

    ifacebasic iface;
    iface.show();

    return a.exec();
}
