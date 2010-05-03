#include <QtGui>
#include <QApplication>

#include "ifacebasic.h"
#include "global.h"

QString ErrorString;
//ifacebasic *iface;
// glocal declaration of window form

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
//    iface = new ifacebasic();
    QFileInfo info(argv[0]);
    op.LisemDir = info.path() + "/";

       QApplication::setStyle("Cleanlooks");
   //  QApplication::setStyle("WindowsVista");

    ifacebasic iface;
    iface.show();


    return a.exec();
}
