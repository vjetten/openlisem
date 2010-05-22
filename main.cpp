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
    op.LisemDir = QCoreApplication::applicationDirPath() + QDir::separator();

  //       QApplication::setStyle("Cleanlooks");
    QApplication::setStyle("WindowsVista");
     //a.setStyleSheet("QLabel# { background-color: yellow }");
    //a.addLibraryPath("./libs/debug");

    ifacebasic iface;
    iface.show();


    return a.exec();
}
