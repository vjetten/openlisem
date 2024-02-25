#include "lisemqt.h"
//#include "lismpeg.h"
//#include "ui_lismpeg.h"


//---------------------------------------------------------------------------
lismpeg::lismpeg(QWidget *parent) :
        QDialog(parent)
{
    setupUi(this);
    this->setWindowTitle("Create mpeg from screenshots");
   // this->setWindowFlags(Qt::WindowSystemMenuHint | Qt::WindowTitleHint);
}
//---------------------------------------------------------------------------
lismpeg::~lismpeg()
{
  //  delete ui;
}
//---------------------------------------------------------------------------
