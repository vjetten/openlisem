#ifndef LISMPEG_H
#define LISMPEG_H

#include <QDialog>
#include "ui_lismpeg.h"
//namespace Ui {
//class lsmpeg;
//}

class lismpeg : public QDialog, private Ui::lismpeg
{
    Q_OBJECT

public:
    explicit lismpeg(QWidget *parent = nullptr);
    ~lismpeg();



//public slots:
//     bool findCondaDir();

//    private:
//        Ui::lismpeg *ui;
};

#endif // LISMPEG_H
