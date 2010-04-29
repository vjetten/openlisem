#ifndef ifacebasic_H
#define ifacebasic_H

#include <QtGui/QWidget>
#include <QtGui>
#include <QFileDialog>

#include "ui_ifacebasic.h"
#include "model.h"

struct output{
    int runstep;
    int maxstep;
    double t;
    double time, maxtime;

    double MB, Qtot, Qtotmm, Qpeak, IntercTot, WaterVolTot, InfilTot, RainTot, SurfStorTot, InfilKWTot;
    double MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedVolTot;
    double ChannelVolTot, ChannelSedTot, ChannelDepTot;

    bool SwitchErosion;
    bool SwitchIncludeChannel;
    QString runfilename;
    QString LisemDir;
};


class ifacebasic : public QWidget, public Ui::ifacebasicClass
{
    Q_OBJECT

public:
    ifacebasic(QWidget *parent = 0);
    ~ifacebasic();

private slots:
    void on_runButton_clicked();
    void on_toolButton_ShowRunfile_clicked();
    void Showit(const int step);
    void worldDone(const QString &results);
    void worldDebug(const QString &results);
   // void on_checkChannel_clicked();
   // void on_checkErosion_clicked();
    void on_toolButton_runfilename_clicked();
    void StorePath();
    void GetStorePath();


private:
    TWorld *W;
};

#endif // ifacebasic_H
