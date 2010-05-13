#ifndef ifacebasic_H
#define ifacebasic_H

#include <QtGui/QWidget>
#include <QtGui>
#include <QFileDialog>

#include "ui_ifacebasic.h"
#include "model.h"

struct output{
    int runstep;
    int printstep;
    int maxstep;
    double CatchmentArea, dx, t,
      time, maxtime, EndTime;

    double MB, Qtot, Qtotmm, Qpeak, IntercTotmm, WaterVolTotmm, InfilTotmm,
      RainTotmm, SurfStorTotmm, InfilKWTotmm,
      MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedVolTot,
      ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot,
      RunoffFraction, RainpeakTime, QpeakTime;

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
    void Showit();
    void worldDone(const QString &results);
    void worldDebug(const QString &results);
   // void on_checkChannel_clicked();
   // void on_checkErosion_clicked();
    void on_toolButton_runfilename_clicked();
    void StorePath();
    void GetStorePath();
    void SetStyleUI();

private:
    TWorld *W;
};

#endif // ifacebasic_H
