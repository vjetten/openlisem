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
 \file LisUIModel.cpp
 \brief manage the model runs

  This file maintains the link between the interface and the LISEM model
  with data and graph displayed on the screen
  it creates the model "world" W and activates and kills the model thread
 */


#include "lisemqt.h"
#include "model.h"
#include "global.h"

//---------------------------------------------------------------------------
/** Run the model:
Save the current inyerface as a temporari run file, read by the model
Make the model world and run it
*/
void lisemqt::runmodel()
{
    //NOTE op.runfilename is set in function openRunFile()
    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this,"openLISEM",QString("Load a runfile first!"));
        return;
    }

    label_runfilename->setText(op.runfilename);
    /* TODO if run from commandline this name must exist */

    savefile(QString(op.LisemDir+"openlisemtmp.run"));
    // save the current settings as a runfile that is read by the model
    // in savefile(string) the runfile is updated with all user options and map names

    tabWidget->setCurrentIndex(2);
    //switch to output screen

    initPlot();

    InitOP();
    // wipe the result screen

    //=======================================================================================//
    W = new TWorld();
    // make the model world !!!

    connect(W, SIGNAL(show(void)),this, SLOT(worldShow(void)),Qt::QueuedConnection);
    connect(W, SIGNAL(done(QString)),this, SLOT(worldDone(QString)),Qt::QueuedConnection);
    connect(W, SIGNAL(debug(QString)),this, SLOT(worldDebug(QString)),Qt::QueuedConnection);
    // connect emitted signals from the model thread to the interface routines that handle them

    W->stopRequested = false;
    // stoprequested is used to stop the thread with the interface
    W->waitRequested = false;
    // waitrequested is used to pause the thread with the interface, only on windows machines!
    W->start();
    // start the model thread, executes W->run()
    //=======================================================================================//

}
//---------------------------------------------------------------------------
void lisemqt::pausemodel()
{
    if(W)
    {
        W->waitRequested = !W->waitRequested;
        if (!W->waitRequested)
        {
            label_debug->setText("User continue...");
            W->condition.wakeAll();
        }
    }
}
//---------------------------------------------------------------------------
void lisemqt::stopmodel()
{
    if(W)
        W->stopRequested = true;
}
//---------------------------------------------------------------------------
void lisemqt::worldShow()
{
    // copy the run results from the "output structure op" to the ui labels
    // "op" is filled in the model run each timestep
    // "op" struct is declared in lisUIoutput.h
    // "op" struct is shared everywhere in global.h

    label_dx->setText(QString::number(op.dx,'f',3));
    label_area->setText(QString::number(op.CatchmentArea/10000,'f',3));
    label_time->setText(QString::number(op.time,'f',3));
    label_endtime->setText(QString::number(op.EndTime,'f',3));
    label_runtime->setText(QString::number(op.t,'f',3));
    label_endruntime->setText(QString::number(op.maxtime,'f',3));

    label_MB->setText(QString::number(op.MB,'e',3));
    label_raintot->setText(QString::number(op.RainTotmm,'f',3));
    label_watervoltot->setText(QString::number(op.WaterVolTotmm,'f',3));
    label_qtot->setText(QString::number(op.Qtotmm,'f',3));
    label_infiltot->setText(QString::number(op.InfilTotmm,'f',3));
    label_surfstor->setText(QString::number(op.SurfStormm,'f',3));
    label_interctot->setText(QString::number(op.IntercTotmm,'f',3));
    label_qtotm3->setText(QString::number(op.Qtot,'f',3));
    label_qtotm3sub->setText(QString::number(op.QtotPlot,'f',3));
    label_qpeak->setText(QString::number(op.Qpeak,'f',3));
    label_qpeaktime->setText(QString::number(op.QpeakTime,'f',3));
    label_ppeaktime->setText(QString::number(op.RainpeakTime,'f',3));
    label_discharge->setText(QString::number(op.Q,'f',3));
    label_QPfrac->setText(QString::number((op.RainTotmm > 0 ? op.Qtotmm/op.RainTotmm*100 : 0),'f',3));
    if (checkBuffers->isChecked())
        label_buffervol->setText(QString::number(op.BufferVolTot,'f',3));

    if (!checkNoErosion->isChecked())
    {
        int dig = 2;
        label_MBs->setText(QString::number(op.MBs,'e',dig));
        label_splashdet->setText(QString::number(op.DetTotSplash,'f',dig));
        label_flowdet->setText(QString::number(op.DetTotFlow,'f',dig));
        label_sedvol->setText(QString::number(op.SedTot,'f',dig));
        label_dep->setText(QString::number(op.DepTot,'f',dig));

        label_detch->setText(QString::number(op.ChannelDetTot,'f',dig));
        label_depch->setText(QString::number(op.ChannelDepTot,'f',dig));
        label_sedvolch->setText(QString::number(op.ChannelSedTot,'f',dig));

        label_soilloss->setText(QString::number(op.SoilLossTot,'f',dig));
        label_soillosskgha->setText(QString::number(op.SoilLossTot/(op.CatchmentArea/10000)*1000,'f',dig));

        double SDR = op.DetTotSplash + op.ChannelDetTot + op.DetTotFlow;
        SDR = (SDR > 0? 100*op.SoilLossTot/(SDR) : 0);
        SDR = min(SDR ,100);
        label_SDR->setText(QString::number(SDR,'f',dig));
        if (checkBuffers->isChecked() || checkSedtrap->isChecked())
            label_buffersed->setText(QString::number(op.BufferSedTot,'f',dig));
    }

    progressBar->setMaximum(op.maxstep);
    progressBar->setValue(op.runstep);

    showPlot();

    ShowMap();

    // 6 line text output below hydrographs
    //textGraph->setMaximumBlockCount(7);

    if (checkNoErosion->isChecked())
    {
        if(!checkIncludeTiledrains->isChecked())
            textGraph->appendPlainText(QString("%1 %2 %3").arg(op.time,15,'f',3,' ').arg(op.P,15,'f',3,' ').arg(op.Q,15,'f',3,' '));
        else
            textGraph->appendPlainText(QString("%1 %2 %3 %4").arg(op.time,15,'f',3,' ').arg(op.P,15,'f',3,' ').arg(op.Q,15,'f',3,' ').arg(op.Qtile,15,'f',3,' '));
    }
    else
    {
        if(!checkIncludeTiledrains->isChecked())
            textGraph->appendPlainText(QString("%1 %2 %3 %4 %5").arg(op.time,15,'f',3,' ').arg(op.P,15,'f',3,' ').arg(op.Q,15,'f',3,' ').arg(op.Qs,12,'f',3).arg(op.C,15,'f',3,' '));
        else
            textGraph->appendPlainText(QString("%1 %2 %3 %4 %5 %6").arg(op.time,15,'f',3,' ').arg(op.P,15,'f',3,' ').arg(op.Q,15,'f',3,' ').arg(op.Qs,12,'f',3).arg(op.C,15,'f',3,' ').arg(op.Qtile,15,'f',3,' '));
    }

}
//---------------------------------------------------------------------------
void lisemqt::ShowMap()
{
    //MapPlot->replot();
    //if all calcs are done in the model only replot here?
    //op needs a poiinter to TMMap
}
//---------------------------------------------------------------------------
void lisemqt::worldDone(const QString &results)
{
    label_debug->setText(results);
    if (results.contains("ERROR"))
        QMessageBox::critical(this,QString("openLISEM"), results, QMessageBox::Ok );

    shootScreen();

    // arrive here after model emits done signal
    if (W)
    {
        delete W;
        W=NULL;
    }
    //free the world instance

    killPlot();
    // free the plot discharge bdata

    QFile(QString(op.LisemDir+"openlisemtmp.run")).remove();
    // delete the temp run file
}
//---------------------------------------------------------------------------
// this function is linked to the debug signal emitted from the model world
void lisemqt::worldDebug(const QString &results)
{
    label_debug->setText(results);
    // show messages from the World model on the screen
}
//---------------------------------------------------------------------------
void lisemqt::InitOP()
{
    op.runstep = 0;
    op.printstep = 0;
    op.maxstep = 0;
    op.CatchmentArea = 0;
    op.dx = 0;
    op.t = 0;
    op.time = 0;
    op.maxtime = 0;
    op.EndTime = 0;
    op.BeginTime = 0;
    op.MB = 0;
    op.Qtot = 0;
    op.Qtotmm = 0;
    op.Qpeak = 0;
    op.IntercTotmm = 0;
    op.WaterVolTotmm = 0;
    op.InfilTotmm = 0;
    op.RainTotmm = 0;
    op.SurfStormm = 0;
    op.InfilKWTotmm = 0;
    op.MBs = 0;
    op.DetTot = 0;
    op.DetTotSplash = 0;
    op.DetTotFlow = 0;
    op.DepTot = 0;
    op.SoilLossTot = 0;
    op.SedTot = 0;
    op.ChannelVolTot = 0;
    op.ChannelSedTot = 0;
    op.ChannelDepTot = 0;
    op.ChannelDetTot = 0;
    op.RunoffFraction = 0;
    op.RainpeakTime = 0;
    op.QpeakTime = 0;
    op.Q = 0;
    op.Qs = 0;
    op.C = 0;
    op.P = 0;
    op.BufferVolTot = 0;
    op.BufferSedTot = 0;
    label_dx->setText(QString::number(op.dx,'f',3));
    label_area->setText(QString::number(op.CatchmentArea/10000,'f',3));
    label_time->setText(QString::number(op.time,'f',3));
    label_endtime->setText(QString::number(op.EndTime,'f',3));
    label_runtime->setText(QString::number(op.t,'f',3));
    label_endruntime->setText(QString::number(op.maxtime,'f',3));

    label_MB->setText(QString::number(op.MB,'e',3));
    label_raintot->setText(QString::number(op.RainTotmm,'f',3));
    label_watervoltot->setText(QString::number(op.WaterVolTotmm,'f',3));
    label_qtot->setText(QString::number(op.Qtotmm,'f',3));
    label_infiltot->setText(QString::number(op.InfilTotmm,'f',3));
    label_surfstor->setText(QString::number(op.SurfStormm,'f',3));
    label_interctot->setText(QString::number(op.IntercTotmm,'f',3));
    label_qtotm3->setText(QString::number(op.Qtot,'f',3));
    label_qtotm3->setText(QString::number(op.QtotPlot,'f',3));
    label_qpeak->setText(QString::number(op.Qpeak,'f',3));
    label_qpeaktime->setText(QString::number(op.QpeakTime,'f',3));
    label_ppeaktime->setText(QString::number(op.RainpeakTime,'f',3));
    label_discharge->setText(QString::number(op.Q,'f',3));
    label_QPfrac->setText(QString::number((op.RainTotmm > 0 ? op.Qtotmm/op.RainTotmm*100 : 0),'f',3));
    if (checkBuffers->isChecked())
        label_buffervol->setText(QString::number(op.BufferVolTot,'f',3));

    if (!checkNoErosion->isChecked())
    {
        int dig = 2;
        label_MBs->setText(QString::number(op.MBs,'e',dig));
        label_splashdet->setText(QString::number(op.DetTotSplash,'f',dig));
        label_flowdet->setText(QString::number(op.DetTotFlow,'f',dig));
        label_sedvol->setText(QString::number(op.SedTot,'f',dig));
        label_dep->setText(QString::number(op.DepTot,'f',dig));

        label_detch->setText(QString::number(op.ChannelDetTot,'f',dig));
        label_depch->setText(QString::number(op.ChannelDepTot,'f',dig));
        label_sedvolch->setText(QString::number(op.ChannelSedTot,'f',dig));

        label_soilloss->setText(QString::number(op.SoilLossTot,'f',dig));
        label_soillosskgha->setText(QString::number(0,'f',dig));
        label_SDR->setText(QString::number(0,'f',dig));
        if (checkBuffers->isChecked() || checkSedtrap->isChecked())
            label_buffersed->setText(QString::number(op.BufferSedTot,'f',dig));
    }

    textGraph->clear();

}
//---------------------------------------------------------------------------
