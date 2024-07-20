/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
 \file LisUIPlot.cpp
 \brief plot the hydrographs using the Qwt library

  This file contains all the functions to plot the hydrographs
  it uses the op data for each next step


 */

#include <algorithm>
#include "lisemqt.h"
#include "global.h"

#define DIGITS 3

// defines for dual axes, branch of qwt plot
#define axisYL1 QwtAxisId(QwtAxis::YLeft,0)
#define axisYL2 QwtAxisId(QwtAxis::YLeft,1)
#define axisYR1 QwtAxisId(QwtAxis::YRight,0)
#define axisYR2 QwtAxisId(QwtAxis::YRight,1)
#define axisXB  QwtAxis::XBottom

//---------------------------------------------------------------------------
/// set up discharge plot, graphics part (not data)
/// this is called at the start of lisem, initplot() is called at the start of a model run
void lisemqt::setupPlot()
{
    QColor col;
    QwtText title;
    title.setText("Hydrograph outlet");
    title.setFont(QFont("MS Shell Dlg 2",12));
    HPlot = new QwtPlot(title);
    layout_Plot->insertWidget(0, HPlot, 1);

    // panning with the left mouse button
    (void) new QwtPlotPanner( HPlot->canvas() );

    // zoom in/out with the wheel
    (void) new QwtPlotMagnifier( HPlot->canvas() );

    PGraph = new QwtPlotCurve("Rainfall intensity");
    QGraph = new QwtPlotCurve("Discharge");
    QbGraph = new QwtPlotCurve("Water Height boundary");
    QsGraph = new QwtPlotCurve("Sediment discharge");
    CGraph = new QwtPlotCurve("Concentration");
    QtileGraph = new QwtPlotCurve("Tile drain");

    PGraph->attach(HPlot);
    QGraph->attach(HPlot);

    // order of attaching determines order of display in Legend
    PGraph->setStyle(QwtPlotCurve::Steps);
    QGraph->setStyle(QwtPlotCurve::Lines);

    QPen pen1, pen2, pen3, pen4, pen5;
    pen1.setWidth(2);
    col.setRgb( 0,30,200,255 );
    pen1.setColor(col);
    pen1.setCosmetic(true);

    pen2.setWidth(2);
    pen2.setColor("#808080");
    pen2.setCosmetic(false);

    col.setRgb( 0,180,255,255 );
    pen3.setWidth(2);
    pen3.setColor(col);
    pen3.setCosmetic(true);
    //pen3.setStyle(Qt::DashLine);


    col.setRgb( 220,0,0,255 );
    pen4.setWidth(2);
    pen4.setColor(col);//Qt::red);
    pen4.setCosmetic(false);

    //col.setRgb( 200,0,0,255 ); // darkred
    col.setRgb( 255,155,0,255 ); // darkred
    pen5.setWidth(2);
    pen5.setColor(col);
    pen5.setCosmetic(false);

    //https://forum.qt.io/topic/47823/setting-axis-label-in-qwt-multiaxis-branches-solved/5

    HPlot->setAxesCount(QwtAxis::YLeft, 1);
    HPlot->setAxesCount(QwtAxis::YRight, 1);

    HPlot->setAxisVisible( QwtAxis::YLeft );
    HPlot->setAxisVisible( QwtAxis::YRight );
    HPlot->setAxisVisible( QwtAxis::XBottom );

    QGraph->setPen(pen1);
    QGraph->setAxes(axisXB, QwtAxis::YLeft);

    PGraph->setPen(pen2);
    PGraph->setAxes(axisXB, QwtAxis::YRight);

    QbGraph->setPen(pen3);
    QbGraph->setAxes(axisXB, QwtAxis::YLeft);
    QbGraph->setStyle(QwtPlotCurve::Lines);

    QtileGraph->setPen(pen3);
    QtileGraph->setAxes(axisXB, QwtAxis::YLeft);
    QtileGraph->setStyle(QwtPlotCurve::Lines);

    QsGraph->setPen(pen4);
    QsGraph->setAxes(axisXB, QwtAxis::YRight);

    CGraph->setPen(pen5);
    CGraph->setAxes(axisXB, QwtAxis::YRight);
    QsGraph->setStyle(QwtPlotCurve::Lines);
    CGraph->setStyle(QwtPlotCurve::Lines);

    if (checkEventBased->isChecked())
        HPlot->setAxisTitle(axisXB, "time (min)");
    else
        HPlot->setAxisTitle(axisXB, "time (day)");

    if (checkUnits_ls->isChecked())
        HPlot->setAxisTitle(axisYL1, "Q (l/s)");
    else
        HPlot->setAxisTitle(axisYL1, "Q (m3/s)");
    HPlot->setAxisTitle(axisYR1, "P (mm/h)");

    HPlot->setAxisAutoScale(axisYL1,true);
    HPlot->setAxisAutoScale(axisYL2,true);
    HPlot->setAxisAutoScale(axisYR1,true);
    HPlot->setAxisAutoScale(axisYR2,true);

    HPlot->setAxisScale(axisYL1, 0.0, 1.0 );
    HPlot->setAxisScale(axisYL2, 0.0, 1.0 );
    HPlot->setAxisScale(axisYR1, 0.0, 1.0 );
    HPlot->setAxisScale(axisYR2, 0.0, 1.0 );

    QwtPlotGrid *grid = new QwtPlotGrid();
    col.setRgb( 180,180,180,180 );
    grid->setMajorPen(QPen(col, 0, Qt::DashLine));
    col.setRgb( 210,210,210,180 );
    grid->setMinorPen(QPen(col, 0 , Qt::DotLine));
    grid->attach(HPlot);


    HPlot->replot();
    // draw empty plot
}
//---------------------------------------------------------------------------
void lisemqt::onOutletChanged(int point)
{
    if(!startplot)
    {
        int index= 0;
        int oldindex = OutletIndices.indexOf(outletpoint);

        if(outletpoint == point)
        {
            return;
        }

        if(oldindex == -1)
        {
            outletpoint = 0;
            spinBoxPointtoShow->setValue(1);

            //SetTextHydrographs();
            outletgroup->setTitle(QString("Catchment outlet(s)"));

            showPlot(); // show main plot for point X

            HPlot->setTitle(QString("Catchment outlet(s)"));

        } else {
            if(point > outletpoint)
            {
                index = oldindex + 1;
            } else {
                index = oldindex - 1;
            }

            outletpoint = OutletIndices.at(index);
            spinBoxPointtoShow->setValue(outletpoint);

            //SetTextHydrographs();
            if (outletpoint > 0)
                outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
            else
                outletgroup->setTitle(QString("Catchment outflow"));


            showPlot(); // show main plot for point X

            if (outletpoint > 0)
            {
                HPlot->setTitle(QString("Hydrograph outlet %1").arg(outletpoint));
            }
            else
            {
                HPlot->setTitle(QString("Hydrograph all outflow"));
            }
        }
    }
}

//---------------------------------------------------------------------------
/// initialize graph before plotting at the start of a run
void lisemqt::initPlot()
{
    HPlot->setTitle("Hydrograph Outlet");

   // QbGraph->detach();
    QsGraph->detach();
    CGraph->detach();
    QtileGraph->detach();

    if (checkEventBased->isChecked())
        HPlot->setAxisTitle(axisXB, "time (min)");
    else
        HPlot->setAxisTitle(axisXB, "time (day)");

//    if (checkChannelBaseflow->isChecked()) {
//       QbGraph->attach(HPlot);
//    }

//    if(checkIncludeTiledrains->isChecked()) {
//        QtileGraph->attach(HPlot);
//    }


    if(checkDoErosion->isChecked())
    {
        HPlot->setAxesCount(QwtAxis::YLeft, 2);
        HPlot->setAxesCount(QwtAxis::YRight, 2);

        QsGraph->attach(HPlot);
        CGraph->attach(HPlot);

         QGraph->setAxes(axisXB, axisYL1);
         PGraph->setAxes(axisXB, axisYL2);
        QsGraph->setAxes(axisXB, axisYR1);
         CGraph->setAxes(axisXB, axisYR2);

        if (checkUnits_ls->isChecked())
            HPlot->setAxisTitle(axisYL1, "Q (l/s)");
        else
            HPlot->setAxisTitle(axisYL1, "Q (m3/s)");

        HPlot->setAxisTitle(axisYL2, "P (mm/h)");
        HPlot->setAxisTitle(axisYR1, "Qs (kg/s)");
        HPlot->setAxisTitle(axisYR2, "C (g/l)");
    }
    else
    {
        HPlot->setAxesCount(QwtAxis::YLeft, 1);
        HPlot->setAxesCount(QwtAxis::YRight, 1);

        QGraph->setAxes(axisXB, axisYL1);
        PGraph->setAxes(axisXB, axisYR1);
        QbGraph->setAxes(axisXB, axisYR1);

        if (checkUnits_ls->isChecked())
            HPlot->setAxisTitle(axisYL1, "Q (l/s)");
        else
            HPlot->setAxisTitle(axisYL1, "Q (m3/s)");

        HPlot->setAxisTitle(axisYR1, "P (mm/h)");
        if (checkWaterUserIn->isChecked())
            HPlot->setAxisTitle(axisYR1, "WH (m)");
    }

    // redraw legend with nr of variables
    QwtLegend *legend = new QwtLegend(HPlot);
    legend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
    HPlot->insertLegend(legend, QwtPlot::BottomLegend);
    //legend

}
//---------------------------------------------------------------------------
void lisemqt::showPlot()
{
    if (!checkEventBased->isChecked())
        HPlot->setAxisScale(axisXB, op.BeginTime/1440, op.EndTime/1440);
    else
        HPlot->setAxisScale(axisXB, op.BeginTime, op.EndTime);

    int index = OutletIndices.indexOf(this->outletpoint);

    QGraph->setSamples(op.Time,*op.OutletQ[index]);
    PGraph->setSamples(op.Time,op.Pmm);

    if (checkWaterUserIn->isChecked())
        QbGraph->setSamples(op.Time,*op.Wavein[index]);

    int _j = op.OutletQ[index]->count()-1; // last value index

    // qmax[index] = std::max(qmax[index] , 1.1*op.OutletQ[index]->at(_j));

    for (int i = 0; i < OutletIndices.count(); i++)  {
        qmax[i] = std::max(qmax[i] , 1.1*op.OutletQ[i]->at(_j));

        if (checkDoErosion->isChecked()) {
            qsmax[i] = std::max(qsmax[i] , 1.1*op.OutletQs[i]->at(_j));
            cmax[i] = std::max(cmax[i] , 1.1*op.OutletC[i]->at(_j));
        }
    }

    pmax = std::max(pmax, 2.0*op.Pmm[_j]);

    if(checkDoErosion->isChecked())
    {
        //        qsmax[index] = std::max(qsmax[index] , 1.1*op.OutletQs[index]->at(_j));
        //        cmax[index] = std::max(cmax[index] , 1.1*op.OutletC[index]->at(_j));

        QsGraph->setSamples(op.Time,*op.OutletQs[index]);
        CGraph->setSamples(op.Time,*op.OutletC[index]);

        HPlot->setAxisScale( axisYL1, 0.0, qmax[index] );
        HPlot->setAxisScale( axisYL2, pmax, 0.0 );
        HPlot->setAxisScale( axisYR1, 0.0, qsmax[index] );
        HPlot->setAxisScale( axisYR2, 0.0, cmax[index] );
    } else {
        HPlot->setAxisScale(axisYL1, 0.0, qmax[index] );
        HPlot->setAxisScale(axisYR1, pmax,0.0 );
    }


    if(checkIncludeTiledrains->isChecked())
        QtileGraph->setSamples(op.Time,op.Qtile);

    HPlot->replot();
}
//---------------------------------------------------------------------------
// initializes plot with world data, called in worldshow, done once
void lisemqt::startPlots()
{
    if (!startplot)
        return;

    times.clear();

    qmax.clear();
    qsmax.clear();
    cmax.clear();
    pmax = 1.0;
    // to start the max finding
    for(int i =0; i < op.OutletIndices.length(); i++)
    {
        qmax.append(1.0);
        qsmax.append(1.0);
        cmax.append(1.0);
    }

    OutletIndices.clear();
    OutletLocationX.clear();
    OutletLocationY.clear();

    OutletIndices.append(op.OutletIndices);
    OutletLocationX.append(op.OutletLocationX);
    OutletLocationY.append(op.OutletLocationY);

    outletpoint = op.OutletIndices.at(1);//1;
    spinBoxPointtoShow->setValue(1);
    spinBoxPointtoShow->setMaximum(OutletIndices.at(OutletIndices.length()-1));
    label_hydroCount->setText(QString("Output all (0) or point (1-%1)").arg(OutletIndices.count()-1));

    if (outletpoint > 0)
    {
        outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
        HPlot->setTitle(QString("Hydrograph point %1").arg(outletpoint));
    }
    else
    {
        outletgroup->setTitle(QString("Total domain outflow"));
        HPlot->setTitle(QString("Combined hydrograph domain outflow"));
    }
}
//---------------------------------------------------------------------------
void lisemqt::showOutputData()
{
    // copy the run results from the "output structure op" to the ui labels
    // "op" is filled in the model run each timestep
    // "op" struct is declared in lisUIoutput.h
    // "op" struct is shared everywhere in global.h

    int dig = E_DigitsOut->value();//DIGITS;

    if(op.OutletQ.length() == 0)
    {
        return;
    }

    QString format;
    // if(darkLISEM)
    //     format = QString("<font color=#ffffaa>%2</font>");
    // else
    //     format= QString("<font color=#000000>%2</font>");
    format="%2";

    if (startplot) {
        label_dx->setText(format.arg(QString::number(op._dx,'f',dig)));
        label_area->setText(format.arg(QString::number(op.CatchmentArea/1000000,'f',dig)));
        int days = op.EndTime/1440;
        int mins = long(op.EndTime) % 1440;
        QString ts = QString("%1:%2").arg(days,3,10,QLatin1Char('0')).arg(mins,4,10,QLatin1Char('0'));
        label_endtime->setText(ts);//format.arg(QString::number(op.EndTime,'f',dig)));
    }

    int days = op.time/1440;
    int mins = long(op.time) % 1440;
    QString ts = QString("%1:%2").arg(days,3,10,QLatin1Char('0')).arg(mins,4,10,QLatin1Char('0'));
    label_time->setText(ts);//format.arg(QString::number(op.time,'f',dig)));

    label_runtime->setText(format.arg(QString::number(op.t,'f',dig)));
    label_endruntime->setText(format.arg(QString::number(op.maxtime,'f',dig)));

    // mm output
    label_ETatot->setText(format.arg(QString::number(op.ETaTotmm,'f',dig)));
    label_raintot->setText(format.arg(QString::number(op.RainTotmm,'f',dig)));
    label_watervoltot->setText(format.arg(QString::number(op.WaterVolTotmm,'f',dig)));
    if (checkStormDrains->isChecked()) {
        label_99->setText("Storm drains");
        label_stormdraintot->setText(format.arg(QString::number(op.StormDrainTotmm,'f',dig)));
    } else
        if (doChannelBaseflow) {
            label_99->setText("Added Baseflow");
            label_stormdraintot->setText(format.arg(QString::number(op.BaseFlowTotmm,'f',dig)));
        }

    label_qtot->setText(format.arg(QString::number(op.Qtotmm,'f',dig)));
    label_infiltot->setText(format.arg(QString::number(op.InfilTotmm,'f',dig)));
    label_surfstor->setText(format.arg(QString::number(op.SurfStormm,'f',dig)));
    label_interctot->setText(format.arg(QString::number(op.IntercTotmm+op.IntercHouseTotmm+op.IntercLitterTotmm,'f',dig)));
    if (/*checkOverlandFlow1D->isChecked()*/E_OFWaveType->currentIndex() == 0 && !checkIncludeChannel->isChecked())
        label_floodVolmm->setText(format.arg(QString::number(0,'f',dig)));
    else
        label_floodVolmm->setText(format.arg(QString::number(op.volFloodmm,'f',dig)));

    label_watervolchannel->setText(format.arg(QString::number(op.ChannelVolTotmm,'f',dig)));
    //label_baseflowtot->setText(format.arg(QString::number(op.BaseFlowtotmm,'f',dig)));
    //  label_litterstore->setText(QString::number(op.LitterStorageTotmm,'f',dig));

    // peak time
    label_QPfrac->setText(format.arg(QString::number((op.RainTotmm > 0 ? std::max(0.0,op.Qtotmm-op.BaseFlowTotmm)/op.RainTotmm*100 : 0),'f',dig)));
    label_ppeaktime->setText(format.arg(QString::number(op.RainpeakTime,'f',2)));

    // mass balance
    label_MB->setText(QString::number(op.MB,'e',dig));
    if (op.MB > 0)
        label_MB->setText(" "+label_MB->text());

    int j = OutletIndices.indexOf(this->outletpoint);

    double vv = op.OutletQpeak.at(j);
    if (checkUnits_ls->isChecked()) {
    //if (vv < 10000) {
        label_82->setText("Qpeak (l/s)");
        label_qpeaksub->setText(format.arg(QString::number(vv,'f',2)));
    } else {
       // vv = vv/1000.0;
        label_82->setText("Qpeak (m3/s)");
        label_qpeaksub->setText(format.arg(QString::number(vv,'f',3)));
    }

    label_qpeaktime->setText(format.arg(QString::number(op.OutletQpeaktime.at(j),'f',2)));
    if (op.OutletQtot.at(j) < 10000)
        label_qtotm3sub->setText(format.arg(QString::number(op.OutletQtot.at(j),'f',2)));
    else
        label_qtotm3sub->setText(format.arg(QString::number(op.OutletQtot.at(j),'e',3)));


    int len = op.OutletQ.at(j)->length()-1;

    vv = op.OutletQ.at(j)->at(len);
    if (checkUnits_ls->isChecked()) {
        label_54->setText("Q (l/s)");
        label_dischargesub->setText(format.arg(QString::number(vv,'f',2)));
    } else {
        label_54->setText("Q (m3/s)");
        label_dischargesub->setText(format.arg(QString::number(vv,'f',3)));
    }

    if (checkUnits_ls->isChecked())
        checkBox_OutRunoff->setText("ro - runoff + channel discharge (l/s)");
    else
        checkBox_OutRunoff->setText("ro - runoff + channel discharge (m3/s)");

    label_soillosssub->setEnabled(checkDoErosion->isChecked());
    label_94->setEnabled(checkDoErosion->isChecked());
    label_Qssub->setEnabled(checkDoErosion->isChecked());
    label_105->setEnabled(checkDoErosion->isChecked());

    if(checkDoErosion->isChecked())
    {
        dig = E_DigitsOut->value();//DIGITS;
        label_MBs->setText(QString::number(op.MBs,'e',dig));
        if (op.MBs > 0)
            label_MBs->setText(" "+label_MBs->text());

        label_soillosssub->setText(format.arg(QString::number(op.OutletQstot.at(j),'f',dig)));
        label_Qssub->setText(format.arg(QString::number(op.OutletQs.at(j)->at(len),'f',dig)));
        label_splashdet->setText(format.arg(QString::number(op.DetTotSplash,'f',dig)));
        label_flowdet->setText(format.arg(QString::number(op.DetTotFlow+op.FloodDetTot,'f',dig)));
        label_sedvol->setText(format.arg(QString::number(op.SedTot+op.FloodSedTot,'f',dig)));
        label_dep->setText(format.arg(QString::number(op.DepTot+op.FloodDepTot,'f',dig)));

        label_detch->setText(format.arg(QString::number(op.ChannelDetTot,'f',dig)));
        label_depch->setText(format.arg(QString::number(op.ChannelDepTot,'f',dig)));
        label_sedvolch->setText(format.arg(QString::number(op.ChannelSedTot,'f',dig)));

        label_soilloss->setText(format.arg(QString::number(op.SoilLossTot,'f',dig)));
        label_soillosskgha->setText(format.arg(QString::number(op.SoilLossTot/(op.CatchmentArea/10000)*1000,'f',dig)));

        double SDR = op.DetTotSplash + op.ChannelDetTot + op.DetTotFlow;
        SDR = (SDR > 0? 100*op.SoilLossTot/(SDR) : 0);
        SDR = std::min(SDR ,100.0);
        label_SDR->setText(format.arg(QString::number(SDR,'f',dig)));
    } else {
        QString zero = QString::number(0,'f',E_DigitsOut->value());
        label_MBs->setText(zero);

        label_soillosssub->setText(zero);
        label_Qssub->setText(zero);
        label_splashdet->setText(zero);
        label_flowdet->setText(zero);
        label_sedvol->setText(zero);
        label_dep->setText(zero);

        label_detch->setText(zero);
        label_depch->setText(zero);
        label_sedvolch->setText(zero);

        label_soilloss->setText(zero);
        label_soillosskgha->setText(zero);
        label_SDR->setText(zero);
    }
}

void lisemqt::showOutputDataZero()
{
    if(op.OutletQ.length() == 0)
    {
        return;
    }

    // QString format;
    // // if(darkLISEM)
    // //     format = QString("<font color=#ffffaa>%2</font>");
    // // else
    // //     format= QString("<font color=#000000>%2</font>");
    // format = "%2";

    QString zero = QString::number(0,'f',E_DigitsOut->value());

    label_dx->setText(zero);
    label_area->setText(zero);
    label_endtime->setText(zero);
    label_time->setText(zero);

    label_runtime->setText(zero);
    label_endruntime->setText(zero);

    // mm output
    label_ETatot->setText(zero);
    label_raintot->setText(zero);
    label_watervoltot->setText(zero);
    label_stormdraintot->setText(zero);
    label_qtot->setText(zero);
    label_infiltot->setText(zero);
    label_surfstor->setText(zero);
    label_interctot->setText(zero);
    label_floodVolmm->setText(zero);

    label_watervolchannel->setText(zero);

    // peak time
    label_QPfrac->setText(zero);
    label_ppeaktime->setText(zero);

    // mass balance
    label_MB->setText(zero);
    label_qpeaksub->setText(zero);
    label_qpeaktime->setText(zero);
    label_qtotm3sub->setText(zero);
    label_dischargesub->setText(zero);

    label_MBs->setText(zero);

    label_soillosssub->setText(zero);
    label_Qssub->setText(zero);
    label_splashdet->setText(zero);
    label_flowdet->setText(zero);
    label_sedvol->setText(zero);
    label_dep->setText(zero);

    label_detch->setText(zero);
    label_depch->setText(zero);
    label_sedvolch->setText(zero);

    label_soilloss->setText(zero);
    label_soillosskgha->setText(zero);
    label_SDR->setText(zero);

}
