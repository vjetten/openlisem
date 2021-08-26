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
**  website, information and code: http://lisem.sourceforge.net
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

    PGraph = new QwtPlotCurve("Rainfall");
    QGraph = new QwtPlotCurve("Discharge");
//    if (checkChannelBaseflow->isChecked())
        QbGraph = new QwtPlotCurve("Baseflow");

    if(checkDoErosion->isChecked()) {
        QsGraph = new QwtPlotCurve("Sediment discharge");
        CGraph = new QwtPlotCurve("Concentration");
    }
    if(checkIncludeTiledrains->isChecked())
        QtileGraph = new QwtPlotCurve("Tile drain");

    PGraph->attach(HPlot);
    QGraph->attach(HPlot);

    // order of attaching determines order of display in Legend
    PGraph->setStyle(QwtPlotCurve::Steps);
    QGraph->setStyle(QwtPlotCurve::Lines);

    QPen pen1, pen2, pen3, pen4, pen5;
    pen1.setWidth(2);
    col.setRgb( 60,60,200,255 );
    pen1.setColor(col);
    pen1.setCosmetic(true);

    pen2.setWidth(1);
    pen2.setColor("#000000");
    pen2.setCosmetic(false);

    col.setRgb( 80,160,230,255 );
    pen3.setWidth(2);
    pen3.setColor(Qt::darkCyan);
    pen3.setCosmetic(true);

    pen4.setWidth(2);
    pen4.setColor(Qt::red);
    pen4.setCosmetic(false);

    col.setRgb( 200,0,0,255 ); // darkred
    pen5.setWidth(2);
    pen5.setColor(col);
    pen5.setCosmetic(false);

    axisYL1 = new QwtAxisId(QwtPlot::yLeft,0);
    axisYL2 = new QwtAxisId(QwtPlot::yLeft,1);
    axisYR1 = new QwtAxisId(QwtPlot::yRight,0);
    axisYR2 = new QwtAxisId(QwtPlot::yRight,1);

    HPlot->setAxesCount(QwtPlot::yLeft, 2);
    HPlot->setAxesCount(QwtPlot::yRight, 2);

    QGraph->setPen(pen1);
    QGraph->setAxes(HPlot->xBottom, *axisYL1);

    PGraph->setPen(pen2);
    PGraph->setAxes(HPlot->xBottom, *axisYR1);// HPlot->yRight);

  //  if(checkChannelBaseflow->isChecked()) {
        QbGraph->setPen(pen3);
        QbGraph->setAxes(HPlot->xBottom, *axisYL1);
        QbGraph->setStyle(QwtPlotCurve::Lines);
  //  }
    if(checkIncludeTiledrains->isChecked()) {
        QtileGraph->setPen(pen3);
        QtileGraph->setAxes(HPlot->xBottom, *axisYL1);
        QtileGraph->setStyle(QwtPlotCurve::Lines);
    }
    if(checkDoErosion->isChecked()) {
        QsGraph->setPen(pen4);
        QsGraph->setAxes(HPlot->xBottom, *axisYR1);

        CGraph->setPen(pen5);
        CGraph->setAxes(HPlot->xBottom, *axisYR2);
        QsGraph->setStyle(QwtPlotCurve::Lines);
        CGraph->setStyle(QwtPlotCurve::Lines);
    }

//    PGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
//    QGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
//    QtileGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
//    QsGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
//    CGraph->setRenderHint(QwtPlotItem::RenderAntialiased);

    // set axes
    HPlot->enableAxis(HPlot->yRight,true);
    HPlot->enableAxis(HPlot->yLeft,true);
    HPlot->enableAxis(HPlot->xBottom,true);
    if (checkEventBased->isChecked())
        HPlot->setAxisTitle(HPlot->xBottom, "time (min)");
    else
        HPlot->setAxisTitle(HPlot->xBottom, "time (day)");

    HPlot->setAxisTitle(*axisYL1, "Q (l/s)");
    HPlot->setAxisTitle(*axisYR1, "P (mm/h)");

    HPlot->setAxisAutoScale(*axisYL1,true);
    HPlot->setAxisAutoScale(*axisYL2,true);
    HPlot->setAxisAutoScale(*axisYR1,true);
    HPlot->setAxisAutoScale(*axisYR2,true);

    HPlot->setAxisScale(*axisYL1, 0.0, 1.0 );
    HPlot->setAxisScale(*axisYL2, 0.0, 1.0 );
    HPlot->setAxisScale(*axisYR1, 0.0, 1.0 );
    HPlot->setAxisScale(*axisYR2, 0.0, 1.0 );

    if (darkLISEM)
        HPlot->setCanvasBackground(QBrush("#777777"));
    else
        HPlot->setCanvasBackground(QBrush(Qt::white));   // set gridlines

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

            if (doNewPlot)
                newPlot(false);
            else {
                showPlot(); // show main plot for point X

                HPlot->setTitle(QString("Catchment outlet(s)"));
                //HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
                //outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
            }
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


            if (doNewPlot)
                newPlot(true);
            else {
                showPlot(); // show main plot for point X


                if (outletpoint > 0)
                {
                    HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
                }
                else
                {
                    HPlot->setTitle(QString("Hydrograph all outflow"));
                }
            }
        }
    }
}
//---------------------------------------------------------------------------
/// initialize graph before plotting at the start of a run
void lisemqt::initPlot()
{
    HPlot->setTitle("Hydrograph Outlet");


    if (checkEventBased->isChecked())
        HPlot->setAxisTitle(HPlot->xBottom, "time (min)");
    else
        HPlot->setAxisTitle(HPlot->xBottom, "time (day)");

  //  if (checkChannelBaseflow->isChecked()) {
        QbGraph->attach(HPlot);
  //  }

    if(checkIncludeTiledrains->isChecked()) {
        QtileGraph->attach(HPlot);
    }


    if(checkDoErosion->isChecked())
    {
        QsGraph->attach(HPlot);
        CGraph->attach(HPlot);

        HPlot->setAxesCount(QwtPlot::yLeft, 2);
        HPlot->setAxesCount(QwtPlot::yRight, 2);

        QGraph->setAxes(HPlot->xBottom, *axisYL1);
        PGraph->setAxes(HPlot->xBottom, *axisYL2);
        QsGraph->setAxes(HPlot->xBottom, *axisYR1);
        CGraph->setAxes(HPlot->xBottom, *axisYR2);

        HPlot->setAxisTitle(*axisYL1, "Q (l/s)");
        HPlot->setAxisTitle(*axisYL2, "P (mm/h)");
        HPlot->setAxisTitle(*axisYR1, "Qs (kg/s)");
        HPlot->setAxisTitle(*axisYR2, "C (g/l)");
    }
    else
    {
        //        QsGraph->detach();
        //        CGraph->detach();

        QGraph->setAxes(HPlot->xBottom, *axisYL1);
        PGraph->setAxes(HPlot->xBottom, *axisYR1);
        HPlot->setAxesCount(QwtPlot::yLeft, 1);
        HPlot->setAxesCount(QwtPlot::yRight, 1);
        HPlot->setAxisTitle(*axisYL1, "Q (l/s)");
        HPlot->setAxisTitle(*axisYR1, "P (mm/h)");
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
        HPlot->setAxisScale(HPlot->xBottom, op.BeginTime/1440, op.EndTime/1440);
    else
        HPlot->setAxisScale(HPlot->xBottom, op.BeginTime, op.EndTime);
    int index = OutletIndices.indexOf(this->outletpoint);

    QGraph->setSamples(op.Time,*op.OutletQ[index]);
    PGraph->setSamples(op.Time,op.Pmm);
   // if (checkChannelBaseflow->isChecked())
        QbGraph->setSamples(op.Time,*op.OutletQb[index]);

    int _j = op.OutletQ[index]->count()-1;
//    if (qmax[index] < op.OutletQ[index]->at(_j))
//        qmax[index] = op.OutletQ[index]->at(_j);

    qmax[index] = std::max(qmax[index] , 1.1*op.OutletQ[index]->at(_j));
    pmax = std::max(pmax, 1.1*op.Pmm[_j]);

    if(checkDoErosion->isChecked())
    {
        qsmax[index] = std::max(qsmax[index] , 1.1*op.OutletQs[index]->at(_j));
        cmax[index] = std::max(cmax[index] , 1.1*op.OutletC[index]->at(_j));

        QsGraph->setSamples(op.Time,*op.OutletQs[index]);
        CGraph->setSamples(op.Time,*op.OutletC[index]);

        HPlot->setAxisScale(*axisYL1, 0.0, qmax[index] );
        HPlot->setAxisScale(*axisYL2, 0.0, pmax );
        HPlot->setAxisScale(*axisYR1, 0.0, qsmax[index] );
        HPlot->setAxisScale(*axisYR2, 0.0, cmax[index] );
    } else {
        HPlot->setAxisScale(*axisYL1, 0.0, qmax[index] );
        HPlot->setAxisScale(*axisYR1, 0.0, pmax );
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

    outletpoint = 1;
    spinBoxPointtoShow->setValue(1);
    spinBoxPointtoShow->setMaximum(OutletIndices.at(OutletIndices.length()-1));
    label_hydroCount->setText(QString("Output all (0) or point (1-%1)").arg(OutletIndices.count()-1));

    if (outletpoint > 0)
    {
        outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
        if (doNewPlot)
            PQSchart->setTitle(QString("Hydrograph point %1").arg(outletpoint));
        else
            HPlot->setTitle(QString("Hydrograph point %1").arg(outletpoint));
    }
    else
    {
        outletgroup->setTitle(QString("Total domain outflow"));
        if (doNewPlot)
            PQSchart->setTitle("Combined hydrograph domain outflow");
        else
            HPlot->setTitle(QString("Combined hydrograph domain outflow"));
    }
}

//---------------------------------------------------------------------------
// max 6 line text output below hydrographs
void lisemqt::initOutputData()
{

    //textGraph->setMaximumBlockCount(4);
    textGraph->setWordWrapMode(QTextOption::NoWrap);
    textGraph->setMaximumHeight(80);
    textGraph->clear();

    textGraph->setVisible(false);
    label_headerTextGraph->setVisible(false);
}

//---------------------------------------------------------------------------
// max 6 line text output below hydrographs
void lisemqt::SetTextHydrographs()
{

   // textGraph->clear();
    QStringList SL;

    if(op.OutletQ.length() == 0)
    {
        return;
    }

    int j = OutletIndices.indexOf(this->outletpoint);
    int steps = op.OutletQ.at(0)->size();
    times << op.time;
    for(int i = std::max(0,steps-6); i < steps; i++)
    {
        int days = trunc(times.at(i)/1440.0);
        double mins = times.at(i) - (double)days*1440.0;//long(op.time) % 1440;
        double Pmm = op.Pmm.at(i); //Rainfall.at(i);
        double QPlot = op.OutletQ.at(j)->at(i);
        double ChannelWH = op.OutletChannelWH.at(j)->at(i);
        double Qsplot = checkDoErosion->isChecked() ? op.OutletQs.at(j)->at(i) : 0;
        double Cplot = checkDoErosion->isChecked() ? op.OutletC.at(j)->at(i) : 0;

        QString outS;

        // total outflow
        if (j == 0) {
            label_headerTextGraph->setText(QString("%1 %2 %3 %4 %5 %6")
                                           .arg("Total ")
                                           .arg("Time (day:min)   ")
                                           .arg("   Rain (mm)  ")
                                           .arg("      Q (l/s)    ")
                                           .arg("  Qs (kg/s)     ")
                                           .arg("Conc (g/l)    ")
                                           );

            outS = QString("%1    %2:%3   %4 %5 %6 %7")
                                       .arg(j,5)
                                       .arg(days,3,10,QLatin1Char('0'))
                                       .arg(mins,6,'f',1,'0')//QLatin1Char('0'))
//                    .arg(mins,4,10,QLatin1Char('0'))
                                       .arg(Pmm,15,'f',3,' ')
                                       .arg(QPlot,15,'f',3,' ')
                                       .arg(Qsplot,15,'f',3)
                                       .arg(Cplot,15,'f',3,' ')
                                        ;
        } else {
            // point outflow
            label_headerTextGraph->setText(QString("%1 %2 %3 %4 %5 %6 %7")
                                           .arg("Point ")
                                           .arg("Time (day:min)   ")
                                           .arg("   Rain (mm)  ")
                                           .arg("      Q (l/s)    ")
                                           .arg("  ChanWH (m)   ")
                                           .arg("  Qs (kg/s)     ")
                                           .arg("Conc (g/l)    ")
                                            );

            outS = QString("%1    %2:%3   %4 %5 %6 %7 %8")
                                       .arg(j,5)
                                       .arg(days,3,10,QLatin1Char('0'))
                                       .arg(mins,6,'f',1,'0')//QLatin1Char('0'))
                                       //.arg(mins,4,10,QLatin1Char('0'))
                                       .arg(Pmm,15,'f',3,' ')
                                       .arg(QPlot,15,'f',3,' ')
                                       .arg(ChannelWH,15,'f',3,' ')
                                       .arg(Qsplot,15,'f',3)
                                       .arg(Cplot,15,'f',3,' ')
                                        ;
        }
        SL << outS;
    }
 //   textGraph->appendPlainText(SL.join('\n'));

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
    if(darkLISEM)
        format = QString("<font color=#ffffaa>%2</font>");
    else
        format= QString("<font color=#000000>%2</font>");

    if (startplot) {
        label_dx->setText(format.arg(QString::number(op.dx,'f',dig)));
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
    label_raintot->setText(format.arg(QString::number(op.RainTotmm,'f',dig)));
    label_watervoltot->setText(format.arg(QString::number(op.WaterVolTotmm,'f',dig)));
    label_stormdraintot->setText(format.arg(QString::number(op.StormDrainTotmm,'f',dig)));
    label_qtot->setText(format.arg(QString::number(op.Qtotmm,'f',dig)));
    label_infiltot->setText(format.arg(QString::number(op.InfilTotmm,'f',dig)));
    label_surfstor->setText(format.arg(QString::number(op.SurfStormm,'f',dig)));
    label_interctot->setText(format.arg(QString::number(op.IntercTotmm+op.IntercHouseTotmm+op.IntercLitterTotmm,'f',dig)));
    if (checkOverlandFlow1D->isChecked() && !checkIncludeChannel->isChecked())
        label_floodVolmm->setText(format.arg(QString::number(0,'f',dig)));
    else
        label_floodVolmm->setText(format.arg(QString::number(op.volFloodmm,'f',dig)));

    label_watervolchannel->setText(format.arg(QString::number(op.ChannelVolTotmm,'f',dig)));
    label_baseflowtot->setText(format.arg(QString::number(op.BaseFlowtotmm,'f',dig)));
    //  label_litterstore->setText(QString::number(op.LitterStorageTotmm,'f',dig));

    // peak time
    label_QPfrac->setText(format.arg(QString::number((op.RainTotmm > 0 ? op.Qtotmm/op.RainTotmm*100 : 0),'f',dig)));
    label_ppeaktime->setText(format.arg(QString::number(op.RainpeakTime,'f',2)));

    // mass balance
    label_MB->setText(QString::number(op.MB,'e',dig));
    if (op.MB > 0)
        label_MB->setText(" "+label_MB->text());

    int j = OutletIndices.indexOf(this->outletpoint);
    label_qpeaksub->setText(format.arg(QString::number(op.OutletQpeak.at(j),'f',2)));
    label_qpeaktime->setText(format.arg(QString::number(op.OutletQpeaktime.at(j),'f',2)));
    label_qtotm3sub->setText(format.arg(QString::number(op.OutletQtot.at(j),'f',2)));
    label_dischargesub->setText(format.arg(QString::number(op.OutletQ.at(j)->at(op.OutletQ.at(j)->length()-1),'f',2)));

    label_soillosssub->setEnabled(checkDoErosion->isChecked());
    label_94->setEnabled(checkDoErosion->isChecked());

    if(checkDoErosion->isChecked())
    {
        int dig = E_DigitsOut->value();//DIGITS;
        label_MBs->setText(QString::number(op.MBs,'e',dig));
        if (op.MBs > 0)
            label_MBs->setText(" "+label_MBs->text());

        label_soillosssub->setText(format.arg(QString::number(op.OutletQstot.at(j),'f',dig)));
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
    }
}
