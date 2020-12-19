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
    HPlot = new QwtPlot(title, this);
    layout_Plot->insertWidget(0, HPlot, 1);
    //HPlot->canvas()->setFrameStyle( QFrame::StyledPanel);//QFrame::Box | QFrame::Plain );

    // panning with the left mouse button
    (void) new QwtPlotPanner( HPlot->canvas() );

    // zoom in/out with the wheel
    (void) new QwtPlotMagnifier( HPlot->canvas() );

    PGraph = new QwtPlotCurve("Rainfall");
    QGraph = new QwtPlotCurve("Discharge");
    QsGraph = new QwtPlotCurve("Sediment discharge");
    CGraph = new QwtPlotCurve("Concentration");
    QtileGraph = new QwtPlotCurve("Tile drain");

    PGraph->attach(HPlot);
    QGraph->attach(HPlot);
    // order determines order of display in Legend
    //VJ 101223 changed for qwt 6.0.0

    QPen pen1, pen2, pen3, pen4, pen5;
    pen1.setWidth(2);
    col.setRgb( 60,60,200,255 );
    pen1.setColor(col);
    pen1.setCosmetic(true);
    pen2.setWidth(2);
    pen2.setColor("#000000");
    pen2.setCosmetic(false);
    col.setRgb( 0,160,160,255 );
    pen3.setWidth(2);
    pen3.setColor(col);
    pen3.setCosmetic(true);
    pen4.setWidth(2);
    pen4.setColor(Qt::red);
    pen4.setCosmetic(false);
    col.setRgb( 200,0,0,255 ); // darkred
    pen5.setWidth(2);
    pen5.setColor(col);
    pen5.setCosmetic(false);

    QGraph->setPen(pen1);
    QGraph->setAxes(HPlot->xBottom, HPlot->yLeft);

    PGraph->setPen(pen2);
    PGraph->setAxes(HPlot->xBottom, HPlot->yRight);

    QtileGraph->setPen(pen3);
    QtileGraph->setAxes(HPlot->xBottom, HPlot->yLeft);

    QsGraph->setPen(pen4);
    QsGraph->setAxes(HPlot->xBottom, HPlot->yRight);

    CGraph->setPen(pen5);
    CGraph->setAxes(HPlot->xBottom, HPlot->yRight);

    PGraph->setStyle(QwtPlotCurve::Steps);
    QGraph->setStyle(QwtPlotCurve::Lines);
    QtileGraph->setStyle(QwtPlotCurve::Lines);
    QsGraph->setStyle(QwtPlotCurve::Lines);
    CGraph->setStyle(QwtPlotCurve::Lines);

    PGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
    QGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
    QtileGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
    QsGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
    CGraph->setRenderHint(QwtPlotItem::RenderAntialiased);

    HPlot->setCanvasBackground(QBrush(Qt::white));

    // set axes
    HPlot->enableAxis(HPlot->yRight,true);
    HPlot->enableAxis(HPlot->yLeft,true);
    HPlot->enableAxis(HPlot->xBottom,true);
    HPlot->setAxisTitle(HPlot->xBottom, "time (min)");

    HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s)");
    HPlot->setAxisTitle(HPlot->yRight, "P (mm/h)");
    multiplierRain->setValue(0);

    HPlot->setAxisScale(HPlot->yRight, 0, 1);
    HPlot->setAxisScale(HPlot->yLeft, 0, 1);
    HPlot->setAxisScale(HPlot->xBottom, 0, 100);

    // set gridlines
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

            showPlot();
            SetTextHydrographs();

            HPlot->setTitle(QString("Catchment outlet(s)"));
            outletgroup->setTitle(QString("Catchment outlet(s)"));
            //HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
            //outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));

        }else
        {
            if(point > outletpoint)
            {
                index = oldindex + 1;
            }else
            {
                index = oldindex - 1;
            }

            outletpoint = OutletIndices.at(index);
            spinBoxPointtoShow->setValue(outletpoint);

            showPlot();
            SetTextHydrographs();

            if (outletpoint > 0)
            {
                outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
                HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
            }
            else
            {
                outletgroup->setTitle(QString("Catchment outflow"));
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
    // VJ 110630 show hydrograph for selected output point

    if(checkIncludeTiledrains->isChecked())
        QtileGraph->attach(HPlot);
    else
        QtileGraph->detach();

    if(checkDoErosion->isChecked())
    {
        QsGraph->attach(HPlot);
        CGraph->attach(HPlot);
    }
    else
    {
        QsGraph->detach();
        CGraph->detach();
        //       sQsGraph->detach();
    }

    if(checkDoErosion->isChecked())
    {
        PGraph->setAxes(HPlot->xBottom, HPlot->yLeft);
        QsGraph->setAxes(HPlot->xBottom, HPlot->yRight);
        CGraph->setAxes(HPlot->xBottom, HPlot->yRight);

        HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s) - P (mm/h)");
        HPlot->setAxisTitle(HPlot->yRight, "Qs (kg/s) - C (g/l)");
    }
    else
    {
        PGraph->setAxes(HPlot->xBottom, HPlot->yRight);

        HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s)");
        HPlot->setAxisTitle(HPlot->yRight, "P (mm/h)");
    }


}
//---------------------------------------------------------------------------
/// free data structures graph
void lisemqt::killPlot()
{
    OutletQ.clear();
    OutletQs.clear();
    OutletC.clear();
    OutletChannelWH.clear();

    OutletIndices.clear();
    OutletLocationX.clear();
    OutletLocationY.clear();

    Rainfall.clear();
    OutletQpeak.clear();
    OutletQpeaktime.clear();
    OutletQtot.clear();
    OutletQstot.clear();

    PData.clear();
    TData.clear();
    QData.clear();
    QData1.clear();
    QData2.clear();
    QtileData.clear();
    QsData.clear();
    CData.clear();

}
void lisemqt::GetPlotData()
{   
    QtileData << op.Qtile;
//    PData << op.Pmm*qPow(10.0, multiplierRain->value());
    TData << op.time;

    for(int i = 0; i < OutletIndices.length(); i++)
    {
        OutletQ.at(i)->clear();
        OutletQs.at(i)->clear();
        OutletC.at(i)->clear();
        OutletChannelWH.at(i)->clear();

        OutletQ.at(i)->append(*op.OutletQ.at(i));
        OutletQs.at(i)->append(*op.OutletQs.at(i));
        OutletC.at(i)->append(*op.OutletC.at(i));
        OutletChannelWH.at(i)->append(*op.OutletChannelWH.at(i));
    }

    Rainfall.append(op.Pmm);

    OutletQpeak.clear();
    OutletQpeaktime.clear();
    OutletQpeak.append(op.OutletQpeak);
    OutletQpeaktime.append(op.OutletQpeaktime);
    OutletQtot.clear();
    OutletQstot.clear();
    OutletQtot.append(op.OutletQtot);
    OutletQstot.append(op.OutletQstot);
    timestep = op.timestep;

//    double mf = 0;
//     for(int i = 0; i < 6; i++) {
//         mf = (double)i;
//         if(op.Rainpeak*qPow(10.0, mf) > op.OutletQpeak.at(0))
//             break;
//     }
//     mf -= 1.0;
     PData << op.Pmm*mult;//qPow(10.0, mf);
}

//---------------------------------------------------------------------------
void lisemqt::on_multiplierRain_valueChanged(double)
{
    double mult = qPow(10.0, multiplierRain->value());
    if(multiplierRain->value() > 0)
        HPlot->setAxisTitle(HPlot->yLeft, QString("Q (l/s) - P (%1 mm/h)").arg(1/mult));
    else
        HPlot->setAxisTitle(HPlot->yLeft, QString("Q (l/s) - P (mm/h)"));
}
//---------------------------------------------------------------------------
void lisemqt::showPlot()
{
  //  double mult = qPow(10.0, multiplierRain->value());
    mult = 1;
    int i;
   double mf[6] ={1,10,100,1000,10000,1000000};
     for(i = 0; i < 6; i++) {
         if(op.Rainpeak*mf[i] > op.OutletQpeak.at(0)) {
             break;
         }
     }
     mult = mf[i-1]*0.5;
  //   qDebug() << mult << i << op.Rainpeak*mult << op.OutletQpeak.at(0);

    QData.clear();
    QsData.clear();
    CData.clear();
    PData.clear();

    int index = OutletIndices.indexOf(this->outletpoint);

    for(int i = 0; i < OutletQ.at(index)->length();i++)
    {
        PData << Rainfall.at(i)*mult;
        QData << OutletQ.at(index)->at(i);
        QsData <<OutletQs.at(index)->at(i);
        CData << OutletC.at(index)->at(i);

        qmax.replace(index,OutletQ.at(index)->at(i) > qmax.at(index)? OutletQ.at(index)->at(i) : qmax.at(index));
        qsmax.replace(index,OutletQs.at(index)->at(i) > qsmax.at(index)? OutletQs.at(index)->at(i) : qsmax.at(index));
        cmax.replace(index,OutletC.at(index)->at(i) > cmax.at(index)? OutletC.at(index)->at(i) : cmax.at(index));
    }

    QGraph->setSamples(TData,QData);

    yas = std::max(0.01,qmax.at(index));
    yasP = std::max(yasP, op.Pmm*mult);

    PGraph->setSamples(TData,PData);

    HPlot->setAxisScale(HPlot->yLeft, 0, yas*1.05);

    if(checkDoErosion->isChecked())
    {
        QsGraph->setSamples(TData,QsData);
        CGraph->setSamples(TData,CData);
        y2as = std::max(0.1,std::max(qsmax.at(index), cmax.at(index)));
        HPlot->setAxisScale(HPlot->yRight, 0, y2as*1.05);

        yas = std::max(0.1,std::max(yas, op.Pmm*mult));
    }
    else
    {
        y2as = std::max(0.1,std::max(y2as, op.Pmm*mult));
        HPlot->setAxisScale(HPlot->yRight, 0, y2as*1.05);
    }


    if(checkIncludeTiledrains->isChecked())
        QtileGraph->setSamples(TData,QtileData);

    HPlot->replot();

}
//---------------------------------------------------------------------------
// initializes plot with world data, called in worldshow, done once
void lisemqt::startPlots()
{
    if (!startplot)
        return;

    yas = 0.1;
    yasP = 0;
    y2as = 0.1;

    qmax.clear();
    qsmax.clear();
    cmax.clear();

    killPlot();
    // clear() plot data

    HPlot->setAxisScale(HPlot->xBottom, op.BeginTime, op.EndTime);

    QwtLegend *legend = new QwtLegend(HPlot);
    legend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
    HPlot->insertLegend(legend, QwtPlot::BottomLegend);
    //legend

    OutletIndices.append(op.OutletIndices);
    OutletLocationX.append(op.OutletLocationX);
    OutletLocationY.append(op.OutletLocationY);
    OutletQtot.append(op.OutletQtot);
    OutletQstot.append(op.OutletQstot);

    for(int i =0; i < OutletIndices.length(); i++)
    {
        OutletQ.append(new QList<double>);
        OutletQs.append(new QList<double>);
        OutletC.append(new QList<double>);
        OutletChannelWH.append(new QList<double>);
        qmax.append(0);
        qsmax.append(0);
        cmax.append(0);
    }

    outletpoint = 1;
    spinBoxPointtoShow->setValue(1);
    spinBoxPointtoShow->setMaximum(OutletIndices.at(OutletIndices.length()-1));
    label_hydroCount->setText(QString("Hydrograph Point (0, 1-%1)").arg(OutletIndices.count()-1));

    if (outletpoint > 0)
    {
        outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
        HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
    }
    else
    {
        outletgroup->setTitle(QString("Catchment outflow"));
        HPlot->setTitle(QString("Hydrograph all outflow"));
    }
    // VJ 110630 show hydrograph for selected output point

}

//---------------------------------------------------------------------------
// max 6 line text output below hydrographs
void lisemqt::initOutputData()
{
    //textGraph->setMaximumBlockCount(4);
    textGraph->setWordWrapMode(QTextOption::NoWrap);
    textGraph->setMaximumHeight(80);
    textGraph->clear();
}

//---------------------------------------------------------------------------
// max 6 line text output below hydrographs
void lisemqt::SetTextHydrographs()
{
    textGraph->clear();

    if(OutletQ.length() == 0)
    {
        return;
    }

    int j = OutletIndices.indexOf(this->outletpoint);

    int dig = E_DigitsOut->value(); //DIGITS;

    label_qpeaksub->setText(QString::number(OutletQpeak.at(j),'e',dig));
    label_qpeaktime->setText(QString::number(OutletQpeaktime.at(j),'e',dig));
    label_qtotm3sub->setText(QString::number(OutletQtot.at(j),'e',dig));
    label_dischargesub->setText(QString::number(OutletQ.at(j)->at(OutletQ.at(j)->length()-1),'e',dig));

    if(checkDoErosion->isChecked())
        label_soillosssub->setText(QString::number(OutletQstot.at(j),'f',dig));

    int steps = OutletQ.at(0)->length();
    for(int i = 0; i < steps; i++)
    {
        double time = timestep* i;
        double Pmm = Rainfall.at(i);
        double QPlot = OutletQ.at(j)->at(i);
        double ChannelWH = OutletChannelWH.at(j)->at(i);
        double Qsplot = OutletQs.at(j)->at(i);
        double Cplot = OutletC.at(j)->at(i);

        if(checkDoErosion->isChecked())
        {
            if(!checkIncludeTiledrains->isChecked())
                textGraph->appendPlainText(QString("%1 %2 %3 %4 %5 %6")
                                           .arg(time,15,'f',3,' ')
                                           .arg(Pmm,15,'f',3,' ')
                                           .arg(QPlot,15,'f',3,' ')
                                           .arg(ChannelWH,15,'f',3,' ')
                                           .arg(Qsplot,12,'f',3)
                                           .arg(Cplot,15,'f',3,' '));
            else
                textGraph->appendPlainText(QString("%1 %2 %3 %4 %5 %6 %7")
                                           .arg(time,15,'f',3,' ')
                                           .arg(Pmm,15,'f',3,' ')
                                           .arg(QPlot,15,'f',3,' ')
                                           .arg(ChannelWH,15,'f',3,' ')
                                           .arg(Qsplot,12,'f',3)
                                           .arg(Cplot,15,'f',3,' ')
                                           .arg(0.0,15,'f',3,' '));
        }
        else
        {
            double dummy = 0.0;
            if(!checkIncludeTiledrains->isChecked())
                textGraph->appendPlainText(QString("%1 %2 %3 %4 %5 %6")
                                           .arg(time,15,'f',3,' ')
                                           .arg(Pmm,15,'f',3,' ')
                                           .arg(QPlot,15,'f',3,' ')
                                           .arg(ChannelWH,15,'f',3,' ')
                                           .arg(dummy,15,'f',3,' ')
                                           .arg(dummy,15,'f',3,' ')
                                           );
            else
                textGraph->appendPlainText(QString("%1 %2 %3 %4 %5 %6 %7")
                                           .arg(time,15,'f',3,' ')
                                           .arg(Pmm,15,'f',3,' ')
                                           .arg(QPlot,15,'f',3,' ')
                                           .arg(ChannelWH,15,'f',3,' ')
                                           .arg(0.0,15,'f',3,' ')
                                           .arg(dummy,15,'f',3,' ')
                                           .arg(dummy,15,'f',3,' ')
                                           );
        }
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

    label_dx->setText(QString::number(op.dx,'f',dig));
    label_area->setText(QString::number(op.CatchmentArea/1000000,'f',dig));
    label_time->setText(QString::number(op.time,'f',dig));
    label_endtime->setText(QString::number(op.EndTime,'f',dig));
    label_runtime->setText(QString::number(op.t,'f',dig));
    label_endruntime->setText(QString::number(op.maxtime,'f',dig));

    // mass balance
    label_MB->setText(QString::number(op.MB,'e',dig));
    if (op.MB > 0)
        label_MB->setText(" "+label_MB->text());

    // mm output
    label_raintot->setText(QString::number(op.RainTotmm,'f',dig));

    label_watervoltot->setText(QString::number(op.WaterVolTotmm,'f',dig));
    label_stormdraintot->setText(QString::number(op.StormDrainTotmm,'f',dig));
    label_qtot->setText(QString::number(op.Qtotmm,'f',dig));
    label_infiltot->setText(QString::number(op.InfilTotmm,'f',dig));
    label_surfstor->setText(QString::number(op.SurfStormm,'f',dig));
    label_interctot->setText(QString::number(op.IntercTotmm+op.IntercHouseTotmm+op.IntercLitterTotmm,'f',dig));
    if (checkOverlandFlow1D->isChecked() && !checkIncludeChannel->isChecked())
        label_floodVolmm->setText(QString::number(0,'f',dig));
    else
        label_floodVolmm->setText(QString::number(op.volFloodmm,'f',dig));

    label_watervolchannel->setText(QString::number(op.ChannelVolTotmm,'f',dig));
    label_baseflowtot->setText(QString::number(op.BaseFlowtotmm,'f',dig));
    //  label_litterstore->setText(QString::number(op.LitterStorageTotmm,'f',dig));

    // outlet
    //    label_qtotm3->setText(QString::number(op.Qtot,'f',2));
    //    label_discharge->setText(QString::number(op.Q,'f',2));

    // peak time
    label_QPfrac->setText(QString::number((op.RainTotmm > 0 ? op.Qtotmm/op.RainTotmm*100 : 0),'f',dig));
    label_ppeaktime->setText(QString::number(op.RainpeakTime,'f',dig));

    if(checkDoErosion->isChecked())
    {
        int dig = E_DigitsOut->value();//DIGITS;
        label_MBs->setText(QString::number(op.MBs,'e',dig));
        if (op.MBs > 0)
            label_MBs->setText(" "+label_MBs->text());

        label_splashdet->setText(QString::number(op.DetTotSplash,'f',dig));
        label_flowdet->setText(QString::number(op.DetTotFlow+op.FloodDetTot,'f',dig));
        label_sedvol->setText(QString::number(op.SedTot+op.FloodSedTot,'f',dig));
        label_dep->setText(QString::number(op.DepTot+op.FloodDepTot,'f',dig));

        label_detch->setText(QString::number(op.ChannelDetTot,'f',dig));
        label_depch->setText(QString::number(op.ChannelDepTot,'f',dig));
        label_sedvolch->setText(QString::number(op.ChannelSedTot,'f',dig));

//        label_flooddet->setText(QString::number(op.FloodDetTot,'f',dig));
//        label_flooddep->setText(QString::number(op.FloodDepTot,'f',dig));
//        label_floodsed->setText(QString::number(op.FloodSedTot,'f',dig));

        label_soilloss->setText(QString::number(op.SoilLossTot,'f',dig));
        label_soillosskgha->setText(QString::number(op.SoilLossTot/(op.CatchmentArea/10000)*1000,'f',dig));

        double SDR = op.DetTotSplash + op.ChannelDetTot + op.DetTotFlow;
        SDR = (SDR > 0? 100*op.SoilLossTot/(SDR) : 0);
        SDR = std::min(SDR ,100.0);
        label_SDR->setText(QString::number(SDR,'f',dig));
    }
}
