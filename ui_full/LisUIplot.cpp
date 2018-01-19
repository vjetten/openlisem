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
 \file LisUIPlot.cpp
 \brief plot the hydrographs using the Qwt library

  This file contains all the functions to plot the hydrographs
  it uses the op data for each next step


 */

#include <algorithm>
#include "lisemqt.h"
#include "global.h"

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

    col.setRgb( 60,60,200,255 );
    QGraph->setPen(QPen(col));
    PGraph->setPen(QPen("#000000"));

    PGraph->setAxes(HPlot->xBottom, HPlot->yRight);

    QGraph->setAxes(HPlot->xBottom, HPlot->yLeft);

    QtileGraph->setAxes(HPlot->xBottom, HPlot->yLeft);
    col.setRgb( 0,160,160,255 );
    QtileGraph->setPen(QPen(col));

    QsGraph->setAxes(HPlot->xBottom, HPlot->yRight);
    CGraph->setAxes(HPlot->xBottom, HPlot->yRight);
    QsGraph->setPen(QPen(Qt::red));
    col.setRgb( 200,0,0,255 ); // darkred
    CGraph->setPen(QPen(col));

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

    HPlot->setAxisScale(HPlot->yRight, 0, 1);
    HPlot->setAxisScale(HPlot->yLeft, 0, 1);
    HPlot->setAxisScale(HPlot->xBottom, 0, 100);

    // set gridlines
    QwtPlotGrid *grid = new QwtPlotGrid();
//    grid->enableXmin(true);
 //   grid->enableYmin(true);
    col.setRgb( 180,180,180,180 );
    grid->setMajorPen(QPen(col, 0, Qt::DashLine));
    col.setRgb( 210,210,210,180 );
    grid->setMinorPen(QPen(col, 0 , Qt::DotLine));
    grid->attach(HPlot);

    HPlot->replot();
    // draw empty plot
}
//---------------------------------------------------------------------------
/// set up small discharge plot on mapplot page
/// this is called at the start of lisem, initsmallplot() is called at the start of a model run
/*void lisemqt::setupSmallPlot()
{

    QwtText title;
    title.setText("Hydrograph");
    title.setFont(QFont("MS Shell Dlg 2",8));
    smallPlot = new QwtPlot(title, this);
    smallPlot->setMinimumSize(300,300);
    smallPlot->resize(500,500);
    verticalLayout_6->insertWidget(0, smallPlot, 1);
    smallPlot->canvas()->setFrameStyle( QFrame::StyledPanel);

    sPGraph = new QwtPlotCurve("Rainfall");
    sQGraph = new QwtPlotCurve("Discharge");

    sPGraph->attach(smallPlot);
    sQGraph->attach(smallPlot);
    // order determines order of display in Legend
    //VJ 101223 changed for qwt 6.0.0

    sPGraph->setAxes(smallPlot->xBottom, smallPlot->yRight);
    sQGraph->setAxes(smallPlot->xBottom, smallPlot->yLeft);

    sQsGraph = new QwtPlotCurve("Sediment discharge");
    sQsGraph->setPen(QPen(Qt::red));
    sQsGraph->setStyle(QwtPlotCurve::Lines);
    sQsGraph->setRenderHint(QwtPlotItem::RenderAntialiased);

    QColor col;
    col.setRgb( 60,60,200,255 );
    sQGraph->setPen(QPen(col));
    sPGraph->setPen(QPen("#000000"));

    sPGraph->setStyle(QwtPlotCurve::Steps);
    sQGraph->setStyle(QwtPlotCurve::Lines);

    sPGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
    sQGraph->setRenderHint(QwtPlotItem::RenderAntialiased);

    smallPlot->setCanvasBackground(QBrush(Qt::white));

    // set axes
    smallPlot->enableAxis(smallPlot->yLeft,true);
    smallPlot->enableAxis(smallPlot->xBottom,true);
    smallPlot->enableAxis(smallPlot->yRight,true);

    title.setText("time (min)");
    title.setFont(QFont("MS Shell Dlg 2",8));
    smallPlot->setAxisTitle(smallPlot->xBottom, title);
    title.setText("Q (l/s)");
    smallPlot->setAxisTitle(smallPlot->yLeft, title);
    smallPlot->setAxisScale(smallPlot->yLeft, 0, 1);
    smallPlot->setAxisScale(smallPlot->xBottom, 0, 100);
    title.setText("P (mm/h)");
    smallPlot->setAxisTitle(smallPlot->yRight, title);
    smallPlot->setAxisScale(smallPlot->yRight, 0, 1);

    // set gridlines

    QwtPlotGrid *grid = new QwtPlotGrid();
    col.setRgb( 180,180,180,180 );
    grid->setMajPen(QPen(col, 0, Qt::DotLine));
    grid->attach(smallPlot);

    smallPlot->replot();

}*/
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
            spinBoxPointtoShow->setValue(0);

            showPlot();
            SetTextHydrographs();

            HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
            outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));

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

            HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
            outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
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

    //    if(!checkNoErosion->isChecked())
    if(checkDoErosion->isChecked())
    {
        QsGraph->attach(HPlot);
        CGraph->attach(HPlot);

        //sQsGraph->attach(smallPlot);
    }
    else
    {
        QsGraph->detach();
        CGraph->detach();
 //       sQsGraph->detach();
    }

    //    if(!checkNoErosion->isChecked())
    if(checkDoErosion->isChecked())
    {
        PGraph->setAxes(HPlot->xBottom, HPlot->yLeft);
        QsGraph->setAxes(HPlot->xBottom, HPlot->yRight);
        CGraph->setAxes(HPlot->xBottom, HPlot->yRight);

   //     sPGraph->setAxes(smallPlot->xBottom, smallPlot->yLeft);
   //     sQsGraph->setAxes(smallPlot->xBottom, smallPlot->yRight);

        HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s) - P (mm/h)");
        HPlot->setAxisTitle(HPlot->yRight, "Qs (kg/s) - C (g/l)");
    //    smallPlot->setAxisTitle(smallPlot->yLeft, "Q (l/s) & P (mm/h)");
   //     smallPlot->setAxisTitle(smallPlot->yRight, "Qs (kg/s)");
    }
    else
    {
        PGraph->setAxes(HPlot->xBottom, HPlot->yRight);
 //       sPGraph->setAxes(smallPlot->xBottom, smallPlot->yRight);


        HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s)");
        HPlot->setAxisTitle(HPlot->yRight, "P (mm/h)");
   //    smallPlot->setAxisTitle(smallPlot->yLeft, "Q (l/s)");
   //     smallPlot->setAxisTitle(smallPlot->yRight, "P (mm/h)");
    }

}
//---------------------------------------------------------------------------
/// free data structures graph
void lisemqt::killPlot()
{
    /*for(int i =OutletQ.length()-1; i >-1 ; i--)
    {
        delete OutletQ.at(i);
    }
    for(int i =OutletQs.length()-1; i >-1 ; i--)
    {
        delete OutletQs.at(i);
    }
    for(int i =OutletC.length()-1; i >-1 ; i--)
    {
        delete OutletC.at(i);
    }
    for(int i =OutletChannelWH.length()-1; i >-1 ; i--)
    {
        delete OutletChannelWH.at(i);
    }*/

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
    PData << op.Pmm;
    TData << op.time;

    for(int i =0; i < OutletIndices.length(); i++)
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
}

//---------------------------------------------------------------------------
void lisemqt::showPlot()
{

    QData.clear();
    QsData.clear();
    CData.clear();


    int index = OutletIndices.indexOf(this->outletpoint);

    for(int i = 0; i < OutletQ.at(index)->length();i++)
    {
        QData << OutletQ.at(index)->at(i);
        QsData <<OutletQs.at(index)->at(i);
        CData << OutletC.at(index)->at(i);

        qmax.replace(index,OutletQ.at(index)->at(i) > qmax.at(index)? OutletQ.at(index)->at(i) : qmax.at(index));
        qsmax.replace(index,OutletQs.at(index)->at(i) > qsmax.at(index)? OutletQs.at(index)->at(i) : qsmax.at(index));
        cmax.replace(index,OutletC.at(index)->at(i) > cmax.at(index)? OutletC.at(index)->at(i) : cmax.at(index));
    }

    QGraph->setSamples(TData,QData);

    yas = std::max(0.1,qmax.at(index));
    yasP = std::max(yasP, op.Pmm);

    PGraph->setSamples(TData,PData);

    HPlot->setAxisScale(HPlot->yLeft, 0, yas*1.05);

    if(checkDoErosion->isChecked())
    {
        QsGraph->setSamples(TData,QsData);
        CGraph->setSamples(TData,CData);
        y2as = std::max(0.1,std::max(qsmax.at(index), cmax.at(index)));
        HPlot->setAxisScale(HPlot->yRight, 0, y2as*1.05);

        yas = std::max(0.1,std::max(yas, op.Pmm));
    }
    else
    {
        y2as = std::max(0.1,std::max(y2as, op.Pmm));
        HPlot->setAxisScale(HPlot->yRight, 0, y2as*1.05);
    }


    if(checkIncludeTiledrains->isChecked())
        QtileGraph->setSamples(TData,QtileData);

    HPlot->replot();

}
//---------------------------------------------------------------------------
/*void lisemqt::showSmallPlot()
{

    sQGraph->setSamples(TData,QData);
    sPGraph->setSamples(TData,PData);
    //    if(!checkNoErosion->isChecked())
    if(checkDoErosion->isChecked())
        sQsGraph->setSamples(TData,QsData);

    smallPlot->setTitle(QString("Hydrograph %1").arg(op.outputpointdata));

    smallPlot->setAxisScale(smallPlot->yLeft, 0, yas*1.05);

    smallPlot->setAxisScale(smallPlot->yRight, 0, y2as*1.05);

    smallPlot->replot();


}*/
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

 //   smallPlot->setAxisScale(smallPlot->xBottom, op.BeginTime, op.EndTime);

    QwtLegend *legend = new QwtLegend(HPlot);
    legend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
    HPlot->insertLegend(legend, QwtPlot::BottomLegend);
    //legend
   // QwtLegend *slegend = new QwtLegend(smallPlot);
   // slegend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
  //  smallPlot->insertLegend(slegend, QwtPlot::BottomLegend);
    //legend

 //   label_pointOutput->setText(QString("Hydrograph %1").arg(op.outputpointdata));
    // VJ 110630 show hydrograph for selected output point

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

    outletpoint = 0;
    spinBoxPointtoShow->setValue(0);
    spinBoxPointtoShow->setMaximum(OutletIndices.at(OutletIndices.length()-1));


    if (outletpoint == 0)
    outletgroup->setTitle(QString("Catchment outflow (all)"));
    else
    outletgroup->setTitle(QString("Catchment outlet %1").arg(outletpoint));
    HPlot->setTitle(QString("Hydrograph %1").arg(outletpoint));
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
    //    if(!checkNoErosion->isChecked())
    textGraph->clear();

    if(OutletQ.length() == 0)
    {
        return;
    }
    int j = OutletIndices.indexOf(this->outletpoint);

    label_qpeaksub->setText(QString::number(OutletQpeak.at(j),'f',1));
    label_qpeaktime->setText(QString::number(OutletQpeaktime.at(j),'f',3));
    label_qtotm3sub->setText(QString::number(OutletQtot.at(j),'f',1));
    label_dischargesub->setText(QString::number(OutletQ.at(j)->at(OutletQ.at(j)->length()-1),'f',1));

    if(checkDoErosion->isChecked())
        label_soillosssub->setText(QString::number(OutletQstot.at(j),'f',2));

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
            if(!checkIncludeTiledrains->isChecked())
                textGraph->appendPlainText(QString("%1 %2 %3 %4    --           --")
                                           .arg(time,15,'f',3,' ')
                                           .arg(Pmm,15,'f',3,' ')
                                           .arg(QPlot,15,'f',3,' ')
                                           .arg(ChannelWH,15,'f',3,' '));
            else
                textGraph->appendPlainText(QString("%1 %2 %3 %4 %5     --           --")
                                           .arg(time,15,'f',3,' ')
                                           .arg(Pmm,15,'f',3,' ')
                                           .arg(QPlot,15,'f',3,' ')
                                           .arg(ChannelWH,15,'f',3,' ')
                                           .arg(0.0,15,'f',3,' '));
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

    label_dx->setText(QString::number(op.dx,'f',3));
    label_area->setText(QString::number(op.CatchmentArea/1000000,'f',3));
    label_time->setText(QString::number(op.time,'f',3));
    label_endtime->setText(QString::number(op.EndTime,'f',3));
    label_runtime->setText(QString::number(op.t,'f',3));
    label_endruntime->setText(QString::number(op.maxtime,'f',3));

    // mass balance
    label_MB->setText(QString::number(op.MB,'e',3));

    // mm output
    label_raintot->setText(QString::number(op.RainTotmm,'f',3));
    label_watervoltot->setText(QString::number(op.WaterVolTotmm-op.volFloodmm,'f',3));
    label_qtot->setText(QString::number(op.Qtotmm,'f',3));
    label_infiltot->setText(QString::number(op.InfilTotmm,'f',3));
    label_surfstor->setText(QString::number(op.SurfStormm,'f',3));
    label_interctot->setText(QString::number(op.IntercTotmm+op.IntercHouseTotmm,'f',3));
    label_watervolchannel->setText(QString::number(op.ChannelVolTot,'f',3));
    label_baseflowtot->setText(QString::number(op.BaseFlowtot,'f',3));
    label_floodVolmm->setText(QString::number(op.volFloodmm,'f',3));
    //label_litterstore->setText(QString::number(op.LitterStorageTot,'f',3));

 //   if (op.outputpointnr > 1)
 //   {
 //   }

    // outlet
//    label_qtotm3->setText(QString::number(op.Qtot,'f',2));
//    label_discharge->setText(QString::number(op.Q,'f',2));

    // peak time
    label_QPfrac->setText(QString::number((op.RainTotmm > 0 ? op.Qtotmm/op.RainTotmm*100 : 0),'f',3));
    label_ppeaktime->setText(QString::number(op.RainpeakTime,'f',3));

    if(checkDoErosion->isChecked())
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
        SDR = std::min(SDR ,100.0);
        label_SDR->setText(QString::number(SDR,'f',dig));
    }




}
