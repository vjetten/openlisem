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

#include "lisemqt.h"
#include "global.h"


//---------------------------------------------------------------------------
/// initialize graph plotting
void lisemqt::initPlot()
{
  startplot = true;

  op.outputpointnr = spinBoxPointtoShow->value();
  spinBoxPointtoShow->setEnabled(false);
  HPlot->setTitle(op.outputpointdata);//QString("Hydrograph point %1").arg(op.outputpointnr));
  // VJ 110630 show hydrograph for selected output point
  label_qtotm3sub->setEnabled(op.outputpointnr > 1);

  textGraph->setMaximumBlockCount(6);
  textGraph->setWordWrapMode(QTextOption::NoWrap);
  textGraph->setMaximumHeight(90);
}
//---------------------------------------------------------------------------
/// free data structures graph
void lisemqt::killPlot()
{
  PData.clear();
  TData.clear();
  QData.clear();
  QtileData.clear();
  QsData.clear();
  CData.clear();

  spinBoxPointtoShow->setEnabled(true);

  startplot = true;

}
//---------------------------------------------------------------------------
/// set up discharge plot, graphics part (not data)
void lisemqt::setupPlot()
{
  QColor col;
  QwtText title;
  title.setText("Hydrograph outlet");
  title.setFont(QFont("MS Shell Dlg 2",12));
  HPlot = new QwtPlot(title, this);
  layout_Plot->insertWidget(0, HPlot, 1);
  HPlot->canvas()->setFrameStyle( QFrame::StyledPanel);//QFrame::Box | QFrame::Plain );

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
  PGraph->setAxes(HPlot->xBottom, HPlot->yLeft);
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
  HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s) - P (mm/h)");
  HPlot->setAxisTitle(HPlot->yRight, "Qs (kg/s) - C (g/l)");
  HPlot->setAxisScale(HPlot->yRight, 0, 1);
  HPlot->setAxisScale(HPlot->yLeft, 0, 100);
  HPlot->setAxisScale(HPlot->xBottom, 0, 100);

  // set gridlines
  QwtPlotGrid *grid = new QwtPlotGrid();
  grid->enableXMin(true);
  grid->enableYMin(true);
  col.setRgb( 180,180,180,180 );
  grid->setMajPen(QPen(col, 0, Qt::DashLine));
  col.setRgb( 210,210,210,180 );
  grid->setMinPen(QPen(col, 0 , Qt::DotLine));
  grid->attach(HPlot);


  HPlot->replot();
  // draw empty plot

}
//---------------------------------------------------------------------------
void lisemqt::showPlot()
{

  QData << op.Q;
  QtileData << op.Qtile;
  QsData << op.Qs;
  CData << op.C;
  PData << op.P;
  TData << op.time;

  QGraph->setSamples(TData,QData);
  PGraph->setSamples(TData,PData);
  if(!checkNoErosion->isChecked())
    {
      QsGraph->setSamples(TData,QsData);
      CGraph->setSamples(TData,CData);
      y2as = max(y2as, op.Qs);
      y2as = max(y2as, op.C);
      HPlot->setAxisScale(HPlot->yRight, 0, y2as*1.05);
    }
  if(checkIncludeTiledrains->isChecked())
    QtileGraph->setSamples(TData,QtileData);


  yas = max(yas, op.Q);
  yas = max(yas, op.P);
  HPlot->setAxisScale(HPlot->yLeft, 0, yas*1.05);

  HPlot->replot();

}
//---------------------------------------------------------------------------
/// set up discharge plot, graphics part (not data)
void lisemqt::setupSmallPlot()
{
  QwtText title;
  title.setText("Hydrograph");
  title.setFont(QFont("MS Shell Dlg 2",8));
  smallPlot = new QwtPlot(title, this);
  verticalLayout_6->insertWidget(0, smallPlot, 1);
  smallPlot->canvas()->setFrameStyle( QFrame::StyledPanel);

  sPGraph = new QwtPlotCurve("Rainfall");
  sQGraph = new QwtPlotCurve("Discharge");

  sPGraph->attach(smallPlot);
  sQGraph->attach(smallPlot);
  // order determines order of display in Legend
  //VJ 101223 changed for qwt 6.0.0
  sPGraph->setAxes(smallPlot->xBottom, smallPlot->yLeft);
  sQGraph->setAxes(smallPlot->xBottom, smallPlot->yLeft);

  // do not attach yet
  if(!checkNoErosion->isChecked())
    {
      sQsGraph = new QwtPlotCurve("Sediment discharge");
      sQsGraph->setAxes(smallPlot->xBottom, smallPlot->yRight);
      sQsGraph->setPen(QPen(Qt::red));
      sQsGraph->setStyle(QwtPlotCurve::Lines);
      sQsGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
    }

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
  if(!checkNoErosion->isChecked())
    smallPlot->enableAxis(smallPlot->yRight,true);

  title.setText("time (min)");
  title.setFont(QFont("MS Shell Dlg 2",8));
  smallPlot->setAxisTitle(smallPlot->xBottom, title);
  title.setText("Q (l/s)/P (mm/h)");
  smallPlot->setAxisTitle(smallPlot->yLeft, title);
  smallPlot->setAxisScale(smallPlot->yLeft, 0, 100);
  smallPlot->setAxisScale(smallPlot->xBottom, 0, 100);
  if(!checkNoErosion->isChecked())
    {
      title.setText("Qs (kg/s)");
      smallPlot->setAxisTitle(smallPlot->yRight, title);
      smallPlot->setAxisScale(smallPlot->yRight, 0, 1);
    }

  // set gridlines

  QwtPlotGrid *grid = new QwtPlotGrid();
  col.setRgb( 180,180,180,180 );
  grid->setMajPen(QPen(col, 0, Qt::DotLine));
  grid->attach(smallPlot);

  smallPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::showSmallPlot()
{
  sQGraph->setSamples(TData,QData);
  sPGraph->setSamples(TData,PData);
  if(!checkNoErosion->isChecked())
    sQsGraph->setSamples(TData,QsData);

  smallPlot->setTitle(op.outputpointdata);//QString("Hydrograph point %1").arg(op.outputpointnr));

  smallPlot->setAxisScale(smallPlot->yLeft, 0, yas*1.05);
  if(!checkNoErosion->isChecked())
    smallPlot->setAxisScale(smallPlot->yRight, 0, y2as*1.05);

  smallPlot->replot();

}
//---------------------------------------------------------------------------
void lisemqt::startPlots()
{
  if (startplot)
    {
      startplot = false;

      yas = 0.1;
      y2as = 0.1;

      PData.clear();
      TData.clear();
      QData.clear();
      QtileData.clear();
      QsData.clear();
      CData.clear();

      HPlot->setAxisScale(HPlot->xBottom, op.BeginTime, op.EndTime);

      smallPlot->setAxisScale(smallPlot->xBottom, op.BeginTime, op.EndTime);

      if(checkIncludeTiledrains->isChecked())
        QtileGraph->attach(HPlot);

      if(!checkNoErosion->isChecked())
        {
          QsGraph->attach(HPlot);
          CGraph->attach(HPlot);

          sQsGraph->attach(smallPlot);
        }
      else
        {
          HPlot->enableAxis(HPlot->yRight,false);
          smallPlot->enableAxis(smallPlot->yRight,false);
        }

      if(checkNoErosion->isChecked())
        {
          HPlot->setAxisTitle(HPlot->yRight, "");
          HPlot->setAxisScale(HPlot->yRight, 0, 1);
          smallPlot->setAxisTitle(HPlot->yRight, "");
          smallPlot->setAxisScale(HPlot->yRight, 0, 1);
        }

      QwtLegend *legend = new QwtLegend(HPlot);
      legend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
      HPlot->insertLegend(legend, QwtPlot::BottomLegend);
      //legend
      QwtLegend *slegend = new QwtLegend(smallPlot);
      slegend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
      smallPlot->insertLegend(slegend, QwtPlot::BottomLegend);
      //legend

      label_pointOutput->setText(op.outputpointdata);
      // VJ 110630 show hydrograph for selected output point

      HPlot->setTitle(op.outputpointdata);//QString("Hydrograph point %1").arg(op.outputpointnr));
      // VJ 110630 show hydrograph for selected output point
    }


}

//---------------------------------------------------------------------------
