#include "lisemqt.h"
#include "global.h"


//---------------------------------------------------------------------------
void lisemqt::initPlot()
{
   startplot = true;
   QData = NULL;
   QsData = NULL;
   CData = NULL;
   PData = NULL;
   timeData = NULL;
   //intialize plot stuff for this run
}
//---------------------------------------------------------------------------
// free data structures graph
void lisemqt::killPlot()
{
   delete QData;
   delete QsData;
   delete CData;
   delete PData;
   delete timeData;
   QData = NULL;
   QsData = NULL;
   CData = NULL;
   PData = NULL;
   timeData = NULL;
}
//---------------------------------------------------------------------------
void lisemqt::setupPlot()
{
   textGraph->setMaximumBlockCount(6);
   textGraph->setWordWrapMode(QTextOption::NoWrap);
   textGraph->setMaximumHeight(96);

   QwtText title;
   title.setText("Hydrograph/Sedigraph outlet");
   HPlot = new QwtPlot(title, this);
   // make the plot window
   verticalLayout_6->insertWidget(0, HPlot);

   // panning with the left mouse button
   (void) new QwtPlotPanner( HPlot->canvas() );

   // zoom in/out with the wheel
   (void) new QwtPlotMagnifier( HPlot->canvas() );

   PGraph = new QwtPlotCurve("Rainfall");
   QGraph = new QwtPlotCurve("Discharge");
   QsGraph = new QwtPlotCurve("Sediment discharge");
   CGraph = new QwtPlotCurve("Concentration");

   PGraph->attach(HPlot);
   QGraph->attach(HPlot);
   if(!checkNoErosion->isChecked())
   {
      QsGraph->attach(HPlot);
      CGraph->attach(HPlot);
   }

   // order determines order of display in Legend
   //VJ 101223 changed for qwt 6.0.0
   PGraph->setAxes(HPlot->xBottom, HPlot->yLeft);
   QGraph->setAxes(HPlot->xBottom, HPlot->yLeft);
   QsGraph->setAxes(HPlot->xBottom, HPlot->yRight);
   CGraph->setAxes(HPlot->xBottom, HPlot->yRight);

   QColor col;
   col.setRgb( 200,0,0,255 ); // darkred
   CGraph->setPen(QPen(col));
   QsGraph->setPen(QPen(Qt::red));
   col.setRgb( 60,100,160,255 );
   QGraph->setPen(QPen(col));
   PGraph->setPen(QPen("#000000"));

   PGraph->setStyle(QwtPlotCurve::Steps);
   QGraph->setStyle(QwtPlotCurve::Lines);
   QsGraph->setStyle(QwtPlotCurve::Lines);
   CGraph->setStyle(QwtPlotCurve::Lines);

   PGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
   QGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
   QsGraph->setRenderHint(QwtPlotItem::RenderAntialiased);
   CGraph->setRenderHint(QwtPlotItem::RenderAntialiased);

   QwtLegend *legend = new QwtLegend(HPlot);//this);//widgetGraph);
   legend->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
   HPlot->insertLegend(legend, QwtPlot::BottomLegend);
   //legend

   HPlot->setCanvasBackground(QBrush(Qt::white));

   HPlot->enableAxis(HPlot->yRight,true);
   HPlot->enableAxis(HPlot->yLeft,true);
   HPlot->enableAxis(HPlot->xBottom,true);
   HPlot->setAxisTitle(HPlot->xBottom, "time (min)");
   HPlot->setAxisTitle(HPlot->yLeft, "Q (l/s)/P (mm/h)");
   HPlot->setAxisTitle(HPlot->yRight, "Qs(kg/s)/C(g/l)");
   HPlot->setAxisScale(HPlot->yRight, 0, 1);
   HPlot->setAxisScale(HPlot->yLeft, 0, 100);
   HPlot->setAxisScale(HPlot->xBottom, 0, 100);

   // set axes

   QwtPlotGrid *grid = new QwtPlotGrid();
   grid->enableXMin(true);
   grid->enableYMin(true);
   col.setRgb( 180,180,180,180 );
   grid->setMajPen(QPen(col, 0, Qt::DashLine));
   col.setRgb( 210,210,210,180 );
   grid->setMinPen(QPen(col, 0 , Qt::DotLine));
   grid->attach(HPlot);
   // set gridlines


   HPlot->replot();
   // draw empty plot

   QData = NULL; //discharge
   QsData = NULL;  //sed discharge
   CData = NULL; //conc
   PData = NULL; //rainfall
   timeData = NULL;  //time
   // init data arrays for plot data

}
//---------------------------------------------------------------------------

void lisemqt::showPlot()
{
   // first time do this
   if (startplot)
   {
      startplot = false;

      stepP = 0;
      // VJ 110417 op.runstep did not work properly

      yas = 0.1;
      y2as = 0.1;

      // create the arrays for the curves in the first timestep when the total size is known
      timeData = new double[op.maxstep+2];
      QData = new double[op.maxstep+2];
      QsData = new double[op.maxstep+2];
      CData = new double[op.maxstep+2];
      PData = new double[op.maxstep+2];

      memset(timeData, 0, sizeof(timeData));
      memset(QData, 0, sizeof(QData));
      memset(QsData, 0, sizeof(QsData));
      memset(CData, 0, sizeof(CData));
      memset(PData, 0, sizeof(PData));

      HPlot->setAxisScale(HPlot->xBottom, op.BeginTime, op.EndTime);
   }

   stepP++;
   timeData[stepP] = op.time;
   PData[stepP] = op.P;
   QData[stepP] = op.Q;
   QsData[stepP] = op.Qs;
   CData[stepP] = op.C;

   //qwt 6.0.0:
   QGraph->setRawSamples(timeData,QData,stepP);
   PGraph->setRawSamples(timeData,PData,stepP);
   if(!checkNoErosion->isChecked())
   {
      //qwt 6.0.0:
      QsGraph->setRawSamples(timeData,QsData,stepP);
      CGraph->setRawSamples(timeData,CData,stepP);
   }

   y2as = max(y2as, op.Qs);
   y2as = max(y2as, op.C);
   HPlot->setAxisScale(HPlot->yRight, 0, y2as*1.05);
   yas = max(yas, op.Q);
   yas = max(yas, op.P);
   HPlot->setAxisScale(HPlot->yLeft, 0, yas*1.05);

   HPlot->replot();
   //   HPlot->canvas()->invalidatePaintCache();
   //   HPlot->canvas()->update(canvas()->contentsRect());

}
//---------------------------------------------------------------------------
