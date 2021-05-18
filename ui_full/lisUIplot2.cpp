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

#include "lisemqt.h"
#include "global.h"

//#include <QtCharts/QChartView>
//#include <QtCharts/QLineSeries>
//#include <QtCharts/QAreaSeries>

//QT_CHARTS_USE_NAMESPACE


void lisemqt::SetupnewPlot()
{
    QGraphN = new QLineSeries();
    PGraphN = new QLineSeries();

    QGraphN->setName("Discharge");
    QPen pen("#4444DD");
    pen.setWidth(3);
    QGraphN->setPen(pen);

    PGraphN->setName("Rain intensity");
    pen.setColor("#111111");
    pen.setWidth(3);
    PGraphN->setPen(pen);

    PQSchart = new QChart();

    PQSchart->addSeries(PGraphN);
    PQSchart->addSeries(QGraphN);

    if (checkDoErosion->isChecked()) {
        QsGraphN = new QLineSeries();
        QsGraphN->setName("Sediment discharge");
        pen.setColor("#FF1111");
        pen.setWidth(3);
        QsGraphN->setPen(pen);

        CGraphN = new QLineSeries();
        CGraphN->setName("Sediment concentration");
        pen.setColor("#cc0000");
        pen.setWidth(3);
        CGraphN->setPen(pen);

        PQSchart->addSeries(QsGraphN);
        PQSchart->addSeries(CGraphN);
    }

    PQSchart->setTitle("Simple example");

    axisX = new QValueAxis;
    axisX->setTickCount(10);
    axisX->setTitleText("Time (min)");

    axisYQ = new QValueAxis;
    axisYQ->setTickCount(10);
    axisYQ->setTitleText("Q (l/s)");

    axisYP = new QValueAxis;
    axisYP->setTickCount(10);
    axisYP->setTitleText("P (mm/h)");

    PQSchart->addAxis(axisX, Qt::AlignBottom);
    PQSchart->addAxis(axisYQ, Qt::AlignLeft);
    PQSchart->legend()->setAlignment(Qt::AlignBottom);

   if(checkDoErosion->isChecked()) {
        PQSchart->addAxis(axisYP, Qt::AlignLeft);
        axisYQs = new QValueAxis;
        axisYC = new QValueAxis;
        axisYQs->setTitleText("Q (kg)");
        axisYQs->setTitleText("Q (g/l)");
        QsGraphN->attachAxis(axisX);
        QsGraphN->attachAxis(axisYQs);
        CGraphN->attachAxis(axisX);
        CGraphN->attachAxis(axisYC);
        PQSchart->addAxis(axisYQs, Qt::AlignRight);
        PQSchart->addAxis(axisYC, Qt::AlignRight);
    }
    else
        PQSchart->addAxis(axisYP, Qt::AlignRight);

    PGraphN->attachAxis(axisX);
    PGraphN->attachAxis(axisYP);

    QGraphN->attachAxis(axisX);
    QGraphN->attachAxis(axisYQ);


    chartView = new QChartView(PQSchart);
    //chartView->setRenderHint(QPainter::Antialiasing);
    verticalLayout_2->insertWidget(0, chartView,1);
}


void lisemqt::newPlot()
{

      PQSchart->axes(Qt::Horizontal).first()->setRange(op.BeginTime,op.EndTime);

      int index = OutletIndices.indexOf(this->outletpoint);

    //  dataRain.clear();
    //  dataQ.clear();

   //   for(int i = 0; i < OutletQ[index]->length();i++)
    //  {
      int i = std::max(0,OutletQ[index]->length()-1);
          if(checkDoErosion->isChecked()) {
              QPointF _Qs;
              QPointF _C;
              _Qs.setX(TData.at(i));
              _Qs.setY(OutletQs.at(index)->at(i));
              _C.setX(TData.at(i));
              _C.setY(OutletC.at(index)->at(i));
              qsmax[index] = std::max(qsmax[index], OutletQs[index]->at(i));
              cmax[index] = std::max(cmax[index], OutletC[index]->at(i));
         }
          QPointF _P = QPointF(TData.at(i), Rainfall.at(i));
          QPointF _Q = QPointF(TData.at(i), OutletQ.at(index)->at(i));
//          dataRain << _P;
  //        dataQ << _Q;
          qmax[index] = std::max(qmax[index], OutletQ[index]->at(i));


      //PGraphN->replace(dataRain);
      //QGraphN->replace(dataQ);
          PGraphN->append(_P);
          QGraphN->append(_Q);
      axisYP->setRange(0.0,op.maxRainaxis);
      axisYQ->setRange(0.0,qmax[index]*1.1);

      if(checkDoErosion->isChecked()) {
          axisYQs->setMax(qsmax[index]);
          axisYC->setMax(cmax[index]);
      }

      chartView->update();

}
