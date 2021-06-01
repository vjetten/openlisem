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

void lisemqt::initNewPlot()
{
    QGraphN->clear();
    PGraphN->clear();
    if (checkDoErosion->isChecked()) {
        QsGraphN->clear();
        CGraphN->clear();
    }
}


void lisemqt::setupNewPlot()
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

    PQSchart->setTitle("Hydrograph");

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
    chartView->chart()->setTheme(QChart::ChartThemeDark);
    //chartView->setRenderHint(QPainter::Antialiasing); very slow
    layout_Plot->insertWidget(0, chartView,1);

}


void lisemqt::newPlot(bool refresh)
{    
    int index = OutletIndices.indexOf(this->outletpoint);

    if (refresh) {
        // fill everything with a new index point and replace entire series
        dataRain.clear();
        dataQ.clear();
        qmax.clear();
        qsmax.clear();
        cmax.clear();
        qmax.append(0);
        qsmax.append(0);
        cmax.append(0);

        for(int i = 0; i < op.OutletQ[index]->length()-1;i++)
        {
            QPointF _P = QPointF(op.Time.at(i), op.Pmm.at(i));
            QPointF _Q = QPointF(op.Time.at(i), op.OutletQ[index]->at(i));
            qmax[index] = std::max(qmax[index], op.OutletQ[index]->at(i));
            dataRain << _P;
            dataQ << _Q;

            if(checkDoErosion->isChecked()) {
                QPointF _Qs = QPointF(op.Time.at(i),op.OutletQs.at(index)->at(i));
                QPointF _C = QPointF(op.Time.at(i),op.OutletC.at(index)->at(i));
                qsmax[index] = std::max(qsmax[index], op.OutletQs[index]->at(i));
                cmax[index] = std::max(cmax[index], op.OutletC[index]->at(i));
                dataQs << _Qs;
                dataC << _C;
            }
        //qDebug() << _P << _Q << i << OutletQ[index]->length();
        }
        PGraphN->replace(dataRain);
        QGraphN->replace(dataQ);


        if(checkDoErosion->isChecked()) {
            QsGraphN->replace(dataQs);
            CGraphN->replace(dataC);
        }
    } else {
        //add to the last point

          int i = std::max(0,op.OutletQ[index]->length()-1);

          QPointF _P = QPointF(op.Time.at(i), op.Pmm.at(i));
          QPointF _Q = QPointF(op.Time.at(i), op.OutletQ.at(index)->at(i));
          PGraphN->append(_P);
          QGraphN->append(_Q);

//          qmax[index] = std::max(qmax[index], op.OutletQ[index]->at(i));
//          if(checkDoErosion->isChecked()) {
//              QPointF _Qs = QPointF(op.Time.at(i),op.OutletQs.at(index)->at(i));
//              QPointF _C = QPointF(op.Time.at(i),op.OutletC.at(index)->at(i));
//              qsmax[index] = std::max(qsmax[index], op.OutletQs[index]->at(i));
//              cmax[index] = std::max(cmax[index], op.OutletC[index]->at(i));

//              QsGraphN->append(_Qs);
//              CGraphN->append(_C);
//          }

//          double step = floor(log10(qmax[index]))+1.0;
//          axisYQ->setTickCount(int(step)*2);
//          axisYQ->setRange(0.0,pow(10.0,step));//qmax[index]*1.1);

//          step = floor(log10(op.maxRainaxis))+1.0;
//          axisYP->setTickCount(int(step)*2);
//          axisYP->setRange(0.0,pow(10.0,step));

//          if(checkDoErosion->isChecked()) {
//              axisYQs->setMax(qsmax[index]);
//              axisYC->setMax(cmax[index]);
//          }
      }
    PQSchart->axes(Qt::Horizontal).first()->setRange(op.BeginTime,op.EndTime);
    chartView->update();
}
