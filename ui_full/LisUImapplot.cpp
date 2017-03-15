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

#include <algorithm>
#include "lisemqt.h"
#include "operation.h"
#include "global.h"
#include "LisUImapplot.h"

//---------------------------------------------------------------------------
QwtLinearColorMapVJ::QwtLinearColorMapVJ( const QColor &color1,const QColor &color2, QwtLinearColorMap::Format format ):
    QwtLinearColorMap( format )
{
    thresholdLCM = 0;
    QwtLinearColorMap::setColorInterval( color1, color2 );
}

QwtLinearColorMapVJ::~QwtLinearColorMapVJ()
{
}

void QwtLinearColorMapVJ::setThreshold( double value)
{
    thresholdLCM = value;
}

//---------------------------------------------------------------------------
void lisemqt::ssetAlpha(int v)
{
    //drawMap->setAlpha(v);
    baseMap->setAlpha(v);

    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaChannel(int v)
{
    channelMap->setAlpha(v);
    if (v > 0)
        MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaRoad(int v)
{
    roadMap->setAlpha(v);
    if (v > 0)
        MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaHouse(int v)
{
    houseMap->setAlpha(v);
    if (v > 0)
        MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaBarrier(int v)
{
  flowbarriersMap->setAlpha(v);
  if (v > 0)
    MPlot->replot();
}

//---------------------------------------------------------------------------
// called when a model run is started
void lisemqt::initMapPlot()
{
    maxAxis1 = -1e20;
    maxAxis2 = -1e20;
    maxAxis3 = -1e20;
    maxAxis4 = -1e20;
    maxAxis5 = -1e20;
    pstep = 0;

    //    transparency->setValue(128);  //main data
    //    transparency2->setValue(160); //channels
    //    transparency3->setValue(160); //roads
    //    transparency4->setValue(100); //houses
}

//---------------------------------------------------------------------------
// called at the start of openLisem, creates structures to hold maps
void lisemqt::setupMapPlot()
{
    title.setText("Runoff (l/s)");
    title.setFont(QFont("MS Shell Dlg 2",12));
    MPlot = new QwtPlot(title, this);
    // make the plot window
    //Layout_Map_2

  // verticalLayout_2->addWidget( MPlot, 0);
       maplayout->insertWidget(1, MPlot, 0, 0);
    // put it on screen
    MPlot->canvas()->setFrameStyle( QFrame::StyledPanel);
    MPlot->enableAxis( MPlot->yRight );
    MPlot->setAxisTitle(HPlot->xBottom, "m");
    MPlot->setAxisTitle(HPlot->yLeft, "m");

    // attach plot to widget in UI

    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->setPen( QPen( Qt::DotLine ) );
    grid->attach( MPlot );

    // NOTE the order in which these are attached is the order displayed.
    baseMapDEM = new QwtPlotSpectrogram();
    baseMapDEM->setRenderThreadCount( 0 );
    baseMapDEM->attach( MPlot );
    // dem

    drawMap = new QwtPlotSpectrogram();
    drawMap->setRenderThreadCount( 0 );
    drawMap->attach( MPlot );
    //map for runoff, infil, flood etc

    baseMap = new QwtPlotSpectrogram();
    baseMap->setRenderThreadCount( 0 );
    baseMap->attach( MPlot );
    // shaded relief   

    houseMap = new QwtPlotSpectrogram();
    houseMap->setRenderThreadCount( 0 );
    houseMap->attach( MPlot );
    // building structure map

    roadMap = new QwtPlotSpectrogram();
    roadMap->setRenderThreadCount( 0 );
    roadMap->attach( MPlot );
    // road map

    flowbarriersMap = new QwtPlotSpectrogram();
    flowbarriersMap->setRenderThreadCount( 0 );
    flowbarriersMap->attach( MPlot );

    channelMap = new QwtPlotSpectrogram();
    channelMap->setRenderThreadCount( 0 );
    channelMap->attach( MPlot );
    // channel map

    RD = new QwtMatrixRasterData();
    RDb = new QwtMatrixRasterData();
    RDbb = new QwtMatrixRasterData();
    RDc = new QwtMatrixRasterData();
    RDd = new QwtMatrixRasterData();
    RDe = new QwtMatrixRasterData();
    RDf = new QwtMatrixRasterData();

    // raster data to link to plot

    rightAxis = new QwtScaleWidget();
    rightAxis = MPlot->axisWidget( MPlot->yRight );
    rightAxis->setColorBarEnabled( true );
    rightAxis->setColorBarWidth( 20 );
    // legend to the right of the plot

    mapRescaler = new QwtPlotRescaler( MPlot->canvas() );
    //  mapRescaler->setReferenceAxis( QwtPlot::yLeft );
    //NOT resets the plot all the time after checking another map !!!
    mapRescaler->setAspectRatio( QwtPlot::xBottom, 1.0 );
//    mapRescaler->setAspectRatio( QwtPlot::yLeft, 1.0 );
    mapRescaler->setAspectRatio( QwtPlot::yRight, 0.0 );
    mapRescaler->setAspectRatio( QwtPlot::xTop, 0.0 );
//      mapRescaler->setRescalePolicy( QwtPlotRescaler::Fitting ); // every tmestep fits map to lower boundary, position not maintained
    //  mapRescaler->setRescalePolicy( QwtPlotRescaler::Expanding ); //=DEFAULT ANYWAY
     mapRescaler->setExpandingDirection( QwtPlotRescaler::ExpandUp );
    //mapRescaler->setEnabled( true );
    // rescaling fixed to avoid deformation

    magnifier = new QwtPlotMagnifier( MPlot->canvas() );
    magnifier->setAxisEnabled( MPlot->yRight, false );
    // exclude right axis legend from rescaling
    magnifier->setZoomInKey((int)Qt::Key_Plus, Qt::ShiftModifier);

    panner = new QwtPlotPanner( MPlot->canvas() );
    panner->setAxisEnabled( MPlot->yRight, false );
    // exclude right axis legend from panning

    picker = new MyPicker( MPlot->canvas() );

    maxAxis1 = -1e20;
    maxAxis2 = -1e20;
    maxAxis3 = -1e20;
    maxAxis4 = -1e20;
    pstep = 0;

}
//---------------------------------------------------------------------------
// fill the current raster data structure with new data, called each run step
double lisemqt::fillDrawMapData(cTMap *_M, QwtMatrixRasterData *_RD)//, double type)
{
    double maxV = -1e20;
    mapData.clear();  //QVector double

    if (_M == NULL)
        return (maxV);

    // copy map data into vector for the display structure
    for(int r = _M->nrRows()-1; r >= 0; r--)
        for(int c=0; c < _M->nrCols(); c++)
        {
            if(!pcr::isMV(_M->Drc))
            {
                mapData << _M->Drc;
                maxV = std::max(maxV, _M->Drc);
            }
            else
                mapData << (double)-1e20;
        }

  //  mapData.replace(0, (double)type);
    // highjack position 0,0 with flag to get the variable unit in the cursor in trackerTextF
    //obsolete

    // set intervals for rasterdata, x,y,z min and max
    _RD->setValueMatrix( mapData, _M->nrCols() );
    // set column number to divide vector into rows

    _RD->setInterval( Qt::XAxis, QwtInterval( 0, (double)_M->nrCols()*_M->cellSize(), QwtInterval::ExcludeMaximum ) );
    _RD->setInterval( Qt::YAxis, QwtInterval( 0, (double)_M->nrRows()*_M->cellSize(), QwtInterval::ExcludeMaximum ) );
    // set x/y axis intervals
    return maxV;
}
//---------------------------------------------------------------------------
void lisemqt::showMapb(bool)
{
    showMap();
}
void lisemqt::showMapd(double)
{
    showMap();
}
//---------------------------------------------------------------------------
// show the maps on screen
// the order of showing layers is determined by the order in how they are added to MPlot,
// not how they are done here!
void lisemqt::showMap()
{
    if(op.comboboxset == false)
    {
        op.comboboxset = true;
        for(int i = ColorMapList.length() - 1; i >-1 ; i--)
        {
            delete ColorMapList.at(i);

        }
        ColorMapList.clear();
        DisplayComboBox->clear();
        DisplayComboBox2->clear();
        NameList.clear();
        UnitList.clear();
        SymList.clear();
        LogList.clear();
        picker->NameList.clear();
        picker->UnitList.clear();
        IndexList.clear();
        IndexList1.clear();

        for(int i = 0; i < op.ComboMapsSafe.length(); i++)
        {
            QwtComboColorMap *cm = new QwtComboColorMap(QColor(op.ComboColors.at(i).at(0)),QColor(op.ComboColors.at(i).at(op.ComboColors.at(i).length()-1)),op.ComboColorMap.at(i),op.ComboColors.at(i));

            ColorMapList.append(cm);
            NameList.append(op.ComboMapNames.at(i));
            UnitList.append(op.ComboUnits.at(i));
            SymList.append(op.ComboSymColor.at(i)); // symetric colors
            LogList.append(op.ComboLogaritmic.at(i));  //log display
            picker->NameList.append(op.ComboMapNames.at(i));
            picker->UnitList.append(op.ComboUnits.at(i));
        }
        QStringList S;
        QStringList S1;

        for(int i = 0; i < op.ComboMapsSafe.length(); i++)
        {

            if(op.ComboLists.at(i) == 0)
            {
                S << QString(op.ComboMapNames.at(i) + " (" + op.ComboUnits.at(i) + ")");
                IndexList.append(i);
            }else
            {
                S1 << QString(op.ComboMapNames.at(i) + " (" + op.ComboUnits.at(i) + ")");
                IndexList1.append(i);
            }
        }
        DisplayComboBox->addItems(S);
        DisplayComboBox2->addItems(S1);
        ActiveList = 0;

        DisplayComboBox->setFixedWidth(180);
        DisplayComboBox2->setFixedWidth(180);

        DisplayComboBox2->setCurrentIndex(0);
        DisplayComboBox->setCurrentIndex(0);

        DisplayComboBox->setMaxVisibleItems(IndexList.count());
        DisplayComboBox2->setMaxVisibleItems(IndexList1.count());
    }

    // needed if user clicks while nothing is running:
    if (IndexList.count() == 0)
        return;

    if (ActiveList == 0)
    {
        showComboMap(IndexList.at(DisplayComboBox->currentIndex()));
    }
    else
        if (ActiveList == 1)
    {
        showComboMap(IndexList1.at(DisplayComboBox2->currentIndex()));
    }

    channelMap->setAlpha(checkMapChannels->isChecked() ? transparency2->value() : 0);
    roadMap->setAlpha(checkMapRoads->isChecked() ? transparency3->value() : 0);
    houseMap->setAlpha(checkMapBuildings->isChecked() ? transparency4->value() : 0);
    flowbarriersMap->setAlpha(checkMapFlowBarriers->isChecked() ? transparency5->value() : 0);

    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::showComboMap(int i)
{
    if( i < 0 || i >= op.ComboMapsSafe.length())
    {
        return;
    }

    MPlot->setTitle(op.ComboMapNames.at(i) + " (" + op.ComboUnits.at(i) + ")");

    // fill vector RD with matrix data and find the new max value
    double MaxV = fillDrawMapData(op.ComboMapsSafe.at(i), RD);//, i);
    if (MaxV ==-1e20)
        return;

    double MinV = mapMinimum(*op.ComboMapsSafe.at(i));

    // set stepsize
    if(op.ComboLists.at(i) == 0)
    {
        ComboMinSpinBox->setSingleStep(op.comboStep.at(i));
        ComboMaxSpinBox->setSingleStep(op.comboStep.at(i));
    }
    else
    {
        ComboMinSpinBox2->setSingleStep(op.comboStep.at(i));
        ComboMaxSpinBox2->setSingleStep(op.comboStep.at(i));
    }

    // get spinbox values, can be 0
    double mi = op.userMinV.at(i);
    double ma = op.userMaxV.at(i);

    if (ma == 0)
        ma = MaxV; // use map max when ma = 0
    if (mi == 0)
        mi = MinV;
    if(mi == ma) // because of they are equal nothing is displayed (raincum)
        mi = 0;
//qDebug() << mi << ma << MinV << MaxV;

    if (op.ComboSymColor.at(i)) // symetric coloring for soilloss
        mi = -ma;

    RD->setInterval( Qt::ZAxis, QwtInterval( mi, ma));

    QwtComboColorMap *cm = new QwtComboColorMap(QColor(op.ComboColors.at(i).at(0)),
                                                QColor(op.ComboColors.at(i).at(op.ComboColors.at(i).length()-1)),
                                                op.ComboColorMap.at(i),op.ComboColors.at(i));
    QwtComboColorMap *cmL = new QwtComboColorMap(QColor(op.ComboColors.at(i).at(0)),
                                                 QColor(op.ComboColors.at(i).at(op.ComboColors.at(i).length()-1)),
                                                 op.ComboColorMap.at(i),op.ComboColors.at(i));

    cm->thresholduse = ActiveList == 0;
    cmL->thresholduse = false; // !op.ComboSymColor.at(i);

    cm->thresholdmin = mi;
    cmL->thresholdmin = mi;

    drawMap->setData(RD);
    drawMap->setColorMap(cm);

    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), cmL);

    if(op.ComboLogaritmic.at(i))
    {
        int coef = int(log10(ma));
        mi = (mi == 0 ? std::pow(0.1,3-coef) : mi);
        ma = (ma == 0 ? std::pow(10,coef)    : ma);

        MPlot->setAxisScale( MPlot->yRight, mi, ma );
        MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLog10ScaleEngine() );
    }
    else
    {

        MPlot->setAxisScale( MPlot->yRight, mi, ma);//std::max(mi,0.001)
        MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
    }
}
//---------------------------------------------------------------------------
void lisemqt::showBaseMap()
{
    if (!startplot)
        return;

    double res = fillDrawMapData(op.baseMap, RDb);//, 0);
    if (res == -1e20)
        return;

    baseMap->setAlpha(75);
    baseMap->setColorMap(new colorMapGray());
    RDb->setInterval( Qt::ZAxis, QwtInterval( 0,res));
    baseMap->setData(RDb);
    // setdata sets a pointer to DRb to the private QWT d_data Qvector

    res = fillDrawMapData(op.baseMapDEM, RDbb);//, 7);
    if (res == -1e20)
        return;
    double mindem = mapMinimum(*op.baseMapDEM);

    baseMapDEM->setAlpha(255);
    baseMapDEM->setColorMap(new colorMapElevation());//colorMapGray());//
    RDbb->setInterval( Qt::ZAxis, QwtInterval( mindem,res));
    baseMapDEM->setData(RDbb);

    double nrCols = (double)op.baseMap->nrCols()*op.baseMap->cellSize();
    double nrRows = (double)op.baseMap->nrRows()*op.baseMap->cellSize();
    double dx = std::max(nrCols,nrRows)/20;
    // reset the axes to the correct rows/cols,
    // do only once because resets zooming and panning

    //  MPlot->setAxisAutoScale(MPlot->yRight, false);
    //  MPlot->setAxisAutoScale(MPlot->xBottom, true);
    //  MPlot->setAxisAutoScale(MPlot->yLeft, false);

    MPlot->setAxisScale( MPlot->xBottom, 0.0, nrCols, dx);
    MPlot->setAxisMaxMinor( MPlot->xBottom, 0 );
    MPlot->setAxisScale( MPlot->yLeft, 0.0, nrRows, dx);
    MPlot->setAxisMaxMinor( MPlot->yLeft, 0 );

}
//---------------------------------------------------------------------------
void lisemqt::showChannelMap()
{
    if (startplot)
    {
        double res = fillDrawMapData(op.channelMap, RDc);//,0 );
        if (res ==-1e20)
            return;
        QwtLinearColorMapVJ *pala = new colorMapFlood();
        pala->thresholdLCM = 0.01;

        channelMap->setColorMap(pala);
        RDc->setInterval( Qt::ZAxis, QwtInterval( 0,1.0));
        channelMap->setData(RDc);
    }

    if (checkMapChannels->isChecked())
        channelMap->setAlpha(transparency2->value());
    else
        channelMap->setAlpha(0);
}
//---------------------------------------------------------------------------
void lisemqt::showRoadMap()
{
    if (startplot)
    {
        double res = fillDrawMapData(op.roadMap, RDd);//, 0);
        if (res ==-1e20)
            return;
        RDd->setInterval( Qt::ZAxis, QwtInterval( 0,0.5));
        roadMap->setData(RDd);
    }

    if (checkMapRoads->isChecked())
        roadMap->setAlpha(transparency3->value());
    else
        roadMap->setAlpha(0);

    roadMap->setColorMap(new colorMapRoads3());
}
//---------------------------------------------------------------------------
void lisemqt::showHouseMap()
{
    if (startplot)
    {
        // set intervals for rasterdata, x,y,z min and max
        double res = fillDrawMapData(op.houseMap, RDe);//, 0);
        if (res ==-1e20)
            return;
        RDe->setInterval( Qt::ZAxis, QwtInterval( 0.0 ,res));
        houseMap->setData(RDe);
    }
    if (checkMapBuildings->isChecked())
        houseMap->setAlpha(transparency4->value());
    else
        houseMap->setAlpha(0);

    houseMap->setColorMap(new colorMapHouse());
}
//---------------------------------------------------------------------------
void lisemqt::showFlowBarriersMap()
{
  if (startplot)
    {

      // set intervals for rasterdata, x,y,z min and max
      double res = fillDrawMapData(op.flowbarriersMap, RDf);//, 0);
      if (res ==-1e20)
        return;
      RDf->setInterval( Qt::ZAxis, QwtInterval( 0.0, res));
      flowbarriersMap->setData(RDf);
    }

  if (checkMapFlowBarriers->isChecked())
    flowbarriersMap->setAlpha(transparency5->value());
  else
    flowbarriersMap->setAlpha(0);

  flowbarriersMap->setColorMap(new colorMapFlowBarrier());
}

