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

#include "lisemqt.h"
#include "global.h"
#include "LisUImapplot.h"



//---------------------------------------------------------------------------
void lisemqt::ssetAlpha(int v)
{
    drawMap->setAlpha(v);
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlpha2(int v)
{
    channelMap->setAlpha(v);
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlpha3(int v)
{
    roadMap->setAlpha(v);
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlpha4(int v)
{
    houseMap->setAlpha(v);
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::selectMapType(bool doit)
{
    if (stopplot)
        return;

    if (radioButton_RO->isChecked())    op.drawMapType = 1;
    if (radioButton_INF->isChecked())   op.drawMapType = 2;
    if (radioButton_SL->isChecked())    op.drawMapType = 3;
    if (radioButton_FL->isChecked())    op.drawMapType = 4;

    showMap(); // show map

    //    showChannelMap(); // show channel map

    //    showRoadMap(); // show road map

}
//---------------------------------------------------------------------------
// called when a model run is started
void lisemqt::initMapPlot()
{
    op.drawMapType = 1;
    radioButton_RO->setChecked(true);
    radioButton_INF->setChecked(false);
    if (checkNoErosion->isChecked())
        radioButton_SL->setEnabled(false);
    else
        radioButton_SL->setEnabled(true);
    if (checkChannelFlood->isChecked())
        radioButton_FL->setEnabled(true);
    else
        radioButton_FL->setEnabled(false);


    maxAxis1 = -1e20;
    maxAxis2 = -1e20;
    maxAxis3 = -1e20;
    maxAxis4 = -1e20;
    pstep = 0;
    transparency->setValue(180);  //main data
    transparency2->setValue(128); //channels
    transparency3->setValue(180); //roads
    transparency4->setValue(100); //houses
    // slider setting basemap transparency

}

//---------------------------------------------------------------------------
// called at the start of openLisem, creates structures to hold maps
void lisemqt::setupMapPlot()
{
    op.drawMapType = 1;

    title.setText("Runoff (l/s)");
    title.setFont(QFont("MS Shell Dlg 2",12));
    MPlot = new QwtPlot(title, this);
    // make the plot window
    Layout_Map->insertWidget(0, MPlot, 1);
    // put it on screen
    //MPlot->canvas()->setBorderRadius( 0 );
    MPlot->canvas()->setFrameStyle( QFrame::StyledPanel);
    MPlot->enableAxis( MPlot->yRight );
    MPlot->setAxisTitle(HPlot->xBottom, "m");
    MPlot->setAxisTitle(HPlot->yLeft, "m");

    // attach plot to widget in UI

    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->setPen( QPen( Qt::DotLine ) );
    grid->attach( MPlot );


    baseMap = new QwtPlotSpectrogram();
    baseMap->setRenderThreadCount( 0 );
    baseMap->attach( MPlot );
    // shaded relief base map

    drawMap = new QwtPlotSpectrogram();
    drawMap->setRenderThreadCount( 0 );
    drawMap->attach( MPlot );
    drawMap->setAlpha(180);
    // NOTE the order in which these are attached is the order displayed.

    houseMap = new QwtPlotSpectrogram();
    houseMap->setRenderThreadCount( 0 );
    houseMap->attach( MPlot );
    // building structure map

    roadMap = new QwtPlotSpectrogram();
    roadMap->setRenderThreadCount( 0 );
    roadMap->attach( MPlot );
    // channel map

    channelMap = new QwtPlotSpectrogram();
    channelMap->setRenderThreadCount( 0 );
    channelMap->attach( MPlot );
    // channel map

    RD = new QwtMatrixRasterData();
    RDb = new QwtMatrixRasterData();
    RDc = new QwtMatrixRasterData();
    RDd = new QwtMatrixRasterData();
    RDe = new QwtMatrixRasterData();

    // raster data to link to plot

    rightAxis = new QwtScaleWidget();
    rightAxis = MPlot->axisWidget( MPlot->yRight );
    rightAxis->setColorBarEnabled( true );
    rightAxis->setColorBarWidth( 20 );
    // legend to the right of the plot

    mapRescaler = new QwtPlotRescaler( MPlot->canvas() );
    mapRescaler->setReferenceAxis( QwtPlot::yLeft );
    mapRescaler->setAspectRatio( QwtPlot::xBottom, 1.0 );
    mapRescaler->setAspectRatio( QwtPlot::yRight, 0.0 );
    mapRescaler->setAspectRatio( QwtPlot::xTop, 0.0 );
    mapRescaler->setRescalePolicy( QwtPlotRescaler::Fitting );
    mapRescaler->setExpandingDirection( QwtPlotRescaler::ExpandUp );
    mapRescaler->setEnabled( true );
    // rescaling fixed to avoid deformation

    magnifier = new QwtPlotMagnifier( MPlot->canvas() );
    magnifier->setAxisEnabled( MPlot->yRight, false );

    panner = new QwtPlotPanner( MPlot->canvas() );
    panner->setAxisEnabled( MPlot->yRight, false );

    picker = new MyPicker( MPlot->canvas() );  

    maxAxis1 = -1e20;
    maxAxis2 = -1e20;
    maxAxis3 = -1e20;
    maxAxis4 = -1e20;
    pstep = 0;

}
//---------------------------------------------------------------------------
// fill the current raster data structure with new data, called each run step
double lisemqt::fillDrawMapData(TMMap *_M, QwtMatrixRasterData *_RD)
{
    double maxV = -1e20;
    mapData.clear();  //QVector double

    // copy map data into vector for the display structure
    for(int r = _M->nrRows-1; r >= 0; r--)
        //      for(int r = 0; r < _M->nrRows; r++)
        for(int c=0; c < _M->nrCols; c++)
        {
            if(!IS_MV_REAL8(&_M->Drc))
            {
                mapData << _M->Drc;
                maxV = qMax(maxV, _M->Drc);
            }
            else
                mapData << (double)-1e20;
        }

    // set intervals for rasterdata, x,y,z min and max
    _RD->setValueMatrix( mapData, _M->nrCols );
    // set column number to divide vector into rows
    _RD->setInterval( Qt::XAxis, QwtInterval( 0, (double)_M->nrCols*_M->MH.cellSizeX, QwtInterval::ExcludeMaximum ) );
    _RD->setInterval( Qt::YAxis, QwtInterval( 0, (double)_M->nrRows*_M->MH.cellSizeY, QwtInterval::ExcludeMaximum ) );
    // set x/y axis intervals
    return maxV;
}
//---------------------------------------------------------------------------
// show the maps on screen
// the order of showing layers is determined by the order in how they are added to MPlot,
// not how they are done here!
void lisemqt::showMap()
{
    drawMap->setAlpha(transparency->value());
    if (op.drawMapType == 1) showMap1();
    if (op.drawMapType == 2) showMap2();
    if (op.drawMapType == 3) showMap3();
    if (op.drawMapType == 4) showMap4();

    MPlot->replot();
    // do not do resets panning
    //   mapRescaler->rescale();
}
//---------------------------------------------------------------------------
void lisemqt::showBaseMap()
{
    if (!startplot)
        return;

    double res = fillDrawMapData(op.baseMap, RDb);

    baseMap->setAlpha(255);
    baseMap->setColorMap(new colorMapGray());
    RDb->setInterval( Qt::ZAxis, QwtInterval( 0,res));
    baseMap->setData(RDb);
    // setdata sets a pointer to DRb to the private QWT d_data Qvector

    double nrCols = (double)op.baseMap->nrCols*op.baseMap->MH.cellSizeX;
    double nrRows = (double)op.baseMap->nrRows*op.baseMap->MH.cellSizeY;

    // reset the axes to the correct rows/cols,
    // do only once because resets zooming and panning
    MPlot->setAxisScale( MPlot->xBottom, 0.0, nrCols, nrCols/20);
    MPlot->setAxisMaxMinor( MPlot->xBottom, 0 );
    MPlot->setAxisScale( MPlot->yLeft, 0.0, nrRows, nrRows/20);
    MPlot->setAxisMaxMinor( MPlot->yLeft, 0 );

   // startplot = false;
}
//---------------------------------------------------------------------------
void lisemqt::showChannelMap()
{
    if (!startplot)
        return;
    double res = fillDrawMapData(op.channelMap, RDc);

    channelMap->setAlpha(transparency2->value());

    channelMap->setColorMap(new colorMapFlood());
    RDc->setInterval( Qt::ZAxis, QwtInterval( 0,0.5));
    channelMap->setData(RDc);
}
//---------------------------------------------------------------------------
void lisemqt::showRoadMap()
{
//    if (!startplot)
//        return;
    double res = fillDrawMapData(op.roadMap, RDd);

    roadMap->setAlpha(transparency3->value());

    if (op.drawMapType <= 2)
        roadMap->setColorMap(new colorMapRoads());
    else
        roadMap->setColorMap(new colorMapRoads2());

    RDd->setInterval( Qt::ZAxis, QwtInterval( 0,0.5));
    roadMap->setData(RDd);
}
//---------------------------------------------------------------------------
void lisemqt::showHouseMap()
{
    if (!startplot)
        return;

    // set intervals for rasterdata, x,y,z min and max
    double res = fillDrawMapData(op.houseMap, RDe);

    houseMap->setAlpha(transparency4->value());

    houseMap->setColorMap(new colorMapHouse());

    RDe->setInterval( Qt::ZAxis, QwtInterval( 0.0 ,res));
    houseMap->setData(RDe);
}
//---------------------------------------------------------------------------
void lisemqt::showHouseMapA()
{
    if (!startplot)
        return;

    double maxV = -1;
    QVector<double> mapDataa;

    TMMap *M = new TMMap();
    M->PathName = W->inputDir + "lu5m.map";
    bool res = M->LoadFromFile();
    if (!res)
    {
        ErrorString = "Cannot find map " +M->PathName;
        throw 1;
    }

    mapDataa.clear();

    // copy map data into vector for the display structure
    for(int r = M->nrRows-1; r >= 0; r--)
        for(int c=0; c < M->nrCols; c++)
        {
            if(!IS_MV_REAL8(&M->Drc))
            {
                mapDataa << M->Drc;
                maxV = qMax(maxV, M->Drc);
            }
            else
                mapDataa << (double)-1e20;
        }

    // set intervals for rasterdata, x,y,z min and max
    RDe->setValueMatrix( mapDataa, M->nrCols );

    houseMap->setAlpha(transparency4->value());

    houseMap->setColorMap(new colorMapGray());

    RDe->setInterval( Qt::ZAxis, QwtInterval( 0.0 ,maxV));
    houseMap->setData(RDe);
    MPlot->setAxisScale( MPlot->xBottom, 0.0, M->nrCols, M->nrCols/20);
    MPlot->setAxisScale( MPlot->yLeft, 0.0, M->nrRows, M->nrRows/20);

}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
void lisemqt::showMap1()
{
    MPlot->setTitle("Runoff (l/s)");

    double MaxV = fillDrawMapData(op.DrawMap1, RD);
    MaxV = fillDrawMapData(op.DrawMap1, RD);
    // fill vector and find the new max value

    // set intervals for rasterdata, x,y,z min and max
    maxAxis1 = qMax(maxAxis1, MaxV);
    if (doubleSpinBoxRO->value() > 0)
        maxAxis1 = doubleSpinBoxRO->value();
    /*    else
        maxAxis1 = MaxV*/;
    RD->setInterval( Qt::ZAxis, QwtInterval( 0, qMax(0.1, maxAxis1)));

    drawMap->setData(RD);
    drawMap->setColorMap(new colorMapWaterLog());

    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapWaterLog());

    if (maxAxis1 < 10)
        MPlot->setAxisScale( MPlot->yRight, 0.01, qMax(1.0,maxAxis1));
    else
        if (maxAxis1 < 100)
            MPlot->setAxisScale( MPlot->yRight, 0.1, qMax(10.0,maxAxis1));
        else
            MPlot->setAxisScale( MPlot->yRight, 0.1, qMax(100.0,maxAxis1));

    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLog10ScaleEngine() );
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
void lisemqt::showMap2()
{
    MPlot->setTitle("Infiltration (mm)");

    // fill vector RD with matrix data and find the new max value
    double MaxV = fillDrawMapData(op.DrawMap2, RD);
    MaxV = fillDrawMapData(op.DrawMap2, RD);

    // set the new interval to the new max value
    maxAxis2 = qMax(maxAxis2, MaxV);
    if (doubleSpinBoxINF->value() > 0)
        maxAxis2 = doubleSpinBoxINF->value();
    else
        maxAxis2 = MaxV;
    RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis2));

    // point spectrogram to data
    drawMap->setData(RD);
    drawMap->setColorMap(new colorMapWater());
    //drawMap = QwtPlotSpectrogram

    // set the right axis legend to the new interval
    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapWater());
    MPlot->setAxisScale( MPlot->yRight, 0, maxAxis2);
    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
void lisemqt::showMap3()
{
    MPlot->setTitle("Soil loss (ton/ha)");

    double MaxV = fillDrawMapData(op.DrawMap3, RD);
    MaxV = fillDrawMapData(op.DrawMap3, RD);
    // fill vector and find the new max value

    maxAxis3 = qMax(maxAxis3, MaxV);
    if (doubleSpinBoxSL->value() > 0)
        maxAxis3 = doubleSpinBoxSL->value();
    else
        maxAxis3 = MaxV;
    RD->setInterval( Qt::ZAxis, QwtInterval( -maxAxis3, maxAxis3));

    drawMap->setData(RD);
    drawMap->setColorMap(new colorMapSed());
    //QwtPlotSpectrogram

    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapSed());
    MPlot->setAxisScale( MPlot->yRight, -maxAxis3, maxAxis3);
    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
void lisemqt::showMap4()
{
    MPlot->setTitle("Flood level (m)");

    double MinV = 0;
    double MaxV = fillDrawMapData(op.DrawMap4, RD);
    MaxV = fillDrawMapData(op.DrawMap4, RD);
    // fill vector and find the new max value

    maxAxis4 = qMax(maxAxis4, MaxV);
    if (doubleSpinBoxFL->value() > 0)
        maxAxis4 = doubleSpinBoxFL->value();
    else
        maxAxis4 = MaxV;
    //    if (doubleSpinBoxFLmin->value() > 0)
    //        MinV = doubleSpinBoxFLmin->value();

    RD->setInterval( Qt::ZAxis, QwtInterval( MinV, maxAxis4));

    drawMap->setData(RD);
    if (floodCutoffLevel->value() == 0) drawMap->setColorMap(new colorMapFlood());
    if (floodCutoffLevel->value() == 1) drawMap->setColorMap(new colorMapFlood005());
    if (floodCutoffLevel->value() == 2) drawMap->setColorMap(new colorMapFlood01());
    // draw map

    if (floodCutoffLevel->value() == 0) rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood());
    if (floodCutoffLevel->value() == 1) rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood005());
    if (floodCutoffLevel->value() == 2) rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood01());
    MPlot->setAxisScale( MPlot->yRight, MinV, maxAxis4);
    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
    // draw legend
}
//---------------------------------------------------------------------------
