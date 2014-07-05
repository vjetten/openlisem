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

/*
void lisemqt::on_checkBoxOverlay_stateChanged(int yes)
{
    if (yes == Qt::Unchecked)
    {
//        ssetAlpha2(0);
//        ssetAlpha3(0);
//        ssetAlpha4(0);
    //    channelMap->setAlpha(0);
        houseMap->setAlpha(0);
  //      roadMap->setAlpha(0);
 //       MPlot->replot();
    }
    else
    {
  //      ssetAlpha2(transparency2->sliderPosition ());
   //     ssetAlpha3(transparency3->sliderPosition ());
   //     ssetAlpha4(transparency4->sliderPosition ());

    }
}
*/
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
    if (v > 0)
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlpha3(int v)
{
    roadMap->setAlpha(v);
    if (v > 0)
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlpha4(int v)
{
    houseMap->setAlpha(v);
    if (v > 0)
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::selectMapType(bool /* doit */)
{
    //    if (!startplot)
    //        return;

    op.addWHtohmx = checkAddWHtohmx->isChecked();
    op.displayPcum = checkDisplayPcum->isChecked();

    if (radioButton_RO->isChecked())    op.drawMapType = 1;
    if (radioButton_INF->isChecked())   op.drawMapType = 2;
    if (radioButton_SL->isChecked())    op.drawMapType = 3;
    if (radioButton_FL->isChecked())    op.drawMapType = 4;
    if (radioButton_FLV->isChecked())   op.drawMapType = 5;
    if (radioButton_P->isChecked())     op.drawMapType = 6;

    showMap(); // show map

    //    showChannelMap(); // show channel map

    //    showRoadMap(); // show road map

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
    radioButton_P->setChecked(false);

    transparency->setValue(180);  //main data
    transparency2->setValue(180); //channels
    transparency3->setValue(180); //roads
    transparency4->setValue(100); //houses

    // link with runfile
    //    op.drawMapType = MapDisplayMapSelection;
    //    transparency->setValue(qMin(MapDisplayHydrology, 200));  //main data
    //    transparency2->setValue(qMin(MapDisplayChannels, 200)); //channels
    //    transparency3->setValue(qMin(MapDisplayRoads, 200)); //roads
    //    transparency4->setValue(qMin(MapDisplayBuilding, 200)); //houses

    //    doubleSpinBoxRO->setValue(MapDisplayRunoffMax);
    //    doubleSpinBoxINF->setValue(MapDisplayInfiltrationMax);
    //    doubleSpinBoxSL->setValue(MapDisplaySoillossMax);
    //    doubleSpinBoxFL->setValue(MapDisplayFlooddepthMax);
    //    op.addWHtohmx = MapDisplayIncludeRunoff == 1;
    //    floodCutoffLevel->setValue(qMax(qMin(MapDisplayMinimumDepth,2),0));

    // slider setting basemap transparency
    //    if (p1.compare("Building display")==0)  MapDisplayBuilding = iii;
    //    if (p1.compare("Roads display")==0)     MapDisplayRoads = iii;
    //    if (p1.compare("Channels display")==0)  MapDisplayChannels = iii;
    //    if (p1.compare("Hydrology display")==0) MapDisplayHydrology = iii;
    //    if (p1.compare("Runoff max")==0)        MapDisplayRunoffMax = val;
    //    if (p1.compare("Infiltration max")==0)  MapDisplayInfiltrationMax = val;
    //    if (p1.compare("Soilloss max")==0)      MapDisplaySoillossMax = val;
    //    if (p1.compare("Flooddepth max")==0)    MapDisplayFlooddepthMax = val;
    //    if (p1.compare("Include runoff")==0)    MapDisplayIncludeRunoff = iii;
    //    if (p1.compare("Minimum depth")==0)     MapDisplayMinimumDepth = iii;
    //    if (p1.compare("Screendumps")==0)       MapDisplayScreenDumps = iii;
}
//---------------------------------------------------------------------------
// called at the start of openLisem, creates structures to hold maps
void lisemqt::setupMapPlot()
{
    op.drawMapType = 1;
    //   double alpha = 1;

    title.setText("Runoff (l/s)");
    title.setFont(QFont("MS Shell Dlg 2",12));
    MPlot = new QwtPlot(title, this);
    // make the plot window
    Layout_Map_2->insertWidget(0, MPlot, 1);
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

    baseMapDEM = new QwtPlotSpectrogram();
    baseMapDEM->setRenderThreadCount( 0 );
    baseMapDEM->attach( MPlot );
    // DEM information
    RD = new QwtMatrixRasterData();
    RDb = new QwtMatrixRasterData();
    RDbb = new QwtMatrixRasterData();
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
double lisemqt::fillDrawMapData(TMMap *_M, QwtMatrixRasterData *_RD, double type)
{
    double maxV = -1e20;
    mapData.clear();  //QVector double

    if (_M == NULL)
        return (maxV);

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

    mapData.replace(0, (double)type);
    // highjack position 0,0 with flag to get the variable unit in the cursor in trackerTextF

    // set intervals for rasterdata, x,y,z min and max
    _RD->setValueMatrix( mapData, _M->nrCols );
    // set column number to divide vector into rows

    _RD->setInterval( Qt::XAxis, QwtInterval( 0, (double)_M->nrCols*_M->MH.cellSize, QwtInterval::ExcludeMaximum ) );
    _RD->setInterval( Qt::YAxis, QwtInterval( 0, (double)_M->nrRows*_M->MH.cellSize, QwtInterval::ExcludeMaximum ) );
    // set x/y axis intervals
    return maxV;
}
//---------------------------------------------------------------------------
// show the maps on screen
// the order of showing layers is determined by the order in how they are added to MPlot,
// not how they are done here!
void lisemqt::showMap()
{
    //        if (!startplot)
    //            return;

    drawMap->setAlpha(transparency->value());
    if (op.drawMapType == 1) showMap1();
    if (op.drawMapType == 2) showMap2();
    if (op.drawMapType == 3) showMap3();
    if (op.drawMapType == 4) showMap4();
    if (op.drawMapType == 5) showMap5();
    if (op.drawMapType == 6) showMap6();

    MPlot->replot();
    // do not do resets panning
    //   mapRescaler->rescale();
}
//---------------------------------------------------------------------------
void lisemqt::showBaseMap()
{
    if (!startplot)
        return;

    double res = fillDrawMapData(op.baseMap, RDb, 0);
    if (res == -1e20)
        return;

    baseMap->setAlpha(255);
    baseMap->setColorMap(new colorMapGray());
    RDb->setInterval( Qt::ZAxis, QwtInterval( 0,res));
    baseMap->setData(RDb);
    // setdata sets a pointer to DRb to the private QWT d_data Qvector

    res = fillDrawMapData(op.baseMapDEM, RDbb, 7);
    if (res == -1e20)
        return;
    double mindem = op.baseMapDEM->mapMinimum();
    baseMapDEM->setAlpha(0);
    baseMapDEM->setColorMap(new colorMapGray());
    RDbb->setInterval( Qt::ZAxis, QwtInterval( mindem,res));
    baseMapDEM->setData(RDbb);
    double nrCols = (double)op.baseMap->nrCols*op.baseMap->MH.cellSize;
    double nrRows = (double)op.baseMap->nrRows*op.baseMap->MH.cellSize;
    double dx = max(nrCols,nrRows)/20;
    // reset the axes to the correct rows/cols,
    // do only once because resets zooming and panning

    MPlot->setAxisAutoScale(MPlot->yRight, false);
    MPlot->setAxisAutoScale(MPlot->xBottom, false);
    MPlot->setAxisAutoScale(MPlot->yLeft, false);

    MPlot->setAxisScale( MPlot->xBottom, 0.0, nrCols, dx);
    MPlot->setAxisMaxMinor( MPlot->xBottom, 0 );
    MPlot->setAxisScale( MPlot->yLeft, 0.0, nrRows, dx);
    MPlot->setAxisMaxMinor( MPlot->yLeft, 0 );

    /*
    if (nrRows < nrCols)
    {
        MPlot->setAxisScale( MPlot->xBottom, 0.0, nrCols, nrCols/20);
        MPlot->setAxisMaxMinor( MPlot->xBottom, 0 );
        MPlot->setAxisScale( MPlot->yLeft, 0.0, nrCols, nrCols/20);
        MPlot->setAxisMaxMinor( MPlot->yLeft, 0 );
    }
    else
    {
        MPlot->setAxisScale( MPlot->xBottom, 0.0, nrRows, nrRows/20);
        MPlot->setAxisMaxMinor( MPlot->xBottom, 0 );
        MPlot->setAxisScale( MPlot->yLeft, 0.0, nrRows, nrRows/20);
        MPlot->setAxisMaxMinor( MPlot->yLeft, 0 );
    }
    */
    // startplot = false;
}
//---------------------------------------------------------------------------
void lisemqt::showChannelMap()
{
    if (!startplot)
        return;

    double res = fillDrawMapData(op.channelMap, RDc,0 );
    if (res ==-1e20)
        return;

    channelMap->setAlpha(transparency2->value());

    channelMap->setColorMap(new colorMapFlood());
    RDc->setInterval( Qt::ZAxis, QwtInterval( 0,0.5));

    channelMap->setData(RDc);
}
//---------------------------------------------------------------------------
void lisemqt::showRoadMap()
{
    if (startplot)
    {
        double res = fillDrawMapData(op.roadMap, RDd, 0);
        if (res ==-1e20)
            return;
        RDd->setInterval( Qt::ZAxis, QwtInterval( 0,0.5));
        roadMap->setData(RDd);
    }


    roadMap->setAlpha(transparency3->value());

    if (op.drawMapType <= 2)
        roadMap->setColorMap(new colorMapRoads());
    else
        roadMap->setColorMap(new colorMapRoads2());

}
//---------------------------------------------------------------------------
void lisemqt::showHouseMap()
{
    if (startplot)
    {
        // set intervals for rasterdata, x,y,z min and max
        double res = fillDrawMapData(op.houseMap, RDe, 0);
        if (res ==-1e20)
            return;
        RDe->setInterval( Qt::ZAxis, QwtInterval( 0.0 ,res));
        houseMap->setData(RDe);
    }

    houseMap->setAlpha(transparency4->value());

    //    if (op.drawMapType == 4)
    houseMap->setColorMap(new colorMapHouse());
    //    else
    //        houseMap->setColorMap(new colorMapWhite());
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
void lisemqt::showMap1()
{
    MPlot->setTitle("Runoff (l/s)");

    double MaxV = fillDrawMapData(op.DrawMap1, RD, 1);
    if (MaxV ==-1e20)
        return;
    // fill vector and find the new max value

    // set intervals for rasterdata, x,y,z min and max
    maxAxis1 = qMax(maxAxis1, MaxV);
    if (doubleSpinBoxRO->value() > 0)
        maxAxis1 = doubleSpinBoxRO->value();
    // use userdefined if spinbox not 0

    RD->setInterval( Qt::ZAxis, QwtInterval( 0, qMax(0.1, maxAxis1)));
    // classify data from 0 to,ax

    drawMap->setData(RD);
    drawMap->setColorMap(new colorMapWaterLog());
    // link data to map with a color palette

    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapWaterLog());
    // set the legend with the same palette

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
// INFILTRATION
void lisemqt::showMap2()
{
    MPlot->setTitle("Infiltration (mm)");

    // fill vector RD with matrix data and find the new max value
    double MaxV = fillDrawMapData(op.DrawMap2, RD, 2);
    if (MaxV ==-1e20)
        return;
    // set the new interval to the new max value
    maxAxis2 = qMax(maxAxis2, MaxV);
    if (doubleSpinBoxINF->value() > 0)
        maxAxis2 = doubleSpinBoxINF->value();
    else
        maxAxis2 = MaxV;
    RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis2));

    drawMap->setData(RD);
    // point spectrogram to data
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

    double MaxV = fillDrawMapData(op.DrawMap3, RD, 3);
    if (MaxV ==-1e20)
        return;
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
    double MaxV = fillDrawMapData(op.DrawMap4, RD, 4);
    if (MaxV ==-1e20)
        return;
    // fill vector and find the new max value

    maxAxis4 = qMax(maxAxis4, MaxV);
    if (doubleSpinBoxFL->value() > 0)
        maxAxis4 = doubleSpinBoxFL->value();
    //    else
    //        maxAxis4 = MaxV;
    //    if (doubleSpinBoxFLmin->value() > 0)
    //        MinV = doubleSpinBoxFLmin->value();

    RD->setInterval( Qt::ZAxis, QwtInterval( MinV, maxAxis4));

    drawMap->setData(RD);
    if (checkFloodCutoff->isChecked())
        drawMap->setColorMap(new colorMapFlood005());
    else
        drawMap->setColorMap(new colorMapFlood());

//    if (floodCutoffLevel->value() == 0) drawMap->setColorMap(new colorMapFlood());
//    if (floodCutoffLevel->value() == 1) drawMap->setColorMap(new colorMapFlood005());
//    if (floodCutoffLevel->value() == 2) drawMap->setColorMap(new colorMapFlood01());
    // draw map
    if (checkFloodCutoff->isChecked())
        rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood005());
    else
        rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood());
    //    if (floodCutoffLevel->value() == 0) rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood());
    //    if (floodCutoffLevel->value() == 1) rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood005());
    //    if (floodCutoffLevel->value() == 2) rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood01());
    MPlot->setAxisScale( MPlot->yRight, MinV, maxAxis4);
    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
    // draw legend
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
void lisemqt::showMap5()
{
    MPlot->setTitle("Flood Velocity (m/s)");

    double MaxV = fillDrawMapData(op.DrawMap5, RD, 5);
    if (MaxV ==-1e20)
        return;
    // fill vector and find the new max value

    maxAxis5 = qMax(maxAxis5, MaxV);
    if (doubleSpinBoxSL->value() > 0)
        maxAxis5 = doubleSpinBoxFLV->value();
    else
        maxAxis5 = MaxV;
    RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis5));

    drawMap->setData(RD);
    drawMap->setColorMap(new colorMapFloodV());
    //QwtPlotSpectrogram

    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFloodV());
    MPlot->setAxisScale( MPlot->yRight, 0, maxAxis5);
    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
void lisemqt::showMap6()
{
    if (op.displayPcum)
        MPlot->setTitle("Cumulative Precipitation (mm)");
    else
        MPlot->setTitle("Precipitation (mm)");

    double MaxV = fillDrawMapData(op.DrawMap6, RD, 6);
    if (MaxV ==-1e20)
        return;
    // fill vector and find the new max value

    maxAxis6 = qMax(maxAxis6, MaxV);
    if (doubleSpinBoxP->value() > 0)
        maxAxis6 = doubleSpinBoxP->value();
    else
        maxAxis6 = MaxV;
    RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis6));

    drawMap->setData(RD);
    drawMap->setColorMap(new colorMapP());
    //QwtPlotSpectrogram

    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapP());
    MPlot->setAxisScale( MPlot->yRight, 0, maxAxis6);
    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}//---------------------------------------------------------------------------
