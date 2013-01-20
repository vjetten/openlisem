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
#include "LisUImapcolor.h"

//---------------------------------------------------------------------------
void lisemqt::ssetAlpha(int v)
{
    baseMap->setAlpha(v);
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::selectMapType(bool doit)
{
    if (radioButton_RO->isChecked())    op.drawMapType = 1;
    if (radioButton_INF->isChecked())   op.drawMapType = 2;
    if (radioButton_SL->isChecked())    op.drawMapType = 3;
    if (radioButton_FL->isChecked())    op.drawMapType = 4;
}
//---------------------------------------------------------------------------
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
    transparency->setValue(80);

}

//---------------------------------------------------------------------------
void lisemqt::setupMapPlot()
{
    op.drawMapType = 1;

    title.setText("Runoff (l/s)");
    title.setFont(QFont("MS Shell Dlg 2",12));
    MPlot = new QwtPlot(title, this);
    // make the plot window
    Layout_Map->insertWidget(0, MPlot, 1);
    //MPlot->canvas()->setBorderRadius( 0 );
    MPlot->canvas()->setFrameStyle( QFrame::StyledPanel);

    // attach plot to widget in UI

    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->setPen( QPen( Qt::DotLine ) );
    grid->attach( MPlot );

    drawMap = new QwtPlotSpectrogram();
    drawMap->setRenderThreadCount( 0 );
    drawMap->attach( MPlot );
    // make the raster drawing

    baseMap = new QwtPlotSpectrogram();
    baseMap->setRenderThreadCount( 0 );
    baseMap->attach( MPlot );
    // shaded relief base map

    RD = new QwtMatrixRasterData();
    RDb = new QwtMatrixRasterData();

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
    // rescaling fixed to avoid deformation


    magnifier = new QwtPlotMagnifier( MPlot->canvas() );
    magnifier->setAxisEnabled( MPlot->yRight, false );

    panner = new QwtPlotPanner( MPlot->canvas() );
    panner->setAxisEnabled( MPlot->yRight, false );

    //   (void) new QwtPlotPanner( MPlot->canvas() );
    //   (void) new QwtPlotMagnifier( MPlot->canvas() );

    maxAxis1 = -1e20;
    maxAxis2 = -1e20;
    maxAxis3 = -1e20;
    maxAxis4 = -1e20;
    pstep = 0;

}
//---------------------------------------------------------------------------
void lisemqt::fillDrawMapData(TMMap *_M, QwtMatrixRasterData *_RD)
{
    mapData.clear();
    // copy map data into vector for the display structure
    for(int r = _M->nrRows-1; r >= 0; r--)
        //      for(int r = 0; r < _M->nrRows; r++)
        for(int c=0; c < _M->nrCols; c++)
        {
            if(!IS_MV_REAL8(&_M->Data[r][c]))
            {
                mapData << _M->Data[r][c];
            }
            else
                mapData << (double)-1e20;
        }

    // set intervals for rasterdata, x,y,z min and max
    _RD->setValueMatrix( mapData, _M->nrCols );
    // set column number to divide vector into rows
    _RD->setInterval( Qt::XAxis, QwtInterval( 0, (double)_M->nrCols, QwtInterval::ExcludeMaximum ) );
    _RD->setInterval( Qt::YAxis, QwtInterval( 0, (double)_M->nrRows, QwtInterval::ExcludeMaximum ) );
    // set x/y axis intervals
}

//---------------------------------------------------------------------------
void lisemqt::fillDrawMapBaseData(TMMap *_M, TMMap *_M2)
{
    mapData.clear();
    // copy map data into vector for the display structure
    for(int r = _M->nrRows-1; r >= 0; r--)
        for(int c=0; c < _M->nrCols; c++)
        {
            if(!IS_MV_REAL8(&_M->Data[r][c]))
            {
                if (_M2->Data[r][c] > 0)
                    mapData << (double)-1e20;
                else
                    mapData << _M->Data[r][c];

            }
            else
                mapData << (double)-1e20;
        }

    // set intervals for rasterdata, x,y,z min and max
    RDb->setValueMatrix( mapData, _M->nrCols );
    // set column number to divide vector into rows
    RDb->setInterval( Qt::XAxis, QwtInterval( 0, (double)_M->nrCols, QwtInterval::ExcludeMaximum ) );
    RDb->setInterval( Qt::YAxis, QwtInterval( 0, (double)_M->nrRows, QwtInterval::ExcludeMaximum ) );
    // set x/y axis intervals
}

//---------------------------------------------------------------------------
// the order of showing layers is determined by the order in how they are added to MPlot,
// not how they are done here!
void lisemqt::ShowMap()
{
    double MinV = 0.1;
    fillDrawMapData(op.DrawMap, RD);

    op.DrawMap->ResetMinMax();
    double MaxV = (double)op.DrawMap->MH.maxVal;
    double nrCols = (double)op.DrawMap->nrCols;
    double nrRows = (double)op.DrawMap->nrRows;

    if (op.drawMapType == 1)
    {
        MPlot->setTitle("Runoff (l/s)");
        drawMap->setColorMap(new colorMapWaterLog());
        maxAxis1 = qMax(maxAxis1, MaxV);
        if (doubleSpinBoxRO->value() > 0)
            maxAxis1 = doubleSpinBoxRO->value();
        RD->setInterval( Qt::ZAxis, QwtInterval( 0, qMax(MinV, maxAxis1)));
    }
    else
        if (op.drawMapType == 2)
        {
            MPlot->setTitle("Infiltration (mm)");
            drawMap->setColorMap(new colorMapWater());
            maxAxis2 = qMax(maxAxis2, MaxV);
            if (doubleSpinBoxINF->value() > 0)
                maxAxis2 = doubleSpinBoxINF->value();
            RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis2));
        }
        else
            if (op.drawMapType == 3)
            {
                MPlot->setTitle("Soil loss (ton/ha)");
                drawMap->setColorMap(new colorMapSedB());
                maxAxis3 = qMax(maxAxis3, MaxV);
                if (doubleSpinBoxSL->value() > 0)
                    maxAxis3 = doubleSpinBoxSL->value();
                RD->setInterval( Qt::ZAxis, QwtInterval( -maxAxis3, maxAxis3));
                // use max and -max for sediemnt delivery so that white legend color is in the middle, no activity
                // cyan is deposition and dark orange is erosion
            }
            else
                if (op.drawMapType == 4)
                {
                    MPlot->setTitle("Flood level (m)");
                    drawMap->setColorMap(new colorMapFlood());
                    maxAxis4 = qMax(maxAxis4, MaxV);
                    if (doubleSpinBoxFL->value() > 0)
                        maxAxis4 = doubleSpinBoxFL->value();
                    RD->setInterval( Qt::ZAxis, QwtInterval( 0.001, maxAxis4));
                }

    drawMap->setData(RD);
    // link raster data to drawMap

    // add legend right of axis
    if (op.drawMapType == 1)
    {
        // log scale for runoff
        rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapWaterLog());
        if (maxAxis1 < 100)
            MPlot->setAxisScale( MPlot->yRight, MinV, qMax(1.0,maxAxis1));
        else
            MPlot->setAxisScale( MPlot->yRight, MinV, qMax(10.0,maxAxis1));
        MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLog10ScaleEngine() );
    }
    else
        if (op.drawMapType == 2)
        {
            // lin scale for infiltration
            rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapWater());
            MPlot->setAxisScale( MPlot->yRight, 0.0, maxAxis2);
            MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
        }
        else
            if (op.drawMapType == 3)
            {
                //lin scale with mirrored max and -max for soil loss
                rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapSedB());
                MPlot->setAxisScale( MPlot->yRight, -maxAxis3, maxAxis3);
                MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
            }
            else
                if (op.drawMapType == 4)
                {
                    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), new colorMapFlood());
                    MPlot->setAxisScale( MPlot->yRight, 0.001, maxAxis4);
                    MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
                }
    MPlot->enableAxis( MPlot->yRight );

    MPlot->plotLayout()->setAlignCanvasToScales( true );

    MPlot->setAxisScale( MPlot->xBottom, 0.0, nrCols, nrCols/20);
    MPlot->setAxisMaxMinor( MPlot->xBottom, 0 );
    MPlot->setAxisScale( MPlot->yLeft, 0.0, nrRows, nrRows/20);
    MPlot->setAxisMaxMinor( MPlot->yLeft, 0 );

    mapRescaler->setEnabled( true );
    mapRescaler->setExpandingDirection( QwtPlotRescaler::ExpandUp );

//      ShowBaseMap();
//    drawMap->setAlpha(transparency->value());

    MPlot->replot();
    mapRescaler->rescale();

}
//---------------------------------------------------------------------------
void lisemqt::ShowBaseMap()
{
    if (!startplot)
        return;

    fillDrawMapData(op.baseMap, RDb);
    //fillDrawMapBaseData(op.baseMap, op.DrawMap);

    baseMap->setColorMap(new colorMapGray());
    RDb->setInterval( Qt::ZAxis, QwtInterval( 0,1.0));
    baseMap->setData(RDb);

    baseMap->setAlpha(transparency->value());

    startplot = false;
}
//---------------------------------------------------------------------------
