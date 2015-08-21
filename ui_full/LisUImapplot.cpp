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
  op.displayPcum = checkDisplayPcum->isChecked();

  if (radioButton_RO->isChecked())    op.drawMapType = 1;
  if (radioButton_INF->isChecked())   op.drawMapType = 2;
  if (radioButton_SL->isChecked())    op.drawMapType = 3;
  if (radioButton_FL->isChecked())    op.drawMapType = 4;
  if (radioButton_FLV->isChecked())   op.drawMapType = 5;
  if (radioButton_P->isChecked())     op.drawMapType = 6;
  if (radioButton_FEW->isChecked())   op.drawMapType = 7;

  showMap();
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
  radioButton_SL->setEnabled(!checkNoErosion->isChecked());
  radioButton_FL->setEnabled(checkChannelFlood->isChecked());
  radioButton_P->setChecked(false);

  //    transparency->setValue(128);  //main data
  //    transparency2->setValue(160); //channels
  //    transparency3->setValue(160); //roads
  //    transparency4->setValue(100); //houses

  // link with runfile
  //    op.drawMapType = MapDisplayMapSelection;
  //    transparency->setValue(std::min(MapDisplayHydrology, 200));  //main data
  //    transparency2->setValue(std::min(MapDisplayChannels, 200)); //channels
  //    transparency3->setValue(std::min(MapDisplayRoads, 200)); //roads
  //    transparency4->setValue(std::min(MapDisplayBuilding, 200)); //houses

  //    doubleSpinBoxRO->setValue(MapDisplayRunoffMax);
  //    doubleSpinBoxINF->setValue(MapDisplayInfiltrationMax);
  //    doubleSpinBoxSL->setValue(MapDisplaySoillossMax);
  //    doubleSpinBoxFL->setValue(MapDisplayFlooddepthMax);
  //    floodCutoffLevel->setValue(std::max(std::min(MapDisplayMinimumDepth,2),0));

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

  houseMap = new QwtPlotSpectrogram();
  houseMap->setRenderThreadCount( 0 );
  houseMap->attach( MPlot );
  // building structure map

  roadMap = new QwtPlotSpectrogram();
  roadMap->setRenderThreadCount( 0 );
  roadMap->attach( MPlot );
  // road map

  channelMap = new QwtPlotSpectrogram();
  channelMap->setRenderThreadCount( 0 );
  channelMap->attach( MPlot );
  // channel map

  baseMap = new QwtPlotSpectrogram();
  baseMap->setRenderThreadCount( 0 );
  baseMap->attach( MPlot );
  // shaded relief


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
//  mapRescaler->setReferenceAxis( QwtPlot::yLeft );  //NOT resets the plot after checking another map !!!
  mapRescaler->setAspectRatio( QwtPlot::xBottom, 1.0 );
  mapRescaler->setAspectRatio( QwtPlot::yRight, 0.0 );
  mapRescaler->setAspectRatio( QwtPlot::xTop, 0.0 );
//  mapRescaler->setRescalePolicy( QwtPlotRescaler::Fitting ); // also resets map to lower boundary
//  mapRescaler->setRescalePolicy( QwtPlotRescaler::Expanding ); //=DEFAULT ANYWAY
  mapRescaler->setExpandingDirection( QwtPlotRescaler::ExpandUp );
  //mapRescaler->setEnabled( true );
  // rescaling fixed to avoid deformation

  magnifier = new QwtPlotMagnifier( MPlot->canvas() );
  magnifier->setAxisEnabled( MPlot->yRight, false );
  magnifier->setZoomInKey((int)Qt::Key_Plus, Qt::ShiftModifier);
  // exclude right axis legend from rescaling

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
double lisemqt::fillDrawMapData(cTMap *_M, QwtMatrixRasterData *_RD, double type)
{
  double maxV = -1e20;
  mapData.clear();  //QVector double

  if (_M == NULL)
    return (maxV);

  // copy map data into vector for the display structure
  for(int r = _M->nrRows()-1; r >= 0; r--)
    //      for(int r = 0; r < _M->nrRows(); r++)
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

  mapData.replace(0, (double)type);
  // highjack position 0,0 with flag to get the variable unit in the cursor in trackerTextF

  // set intervals for rasterdata, x,y,z min and max
  _RD->setValueMatrix( mapData, _M->nrCols() );
  // set column number to divide vector into rows

  _RD->setInterval( Qt::XAxis, QwtInterval( 0, (double)_M->nrCols()*_M->cellSize(), QwtInterval::ExcludeMaximum ) );
  _RD->setInterval( Qt::YAxis, QwtInterval( 0, (double)_M->nrRows()*_M->cellSize(), QwtInterval::ExcludeMaximum ) );
  // set x/y axis intervals
  return maxV;
}
//---------------------------------------------------------------------------
// show the maps on screen
// the order of showing layers is determined by the order in how they are added to MPlot,
// not how they are done here!
void lisemqt::showMap()
{
  //drawMap->setAlpha(transparency->value());
  //drawMap->setAlpha(255);
  if (op.drawMapType == 1) showMap1();
  if (op.drawMapType == 2) showMap2();
  if (op.drawMapType == 3) showMap3();
  if (op.drawMapType == 4) showMap4();
  if (op.drawMapType == 5) showMap5();
  if (op.drawMapType == 6) showMap6();
  if (op.drawMapType == 7) showMap7();

  MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::showBaseMap()
{
  if (!startplot)
    return;

  double res = fillDrawMapData(op.baseMap, RDb, 0);
  if (res == -1e20)
    return;

  baseMap->setAlpha(75);
  baseMap->setColorMap(new colorMapGray());
  RDb->setInterval( Qt::ZAxis, QwtInterval( 0,res));
  baseMap->setData(RDb);
  // setdata sets a pointer to DRb to the private QWT d_data Qvector

  res = fillDrawMapData(op.baseMapDEM, RDbb, 7);
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
      double res = fillDrawMapData(op.channelMap, RDc,0 );
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
  //channelMap->setAlpha(transparency2->value());

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

  if (checkMapRoads->isChecked())
    roadMap->setAlpha(transparency3->value());
  else
    roadMap->setAlpha(0);

  //roadMap->setAlpha(transparency3->value());

  if (op.drawMapType == 2)
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
  if (checkMapBuildings->isChecked())
    houseMap->setAlpha(transparency4->value());
  else
    houseMap->setAlpha(0);

  houseMap->setColorMap(new colorMapHouse());
}
//---------------------------------------------------------------------------
// RUNOFF
void lisemqt::showMap1()
{
  MPlot->setTitle("Runoff (l/s)");

  QwtLinearColorMapVJ *pal1a = new colorMapWaterLog();
  QwtLinearColorMapVJ *pal1b = new colorMapWaterLog();

  double MaxV = fillDrawMapData(op.DrawMap1, RD, 1);
  if (MaxV ==-1e20)
    return;

  maxAxis1 = std::max(maxAxis1, MaxV);
  if (doubleSpinBoxRO->value() > 0)
    maxAxis1 = doubleSpinBoxRO->value();

  RD->setInterval( Qt::ZAxis, QwtInterval( 0, std::max(0.01, maxAxis1)));
  // classify data from 0 to,ax

  drawMap->setData(RD);

  pal1a->setThreshold(doubleSpinBoxROmin->value());
  pal1b->setThreshold(doubleSpinBoxROmin->value());

  drawMap->setColorMap(pal1a);
  // link data to map with a color palette

  rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), pal1b);
  // set the legend with the same palette

  if (maxAxis1 < 10)
    MPlot->setAxisScale( MPlot->yRight, 0.001, std::max(1.0,maxAxis1));
  else
    if (maxAxis1 < 100)
      {
      MPlot->setAxisScale( MPlot->yRight, 0.01, std::max(10.0,maxAxis1));
      doubleSpinBoxROmin->setValue(std::max(0.01,doubleSpinBoxROmin->value() ));
      }
    else
      if (maxAxis1 < 1000)
        {
          MPlot->setAxisScale( MPlot->yRight, 0.1, std::max(100.0,maxAxis1));
          doubleSpinBoxROmin->setValue(std::max(0.1,doubleSpinBoxROmin->value() ));
        }
      else
        {
          MPlot->setAxisScale( MPlot->yRight, 1, std::max(1000.0,maxAxis1));
          doubleSpinBoxROmin->setValue(std::max(1.0,doubleSpinBoxROmin->value() ));
        }

  MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLog10ScaleEngine() );
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
// INFILTRATION
void lisemqt::showMap2()
{
  MPlot->setTitle("Infiltration (mm)");

  pal2a = new colorMapWater();
  pal2b = new colorMapWater();

  // fill vector RD with matrix data and find the new max value
  double MaxV = fillDrawMapData(op.DrawMap2, RD, 2);
  if (MaxV ==-1e20)
    return;

  maxAxis2 = std::max(maxAxis2, MaxV);
  if (doubleSpinBoxINF->value() > 0)
    maxAxis2 = doubleSpinBoxINF->value();
  else
    maxAxis2 = MaxV;

  RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis2));

  drawMap->setData(RD);
  drawMap->setColorMap(pal2a);

  rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), pal2b);

  MPlot->setAxisScale( MPlot->yRight, 0, maxAxis2);
  MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
// draw a map, RD (QVector) and mapData (QwtPlotSpectrogram) are reused
// SOILLOSS
void lisemqt::showMap3()
{
  MPlot->setTitle("Soil loss (ton/ha)");

  pal3a = new colorMapSed();
  pal3b = new colorMapSed();

  double MaxV = fillDrawMapData(op.DrawMap3, RD, 3);
  if (MaxV ==-1e20)
    return;

  maxAxis3 = std::max(maxAxis3, MaxV);
  if (doubleSpinBoxSL->value() > 0)
    maxAxis3 = doubleSpinBoxSL->value();
  else
    maxAxis3 = MaxV;

  RD->setInterval( Qt::ZAxis, QwtInterval( -maxAxis3, maxAxis3));

  drawMap->setData(RD);
  drawMap->setColorMap(pal3a);

  rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), pal3b);

  MPlot->setAxisScale( MPlot->yRight, -maxAxis3, maxAxis3);
  MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
// FLOOD LEVEL
void lisemqt::showMap4()
{
  MPlot->setTitle("Flood level (m)");

  pal4a = new colorMapFlood();
  pal4b = new colorMapFlood();

  double MaxV = fillDrawMapData(op.DrawMap4, RD, 4);
  if (MaxV ==-1e20)
    return;

  maxAxis4 = std::max(maxAxis4, MaxV);
  if (doubleSpinBoxFL->value() > 0)
    maxAxis4 = doubleSpinBoxFL->value();
  else
    maxAxis4 = MaxV;

  RD->setInterval( Qt::ZAxis, QwtInterval( 0.000, maxAxis4));

  drawMap->setData(RD);

  pal4a->setThreshold(std::max(0.001,doubleSpinBoxFLmin->value()));
  pal4b->setThreshold(doubleSpinBoxFLmin->value());

  drawMap->setColorMap(pal4a);
  rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), pal4b);

  MPlot->setAxisScale( MPlot->yRight, 0, maxAxis4);
  MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
void lisemqt::showMap5()
{
  MPlot->setTitle("Flood Velocity (m/s)");

  pal5a = new colorMapFloodV();
  pal5b = new colorMapFloodV();

  double MaxV = fillDrawMapData(op.DrawMap5, RD, 5);
  if (MaxV ==-1e20)
    return;

  maxAxis5 = std::max(maxAxis5, MaxV);
  if (doubleSpinBoxFLV->value() > 0)
    maxAxis5 = doubleSpinBoxFLV->value();
  else
    maxAxis5 = MaxV;

  RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis5));

  drawMap->setData(RD);

  pal5a->setThreshold(std::max(doubleSpinBoxFLVmin->value(),0.001));
  pal5b->setThreshold(std::max(doubleSpinBoxFLVmin->value(),0.001));

  drawMap->setColorMap(pal5a);
  rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), pal5b);

  MPlot->setAxisScale( MPlot->yRight, 0, maxAxis5);
  MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
// PRECIPITATION
void lisemqt::showMap6()
{

  if (op.displayPcum)
    MPlot->setTitle("Cumulative Precipitation (mm)");
  else
    MPlot->setTitle("Precipitation (mm)");

  pal6a = new colorMapP();
  pal6b = new colorMapP();

  double MaxV = fillDrawMapData(op.DrawMap6, RD, 6);
  if (MaxV ==-1e20)
    return;
  // fill vector and find the new max value

  maxAxis6 = std::max(maxAxis6, MaxV);
  if (doubleSpinBoxP->value() > 0)
    maxAxis6 = doubleSpinBoxP->value();
  else
    maxAxis6 = MaxV;

  RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis6));

  drawMap->setData(RD);
  drawMap->setColorMap(pal6a);

  rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), pal6b);
  MPlot->setAxisScale( MPlot->yRight, 0, maxAxis6);
  MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );

}
//---------------------------------------------------------------------------
// EARLY WARNING
void lisemqt::showMap7()
{
  MPlot->setTitle("Moment of inundation after start of rainfall (min)");

  pal7a = new colorMapFEW();
  pal7b = new colorMapFEW();

  double MaxV = fillDrawMapData(op.DrawMap7, RD, 7);
  if (MaxV ==-1e20)
    return;
  // fill vector and find the new max value

  maxAxis7 = std::max(maxAxis7, MaxV);
  if (doubleSpinBoxFEW->value() > 0)
    maxAxis7 = doubleSpinBoxFEW->value();

  RD->setInterval( Qt::ZAxis, QwtInterval( 0, maxAxis7));

  drawMap->setData(RD);

  pal7a->setThreshold(0.001);
  pal7b->setThreshold(0.001);

  drawMap->setColorMap(pal7a);
  rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), pal7b);

  MPlot->setAxisScale( MPlot->yRight, 0, maxAxis7);
  MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
}
//---------------------------------------------------------------------------
