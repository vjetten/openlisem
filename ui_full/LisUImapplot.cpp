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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#include <algorithm>
#include "lisemqt.h"
#include "operation.h"
#include "global.h"
#include "LisUImapplot.h"

//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+_dy[ldd]==rTo && cFrom+_dx[ldd]==cTo )

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
void lisemqt::onImageToggled(bool b)
{
    baseMapImage->setAlpha(b? 255 : 0);
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlpha(int v)
{
    baseMap->setAlpha(v);

    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaMap(int v)
{
    drawMap->setAlpha(v);
    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaChannelOutlet(int v)
{
   // showChannelVector(false);
    showChannelVector(true);
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaChannel(int v)
{
   // showChannelVector(false);
    showChannelVector(true);
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaRoad(int v)
{
    roadMap->setAlpha(v);
    if (v > 0 && checkMapRoads->isChecked())
        MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaHouse(int v)
{
    houseMap->setAlpha(v);
    if (v > 0 && checkMapBuildings->isChecked())
        MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaHardSurface(int v)
{
    bool doit = (checkHouses->isChecked() || checkHardsurface->isChecked() || checkRoadsystem->isChecked());
    if (v > 0 && checkHouses->isChecked())
        houseMap->setAlpha(v);
    if (v > 0 && checkHardsurface->isChecked())
        hardsurfMap->setAlpha(v);
    if (v > 0 && checkRoadsystem->isChecked())
        roadMap->setAlpha(v);
    if (v > 0 && doit)
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

    // maps structures do not have to be deleted. every run has the same maps, and the content is filled ynamically
    // the rivers and culverts etc are really recreated and need to be deleted
    showChannelVector(false);

    rivers.clear();
    culverts.clear();
    //outlets.clear();
    //obspoints.clear();

}

void lisemqt::changeSize()
{
    double h = MPlot->height();
    double w = MPlot->width();

    MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrRows*op._dx*w/h);//, op._dx*10);
    MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrRows*op._dx);//, op._dx*10);

    MPlot->replot();

// from mapedit:
//    double h = MPlot->height();
//    double w = MPlot->width();

//    if(_nrCols >= _nrRows ) {
//        MPlot->setAxisScale( MPlot->xBottom, _llx, _llx+_nrCols*_dx*w/h, _dx*10);
//        MPlot->setAxisScale( MPlot->yLeft, _lly, _lly+_nrCols*_dx, _dx*10);
//    } else {
//        MPlot->setAxisScale( MPlot->xBottom, _llx, _llx+_nrRows*_dx*w/h,_dx*10);
//        MPlot->setAxisScale( MPlot->yLeft, _lly, _lly+_nrRows*_dx, _dx*10);
//    }
//    MPlot->replot();



}

//---------------------------------------------------------------------------
// called at the start of openLisem, creates structures to hold maps
void lisemqt::setupMapPlot()
{
    title.setText("Runoff (l/s)");
    title.setFont(QFont("MS Shell Dlg 2",12));

    MPlot = new QwtPlot(title, this);
    // make the plot window
    maplayout->insertWidget(1, MPlot, 0);

  //  MPlot->setStyleSheet(QString("* { background-color: %1 }").arg("#555555"));
    // put it on screen
    MPlot->enableAxis( MPlot->yRight );
    MPlot->setAxisTitle(MPlot->xBottom, "m");
    MPlot->setAxisTitle(MPlot->yLeft, "m");
    MPlot->setAxisLabelRotation(MPlot->yLeft, 270);
    MPlot->setAxisLabelAlignment(MPlot->yLeft, Qt::AlignVCenter);

    // attach plot to widget in UI

    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->setPen( QPen( Qt::DotLine ) );
    grid->attach( MPlot );

    // NOTE the order in which these are attached is the order displayed.
    // 0
    baseMapDEM = new QwtPlotSpectrogram();
    baseMapDEM->setRenderThreadCount( 0 );
    baseMapDEM->attach( MPlot );
    // dem
    // 1
    baseMapImage = new QwtPlotSpectrogram();
    baseMapImage->setRenderThreadCount( 0 );
    baseMapImage->attach( MPlot );
    //image

    // 2
    baseMap = new QwtPlotSpectrogram();
    baseMap->setRenderThreadCount( 0);
    baseMap->attach( MPlot );
    // shaded relief
    // 3
    roadMap = new QwtPlotSpectrogram();
    roadMap->setRenderThreadCount( 0 );
    roadMap->attach( MPlot );
    // road map

    // 4
    hardsurfMap = new QwtPlotSpectrogram();
    hardsurfMap->setRenderThreadCount( 0 );
    hardsurfMap->attach( MPlot );

    //5 data
    drawMap = new QwtPlotSpectrogram();
    drawMap->setRenderThreadCount( 0 );
    drawMap->attach( MPlot );
    //map for runoff, infil, flood etc

    // // 5
    // roadMap = new QwtPlotSpectrogram();
    // roadMap->setRenderThreadCount( 0 );
    // roadMap->attach( MPlot );
    // // road map

    // // 6
    // hardsurfMap = new QwtPlotSpectrogram();
    // hardsurfMap->setRenderThreadCount( 0 );
    // hardsurfMap->attach( MPlot );

    // 4
    houseMap = new QwtPlotSpectrogram();
    houseMap->setRenderThreadCount( 0 );
    houseMap->attach( MPlot );
    // building structure map

    //7
    outletMap = new QwtPlotSpectrogram();
    outletMap->setRenderThreadCount( 0 );
    outletMap->attach( MPlot );
    // outlet map used for outlet number when hovering (?)

    //8
    contourDEM = new QwtPlotSpectrogram();
    contourDEM->setRenderThreadCount( 0 );
    contourDEM->attach( MPlot );
    // contours

    RD = new QwtMatrixRasterData();
    RDb = new QwtMatrixRasterData();
    RDbb = new QwtMatrixRasterData();
    RDc = new QwtMatrixRasterData();
    RDd = new QwtMatrixRasterData();
    RDe = new QwtMatrixRasterData();
    RDf = new QwtMatrixRasterData();
    RImage = new QwtMatrixRasterData();

    // raster data to link to plot

    rightAxis = new QwtScaleWidget();
    rightAxis = MPlot->axisWidget( MPlot->yRight );
    rightAxis->setColorBarEnabled( true );
    rightAxis->setColorBarWidth( 20 );
    // legend to the right of the plot

    magnifier = new QwtPlotMagnifier( MPlot->canvas() );
    magnifier->setAxisEnabled( MPlot->yRight, false );
    // exclude right axis legend from rescaling
    magnifier->setZoomInKey(Qt::Key_Plus, Qt::ShiftModifier);
    magnifier->setZoomOutKey(Qt::Key_Minus, Qt::NoModifier );
    magnifier->setZoomInKey(Qt::Key_Plus, Qt::KeypadModifier);
    magnifier->setZoomOutKey(Qt::Key_Minus, Qt::KeypadModifier);

//    zoomer = new QwtPlotZoomer( MPlot->canvas() );
//    zoomer->setKeyPattern( QwtEventPattern::KeyRedo, Qt::Key_I, Qt::ShiftModifier );
//    zoomer->setKeyPattern( QwtEventPattern::KeyUndo, Qt::Key_O, Qt::ShiftModifier );
//    zoomer->setKeyPattern( QwtEventPattern::KeyHome, Qt::Key_Home );

    panner = new QwtPlotPanner( MPlot->canvas() );
    panner->setAxisEnabled( MPlot->yRight, false );
    // exclude right axis legend from panning

    picker = new MyPicker( (QwtPlotCanvas *) MPlot->canvas() );
    picker->setEnabled(true);

    mapRescaler = new QwtPlotRescaler( MPlot->canvas() );
 //   mapRescaler->setReferenceAxis( QwtPlot::xBottom );
    mapRescaler->setAspectRatio( QwtPlot::xBottom, 1.0 );
    mapRescaler->setAspectRatio( QwtPlot::yLeft, 1.0 );
    mapRescaler->setAspectRatio( QwtPlot::yRight, 0.0 );
    mapRescaler->setAspectRatio( QwtPlot::xTop, 1.0 );
    mapRescaler->setExpandingDirection( QwtPlotRescaler::ExpandUp );
    //mapRescaler->setRescalePolicy( QwtPlotRescaler::Fixed );

    MPlot->replot();

    maxAxis1 = -1e20;
    maxAxis2 = -1e20;
    maxAxis3 = -1e20;
    maxAxis4 = -1e20;
    pstep = 0;

}
//---------------------------------------------------------------------------
// fill the current raster data structure with new data, called each run step
double lisemqt::fillDrawMapData(cTMap *_M, double scale, QwtMatrixRasterData *_RD, double *minv, double *maxv){
    double maxV = -1e20;
    double minV = 1e20;
    mapData.clear();  //QVector double
    double sum = 0;

    if (_M == nullptr)
        return (maxV);

    // copy map data into vector for the display structure

    for(int r = _M->nrRows()-1; r >= 0; r--)
        for(int c=0; c < _M->nrCols(); c++)
        {
            if(!pcr::isMV(_M->Drc))
            {
                double v =_M->Drc*scale;
                mapData << v;
                maxV = std::max(maxV, v);
                minV = std::min(minV, v);
                sum += v;
            }
            else
                mapData << (double)-1e20;
        }

    // set intervals for rasterdata, x,y,z min and max
    _RD->setValueMatrix( mapData, _M->nrCols() );
    // set column number to divide vector into rows
    _RD->setInterval( Qt::XAxis, QwtInterval( op._llx,op._llx+(double)op._nrCols*op._dx, QwtInterval::ExcludeMaximum ) );
    _RD->setInterval( Qt::YAxis, QwtInterval( op._lly,op._lly+(double)op._nrRows*op._dx, QwtInterval::ExcludeMaximum ) );
    // set x/y axis intervals
    //qDebug() << sum << maxV;
    if (sum == 0) maxV = -1e-20;
    *maxv = maxV;
    *minv = minV;
    return maxV;
}
//---------------------------------------------------------------------------
// fill the current raster data structure with new data, called each run step
double lisemqt::fillDrawMapDataRGB(cTMap * base, cTRGBMap *_M, QwtMatrixRasterData *_RD)//, double type)
{
    double maxV = -1e20;
    RGBData.clear();  //QVector double

    if (_M == nullptr)
        return (maxV);

    // copy map data into vector for the display structure
    for(int r = _M->nrRows()-1; r >= 0; r--)
        for(int c=0; c < _M->nrCols(); c++)
        {

            if(true)// !pcr::isMV(_M->dataR[r][c]))
            {
                double value = 0;
                char * valuechar = ((char*)(&value));
                valuechar[0] = _M->dataR[r][c];//*rc;//(*base->data[r][c]);
                if(_M->bands > 1)
                {
                    valuechar[1] = _M->dataG[r][c];//*rg;
                    valuechar[2] = _M->dataB[r][c];//*rb;
                 //   valuechar[3] = base->data[r][c]*255;
                }else
                {
                    valuechar[1] = _M->dataR[r][c];
                    valuechar[2] = _M->dataR[r][c];
                  //  valuechar[3] = base->data[r][c]*255;
                }

                RGBData << value;
                maxV = std::max(maxV, 1.0);
            }
            else
            {
                RGBData << (double)-1e20;
            }
        }

    // set intervals for rasterdata, x,y,z min and max
    _RD->setValueMatrix( RGBData, _M->nrCols() );
    // set column number to divide vector into rows

    _RD->setInterval( Qt::XAxis, QwtInterval( op._llx,op._llx + (double)_M->nrCols()*_M->cellSize(), QwtInterval::ExcludeMaximum ) );
    _RD->setInterval( Qt::YAxis, QwtInterval( op._lly,op._lly + (double)_M->nrRows()*_M->cellSize(), QwtInterval::ExcludeMaximum ) );
    // set x/y axis intervals
    return maxV;
}
//---------------------------------------------------------------------------
void lisemqt::showMapb(bool)
{
    showMap();
}
//---------------------------------------------------------------------------
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
    //initialize and plot at the start
    if(op.comboboxset == false)
    {
        op.comboboxset = true;
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

        for(int i = 0; i < op.ComboMaps.length(); i++)
        {
            QwtComboColorMap *cm1 = new QwtComboColorMap(QColor(op.ComboColors.at(i).at(0)),QColor(op.ComboColors.at(i).at(op.ComboColors.at(i).length()-1)),op.ComboColorMap.at(i),op.ComboColors.at(i));
            QwtComboColorMap *cm2 = new QwtComboColorMap(QColor(op.ComboColors.at(i).at(0)),QColor(op.ComboColors.at(i).at(op.ComboColors.at(i).length()-1)),op.ComboColorMap.at(i),op.ComboColors.at(i));
            cmMap.append(cm1);
            cmLeg.append(cm2);

            NameList.append(op.ComboMapNames.at(i));
            UnitList.append(op.ComboUnits.at(i));
            SymList.append(op.ComboSymColor.at(i)); // symetric colors
            LogList.append(op.ComboLogaritmic.at(i));  //log display
            picker->NameList.append(op.ComboMapNames.at(i));
            picker->UnitList.append(op.ComboUnits.at(i));
        }
        QStringList S;
        QStringList S1;

        for(int i = 0; i < op.ComboMaps.length(); i++)
        {

            if(op.ComboLists.at(i) == 0) {
                S << QString(op.ComboMapNames.at(i) + " (" + op.ComboUnits.at(i) + ")");
                IndexList.append(i);
            } else {
                S1 << QString(op.ComboMapNames.at(i) + " (" + op.ComboUnits.at(i) + ")");
                IndexList1.append(i);
            }

        }
        DisplayComboBox->addItems(S);
        DisplayComboBox2->addItems(S1);
        ActiveList = 0;


        DisplayComboBox->setSizeAdjustPolicy(QComboBox::AdjustToContents);
        DisplayComboBox2->setSizeAdjustPolicy(QComboBox::AdjustToContents);

        DisplayComboBox2->setCurrentIndex(0);
        DisplayComboBox->setCurrentIndex(0);

        DisplayComboBox->setMaxVisibleItems(IndexList.count());
        DisplayComboBox2->setMaxVisibleItems(IndexList1.count());

        MPlot->replot();

    }

    // needed if user clicks while nothing is running:
    if (IndexList.count() == 0)
        return;

    if (tabWidget_out->currentIndex() < 1)
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

    roadMap->setAlpha(checkMapRoads->isChecked() ? transparencyHardSurface->value() : 0);
    houseMap->setAlpha(checkMapBuildings->isChecked() ? transparencyHardSurface->value() : 0);
    hardsurfMap->setAlpha(checkMapHardSurface->isChecked() ? transparencyHardSurface->value() : 0);

    // imageMap->setAlpha(0);  // flow barriers for now not used, sat image instead
    if (checksatImage->isChecked()){
        baseMapImage->setAlpha(checkMapImage->isChecked() ? 255 /*transparencyImage->value()*/ : 0);
    }

    baseMapDEM->setAlpha(checkMapImage->isChecked() ? 0 : 255);

    if (nrcontourlevels->value() > 0)
    {
        contourLevels.clear();
        for ( double level = contourmin; level < contourmax; level += (contourmax-contourmin)/nrcontourlevels->value() )
            contourLevels += level;
        contourDEM->setContourLevels( contourLevels );
        contourDEM->setDisplayMode( QwtPlotSpectrogram::ContourMode, nrcontourlevels->value() > 0 );
    }

    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::showComboMap(int i)
{
    if( i < 0 || i >= op.ComboMaps.length())
        return;
   // qDebug() << i << op.ComboMapNames.at(i);
    MPlot->setTitle(op.ComboMapNames.at(i) + " (" + op.ComboUnits.at(i) + ")");
    // fill vector RD with matrix data and find the new max value
    double MinV;
    double MaxV;
    double res = fillDrawMapData(op.ComboMaps.at(i), op.ComboScaling.at(i), RD, &MinV, &MaxV);
    if (res <=-1e20)
        return;
    //double MinV = mapMinimum(*op.ComboMaps.at(i));
    if ( i < 20 && MaxV < 1e-10)
       return;
    if ( i > 20 && MaxV < 1e-10 && MinV > -1e-10)
       return;


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
    bool domin = true;

    if (ma == 0)
        ma = MaxV; // use map max when ma = 0
    if (mi == 0) {
        mi = MinV;
      //  domin = false;
    }
    if(mi == ma) // because of they are equal nothing is displayed (raincum)
        mi = 0;
    if(mi > ma)
        ma = mi;
    //   qDebug() << mi << ma << MinV << MaxV;

    if (op.ComboSymColor.at(i)) // symetric coloring for soilloss
    {
        if (ma == 0) mi = 0;
        mi = -ma;
        if (ComboMaxSpinBox2->value() > 0)
            ComboMinSpinBox2->setValue(mi);
        else
            ComboMinSpinBox2->setValue(0);

    }
   // qDebug() << mi << ma << MinV << MaxV;


    RD->setInterval( Qt::ZAxis, QwtInterval( mi, ma));

    QwtComboColorMap *cm = new QwtComboColorMap(QColor(op.ComboColors.at(i).at(0)),
                                                QColor(op.ComboColors.at(i).at(op.ComboColors.at(i).length()-1)),
                                                op.ComboColorMap.at(i),op.ComboColors.at(i));
    QwtComboColorMap *cmL = new QwtComboColorMap(QColor(op.ComboColors.at(i).at(0)),
                                                 QColor(op.ComboColors.at(i).at(op.ComboColors.at(i).length()-1)),
                                                 op.ComboColorMap.at(i),op.ComboColors.at(i));
//    cm->setMode(cm->FixedColors);
//    cmL->setMode(cm->FixedColors);
    cm->thresholduse = domin;
    cmL->thresholduse = true;
    cm->thresholdmin = mi;
    cmL->thresholdmin = mi;
    if (op.ComboSymColor.at(i)) {
        cm->thresholdmin = MinV;
        cmL->thresholdmin = mi;
    }
//    cmMap.at(i)->thresholduse = domin;
//    cmLeg.at(i)->thresholduse = true;
//    cmMap.at(i)->thresholdmin = mi;
//    cmLeg.at(i)->thresholdmin = mi;
//    if (op.ComboSymColor.at(i)) {
//        cmMap.at(i)->thresholdmin = MinV;
//        cmLeg.at(i)->thresholdmin = mi;
//    }

    drawMap->setData(RD);
    drawMap->setColorMap(cm);//Map.at(i));
    drawMap->setAlpha(transparencyMap->value());

    rightAxis->setColorMap( drawMap->data()->interval( Qt::ZAxis ), cmL);//eg.at(i));

    if(op.ComboLogaritmic.at(i))
    {
        int coef = int(log10(ma));
        mi = (mi == 0 ? std::pow(0.1,3-coef) : mi);
        ma = (ma == 0 ? std::pow(10,coef)    : ma);

        MPlot->setAxisScale( MPlot->yRight, mi, ma );
        MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLogScaleEngine() );
    }
    else
    {

        MPlot->setAxisScale( MPlot->yRight, mi, ma);
        MPlot->setAxisScaleEngine( MPlot->yRight, new QwtLinearScaleEngine() );
    }

}
//---------------------------------------------------------------------------
void lisemqt::showBaseMap()
{
    if (!startplot)
        return;
    // only done once

    double m1, m2;
    double res = fillDrawMapData(op.baseMap, 1.0, RDb, &m1, &m2);
    if (res == -1e20)
        return;

    baseMap->setAlpha(transparency->value());
    baseMap->setColorMap(new colorMapGray());
    RDb->setInterval( Qt::ZAxis, QwtInterval( 0, m2));
    baseMap->setData(RDb);
    // setdata sets a pointer to DRb to the private QWT d_data Qvector

    res = fillDrawMapData(op.baseMapDEM, 1.0, RDbb, &m1, &m2);
    if (res == -1e20)
        return;
    //double mindem = mapMinimum(*op.baseMapDEM);
    contourmax = m2;//mapMaximum(*op.baseMapDEM);
    contourmin = m1;//mindem;

    baseMapDEM->setAlpha(checkMapImage->isChecked() ? 0 : 255);
    baseMapDEM->setColorMap(new colorMapElevation());
    RDbb->setInterval( Qt::ZAxis, QwtInterval( m1, m2));
    baseMapDEM->setData(RDbb);

    contourDEM->setAlpha(255);
    contourDEM->setColorMap(new colorMapTransparent());
    RDbb->setInterval( Qt::ZAxis, QwtInterval( m1, m2));
    contourDEM->setData(RDbb);

    tabWidget_out->setCurrentIndex(1);

    // fit into screen
    changeSize();

}
//---------------------------------------------------------------------------
void lisemqt::showChannelVector(bool yes)
{
    if (!checkIncludeChannel->isChecked())
        return;

    if(rivers.isEmpty())
        return;

    if (yes) {
        // attach everything and use new pensizes

        QPen pen1;
        pen1.setWidth(spinChannelSize->value());
        pen1.setColor(QColor("#000000"));
        pen1.setCosmetic(true);

        for (int i = 0; i < rivers.length(); i++) {
            rivers[i]->setPen(pen1);
            rivers[i]->attach( MPlot );
            rivers[i]->setAxes(MPlot->xBottom, MPlot->yLeft);
        }

        QPen pen2;
        pen2.setWidth(spinChannelSize->value()+1);
        pen2.setColor(QColor("#FFFFFF"));
        pen2.setCosmetic(true);

        // culverts get white channel color
        for (int i = 0; i < culverts.length(); i++) {
            culverts[i]->setPen(pen2);
            culverts[i]->attach( MPlot );
            culverts[i]->setAxes(MPlot->xBottom, MPlot->yLeft);
        }

        int dxi = spinCulvertSize->value();
        outlets.setSymbol(new QwtSymbol( QwtSymbol::Ellipse, Qt::white, QPen( Qt::black ), QSize( dxi,dxi )));
        outlets.attach( MPlot );

        obspoints.setSymbol(new QwtSymbol( QwtSymbol::Ellipse, Qt::cyan, QPen( Qt::black ), QSize( dxi,dxi )));
        obspoints.attach( MPlot );

    } else {
        // detach everything

        obspoints.detach();

        outlets.detach();

        if (culverts.length() > 0 && !culverts.isEmpty()) {
            for (int i = 0; i < culverts.length(); i++)
                culverts[i]->detach();
        }

        if (rivers.length() > 0 && !rivers.isEmpty()) {
            for (int i = 0; i < rivers.length(); i++)
                rivers[i]->detach();
        }

    }

    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::showChannelVectorNew()
{
    if (!checkIncludeChannel->isChecked())
        return;

    // fill the line and dot structures once at startplot
    // draw once with showChannelVector(true);
    if (startplot) {

        // Channel network
        QVector <double> X;
        QVector <double> Y;
        checkMapChannels->setChecked(true);
        int _dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
        int _dy[10] = {0, -1,-1,-1,  0, 0, 0,  1, 1, 1};

        double xend, yend;
        double dx = op.channelMap->cellSize();
        double _nrRows = op.channelMap->nrRows();
        double cy = op.channelMap->north()-_nrRows*dx;
        double cx = op.channelMap->west();

        for(long i_ =  0; i_ < op.lddch_.size(); i_++) {

            int r = _nrRows-op.lddch_[i_].r-1;
            int c = op.lddch_[i_].c;
            int ldd = std::abs(op.lddch_[i_].ldd);

            xend = cx+c*dx + 0.5*dx;
            yend = cy+r*dx + 0.5*dx;
            if (op.lddch_[i_].nr == 0) {
                Y << yend;
                X << xend;
            }
            Y << yend + _dy[ldd]*dx;
            X << xend + _dx[ldd]*dx;

            if (i_ < op.lddch_.size()-1 && op.lddch_[i_+1].nr == 0) {
                Xa.push_back(X);
                Ya.push_back(Y);
                X.clear();
                Y.clear();
            }
        }
        Xa.push_back(X);
        Ya.push_back(Y);
        X.clear();
        Y.clear();

        for (int i = 0; i < Xa.length(); i++) {
            QwtPlotCurve *river = new QwtPlotCurve();
            river->setSamples(Xa.at(i),Ya.at(i));
            rivers << river;
        }

        // culvert parts of channel network
        if (checkChannelCulverts->isChecked()) {
            bool first = true;
            int count = 0;
            for(long i_ =  0; i_ < op.lddch_.size(); i_++) {

                if(op.lddch_[i_].ldd < 0) {
                    int r,c,ldd;

                    if (first) {
                        r = _nrRows-op.lddch_[i_].r-1;
                        c = op.lddch_[i_].c;
                        xend = cx+c*dx + 0.5*dx;
                        yend = cy+r*dx + 0.5*dx;
                        Y << yend;
                        X << xend;
                        first = false;
                    }
                    r = _nrRows-op.lddch_[i_].r-1;
                    c = op.lddch_[i_].c;
                    xend = cx+c*dx + 0.5*dx;
                    yend = cy+r*dx + 0.5*dx;
                    Y << yend;
                    X << xend;
                    count++;

                    if (i_ < op.lddch_.size() -1 && (op.lddch_[i_+1].ldd > 0)) {
                        if (count == 1) {
                            ldd = std::abs(op.lddch_[i_].ldd);
                            Y << yend + _dy[ldd]*dx;
                            X << xend + _dx[ldd]*dx;
                        }
                        Xc.push_back(X);
                        Yc.push_back(Y);
                        X.clear();
                        Y.clear();
                        first = true;
                        count = 0;
                    }
                }
            }

            Xc.push_back(X);
            Yc.push_back(Y);
            X.clear();
            Y.clear();

            for (int i = 0; i < Xc.length(); i++) {
                QwtPlotCurve *culvert = new QwtPlotCurve();
                culvert->setSamples(Xc.at(i),Yc.at(i));
                culverts << culvert;
            }
        }

        // dot size
        int dxi = MPlot->invTransform(MPlot->xBottom,dx*1.5);
        dxi = dxi - MPlot->invTransform(MPlot->xBottom,dx);        
        dxi = std::min(9,dxi);
        spinCulvertSize->setValue(dxi);

        // points in outlet.map
        outlets.setSymbol(new QwtSymbol( QwtSymbol::Ellipse, Qt::white, QPen( Qt::black ), QSize( dxi,dxi )));
        outlets.setPen( Qt::black );
        outlets.setStyle( QwtPlotCurve::NoCurve );
        outlets.setAxes(MPlot->xBottom, MPlot->yLeft);
        outlets.setSamples(op.EndPointX,op.EndPointY);

        // points in outpoint.map
        obspoints.setSymbol(new QwtSymbol( QwtSymbol::Ellipse, Qt::cyan, QPen( Qt::black ), QSize( dxi,dxi )));
        obspoints.setPen( Qt::black );
        obspoints.setStyle( QwtPlotCurve::NoCurve );
        obspoints.setAxes(MPlot->xBottom, MPlot->yLeft);
        obspoints.setSamples(op.ObsPointX,op.ObsPointY);

        // clear all structures here for the next run of a different area
        Xa.clear();
        Ya.clear();
        Xc.clear();
        Yc.clear();
        op.ObsPointX.clear();
        op.ObsPointY.clear();
        op.EndPointX.clear();
        op.EndPointY.clear();

        // attach everything to MPlot and draw
        showChannelVector(true);

    } // startplot
}

void lisemqt::getOutletMap()
{
    if (!checkIncludeChannel->isChecked())
        return;

    if (startplot)
    {
        double m1, m2;
        double res = fillDrawMapData(op.outletMap, 1.0, RDc, &m1, &m2);
        if (res ==-1e20)
            return;
        RDc->setInterval( Qt::ZAxis, QwtInterval( 0,1.0));
        outletMap->setData(RDc);
        outletMap->setAlpha(0);
    }

}
//---------------------------------------------------------------------------
void lisemqt::showRoadMap()
{
    if (startplot)
    {
        double m1, m2;
        double res = fillDrawMapData(op.roadMap,1.0, RDd, &m1, &m2);
        if (res ==-1e20)
            return;
        RDd->setInterval( Qt::ZAxis, QwtInterval( 0,0.5));
        roadMap->setData(RDd);
    }

    if (checkMapRoads->isChecked())
        roadMap->setAlpha(transparencyHardSurface->value());
    else
        roadMap->setAlpha(0);

    roadMap->setColorMap(new colorMapRoads2());
}
//---------------------------------------------------------------------------
void lisemqt::showHouseMap()
{
    if (startplot)
    {
        double m1, m2;
        double res = fillDrawMapData(op.houseMap,1.0, RDe, &m1, &m2);

        if (res ==-1e20)
            return;

        RDe->setInterval( Qt::ZAxis, QwtInterval( m1, m2));
        houseMap->setData(RDe);
    }

    if (checkMapBuildings->isChecked())
        houseMap->setAlpha(transparencyHardSurface->value());
    else
        houseMap->setAlpha(0);
    houseMap->setColorMap(new colorMapHouse());
}
//---------------------------------------------------------------------------
void lisemqt::showHardSurfaceMap()
{

    if (startplot)
    {
        double m1, m2;
        // set intervals for rasterdata, x,y,z min and max
        double res = fillDrawMapData(op.hardsurfaceMap,1.0, RDf, &m1, &m2);
        if (res ==-1e20)
            return;

        RDf->setInterval( Qt::ZAxis, QwtInterval( m1, m2));
        hardsurfMap->setData(RDf);
    }

    if (checkHardsurface->isChecked())
        hardsurfMap->setAlpha(transparencyHardSurface->value());
    else
        hardsurfMap->setAlpha(0);
    hardsurfMap->setColorMap(new colorMapRoads());
}
//---------------------------------------------------------------------------
void lisemqt::showImageMap()
{
    if (startplot && checksatImage->isChecked())
    {
        // set intervals for rasterdata, x,y,z min and max
//        double res = fillDrawMapDataRGB(op.baseMapDEM,op.Image, RImage);
        double res = fillDrawMapDataRGB(op.baseMap,op.Image, RImage);
        RImage->setInterval( Qt::ZAxis, QwtInterval( 0.0, 1.0));
        baseMapImage->setData(RImage);
    }

    if (checksatImage->isChecked())
        baseMapImage->setAlpha(255);//transparencyImage->value());
    else
        baseMapImage->setAlpha(0);

    baseMapImage->setColorMap(new colorMapRGB());
}
//---------------------------------------------------------------------------
