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
    if (spinChannelSize->value() > 0)
        hideChannelVector(true);
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaChannel(int v)
{
    hideChannelVector(false);
    hideChannelVector(v > 0);
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
    //    transvalue = ((double) v)/256.0;
    //    doHouse = true;
    //    showHouseMap();
    if (v > 0 && checkMapBuildings->isChecked())
        MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::ssetAlphaBarrier(int v)
{
    imageMap->setAlpha(v);
    if (v > 0 && checkMapImage->isChecked())
        MPlot->replot();
    // flowbarrier is now sat image, highjacked!
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
}

void lisemqt::changeSize()
{
    int i = 0;

    MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrCols*op._dx, op._dx*100);
    MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrRows*op._dx, op._dx*100);

    if(op._nrCols/op._nrRows > Masp) {
        MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrCols*op._dx, op._dx*100);
        MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrCols*op._dx, op._dx*100);
        i = 1;
    } else
        if(op._nrRows > op._nrCols)
        {
            MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrRows*op._dx*Masp, op._dx*100);
            MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrRows*op._dx*Masp, op._dx*100);
            i = 2;
        }

    MPlot->replot();
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
    baseMapDEM->setRenderThreadCount( 2 );
    baseMapDEM->attach( MPlot );
    // dem
    // 1
    baseMapImage = new QwtPlotSpectrogram();
    baseMapImage->setRenderThreadCount( 2 );
    baseMapImage->attach( MPlot );
    //image

    // 2
    baseMap = new QwtPlotSpectrogram();
    baseMap->setRenderThreadCount( 2 );
    baseMap->attach( MPlot );
    // shaded relief

    // 3 data
    drawMap = new QwtPlotSpectrogram();
    drawMap->setRenderThreadCount( 2 );
    drawMap->attach( MPlot );
    //map for runoff, infil, flood etc

    // 4
    houseMap = new QwtPlotSpectrogram();
    houseMap->setRenderThreadCount( 2 );
    houseMap->attach( MPlot );
    // building structure map

    // 5
    roadMap = new QwtPlotSpectrogram();
    roadMap->setRenderThreadCount( 2 );
    roadMap->attach( MPlot );
    // road map

    // 6
    //imageMap = new QwtPlotSpectrogram();
    //imageMap->setRenderThreadCount( 0 );
    //imageMap->attach( MPlot );


    //7
    outletMap = new QwtPlotSpectrogram();
    outletMap->setRenderThreadCount( 2 );
    outletMap->attach( MPlot );
    // channel map-

    //8
    contourDEM = new QwtPlotSpectrogram();
    contourDEM->setRenderThreadCount( 2 );
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


    panner = new QwtPlotPanner( MPlot->canvas() );
    panner->setAxisEnabled( MPlot->yRight, false );
    // exclude right axis legend from panning

    picker = new MyPicker( (QwtPlotCanvas *) MPlot->canvas() );
    picker->setEnabled(false);

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

//                char r = (_M->dataR[r][c]);
//                char g = (_M->dataG[r][c]);
//                char b = (_M->dataB[r][c]);


//                char rc[1];
//                char rg[1];
//                char rb[1];
//                itoa(r,rc,10);
//                itoa(g,rg,10);
//                itoa(b,rb,10);
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

    //channelMap->setAlpha(checkMapChannels->isChecked() ? transparencyChannel->value() : 0);
    roadMap->setAlpha(checkMapRoads->isChecked() ? transparencyRoad->value() : 0);
    houseMap->setAlpha(checkMapBuildings->isChecked() ? transparencyHouse->value() : 0);

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
    double res = fillDrawMapData(op.ComboMapsSafe.at(i), op.ComboScaling.at(i), RD, &MinV, &MaxV);
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

        mi = -ma;
        if (ComboMaxSpinBox2->value() > 0)
            ComboMinSpinBox2->setValue(mi);

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

    // reset the axes to the correct rows/cols,
    // do only once because resets zooming and panning

    // fit into screen first time
    tabWidget_out->setCurrentIndex(1);
 //   MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrCols*op._dx, op._dx*10);
 //   MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrRows*op._dx, op._dx*10);


    //    int h = MPlot->height();
    //    int w = MPlot->width();
    //    double asp = (double)w/(double)h;
    //    int i = 0;

    //tabWidget_out->setCurrentIndex(1);

    changeSize();

//    MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrCols*op._dx, op._dx*10);
//    MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrRows*op._dx, op._dx*10);

//    if(op._nrCols/op._nrRows > asp) {
//        MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrCols*op._dx, op._dx*10);
//        MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrCols*op._dx, op._dx*10);
//        i = 1;
//    } else
//        if(op._nrRows > op._nrCols)
//        {
//            MPlot->setAxisScale( MPlot->xBottom, op._llx, op._llx+(double)op._nrRows*op._dx*asp, op._dx*10);
//            MPlot->setAxisScale( MPlot->yLeft, op._lly, op._lly+(double)op._nrRows*op._dx*asp, op._dx*10);
//            i = 2;
//        }

 //   tabWidget_out->setCurrentIndex(0);

}
//---------------------------------------------------------------------------
void lisemqt::hideChannelVector(bool yes)
{
    if (!checkIncludeChannel->isChecked())
        return;

    if(rivers.isEmpty())
        return;

    if (!yes) {

        for (int i = 0; i < rivers.length(); i++)
            rivers[i]->detach();

        if (checkChannelCulverts->isChecked() && !culverts.isEmpty()) {
        if (culverts.length() > 0)
            for (int i = 0; i < culverts.length(); i++)
                culverts[i]->detach();
        }
        if (!obspoints.isEmpty()) {
        if (obspoints.length() > 0)
            for (int i = 0; i < obspoints.length(); i++)
                obspoints[i]->detach();
        }

        for (int i = 0; i < outlets.length(); i++)
            outlets[i]->detach();
    } else {
        QPen pen1;
        pen1.setWidth(spinChannelSize->value());//showRiverSize->value());
        pen1.setColor(QColor("#000000"));
        pen1.setCosmetic(true);

        for (int i = 0; i < rivers.length(); i++) {
            rivers[i]->setPen(pen1);
            rivers[i]->attach( MPlot );
            rivers[i]->setAxes(MPlot->xBottom, MPlot->yLeft);
        }

        for (int i = 0; i < obspoints.length(); i++) {
            int dxi = spinCulvertSize->value();

            QwtSymbol *bluedot = new QwtSymbol( QwtSymbol::Ellipse, Qt::cyan, QPen( Qt::black ), QSize( dxi,dxi ));

            obspoints[i]->setSymbol(bluedot);
            obspoints[i]->setStyle( QwtPlotCurve::NoCurve );
            obspoints[i]->attach( MPlot );
            obspoints[i]->setAxes(MPlot->xBottom, MPlot->yLeft);
        }

        for (int i = 0; i < outlets.length(); i++) {
            int dxi = spinCulvertSize->value();

            QwtSymbol *whitedot = new QwtSymbol( QwtSymbol::Ellipse, Qt::white, QPen( Qt::black ), QSize( dxi,dxi ));

            outlets[i]->setSymbol(whitedot);
            outlets[i]->setStyle( QwtPlotCurve::NoCurve );
            outlets[i]->attach( MPlot );
            outlets[i]->setAxes(MPlot->xBottom, MPlot->yLeft);
        }
    }

    MPlot->replot();
}
//---------------------------------------------------------------------------
void lisemqt::showChannelVectorNew()
{
    if (!checkIncludeChannel->isChecked())
        return;

    if (startplot) {
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

        for(long i_ =  0; i_ < op.lddch_.size(); i_++)
        {
            int r = _nrRows-op.lddch_[i_].r-1;
            int c = op.lddch_[i_].c;
            double ldd = op.lddch_[i_].ldd;

            xend = cx+c*dx + 0.5*dx;
            yend = cy+r*dx + 0.5*dx;
            if (op.lddch_[i_].nr == 0) {
                Y << yend;
                X << xend;
            }
            Y << yend + _dy[(int)ldd]*dx;
            X << xend + _dx[(int)ldd]*dx;

            if (i_ < op.lddch_.size() -1 && op.lddch_[i_+1].nr == 0) {
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

        QPen pen1;
        pen1.setWidth(spinChannelSize->value());
        pen1.setColor(QColor("#000000"));
        pen1.setCosmetic(false);
        for (int i = 0; i < Xa.length(); i++) {
            rivera = new QwtPlotCurve();
            rivers << rivera;

            rivera->setPen(pen1);
            rivera->attach( MPlot );
            rivera->setAxes(MPlot->xBottom, MPlot->yLeft);
            rivera->setSamples(Xa.at(i),Ya.at(i));
        }

        int dxi = (int) (op.channelMap->cellSize()*0.5);
        dxi = std::min(5,dxi);
        spinCulvertSize->setValue(dxi);

        QwtPlotCurve *outlet = new QwtPlotCurve();
        QwtSymbol *whitedot = new QwtSymbol( QwtSymbol::Ellipse, Qt::white, QPen( Qt::black ), QSize( dxi,dxi ));
        outlets << outlet;
        outlet->setSymbol(whitedot);
        outlet->setStyle( QwtPlotCurve::NoCurve );
        outlet->attach( MPlot );
        outlet->setAxes(MPlot->xBottom, MPlot->yLeft);
        outlet->setSamples(op.EndPointX,op.EndPointY);

        // clear here for the next run of a different area
        op.EndPointX.clear();
        op.EndPointY.clear();
        Xa.clear();
        Ya.clear();

        QwtSymbol *greendot = new QwtSymbol( QwtSymbol::Ellipse, Qt::green,QPen( Qt::black ), QSize( dxi, dxi ) );

        if (checkChannelCulverts->isChecked()) {
            culvert = new QwtPlotCurve();
            culverts << culvert;
            culvert->setSymbol(greendot);
            culvert->setStyle( QwtPlotCurve::NoCurve );
            culvert->attach( MPlot );
            culvert->setAxes(MPlot->xBottom, MPlot->yLeft);
            culvert->setSamples(op.CulvertX,op.CulvertY);
        }

        QwtSymbol *bluedot = new QwtSymbol( QwtSymbol::Ellipse, Qt::cyan,QPen( Qt::black ), QSize( dxi, dxi ) );

        obspoint = new QwtPlotCurve();
        obspoints << obspoint;
        obspoint->setSymbol(bluedot);
        obspoint->setStyle( QwtPlotCurve::NoCurve );
        obspoint->attach( MPlot );
        obspoint->setAxes(MPlot->xBottom, MPlot->yLeft);
        obspoint->setSamples(op.ObsPointX,op.ObsPointY);

        op.CulvertX.clear();
        op.CulvertY.clear();
    }
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
        roadMap->setAlpha(transparencyRoad->value());
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
        doHouse = false;
    }
    if (checkMapBuildings->isChecked())
        houseMap->setAlpha(transparencyHouse->value());
    else
        houseMap->setAlpha(0);
    houseMap->setColorMap(new colorMapHouse());
}
//---------------------------------------------------------------------------
// NOT USED FOR NOW
void lisemqt::showFlowBarriersMap()
{

//    if (startplot)
//    {
//        double m1, m2;
//        // set intervals for rasterdata, x,y,z min and max
//        double res = fillDrawMapData(op.flowbarriersMap,1.0, RDf, &m1, &m2);
//        if (res ==-1e20)
//            return;
//        RDf->setInterval( Qt::ZAxis, QwtInterval( m1, m2));
//        imageMap->setData(RDf);
//    }

//    if (checkMapImage->isChecked())
//        imageMap->setAlpha(transparencyImage->value());
//    else
//        imageMap->setAlpha(0);

//    imageMap->setColorMap(new colorMapImage());

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
