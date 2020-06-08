#ifndef LISUIMAPCOLOR_H
#define LISUIMAPCOLOR_H

#include "lisemqt.h"

//---------------------------------------------------------------------------
#define BGc "#eeeeee" // background grey for missing value in maps

//http://www.color-hex.com/color/
//http://www.colorschemer.com/schemes/index.php?tag=nature

//---------------------------------------------------------------------------
/// Shows value on cursor in map window
class MyPicker: public QwtPlotPicker
{
public:

    QList<QString> NameList;
    QList<QString> UnitList;

    MyPicker( QwtPlotCanvas *canvas ):
        QwtPlotPicker( canvas )
    {
        setTrackerMode( AlwaysOn );
    }

    virtual QwtText trackerTextF( const QPointF &pos ) const
    {
        QColor bg( Qt::white );
        bg.setAlpha( 100 );
        QString txt = "";

        //layer 0 is dem, layer 1 is shade, layer 3 is thematic
        QwtPlotItemList list = plot()->itemList(QwtPlotItem::Rtti_PlotSpectrogram);
        QwtPlotSpectrogram * sp2 = static_cast<QwtPlotSpectrogram *> (list.at(3));
        QwtPlotSpectrogram * sp0 = static_cast<QwtPlotSpectrogram *> (list.at(0));
//        // elevation info
        QwtPlotSpectrogram * sp3 = static_cast<QwtPlotSpectrogram *> (list.at(7));
//        // outlet info

        if (sp2->data() == nullptr)
            return QwtText(txt);
        double z2 = sp2->data()->value(pos.x(), pos.y());
        double z0 = sp0->data()->value(pos.x(), pos.y());
        int z3 = 0;
        if (sp3->data() != nullptr)
            z3 = (int) sp3->data()->value(pos.x(), pos.y());

        if (z2 > -1e10)
        {
            int dig = 2;
            if (z0 < 1.0)
                dig = 4;
            if (fabs(z2) < 1.0)
                txt = (QString("%1 ")/* + unit*/ + QString(" [%2m]")).arg(z2,0,'g',3).arg(z0,0,'f',dig);
            else
                txt = (QString("%1 ")/* + unit */+ QString(" [%2m]")).arg(z2,0,'f',3).arg(z0,0,'f',dig);
        }
        if (z3 > 0) txt = txt + QString(" (Outlet %1)").arg(z3);

        QwtText text = QwtText(txt);
        text.setColor(Qt::black);
        text.setBackgroundBrush( QBrush( bg ) );
        return text;
    }
};

//---------------------------------------------------------------------------
// class derived from QwtLinearColorMap to enable transparency thresholds
class QwtLinearColorMapVJ: public QwtLinearColorMap
{
public:

    QwtLinearColorMapVJ( const QColor &from, const QColor &to,
        QwtLinearColorMap::Format = QwtColorMap::RGB );

    virtual ~QwtLinearColorMapVJ();

    double thresholdLCM;

    void setThreshold(double v);

};
//---------------------------------------------------------------------------
// class derived from QwtLinearColorMap to enable transparency thresholds
class QwtComboColorMap: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {

        if(thresholduse)
        {
            if ( value <= thresholdmin )
            {
                return qRgba( 0, 0, 0, 0 );
            }
        }else
        {
            /*if ( value <= interval.minValue() )
            {
                return qRgba( 0, 0, 0, 0 );
            }*/
        }
        return QwtLinearColorMap::rgb( interval, value );
    }

public:

    QwtComboColorMap( const QColor &from, const QColor &to,
                      QList<double> map,
                      QList<QString> colors,
        QwtLinearColorMap::Format = QwtColorMap::RGB ):QwtLinearColorMap( from, to)
    {
       //     thresholdmin = 0;
       //     thresholduse = false;

            setColorInterval(from,to);

            for(int j = 1; j < colors.length()-1; j++)
            {
                addColorStop(map.at(j),QColor(colors.at(j)));
            }
    }

    virtual ~QwtComboColorMap()
    {

    }

    double thresholdmin;
    bool thresholduse;
};
//---------------------------------------------------------------------------
class colorMapRGB: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {

        char * valuechar = ((char*)(&value));
        return qRgba( (int) valuechar[0], (int) valuechar[1], (int) valuechar[2], 255 );
    }
public:
    colorMapRGB():
        QwtLinearColorMap( QColor("#555555"),QColor("#ffffff"))
    {

    }
};
//---------------------------------------------------------------------------
/// for contour layer map, transparent
class colorMapTransparent: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < 1e20 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapTransparent():
          QwtLinearColorMap( QColor("#AAAAAA"), QColor("#222222"))
    { }
};
//---------------------------------------------------------------------------
/// House color legend
class colorMapHouse: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < 0.05 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapHouse():
//        QwtLinearColorMap( QColor("#333300"), QColor("#ada399"))
 //     QwtLinearColorMap( QColor("#c06969"), QColor("#421010"))
          QwtLinearColorMap( QColor("#AAAAAA"), QColor("#222222"))
    {

    }
};
//---------------------------------------------------------------------------
/// flow barrier map
class colorMapFlowBarrier: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < 0.05 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFlowBarrier():
        QwtLinearColorMap( QColor("#000000"), QColor("#FFFF00"))
    {

    }
};


//---------------------------------------------------------------------------
///  relief map display
///http://www.colorschemer.com/schemes/viewscheme.php?id=81
class colorMapElevation: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );
        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );
        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapElevation():
//        QwtLinearColorMap( QColor("#B17142"),QColor("#FFDFC7"))
//              QwtLinearColorMap( QColor("#8d5524"),QColor("#ffdbac").lighter())
 //                           QwtLinearColorMap( QColor("#D7191C"),QColor("#ffffbf"))
  QwtLinearColorMapVJ( QColor(141,116,94),QColor(255,251,244)) //skincolor
  //     QwtLinearColorMap( QColor(188,170,39),QColor(20,129,2))
   //     QwtLinearColorMap( QColor("#279400").lighter(),QColor("#F3EEC4").lighter())//kurt arcgis
    {
//        addColorStop(0.000,QColor("#B17142"));
//        addColorStop(0.250,QColor("#C48C63"));
//        addColorStop(0.500,QColor("#D8A885"));
//        addColorStop(0.750,QColor("#ECC4A6"));
//        addColorStop(1.000,QColor("#FFDFC7"));

//        addColorStop(0.250,QColor("#c68642"));
//        addColorStop(0.500,QColor("#e0ac69"));
//        addColorStop(0.750,QColor("#ECC4A6"));

//      addColorStop(0.500,QColor("#fdae61"));


        addColorStop(0.250,QColor(198,164,136)); // "skincolor"
        addColorStop(0.500,QColor(224,201,173));
        addColorStop(0.750,QColor(236,226,214));

//addColorStop(0.250,QColor(119,14,11));
//addColorStop(0.500,QColor(245,251,91));
//addColorStop(0.750,QColor(20,129,2));

//        addColorStop(0.5,QColor("#EAD49A"));

    }
};
//---------------------------------------------------------------------------
/// Gray scale legend for shaded relief map display
class colorMapGray: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );
        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapGray():
        QwtLinearColorMapVJ( QColor("#555555"),QColor("#ffffff"))
    {   }
};
//---------------------------------------------------------------------------
// /// Gray scale legend for shaded relief map display
//class colorMapGrayGap: public QwtLinearColorMapVJ
//{
//    virtual QRgb rgb( const QwtInterval &interval, double value ) const
//    {
//        if ( value <= thresholdLCM )
//            return qRgba( 0, 0, 0, 0 );
//        if ( value > 0.3 && value < 0.7)
//            return qRgba( 0, 0, 0, 0 );
//        return QwtLinearColorMap::rgb( interval, value );
//    }
//public:
//    colorMapGrayGap():
//        QwtLinearColorMapVJ( QColor("#111111"),QColor("#ffffff"))
//    {
////        addColorStop(0.000,QColor("#111111"));
//    }
//};
//---------------------------------------------------------------------------
/// Dark yellow legend for maps for road map overlay
class colorMapRoads3: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapRoads3():
        QwtLinearColorMapVJ( QColor("#ffffff"), QColor("#ffffff"))
    {
        addColorStop(0.0, QColor("#ffffff"));
    }
};
//---------------------------------------------------------------------------
/// Dark yellow legend for maps for road map overlay
class colorMapRoads: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapRoads():
        QwtLinearColorMapVJ( QColor("#277e00"), QColor("#277e00"))
    {
        addColorStop(0.0, QColor("#277e00"));
    }
};
//---------------------------------------------------------------------------
/// Green legend for maps for road map overlay
class colorMapRoads2: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapRoads2():
        QwtLinearColorMapVJ( QColor("#f9fb44"), QColor("#f9fb44")  )
    {
        addColorStop( 0.0, QColor("#f9fb44") );
    }
};
//---------------------------------------------------------------------------
/// Logarithmic light bluie to blue to red legend for runoff map display
class colorMapWaterLog: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );
        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }

public:
  colorMapWaterLog():
    QwtLinearColorMapVJ( QColor("#8c8cff"), QColor("#ff3300"))

  {
      addColorStop( 0.0005,QColor("#8080FF"));
      addColorStop( 0.01, QColor("#4040ff") );
      addColorStop( 0.05, QColor("#0000FF"));
      addColorStop( 0.1, QColor("#00006F"));
      addColorStop( 0.9, QColor("#FF0000"));
  }
};
//---------------------------------------------------------------------------
/// Linear Yellow to blue legend for infil map display
class colorMapWater: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );
        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapWater():
        QwtLinearColorMapVJ( QColor("#f6f666"),QColor("#0000AA"))
    {
        addColorStop( 0.1, QColor("#FFFF55") );
        addColorStop( 0.4, QColor("#8080FF") );
        addColorStop( 0.9, Qt::blue );
    }
};
//---------------------------------------------------------------------------
/// Transparent  light to dark blue legend for flood display
class colorMapFlood: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFlood():
        //#8c8cff
        QwtLinearColorMapVJ( QColor("#5477ff"), QColor("#001462"))
    {
        addColorStop(0.500,QColor("#0023b1"));
    }
};
//---------------------------------------------------------------------------
/// Cyan to red legend for sediment display
class colorMapSed: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

                if ( value < thresholdLCM && value > -thresholdLCM)
                    return qRgba( 0, 0, 0, 0 );
        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapSed():
        QwtLinearColorMapVJ( Qt::darkCyan,Qt::red)
    {
        addColorStop( 0.0, QColor("#457A60"));//darkCyan );
        addColorStop( 0.3, QColor("#96B547"));//cyan );
        addColorStop( 0.5, Qt::white );
        addColorStop( 0.7, Qt::yellow);
    }
};
//---------------------------------------------------------------------------
/// Blue to red legend for sediment display
class colorMapSedB: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value < thresholdLCM && value > -thresholdLCM)
            return qRgba( 0, 0, 0, 0 );
        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapSedB():
        QwtLinearColorMapVJ( Qt::blue,Qt::red)
    {
        addColorStop( 0.0, Qt::blue );
        addColorStop( 0.3, Qt::cyan );
        addColorStop( 0.5, Qt::white );
        addColorStop( 0.7, Qt::yellow);
    }
};
//---------------------------------------------------------------------------
class colorMapFloodV: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value < thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFloodV():
        QwtLinearColorMapVJ( QColor("#007300"),QColor("#A60000"))
    {
        addColorStop( 0.0, Qt::green);
        addColorStop( 0.25, Qt::yellow);
        addColorStop( 0.75, Qt::red);
    }
};
//---------------------------------------------------------------------------
class colorMapP: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value < thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapP():
        QwtLinearColorMapVJ( QColor(BGc),Qt::red)
    {
        addColorStop(0.0,QColor("#8888FF"));
        addColorStop( 0.5, Qt::blue);
    }
};
//---------------------------------------------------------------------------
class colorMapFEW: public QwtLinearColorMapVJ
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
//        if ( value < -1e19 )
//            return qRgba( 228, 228, 228, 255 );

        if ( value < thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFEW():
        QwtLinearColorMapVJ( QColor("#A60000"),QColor("#007300"))
    {
        addColorStop( 0.25, Qt::red);
        addColorStop( 0.5, Qt::yellow);
        addColorStop( 0.75, Qt::green);
    }
};

#endif // LISUIMAPCOLOR_H

