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
        QString unit = "";

        //layer 0 is dem, layer 1 is shade, layer 3 is thematic
        QwtPlotItemList list = plot()->itemList(QwtPlotItem::Rtti_PlotSpectrogram);
        QwtPlotSpectrogram * sp2 = static_cast<QwtPlotSpectrogram *> (list.at(1));
        QwtPlotSpectrogram * sp0 = static_cast<QwtPlotSpectrogram *> (list.at(0));
        // elevation info

        if (sp2->data() == NULL)
            return QwtText(txt);
        double z2 = sp2->data()->value(pos.x(), pos.y());
        double z0 = sp0->data()->value(pos.x(), pos.y());

        if (z2 > -1e10)
        {
            if (sp2->data()->value(0,0) == 1) txt = QString("%1 l/s [%2m]").arg(z2,0,'f',1).arg(z0,0,'f',1);
            if (sp2->data()->value(0,0) == 2) txt = QString("%1 mm [%2m]").arg(z2,0,'f',1).arg(z0,0,'f',1);
            if (sp2->data()->value(0,0) == 3) txt = QString("%1 %2[%3m]").arg(z2,0,'f',1).arg(unit).arg(z0,0,'f',1);
            if (sp2->data()->value(0,0) == 4) txt = QString("%1 m [%2m]").arg(z2,0,'f',2).arg(z0,0,'f',1);
            if (sp2->data()->value(0,0) == 5) txt = QString("%1 m/s [%2m]").arg(z2,0,'f',2).arg(z0,0,'f',1);
            if (sp2->data()->value(0,0) == 6) txt = QString("%1 mm [%2m]").arg(z2,0,'f',3).arg(z0,0,'f',1);
            if (sp2->data()->value(0,0) == 7) txt = QString("%1 min [%2m]").arg(z2,0,'f',3).arg(z0,0,'f',1);
        }

        QwtText text = QwtText(txt);
        text.setColor(Qt::black);
        text.setBackgroundBrush( QBrush( bg ) );
        return text;
    }
};
//---------------------------------------------------------------------------
/// House color legend
class colorMapHouse: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value < 0.05 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapHouse():
//        QwtLinearColorMap( QColor("#333300"), QColor("#ada399"))
 //     QwtLinearColorMap( QColor("#c06969"), QColor("#421010"))
          QwtLinearColorMap( QColor("#777777"), QColor("#222222"))
    {

    }
};
//---------------------------------------------------------------------------
///  relief map display
///http://www.colorschemer.com/schemes/viewscheme.php?id=81
class colorMapElevation: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapElevation():
//        QwtLinearColorMap( QColor("#B17142"),QColor("#FFDFC7"))
//              QwtLinearColorMap( QColor("#8d5524"),QColor("#ffdbac").lighter())
//                          QwtLinearColorMap( QColor("#D7191C"),QColor("#ffffbf"))
    QwtLinearColorMap( QColor(141,116,94),QColor(255,251,244))
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



    }
};//---------------------------------------------------------------------------
/// Gray scale legend for shaded relief map display
class colorMapGray: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapGray():
        QwtLinearColorMap( QColor("#555555"),QColor("#ffffff"))
    {
        addColorStop(0.000,QColor("#555555"));
    }
};
//---------------------------------------------------------------------------
/// Dark yellow legend for maps for road map overlay
class colorMapRoads: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapRoads():
        QwtLinearColorMap( QColor("#277e00"), QColor("#277e00"))//00cc00"))//888800"))
    {
        addColorStop(0.0, QColor("#277e00"));//QColor("#cc8800"));
    }
};
//---------------------------------------------------------------------------
/// Green legend for maps for road map overlay
class colorMapRoads2: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapRoads2():
        QwtLinearColorMap( Qt::yellow, QColor("#eecc00")  )
    {
        addColorStop( 0.0, Qt::yellow );
    }
};
//---------------------------------------------------------------------------
/// Logarithmic Yellow to blue legend for runoff map display
class colorMapWaterLog: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );
        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }

public:
//    colorMapWaterLog():
//        QwtLinearColorMap( QColor("#f6f633"), QColor("#000080"))
//    {
//        addColorStop( 0.0, QColor("#f6f633").lighter(125));//Qt::yellow );
//        addColorStop( 0.003,QColor("#8080FF"));
//        addColorStop( 0.03, QColor("#4040ff") );
//        addColorStop( 0.5, QColor("#0000FF"));
//    }
  colorMapWaterLog():
//      QwtLinearColorMap( QColor("#f6f666"), QColor("#ff3300"))
    QwtLinearColorMap( QColor("#8c8cff"), QColor("#ff3300"))

  {
     // addColorStop( 0.0, QColor("#f6f633"));
      addColorStop( 0.0005,QColor("#8080FF"));
      addColorStop( 0.01, QColor("#4040ff") );
      addColorStop( 0.05, QColor("#0000FF"));
      addColorStop( 0.1, QColor("#00006F"));
      addColorStop( 0.9, QColor("#FF0000"));
  }
};
//---------------------------------------------------------------------------
/// Linear Yellow to blue legend for infil map display
class colorMapWater: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );
        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapWater():
        QwtLinearColorMap( QColor("#f6f666"),QColor("#0000AA"))
    {
     //   addColorStop( 0.0, QColor("#FFFF00").lighter(125) );
        addColorStop( 0.1, QColor("#FFFF55") );
        addColorStop( 0.4, QColor("#8080FF") );
        addColorStop( 0.9, Qt::blue );
    }
};
//---------------------------------------------------------------------------
/// Transparent  light to dark blue legend for flood display
class colorMapFlood: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value <= thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFlood():
        QwtLinearColorMap( QColor(Qt::blue).lighter(150) ,  QColor(Qt::blue).darker(300))
    {
        addColorStop(0.00,QColor(Qt::blue).lighter(150));
        addColorStop(0.500,Qt::blue);
        addColorStop(1.0,QColor(Qt::blue).darker(300));
    }
};
//---------------------------------------------------------------------------
/// Cyan to red legend for sediment display
class colorMapSed: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value < 0.01 && value > -0.01)//thresholdLCM )
            return qRgba( 0, 0, 0, 0 );
        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapSed():
        QwtLinearColorMap( Qt::darkCyan,Qt::red)
    {
        addColorStop( 0.0, Qt::darkCyan );
        addColorStop( 0.3, Qt::cyan );
        addColorStop( 0.5, QColor("#fee691"));//Qt::white );
        addColorStop( 0.7, Qt::yellow);
    }
};
//---------------------------------------------------------------------------
/// Blue to red legend for sediment display
class colorMapSedB: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

//        if ( value < thresholdLCM )
//            return qRgba( 0, 0, 0, 0 );
        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapSedB():
        QwtLinearColorMap( Qt::blue,Qt::red)
    {
        addColorStop( 0.0, Qt::blue );
        addColorStop( 0.3, Qt::cyan );
        addColorStop( 0.5, Qt::white );
        addColorStop( 0.7, Qt::yellow);
    }
};
//---------------------------------------------------------------------------
class colorMapFloodV: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value < thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFloodV():
        QwtLinearColorMap( QColor("#009900"),QColor("#BF0000"))
    {
        addColorStop( 0.0, QColor("#009900"));
        addColorStop( 0.25, Qt::yellow);
        addColorStop( 0.75, Qt::red);
    }
};
//---------------------------------------------------------------------------
class colorMapP: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value < thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapP():
        QwtLinearColorMap( QColor(BGc),Qt::red)
    {
        addColorStop(0.0,QColor("#8888FF"));
        addColorStop( 0.5, Qt::blue);
    }
};
//---------------------------------------------------------------------------
class colorMapFEW: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value < thresholdLCM )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFEW():
        QwtLinearColorMap( Qt::darkRed,Qt::darkGreen)
    {
        addColorStop( 0.25, Qt::red);
        addColorStop( 0.5, Qt::yellow);
        addColorStop( 0.75, Qt::green);
    }
};
//---------------------------------------------------------------------------



#endif // LISUIMAPCOLOR_H

