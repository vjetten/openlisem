#ifndef LISUIMAPCOLOR_H
#define LISUIMAPCOLOR_H

#include "lisemqt.h"

//---------------------------------------------------------------------------
#define BGc "#eeeeee" // background grey for missing value in maps

//http://www.color-hex.com/color/


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

        QwtPlotItemList list = plot()->itemList(QwtPlotItem::Rtti_PlotSpectrogram);
        QwtPlotSpectrogram * sp = static_cast<QwtPlotSpectrogram *> (list.at(1));
        QwtPlotSpectrogram * sp0 = static_cast<QwtPlotSpectrogram *> (list.at(list.count()-1));
        // elevation info

        if (sp->data() == NULL)
            return QwtText(txt);
        double z = sp->data()->value(pos.x(), pos.y());
        double z0 = sp0->data()->value(pos.x(), pos.y());

        //        if (checkUnits_tonha->isChecked()) unit = "ton/ha";// ton/ha
        //        if (checkUnits_kgcell->isChecked()) unit = "kg/cell"; // in kg/cell
        //        if (checkUnits_kgm2->isChecked()) unit = "kg/m2"; // in kg/m2
        if (z > -1e10)
        {
            if (sp->data()->value(0,0) == 1) txt = QString("%1 l/s [%2m]").arg(z,0,'f',1).arg(z0,0,'f',1);
            if (sp->data()->value(0,0) == 2) txt = QString("%1 mm [%2m]").arg(z,0,'f',1).arg(z0,0,'f',1);
            if (sp->data()->value(0,0) == 3)
                txt = QString("%1 %2[%3m]").arg(z,0,'f',1).arg(unit).arg(z0,0,'f',1);
            if (sp->data()->value(0,0) == 4) txt = QString("%1 m [%2m]").arg(z,0,'f',2).arg(z0,0,'f',1);
            if (sp->data()->value(0,0) == 5) txt = QString("%1 m/s [%2m]").arg(z,0,'f',2).arg(z0,0,'f',1);
            if (sp->data()->value(0,0) == 6) txt = QString("%1 mm [%2m]").arg(z,0,'f',3).arg(z0,0,'f',1);
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

        if ( value < 0.1 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapHouse():
        QwtLinearColorMap( QColor(BGc), QColor("#331900"))
    {
        addColorStop(0, QColor("#cc6600"));
    }
};
//---------------------------------------------------------------------------
/// Gray scale legend for shaded relief map display
class colorMapGray: public QwtLinearColorMap
{
public:
    colorMapGray():
        QwtLinearColorMap( QColor(BGc),Qt::white  )
    {
        addColorStop(0, QColor("#111111"));
    }
};
//---------------------------------------------------------------------------
class colorMapWhite: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < -1e19 )
            return qRgba( 228, 228, 228, 255 );

        if ( value < 0.1 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapWhite():
        QwtLinearColorMap( QColor(BGc),QColor("#111111")  )
    {
        addColorStop(0, Qt::white );
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
        QwtLinearColorMap( QColor(BGc), QColor("#277e00"))//00cc00"))//888800"))
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
        QwtLinearColorMap( QColor(BGc), QColor("#eecc00")  )
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
    colorMapWaterLog():
        QwtLinearColorMap( QColor(BGc), QColor("#000080"))
    {
        addColorStop( 0.0, QColor("#f6f633"));//Qt::yellow );
        addColorStop( 0.003,QColor("#8080FF"));
        addColorStop( 0.03, QColor("#4040ff") );
        addColorStop( 0.5, QColor("#0000FF"));
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
        QwtLinearColorMap( QColor(BGc),QColor("#0000AA"))
    {
        addColorStop( 0.0, Qt::yellow );
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
        addColorStop(1.0,QColor(Qt::blue).darker(400));
    }
};
//---------------------------------------------------------------------------
/// Cyan to red legend for sediment display
class colorMapSed: public QwtLinearColorMap
{
    //    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    //    {
    //        if ( value > -thresholdCLM && value < thresholdCLM )
    //            return qRgba( 0, 0, 0, 0 );

    //        return QwtLinearColorMap::rgb( interval, value );
    //    }
public:
    colorMapSed():
        QwtLinearColorMap( QColor(BGc),Qt::red)//QColor("#903000") )//QColor("#cc3000"));//Qt::darkYellow);
    {
        addColorStop( 0.0, Qt::darkCyan );//QColor("#108030"));
        addColorStop( 0.3, Qt::cyan );//QColor("#30ffcc"));
        addColorStop( 0.5, Qt::white );
        addColorStop( 0.7, Qt::yellow);
    }
};
//---------------------------------------------------------------------------
/// Blue to red legend for sediment display
class colorMapSedB: public QwtLinearColorMap
{
    //    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    //    {
    //        if ( value == 0)//< 1 && value > -1 )
    //            return qRgba( 0, 0, 0, 0 );

    //        return QwtLinearColorMap::rgb( interval, value );
    //    }
public:
    colorMapSedB():
        QwtLinearColorMap( QColor(BGc),Qt::red)
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
        QwtLinearColorMap( QColor(BGc),QColor("#BF0000"))
    {
        addColorStop( 0.0, QColor("#006600"));
        addColorStop( 0.4, Qt::yellow);
        addColorStop( 0.8, Qt::red);
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



#endif // LISUIMAPCOLOR_H

