#ifndef LISUIMAPCOLOR_H
#define LISUIMAPCOLOR_H

#include "lisemqt.h"

//---------------------------------------------------------------------------
#define BGc "#eeeeee" // background grey for missing value in maps

class alphaLinearColorMap: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value == 0 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
};

class colorMapGray: public QwtLinearColorMap
{
public:
    colorMapGray():
        QwtLinearColorMap( QColor(BGc), Qt::white  )
    {
        addColorStop(0, QColor("#111111"));
    }
};

class colorMapYellow: public QwtLinearColorMap
{
public:
    colorMapYellow():
        QwtLinearColorMap( QColor(BGc), Qt::white  )
    {
        addColorStop(0.0, Qt::darkYellow);
        addColorStop(0.5, Qt::yellow);
    }
};

class colorMapWaterLog: public QwtLinearColorMap
{
//    virtual QRgb rgb( const QwtInterval &interval, double value ) const
//    {
//        if ( value == 0 )
//            return qRgba( 0, 0, 0, 0 );

//        return QwtLinearColorMap::rgb( interval, value );
//    }
public:
    colorMapWaterLog():
        QwtLinearColorMap( QColor(BGc), QColor(0,0,128))//Qt::darkBlue )
    {
        addColorStop( 0.0, Qt::yellow );
        addColorStop( 0.05, QColor(128,128,255));
        addColorStop( 0.1, QColor(64,64,255) );
        addColorStop( 0.5, QColor(0,0,255));//Qt::blue );
    }
};

class colorMapWater: public QwtLinearColorMap
{
//    virtual QRgb rgb( const QwtInterval &interval, double value ) const
//    {
//        if ( value == 0 )
//            return qRgba( 0, 0, 0, 0 );

//        return QwtLinearColorMap::rgb( interval, value );
//    }
public:
    colorMapWater():
        QwtLinearColorMap( QColor(BGc), Qt::darkBlue  )
    {
        addColorStop( 0.0, Qt::yellow );
        addColorStop( 0.1, QColor("#FFFF55") );
        addColorStop( 0.4, QColor("#8080FF") );
        addColorStop( 0.9, Qt::blue );
    }
};

//http://www.color-hex.com/color/8080ff
class colorMapFlood: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value == 0 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFlood():
        QwtLinearColorMap( QColor(BGc), QColor("#0e0e4c"))
    {
        addColorStop( 0.000, QColor("#6666d2"));
        addColorStop( 0.125, QColor("#3a3ac5"));
        addColorStop( 0.250, QColor("#2525bf"));
        addColorStop( 0.375, QColor("#2121ab"));
        addColorStop( 0.500, QColor("#1d1d98"));
        addColorStop( 0.625, QColor("#191985"));
        addColorStop( 0.750, QColor("#161672"));
        addColorStop( 0.875, QColor("#12125f"));
    }
};

//http://www.color-hex.com/color/2525bf
//#ffffff
//#d3d3f2
//#bdbdeb
//#a7a7e5
//#9292df
//#7c7cd8
//#6666d2
//#5050cb

//addColorStop( 0.000, QColor("#3a3ac5"));
//addColorStop( 0.125, QColor("#2525bf"));
//addColorStop( 0.250, QColor("#2121ab"));
//addColorStop( 0.375, QColor("#1d1d98"));
//addColorStop( 0.500, QColor("#191985"));
//addColorStop( 0.625, QColor("#161672"));
//addColorStop( 0.750, QColor("#12125f"));

//#0e0e4c
//#0b0b39
//#070726
//#030313
//#000000


class colorMapSed: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value == 0 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapSed():
        QwtLinearColorMap( QColor(BGc),Qt::red)//QColor("#903000") )//QColor("#cc3000"));//Qt::darkYellow);
    {
        addColorStop( 0.0, Qt::darkCyan );//QColor("#108030"));
        addColorStop( 0.3, Qt::cyan );//QColor("#30ffcc"));
        addColorStop( 0.5, Qt::transparent);//Qt::white );
        addColorStop( 0.7, Qt::yellow);
    }
};

class colorMapSedB: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value == 0 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
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

#endif // LISUIMAPCOLOR_H
