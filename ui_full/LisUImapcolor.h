#ifndef LISUIMAPCOLOR_H
#define LISUIMAPCOLOR_H

#include "lisemqt.h"

//---------------------------------------------------------------------------
#define BGc "#eeeeee" // background grey for missing value in maps


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
        QwtLinearColorMap( QColor(BGc), QColor("#0000AA"))
    {
        addColorStop( 0.0, Qt::yellow );
        addColorStop( 0.05, QColor(128,128,255));
        addColorStop( 0.1, QColor(64,64,255) );
        addColorStop( 0.5, QColor("#0000FF"));
    }
};

class colorMapWater: public QwtLinearColorMap
{
public:
    colorMapWater():
        QwtLinearColorMap( QColor(BGc),QColor("#0000AA"))
    {
        addColorStop( 0.0, Qt::yellow );
////        addColorStop(0.000,QColor("#6565FF"));
//        addColorStop(0.125,QColor("#4B4BFF"));
//        addColorStop(0.250,QColor("#3333FF"));
//        addColorStop(0.375,QColor("#1919FF"));
//        addColorStop(0.500,QColor("#0000FE"));
//        addColorStop(0.625,QColor("#0000E4"));
//        addColorStop(0.750,QColor("#0000CC"));
//        addColorStop(0.875,QColor("#0000B2"));
        addColorStop( 0.1, QColor("#FFFF55") );
        addColorStop( 0.4, QColor("#8080FF") );
        addColorStop( 0.9, Qt::blue );
    }
};

class colorMapGreen: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value == 0 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapGreen():
        QwtLinearColorMap( QColor(BGc), QColor("#008800")  )
    {
        addColorStop( 0.0, Qt::darkYellow );
    }
};

//http://www.color-hex.com/color/
class colorMapFlood: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value < 0.0001 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
public:
    colorMapFlood():
        QwtLinearColorMap( QColor(BGc),  QColor("#000098"))
    {
        addColorStop(0.000,QColor("#6565FF"));
        addColorStop(0.125,QColor("#4B4BFF"));
        addColorStop(0.250,QColor("#3333FF"));
        addColorStop(0.375,QColor("#1919FF"));
        addColorStop(0.500,QColor("#0000FE"));
        addColorStop(0.625,QColor("#0000E4"));
        addColorStop(0.750,QColor("#0000CC"));
        addColorStop(0.875,QColor("#0000B2"));
    }
};


class colorMapSed: public QwtLinearColorMap
{
    virtual QRgb rgb( const QwtInterval &interval, double value ) const
    {
        if ( value > -0.1 && value < 0.1 )
            return qRgba( 0, 0, 0, 0 );

        return QwtLinearColorMap::rgb( interval, value );
    }
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

#endif // LISUIMAPCOLOR_H

//class alphaLinearColorMap: public QwtLinearColorMap
//{
//    virtual QRgb rgb( const QwtInterval &interval, double value ) const
//    {
//        if ( value == 0 )
//            return qRgba( 0, 0, 0, 0 );

//        return QwtLinearColorMap::rgb( interval, value );
//    }
//};
