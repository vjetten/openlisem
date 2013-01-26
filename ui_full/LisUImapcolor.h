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
    }
};

class colorMapWaterLog: public QwtLinearColorMap
{
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
public:
    colorMapWater():
        QwtLinearColorMap( QColor(BGc), Qt::darkBlue  )
    {
//        addColorStop( 0.0, qRgba(0, 0, 0, 0));//Qt::transparent);
//        addColorStop( 0.0, qRgba( 191, 191, 0, 0));
        addColorStop( 0.0, Qt::yellow );
        addColorStop( 0.1, QColor("#FFFF55") );
        addColorStop( 0.4, QColor("#8080FF") );
        addColorStop( 0.9, Qt::blue );
    }
};

//http://www.color-hex.com/color/8080ff
class colorMapFlood: public QwtLinearColorMap
{
public:
    colorMapFlood():
        QwtLinearColorMap( QColor(BGc), QColor("#000080"))
    {
//        addColorStop( 0.0, qRgba( 191, 191, 254, 0));//Qt::transparent);
        addColorStop( 0.0, QColor("#bfbffe"));
        addColorStop( 0.2, QColor("#8080ff"));
        addColorStop( 0.6, QColor("#0000ff"));
    }
};

class colorMapSed: public QwtLinearColorMap
{
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
