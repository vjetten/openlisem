#ifndef CSFRGBMAP_H
#define CSFRGBMAP_H

#include <QString>
#include "masked_raster.h"

class cTRGBMap
{
public:
    cTRGBMap();

    cTRGBMap (MaskedRaster<char>&& dataR,
           QString const& projection,
           QString const& mapName);

    cTRGBMap (MaskedRaster<char>&& dataR,
              MaskedRaster<char>&& dataG,
           QString const& projection,
           QString const& mapName);

    cTRGBMap (MaskedRaster<char>&& data,
              MaskedRaster<char>&& dataG,
              MaskedRaster<char>&& dataB,
           QString const& projection,
           QString const& mapName);

    //! The actual raster.
    MaskedRaster<char> dataR; //BAND 1

    //! The actual raster.
    MaskedRaster<char> dataG; //BAND 2

    //! The actual raster.
    MaskedRaster<char> dataB; //BAND 3

    int bands  = 0;

    int            nrRows              () const;

    int            nrCols              () const;

    double         north               () const;

    double         west                () const;

    double         cellSize            () const;

    QString const& projection          () const;

    QString const& mapName             () const;

private:

    //! Projection string as WKT string. Possibly empty.
    QString        _projection;

    QString        _mapName;



};

#endif // CSFRGBMAP_H
