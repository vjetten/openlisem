#include "csfrgbmap.h"


cTRGBMap::cTRGBMap()
{



}


cTRGBMap::cTRGBMap(
    MaskedRaster<char>&& dataR,
    QString const& projection,
    QString const& mapName)

    : dataR(std::forward<MaskedRaster<char>>(dataR)),
      _projection(projection),
      _mapName(mapName)

{
    bands = 1;
}

cTRGBMap::cTRGBMap(
    MaskedRaster<char>&& dataR,
    MaskedRaster<char>&& dataG,
    QString const& projection,
    QString const& mapName)

    : dataR(std::forward<MaskedRaster<char>>(dataR)),
      dataG(std::forward<MaskedRaster<char>>(dataG)),
      _projection(projection),
      _mapName(mapName)

{
    bands = 2;
}

cTRGBMap::cTRGBMap(
    MaskedRaster<char>&& dataR,
    MaskedRaster<char>&& dataG,
    MaskedRaster<char>&& dataB,
    QString const& projection,
    QString const& mapName)

    : dataR(std::forward<MaskedRaster<char>>(dataR)),
      dataG(std::forward<MaskedRaster<char>>(dataG)),
      dataB(std::forward<MaskedRaster<char>>(dataB)),
      _projection(projection),
      _mapName(mapName)

{
    bands = 3;
}


int cTRGBMap::nrRows() const
{
    return static_cast<int>(dataR.nr_rows());
}


int cTRGBMap::nrCols() const
{
    return static_cast<int>(dataR.nr_cols());
}


double cTRGBMap::north() const
{
    return dataR.north();
}


double cTRGBMap::west() const
{
    return dataR.west();
}


double cTRGBMap::cellSize() const
{
    return dataR.cell_size();
}


QString const& cTRGBMap::projection() const
{
    return _projection;
}


QString const& cTRGBMap::mapName() const
{
    return _mapName;
}
