/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/
#include "CsfRGBMap.h"


cTRGBMap::cTRGBMap()
{



}


cTRGBMap::cTRGBMap(
    MaskedRaster<char>&& dataR,
//    MaskedRaster<char>&& dataA,
    QString const& projection,
    QString const& mapName)

    : dataR(std::forward<MaskedRaster<char>>(dataR)),
      //dataA(std::forward<MaskedRaster<char>>(dataA)),
      _projection(projection),
      _mapName(mapName)

{
    bands = 1;
}

cTRGBMap::cTRGBMap(
    MaskedRaster<char>&& dataR,
    MaskedRaster<char>&& dataG,
   // MaskedRaster<char>&& dataA,
    QString const& projection,
    QString const& mapName)

    : dataR(std::forward<MaskedRaster<char>>(dataR)),
      dataG(std::forward<MaskedRaster<char>>(dataG)),
    //  dataA(std::forward<MaskedRaster<char>>(dataA)),
      _projection(projection),
      _mapName(mapName)

{
    bands = 2;
}

cTRGBMap::cTRGBMap(
    MaskedRaster<char>&& dataR,
    MaskedRaster<char>&& dataG,
    MaskedRaster<char>&& dataB,
  //  MaskedRaster<char>&& dataA,
    QString const& projection,
    QString const& mapName)

    : dataR(std::forward<MaskedRaster<char>>(dataR)),
      dataG(std::forward<MaskedRaster<char>>(dataG)),
      dataB(std::forward<MaskedRaster<char>>(dataB)),
  //    dataA(std::forward<MaskedRaster<char>>(dataA)),
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
