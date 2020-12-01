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
#include <cassert>
#include "CsfMap.h"
#include "error.h"
#include <gdal.h>
/*!
    @brief      Constructor.
    @param      data Properties of the raster and the cell values.
    @param      projection Projection as a WKT string or an empty string.
    @param      mapName Map name.
*/
cTMap::cTMap(
    MaskedRaster<double>&& data,
    QString const& projection,
    QString const& mapName)

    : data(std::forward<MaskedRaster<double>>(data)),
      _projection(projection),
      _mapName(mapName)

{
}


int cTMap::nrRows() const
{
    return static_cast<int>(data.nr_rows());
}


int cTMap::nrCols() const
{
    return static_cast<int>(data.nr_cols());
}


double cTMap::north() const
{
    return data.north();
}


double cTMap::west() const
{
    return data.west();
}


double cTMap::cellSize() const
{
    return data.cell_size();
}


QString const& cTMap::projection() const
{
    return _projection;
}


QString const& cTMap::mapName() const
{
    return _mapName;
}


void cTMap::setAllMV()
{
    data.set_all_mv();
}


// make a new map according to dup as a mask and filled with value
void cTMap::MakeMap(
    cTMap *dup,
    REAL8 value)
{
  if (dup == nullptr)
    return;

  data = MaskedRaster<REAL8>(dup->nrRows(), dup->nrCols(), dup->north(),
      dup->west(), dup->cellSize());

  data.set_all_mv();

  for(int r=0; r < nrRows(); r++)
    for(int c=0; c < nrCols(); c++)
      if (!pcr::isMV(dup->data[r][c]))
        {
          data[r][c] = value;
        }
}
