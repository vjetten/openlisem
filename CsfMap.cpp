/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
 \file CsfMap.cpp
 \brief file operations  class for PCRaster maps.

  provide basic functionality to read and write PCRaster CSF maps, \n
  can be altered to link to other file formats.

 functions: \n
 - void cTMap::KillMap(void)  \n
 - void cTMap::GetMapHeader(QString Name)  \n
 - void cTMap::CreateMap(QString Name) \n
 - bool cTMap::LoadFromFile() \n
 - void cTMap::MakeMap(cTMap *dup, REAL8 value) \n
 - void cTMap::ResetMinMax(void) \n
 - void cTMap::WriteMap(QString Name) \n
 - void cTMap::WriteMapSeries(QString Dir, QString Name, int count) \n
*/

#include <memory>
#include <cassert>
#include "CsfMap.h"
#include "error.h"


cTMap::cTMap(
    MaskedRaster<double>&& Data,
    QString const& mapName)

    : Data(std::forward<MaskedRaster<double>>(Data)),
      _MapName(mapName)

{
}


int cTMap::nrRows() const
{
    return static_cast<int>(Data.nr_rows());
}


int cTMap::nrCols() const
{
    return static_cast<int>(Data.nr_cols());
}


double cTMap::north() const
{
    return Data.north();
}


double cTMap::west() const
{
    return Data.west();
}


double cTMap::cellSize() const
{
    return Data.cell_size();
}


QString const& cTMap::MapName() const
{
    return _MapName;
}


void cTMap::setMV()
{
    Data.set_all_mv();
}


// make a new map according to dup as a mask and filled with value
void cTMap::MakeMap(
    cTMap *dup,
    REAL8 value)
{
  if (dup == NULL)
    return;

  Data = MaskedRaster<REAL8>(dup->nrRows(), dup->nrCols(), dup->north(),
      dup->west(), dup->cellSize());

  Data.set_all_mv();

  for(int r=0; r < nrRows(); r++)
    for(int c=0; c < nrCols(); c++)
      if (!pcr::isMV(dup->Data[r][c]))
        {
          Data[r][c] = value;
        }
}
