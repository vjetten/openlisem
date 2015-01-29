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


#include "CsfMap.h"
#include "error.h"
#include <QtGui>

//---------------------------------------------------------------------------
cTMap::cTMap()
{
  Created = false;
  nrRows = 0;
  nrCols = 0;
  MapName = "";
  PathName = "";
}
//---------------------------------------------------------------------------
cTMap::~cTMap()
{
  KillMap();
}
//---------------------------------------------------------------------------
void cTMap::KillMap()
{
  Data = MaskedRaster<REAL8>();
  Created = false;
}
//---------------------------------------------------------------------------
void cTMap::GetMapHeader(QString Name)
{
  MAP *m = Mopen(Name.toAscii().constData(), M_READ);
  if (m == NULL)
    Error(QString("Map %1 cannot be opened.").arg(Name));
  //   {
  //      Error(QString("Map %1 cannot be opened.").arg(Name));
  //      throw 1;
  //   }

  MH = m->raster;
  projection = m->main.projection;

  // easier to separate these
  nrRows = (int)MH.nrRows;
  nrCols = (int)MH.nrCols;

  Mclose(m);
}
//---------------------------------------------------------------------------
/// make an empty map structure
void cTMap::CreateMap(QString Name)
{
  // if exist delete main map data structure and set Created to false
  if (Created)
    KillMap();

  // now get the header for nrrows and nrcols
  GetMapHeader(Name);

  Data = MaskedRaster<REAL8>(nrRows, nrCols);

  Created = true;
}
//---------------------------------------------------------------------------
// make an empty a map and load it from disk
bool cTMap::LoadFromFile()
{
  MAP *m;
  QFileInfo fi(PathName);

  if (!fi.exists())
    Error(QString("Map %1 does not exist.").arg(PathName));
  //   {
  //      Error(QString("Map %1 does not exist.").arg(PathName));
  //      throw 1;
  //   }

  // make map structure
  CreateMap(PathName);

  if (!Created)
    return(false);

  MapName = PathName;

  m = Mopen(MapName.toAscii().constData(), M_READ);

  if (!m)
    return(false);

  RuseAs(m, CR_REAL8); //RgetCellRepr(m));
  for(int r=0; r < nrRows; r++)
    RgetSomeCells(m, (UINT4)r*nrCols, (UINT4)nrCols, Data[r]);

  if (RgetCellRepr(m) == CR_REAL8)
    ResetMinMax();

  Mclose(m);
  return(true);
}
//---------------------------------------------------------------------------
// make a new map according to dup as a mask and filled with value
void cTMap::MakeMap(cTMap *dup, REAL8 value)
{
  if (dup == NULL)
    return;
  //MH = dup->MH;
  MH.minVal = value;
  MH.maxVal = value;
  nrRows = dup->nrRows;
  nrCols = dup->nrCols;
  MH.nrRows = nrRows;
  MH.nrCols = nrCols;
  MH.valueScale = dup->MH.valueScale;
  MH.cellRepr = dup->MH.cellRepr;
  MH.xUL = dup->MH.xUL;
  MH.yUL = dup->MH.yUL;
  MH.cellSize  = dup->MH.cellSize;
  //MH.cellSizeY  = dup->MH.cellSizeX;
  MH.angle = dup->MH.angle;
  projection = dup->projection;

  Data = MaskedRaster<REAL8>(nrRows, nrCols);

  for(int r = 0; r < nrRows; r++)
    SetMemMV(Data[r],nrCols,CR_REAL8);

  for(int r=0; r < nrRows; r++)
    for(int c=0; c < nrCols; c++)
      if (!IS_MV_REAL8(&dup->Data[r][c]))
        {
          Data[r][c] = value;
        }

  Created = true;
}
//---------------------------------------------------------------------------
void cTMap::ResetMinMax(void)
{
  REAL8 minv = 1e30, maxv = -1e30;

  if (!Created)
    return;

  for(int r=0; r < nrRows; r++)
    for(int c=0; c < nrCols; c++)
      if (!IS_MV_REAL8(&Data[r][c]))
        {
          if (maxv < Data[r][c]) maxv = Data[r][c];
          if (minv > Data[r][c]) minv = Data[r][c];
        }

  MH.minVal = minv;
  MH.maxVal = maxv;
}
//---------------------------------------------------------------------------
/// write a map to disk
void cTMap::WriteMap(QString Name)
{
  MAP *out;
  long r, c;
  REAL4 *Dt;

  if (!Created)
    return;

  if (Name.isEmpty())
    {
      ErrorString = "Cannot write file, file name empty";
      throw 1;
    }

  ResetMinMax();

  Dt = new REAL4[nrCols];
  // make an array for output

  MH.cellRepr = CR_REAL4;
  out = Rcreate(Name.toAscii().constData(),nrRows, nrCols, (CSF_CR)MH.cellRepr, VS_SCALAR,
                (CSF_PT)projection, MH.xUL, MH.yUL, MH.angle, MH.cellSize);
  RuseAs(out, CR_REAL4);

  for(r=0; r < nrRows; r++)
    {
      for(c=0; c < nrCols; c++)
        Dt[c] = (REAL4)Data[r][c];


      if (RputRow(out, r, Dt) != (UINT4)nrCols)
        {
          ErrorString = "rputrow write error with" + Name;
          throw 1;
        }
    }

  delete Dt;

  Mclose(out);

}
//---------------------------------------------------------------------------
/// makes mapname if (name.map) or mapseries (name0000.001 to name0009.999)
void cTMap::WriteMapSeries(QString Dir, QString Name, int count)
{
  QString path;
  QFileInfo fi(Name);

  if (Name.indexOf(".") < 0)
    {
      QString nam, dig;

      nam = Name + "00000000";

      nam.remove(7, 10);
      dig = QString("%1").arg(count,4,10,QLatin1Char('0'));
      dig.insert(1,".");
      Name = nam + dig;

    }
  path = Dir + Name;
  WriteMap(path);
}
//---------------------------------------------------------------------------


