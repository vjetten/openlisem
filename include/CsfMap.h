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
  \file CsfMap.h
  \brief file operations  class for PCRaster maps.
  */


#ifndef CsfMapH
#define CsfMapH

#include "csf.h"

#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))

//---------------------------------------------------------------------------

/** class to deal with CSF map construction, reading and writing etc.
   Reading and writing of maps of other GIS systems can be added here in the future.
*/
class cTMap
{
protected:

public:
    CSF_RASTER_HEADER MH; ///PCRaster map header
    UINT2 projection;
    REAL8 **Data;

    QString MapName;
    QString PathName;
    int nrRows, nrCols;
    bool Created;
    
    void KillMap();
    void GetMapHeader(QString Name);
    void CreateMap(QString Name);
    void MakeMap(cTMap *dup, REAL8 value);
    void WriteMap(QString Name);
    void WriteMapSeries(QString Dir, QString Name, int count);
    bool LoadFromFile();
    void ResetMinMax(void);


    cTMap();
    ~cTMap();
};

#endif

/* for info:
typedef struct CSF_RASTER_HEADER
{
	 UINT2    valueScale;
	 UINT2    cellRepr;
	 CSF_VAR_TYPE minVal;
	 CSF_VAR_TYPE maxVal;
	 REAL8    xUL;
	 REAL8    yUL;
	 UINT4    nrRows;
	 UINT4    nrCols;
	 REAL8    cellSizeX;
	 REAL8    cellSizeY;
	 REAL8    angle;
} CSF_RASTER_HEADER;
*/

