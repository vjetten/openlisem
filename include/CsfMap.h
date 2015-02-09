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

#include <QString>
#include "masked_raster.h"

//---------------------------------------------------------------------------
/// CSF map construction, reading, writing series etc.
/** class to deal with CSF map construction, reading and writing etc.
   Reading and writing of maps of other GIS systems can be added here in the future.
*/
class cTMap
{

public:

    //! The actual raster.
    MaskedRaster<double> Data;

                   cTMap               ()=default;

                   cTMap               (MaskedRaster<double>&& Data,
                                        QString const& mapName);

                   cTMap               (cTMap const& other)=delete;

                   cTMap               (cTMap&& other)=default;

                   ~cTMap              ()=default;

    cTMap&         operator=           (cTMap const& other)=delete;

    cTMap&         operator=           (cTMap&& other)=default;

    int            nrRows              () const;

    int            nrCols              () const;

    double         north               () const;

    double         west                () const;

    double         cellSize            () const;

    QString const& MapName             () const;

    void           setMV               ();

    void           MakeMap             (cTMap *dup,
                                        REAL8 value);

private:

    QString        _MapName;

};

#endif
