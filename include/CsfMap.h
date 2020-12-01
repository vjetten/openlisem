/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
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
#ifndef CsfMapH
#define CsfMapH
#include <QString>
#include "masked_raster.h"


/*!
    @brief      A cTMap contains all relevant information about a raster.
    @todo       The data member must be made private.

    cTMap instances contain raster data, projection information and a map name.
    I/O of cTMap instances is handles by functions defined in the io module.
*/
class cTMap
{

public:

    //! The actual raster.
    MaskedRaster<double> data;

                   cTMap               ()=default;

                   cTMap               (MaskedRaster<double>&& data,
                                        QString const& projection,
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

    QString const& projection          () const;

    QString const& mapName             () const;

    void           setAllMV            ();

    void           MakeMap             (cTMap *dup,
                                        REAL8 value);

private:

    //! Projection string as WKT string. Possibly empty.
    QString        _projection;

    QString        _mapName;

};

#endif
