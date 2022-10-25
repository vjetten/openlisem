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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/
#ifndef CSFRGBMAP_H
#define CSFRGBMAP_H

#include <QString>
#include "masked_raster.h"

class cTRGBMap
{
public:
    cTRGBMap();

    cTRGBMap (MaskedRaster<char>&& dataR,
         //     MaskedRaster<char>&& dataA,
           QString const& projection,
           QString const& mapName);

    cTRGBMap (MaskedRaster<char>&& dataR,
              MaskedRaster<char>&& dataG,
       //       MaskedRaster<char>&& dataA,
           QString const& projection,
           QString const& mapName);

    cTRGBMap (MaskedRaster<char>&& dataR,
              MaskedRaster<char>&& dataG,
              MaskedRaster<char>&& dataB,
       //       MaskedRaster<char>&& dataA,
           QString const& projection,
           QString const& mapName);

    //! The actual raster.
    MaskedRaster<char> dataR; //BAND 1

    //! The actual raster.
    MaskedRaster<char> dataG; //BAND 2

    //! The actual raster.
    MaskedRaster<char> dataB; //BAND 3

//    MaskedRaster<char> dataA; //shade

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
