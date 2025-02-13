/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#pragma once

#include <QString>


#define ADD 0
#define SUB 1
#define MUL 2
#define DIV 3
#define POW 4
#define MIN 5  //VJ 041120 added this functionality
#define MAX 6
#define LARGER 7
#define SMALLER 8
#define LARGEREQUAL 9
#define SMALLEREQUAL 10
#define HIGHER 11
#define LOWER 12


class cTMap;

void               copy                (cTMap& raster,
                                        cTMap const& other);

QList <int>        countUnits          (cTMap const& raster);

void               fill                (cTMap& raster,
                                        double value);

double             mapTotal            (cTMap const& raster);

double             mapAverage          (cTMap const& raster);

double             mapMinimum          (cTMap const& raster);

double             mapMaximum          (cTMap const& raster);

double             getWindowAverage    (cTMap const& raster,
                                        int r,
                                        int c,
                                        bool center);

void               cover               (cTMap& raster,
                                        cTMap const& value1,
                                        double value2);

void               calcValue           (cTMap& raster,
                                        double value,
                                        int oper);

void               calcMap             (cTMap& raster,
                                        cTMap const& value,
                                        int oper);

void               calc2Maps           (cTMap& raster,
                                        cTMap const& value1,
                                        cTMap const& value2,
                                        int oper);

void               calcMapValue        (cTMap& raster,
                                        cTMap const& value1,
                                        double value2,
                                        int oper);

void               checkMap            (cTMap const& raster,
                                        int oper,
                                        double value,
                                        QString SS);
