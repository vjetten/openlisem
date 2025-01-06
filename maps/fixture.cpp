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
#include "fixture.h"
#include <gdal_priv.h>


/*!
    @brief      Initializes the Lisem runtime environment.

    All GDAL I/O drivers are registered and GDAL is configured not to throw
    exceptions, but to return error codes.
*/
Fixture::Fixture()
{
    // Register all GDAL drivers.
    GDALAllRegister();

    // GDAL mustn't throw in case of an error.
    CPLSetErrorHandler(CPLQuietErrorHandler);
}


Fixture::~Fixture()
{
}
