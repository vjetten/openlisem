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
