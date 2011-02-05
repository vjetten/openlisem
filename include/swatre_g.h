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
  \file swatre_g.h
  \brief SWATRE global structure for soil profile information and swatre options

SWATRE: structure for soil profile information and swatre options.\n
This model is declared for normal soils, crustd soils, compacted soils and grass strips. \n
*/

#ifndef SWATRE_G_H
#define SWATRE_G_H

/// SWATRE: structure for soil profile information and swatre options.
typedef struct SOIL_MODEL {
   struct PIXEL_INFO  *pixel; //defined in swatre_p.h
  // double precision;
   double minDt;
  // double calibrationfactor;
  // bool geometric;
  // bool swatreBottomClosed;
  // long nrCells;
} SOIL_MODEL;
//---------------------------------------------------------------------------

#endif // SWATRE_G_H
