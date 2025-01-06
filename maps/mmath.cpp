
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
  \file mmath.cpp
  \brief basic map mathematics on two maps or maps and variables

functions: \n
- void CTMap::fill(double value)   \n
- void CTMap::calcV(double v, int oper)   \n
- void CTMap::calc(cTMap *m, int oper)   \n
- void CTMap::calc2(cTMap *m1, cTMap *m2, int oper)   \n
- void CTMap::calc2V(cTMap *m1, double V, int oper)   \n
- void CTMap::copy(cTMap *m)   \n
- void CTMap::cover(double v)   \n
- void CTMap::setMV()   \n
- double CTMap::mapTotal()   \n
operations for 'oper' are ADD, SUB, MUL, DIV, POW, MIN, MAX
*/

#include <algorithm>
#include "model.h"
#include "operation.h"


