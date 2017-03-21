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
 \file lisUnifiedFlowFluxLimiter.cpp
 \brief contains the flux limiting functions
 These functions are crucial in limiting the flux so that no oscillations
 emerge, while allowing for wave-like behavior and other large flow volumes

 MinMod works best on sloped surfaces, where pressure terms are relatively low
 In case of large scale flooding, use the HLL or HLL2 flux limiters

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"


double TWorld::FluxLimiter2D()
{



    return 0.0;
}


double TWorld::FluxLimiter1D()
{



return 0.0;
}


double TWorld::FL_MinMod()
{


    return 0.0;
}


double TWorld::FL_HLL2()
{


    return 0.0;
}
