/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file lookup.cpp
  \brief SWATRE: computes Theta from head, K from head or Diff Moist Cap from head.

  functions:\n
- double HNode(double theta, const  HORIZON *hor)\n
- double TheNode(double head, const  HORIZON *hor)\n
- double HcoNode(double head, const HORIZON *hor, double calib, double SEC)\n
- double DmcNode(double head, const  HORIZON *hor) \n

*/

#include <algorithm>
#include "swatre_p.h"
#include "swatresoillut.h"
#include "swatreLookup.h"

//-----------------------------------------------------------------------------------
/// head from theta
double HNode(
        double theta,           // current theta value of this node
        const  HORIZON *hor)    // parameters of horizon this node belongs to
{
    return LUT_LinIntPol(hor->lut,H_COL, theta,THETA_COL);
}
//-----------------------------------------------------------------------------------
/// theta from head
double TheNode(
        double head,           // current head value of this node
        const  HORIZON *hor)   // parameters of horizon this node belongs to
{
    //head = std::min( head, -1e-10);  //max or min ???? org was max! but head is < 0 !
    //VJ 110825 better to comment out this line, not useful
    if (head >= -1.0E-2)
       return LUT_Highest(hor->lut, THETA_COL);
    return LUT_LinIntPol(hor->lut, THETA_COL, head, H_COL);
}
//-----------------------------------------------------------------------------------
/// hydraulic conductivity from head
double HcoNode(
        double head,
        const HORIZON *hor,
        double calib)
{
    if (head >= -1.0E-2)
        return (LUT_Highest(hor->lut, K_COL)*calib/86400.0);
    // table is in cm/day, funcion return cm/sec
    return (LUT_LinIntPol(hor->lut, K_COL, head, H_COL)/86400.0);
}
//-----------------------------------------------------------------------------------
/// Differential Moisture Capacity from head
double DmcNode(
        double head,           // current head value of this node
        const  HORIZON *hor)   // parameters of horizon this node belongs to
{
    int i;         // index in LUT where dmch[i] <= head <= dmch[i+1]
    const LUT *l;  // lut of this horizon

    // dit gaat niet goed als profiel van verzadigd naar onversazigd switht:
    //if (head >= 0) return 0;

    if (head >= -1.0E-2)
       return LUT_Highest(hor->lut, DMCC_COL);
    //       return LUT_LinIntPol(hor->lut, DMCC_COL, head, DMCH_COL);

    l = hor->lut;
    i = LUT_Index_LE(l, head, DMCH_COL);
    i = std::min(LUT_nrRows(l)-2, i);
    i = std::max(i, 0);

    return LUT_ValueAt(l, DMCC_COL, i) +
            (head - LUT_ValueAt(l, DMCH_COL, i)) *
            (LUT_ValueAt(l,DMCC_COL, i+1)-LUT_ValueAt(l,DMCC_COL, i))/
            (LUT_ValueAt(l,DMCH_COL, i+1)-LUT_ValueAt(l,DMCH_COL, i));
}
//-----------------------------------------------------------------------------------
