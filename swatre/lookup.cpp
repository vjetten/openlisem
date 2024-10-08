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
// #include "swatre_p.h"
// #include "swatresoillut.h"
// #include "swatreLookup.h"

#include "model.h"


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
    LUT *l = hor->lut;

    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);
    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[THETA_COL][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[THETA_COL][l->nrRows-1]);
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);

        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        double lTh = l->hydro[THETA_COL][lowerIndex];
        double uTh = l->hydro[THETA_COL][upperIndex];
        return (lTh + f*(uTh-lTh));
    }

    //head = std::min( head, -1e-10);  //max or min ???? org was max! but head is < 0 !
    //VJ 110825 better to comment out this line, not useful
    // if (head >= -1.0E-2)
    //    return LUT_Highest(hor->lut, THETA_COL);
    // return LUT_LinIntPol(hor->lut, THETA_COL, head, H_COL);
}
//-----------------------------------------------------------------------------------
/// hydraulic conductivity from head
double HcoNode(double head,const HORIZON *hor,double calib)
{
    LUT *l = hor->lut;

    auto it = std::lower_bound(l->hydro[H_COL].begin(), l->hydro[H_COL].end(), head);
    if (it == l->hydro[H_COL].begin()) {
        return(l->hydro[K_COL][0]);
    } else if (it == l->hydro[H_COL].end()) {
        return(l->hydro[K_COL][l->nrRows-1]);
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);
        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        double lK = l->hydro[K_COL][lowerIndex];
        double uK = l->hydro[K_COL][upperIndex];
       // qDebug() << head << lH << uH << lK << uK << (lK+f*(uK-lK));
        return (lK+f*(uK-lK));
    }
}
//-----------------------------------------------------------------------------------
/// Differential Moisture Capacity from head
double DmcNode(
        double head,           // current head value of this node
        const  HORIZON *hor)   // parameters of horizon this node belongs to
{   
    // dit gaat niet goed als profiel van verzadigd naar onverzadigd switht:
    //if (head >= 0) return 0;

    LUT *l = hor->lut;

    auto it = std::lower_bound(l->hydro[DMCH_COL].begin(), l->hydro[DMCH_COL].end(), head);

    if (it == l->hydro[DMCH_COL].begin()) {
        return(l->hydro[DMCC_COL][0]);
    } else if (it == l->hydro[DMCH_COL].end()) {
        return(l->hydro[DMCC_COL][l->nrRows-1]);
    } else {
        double lH = *(it - 1);
        double uH = *it;
        double f = (head-lH)/(uH-lH);
        int lowerIndex = std::distance(l->hydro[H_COL].begin(), it - 1); // Index of lower bound
        int upperIndex = std::distance(l->hydro[H_COL].begin(), it);     // Index of upper bound

        double lcH = l->hydro[DMCH_COL][lowerIndex];
        double ucH = l->hydro[DMCH_COL][upperIndex];
        double lC = l->hydro[DMCC_COL][lowerIndex];
        double uC = l->hydro[DMCC_COL][upperIndex];

        return (lC + (head-lcH)*(uC-lC)/(ucH-lcH));
    }

/*
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
            */
}
//-----------------------------------------------------------------------------------
