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
  \file lisPesticide_MC.cpp
  \brief Transport and partitioning of pesticides with Euler forward method

functions: \n
- void TWorld::SimplePestCalc() \n
-
*/

#include "model.h"

// check if cell From flows to To
//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )


//---------------------------------------------------------------------------
/**
 * @fn double TWorld::simplePestCalc(double Crw_old, double Cmw_old, double Kfilm, double Qinf, double zm, double kr, double Kd,
                              double Crs_old, double Cms_old, double Ez, double Me, double A)
 * @brief Simple calculation of pesticide concentrations in a cell in both water and sediment
 * Simple calculation of pesticide concentrations in a cell in both water and sediment,
 * j = time and i = place: j1i1 is the new output, j1i is the new flux at the upstream 'entrance' flowing into the gridcell
 * @param Qj1i1 : result kin wave for this cell ;j = time, i = place
 * @param Qj1i : sum of all upstreamwater from kin wave
 * @param  Sj1i : sum of all upstream sediment
 * @param  dt : timestep
 * @param  vol : current volume of water in cell
 * @param  sed : current mass of sediment in cell
 * @return sediment outflow in next timestep
 *
 */
double TWorld::simplePestConc(double Crw_old, double Cmw_old, double Kfilm, double Qinf, double zm, double kr, double Kd,
                              double Crs_old, double Cms_old, double Ez, double Me, double A)
{
    double Crw_n = 0;
    double Cinf_n = 0;
    double Cmw_n = 0;
    double Cms_n = 0;
    double Crs_n = 0;
    double rho = 2650;
    FOR_ROW_COL_MV {
    double pore = ThetaS1 ->Drc;

    Crw_n = Crw_old + _dt * (Kfilm * (Cmw_old - Crw_old));
    Cmw_n = Cmw_old + _dt * ((Kfilm * (Crw_old - Cmw_old) + Qinf * (Crw_old - Cmw_old) - zm * rho * kr * (Kd * Cmw_old - Cms_old))/ pore * zm);
    }

    if (Ez > 0 ) {
        Cms_n = Cms_old + _dt * (kr * (Kd * Cmw_old - Cms_old));
        Crs_n = Crs_old + _dt * (((Crs_old * Me) + (Cms_old * Ez * rho * A))/(Me + (Ez * rho * A)));
    } else {
        Cms_n = Cms_old + _dt * (kr * (Kd * Cmw_old - Cms_old) + (((zm + Ez * Cms_old) - Crs_old * Ez)/ zm));
        Crs_n = Crs_old;
    }

    Cinf_n = 0.5 * (Crw_old + Crw_n);
}


