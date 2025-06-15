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


#include "model.h"


//---------------------------------------------------------------------------
double TWorld::pipeThetafroma(int r, int c, double a)
{
    if (a < 1) {
        double theta_next;
        double theta = PI;
        double tol = 1e-6;
        // get angle theta from a
        for (int j = 0; j < 50; j++ ) {
           double f = (theta - sin(theta)) / (2 * PI) - a;
           double df = (1 - cos(theta)) / (2 * PI);
           theta_next = theta - f / df;
           if (abs(theta_next - theta) < tol)
               break;
           theta = theta_next;
        }
        return(theta_next);
    } else {
        return(2.*PI);
    }
}
//---------------------------------------------------------------------------
void TWorld::chanHandPCirc(int r, int c)//, double Area)
{
    // is always a culvert
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
    double a = Area/ChannelMaxArea->Drc;
    double theta = pipeThetafroma(r,c,a);
    ChannelPerimeter->Drc = ChannelDiameter->Drc/2.0*theta;
    ChannelWH->Drc = 0.5*ChannelDiameter->Drc*(1-cos(theta/2.0));
}
//---------------------------------------------------------------------------
void TWorld::chanHandPTrap(int r, int c)//, double Area)
{
    // if not drowned use abc rule, else assume a rectangular section above trapezium channel
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
    if (Area < ChannelMaxArea->Drc) {
        ChannelWH->Drc = (-ChannelWidthB->Drc+std::sqrt(ChannelWidthB->Drc*ChannelWidthB->Drc
                            -4*ChannelSide->Drc*Area))/(2*ChannelSide->Drc);
        ChannelPerimeter->Drc = ChannelWidthB->Drc+2*ChannelWH->Drc/ChannelCos->Drc;
                //*std::sqrt(1+ChannelSide->Drc*ChannelSide->Drc);
    } else {
        // drowned
        ChannelWH->Drc = ChannelDepth->Drc + (Area-ChannelMaxArea->Drc)/ChannelWidth->Drc;
        ChannelPerimeter->Drc = ChannelWidthB->Drc+2*ChannelDepth->Drc/ChannelCos->Drc;//*std::sqrt(1+ChannelSide->Drc*ChannelSide->Drc);
    }
}
//---------------------------------------------------------------------------
void TWorld::chanHandPTri(int r, int c)//, double Area)
{
    // if drowned assume a rectangular section above triangular channel
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
    if (Area < ChannelMaxArea->Drc) {
        ChannelWH->Drc = std::sqrt(Area/ChannelSide->Drc);
        ChannelPerimeter->Drc = 2*ChannelWH->Drc/ChannelCos->Drc;//*std::sqrt(1+ChannelSide->Drc*ChannelSide->Drc);
    } else {
        ChannelWH->Drc = ChannelDepth->Drc + (Area-ChannelMaxArea->Drc)/ChannelWidth->Drc;
        ChannelPerimeter->Drc = 2*ChannelDepth->Drc/ChannelCos->Drc;//*std::sqrt(1+ChannelSide->Drc*ChannelSide->Drc);
    }
}
//---------------------------------------------------------------------------
void TWorld::chanHandPRect(int r, int c)//, double Area)
{
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
    ChannelWH->Drc = Area/ChannelWidth->Drc;
    ChannelPerimeter->Drc = ChannelWidth->Drc+2*ChannelWH->Drc;
}

