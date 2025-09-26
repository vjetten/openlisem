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
    const double TOLERANCE = 1e-6;
    const int MAX_ITER = 100;

    if (a <= 0.0) return 0.0;
    if (a >= 1.0) return M_PI * 2;

    double theta = M_PI; // Initial guess
    for (int i = 0; i < MAX_ITER; ++i) {
        double f = (theta - std::sin(theta)) / (2.0 * M_PI) - a;
        double df = (1.0 - std::cos(theta)) / (2.0 * M_PI);
        double delta = f / df;
        theta -= delta;
        if (std::abs(delta) < TOLERANCE)
            break;
    }

    return theta;
}

//---------------------------------------------------------------------------
void TWorld::chanHandPCirc(int r, int c)//, double Area)
{
    // is always a culvert
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
    double a = Area/ChannelMaxArea->Drc;
    double theta = pipeThetafroma(r,c,a);
    ChannelPerimeter->Drc = 0.5*ChannelDiameter->Drc*theta;
    ChannelWH->Drc = 0.5*ChannelDiameter->Drc*(1-cos(theta/2.0));
}
//---------------------------------------------------------------------------
void TWorld::chanHandPTrap(int r, int c)//, double Area)
{
    // if not drowned use abc rule, else assume a rectangular section above trapezium channel
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;

    if (ChannelSide->Drc == 0) {
        chanHandPRect(r,c);
        return;
    }
    // A = h*(wb + m*h)
    // A=wb​h+mh2 -> mh2+wb​h−A=0
    // P=wb​+2*sqrt(h^2+(mh)^2) =wb​+2h*sqrt(1+m2)

    if (Area < ChannelMaxArea->Drc) {
        double B = ChannelWidthB->Drc;
        ChannelWH->Drc = (-B + std::sqrt(B*B + 4*ChannelSide->Drc*Area))/(2.0*ChannelSide->Drc);
    } else {
        // drowned
        double dh = (ChannelWaterVol->Drc-ChannelMaxArea->Drc*ChannelDX->Drc)/(ChannelWidth->Drc*DX->Drc);
        ChannelWH->Drc = ChannelDepth->Drc + dh;//(Area-ChannelMaxArea->Drc)/ChannelWidth->Drc;
    }
    ChannelPerimeter->Drc = ChannelWidthB->Drc+2*ChannelWH->Drc*std::sqrt(1+ChannelSide->Drc*ChannelSide->Drc);
}
//---------------------------------------------------------------------------
void TWorld::chanHandPTria(int r, int c)//, double Area)
{
    if (ChannelSide->Drc == 0) {
        chanHandPRect(r,c);
        return;
    }
    // if drowned assume a rectangular section above triangular channel
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
    if (Area < ChannelMaxArea->Drc) {
        ChannelWH->Drc = std::sqrt(Area/ChannelSide->Drc);
        ChannelPerimeter->Drc = 2*ChannelWH->Drc*std::sqrt(1+ChannelSide->Drc*ChannelSide->Drc);
    } else {
        ChannelWH->Drc = ChannelDepth->Drc + (Area-ChannelMaxArea->Drc)/ChannelWidth->Drc;
        ChannelPerimeter->Drc = 2*ChannelWH->Drc*std::sqrt(1+ChannelSide->Drc*ChannelSide->Drc);
    }
}
//---------------------------------------------------------------------------
void TWorld::chanHandPRect(int r, int c)//, double Area)
{
    double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
    ChannelWH->Drc = Area/ChannelWidth->Drc;
    ChannelPerimeter->Drc = ChannelWidth->Drc+2*ChannelWH->Drc;
}

