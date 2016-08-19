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
 \file lisUnifiedFlowBoundary.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

double TWorld::UF_BoundaryFlux2D(double dt, double cellx, double celly, double f, double s, double fu, double fv, double su, double sv, double slopeX, double slopeY, double NN, int dr, int dc )
{
    double Perim;
    const double beta = 0.6;
    const double _23 = 2.0/3.0;
    double beta1 = 1/beta;
    double h = (f + s)/(cellx*celly);
    double grad = sqrt(slopeX*slopeX + slopeY*slopeY);
    double R;
    double A;
    double Q;
    double V;

    // avg WH from soil surface and roads, over width FlowWidth
    Perim = 2*h+celly;

    if (Perim > 0)
        R = h*celly/Perim;
    else
        R = 0;

    A = pow(NN/sqrt(grad) * pow(Perim, _23),beta);

    if (A > 0)
        Q = pow((celly*h)/A, beta1);
    else
        Q = 0;

    V = pow(R, _23)*sqrt(grad)/NN;

    return std::min((f+s) * UF_Courant,dt * Q);
}

double TWorld::UF_BoundaryFlux1D(double dt, double width, double f, double s, double fu, double su, double slope,double NN, bool front )
{
    double Perim;
    const double beta = 0.6;
    const double _23 = 2.0/3.0;
    double beta1 = 1/beta;
    double h = (f + s)/(width*_dx);
    double grad = std::fabs(slope);
    double R;
    double A;
    double Q;
    double V;

    // avg WH from soil surface and roads, over width FlowWidth
    Perim = 2*h+width;

    if (Perim > 0)
        R = h*width/Perim;
    else
        R = 0;

    A = pow(NN/sqrt(grad) * pow(Perim, _23),beta);

    if (A > 0)
        Q = pow((width*h)/A, beta1);
    else
        Q = 0;

    V = pow(R, _23)*sqrt(grad)/NN;

    return std::min((f+s) * UF_Courant,dt * Q);


}
