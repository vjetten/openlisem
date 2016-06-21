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
 \file lisUnifiedFlow.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

//common functions
double TWorld::UF_DragDistribution(double ffraction, double sfraction, double fvel, double svel)
{


}
double TWorld::UF_DragPower(double ffraction, double sfraction, double fvel, double svel)
{


}

double TWorld::UF_QuasiReynolds(double density, double viscosity, double fraction)
{
    return (UF_GravitySqrt * sqrt(_dx) * _dx * density) /( viscosity * fraction);
}

double TWorld::UF_DragCoefficient(double ffraction, double sfraction, double gamma,double viscosity, double rocksize, double density)
{
    if(sfraction == 0)
    {
        return 0;
    }
    double F = (gamma/180.0) * pow(ffraction/sfraction,3) * UF_Reynolds(density,viscosity,ffraction,sfraction,rocksize);
    double G = pow(ffraction,3.5 - 1.0);
    double P = UF_P(rocksize,ffraction, viscosity, sfraction, density);
    return ffraction * sfraction * (1-gamma)/pow((UF_Aspect * UF_TerminalVelocity(rocksize,ffraction,viscosity,sfraction,density)*( P * F + (1-P) * G)),UF_j);
}

double TWorld::UF_Reynolds(double density, double viscosity, double ffraction,double sfraction, double rocksize)
{
    return density *  rocksize * UF_TerminalVelocity(rocksize,ffraction,viscosity,sfraction,density)/viscosity;
}

double TWorld::UF_VirtualMassCoeff(double ffraction, double sfraction)
{
    return 0.5 * (1.0 + 2.0*sfraction)/ffraction;

}

double TWorld::UF_P(double rocksize, double ffraction, double viscosity,double sfraction, double density)
{
    return 0.5;
}


double TWorld::UF_TerminalVelocity(double rocksize, double ffraction, double viscosity,double sfraction, double density)
{
    return 1.5;
}

