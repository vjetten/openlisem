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
    double beta1 = 1.0/beta;
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

    return std::min((f+s) * UF_Courant,0.1*dt * Q);
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

void TWorld::UF_ForcedConditions(int thread, cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                       cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                       cTMap * _fu1D,cTMap * _s1D,
                       cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                       cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                       cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                       cTMap * _su2D,cTMap * _sv2D)
{
    if(SwitchUFForced)
    {
        FOR_ROW_COL_UF2DMT
        {
            if(!(UF2D_ForcedFVolume->Drc < 0))
            {
                double vdif = (UF2D_ForcedFVolume->Drc - _f2D->Drc);
                UF_InitializedF += std::fabs(std::max(0.0,vdif));
                UF2D_foutflow += std::fabs(std::min(0.0,vdif));
                _f2D->Drc = UF2D_ForcedFVolume->Drc;

            }
            if(UF_SOLIDPHASE)
            {
                if(!(UF2D_ForcedSVolume->Drc < 0))
                {
                    double vdif = (UF2D_ForcedSVolume->Drc - _s2D->Drc);
                    UF_InitializedS += std::fabs(std::max(0.0,vdif));
                    UF2D_soutflow += std::fabs(std::min(0.0,vdif));
                    _s2D->Drc = UF2D_ForcedSVolume->Drc;

                    if(_s2D->Drc > UF_VERY_SMALL)
                    {
                        _d2D->Drc = UF2D_ForcedSDensity->Drc;
                        _ifa2D->Drc = UF2D_ForcedSIFA->Drc;
                        _rocksize2D->Drc = UF2D_ForcedSRocksize->Drc;
                    }
                }
            }
        }}}
    }
}



void TWorld::UF_Initial( cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                       cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                       cTMap * _fu1D,cTMap * _s1D,
                       cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                       cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                       cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                       cTMap * _su2D,cTMap * _sv2D)
{
    //add initial volume of fluid and solid to the dynamic field
    //this also sets the respective properties
    if(!SwitchUFInitial)
    {
        return;
    }

    FOR_ROW_COL_UF2D
    {
        if(UF2D_Initialized->Drc == 0 && !(UF2D_InitialTime->Drc < 0))
        {

            if(!(UF2D_InitialTime->Drc > (time/ 60.0)))
            {
                UF2D_Initialized->Drc = 1;

                if(!(UF2D_InitialFVolume->Drc < 0))
                {
                    _f2D->Drc = UF2D_InitialFVolume->Drc;
                    UF_InitializedF += UF2D_InitialFVolume->Drc;
                }
                if(UF_SOLIDPHASE)
                {
                    if(!(UF2D_InitialSVolume->Drc < 0))
                    {
                        _s2D->Drc = UF2D_InitialSVolume->Drc;
                        UF_InitializedS += UF2D_InitialSVolume->Drc;
                        if(_s2D->Drc > UF_VERY_SMALL)
                        {
                            _d2D->Drc = UF2D_InitialSDensity->Drc;
                            _ifa2D->Drc = UF2D_InitialSIFA->Drc;
                            _rocksize2D->Drc = UF2D_InitialSRocksize->Drc;
                        }
                    }
                }

            }


        }
    }
}
