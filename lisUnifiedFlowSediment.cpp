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

void TWorld::UF2D_FluidSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_f,cTMap *out_visc)
{
    FOR_ROW_COL_UF2D
    {
        out_f->Drc = _f->Drc;
        out_visc->Drc = _f->Drc;
    }

}

void TWorld::UF2D_SolidSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_s,cTMap * out_d,cTMap *  out_ifa,cTMap *  out_rocksize)
{
    FOR_ROW_COL_UF2D
    {
        out_s->Drc = _s->Drc;
        out_d->Drc = _d->Drc;
        out_ifa->Drc = _ifa->Drc;
        out_rocksize->Drc = _rocksize->Drc;
    }
}


void TWorld::UF1D_FluidSource(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_f,cTMap * out_visc)
{
    FOR_ROW_COL_UF1D
    {
        out_f->Drc = _f->Drc;
        out_visc->Drc = _f->Drc;
    }
}
void TWorld::UF1D_SolidSource(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_s,cTMap * out_d,cTMap *  out_ifa,cTMap *  out_rocksize)
{
    FOR_ROW_COL_UF1D
    {
        out_s->Drc = _s->Drc;
        out_d->Drc = _d->Drc;
        out_ifa->Drc = _ifa->Drc;
        out_rocksize->Drc = _rocksize->Drc;
    }
}



//transport capacity
double TWorld::UnifiedFlowTransportcapacity(double _surface, double _f, double _visc, double _s, double _d,double _v)
{


}


//active entrainment
double TWorld::UnifiedFlowActiveEntrainment(double _surface, double _f, double _visc, double _s, double _d,double _v)
{


}

void TWorld::UF_FlowDetachment(double dt)
{




}

void TWorld::UF_FlowEntrainment(double dt)
{





}

double TWorld::UF_SoilTake(int r, int c, int d, double potential)
{



}

void TWorld::UF_SoilAdd(int r, int c, int d, double mass)
{


}
