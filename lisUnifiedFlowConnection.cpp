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

//connection
void TWorld::UF2D1D_Connection(cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                               cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                               cTMap * _fu1D,cTMap * _s1D,
                               cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                               cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                               cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                               cTMap * _su2D,cTMap * _sv2D)
{

    FOR_ROW_COL_UF2D
    {
        if(!UF_OUTORMV(_ldd,r,c))
        {
            double fV = sqrt(_fu2D->Drc * _fu2D->Drc + _fv2D->Drc * _fv2D->Drc);
            double fractionf = std::max(0.0,std::min(1.0, _dt*fV/std::max(0.01*_dx,0.5*(_dx - _lddw->Drc))));
            double sV = sqrt(_su2D->Drc * _su2D->Drc + _sv2D->Drc * _sv2D->Drc);
            double fractions = std::max(0.0,std::min(1.0, _dt*sV/std::max(0.01*_dx,0.5*(_dx - _lddw->Drc))));

            double vf = _f1D->Drc + _f2D->Drc * (1-fractionf);
            double vs = _s1D->Drc + _s2D->Drc * (1-fractionf);

            _fu1D->Drc = vf > 0? (_fu1D->Drc *_f1D->Drc + fV * _f2D->Drc * (1-fractionf))/vf: 0.0;
            _su1D->Drc = vs > 0? (_su1D->Drc *_s1D->Drc + sV * _s2D->Drc * (1-fractions))/vs: 0.0;
            _visc1D->Drc = vf > 0? (_visc1D->Drc *_f1D->Drc + _visc2D->Drc * _f2D->Drc * (1-fractionf))/vf: 0.0;
            _d1D->Drc = vf > 0? (_d1D->Drc *_f1D->Drc + _d2D->Drc * _f2D->Drc * (1-fractionf))/vf: 0.0;
            _ifa1D->Drc = vf > 0? (_ifa1D->Drc *_f1D->Drc + _ifa2D->Drc * _f2D->Drc * (1-fractionf))/vf: 0.0;
            _rocksize1D->Drc = vf > 0? (_rocksize1D->Drc *_f1D->Drc + _rocksize2D->Drc * _f2D->Drc * (1-fractionf))/vf: 0.0;

            _f1D->Drc = vf;
            _s1D->Drc = vs;

            _f2D->Drc = _f2D->Drc * (1-fractionf);
            _s2D->Drc = _s2D->Drc * (1-fractionf);
        }
    }

}

void TWorld::UFDEMLDD_Connection(cTMap * dt,cTMap * RemovedMaterial1D, cTMap * RemovedMaterial2D, cTMap * out_DEM,cTMap * out_LDD)
{


}
