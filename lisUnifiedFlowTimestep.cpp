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


////Timestep functions
double TWorld::UF_InitTimeStep(cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                               cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                               cTMap * _fu1D,cTMap * _s1D,
                               cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                               cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                               cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                               cTMap * _su2D,cTMap * _sv2D,cTMap * out_t2d,cTMap * out_dt2d, cTMap * out_dtstep2d,cTMap * out_t1d,cTMap * out_dt1d, cTMap * out_dtstep1d)
{
    FOR_ROW_COL_UF2D
    {
        out_t2d->Drc = 0;
        out_dt2d->Drc = 0;
        out_dtstep2d->Drc = 0;
    }

    FOR_ROW_COL_UF1D
    {
        out_t1d->Drc = 0;
        out_dt1d->Drc = 0;
        out_dtstep1d->Drc = 0;
    }


}

double TWorld::UF_TimeStep(double t, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                           cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                           cTMap * _fu1D,cTMap * _s1D,
                           cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                           cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                           cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                           cTMap * _su2D,cTMap * _sv2D,cTMap * out_t2d,cTMap * out_dt2d, cTMap * out_dtstep2d,cTMap * out_t1d,cTMap * out_dt1d, cTMap * out_dtstep1d)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
    UF_DTMIN = _dt;

    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            UF2D_CellR->Drc = -1;
            UF2D_CellC->Drc = -1;
            UF1D_CellR->Drc = -1;
            UF1D_CellC->Drc = -1;
        }
    }

    bool first = true;
    FOR_ROW_COL_UF2D
    {
        double a = (_dx*_dx);
        double vol = _f2D->Drc;
        double svol = 0;
        double svu = 0;
        double svv = 0;
        if(UF_SOLIDPHASE)
        {
            vol += _s2D->Drc;
            svol = _s2D->Drc;
            svu = _su2D->Drc;
            svv = _sv2D->Drc;
        }
        if(vol > UF_VERY_SMALL)
        {
            //empirical estimation
            double vest = 10.0 * ((svol/(a) + _f2D->Drc/(a)) * (svol/(a) + _f2D->Drc/(a)));

            //actual max velocity
            double v = std::max(std::max(std::max(std::fabs(_fu2D->Drc),std::fabs(_fv2D->Drc)),std::fabs(svu)),std::fabs(svv));

            double kinterm = std::min(1.0,std::pow((svol/(a) + _f2D->Drc/(a)),UF_KINEMATIC_TIMESTEP_POWER));
            //actual accaleration
            double dvf = kinterm * std::max(std::fabs(UF2D_fay1->Drc),std::max(std::fabs(UF2D_fay2->Drc),std::max(std::fabs(UF2D_fax2->Drc),std::max(std::fabs(UF2D_fax1->Drc),std::max(std::fabs(UF2D_fax->Drc),std::fabs(UF2D_fay->Drc))))));
            double dvs = 0;
            if(UF_SOLIDPHASE)
            {
                dvs = kinterm * std::max(std::fabs(UF2D_say1->Drc),std::max(std::fabs(UF2D_say2->Drc),std::max(std::fabs(UF2D_sax2->Drc),std::max(std::fabs(UF2D_sax1->Drc),std::max(std::fabs(UF2D_sax->Drc),std::fabs(UF2D_say->Drc))))));
            }

            //based on velocity and accaleration provide needed timestep
            out_dt2d->Drc = ((v + dvf + vest)>0? std::min(_dt,std::max(UF2D_MinimumDT,(UF_Courant *_dx/(std::max(v,std::max(dvf,std::max(dvs,vest)))))) ): _dt);
            /*if(first)
            {
                    out_dt2d->Drc = 0;
                    first = false;
            }*/

        }else
        {
            out_dt2d->Drc = _dt;
        }

        UF_DTMIN = std::max(UF2D_MinimumDT,std::min(out_dt2d->Drc,UF_DTMIN));
        //out_t2d->Drc = 0;
        //out_dtstep2d->Drc = 0;
    }


    FOR_ROW_COL_UF1D
    {
        if(!UF_OUTORMV(_ldd,r,c))
        {
            double a = (_dx*_lddw->Drc);
            double vol = _f1D->Drc;
            double svol = 0;
            double svu = 0;
            if(UF_SOLIDPHASE)
            {
                vol += _s1D->Drc;
                svol = _s1D->Drc;
                svu = _su1D->Drc;
            }

            if(vol > UF_VERY_SMALL)
            {
                double vest = ((svol/(a) + _f1D->Drc/(a)) * (svol/(a) + _f1D->Drc/(a)));

                double v = std::max(std::fabs(_fu1D->Drc),std::fabs(svol));

                double kinterm = std::min(1.0,std::pow((svol/(a) + _f1D->Drc/(a)),2.0));
                //actual accaleration
                double dvf = kinterm *std::max(std::fabs(UF1D_fa1->Drc),std::max(std::fabs(UF1D_fa2->Drc),std::fabs(UF1D_fa->Drc)));
                double dvs = 0;
                if(UF_SOLIDPHASE)
                {
                    dvs =kinterm *std::max(std::fabs(UF1D_sa1->Drc),std::max(std::fabs(UF1D_sa2->Drc),std::fabs(UF1D_sa->Drc)));
                }

                if(_ldd->Drc == 5)
                {
                    out_dt1d->Drc = (0.25 * UF_Courant *_dx/(std::max(v,std::max(dvf,std::max(dvs,vest)))));
                }else
                {
                    out_dt1d->Drc = 0.25 *(UF_Courant *_dx/(std::max(v,std::max(dvf,std::max(dvs,vest)))));
                }

            }else
            {
                out_dt1d->Drc = _dt;
            }

            UF_DTMIN = std::max(UF1D_MinimumDT,std::min(out_dt1d->Drc,UF_DTMIN));


        }
    }

    FOR_ROW_COL_UF1D
    {
        if(!UF_OUTORMV(_ldd,r,c))
        {
            if(_ldd->Drc == 5)
            {
                out_dt1d->Drc = UF_DTMIN;
            }else
            {
                out_dt1d->Drc = std::min(UF_DTMIN* UF1D_OutletDistance->Drc,out_dt1d->Drc);
            }
        }
    }


    if(!SpatiallyDynamicTimestep)
    {
        FOR_ROW_COL_UF2D
        {
            out_dt2d->Drc = UF_DTMIN;
        }

        FOR_ROW_COL_UF1D
        {
            if(!UF_OUTORMV(_ldd,r,c))
            {
                out_dt1d->Drc = UF_DTMIN;
            }
        }
    }


    FOR_ROW_COL_UF2D
    {
        UF_t1->Drc = std::max(1.0,std::floor(std::min(out_dt2d->Drc,_dt)/UF_DTMIN));
    }

    FOR_ROW_COL_UF1D
    {
        if(!UF_OUTORMV(_ldd,r,c))
        {
            UF_t2->Drc = std::max(1.0,std::floor(std::min(out_dt1d->Drc,_dt)/UF_DTMIN));
        }
    }

    FOR_ROW_COL_UF2D
    {
        double dt = UF_t1->Drc;
        double dt1 = !UF_OUTORMV(_dem,r+1,c)? UF_t1->data[r+1][c]:dt;
        double dt2 = !UF_OUTORMV(_dem,r-1,c)? UF_t1->data[r-1][c]:dt;
        double dt3 = !UF_OUTORMV(_dem,r,c+1)? UF_t1->data[r][c+1]:dt;
        double dt4 = !UF_OUTORMV(_dem,r,c-1)? UF_t1->data[r][c-1]:dt;

        double dt5 = !UF_OUTORMV(_dem,r+2,c)? UF_t1->data[r+2][c]:dt;
        double dt6 = !UF_OUTORMV(_dem,r-2,c)? UF_t1->data[r-2][c]:dt;
        double dt7 = !UF_OUTORMV(_dem,r,c+2)? UF_t1->data[r][c+2]:dt;
        double dt8 = !UF_OUTORMV(_dem,r,c-2)? UF_t1->data[r][c-2]:dt;
        double dt9 = !UF_OUTORMV(_dem,r+1,c+1)? UF_t1->data[r+1][c+1]:dt;
        double dt10 = !UF_OUTORMV(_dem,r-1,c+1)? UF_t1->data[r-1][c+1]:dt;
        double dt11 = !UF_OUTORMV(_dem,r+1,c-1)? UF_t1->data[r+1][c-1]:dt;
        double dt12 = !UF_OUTORMV(_dem,r-1,c-1)? UF_t1->data[r-1][c-1]:dt;

        out_dtstep2d->Drc = dt;

        out_dtstep2d->Drc = std::min(dt,std::min(dt1,std::min(dt2,std::min(dt3,dt4))));
        out_dtstep2d->Drc = std::min(out_dtstep2d->Drc,std::min(2.0 * dt9,std::min(2.0 * dt10,std::min(2.0 *dt11,2.0*dt12))));
        out_dtstep2d->Drc = std::min(out_dtstep2d->Drc,std::min(4.0 * dt5,std::min(4.0 * dt6,std::min(4.0 *dt7,2.0*dt8))));
    }

    FOR_ROW_COL_UF1D
    {
        double dt = UF_t1->Drc;
        double dt1 = !UF_OUTORMV(_dem,r+1,c)? UF_t1->data[r+1][c]:dt;
        double dt2 = !UF_OUTORMV(_dem,r-1,c)? UF_t1->data[r-1][c]:dt;
        double dt3 = !UF_OUTORMV(_dem,r,c+1)? UF_t1->data[r][c+1]:dt;
        double dt4 = !UF_OUTORMV(_dem,r,c-1)? UF_t1->data[r][c-1]:dt;

        if(!UF_OUTORMV(_ldd,r,c))
        {
            double dt = UF_t2->Drc;
            for (int i=1;i<=9;i++)
            {
                int r2, c2, ldd = 0;
                if (i==5)  // Skip current cell
                    continue;
                r2 = r+dy[i];
                c2 = c+dx[i];
                if (!UF_OUTORMV(_ldd,r2,c2))
                    ldd = (int) _ldd->data[r2][c2];
                else
                    continue;
                if (!UF_OUTORMV(_ldd,r2,c2) &&
                        FLOWS_TO(ldd, r2,c2,r,c))
                {
                    dt = std::min(dt,UF_t2->data[r2][c2]);
                }
            }
            out_dtstep1d->Drc = std::min(dt,std::min(dt1,std::min(dt2,std::min(dt3,dt4))));
        }
    }
    int rc = 0;
    int cc = 0;

    int count = 0;
    FOR_ROW_COL_UF2D
    {
        if(!(t + UF_DTMIN < _dt ))
        {
            out_dt2d->Drc = _dt - out_t2d->Drc;


        }else if(false)//long timesteps first
        {
            if(!(out_t2d->Drc > t ))
            {
                out_dt2d->Drc = out_dtstep2d->Drc* UF_DTMIN;
                out_dt2d->Drc = std::max(0.0,std::min(out_dt2d->Drc,_dt - out_t2d->Drc));
            }else
            {
                out_dt2d->Drc = 0;
            }

        }else if(true)//short timesteps first
        {
           //this creates weird atrifacts
           //it does handle movements above courant criteria better though!, but just set courant factor lower
           if(!(out_t2d->Drc + out_dtstep2d->Drc* UF_DTMIN > t))
           {
               out_dt2d->Drc = out_dtstep2d->Drc* UF_DTMIN + out_dtstep2d->Drc* UF_DTMIN;
               out_dt2d->Drc = std::max(0.0,std::min(out_dt2d->Drc,_dt - out_t2d->Drc));
           }else
            {
                out_dt2d->Drc = 0;
            }
        }

        out_t2d->Drc += out_dt2d->Drc;

        if(out_dt2d->Drc > UF_VERY_SMALL)
        {
            UF2D_CellR->data[rc][cc] = r;
            UF2D_CellC->data[rc][cc] = c;

            count ++;

            cc ++;
            if(cc == _nrCols)
            {
                rc ++;
                cc = 0;
            }

            if(INSIDE(rc,cc))
            {
                UF2D_CellR->data[rc][cc] = -1;
                UF2D_CellC->data[rc][cc] = -1;
            }
        }
    }

    rc = 0;
    cc = 0;
    FOR_ROW_COL_UF1D
    {
        if(!UF_OUTORMV(_ldd,r,c))
        {

            if(!(t + UF_DTMIN < _dt ))
            {
                out_dt1d->Drc = _dt - out_t1d->Drc;
            }else if(!(out_t1d->Drc > t))
            {
                out_dt1d->Drc = out_dtstep1d->Drc* UF_DTMIN;
                out_dt1d->Drc = std::min(out_dt1d->Drc,_dt - out_t1d->Drc);
            }else
            {
                out_dt1d->Drc = 0;
            }

            out_t1d->Drc += out_dt1d->Drc;
            if(out_dt1d->Drc > UF_VERY_SMALL)
            {

                UF1D_CellR->data[rc][cc] = r;
                UF1D_CellC->data[rc][cc] = c;

                cc ++;
                if(cc == _nrCols)
                {
                    rc ++;
                    cc = 0;
                }
                if(INSIDE(rc,cc))
                {
                    UF1D_CellR->data[rc][cc] = -1;
                    UF1D_CellC->data[rc][cc] = -1;
                }

            }


        }
    }
    return std::min(_dt- t,UF_DTMIN);
}
