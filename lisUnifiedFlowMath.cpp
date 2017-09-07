
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


double TWorld::UF_DEMACCES(cTMap * dem, cTMap * h,int r, int c, int dr, int dc)
{
    return 1.0;

    double dem1 = dem->data[r][c];
    double h1 = h->data[r][c];

    if(OUTORMV(r + dr,c+dc))
    {
        return 0;
    }

    double dem2 = dem->data[r+dr][c+dc];

    if(h1 < UF_VERY_SMALL)
    {
        return 0.0;
    }

    return std::max(0.0,std::min(1.0,1.0 - std::min(std::max(0.0,(dem2 - dem1)),h1)/h1));
}

double TWorld::UF_5CellAverage(cTMap * m,int r, int c)
{
    int count = 0;
    double total = 0;
    if(!OUTORMV(r,c))
    {
        total += m->data[r][c];
        count++;
    }
    if(!OUTORMV(r-1,c))
    {
        total += m->data[r-1][c];
        count++;
    }
    if(!OUTORMV(r+1,c))
    {
        total += m->data[r+1][c];
        count++;
    }
    if(!OUTORMV(r,c-1))
    {
        total += m->data[r][c-1];
        count++;
    }
    if(!OUTORMV(r,c+1))
    {
        total += m->data[r][c+1];
        count++;
    }

    return count > 0? total/((float)count) : 0.0;

}

void TWorld::UF_DEMLDDAnalysis(cTMap * _dem, cTMap * _ldd,cTMap * _lddw,
                               cTMap * _lddh,cTMap * _f1D,cTMap * _s1D,cTMap * _f2D,cTMap * _s2D)
{

    FOR_ROW_COL_UF2D
    {
        UF2D_SlopeX->Drc = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X);
        UF2D_SlopeY->Drc = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y);
    }


}
double TWorld::UF2D_Derivative(cTMap * _dem, cTMap * _in, int r, int c, int direction, int calculationside, bool useflowbarriers)
{
    if(UF_OUTORMV(_dem,r,c))
    {
        return 0;
    }
    if(direction == UF_DIRECTION_X)
    {
        if(calculationside == UF_DERIVATIVE_LR)
        {

            if(useflowbarriers)
            {
                double dx1 = !(UF_OUTORMV(_dem,r,c+1))? std::max(_in->data[r][c+1],GetFlowBarrierHeight(r,c,0,1)) -_in->Drc :0.0;
                double dx2 = !(UF_OUTORMV(_dem,r,c-1))? _in->Drc - std::max(_in->data[r][c-1],GetFlowBarrierHeight(r,c,0,-1)) :0.0;
                return (dx1 + dx2)/(2.0*_dx);
            }else
            {
                double dx1 = !(UF_OUTORMV(_dem,r,c+1))? _in->data[r][c+1] -_in->Drc :0.0;
                double dx2 = !(UF_OUTORMV(_dem,r,c-1))? _in->Drc - _in->data[r][c-1] :0.0;
                return (dx1 + dx2)/(2.0*_dx);
            }
        }else if(calculationside == UF_DERIVATIVE_L)
        {
            if(useflowbarriers)
            {
                return (!UF_OUTORMV(_dem,r,c-1))? (_in->Drc - std::max(_in->data[r][c-1],GetFlowBarrierHeight(r,c,0,-1)))/_dx :0.0;
            }else
            {
                return (!UF_OUTORMV(_dem,r,c-1))? (_in->Drc -_in->data[r][c-1])/_dx :0.0;
            }
        }else if(calculationside == UF_DERIVATIVE_R)
        {
            if(useflowbarriers)
            {
                return (!UF_OUTORMV(_dem,r,c+1))? (std::max(_in->data[r][c+1],GetFlowBarrierHeight(r,c,0,1)) -_in->Drc)/_dx :0.0;
            }else
            {
                return (!UF_OUTORMV(_dem,r,c+1))? (_in->data[r][c+1] -_in->Drc)/_dx :0.0;
            }
        }
    }
    if(direction == UF_DIRECTION_Y)
    {
        if(calculationside == UF_DERIVATIVE_LR)
        {
            if(useflowbarriers)
            {
                double dy1 = (!UF_OUTORMV(_dem,r+1,c))? std::max(GetFlowBarrierHeight(r,c,1,0),_in->data[r+1][c]) -_in->Drc :0.0;
                double dy2 = (!UF_OUTORMV(_dem,r-1,c))? _in->Drc - std::max(GetFlowBarrierHeight(r,c,-1,0),_in->data[r-1][c]) :0.0;
                return (dy1 + dy2)/(2.0*_dx);
            }else
            {
                double dy1 = (!UF_OUTORMV(_dem,r+1,c))? _in->data[r+1][c] -_in->Drc :0.0;
                double dy2 = (!UF_OUTORMV(_dem,r-1,c))? _in->Drc - _in->data[r-1][c] :0.0;
                return (dy1 + dy2)/(2.0*_dx);
            }
        }else if(calculationside == UF_DERIVATIVE_L)
        {
            if(useflowbarriers)
            {
                return (!UF_OUTORMV(_dem,r-1,c))? (_in->Drc - std::max(GetFlowBarrierHeight(r,c,-1,0),_in->data[r-1][c]))/_dx :0.0;
            }else
            {
                return (!UF_OUTORMV(_dem,r-1,c))? (_in->Drc - _in->data[r-1][c])/_dx :0.0;
            }
        }else if(calculationside == UF_DERIVATIVE_R)
        {
            if(useflowbarriers)
            {
                return (!UF_OUTORMV(_dem,r+1,c))? (std::max(GetFlowBarrierHeight(r,c,1,0),_in->data[r+1][c] -_in->Drc))/_dx :0.0;
            }else
            {
                return (!UF_OUTORMV(_dem,r+1,c))? (_in->data[r+1][c] -_in->Drc)/_dx :0.0;
            }
        }
    }
    return 0;

}

double TWorld::UF2D_Derivative_scaled(cTMap * _dem, cTMap * _in, int r, int c, int direction,double scale, int calculationside)
{
    if(UF_OUTORMV(_dem,r,c))
    {
        return 0;
    }
    if(direction == UF_DIRECTION_X)
    {
        if(calculationside == UF_DERIVATIVE_LR)
        {

            double dx1 = !(UF_OUTORMV(_dem,r,c+1))? scale *_in->data[r][c+1] -scale *_in->Drc :0.0;
            double dx2 = !(UF_OUTORMV(_dem,r,c-1))? scale *_in->Drc - scale *_in->data[r][c-1] :0.0;
            return (dx1 + dx2)/(2.0*_dx);

        }else if(calculationside == UF_DERIVATIVE_L)
        {
            return (!UF_OUTORMV(_dem,r,c-1))? (scale *_in->Drc -scale *_in->data[r][c-1])/_dx :0.0;

        }else if(calculationside == UF_DERIVATIVE_R)
        {
            return (!UF_OUTORMV(_dem,r,c+1))? (scale *_in->data[r][c+1] -scale *_in->Drc)/_dx :0.0;

        }
    }
    if(direction == UF_DIRECTION_Y)
    {
        if(calculationside == UF_DERIVATIVE_LR)
        {
            double dy1 = (!UF_OUTORMV(_dem,r+1,c))? scale *_in->data[r+1][c] -scale *_in->Drc :0.0;
            double dy2 = (!UF_OUTORMV(_dem,r-1,c))? scale *_in->Drc - scale *_in->data[r-1][c] :0.0;
            return (dy1 + dy2)/(2.0*_dx);

        }else if(calculationside == UF_DERIVATIVE_L)
        {
            return (!UF_OUTORMV(_dem,r-1,c))? (scale *_in->Drc - scale *_in->data[r-1][c])/_dx :0.0;

        }else if(calculationside == UF_DERIVATIVE_R)
        {
            return (!UF_OUTORMV(_dem,r+1,c))? (scale *_in->data[r+1][c] -scale *_in->Drc)/_dx :0.0;
        }
    }
    return 0;

}


double TWorld::UF2D_Derivative2(cTMap * _dem, cTMap * _in, int r, int c, int direction, int calculationside)
{
    if(UF_OUTORMV(_dem,r,c))
    {
        return 0;
    }
    if(direction == UF_DIRECTION_X)
    {
        double x = _in->Drc;
        double x1 = !UF_OUTORMV(_dem,r,c+1)? _in->data[r][c+1] : x;
        double x2 = !UF_OUTORMV(_dem,r,c-1)? _in->data[r][c-1] : x;

        return (x1 - 2.0 * x + x2)/(_dx*_dx);
    }
    if(direction == UF_DIRECTION_Y)
    {
        double y = _in->Drc;
        double y1 = !UF_OUTORMV(_dem,r+1,c)? _in->data[r+1][c] : y;
        double y2 = !UF_OUTORMV(_dem,r-1,c)? _in->data[r-1][c] : y;

        return (y1 - 2.0 * y + y2)/(_dx*_dx);
    }
    if(direction == UF_DIRECTION_XY)
    {
        double xy = _in->Drc;
        double xy1 = !UF_OUTORMV(_dem,r+1,c+1)? _in->data[r+1][c+1] : xy;
        double xy2 = !UF_OUTORMV(_dem,r+1,c-1)? _in->data[r+1][c-1] : xy;
        double xy3 = !UF_OUTORMV(_dem,r-1,c+1)? _in->data[r-1][c+1] : xy;
        double xy4 = !UF_OUTORMV(_dem,r-1,c-1)? _in->data[r-1][c-1] : xy;

        return (xy1 - xy2 - xy3 + xy4)/(4.0*_dx*_dx);
    }
    return 0;

}

double TWorld::UF1D_Derivative(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c, bool minmod, int calculationside)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    if(UF_OUTORMV(_ldd,r,c))
    {
        return 0;
    }
    double x = _in->Drc;
    double x1 = 0;
    double x2 = 0;

    //front cell
    int lddself = (int) _ldd->data[r][c];
    if(!(lddself == 5))
    {
        int r2 = r+dy[lddself];
        int c2 = c+dx[lddself];

        if(!UF_OUTORMV(_ldd,r2,c2)){
            x1 =  _in->data[r2][c2];
        }
    }
    //back cells
    double totalwidth = 0;
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
            totalwidth += _lddw->data[r2][c2];
        }
    }
    if(totalwidth > 0)
    {
        x2 = 0;
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
                if(!UF_OUTORMV(_ldd,r2,c2)){
                    x2 += _in->data[r2][c2] * _lddw->data[r2][c2]/totalwidth;

                }
            }
        }
    }else
    {
        x2 = x;
    }


    double dx1 = (x1 - x)/_dx;
    double dx2 = (x - x2)/_dx;



    if(calculationside == UF_DERIVATIVE_R)
    {
        return dx1;
    }else if(calculationside == UF_DERIVATIVE_L)
    {
        return dx2;
    }else
    {
        return (minmod? (0.5 * UF_MinMod(dx1,dx2)) : ((dx1 + dx2)/(2.0)));
    }

}

double TWorld::UF1D_Derivative_scaled(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c, double scale, bool minmod, int calculationside)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    if(UF_OUTORMV(_ldd,r,c))
    {
        return 0;
    }
    double x = scale*_in->Drc;
    double x1 = 0;
    double x2 = 0;

    //front cell
    int lddself = (int) _ldd->data[r][c];
    if(!(lddself == 5))
    {
        int r2 = r+dy[lddself];
        int c2 = c+dx[lddself];

        if(!UF_OUTORMV(_ldd,r2,c2)){
            x1 =  scale*_in->data[r2][c2];
        }
    }
    //back cells
    double totalwidth = 0;
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
            totalwidth += _lddw->data[r2][c2];
        }
    }
    if(totalwidth > 0)
    {
        x2 = 0;
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
                if(!UF_OUTORMV(_ldd,r2,c2)){
                    x2 += scale*_in->data[r2][c2] * _lddw->data[r2][c2]/totalwidth;

                }
            }
        }
    }else
    {
        x2 = x;
    }


    double dx1 = (x1 - x)/_dx;
    double dx2 = (x - x2)/_dx;



    if(calculationside == UF_DERIVATIVE_R)
    {
        return dx1;
    }else if(calculationside == UF_DERIVATIVE_L)
    {
        return dx2;
    }else
    {
        return (minmod? (0.5 * UF_MinMod(dx1,dx2)) : ((dx1 + dx2)/(2.0)));
    }

}

double TWorld::UF1D_Derivative2(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c, int calculationside)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
    if(UF_OUTORMV(_ldd,r,c))
    {
        return 0;
    }
    double x = _in->Drc;
    double x1 = x;
    double x2 = x;

    //front cell
    int lddself = (int) _ldd->data[r][c];
    if(!(lddself == 5))
    {
        int r2 = r+dy[lddself];
        int c2 = c+dx[lddself];

        if(!UF_OUTORMV(_ldd,r2,c2)){
            x1 =  _in->data[r2][c2];
        }
    }
    //back cells
    double totalwidth = 0;
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
            totalwidth += _lddw->data[r2][c2];
        }
    }
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
            if(!UF_OUTORMV(_ldd,r2,c2)){
                x2 += _in->data[r2][c2] * _lddw->data[r2][c2]/totalwidth;
            }
        }
    }

    return (x1 - 2.0 * x + x2)/(_dx*_dx);

}



bool TWorld::UF_OUTORMV(cTMap * mask, int r, int c)
{
    if(r>=0 && r<_nrRows && c>=0 && c<_nrCols)
    {
        if(!pcr::isMV(mask->data[r][c]))
        {
            return false;
        }
    }
    return true;

}

bool TWorld::UF_LDDOUT(cTMap * _ldd, int r, int c, bool front)
{
    if(!INSIDE(r,c))
    {
        return true;
    }

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
    if(front)
    {
        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            if(lddself > 9 || lddself < 0)
            {
                qDebug() << "error ldd" << lddself;
            }
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){
                return false;
            }
        }
    }
    else
    {
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
            if (!UF_OUTORMV(_ldd,r2,c2))
            {
                if(FLOWS_TO(ldd, r2,c2,r,c))
                {
                    return false;
                }
            }
        }


    }
    return true;



}

bool TWorld::UF_NOTIME(cTMap * mask,cTMap * dt, int r, int c)
{
    if(r>=0 && r<_nrRows && c>=0 && c<_nrCols)
    {
        if(pcr::isMV(mask->data[r][c]))
        {
            return false;
        }else
        {
            if(dt->data[r][c] == 0)
            {
                return true;
            }
        }
    }
    return false;

}

double TWorld::UF_MinMod(double a, double b)
{
    double rec = 0;
    if (a >= 0 && b >= 0)
      rec = std::min(a, b);
    else
      if (a <= 0 && b <= 0)
        rec = std::max(a, b);
    return rec;

}


//to calculate stored mass
void TWorld::UF2D_Stored_mass(cTMap * dt, cTMap* _dem, cTMap *_f,cTMap * _s, cTMap * out_f, cTMap * out_s)
{
    cTMap * hf = UF_t1;
    cTMap * hs = UF_t2;
    cTMap * h = UF_t3;

    FOR_ROW_COL_UF2D
    {
        out_f->Drc = 0;//_f->Drc;
        out_s->Drc = 0;//_s->Drc;
        double area = DX->Drc * _dx;
        UF_t1->Drc = _f->Drc/area;
        UF_t2->Drc =_s->Drc/area;
        UF_t3->Drc = UF_t1->Drc + UF_t2->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {
        double hself = h->Drc + _dem->Drc;
        double hdem = _dem->Drc;
        double hdemw = h->Drc + _dem->Drc;

        for(int i = 0; i < 4; i++)
        {
            int dx[4] = {0, 1, -1, 0};
            int dy[4] = {1, 0, 0, -1};

            int r2 = r + dy[i];
            int c2 = c + dx[i];

            if(UF_OUTORMV(_dem,r2,c2))
            {
                continue;
            }
            if((_dem->data[r2][c2] + h->data[r2][c2]) < hdemw)
            {
                hdemw = _dem->Drc +h->data[r2][c2];
            }
            if(_dem->data[r2][c2] < hdem)
            {
                hdem = _dem->data[r2][c2];
            }
        }

        double hstore = std::min(hself,std::max(0.0,(_dem->Drc -hdem)));
        if(h->Drc > 0)
        {
            double fraction = hstore/h->Drc;
            out_f->Drc = _f->Drc * fraction;
            out_s->Drc = _s->Drc * fraction;
        }
    }}}

}

//to calculate stored mass
void TWorld::UF1D_Stored_mass(cTMap * dt, cTMap * _ldd,cTMap * _lddw, cTMap *_f,cTMap * _s, cTMap * out_f, cTMap * out_s)
{


}

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::upstream(cTMap *_LDD, cTMap *_M, cTMap *out)
 * @brief Returns the sum of all values upstream
 *
 * Returns the sum of all values upstream using
 * the local drainage direction map (LDD)
 *
 * @param _LDD : Local Drainage Direction map
 * @param _M : Material map, can be any substance
 * @param out : Output map, sum of all upstream material
 *
 * @see LDD
 */
void TWorld::upstream(cTMap *_LDD, cTMap *_M, cTMap *out)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_MV
    {
        double tot = 0;
        for (int i=1; i<=9; i++)
        {
            // this is the current cell
            if (i==5)
                continue;

            // look around in 8 directions
            int row = r+dy[i];
            int col = c+dx[i];
            int ldd = 0;

            if (INSIDE(row, col) && !pcr::isMV(_LDD->data[row][col]))
                ldd = (int) _LDD->Drc;
            else
                continue;

            // if no MVs and row,col flows to central cell r,c
            if (  //INSIDE(row, col) &&
                  // !pcr::isMV(_LDD->data[row][col]) &&
                  FLOWS_TO(ldd, row, col, r, c)
                  )
            {
                tot += _M->data[row][col];
            }
        }
        out->Drc = tot;
    }
}

double TWorld::UF2D_MaxFlux(cTMap * _dem,cTMap * _f, int r, int c, int dr, int dc)
{
    double mf = std::max(0.0,(_f->Drc/(_dx*_dx) + _dem->Drc) - (std::max(_f->data[r + dr][c+dc]/(_dx*_dx),GetFlowBarrierHeight(r,c,dr,dc)) + _dem->data[r + dr][c+dc]));
    if(std::isnan(mf))
    {
        return 0;
    }
    return mf;
}

double TWorld::UF1D_Value(cTMap * _ldd,cTMap * _lddw,int r, int c,bool front, cTMap * _in)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    if(UF_OUTORMV(_ldd,r,c))
    {
        return 0;
    }
    double x = _in->Drc;
    double x1 = 0;
    double x2 = 0;

    if(front)
    {

        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){
                x1 =  _in->data[r2][c2];
            }
        }

        return x1;
    }else
    {

        //back cells
        double totalwidth = 0;
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
                totalwidth += _lddw->data[r2][c2];
            }
        }
        if(totalwidth > 0)
        {
            x2 = 0;
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
                    if(!UF_OUTORMV(_ldd,r2,c2)){
                        x2 += _in->data[r2][c2] * _lddw->data[r2][c2]/totalwidth;
                    }
                }
            }
        }else
        {
            x2 = x;
        }

        return x2;
    }


}

void TWorld::UF1D_AddValue(cTMap * _ldd,cTMap * _lddw,int r, int c, bool front, cTMap * _in, double add, bool zeromin)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    if(UF_OUTORMV(_ldd,r,c))
    {
        return;
    }

    if(front)
    {

        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){
                _in->data[r2][c2] += add;
                if(zeromin)
                {
                    _in->data[r2][c2] = std::max(_in->data[r2][c2],0.0);
                }
            }
        }

    }else
    {

        //back cells
        double totalwidth = 0;
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
                totalwidth += _lddw->data[r2][c2];
            }
        }
        if(totalwidth > 0)
        {
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
                    if(!UF_OUTORMV(_ldd,r2,c2)){
                        _in->data[r2][c2] += add* _lddw->data[r2][c2]/totalwidth;;
                        if(zeromin)
                        {
                            _in->data[r2][c2] = std::max(_in->data[r2][c2],0.0);
                        }
                    }
                }
            }
        }
    }




}