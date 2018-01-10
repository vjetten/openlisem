
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

    double hdemdif = 1.0 - std::max(0.0,std::min(1.0,std::min((dem2 - dem1),h1)/h1));

    return hdemdif;
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
    {
      rec = std::min(a, b);
    }else if (a <= 0 && b <= 0)
    {
        rec = std::max(a, b);
    }
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

// poly34.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
// Thanks to Alexandr Rakhmanin <rakhmanin (at) gmail.com>
// public domain
//

#define	TwoPi  6.28318530717958648
const double eps=1e-14;

//=============================================================================
// _root3, root3 from http://prografix.narod.ru
//=============================================================================
static double _root3 ( double x )
{
    double s = 1.;
    while ( x < 1. )
    {
        x *= 8.;
        s *= 0.5;
    }
    while ( x > 8. )
    {
        x *= 0.125;
        s *= 2.;
    }
    double r = 1.5;
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    return r * s;
}

double TWorld::root3 ( double x )
{
    if ( x > 0 ) return _root3 ( x ); else
    if ( x < 0 ) return-_root3 (-x ); else
    return 0.;
}


// x - array of size 2
// return 2: 2 real roots x[0], x[1]
// return 0: pair of complex roots: x[0]�i*x[1]
int   TWorld::SolveP2(double *x, double a, double b) {			// solve equation x^2 + a*x + b = 0
    double D = 0.25*a*a - b;
    if (D >= 0) {
        D = sqrt(D);
        x[0] = 0.5*a + D;
        x[1] = 0.5*a - D;
        return 2;
    }
    x[0] = 0.5*a;
    x[1] = sqrt(-D);
    return 0;
}
//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] � i*x[2], return 1
int TWorld::SolveP3(double *x,double a,double b,double c) {	// solve cubic equation x^3 + a*x^2 + b*x + c = 0
    /*double a2 = a*a;
    double q  = (a2 - 3.0*b)/9.0;
    double r  = (a*(2.0*a2-9.0*b) + 27.0*c)/54.0;
    // equation x^3 + q*x + r = 0
    double r2 = r*r;
    double q3 = q*q*q;
    double A,B;
    if (r2 <= (q3 + eps)) {//<<-- FIXED!
        double t=r/sqrt(q3);
        if( t<-1.0) t=-1.0;
        if( t> 1.0) t= 1.0;
        t=acos(t);
        a/=3; q=-2.0*sqrt(q);
        x[0]=q*cos(t/3)-a;
        x[1]=q*cos((t+TwoPi)/3.0)-a;
        x[2]=q*cos((t-TwoPi)/3.0)-a;
        return(3);
    } else {
        //A =-pow(fabs(r)+sqrt(r2-q3),1./3);
        A =-root3(fabs(r)+sqrt(r2-q3));
        if( r<0 ) A=-A;
        B = A==0? 0 : B=q/A;

        a/=3.0;
        x[0] =(A+B)-a;
        x[1] =-0.5*(A+B)-a;
        x[2] = 0.5*sqrt(3.)*(A-B);
        if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
        return(1);
    }*/

    FLOAT a1 = a;
        FLOAT a2 = b;
        FLOAT a3 = c;

        double_t Q = (a1 * a1 - 3 * a2) / 9;
        double_t R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
        double_t Qcubed = Q * Q * Q;
        double_t d = Qcubed - R * R;

        /* Three real roots */
        if (d >= 0) {
            double_t theta = acos(R / sqrt(Qcubed));
            double_t sqrtQ = sqrt(Q);
            x[0] = -2 * sqrtQ * cos( theta             / 3) - a1 / 3;
            x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
            x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
            return (3);
        }

        /* One real root */
        else {
            double_t e = pow(sqrt(-d) + fabs(R), 1. / 3.);
            if (R > 0)
                e = -e;
            x[0] = (e + Q / e) - a1 / 3.;
            return (1);
        }

}// SolveP3(double *x,double a,double b,double c) {
//---------------------------------------------------------------------------
// a>=0!
void  TWorld::CSqrt( double x, double y, double &a, double &b) // returns:  a+i*s = sqrt(x+i*y)
{
    double r  = sqrt(x*x+y*y);
    if( y==0 ) {
        r = sqrt(r);
        if(x>=0) { a=r; b=0; } else { a=0; b=r; }
    } else {		// y != 0
        a = sqrt(0.5*(x+r));
        b = 0.5*y/a;
    }
}
//---------------------------------------------------------------------------
int   TWorld::SolveP4Bi(double *x, double b, double d)	// solve equation x^4 + b*x^2 + d = 0
{
    double D = b*b-4*d;
    if( D>=0 )
    {
        double sD = sqrt(D);
        double x1 = (-b+sD)/2;
        double x2 = (-b-sD)/2;	// x2 <= x1
        if( x2>=0 )				// 0 <= x2 <= x1, 4 real roots
        {
            double sx1 = sqrt(x1);
            double sx2 = sqrt(x2);
            x[0] = -sx1;
            x[1] =  sx1;
            x[2] = -sx2;
            x[3] =  sx2;
            return 4;
        }
        if( x1 < 0 )				// x2 <= x1 < 0, two pair of imaginary roots
        {
            double sx1 = sqrt(-x1);
            double sx2 = sqrt(-x2);
            x[0] =    0;
            x[1] =  sx1;
            x[2] =    0;
            x[3] =  sx2;
            return 0;
        }
        // now x2 < 0 <= x1 , two real roots and one pair of imginary root
            double sx1 = sqrt( x1);
            double sx2 = sqrt(-x2);
            x[0] = -sx1;
            x[1] =  sx1;
            x[2] =    0;
            x[3] =  sx2;
            return 2;
    } else { // if( D < 0 ), two pair of compex roots
        double sD2 = 0.5*sqrt(-D);
        CSqrt(-0.5*b, sD2, x[0],x[1]);
        CSqrt(-0.5*b,-sD2, x[2],x[3]);
        return 0;
    } // if( D>=0 )
} // SolveP4Bi(double *x, double b, double d)	// solve equation x^4 + b*x^2 d
//---------------------------------------------------------------------------
#define SWAPdb(a,b) { t=b; b=a; a=t; }
static void  dblSort3( double &a, double &b, double &c) // make: a <= b <= c
{
    double t;
    if( a>b ) SWAPdb(a,b);	// now a<=b
    if( c<b ) {
        SWAPdb(b,c);			// now a<=b, b<=c
        if( a>b ) SWAPdb(a,b);// now a<=b
    }
}
//---------------------------------------------------------------------------
int   TWorld::SolveP4De(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
{
    //if( c==0 ) return SolveP4Bi(x,b,d); // After that, c!=0
    if( fabs(c)<1e-14*(fabs(b)+fabs(d)) ) return SolveP4Bi(x,b,d); // After that, c!=0

    int res3 = SolveP3( x, 2*b, b*b-4*d, -c*c);	// solve resolvent
    // by Viet theorem:  x1*x2*x3=-c*c not equals to 0, so x1!=0, x2!=0, x3!=0
    if( res3>1 )	// 3 real roots,
    {
        dblSort3(x[0], x[1], x[2]);	// sort roots to x[0] <= x[1] <= x[2]
        // Note: x[0]*x[1]*x[2]= c*c > 0
        if( x[0] > 0) // all roots are positive
        {
            double sz1 = sqrt(x[0]);
            double sz2 = sqrt(x[1]);
            double sz3 = sqrt(x[2]);
            // Note: sz1*sz2*sz3= -c (and not equal to 0)
            if( c>0 )
            {
                x[0] = (-sz1 -sz2 -sz3)/2;
                x[1] = (-sz1 +sz2 +sz3)/2;
                x[2] = (+sz1 -sz2 +sz3)/2;
                x[3] = (+sz1 +sz2 -sz3)/2;
                return 4;
            }
            // now: c<0
            x[0] = (-sz1 -sz2 +sz3)/2;
            x[1] = (-sz1 +sz2 -sz3)/2;
            x[2] = (+sz1 -sz2 -sz3)/2;
            x[3] = (+sz1 +sz2 +sz3)/2;
            return 4;
        } // if( x[0] > 0) // all roots are positive
        // now x[0] <= x[1] < 0, x[2] > 0
        // two pair of comlex roots
        double sz1 = sqrt(-x[0]);
        double sz2 = sqrt(-x[1]);
        double sz3 = sqrt( x[2]);

        if( c>0 )	// sign = -1
        {
            x[0] = -sz3/2;
            x[1] = ( sz1 -sz2)/2;		// x[0]�i*x[1]
            x[2] =  sz3/2;
            x[3] = (-sz1 -sz2)/2;		// x[2]�i*x[3]
            return 0;
        }
        // now: c<0 , sign = +1
        x[0] =   sz3/2;
        x[1] = (-sz1 +sz2)/2;
        x[2] =  -sz3/2;
        x[3] = ( sz1 +sz2)/2;
        return 0;
    } // if( res3>1 )	// 3 real roots,
    // now resoventa have 1 real and pair of compex roots
    // x[0] - real root, and x[0]>0,
    // x[1]�i*x[2] - complex roots,
    // x[0] must be >=0. But one times x[0]=~ 1e-17, so:
    if (x[0] < 0) x[0] = 0;
    double sz1 = sqrt(x[0]);
    double szr, szi;
    CSqrt(x[1], x[2], szr, szi);  // (szr+i*szi)^2 = x[1]+i*x[2]
    if( c>0 )	// sign = -1
    {
        x[0] = -sz1/2-szr;			// 1st real root
        x[1] = -sz1/2+szr;			// 2nd real root
        x[2] = sz1/2;
        x[3] = szi;
        return 2;
    }
    // now: c<0 , sign = +1
    x[0] = sz1/2-szr;			// 1st real root
    x[1] = sz1/2+szr;			// 2nd real root
    x[2] = -sz1/2;
    x[3] = szi;
    return 2;
} // SolveP4De(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
//-----------------------------------------------------------------------------
double TWorld::N4Step(double x, double a,double b,double c,double d)	// one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
{
    double fxs= ((4*x+3*a)*x+2*b)*x+c;	// f'(x)
    if (fxs == 0) return x;	//return 1e99; <<-- FIXED!
    double fx = (((x+a)*x+b)*x+c)*x+d;	// f(x)
    return x - fx/fxs;
}
//-----------------------------------------------------------------------------
// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
// return 2: 2 real roots x[0], x[1] and complex x[2]�i*x[3],
// return 0: two pair of complex roots: x[0]�i*x[1],  x[2]�i*x[3],
int   TWorld::SolveP4(double *x,double a,double b,double c,double d) {	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d by Dekart-Euler method
    // move to a=0:
    double d1 = d + 0.25*a*( 0.25*b*a - 3./64*a*a*a - c);
    double c1 = c + 0.5*a*(0.25*a*a - b);
    double b1 = b - 0.375*a*a;
    int res = SolveP4De( x, b1, c1, d1);
    if( res==4) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; x[3]-= a/4; }
    else if (res==2) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; }
    else             { x[0]-= a/4; x[2]-= a/4; }
    // one Newton step for each real root:
    if( res>0 )
    {
        x[0] = N4Step(x[0], a,b,c,d);
        x[1] = N4Step(x[1], a,b,c,d);
    }
    if( res>2 )
    {
        x[2] = N4Step(x[2], a,b,c,d);
        x[3] = N4Step(x[3], a,b,c,d);
    }
    return res;
}
//-----------------------------------------------------------------------------
#define F5(t) (((((t+a)*t+b)*t+c)*t+d)*t+e)
//-----------------------------------------------------------------------------
double TWorld::SolveP5_1(double a,double b,double c,double d,double e)	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
    int cnt;
    if( fabs(e)<eps ) return 0;

    double brd =  fabs(a);			// brd - border of real roots
    if( fabs(b)>brd ) brd = fabs(b);
    if( fabs(c)>brd ) brd = fabs(c);
    if( fabs(d)>brd ) brd = fabs(d);
    if( fabs(e)>brd ) brd = fabs(e);
    brd++;							// brd - border of real roots

    double x0, f0;					// less than root
    double x1, f1;					// greater than root
    double x2, f2, f2s;				// next values, f(x2), f'(x2)
    double dx;

    if( e<0 ) { x0 =   0; x1 = brd; f0=e; f1=F5(x1); x2 = 0.01*brd; }	// positive root
    else	  { x0 =-brd; x1 =   0; f0=F5(x0); f1=e; x2 =-0.01*brd; }	// negative root

    if( fabs(f0)<eps ) return x0;
    if( fabs(f1)<eps ) return x1;

    // now x0<x1, f(x0)<0, f(x1)>0
    // Firstly 10 bisections
    for( cnt=0; cnt<10; cnt++)
    {
        x2 = (x0 + x1) / 2;					// next point
        //x2 = x0 - f0*(x1 - x0) / (f1 - f0);		// next point
        f2 = F5(x2);				// f(x2)
        if( fabs(f2)<eps ) return x2;
        if( f2>0 ) { x1=x2; f1=f2; }
        else       { x0=x2; f0=f2; }
    }

    // At each step:
    // x0<x1, f(x0)<0, f(x1)>0.
    // x2 - next value
    // we hope that x0 < x2 < x1, but not necessarily
    do {
        if(cnt++>50) break;
        if( x2<=x0 || x2>= x1 ) x2 = (x0 + x1)/2;	// now  x0 < x2 < x1
        f2 = F5(x2);								// f(x2)
        if( fabs(f2)<eps ) return x2;
        if( f2>0 ) { x1=x2; f1=f2; }
        else       { x0=x2; f0=f2; }
        f2s= (((5*x2+4*a)*x2+3*b)*x2+2*c)*x2+d;		// f'(x2)
        if( fabs(f2s)<eps ) { x2=1e99; continue; }
        dx = f2/f2s;
        x2 -= dx;
    } while(fabs(dx)>eps);
    return x2;
} // SolveP5_1(double a,double b,double c,double d,double e)	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------
int   TWorld::SolveP5(double *x,double a,double b,double c,double d,double e)	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
    double r = x[0] = SolveP5_1(a,b,c,d,e);
    double a1 = a+r, b1=b+r*a1, c1=c+r*b1, d1=d+r*c1;
    return 1+SolveP4(x+1, a1,b1,c1,d1);
} // SolveP5(double *x,double a,double b,double c,double d,double e)	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------
