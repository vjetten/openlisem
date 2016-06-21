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


void TWorld::UF_DEMLDDAnalysis(cTMap * _dem, cTMap * _ldd,cTMap * _lddw,
                               cTMap * _lddh,cTMap * _f1D,cTMap * _s1D,cTMap * _f2D,cTMap * _s2D)
{
    FOR_ROW_COL_UF2D
    {
        UF2D_SlopeX->Drc = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X);
        UF2D_SlopeY->Drc = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y);
    }


}

void TWorld::UF2D_Derivative(cTMap * _dem, cTMap * _in, cTMap * _out , int direction)
{

}

void TWorld::UF2D_Derivative2(cTMap * _dem, cTMap * _in, cTMap * _out , int direction)
{

}

double TWorld::UF2D_Derivative(cTMap * _dem, cTMap * _in, int r, int c, int direction)
{
    if(UF_OUTORMV(_dem,r,c))
    {
        return 0;
    }
    if(direction == UF_DIRECTION_X)
    {
        double dx1 = !UF_OUTORMV(_dem,r,c+1)? _in->data[r][c+1] -_in->Drc :0.0;
        double dx2 = !UF_OUTORMV(_dem,r,c-1)? _in->Drc - _in->data[r][c-1] :0.0;

        return (dx1 + dx2)/(2.0*_dx);
    }
    if(direction == UF_DIRECTION_Y)
    {
        double dy1 = !UF_OUTORMV(_dem,r+1,c)? _in->data[r+1][c] -_in->Drc :0.0;
        double dy2 = !UF_OUTORMV(_dem,r-1,c)? _in->Drc - _in->data[r-1][c] :0.0;

        return (dy1 + dy2)/(2.0*_dx);
    }
    return 0;

}

double TWorld::UF2D_Derivative2(cTMap * _dem, cTMap * _in, int r, int c, int direction)
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

        return (x1 - 2.0 * x + x1)/(_dx*_dx);
    }
    if(direction == UF_DIRECTION_Y)
    {
        double y = _in->Drc;
        double y1 = !UF_OUTORMV(_dem,r+1,c)? _in->data[r+1][c] : y;
        double y2 = !UF_OUTORMV(_dem,r-1,c)? _in->data[r-1][c] : y;

        return (y1 - 2.0 * y + y1)/(_dx*_dx);
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

void TWorld::UF1D_Derivative(cTMap * _ldd, cTMap * _in, cTMap * _out)
{


}

void TWorld::UF1D_Derivative2(cTMap * _ldd, cTMap * _in, cTMap * _out)
{

}

double TWorld::UF1D_Derivative(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c)
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
    if(!lddself == 5)
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

    return (x1 - x2)/(2.0*_dx);

}

double TWorld::UF1D_Derivative2(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c)
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
    if(!lddself == 5)
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

    return (x1 - 2.0 * x + x1)/(_dx*_dx);

}

void TWorld::UF_MUSCLE_1(cTMap * in)
{

}

void TWorld::UF_MUSCLE_2(cTMap * in)
{

}
void TWorld::UF_MUSCLE_operate(cTMap * in_1, cTMap * in_2,int operation)
{

}
void TWorld::UF_MUSCLE_operate(double in_1, cTMap * in_2,int operation)
{

}

/*void TWorld::UF_SWAP(cTMap *(* from), cTMap *(* to))
{
    cTMap ** temp = to;
    *to = *from;
    *from = *temp;
}*/

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

double TWorld::UF_MinMod(double a, double b)
{
    if( b > 0)
    {
        return std::min(a,b);
    }else
    {
        return std::max(a,b);
    }

}

