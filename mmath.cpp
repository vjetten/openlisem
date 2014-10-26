
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
  \file mmath.cpp
  \brief basic map mathematics on two maps or maps and variables

functions: \n
- void CTMap::fill(double value)   \n
- void CTMap::calcV(double v, int oper)   \n
- void CTMap::calc(cTMap *m, int oper)   \n
- void CTMap::calc2(cTMap *m1, cTMap *m2, int oper)   \n
- void CTMap::calc2V(cTMap *m1, double V, int oper)   \n
- void CTMap::copy(cTMap *m)   \n
- void CTMap::cover(double v)   \n
- void CTMap::setMV()   \n
- double CTMap::mapTotal()   \n
operations for 'oper' are ADD, SUB, MUL, DIV, POW, MIN, MAX
*/


#include "model.h"

/*
//__fastcall ---------------------------------------------------------------------------
CTMap::CTMap()
    : cTMap()
{

}
//__fastcall ---------------------------------------------------------------------------
CTMap::~CTMap()
{

}
*/
//---------------------------------------------------------------------------
void CTMap::fill(double value)
{
    int r, c;

    for (r = 0; r < nrRows; r++)
        for (c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                Data[r][c] = value;
            }
}
//---------------------------------------------------------------------------
void CTMap::cover(cTMap *M, double value)
{
    int r, c;

    for (r = 0; r < nrRows; r++)
        for (c = 0; c < nrCols; c++)
            if (IS_MV_REAL8(&Data[r][c]) && !IS_MV_REAL8(&M->Data[r][c]))
            {
                Data[r][c] = value;
            }
    //     else
    //       SET_MV_REAL4(&Data[r][c]);
}
//---------------------------------------------------------------------------
void CTMap::calcMap(cTMap *M, int oper)
{
    for (int r = 0; r < M->nrRows; r++)
        for (int c = 0; c < M->nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if (!IS_MV_REAL8(&M->Data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: Data[r][c] += M->Data[r][c]; break;
                    case SUB: Data[r][c] -= M->Data[r][c]; break;
                    case MUL: Data[r][c] *= M->Data[r][c]; break;
                    case DIV: if (M->Data[r][c] > 0) Data[r][c] /= M->Data[r][c];
                        else SET_MV_REAL4(&Data[r][c]); break;
                    case POW: Data[r][c] = powl(Data[r][c],M->Data[r][c]); break;
                    case MIN: Data[r][c] = _min(M->Data[r][c], Data[r][c]); break; //VJ 110420 new
                    case MAX: Data[r][c] = _max(M->Data[r][c], Data[r][c]); break;
                    }
                }
                else
                    SET_MV_REAL4(&Data[r][c]);
            }
}
//---------------------------------------------------------------------------
void CTMap::calc2Maps(cTMap *M1, cTMap *M2, int oper)
{
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if (!IS_MV_REAL8(&M1->Data[r][c]) && !IS_MV_REAL8(&M2->Data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: Data[r][c] = M1->Data[r][c] + M2->Data[r][c]; break;
                    case SUB: Data[r][c] = M1->Data[r][c] - M2->Data[r][c]; break;
                    case MUL: Data[r][c] = M1->Data[r][c] * M2->Data[r][c]; break;
                    case DIV: if (M2->Data[r][c] > 0) Data[r][c] = M1->Data[r][c] / M2->Data[r][c];
                        else SET_MV_REAL4(&Data[r][c]); break;
                    case POW: Data[r][c] = pow(M1->Data[r][c], M2->Data[r][c]); break;
                    case MIN: Data[r][c] = _min(M1->Data[r][c], M2->Data[r][c]); break; //VJ 110420 new
                    case MAX: Data[r][c] = _max(M1->Data[r][c], M2->Data[r][c]); break;
                    }
                }
                else
                    SET_MV_REAL4(&Data[r][c]);
            }
}
//---------------------------------------------------------------------------
void CTMap::calcValue(double V, int oper)
{
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                switch (oper)
                {
                case ADD: Data[r][c] += V; break;
                case SUB: Data[r][c] -= V; break;
                case MUL: Data[r][c] *= V; break;
                case DIV: if (V > 0) Data[r][c] /= V;
                    else SET_MV_REAL4(&Data[r][c]); break;
                case POW: Data[r][c] = pow(Data[r][c],V); break;
                case MIN: Data[r][c] = _min(Data[r][c],V); break;//VJ 110420 new
                case MAX: Data[r][c] = _max(Data[r][c],V); break;
                }
            }
}
//---------------------------------------------------------------------------
double CTMap::mapTotal()
{
    double total = 0;
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                total = total + Data[r][c];
            }
    return (total);

}
//---------------------------------------------------------------------------
double CTMap::mapAverage()
{
    double total = 0;
    double nrcells = 0;
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                total = total + Data[r][c];
                nrcells+=1;
            }
    return (total/nrcells);
}
//---------------------------------------------------------------------------
/*
// replaces value inside the areas with the average and retains the original values outside
void CTMap::areaAverage(CTMap *area)
{
    QList <UNIT_LIST> aList;
    QList <double> data;
    int i;

    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if(!data.contains(area->Data[r][c]))
                    data.append(area->Data[r][c]);
            }
    qSort(data);

    for (i = 0; i < data.count(); i++)
    {
        UNIT_LIST ul;
        ul.nr = 0;
        ul.var0 = data[i];
        ul.var1 = 0;
        ul.var2 = 0;
        ul.var3 = 0;
        aList.append(ul);
    }
    // count, sort and initialize diff units

    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&area->Data[r][c]))
            {
                for (i = 0; i < aList.count(); i++)
                {
                    if(area->Data[r][c] == aList[i].var0)
                    {
                        aList[i].var1 += 1.0;
                        aList[i].var2 += Data[r][c];
                    }
                }
            }
    // count the sums within each area

    for (i = 0; i < aList.count(); i++)
        if(aList[i].var1 > 0)
            aList[i].var3 = aList[i].var2/aList[i].var1;
        else
            aList[i].var3 = 0;
    // calculate the average values for each area

    // for (i = 0; i < aList.count(); i++)
    // qDebug() << aList[i].area <<  aList[i].totsl; //aList[i].totdet << aList[i].totdep  <<

    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&area->Data[r][c]))
            {
                for (i = 0; i < aList.count(); i++)
                {
                    if(area->Data[r][c] == aList[i].var0)
                        Data[r][c] = aList[i].var3;
                }
            }

    return;

}
*/
//---------------------------------------------------------------------------
double CTMap::mapMaximum()
{
    double total = -1e20;
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if (total < Data[r][c])
                    total = Data[r][c];
            }
    return (total);

}
//---------------------------------------------------------------------------
double CTMap::mapMinimum()
{
    double total = +1e20;
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if (total > Data[r][c])
                    total = Data[r][c];
            }
    return (total);

}
//---------------------------------------------------------------------------
void CTMap::copy(cTMap *M)
{
    for (int r = 0; r < M->nrRows; r++)
        for (int c = 0; c < M->nrCols; c++)
        {
            if (!IS_MV_REAL8(&M->Data[r][c])&& !IS_MV_REAL8(&Data[r][c]))
            {
                Data[r][c] = M->Data[r][c];
            }
            else
                SET_MV_REAL4(&Data[r][c]);
        }
}
//---------------------------------------------------------------------------
/*
void CTMap::shift(cTMap *M, int dr, int dc)
{
    for (int r = 0; r < M->nrRows; r++)
        for (int c = 0; c < M->nrCols; c++)
        {
            if (!IS_MV_REAL8(&M->Data[r][c])&& !IS_MV_REAL8(&Data[r][c]))
            {
                Data[r][c] = M->Data[r][c];
                if (!IS_MV_REAL8(&Data[r+dr][c+dc]))
                    Data[r][c] = M->Data[r+dr][c+dc];
            }
            else
                SET_MV_REAL4(&Data[r][c]);
        }
}
*/
//---------------------------------------------------------------------------
void CTMap::setMV()
{
    for(int r = 0; r < nrRows; r++)
        SetMemMV(Data[r],nrCols,CR_REAL8);
}
//---------------------------------------------------------------------------
void CTMap::calcMapValue(cTMap *M1, double V, int oper)
{
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if (!IS_MV_REAL8(&M1->Data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: Data[r][c] = M1->Data[r][c] + V; break;
                    case SUB: Data[r][c] = M1->Data[r][c] - V; break;
                    case MUL: Data[r][c] = M1->Data[r][c] * V; break;
                    case DIV: if (V > 0) Data[r][c] = M1->Data[r][c] / V;
                        else SET_MV_REAL4(&Data[r][c]); break;
                    case POW: Data[r][c] = pow(M1->Data[r][c],V); break;
                    case MIN: Data[r][c] = _min(M1->Data[r][c],V); break;//VJ 110420 new
                    case MAX: Data[r][c] = _max(M1->Data[r][c],V); break;
                    }
                }
                else
                    SET_MV_REAL4(&Data[r][c]);
            }
}
//---------------------------------------------------------------------------
void CTMap::checkMap(int oper, double V, QString SS)
{
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if (oper == LARGER && Data[r][c] > V)
                {
                    ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger than %4.\n").arg(r).arg(c).arg(MapName).arg(V) + SS;
                    throw 1;
                }
                else
                    if (oper == SMALLER && Data[r][c] < V)
                    {
                        ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller than %4.\n").arg(r).arg(c).arg(MapName).arg(V) + SS;
                        throw 1;
                    }
                    else
                        if (oper == LARGEREQUAL && Data[r][c] >= V)
                        {
                            ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger or equal than %4.\n").arg(r).arg(c).arg(MapName).arg(V) + SS;
                            throw 1;
                        }
                        else
                            if (oper == SMALLEREQUAL && Data[r][c] <= V)
                            {
                                ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller or equal than %4.\n").arg(r).arg(c).arg(MapName).arg(V) + SS;
                                throw 1;
                            }
            }
}
//---------------------------------------------------------------------------
int CTMap::countUnits()
{
    QList <long> list;
    for (int r = 0; r < nrRows; r++)
        for (int c = 0; c < nrCols; c++)
            if (!IS_MV_REAL8(&Data[r][c]))
            {
                if (!list.contains((long)Data[r][c]))
                    list.append((long)Data[r][c]);
            }
    return(list.count());
}
//---------------------------------------------------------------------------
// gives back the average of surrounding cells
double CTMap::getWindowAverage(int r, int c)
{
  double i = 0;
  double sum = 0, avg = 0;
  if (!IS_MV_REAL8(&Data[r-1][c-1]) && Data[r-1][c-1]> 0) { sum += Data[r-1][c-1]; i+=1.0;}
  if (!IS_MV_REAL8(&Data[r-1][c  ]) && Data[r-1][c  ]> 0) { sum += Data[r-1][c  ]; i+=1.0;}
  if (!IS_MV_REAL8(&Data[r-1][c+1]) && Data[r-1][c+1]> 0) { sum += Data[r-1][c+1]; i+=1.0;}
  if (!IS_MV_REAL8(&Data[r  ][c-1]) && Data[r  ][c-1]> 0) { sum += Data[r  ][c-1]; i+=1.0;}
  if (!IS_MV_REAL8(&Data[r  ][c+1]) && Data[r  ][c+1]> 0) { sum += Data[r  ][c+1]; i+=1.0;}
  if (!IS_MV_REAL8(&Data[r+1][c-1]) && Data[r+1][c-1]> 0) { sum += Data[r+1][c-1]; i+=1.0;}
  if (!IS_MV_REAL8(&Data[r+1][c  ]) && Data[r+1][c  ]> 0) { sum += Data[r+1][c  ]; i+=1.0;}
  if (!IS_MV_REAL8(&Data[r+1][c+1]) && Data[r+1][c+1]> 0) { sum += Data[r+1][c+1]; i+=1.0;}
  avg = (i > 0 ? sum / i : 0);

  return(avg);
}
//---------------------------------------------------------------------------
