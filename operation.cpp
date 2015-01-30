#include "operation.h"
#include <algorithm>
#include <QList>
#include "CsfMap.h"
#include "error.h"


void copy(
    cTMap& raster,
    cTMap const& other)
{
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
        {
            if (!IS_MV_REAL8(&other.Data[r][c])&& !IS_MV_REAL8(&raster.Data[r][c]))
            {
                raster.Data[r][c] = other.Data[r][c];
            }
            else
                SET_MV_REAL4(&raster.Data[r][c]);
        }
}


int countUnits(
    cTMap const& raster)
{
    QList <long> list;
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                if (!list.contains((long)raster.Data[r][c]))
                    list.append((long)raster.Data[r][c]);
            }
    return(list.count());
}


void fill(
    cTMap& raster,
    double value)
{
    int r, c;

    for (r = 0; r < raster.nrRows(); r++)
        for (c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                raster.Data[r][c] = value;
            }
}


/*
// replaces value inside the areas with the average and retains the original values outside
void CTMap::areaAverage(CTMap *area)
{
    QList <UNIT_LIST> aList;
    QList <double> data;
    int i;

    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
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

    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
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

    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
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


double mapTotal(
    cTMap const& raster)
{
    double total = 0;
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                total = total + raster.Data[r][c];
            }
    return (total);
}


double mapAverage(
    cTMap const& raster)
{
    double total = 0;
    double nrcells = 0;
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                total = total + raster.Data[r][c];
                nrcells+=1;
            }
    return (total/nrcells);
}


double mapMinimum(
    cTMap const& raster)
{
    double total = +1e20;
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                if (total > raster.Data[r][c])
                    total = raster.Data[r][c];
            }
    return (total);
}


double mapMaximum(
    cTMap const& raster)
{
    double total = -1e20;
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                if (total < raster.Data[r][c])
                    total = raster.Data[r][c];
            }
    return (total);
}


double getWindowAverage(
    cTMap const& raster,
    int r,
    int c)
{
  double i = 0;
  double sum = 0, avg = 0;
  if (!IS_MV_REAL8(&raster.Data[r-1][c-1]) && raster.Data[r-1][c-1]> 0) { sum += raster.Data[r-1][c-1]; i+=1.0;}
  if (!IS_MV_REAL8(&raster.Data[r-1][c  ]) && raster.Data[r-1][c  ]> 0) { sum += raster.Data[r-1][c  ]; i+=1.0;}
  if (!IS_MV_REAL8(&raster.Data[r-1][c+1]) && raster.Data[r-1][c+1]> 0) { sum += raster.Data[r-1][c+1]; i+=1.0;}
  if (!IS_MV_REAL8(&raster.Data[r  ][c-1]) && raster.Data[r  ][c-1]> 0) { sum += raster.Data[r  ][c-1]; i+=1.0;}
  if (!IS_MV_REAL8(&raster.Data[r  ][c+1]) && raster.Data[r  ][c+1]> 0) { sum += raster.Data[r  ][c+1]; i+=1.0;}
  if (!IS_MV_REAL8(&raster.Data[r+1][c-1]) && raster.Data[r+1][c-1]> 0) { sum += raster.Data[r+1][c-1]; i+=1.0;}
  if (!IS_MV_REAL8(&raster.Data[r+1][c  ]) && raster.Data[r+1][c  ]> 0) { sum += raster.Data[r+1][c  ]; i+=1.0;}
  if (!IS_MV_REAL8(&raster.Data[r+1][c+1]) && raster.Data[r+1][c+1]> 0) { sum += raster.Data[r+1][c+1]; i+=1.0;}
  avg = (i > 0 ? sum / i : 0);

  return(avg);
}


void cover(
    cTMap& raster,
    cTMap const& value1,
    double value2)
{
    int r, c;

    for (r = 0; r < raster.nrRows(); r++)
        for (c = 0; c < raster.nrCols(); c++)
            if (IS_MV_REAL8(&raster.Data[r][c]) && !IS_MV_REAL8(&value1.Data[r][c]))
            {
                raster.Data[r][c] = value2;
            }
    //     else
    //       SET_MV_REAL4(&raster.Data[r][c]);
}


void calcValue(
    cTMap& raster,
    double value,
    int oper)
{
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                switch (oper)
                {
                case ADD: raster.Data[r][c] += value; break;
                case SUB: raster.Data[r][c] -= value; break;
                case MUL: raster.Data[r][c] *= value; break;
                case DIV: if (value > 0) raster.Data[r][c] /= value;
                    else SET_MV_REAL4(&raster.Data[r][c]); break;
                case POW: raster.Data[r][c] = pow(raster.Data[r][c],value); break;
                case MIN: raster.Data[r][c] = std::min(raster.Data[r][c],value); break;//VJ 110420 new
                case MAX: raster.Data[r][c] = std::max(raster.Data[r][c],value); break;
                }
            }
}


void calcMap(
    cTMap& raster,
    cTMap const& value,
    int oper)
{
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                if (!IS_MV_REAL8(&value.Data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: raster.Data[r][c] += value.Data[r][c]; break;
                    case SUB: raster.Data[r][c] -= value.Data[r][c]; break;
                    case MUL: raster.Data[r][c] *= value.Data[r][c]; break;
                    case DIV: if (value.Data[r][c] > 0) raster.Data[r][c] /= value.Data[r][c];
                        else SET_MV_REAL4(&raster.Data[r][c]); break;
                    case POW: raster.Data[r][c] = powl(raster.Data[r][c],value.Data[r][c]); break;
                    case MIN: raster.Data[r][c] = std::min(value.Data[r][c], raster.Data[r][c]); break; //VJ 110420 new
                    case MAX: raster.Data[r][c] = std::max(value.Data[r][c], raster.Data[r][c]); break;
                    }
                }
                else
                    SET_MV_REAL4(&raster.Data[r][c]);
            }
}


void calc2Maps(
    cTMap& raster,
    cTMap const& value1,
    cTMap const& value2,
    int oper)
{
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                if (!IS_MV_REAL8(&value1.Data[r][c]) && !IS_MV_REAL8(&value2.Data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: raster.Data[r][c] = value1.Data[r][c] + value2.Data[r][c]; break;
                    case SUB: raster.Data[r][c] = value1.Data[r][c] - value2.Data[r][c]; break;
                    case MUL: raster.Data[r][c] = value1.Data[r][c] * value2.Data[r][c]; break;
                    case DIV: if (value2.Data[r][c] > 0) raster.Data[r][c] = value1.Data[r][c] / value2.Data[r][c];
                        else SET_MV_REAL4(&raster.Data[r][c]); break;
                    case POW: raster.Data[r][c] = pow(value1.Data[r][c], value2.Data[r][c]); break;
                    case MIN: raster.Data[r][c] = std::min(value1.Data[r][c], value2.Data[r][c]); break; //VJ 110420 new
                    case MAX: raster.Data[r][c] = std::max(value1.Data[r][c], value2.Data[r][c]); break;
                    }
                }
                else
                    SET_MV_REAL4(&raster.Data[r][c]);
            }
}


void calcMapValue(
    cTMap& raster,
    cTMap const& value1,
    double value2,
    int oper)
{
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                if (!IS_MV_REAL8(&value1.Data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: raster.Data[r][c] = value1.Data[r][c] + value2; break;
                    case SUB: raster.Data[r][c] = value1.Data[r][c] - value2; break;
                    case MUL: raster.Data[r][c] = value1.Data[r][c] * value2; break;
                    case DIV: if (value2 > 0) raster.Data[r][c] = value1.Data[r][c] / value2;
                        else SET_MV_REAL4(&raster.Data[r][c]); break;
                    case POW: raster.Data[r][c] = pow(value1.Data[r][c],value2); break;
                    case MIN: raster.Data[r][c] = std::min(value1.Data[r][c],value2); break;//VJ 110420 new
                    case MAX: raster.Data[r][c] = std::max(value1.Data[r][c],value2); break;
                    }
                }
                else
                    SET_MV_REAL4(&raster.Data[r][c]);
            }
}


void checkMap(
    cTMap const& raster,
    int oper,
    double value,
    QString SS)
{
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!IS_MV_REAL8(&raster.Data[r][c]))
            {
                if (oper == LARGER && raster.Data[r][c] > value)
                {
                    ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger than %4.\n").arg(r).arg(c).arg(raster.MapName()).arg(value) + SS;
                    throw 1;
                }
                else
                    if (oper == SMALLER && raster.Data[r][c] < value)
                    {
                        ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller than %4.\n").arg(r).arg(c).arg(raster.MapName()).arg(value) + SS;
                        throw 1;
                    }
                    else
                        if (oper == LARGEREQUAL && raster.Data[r][c] >= value)
                        {
                            ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger or equal than %4.\n").arg(r).arg(c).arg(raster.MapName()).arg(value) + SS;
                            throw 1;
                        }
                        else
                            if (oper == SMALLEREQUAL && raster.Data[r][c] <= value)
                            {
                                ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller or equal than %4.\n").arg(r).arg(c).arg(raster.MapName()).arg(value) + SS;
                                throw 1;
                            }
            }
}
