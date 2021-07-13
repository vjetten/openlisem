#include "operation.h"
#include <algorithm>
#include <QList>
#include <qmath.h>
#include "CsfMap.h"
#include "error.h"


void copy(
    cTMap& raster,
    cTMap const& other)
{

    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
        {
            if (!pcr::isMV(other.data[r][c])&& !pcr::isMV(raster.data[r][c]))
            {
                raster.data[r][c] = other.data[r][c];
            }
            else
                pcr::setMV(raster.data[r][c]);
        }
}


int countUnits(cTMap const& raster)
{
    QList <long> list;
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!pcr::isMV(raster.data[r][c]))
            {
                if (!list.contains((long)raster.data[r][c]))
                    list.append((long)raster.data[r][c]);
            }

    return(list.count());
}


void fill(
    cTMap& raster,
    double value)
{
    int r, c;
//#pragma omp parallel for collapse(2)
    for (r = 0; r < raster.nrRows(); r++)
        for (c = 0; c < raster.nrCols(); c++)
            if (!pcr::isMV(raster.data[r][c]))
            {
                raster.data[r][c] = value;
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
            if (!pcr::isMV(data[r][c]))
            {
                if(!data.contains(area->data[r][c]))
                    data.append(area->data[r][c]);
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
            if (!pcr::isMV(area->data[r][c]))
            {
                for (i = 0; i < aList.count(); i++)
                {
                    if(area->data[r][c] == aList[i].var0)
                    {
                        aList[i].var1 += 1.0;
                        aList[i].var2 += data[r][c];
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
            if (!pcr::isMV(area->data[r][c]))
            {
                for (i = 0; i < aList.count(); i++)
                {
                    if(area->data[r][c] == aList[i].var0)
                        data[r][c] = aList[i].var3;
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
            if (!pcr::isMV(raster.data[r][c]))
            {
                total = total + raster.data[r][c];
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
            if (!pcr::isMV(raster.data[r][c]))
            {
                total = total + raster.data[r][c];
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
            if (!pcr::isMV(raster.data[r][c]))
            {
                if (total > raster.data[r][c])
                    total = raster.data[r][c];
            }
    return (total);
}


double mapMaximum(
    cTMap const& raster)
{
    double total = -1e20;
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!pcr::isMV(raster.data[r][c]))
            {
                if (total < raster.data[r][c])
                    total = raster.data[r][c];
            }
    return (total);
}


double getWindowAverage(
    cTMap const& raster,
    int r,
    int c,
    bool center)
{
  double i = 0;
  double sum = 0, avg = 0;
  if (r > 0 && c > 0 && !pcr::isMV(raster.data[r-1][c-1]) && raster.data[r-1][c-1]> 0) { sum += raster.data[r-1][c-1]; i+=1.0;}
  if (r > 0 && !pcr::isMV(raster.data[r-1][c  ]) && raster.data[r-1][c  ]> 0) { sum += raster.data[r-1][c  ]; i+=1.0;}
  if (r > 0 && c < raster.nrCols()-1 && !pcr::isMV(raster.data[r-1][c+1]) && raster.data[r-1][c+1]> 0) { sum += raster.data[r-1][c+1]; i+=1.0;}
  if (c < raster.nrCols()-1 && !pcr::isMV(raster.data[r  ][c-1]) && raster.data[r  ][c-1]> 0) { sum += raster.data[r  ][c-1]; i+=1.0;}
  if (center && !pcr::isMV(raster.data[r][c]) && raster.data[r][c]> 0) { sum += raster.data[r  ][c]; i+=1.0;}
  if (c < raster.nrCols()-1 && !pcr::isMV(raster.data[r  ][c+1]) && raster.data[r  ][c+1]> 0) { sum += raster.data[r  ][c+1]; i+=1.0;}
  if (r < raster.nrRows()-1 && c > 0 && !pcr::isMV(raster.data[r+1][c-1]) && raster.data[r+1][c-1]> 0) { sum += raster.data[r+1][c-1]; i+=1.0;}
  if (r < raster.nrRows()-1 && !pcr::isMV(raster.data[r+1][c  ]) && raster.data[r+1][c  ]> 0) { sum += raster.data[r+1][c  ]; i+=1.0;}
  if (r < raster.nrRows()-1 && c < raster.nrCols()-1 && !pcr::isMV(raster.data[r+1][c+1]) && raster.data[r+1][c+1]> 0) { sum += raster.data[r+1][c+1]; i+=1.0;}
  avg = (i > 0 ? sum / i : 0);

  return(avg);
}


void cover(
    cTMap& raster,
    cTMap const& value1,
    double value2)
{
    int r, c;
//#pragma omp parallel for collapse(2)
    for (r = 0; r < raster.nrRows(); r++)
        for (c = 0; c < raster.nrCols(); c++)
            if (pcr::isMV(raster.data[r][c]) && !pcr::isMV(value1.data[r][c]))
            {
                raster.data[r][c] = value2;
            }
    //     else
    //       pcr::setMV(raster.data[r][c]);
}


void calcValue(
    cTMap& raster,
    double value,
    int oper)
{
//    #pragma omp parallel for collapse(2)
    for (int r = 0; r < raster.nrRows(); r++)
        for (int c = 0; c < raster.nrCols(); c++)
            if (!pcr::isMV(raster.data[r][c]))
            {
                switch (oper)
                {
                case ADD: raster.data[r][c] += value; break;
                case SUB: raster.data[r][c] -= value; break;
                case MUL: raster.data[r][c] *= value; break;
                case DIV: if (value > 0) raster.data[r][c] /= value;
                    else pcr::setMV(raster.data[r][c]); break;
                case POW: raster.data[r][c] = std::pow(raster.data[r][c],value); break;
                case MIN: raster.data[r][c] = std::min(raster.data[r][c],value); break;//VJ 110420 new
                case MAX: raster.data[r][c] = std::max(raster.data[r][c],value); break;
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
            if (!pcr::isMV(raster.data[r][c]))
            {
                if (!pcr::isMV(value.data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: raster.data[r][c] += value.data[r][c]; break;
                    case SUB: raster.data[r][c] -= value.data[r][c]; break;
                    case MUL: raster.data[r][c] *= value.data[r][c]; break;
                    case DIV: if (value.data[r][c] > 0) raster.data[r][c] /= value.data[r][c];
                        else pcr::setMV(raster.data[r][c]); break;
                    case POW: raster.data[r][c] = powl(raster.data[r][c],value.data[r][c]); break;
                    case MIN: raster.data[r][c] = std::min(value.data[r][c], raster.data[r][c]); break; //VJ 110420 new
                    case MAX: raster.data[r][c] = std::max(value.data[r][c], raster.data[r][c]); break;
                    }
                }
                else
                    pcr::setMV(raster.data[r][c]);
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
            if (!pcr::isMV(raster.data[r][c]))
            {
                if (!pcr::isMV(value1.data[r][c]) && !pcr::isMV(value2.data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: raster.data[r][c] = value1.data[r][c] + value2.data[r][c]; break;
                    case SUB: raster.data[r][c] = value1.data[r][c] - value2.data[r][c]; break;
                    case MUL: raster.data[r][c] = value1.data[r][c] * value2.data[r][c]; break;
                    case DIV: if (value2.data[r][c] > 0) raster.data[r][c] = value1.data[r][c] / value2.data[r][c];
                        else pcr::setMV(raster.data[r][c]); break;
                    case POW: raster.data[r][c] = pow(value1.data[r][c], value2.data[r][c]); break;
                    case MIN: raster.data[r][c] = std::min(value1.data[r][c], value2.data[r][c]); break; //VJ 110420 new
                    case MAX: raster.data[r][c] = std::max(value1.data[r][c], value2.data[r][c]); break;
                    }
                }
                else
                    pcr::setMV(raster.data[r][c]);
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
            if (!pcr::isMV(raster.data[r][c]))
            {
                if (!pcr::isMV(value1.data[r][c]))
                {
                    switch (oper)
                    {
                    case ADD: raster.data[r][c] = value1.data[r][c] + value2; break;
                    case SUB: raster.data[r][c] = value1.data[r][c] - value2; break;
                    case MUL: raster.data[r][c] = value1.data[r][c] * value2; break;
                    case DIV: if (value2 > 0) raster.data[r][c] = value1.data[r][c] / value2;
                        else pcr::setMV(raster.data[r][c]); break;
                    case POW: raster.data[r][c] = pow(value1.data[r][c],value2); break;
                    case MIN: raster.data[r][c] = std::min(value1.data[r][c],value2); break;//VJ 110420 new
                    case MAX: raster.data[r][c] = std::max(value1.data[r][c],value2); break;
                    }
                }
                else
                    pcr::setMV(raster.data[r][c]);
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
            if (!pcr::isMV(raster.data[r][c]))
            {
                if (oper == LARGER && raster.data[r][c] > value)
                {
                    ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger than %4.\n").arg(r).arg(c).arg(raster.mapName()).arg(value) + SS;
                    throw 1;
                }
                else
                    if (oper == SMALLER && raster.data[r][c] < value)
                    {
                        ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller than %4.\n").arg(r).arg(c).arg(raster.mapName()).arg(value) + SS;
                        throw 1;
                    }
                    else
                        if (oper == LARGEREQUAL && raster.data[r][c] >= value)
                        {
                            ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger or equal than %4.\n").arg(r).arg(c).arg(raster.mapName()).arg(value) + SS;
                            throw 1;
                        }
                        else
                            if (oper == SMALLEREQUAL && raster.data[r][c] <= value)
                            {
                                ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller or equal than %4.\n").arg(r).arg(c).arg(raster.mapName()).arg(value) + SS;
                                throw 1;
                            }
            }
}
