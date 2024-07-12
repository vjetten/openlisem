/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/


#include <algorithm>
#include <qstring.h>
#include "io.h"
#include "lisemqt.h"
#include "global.h"

#include "model.h"
#include "operation.h"
#include "CsfRGBMap.h"

//---------------------------------------------------------------------------
/** \n void TWorld::InitMapList(void)
*
*/
void TWorld::InitMapList(void)
{
    maplistnr = 0;
    for (int i = 0; i < NUMNAMES; i++)
    {
        maplistCTMap[i].m = nullptr;
    }
}
//---------------------------------------------------------------------------
cTMap *TWorld::NewMap(double value)
{
    cTMap *_M = new cTMap();

    _M->MakeMap(LDD, value);
    // changed to LDD instead of Mask

    if (_M)
    {
        maplistCTMap[maplistnr].m = _M;
        maplistnr++;
    } else
        qDebug() << "no more space";

    return(_M);
}
//---------------------------------------------------------------------------
//NOT USED
cTMap *TWorld::ReadFullMap(QString name)
{
    cTMap *_M = new cTMap(readRaster(name));

    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
            if (pcr::isMV(_M->Drc))
            {
                _M->Drc = 0;
            }

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    return(_M);

}
//---------------------------------------------------------------------------
// read a map from disk
cTMap *TWorld::ReadMap(cTMap *Mask, QString name)
{
    cTMap *_M = new cTMap(readRaster(name));

    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
            if (!pcr::isMV(Mask->Drc) && pcr::isMV(_M->Drc))
            {
                QString sr, sc;
                sr.setNum(r); sc.setNum(c);
                ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+name+".\n \
                                                                                                     This is a cell with missing values where a flow network esists (either LDD, Channel LDD, tile drain LDD).";
                                                                                                     throw 1;
            }

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    return(_M);

}
//---------------------------------------------------------------------------
void TWorld::DestroyData(void)
{
    if (op.nrRunsDone <= 1)
        return;
qDebug() << "destroy" << op.nrMapsCreated;
    DEBUG("clear all maps");
    if (op.nrMapsCreated > 0) {
        for (int i = 0; i < op.nrMapsCreated; i++)
        {
            if (maplistCTMap[i].m != nullptr)
            {
                delete maplistCTMap[i].m;
                maplistCTMap[i].m = nullptr;
            }
        }
    }
qDebug() << "rain";
DEBUG("clear meteo structures");

    // clear() calls the destruction of all elements in the sturcture
    //if (SwitchRainfall) {
        RainfallSeries.clear();
        RainfallSeriesMaps.clear();
        calibRainfallinFile = false;
    //}
    //if (SwitchSnowmelt) {
        SnowmeltSeries.clear();
        SnowmeltSeriesMaps.clear();
    //}
    //if (SwitchIncludeET) {
        ETSeries.clear();
        ETSeriesMaps.clear();
    //}

    if (InfilMethod == INFIL_SWATRE && initSwatreStructure)
    {
        DEBUG("clear swatre structure");

        FreeSwatreInfo();
        if (SwatreSoilModel)
            CloseSwatre(SwatreSoilModel);
        if (SwatreSoilModelCrust)
            CloseSwatre(SwatreSoilModelCrust);
        if (SwatreSoilModelCompact)
            CloseSwatre(SwatreSoilModelCompact);
        if (SwatreSoilModelGrass)
            CloseSwatre(SwatreSoilModelGrass);
    }
    qDebug() << "network";

    DEBUG("clear network structures");
    cr_.clear();
    crch_.clear();
    crldd5_.clear();
    crlddch5_.clear();

    for(int i_ = 0; i_ < crlinkedldd_.size(); i_++){
        if(crlinkedldd_[i_].inn)
            free(crlinkedldd_[i_].inn);
    }
    crlinkedldd_.clear();

    for(int i_ = 0; i_ < crlinkedlddch_.size(); i_++){
        if(crlinkedlddch_[i_].inn)
            free(crlinkedlddch_[i_].inn);
    }
    crlinkedlddch_.clear();

}
//---------------------------------------------------------------------------
/// separate networks need their own InitMask: LDD, ChannelLDD, TileLDD
cTMap *TWorld::InitMask(QString name)
{
    // read map and make a mask map

    cTMap *_M = new cTMap(readRaster(name));

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    _dx = _M->cellSize()*1.0000000;
    _nrRows = _M->nrRows();
    _nrCols = _M->nrCols();
    _llx =  _M->west();
    _lly = _M->north() - (double)_nrRows * _dx;

    return(_M);

}
//---------------------------------------------------------------------------
/// separate networks need their own InitMask: LDD, ChannelLDD, TileLDD
cTMap *TWorld::InitMaskChannel(QString name)
{

    cTMap *_M = new cTMap(readRaster(name));

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    return(_M);

}
//---------------------------------------------------------------------------
/// separate networks need their own InitMask: LDD, ChannelLDD, TileLDD
cTMap *TWorld::InitMaskTiledrain(QString name)
{

    cTMap *_M = new cTMap(readRaster(name));

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    return(_M);

}
//---------------------------------------------------------------------------
void TWorld::Fill(cTMap &M, double value)
{
#pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        M.Drc = value;
    }}
}
//---------------------------------------------------------------------------
void TWorld::Copy(cTMap &M, cTMap &M1)
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        M1.Drc = M.Drc;
    }}
}
//---------------------------------------------------------------------------
double TWorld::MapTotal(cTMap &M)
{
    double total = 0;
#pragma omp parallel for reduction(+:total) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (!pcr::isMV(M.Drc))
            total = total + M.Drc;
    }}
return (total);
}
//---------------------------------------------------------------------------
void TWorld::Average3x3(cTMap &M, cTMap &mask, bool only)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tm->Drc = M.Drc;
    }}

    FOR_ROW_COL_MV_L {
        double tot = 0;
        double cnt = 0;
        for (int i = 1; i <= 9; i++)
        {
            int rr = r+dy[i];
            int cr = c+dx[i];

            if (INSIDE(rr, cr) && !pcr::isMV(mask.Drcr)) {
                if (only && M.Drcr == 0)
                    continue;
                tot = tot + tm->Drcr;
                cnt += 1.0;
                if (i == 5) {
                    tot = tot + tm->Drcr;
                    cnt += 1.0;
                }
            }
        }
        M.Drc = cnt > 0 ? tot/cnt : tm->Drc;
        if (pcr::isMV(mask.Drc))
            M.Drc = tm->Drc;
    }}
}
//---------------------------------------------------------------------------
void TWorld::Average2x2(cTMap &M, cTMap &mask)
{
    int dx[10] = {0, -1, 1, -1,  1};
    int dy[10] = {0,  1, 1, -1, -1};
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tm->Drc = M.Drc;
    }}

    double f = 0.5;
    FOR_ROW_COL_MV_L {
        double tot = 0;
        double cnt = 0;
        for (int i = 0; i <= 5; i++)
        {
            int rr = r+dy[i];
            int cr = c+dx[i];

            if (INSIDE(rr, cr) && !pcr::isMV(mask.Drcr)) {
                tot = tot + tm->Drcr;
                cnt += 1.0;
            }
        }
        M.Drc = cnt > 0 ? tot/cnt : tm->Drc;
    }}
}
//---------------------------------------------------------------------------
//NOT USED
double TWorld::LogNormalDist(double d50,double s, double d)
{
    double dev = log(1.0 + s/d50);
    double dev2 = (log(d)  - log(d50));
    return (1.0/(d *sqrt(2.0*3.14159) * log(1.0 + s/d50)))*exp(-dev2*dev2)/(4*dev*dev);

}
