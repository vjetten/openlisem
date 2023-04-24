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
#include "operation.h"
#include "model.h"

//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

#define do_row(r, c, w2, n)   if(n == 1) extendRow(r,c,1,w2);\
else for (int i = 1; i < n; i++) extendRow(r,c,i,(i == n-1) ? w2 : adx);tma->data[r][c]=2

#define do_col(r,c, w2, n)   if(n == 1) extendCol(r,c,1,w2);\
else for (int i = 1; i < n; i++) extendCol(r,c,i,(i == n-1) ? w2 : adx);tma->data[r][c]=3



/// when channelwidth > _dx distributes excess width over neighbouring cells
/**
* returns ChannelWidthExtended, ChannelDepthExtended,
* ChannelMaskExtended : mask = 1 for all cells in the extended channel
* ChannelSourceXExtended : row of the channel cell the extended cells belong to
* ChannelSourceXExtended : col of the channel cell the extended cells belong to
* simple assumption: extension is always in X or Y, not diagonally. Diagonal LDD cells
* are just a direction of flow
*/

void TWorld::doExtendRow(int r, int c, int n,  double w2, double adx)
{
    ExtendCH.chCol = c;
    ExtendCH.chRow = r;
    ExtendCH.isExtended = true;
    ExtendCH.childCol.clear();
    ExtendCH.childRow.clear();
    if (n == 1) {
        extendRow(r,c,1,w2);
        ExtendCH.childCol << c;
        ExtendCH.childRow << r-1;
        ExtendCH.childCol << c;
        ExtendCH.childRow << r+1;
    }
    else
        for (int i = 1; i < n; i++) {
           // w = ;
            extendRow(r,c,i,(i == n-1) ? w2 : adx);
            ExtendCH.childCol << c;
            ExtendCH.childRow << r-i;
            ExtendCH.childCol << c;
            ExtendCH.childRow << r+i;
        }
    ExtChannel << ExtendCH;
    tma->data[r][c]=2;
}


void TWorld::doExtendCol(int r, int c, int n, double w2, double adx)
{
    ExtendCH.chCol = c;
    ExtendCH.chRow = r;
    ExtendCH.isExtended = true;
    ExtendCH.childCol.clear();
    ExtendCH.childRow.clear();
    if (n == 1) {
        extendCol(r,c,1,w2);
        ExtendCH.childCol << c-1;
        ExtendCH.childRow << r;
        ExtendCH.childCol << c+1;
        ExtendCH.childRow << r;
    }
    else
        for (int i = 1; i < n; i++) {
           // w = (i == n-1) ? w2 : adx;
            extendCol(r,c,i,(i == n-1) ? w2 : adx);
            ExtendCH.childCol << c-i;
            ExtendCH.childRow << r;
            ExtendCH.childCol << c+i;
            ExtendCH.childRow << r;
        }
    tma->data[r][c]=3;
    ExtChannel << ExtendCH;
}


void TWorld::extendRow(int r, int c, int i, double w)
{

    if (notMVIn(r-1,c) && tma->data[r-i][c] <= 0){
        ChannelWidthExtended->data[r-i][c] = std::max(ChannelWidthExtended->data[r-i][c],w);
        if (ChannelWidthExtended->data[r-i][c] > 0) tma->data[r-i][c] = 1;
        ChannelSourceYExtended->data[r-i][c] = r;
        ChannelSourceXExtended->data[r-i][c] = c;
    }
    if (notMVIn(r+1,c) && tma->data[r+i][c] <= 0){
        ChannelWidthExtended->data[r+i][c] = std::max(ChannelWidthExtended->data[r+i][c],w);
        if (ChannelWidthExtended->data[r+i][c] > 0) tma->data[r+i][c] = 1;
        ChannelSourceYExtended->data[r+i][c] = r;
        ChannelSourceXExtended->data[r+i][c] = c;
    }
}

void TWorld::extendCol(int r, int c, int i, double w)
{
    if (notMVIn(r,c-1) && tma->data[r][c-i] <= 0){
        ChannelWidthExtended->data[r][c-i] = std::max(ChannelWidthExtended->data[r][c-i],w);
        if (ChannelWidthExtended->data[r][c-i] > 0) tma->data[r][c-i] = 1;
        ChannelSourceYExtended->data[r][c-i] = r;
        ChannelSourceXExtended->data[r][c-i] = c;
    }
    if (notMVIn(r,c+1)  && tma->data[r][c+i] <= 0){
        ChannelWidthExtended->data[r][c+i] = std::max(ChannelWidthExtended->data[r][c+i],w);
        if (ChannelWidthExtended->data[r][c+i] > 0) tma->data[r][c+i] = 1;
        ChannelSourceYExtended->data[r][c+i] = r;
        ChannelSourceXExtended->data[r][c+i] = c;
    }
}


bool TWorld::ExtendChannelNew()
{
    ChannelDepthExtended = NewMap(0.0);
    ChannelWidthExtended = NewMap(0.0);
    ChannelMaskExtended = NewMap(0.0);
    ChannelSourceXExtended = NewMap(0.0);
    ChannelSourceYExtended = NewMap(0.0);
    ChannelWHExtended = NewMap(0.0);
    ChannelVolExtended = NewMap(0.0);

    ExtChannel.clear();

    if(!SwitchIncludeChannel)
    {
        return false;
    }

    double ww = 0;
    FOR_ROW_COL_MV_CH {
        ww = std::max(ChannelWidthMax->Drc, ww);
        ChannelWidthExtended->Drc = ChannelWidth->Drc;
        ChannelDepthExtended->Drc = ChannelDepth->Drc;
        ChannelMaskExtended->Drc = 1;
        ChannelSourceXExtended->Drc = c;
        ChannelSourceYExtended->Drc = r;
        ExtendCH.chCol = c;
        ExtendCH.chRow = r;
        ExtendCH.isExtended = false;
        ExtChannel << ExtendCH;
    }
    if (ww < _dx)
        return false; // there is no extended channel


    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1,  0,  1};
    int dy[10] = {0,  1, 1, 1,  0, 0, 0, -1, -1, -1};
    Fill(*tma, 0); // flag, -1 channel cell is not done, 2 row direction, 3 col direction
    double adx = 0.95*_dx;
    // do first batch horizontal, vertical and diagonsl
    FOR_ROW_COL_MV_CH
    {
        double w = ChannelWidthMax->Drc;

        if (w < adx) {
            ChannelWidthExtended->Drc = w;
            ChannelDepthExtended->Drc = ChannelDepth->Drc;
            ChannelMaskExtended->Drc = 1;
            ChannelSourceXExtended->Drc = c;
            ChannelSourceYExtended->Drc = r;
            ExtendCH.chCol = c;
            ExtendCH.chRow = r;
            ExtendCH.isExtended = false;
            ExtChannel << ExtendCH;
            tma->Drc = 1;
            continue;
        }

        int n = 1;
        while ((w-n*adx)/adx > 2)
                    n++;
        double dif = w - (n*adx);
        if (n % 2 == 0){ dif -= adx; n++; }
        double w2 = dif/2.0;

        ChannelWidthExtended->Drc = adx;
        tma->Drc = -1;

        if (LDDChannel->Drc == 2 || LDDChannel->Drc == 8) {
           doExtendCol(r,c,n,w2,adx); //tma = 3
           // do_col(r, c,  w2, n);
        }
        if (LDDChannel->Drc == 4 || LDDChannel->Drc == 6) {
            doExtendRow(r,c,n,w2,adx); //tma = 2
          //  do_row(r, c,  w2, n);
        }

        // check inflow cells are horizontal or vertical
        QList <int> lddfrom;
        // do diagonal cells and pits, 1, 3, 5, 7, 9
        if ((int)LDDChannel->Drc % 2 > 0) {

            int ldd0 = (int) LDDChannel->Drc;
            int lddto = (int) LDDChannel->data[r+dy[ldd0]][c+dx[ldd0]];
            lddfrom.clear();
            for (int i=1; i<=9; i++)
            {
                int ldd = 0;
                if (i==5)
                    continue;
                int rr = r+dy[i];
                int cr = c+dx[i];

                if (notMVIn(rr,cr) && !pcr::isMV(LDDChannel->Drcr))
                    ldd = (int) LDDChannel->Drcr;
                else
                    continue;

                if (FLOWS_TO(ldd, rr, cr, r, c))
                   if (ldd % 2 == 0) lddfrom << ldd;
            }
            int lddfromi = 0;
            for (int j = 0; j < lddfrom.length(); j++)
                 lddfromi = lddfrom.at(j);

            if (lddfromi == 2 || lddto == 2 || lddfromi == 8 || lddto == 8) {
                doExtendCol(r,c,n,w2,adx);
            }

            if (lddfromi == 4 || lddto == 4 || lddfromi == 6 || lddto == 6) {
                doExtendRow(r,c,n,w2,adx);
            }
        }
    }

    // loop through all channel cells if there are still cells not assigned
    // flagged -1. If so, give these the direction of the nearest horizontal or vertical
    // cell flagged 3 or 2
    int flag = 0;
    FOR_ROW_COL_MV_CH
    {
        if (tma->Drc == -1)
            flag = -1;
    }

    while (flag < 0) {
        FOR_ROW_COL_MV_CH
        {
            if (tma->Drc == -1) {
                double w = ChannelWidthMax->Drc;

                int n = 1;
                while ((w-n*adx)/adx > 2)
                    n++;
                double dif = w - (n*adx);
                if (n % 2 == 0){ dif -= adx; n++; }
                double w2 = dif/2.0;

                int ldd0 = (int) LDDChannel->Drc;
                int lddto = (int) LDDChannel->data[r+dy[ldd0]][c+dx[ldd0]];
                QList <int> lddfrom;
                lddfrom.clear();
                for (int i=1; i<=9; i++)
                {
                    int ldd = 0;
                    if (i==5)
                        continue;
                    int rr = r+dy[i];
                    int cr = c+dx[i];

                    if (notMVIn(rr,cr) && !pcr::isMV(LDDChannel->Drcr))
                        ldd = (int) LDDChannel->Drcr;
                    else
                        continue;

                    if (FLOWS_TO(ldd, rr, cr, r, c))
                        lddfrom << ldd;
                }
                int m = 0;
                if (tma->data[r+dy[ldd0]][c+dx[ldd0]] == 3) m++;
                if (tma->data[r+dy[lddto]][c+dx[lddto]] == 3) m++;
                if (tma->data[r+dy[lddfrom.at(0)]][c+dx[lddfrom.at(0)]] == 3) m++;

                if (m >= 2) {
                    doExtendCol(r,c,n,w2,adx);
                } else {
                    doExtendRow(r,c,n,w2,adx);
                }
            }
        }
        flag = 0;
        FOR_ROW_COL_MV_CH
        {
            if (tma->Drc < 0)
                flag = -1;
        }
    }

    //make sure lddchannel cells are included
    FOR_ROW_COL_MV {
        if (ChannelWidthExtended->Drc > 0)
            ChannelMaskExtended->Drc = 1;
        if (LDDChannel->Drc > 0) {
            ChannelSourceXExtended->Drc = c;
            ChannelSourceYExtended->Drc = r;
        }
    }

    //report(*tma,"tma.map");
    //report (*ChannelWidthExtended,"chwe.map");
    copy(*ChannelDepthExtended, *ChannelDepth);
   // report (*ChannelSourceXExtended,"chmx.map");
   // report (*ChannelSourceYExtended,"chmy.map");

    return true;
}

void TWorld::ExtendChannel()
{
    ExtendChannelNew();
    return;


    // channelwidth does add up, channeldepth is altered
    ChannelDepthExtended = NewMap(0.0);
    ChannelWidthExtended = NewMap(0.0);
    ChannelMaskExtended = NewMap(0.0);
   // ChannelFlowWidth = NewMap(0.0);


    if(!SwitchIncludeChannel)
    {
        return;
    }

    copy(*ChannelWidthExtended, *ChannelWidthMax);
    copy(*ChannelDepthExtended, *ChannelDepth);

    ChannelNeighborsExtended = NewMap(0.0);
    ChannelSourceXExtended = NewMap(0.0);
    ChannelSourceYExtended = NewMap(0.0);


    //Channel cells are part of the extended channel
    FOR_ROW_COL_MV_CH
    {
        ChannelMaskExtended->Drc = (ChannelWidthMax->Drc > 0 ? 1.0 : 0.0);
    }

    int dxl[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dyl[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    //set width etc for channel cells, for these we know they are within the extende channel
    double MaxWidth = 0;
    FOR_ROW_COL_MV_CH
    {
        if(ChannelWidthMax->Drc > MaxWidth)
        {
            MaxWidth = ChannelWidthMax->Drc;
        }
        ChannelMaskExtended->Drc = 1;
        ChannelDepthExtended->Drc = ChannelDepth->Drc;
        ChannelWidthExtended->Drc = std::min(_dx,ChannelWidthMax->Drc);
        ChannelNeighborsExtended->Drc = -1;
    }

    //for iteration
    int maxdistance = ((int) ((std::max(0.0,MaxWidth-_dx)/2.0)/_dx)) + 1;

    FOR_ROW_COL_MV
    {
        if(!pcr::isMV(LDDChannel->Drc))
        {
            ChannelSourceXExtended->Drc = c;
            ChannelSourceYExtended->Drc = r;
            continue;
        }

        bool found = false;
        double found_distance = 9999999.9;
        int i = 1;
        while(i < maxdistance + 1 && !(found_distance < double(i) * _dx))
        {
            for(int r2 = r - i; r2 < r + i + 1; r2++)
            {
                int c2 = c - i;
                if(notMVIn(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            double width =  std::min(_dx,std::max(0.0, ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += 0.5*width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                            }
                        }
                    }
                }
            }
            for(int r2 = r - i; r2 < r + i + 1; r2++)
            {
                int c2 = c+i;

                if(notMVIn(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            double width =  std::min(_dx,std::max(0.0, ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += 0.5*width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                            }
                        }
                    }
                }
            }
            for(int c2 = c - i + 1; c2 < c + i; c2++)
            {
                int r2 = r - i;
                if(notMVIn(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            double width =  std::min(_dx,std::max(0.0, ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += 0.5*width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;

                            }
                        }
                    }
                }
            }
            for(int c2 = c - i + 1; c2 < c + i; c2++)
            {
                int r2 = r+i;
                if(notMVIn(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            double width =  std::min(_dx,std::max(0.0, ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += 0.5*width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                            }
                        }
                    }
                }
            }

            i++;
        }

        if(found)
        {
            if(ChannelWidthExtended->Drc > 0 && ChannelNeighborsExtended->Drc > 0) //VJ!!!! ChannelWidthExtended without Drc
            {
                ChannelDepthExtended->Drc /= ChannelWidthExtended->Drc;
                ChannelWidthExtended->Drc /= ChannelNeighborsExtended->Drc;
            }
        }
    }


//    report (*ChannelWidthExtended,"chma.map");
//    report (*ChannelDepthExtended,"chda.map");

//    report (*ChannelSourceXExtended,"chmxa.map");
//    report (*ChannelSourceYExtended,"chmya.map");
    return;
}

bool TWorld::IsExtendedChannel(int r, int c, int dr, int dc)
{
    if(notMVIn(r+dr,c+dc))
    {
        return ChannelMaskExtended->data[r+dr][c+dc] == 1;
    }else
    {
        return true;
    }

}

//---------------------------------------------------------------------------
// Distributes a certain value over the actual channel width (used for display stuff)
//bool do_not_divide = false,bool proportional = true)
void TWorld::DistributeOverExtendedChannel(cTMap * _In, cTMap * _Out)
{
    /*
    if(!SwitchIncludeChannel)
        return;
    double intot = MapTotal(*_In);
    Fill(*_Out, 0.0);

    if(!SwitchChannelExtended) {
#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            _Out->Drc = _In->Drc;
        }}
        return;
    }

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if(ChannelMaskExtended->Drc == 1) {
            int rr = (int)ChannelSourceYExtended->Drc;
            int cr = (int)ChannelSourceXExtended->Drc;
            double div = ChannelWidthExtended->Drc/ChannelWidthMax->Drcr;
            _Out->Drc = _In->Drcr * div;
        //    qDebug() << rr << cr << div << ChannelWidthExtended->Drc << ChannelWidthMax->Drcr;
        }
    }}

    double outtot = mapTotal(*_Out);
    if(outtot > 0) {
#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            if(ChannelMaskExtended->Drc == 1) {
                _Out->Drc *= intot/outtot;
            }
        }}
    }
    outtot = MapTotal(*_Out);
    //qDebug() << intot << outtot;
    */
}
