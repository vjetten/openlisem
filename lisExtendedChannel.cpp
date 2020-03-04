
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

#include <algorithm>
#include "operation.h"
#include "model.h"


void TWorld::ExtendChannel()
{

    ChannelDepthExtended = NewMap(0.0);
    ChannelWidthExtended = NewMap(0.0);
    ChannelMaskExtended = NewMap(0.0);
    ChannelFlowWidth = NewMap(0.0);

    if(!SwitchIncludeChannel)
    {
        return;
    }

    copy(*ChannelWidthExtended, *ChannelWidthMax);
    copy(*ChannelDepthExtended, *ChannelDepth);

    ChannelNeighborsExtended = NewMap(0.0);
    ChannelSourceXExtended = NewMap(0.0);
    ChannelSourceYExtended = NewMap(0.0);
    ChannelBoundaryExtended = NewMap(0.0);
    ChannelBoundaryLExtended = NewMap(0.0);
    ChannelBoundaryRExtended = NewMap(0.0);


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

    //for checking wheter channel lies within cell
    //double maxd = (std::max(0.0,MaxWidth-_dx)/2.0);

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
                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);

                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
                            }
                        }
                    }
                }
            }
            for(int r2 = r - i; r2 < r + i + 1; r2++)
            {
                int c2 = c+i;

                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);
                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
                            }
                        }
                    }
                }
            }
            for(int c2 = c - i + 1; c2 < c + i; c2++)
            {
                int r2 = r - i;
                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //`double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);

                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
                            }
                        }
                    }
                }
            }
            for(int c2 = c - i + 1; c2 < c + i; c2++)
            {
                int r2 = r+i;
                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);

                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
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

    FOR_ROW_COL_MV
    {
        if(ChannelMaskExtended->Drc == 1)
        {
            if(!pcr::isMV(LDDChannel->Drc) && !(ChannelWidthMax->Drc > _dx))
            {

                    ChannelBoundaryExtended->Drc = 1;
                    ChannelBoundaryRExtended->Drc = 1;
                    ChannelBoundaryLExtended->Drc = 1;
            }

            double nn = 0;

            if(IsExtendedChannel(r,c,1,0)) nn+=1.0;
            if(IsExtendedChannel(r,c,-1,0)) nn+=1.0;
            if(IsExtendedChannel(r,c,0,1)) nn+=1.0;
            if(IsExtendedChannel(r,c,0,-1)) nn+=1.0;
            if(nn == 4){
                ChannelBoundaryExtended->Drc = 0;
                ChannelBoundaryRExtended->Drc = 0;
                ChannelBoundaryLExtended->Drc = 0;
            }
//            if(!ChannelBoundaryExtended->Drc == 1 && nn < 4)
                if(ChannelBoundaryExtended->Drc != 1 && nn < 4)
            {
                ChannelBoundaryExtended->Drc = 1;
                ChannelBoundaryRExtended->Drc = 1;
                ChannelBoundaryLExtended->Drc = 1;
            }
            if(IsExtendedChannel(r,c,-1,1)) nn+=1.0;
            if(IsExtendedChannel(r,c,-1,-1)) nn+=1.0;
            if(IsExtendedChannel(r,c,1,1)) nn+=1.0;
            if(IsExtendedChannel(r,c,1,-1)) nn+=1.0;
        }
    }


    return;
}

bool TWorld::IsExtendedChannel(int r, int c, int dr, int dc)
{
    if(!OUTORMV(r+dr,c+dc))
    {
        return ChannelMaskExtended->data[r+dr][c+dc] == 1;
    }else
    {
        return true;
    }

}

//---------------------------------------------------------------------------
// Distributes a certain value over the actual channel width (used for display stuff)

void TWorld::DistributeOverExtendedChannel(cTMap * _In, cTMap * _Out, bool do_not_divide,bool proportional)
{
    double totala=0;
    double totalb=0;
    FOR_ROW_COL_MV
    {
        if(!pcr::isMV(LDDChannel->Drc))
        {
            totala += _In->Drc;
        }
        if(ChannelMaskExtended->Drc == 1)
        {
            double ow = ChannelWidthMax->data[(int)ChannelSourceYExtended->Drc][(int)ChannelSourceXExtended->Drc];
            if(ow> 0 )
            {
                double div = do_not_divide? (proportional? ChannelWidthExtended->Drc/ _dx : 1.0): (ChannelWidthExtended->Drc / ow);
                _Out->Drc = _In->data[(int)ChannelSourceYExtended->Drc][(int)ChannelSourceXExtended->Drc] * div;
            }
            totalb += _Out->Drc;
        }else
        {
            _Out->Drc  = 0.0;
        }
    }
    if(totalb > 0)
    {
        FOR_ROW_COL_MV
        {
            if(ChannelMaskExtended->Drc == 1)
            {
                _Out->Drc *= totala/totalb;
            }
        }
    }
}
