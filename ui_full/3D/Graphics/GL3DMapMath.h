#ifndef GL3DMAPMATH_H
#define GL3DMAPMATH_H

#include "3D/GL3DWidget.h"

class MapMath
{

public:

    inline static void SlopeMap(cTMap * in_dem, cTMap * out_slope_x,cTMap * out_slope_y)
    {

        FOR_ROW_COL_MV(in_dem,in_dem)
        {
            double slope = 0;
            double Dhx = 0;
            double Dhy = 0;

            //DEM
            double dem = in_dem->data[r][c];

            double demx1 = (INSIDE(in_dem,r,c+1))? (!MV(in_dem,r,c+1)? in_dem->data[r][c+1] : dem):dem;
            double demx2 = (INSIDE(in_dem,r,c-1))? (!MV(in_dem,r,c-1)? in_dem->data[r][c-1] : dem):dem;

            double demy1 = (INSIDE(in_dem,r+1,c))? (!MV(in_dem,r+1,c)? in_dem->data[r+1][c] : dem):dem;
            double demy2 = (INSIDE(in_dem,r-1,c))? (!MV(in_dem,r-1,c)? in_dem->data[r-1][c] : dem):dem;

            out_slope_x->data[r][c] = (demx1 -demx2)/in_dem->cellSize();
            out_slope_y->data[r][c] = (demy1 -demy2)/in_dem->cellSize();

        }
    }

    inline static void FillDem(cTMap * in_dem, cTMap * out_dem,cTMap * temp, int iter)
    {
        out_dem->setAllMV();
        temp->setAllMV();
        FOR_ROW_COL_MV(in_dem,in_dem)
        {
            out_dem->Drc = in_dem->Drc;
            temp->Drc = in_dem->Drc;
        }

        for(int i = 0; i < iter; i++)
        {

            FOR_ROW_COL_MV(out_dem,out_dem)
            {
                temp->Drc = out_dem->Drc;
            }

            FOR_ROW_COL_MV(temp,temp)
            {

                int dx[8] = {0, 1, 0, -1,1, -1, 1, -1};
                int dy[8] = {1, 0, -1, 0,1, -1, -1, 1};

                for(int l = 0; l < 4; l++)
                {
                    int r2 = r + dy[l];
                    int c2 = c + dx[l];

                    if(INSIDE(temp,r2,c2))
                    {
                        if(MV(temp,r2,c2))
                        {
                            out_dem->data[r2][c2] = out_dem->Drc;
                        }

                    }
                }

            }
        }

    }

};


#endif // GL3DMAPMATH_H
