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

    inline static void FlowSlopeMap(cTMap * in_dem, cTMap * in_flow, cTMap * out_slope_x,cTMap * out_slope_y)
    {

        FOR_ROW_COL_MV(in_dem,in_dem)
        {
            double slope = 0;
            double Dhx = 0;
            double Dhy = 0;

            //DEM
            double dem = in_dem->data[r][c] + in_flow->data[r][c];

            double demx1 = (INSIDE(in_dem,r,c+1))? (!MV(in_dem,r,c+1)? in_dem->data[r][c+1]+in_flow->data[r][c+1]: dem):dem;
            double demx2 = (INSIDE(in_dem,r,c-1))? (!MV(in_dem,r,c-1)? in_dem->data[r][c-1]+in_flow->data[r][c-1]: dem):dem;

            double demy1 = (INSIDE(in_dem,r+1,c))? (!MV(in_dem,r+1,c)? in_dem->data[r+1][c]+in_flow->data[r+1][c]: dem):dem;
            double demy2 = (INSIDE(in_dem,r-1,c))? (!MV(in_dem,r-1,c)? in_dem->data[r-1][c]+in_flow->data[r-1][c]: dem):dem;

            out_slope_x->data[r][c] = (demx1 -demx2)/in_dem->cellSize();
            out_slope_y->data[r][c] = (demy1 -demy2)/in_dem->cellSize();

        }
    }


    inline static bool GetMVAt(cTMap * map, double x, double z)
    {
        double cs = map->cellSize();
        int cn = std::floor(x/cs);
        int rn = std::floor(z/cs);
        return pcr::isMV(map->data[rn][cn]);

    }

    //bilinear interpolation, mimics opegl
    inline static double GetValueAt(cTMap * map, double x, double z)
    {
        double cs = map->cellSize();
        x = x - 1.0 * cs;
        z = z - 1.0 * cs;

        int cn = std::floor((x)/cs);
        int rn = std::floor((z)/cs);

        int sgnx = 1;
        int sgny = 1;

        double wx =1.0-fabs((x - (cn * cs))/cs);
        double wy =1.0-fabs((z - (rn * cs))/cs);

        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
        double v4 = 0;

        double w1 = 0;
        double w2 = 0;
        double w3 = 0;
        double w4 = 0;

        if(!OUTOFMAP(map,rn,cn))
        {
            if(!MV(map,rn,cn))
            {
                w1 = (wx) * (wy);
                v1 = map->data[rn][cn];
            }
        }
        if(!OUTOFMAP(map,rn + sgnx *1,cn))
        {
            if(!MV(map,rn + sgnx *1,cn))
            {
                w2 = (wx) * (1.0-wy);
                v2 =map->data[rn + sgnx *1][cn];
            }
        }
        if(!OUTOFMAP(map,rn,cn + sgny *1))
        {
            if(!MV(map,rn,cn + sgny *1))
            {
                w3 = (1.0-wx) * (wy);
                v3 =map->data[rn][cn + sgny *1];
            }
        }
        if(!OUTOFMAP(map,rn + sgnx *1,cn + sgny *1))
        {
            if(!MV(map,rn + sgnx *1,cn + sgny *1))
            {
                w4 = (1.0-wx) * (1.0-wy);
                v4 =map->data[rn + sgnx *1][cn + sgny *1];
            }
        }

        return (w1+w2+w3+w4) > 0? ((w1*v1 + w2*v2 + w3*v3 + w4 * v4)/(w1+w2+w3+w4)): 0.0;


    }

};


#endif // GL3DMAPMATH_H
