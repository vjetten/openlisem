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
 \file lisUnifiedFlowMuscle.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

void TWorld::UF2D_MUSCLE(cTMap * _dem,cTMap * dt, cTMap * _in, int target, int dir)
{
     double delta_x1, delta_y1;
     double delta_x2, delta_y2;
     double dx, dy;

     cTMap * outx1 =0;
     cTMap * outy1 =0;
     cTMap * outx2 =0;
     cTMap * outy2 =0;

     if(target == UF_MUSCLE_TARGET_OUT)
     {
         outx1 = UF2D_MUSCLE_OUT_x1;
         outx2 = UF2D_MUSCLE_OUT_x2;
         outy1 = UF2D_MUSCLE_OUT_y1;
         outy2 = UF2D_MUSCLE_OUT_y2;

     }else if(target == UF_MUSCLE_TARGET_IN1)
     {
         outx1 = UF2D_MUSCLE_1_x1;
         outx2 = UF2D_MUSCLE_1_x2;
         outy1 = UF2D_MUSCLE_1_y1;
         outy2 = UF2D_MUSCLE_1_y2;

     }else if(target == UF_MUSCLE_TARGET_IN2)
     {
         outx1 = UF2D_MUSCLE_2_x1;
         outx2 = UF2D_MUSCLE_2_x2;
         outy1 = UF2D_MUSCLE_2_y1;
         outy2 = UF2D_MUSCLE_2_y2;
     }else if(target == UF_MUSCLE_TARGET_IN3)
     {
         outx1 = UF2D_MUSCLE_3_x1;
         outx2 = UF2D_MUSCLE_3_x2;
         outy1 = UF2D_MUSCLE_3_y1;
         outy2 = UF2D_MUSCLE_3_y2;
     }

     if(dt != 0)
     {
         if(dir == UF_DIRECTION_X || dir == UF_DIRECTION_XY)
         {
             FOR_ROW_COL_UF2D_DT
             {
                 delta_x1 = (UF_OUTORMV(_dem,r,c-1)? 0.0: _in->Drc - _in->data[r][c-1])/_dx ;
                 delta_x2 = (UF_OUTORMV(_dem,r,c+1)? 0.0: _in->data[r][c+1] - _in->Drc)/_dx;
                 dx   = 0.5*UF_MinMod(delta_x1, delta_x2);

                 outx1->Drc = _in->Drc-0.5 * dx;
                 outx2->Drc = _in->Drc+0.5 * dx;

             }}}
         }
         if(dir == UF_DIRECTION_Y || dir == UF_DIRECTION_XY)
         {
             FOR_ROW_COL_UF2D_DT
             {
                 delta_y1 = (UF_OUTORMV(_dem,r-1,c)? 0.0: _in->Drc - _in->data[r-1][c])/_dx;
                 delta_y2 = (UF_OUTORMV(_dem,r+1,c)? 0.0: _in->data[r+1][c] - _in->Drc)/_dx;
                 dy   = 0.5*UF_MinMod(delta_y1, delta_y2);

                 outy1->Drc = _in->Drc-0.5 * dy;
                 outy2->Drc = _in->Drc+0.5 * dy;
            }}}
         }
    }else
    {
        if(dir == UF_DIRECTION_X || dir == UF_DIRECTION_XY)
        {
            FOR_ROW_COL_UF2D
            {
                delta_x1 = (UF_OUTORMV(_dem,r,c-1)? 0.0: _in->Drc - _in->data[r][c-1])/_dx ;
                delta_x2 = (UF_OUTORMV(_dem,r,c+1)? 0.0: _in->data[r][c+1] - _in->Drc)/_dx;
                dx   = 0.5*UF_MinMod(delta_x1, delta_x2);

                outx1->Drc = _in->Drc-0.5 * dx;
                outx2->Drc = _in->Drc+0.5 * dx;
            }
        }
        if(dir == UF_DIRECTION_Y || dir == UF_DIRECTION_XY)
        {
            FOR_ROW_COL_UF2D
            {
                delta_y1 = (UF_OUTORMV(_dem,r-1,c)? 0.0: _in->Drc - _in->data[r-1][c])/_dx;
                delta_y2 = (UF_OUTORMV(_dem,r+1,c)? 0.0: _in->data[r+1][c] - _in->Drc)/_dx;
                dy   = 0.5*UF_MinMod(delta_y1, delta_y2);

                outy1->Drc = _in->Drc-0.5 * dy;
                outy2->Drc = _in->Drc+0.5 * dy;
           }
        }
    }
}

void TWorld::UF2D_MUSCLE_operate(cTMap * _dem,cTMap * dt,int in_map1, int in_map2,int operation, int target)
{
    cTMap * outx1 =0;
    cTMap * outy1 =0;
    cTMap * outx2 =0;
    cTMap * outy2 =0;


    if(target == UF_MUSCLE_TARGET_OUT)
    {
        outx1 = UF2D_MUSCLE_OUT_x1;
        outx2 = UF2D_MUSCLE_OUT_x2;
        outy1 = UF2D_MUSCLE_OUT_y1;
        outy2 = UF2D_MUSCLE_OUT_y2;

    }else if(target == UF_MUSCLE_TARGET_IN1)
    {
        outx1 = UF2D_MUSCLE_1_x1;
        outx2 = UF2D_MUSCLE_1_x2;
        outy1 = UF2D_MUSCLE_1_y1;
        outy2 = UF2D_MUSCLE_1_y2;

    }else if(target == UF_MUSCLE_TARGET_IN2)
    {
        outx1 = UF2D_MUSCLE_2_x1;
        outx2 = UF2D_MUSCLE_2_x2;
        outy1 = UF2D_MUSCLE_2_y1;
        outy2 = UF2D_MUSCLE_2_y2;
    }else if(target == UF_MUSCLE_TARGET_IN3)
    {
        outx1 = UF2D_MUSCLE_3_x1;
        outx2 = UF2D_MUSCLE_3_x2;
        outy1 = UF2D_MUSCLE_3_y1;
        outy2 = UF2D_MUSCLE_3_y2;
    }

    UF2D_MUSCLE_operate(_dem,dt,in_map1,in_map2,operation,outx1,outx2,outy1,outy2);


}

void TWorld::UF2D_MUSCLE_operate(cTMap * _dem,cTMap * dt,int in_map1, int in_map2x,int in_map2y,int operation, int target)
{
    cTMap * outx1 =0;
    cTMap * outy1 =0;
    cTMap * outx2 =0;
    cTMap * outy2 =0;


    if(target == UF_MUSCLE_TARGET_OUT)
    {
        outx1 = UF2D_MUSCLE_OUT_x1;
        outx2 = UF2D_MUSCLE_OUT_x2;
        outy1 = UF2D_MUSCLE_OUT_y1;
        outy2 = UF2D_MUSCLE_OUT_y2;

    }else if(target == UF_MUSCLE_TARGET_IN1)
    {
        outx1 = UF2D_MUSCLE_1_x1;
        outx2 = UF2D_MUSCLE_1_x2;
        outy1 = UF2D_MUSCLE_1_y1;
        outy2 = UF2D_MUSCLE_1_y2;

    }else if(target == UF_MUSCLE_TARGET_IN2)
    {
        outx1 = UF2D_MUSCLE_2_x1;
        outx2 = UF2D_MUSCLE_2_x2;
        outy1 = UF2D_MUSCLE_2_y1;
        outy2 = UF2D_MUSCLE_2_y2;
    }else if(target == UF_MUSCLE_TARGET_IN3)
    {
        outx1 = UF2D_MUSCLE_3_x1;
        outx2 = UF2D_MUSCLE_3_x2;
        outy1 = UF2D_MUSCLE_3_y1;
        outy2 = UF2D_MUSCLE_3_y2;
    }

    UF2D_MUSCLE_operate(_dem,dt,in_map1,in_map2x,in_map2y,operation,outx1,outx2,outy1,outy2);


}

void TWorld::UF2D_MUSCLE_operate(cTMap * _dem,cTMap * dt,int in_map,int in_map2, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2)
{
    UF2D_MUSCLE_operate(_dem, dt,in_map, in_map2,in_map2, operation, outx1, outx2, outy1, outy2);

}

void TWorld::UF2D_MUSCLE_operate(cTMap * _dem,cTMap * dt,int in_map1, int in_map2x,int in_map2y, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2)
{

    cTMap * in1x1 =0;
    cTMap * in1y1 =0;
    cTMap * in1x2 =0;
    cTMap * in1y2 =0;

    cTMap * in2x1x =0;
    cTMap * in2y1x =0;
    cTMap * in2x2x =0;
    cTMap * in2y2x =0;

    cTMap * in2x1y =0;
    cTMap * in2y1y =0;
    cTMap * in2x2y =0;
    cTMap * in2y2y =0;

    if(in_map1 == UF_MUSCLE_TARGET_OUT)
    {
        in1x1 = UF2D_MUSCLE_OUT_x1;
        in1x2 = UF2D_MUSCLE_OUT_x2;
        in1y1 = UF2D_MUSCLE_OUT_y1;
        in1y2 = UF2D_MUSCLE_OUT_y2;

    }else if(in_map1 == UF_MUSCLE_TARGET_IN1)
    {
        in1x1 = UF2D_MUSCLE_1_x1;
        in1x2 = UF2D_MUSCLE_1_x2;
        in1y1 = UF2D_MUSCLE_1_y1;
        in1y2 = UF2D_MUSCLE_1_y2;

    }else if(in_map1 == UF_MUSCLE_TARGET_IN2)
    {
        in1x1 = UF2D_MUSCLE_2_x1;
        in1x2 = UF2D_MUSCLE_2_x2;
        in1y1 = UF2D_MUSCLE_2_y1;
        in1y2 = UF2D_MUSCLE_2_y2;
    }else if(in_map1 == UF_MUSCLE_TARGET_IN3)
    {
        in1x1 = UF2D_MUSCLE_3_x1;
        in1x2 = UF2D_MUSCLE_3_x2;
        in1y1 = UF2D_MUSCLE_3_y1;
        in1y2 = UF2D_MUSCLE_3_y2;
    }

    if(in_map2x == UF_MUSCLE_TARGET_OUT)
    {
        in2x1x = UF2D_MUSCLE_OUT_x1;
        in2x2x = UF2D_MUSCLE_OUT_x2;
        in2y1x = UF2D_MUSCLE_OUT_y1;
        in2y2x = UF2D_MUSCLE_OUT_y2;

    }else if(in_map2x == UF_MUSCLE_TARGET_IN1)
    {
        in2x1x = UF2D_MUSCLE_1_x1;
        in2x2x = UF2D_MUSCLE_1_x2;
        in2y1x = UF2D_MUSCLE_1_y1;
        in2y2x = UF2D_MUSCLE_1_y2;

    }else if(in_map2x == UF_MUSCLE_TARGET_IN2)
    {
        in2x1x = UF2D_MUSCLE_2_x1;
        in2x2x = UF2D_MUSCLE_2_x2;
        in2y1x = UF2D_MUSCLE_2_y1;
        in2y2x = UF2D_MUSCLE_2_y2;
    }else if(in_map2x == UF_MUSCLE_TARGET_IN3)
    {
        in2x1x = UF2D_MUSCLE_3_x1;
        in2x2x = UF2D_MUSCLE_3_x2;
        in2y1x = UF2D_MUSCLE_3_y1;
        in2y2x = UF2D_MUSCLE_3_y2;
    }

    if(in_map2y == UF_MUSCLE_TARGET_OUT)
    {
        in2x1y = UF2D_MUSCLE_OUT_x1;
        in2x2y = UF2D_MUSCLE_OUT_x2;
        in2y1y = UF2D_MUSCLE_OUT_y1;
        in2y2y = UF2D_MUSCLE_OUT_y2;

    }else if(in_map2y == UF_MUSCLE_TARGET_IN1)
    {
        in2x1y = UF2D_MUSCLE_1_x1;
        in2x2y = UF2D_MUSCLE_1_x2;
        in2y1y = UF2D_MUSCLE_1_y1;
        in2y2y = UF2D_MUSCLE_1_y2;

    }else if(in_map2y == UF_MUSCLE_TARGET_IN2)
    {
        in2x1y = UF2D_MUSCLE_2_x1;
        in2x2y = UF2D_MUSCLE_2_x2;
        in2y1y = UF2D_MUSCLE_2_y1;
        in2y2y = UF2D_MUSCLE_2_y2;
    }else if(in_map2y == UF_MUSCLE_TARGET_IN3)
    {
        in2x1y = UF2D_MUSCLE_3_x1;
        in2x2y = UF2D_MUSCLE_3_x2;
        in2y1y = UF2D_MUSCLE_3_y1;
        in2y2y = UF2D_MUSCLE_3_y2;
    }


    if(dt != 0)
    {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc*in2x1x->Drc;
                    outx2->Drc = in1x2->Drc*in2x2x->Drc;

                    outy1->Drc = in1y1->Drc*in2y1y->Drc;
                    outy2->Drc = in1y2->Drc*in2y2y->Drc;

                }}}
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc/in2x1x->Drc;
                    outx2->Drc = in1x2->Drc/in2x2x->Drc;

                    outy1->Drc = in1y1->Drc/in2y1y->Drc;
                    outy2->Drc = in1y2->Drc/in2y2y->Drc;
                }}}
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc+in2x1x->Drc;
                    outx2->Drc = in1x2->Drc+in2x2x->Drc;

                    outy1->Drc = in1y1->Drc+in2y1y->Drc;
                    outy2->Drc = in1y2->Drc+in2y2y->Drc;
                }}}
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc-in2x1x->Drc;
                    outx2->Drc = in1x2->Drc-in2x2x->Drc;

                    outy1->Drc = in1y1->Drc-in2y1y->Drc;
                    outy2->Drc = in1y2->Drc-in2y2y->Drc;
                }}}
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = pow(in1x1->Drc,in2x1x->Drc);
                    outx2->Drc = pow(in1x2->Drc,in2x2x->Drc);

                    outy1->Drc = pow(in1y1->Drc,in2y1y->Drc);
                    outy2->Drc = pow(in1y2->Drc,in2y2y->Drc);
                }}}
                break;
        }

   }else
   {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc*in2x1x->Drc;
                    outx2->Drc = in1x2->Drc*in2x2x->Drc;

                    outy1->Drc = in1y1->Drc*in2y1y->Drc;
                    outy2->Drc = in1y2->Drc*in2y2y->Drc;
                }
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc/in2x1x->Drc;
                    outx2->Drc = in1x2->Drc/in2x2x->Drc;

                    outy1->Drc = in1y1->Drc/in2y1y->Drc;
                    outy2->Drc = in1y2->Drc/in2y2y->Drc;
                }
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc+in2x1x->Drc;
                    outx2->Drc = in1x2->Drc+in2x2x->Drc;

                    outy1->Drc = in1y1->Drc+in2y1y->Drc;
                    outy2->Drc = in1y2->Drc+in2y2y->Drc;
                }
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc-in2x1x->Drc;
                    outx2->Drc = in1x2->Drc-in2x2x->Drc;

                    outy1->Drc = in1y1->Drc-in2y1y->Drc;
                    outy2->Drc = in1y2->Drc-in2y2y->Drc;
                }
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = pow(in1x1->Drc,in2x1x->Drc);
                    outx2->Drc = pow(in1x2->Drc,in2x2x->Drc);

                    outy1->Drc = pow(in1y1->Drc,in2y1y->Drc);
                    outy2->Drc = pow(in1y2->Drc,in2y2y->Drc);
                }
                break;
        }
   }



}


void TWorld::UF2D_MUSCLE_operate(cTMap * _dem,cTMap * dt,int in_map,double in_double, int operation, int target)
{

    cTMap * outx1 =0;
    cTMap * outy1 =0;
    cTMap * outx2 =0;
    cTMap * outy2 =0;


    if(target == UF_MUSCLE_TARGET_OUT)
    {
        outx1 = UF2D_MUSCLE_OUT_x1;
        outx2 = UF2D_MUSCLE_OUT_x2;
        outy1 = UF2D_MUSCLE_OUT_y1;
        outy2 = UF2D_MUSCLE_OUT_y2;

    }else if(target == UF_MUSCLE_TARGET_IN1)
    {
        outx1 = UF2D_MUSCLE_1_x1;
        outx2 = UF2D_MUSCLE_1_x2;
        outy1 = UF2D_MUSCLE_1_y1;
        outy2 = UF2D_MUSCLE_1_y2;

    }else if(target == UF_MUSCLE_TARGET_IN2)
    {
        outx1 = UF2D_MUSCLE_2_x1;
        outx2 = UF2D_MUSCLE_2_x2;
        outy1 = UF2D_MUSCLE_2_y1;
        outy2 = UF2D_MUSCLE_2_y2;
    }else if(target == UF_MUSCLE_TARGET_IN3)
    {
        outx1 = UF2D_MUSCLE_3_x1;
        outx2 = UF2D_MUSCLE_3_x2;
        outy1 = UF2D_MUSCLE_3_y1;
        outy2 = UF2D_MUSCLE_3_y2;
    }

    UF2D_MUSCLE_operate(_dem,dt,in_map,in_double,operation,outx1,outx2,outy1,outy2);


}
void TWorld::UF2D_MUSCLE_operate(cTMap * _dem,cTMap * dt,int in_map,double in_double, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2)
{
    cTMap * in1x1 =0;
    cTMap * in1y1 =0;
    cTMap * in1x2 =0;
    cTMap * in1y2 =0;


    if(in_map== UF_MUSCLE_TARGET_OUT)
    {
        in1x1 = UF2D_MUSCLE_OUT_x1;
        in1x2 = UF2D_MUSCLE_OUT_x2;
        in1y1 = UF2D_MUSCLE_OUT_y1;
        in1y2 = UF2D_MUSCLE_OUT_y2;

    }else if(in_map == UF_MUSCLE_TARGET_IN1)
    {
        in1x1 = UF2D_MUSCLE_1_x1;
        in1x2 = UF2D_MUSCLE_1_x2;
        in1y1 = UF2D_MUSCLE_1_y1;
        in1y2 = UF2D_MUSCLE_1_y2;

    }else if(in_map == UF_MUSCLE_TARGET_IN2)
    {
        in1x1 = UF2D_MUSCLE_2_x1;
        in1x2 = UF2D_MUSCLE_2_x2;
        in1y1 = UF2D_MUSCLE_2_y1;
        in1y2 = UF2D_MUSCLE_2_y2;
    }else if(in_map == UF_MUSCLE_TARGET_IN3)
    {
        in1x1 = UF2D_MUSCLE_3_x1;
        in1x2 = UF2D_MUSCLE_3_x2;
        in1y1 = UF2D_MUSCLE_3_y1;
        in1y2 = UF2D_MUSCLE_3_y2;
    }

    if(dt != 0)
    {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc*in_double;
                    outx2->Drc = in1x2->Drc*in_double;

                    outy1->Drc = in1y1->Drc*in_double;
                    outy2->Drc = in1y2->Drc*in_double;
                }}}
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc/in_double;
                    outx2->Drc = in1x2->Drc/in_double;

                    outy1->Drc = in1y1->Drc/in_double;
                    outy2->Drc = in1y2->Drc/in_double;
                }}}
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc+in_double;
                    outx2->Drc = in1x2->Drc+in_double;

                    outy1->Drc = in1y1->Drc+in_double;
                    outy2->Drc = in1y2->Drc+in_double;
                }}}
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = in1x1->Drc-in_double;
                    outx2->Drc = in1x2->Drc-in_double;

                    outy1->Drc = in1y1->Drc-in_double;
                    outy2->Drc = in1y2->Drc-in_double;
                }}}
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF2D_DT
                {
                    outx1->Drc = pow(in1x1->Drc,in_double);
                    outx2->Drc = pow(in1x2->Drc,in_double);

                    outy1->Drc = pow(in1y1->Drc,in_double);
                    outy2->Drc = pow(in1y2->Drc,in_double);
                }}}
                break;
        }

   }else
   {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc*in_double;
                    outx2->Drc = in1x2->Drc*in_double;

                    outy1->Drc = in1y1->Drc*in_double;
                    outy2->Drc = in1y2->Drc*in_double;
                }
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc/in_double;
                    outx2->Drc = in1x2->Drc/in_double;

                    outy1->Drc = in1y1->Drc/in_double;
                    outy2->Drc = in1y2->Drc/in_double;
                }
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc+in_double;
                    outx2->Drc = in1x2->Drc+in_double;

                    outy1->Drc = in1y1->Drc+in_double;
                    outy2->Drc = in1y2->Drc+in_double;
                }
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = in1x1->Drc-in_double;
                    outx2->Drc = in1x2->Drc-in_double;

                    outy1->Drc = in1y1->Drc-in_double;
                    outy2->Drc = in1y2->Drc-in_double;
                }
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF2D
                {
                    outx1->Drc = pow(in1x1->Drc,in_double);
                    outx2->Drc = pow(in1x2->Drc,in_double);

                    outy1->Drc = pow(in1y1->Drc,in_double);
                    outy2->Drc = pow(in1y2->Drc,in_double);
                }
                break;
        }
   }
}







//////////////////////////////////////////////////////////////////////////////////





void TWorld::UF1D_MUSCLE(cTMap * _ldd,cTMap * _lddw,cTMap * dt,cTMap * _in, int target)
{
    double dx;

    cTMap * outx1 =0;
    cTMap * outx2 =0;

    if(target == UF_MUSCLE_TARGET_OUT)
    {
        outx1 = UF1D_MUSCLE_OUT_x1;
        outx2 = UF1D_MUSCLE_OUT_x2;

    }else if(target == UF_MUSCLE_TARGET_IN1)
    {
        outx1 = UF1D_MUSCLE_1_x1;
        outx2 = UF1D_MUSCLE_1_x2;

    }else if(target == UF_MUSCLE_TARGET_IN2)
    {
        outx1 = UF1D_MUSCLE_2_x1;
        outx2 = UF1D_MUSCLE_2_x2;
    }else if(target == UF_MUSCLE_TARGET_IN3)
    {
        outx1 = UF1D_MUSCLE_3_x1;
        outx2 = UF1D_MUSCLE_3_x2;
    }

    if(dt != 0)
    {
        FOR_ROW_COL_UF1D_DT
        {
            dx = UF1D_Derivative(_ldd,_lddw,_in,r,c,true);

            outx1->Drc = _in->Drc+dx;
            outx2->Drc = _in->Drc-dx;
        }
   }else
   {
        FOR_ROW_COL_UF1D
        {
            dx = UF1D_Derivative(_ldd,_lddw,_in,r,c,true);

            outx1->Drc = _in->Drc+dx;
            outx2->Drc = _in->Drc-dx;
        }
   }

}


void TWorld::UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap * dt,int in_map1, int in_map2,int operation, int target)
{
    cTMap * outx1 =0;
    cTMap * outx2 =0;

    if(target == UF_MUSCLE_TARGET_OUT)
    {
        outx1 = UF1D_MUSCLE_OUT_x1;
        outx2 = UF1D_MUSCLE_OUT_x2;

    }else if(target == UF_MUSCLE_TARGET_IN1)
    {
        outx1 = UF1D_MUSCLE_1_x1;
        outx2 = UF1D_MUSCLE_1_x2;

    }else if(target == UF_MUSCLE_TARGET_IN2)
    {
        outx1 = UF1D_MUSCLE_2_x1;
        outx2 = UF1D_MUSCLE_2_x2;
    }else if(target == UF_MUSCLE_TARGET_IN3)
    {
        outx1 = UF1D_MUSCLE_3_x1;
        outx2 = UF1D_MUSCLE_3_x2;
    }

    UF1D_MUSCLE_operate(_ldd,_lddw, dt,in_map1,in_map2, operation, outx1, outx2);



}
void TWorld::UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap * dt,int in_map1,int in_map2, int operation, cTMap * outx1, cTMap * outx2)
{
    cTMap * in1x1 =0;
    cTMap * in1x2 =0;

    cTMap * in2x1 =0;
    cTMap * in2x2 =0;


    if(in_map1 == UF_MUSCLE_TARGET_OUT)
    {
        in1x1 = UF1D_MUSCLE_OUT_x1;
        in1x2 = UF1D_MUSCLE_OUT_x2;

    }else if(in_map1 == UF_MUSCLE_TARGET_IN1)
    {
        in1x1 = UF1D_MUSCLE_1_x1;
        in1x2 = UF1D_MUSCLE_1_x2;

    }else if(in_map1 == UF_MUSCLE_TARGET_IN2)
    {
        in1x1 = UF1D_MUSCLE_2_x1;
        in1x2 = UF1D_MUSCLE_2_x2;
    }else if(in_map1 == UF_MUSCLE_TARGET_IN3)
    {
        in1x1 = UF1D_MUSCLE_3_x1;
        in1x2 = UF1D_MUSCLE_3_x2;
    }

    if(in_map2 == UF_MUSCLE_TARGET_OUT)
    {
        in2x1 = UF1D_MUSCLE_OUT_x1;
        in2x2 = UF1D_MUSCLE_OUT_x2;

    }else if(in_map2 == UF_MUSCLE_TARGET_IN1)
    {
        in2x1 = UF1D_MUSCLE_1_x1;
        in2x2 = UF1D_MUSCLE_1_x2;

    }else if(in_map2 == UF_MUSCLE_TARGET_IN2)
    {
        in2x1 = UF1D_MUSCLE_2_x1;
        in2x2 = UF1D_MUSCLE_2_x2;
    }else if(in_map2 == UF_MUSCLE_TARGET_IN3)
    {
        in2x1 = UF1D_MUSCLE_3_x1;
        in2x2 = UF1D_MUSCLE_3_x2;
    }

    if(dt != 0)
    {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc*in2x1->Drc;
                    outx2->Drc = in1x2->Drc*in2x2->Drc;
                }
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc/in2x1->Drc;
                    outx2->Drc = in1x2->Drc/in2x2->Drc;
                }
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc+in2x1->Drc;
                    outx2->Drc = in1x2->Drc+in2x2->Drc;
                }
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc-in2x1->Drc;
                    outx2->Drc = in1x2->Drc-in2x2->Drc;
                }
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = pow(in1x1->Drc,in2x1->Drc);
                    outx2->Drc = pow(in1x2->Drc,in2x2->Drc);
                }
                break;
        }

   }else
   {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc*in2x1->Drc;
                    outx2->Drc = in1x2->Drc*in2x2->Drc;
                }
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc/in2x1->Drc;
                    outx2->Drc = in1x2->Drc/in2x2->Drc;
                }
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc+in2x1->Drc;
                    outx2->Drc = in1x2->Drc+in2x2->Drc;
                }
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc-in2x1->Drc;
                    outx2->Drc = in1x2->Drc-in2x2->Drc;
                }
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = pow(in1x1->Drc,in2x1->Drc);
                    outx2->Drc = pow(in1x2->Drc,in2x2->Drc);
                }
                break;
        }
   }

}

void TWorld::UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap * dt,int in_map,double in_double, int operation, int target)
{
    cTMap * outx1 =0;
    cTMap * outx2 =0;


    if(target == UF_MUSCLE_TARGET_OUT)
    {
        outx1 = UF1D_MUSCLE_OUT_x1;
        outx2 = UF1D_MUSCLE_OUT_x2;

    }else if(target == UF_MUSCLE_TARGET_IN1)
    {
        outx1 = UF1D_MUSCLE_1_x1;
        outx2 = UF1D_MUSCLE_1_x2;

    }else if(target == UF_MUSCLE_TARGET_IN2)
    {
        outx1 = UF1D_MUSCLE_2_x1;
        outx2 = UF1D_MUSCLE_2_x2;
    }else if(target == UF_MUSCLE_TARGET_IN3)
    {
        outx1 = UF1D_MUSCLE_3_x1;
        outx2 = UF1D_MUSCLE_3_x2;
    }

    UF1D_MUSCLE_operate(_ldd,_lddw,dt,in_map,in_double,operation,outx1,outx2);
}

void TWorld::UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap * dt,int in_map,double in_double, int operation, cTMap * outx1, cTMap * outx2)
{

    cTMap * in1x1 =0;
    cTMap * in1x2 =0;



    if(in_map== UF_MUSCLE_TARGET_OUT)
    {
        in1x1 = UF1D_MUSCLE_OUT_x1;
        in1x2 = UF1D_MUSCLE_OUT_x2;

    }else if(in_map == UF_MUSCLE_TARGET_IN1)
    {
        in1x1 = UF1D_MUSCLE_1_x1;
        in1x2 = UF1D_MUSCLE_1_x2;

    }else if(in_map == UF_MUSCLE_TARGET_IN2)
    {
        in1x1 = UF1D_MUSCLE_2_x1;
        in1x2 = UF1D_MUSCLE_2_x2;
    }else if(in_map == UF_MUSCLE_TARGET_IN3)
    {
        in1x1 = UF1D_MUSCLE_3_x1;
        in1x2 = UF1D_MUSCLE_3_x2;
    }


    if(dt != 0)
    {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc*in_double;
                    outx2->Drc = in1x2->Drc*in_double;

                }
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc/in_double;
                    outx2->Drc = in1x2->Drc/in_double;

                }
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc+in_double;
                    outx2->Drc = in1x2->Drc+in_double;

                }
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = in1x1->Drc-in_double;
                    outx2->Drc = in1x2->Drc-in_double;

                }
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF1D_DT
                {
                    outx1->Drc = pow(in1x1->Drc,in_double);
                    outx2->Drc = pow(in1x2->Drc,in_double);

                }
                break;
        }

   }else
   {
        switch(operation)
        {
            case UF_MUSCLE_mult:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc*in_double;
                    outx2->Drc = in1x2->Drc*in_double;

                }
                break;
            case UF_MUSCLE_div:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc/in_double;
                    outx2->Drc = in1x2->Drc/in_double;

                }
                break;
            case UF_MUSCLE_add:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc+in_double;
                    outx2->Drc = in1x2->Drc+in_double;

                }
                break;
            case UF_MUSCLE_sub:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = in1x1->Drc-in_double;
                    outx2->Drc = in1x2->Drc-in_double;

                }
                break;
            case UF_MUSCLE_pow:
                FOR_ROW_COL_UF1D
                {
                    outx1->Drc = pow(in1x1->Drc,in_double);
                    outx2->Drc = pow(in1x2->Drc,in_double);

                }
                break;
        }
   }
}

