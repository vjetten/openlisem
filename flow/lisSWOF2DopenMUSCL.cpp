/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

// This particular code uses much of the principles of the FullSWOF model
// https://arxiv.org/abs/1204.3210
// https://sourcesup.renater.fr/frs/?group_id=895&release_id=3901#fullswof_2d-_1.10.00-title-content
// the scheme is made suited for parallel processing
// LICENCE: http://cecill.info/licences/Licence_CeCILL_V2-en.html

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

//#define LIMIT(V,L) (V < 0.0 ? -1.0 : 1.0)*qMin(L,fabs(V))
//#define SIGN(V)(V < 0 ? -1.0 : 1.0)

//----------------------------------------------------------------------------------------
double TWorld::fullSWOF2openMUSCL(cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    double timesum = 0;
    double dt_max = qMin(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;
    sumh = getMass(h);

    //F_MaxIter = 10000;
    Fill(*tmd,0);
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (h->Drc > F_minWH)
            tmd->Drc = 1; // flag which cells have to be calculated
    }}

    do {

        //if (SwitchErosion)
        //sumS = getMassSed(SSFlood, 0);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            FloodDT->Drc = dt_max;
            //activeCells->Drc = 0;
            tma->Drc = h->Drc;
            tmb->Drc = u->Drc;
            tmc->Drc = v->Drc;
            // save the values at the start of the run for MUSCL
        }}

        dt_req_min = doSWOFMUSCLdt(dt_max, timesum, h, u, v, z);
        // do MUSCL (optional), Riemann etc, get back smallest dt
        // in the original code this is split in reconstruction/MUSCL and maincalcflux

        doSWOFStV(dt_req_min, h, u, v);
        // Saint-Venant calculations for new h, u, v
        // called maincalcscheme in fullSWOF

        //until here is first order ! just one calculation

        // 2nd order, with avg according to Heun, according to fullswof hean should allways be done!
        int step = 0;
        double dt1;

        if (SwitchMUSCL) {
            do {
                step++;
                dt1 = dt_req_min;

                dt_req_min = doSWOFMUSCLdt(dt1, timesum, h, u, v, z);

            } while (dt1 > dt_req_min && step < 5);

            doSWOFStV(dt_req_min, h, u, v);

            //Heun, see SWOF doc
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                double havg = 0.5*(tma->Drc + h->Drc); // avg original before loops and second estimation
                if (havg >= he_ca){
                    double q1 = 0.5*(tma->Drc*tmb->Drc + h->Drc*u->Drc);
                    u->Drc = q1/havg;
                    double q2 = 0.5*(tma->Drc*tmc->Drc + h->Drc*v->Drc);
                    v->Drc = q2/havg;
                    h->Drc = havg;
                } else {
                    h->Drc = 0.0;
                    u->Drc = 0.0;
                    v->Drc = 0.0;
                }
            }}
        } // MUSCL

        if (SwitchErosion && !SwitchErosionOutsideLoop) {
            SWOFSediment(dt_req_min, h, FlowWidth, u,v);
        }

        if (Switch2DDiagonalFlow) {
            SWOFDiagonalFlowNew(dt_req_min, h, u, v);
        }

        timesum += dt_req_min;
        count++; // nr loops

        stop = timesum > _dt-0.001;
        if(count > F_MaxIter)
        stop = true;

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            tmd->Drc = 0;
            if (h->Drc > F_minWH && qSqrt(u->Drc*u->Drc+v->Drc*v->Drc) > F_minWH)
                tmd->Drc = 1;
        }}

    } while (!stop);

    // small mass balance corrections within 2d flow
    correctMassBalance(sumh, h);

    if (SwitchErosion && SwitchErosionOutsideLoop) {
        SWOFSediment(_dt, h, FlowWidth, u,v);
    }

    if (FlowBoundaryType > 0) {
        Boundary2Ddyn(_dt, h, u, v);
    }

    //floodCount(h);

    iter_n = qMax(1,count);
    return(count > 0 ? _dt/count : _dt);

}
//------------------------------------------------------------------------------------------------------
double TWorld::doSWOFMUSCLdt(double dt, double timesum, cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    // boundary
    double factor = exp(-0.005*_dx); // sort of cell size dpendent, if large cells, farther away so more dip
    double factor2 = factor;//pow(factor,0.667); // manning reduction V=h^2/3

   // Fill(*tmd,0);
    // map edges are zero, avoid domain touching the edges
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        //if (h->Drc > he_ca)
        //     tmd->Drc = 1;

        if (c > 0 && !MV(r,c-1)        )  tmd->data[r][c-1] = 1;
        if (c < _nrCols-1 && !MV(r,c+1))  tmd->data[r][c+1] = 1;
        if (r > 0 && !MV(r-1,c)        )  tmd->data[r-1][c] = 1;
        if (r < _nrRows-1 && !MV(r+1,c))  tmd->data[r+1][c] = 1;

        if (r == 0 || r == _nrRows-1 || c == 0 || c == _nrCols-1)
            tmd->Drc = 0;
        if (DomainEdge->Drc > 0 && FlowBoundary->Drc == 0)
            tmd->Drc = 0;
    }}

    //do all flow and state calculations
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (tmd->Drc == 1) {
            double dx = _dx; // do not do channeladj because the channelflood function does this already
            double dy = _dx;
            double H, Z, U, V;
            bool bc1, bc2, br1, br2;
            double z_x1, z_x2, z_y1, z_y2;
            double h_x1, h_x2, h_y1, h_y2;
            double u_x1, u_x2, u_y1, u_y2;
            double v_x1, v_x2, v_y1, v_y2;
            double dz_x1, dz_x2, dz_y1, dz_y2;
            double fb_x1=0,fb_x2=0,fb_y1=0,fb_y2=0;
            double delzcx=0, delzcy=0;
            double hx1r, hxl, hxr, hx2l;  // c-1, c, c+1
            double hy1d, hyu, hyd, hy2u;  // r-1, r, r+1
            double ux1r, uxl, uxr, ux2l;
            double uy1d, uyu, uyd, uy2u;
            double vx1r, vxl, vxr, vx2l;
            double vy1d, vyu, vyd, vy2u;
            double h_x1r, h_xl, h_xr, h_x2l;
            double h_y1d, h_yu, h_yd, h_y2u;
            vec4 hll_x1;
            vec4 hll_x2;
            vec4 hll_y1;
            vec4 hll_y2;

            bc1 = c > 0 && !MV(r,c-1)        ;
            bc2 = c < _nrCols-1 && !MV(r,c+1);
            br1 = r > 0 && !MV(r-1,c)        ;
            br2 = r < _nrRows-1 && !MV(r+1,c);

            // get value for all 5 cells center, up, down, left, right
            //if MV the cell gets the center cell values
            Z = z->Drc;
            H = h->Drc;
            U = u->Drc;
            V = v->Drc;
            if (bc1) {
                z_x1 = z->data[r][c-1];
                h_x1 = h->data[r][c-1];
                u_x1 = u->data[r][c-1];
                v_x1 = v->data[r][c-1];
            } else {
                z_x1 = Z;
                h_x1 = H;
                u_x1 = U;
                v_x1 = V;
            }

            if (bc2) {
                z_x2 = z->data[r][c+1];
                h_x2 = h->data[r][c+1];
                u_x2 = u->data[r][c+1];
                v_x2 = v->data[r][c+1];
            } else {
                z_x2 = Z;
                h_x2 = H;
                u_x2 = U;
                v_x2 = V;
            }
            if (br1) {
                z_y1 = z->data[r-1][c];
                h_y1 = h->data[r-1][c];
                u_y1 = u->data[r-1][c];
                v_y1 = v->data[r-1][c];
            } else {
                z_y1 = Z;
                h_y1 = H;
                u_y1 = U;
                v_y1 = V;
            }
            if (br2) {
                z_y2 = z->data[r+1][c];
                h_y2 = h->data[r+1][c];
                u_y2 = u->data[r+1][c];
                v_y2 = v->data[r+1][c];
            } else {
                z_y2 = Z;
                h_y2 = H;
                u_y2 = U;
                v_y2 = V;
            }
            double Hc = factor*H;
            double Uc = factor2*U;
            double Vc = factor2*V;

            if (FlowBoundary->Drc > 0) {
                //if left does not exist and right exist estimate gradient
                if (c > 0 && MV(r,c-1) && !MV(r,c+1)) {
                    if (h_x2+z_x2 > H+Z) {
                        h_x1 = Hc;
                        u_x1 = Uc;
                        v_x1 = Vc;
                    }
                }
                if (c < _nrCols-1 && MV(r,c+1) && !MV(r,c-1)) {
                    if (h_x1+z_x1 > H+Z){
                        h_x2 = Hc;
                        u_x2 = Uc;
                        v_x2 = Vc;
                    }
                }
                if (r > 0 && MV(r-1,c) && !MV(r+1,c)) {
                    if (h_y2+z_y2 > H+Z) {
                        h_y1 = Hc;
                        u_y1 = Uc;
                        v_y1 = Vc;
                    }
                }
                if (r < _nrRows-1 && MV(r+1,c) && !MV(r-1,c)) {
                    if (h_y1+z_y1 > H+Z) {
                        h_y2 = Hc;
                        u_y2 = Uc;
                        v_y2 = Vc;
                    }
                }
            }

            dz_x1 = (Z - z_x1);
            dz_x2 = (z_x2 - Z);
            dz_y1 = (Z - z_y1);
            dz_y2 = (z_y2 - Z);

            if (SwitchFlowBarriers) {
                fb_x1 = bc1 ? qMax(FlowBarrierW->Drc, FlowBarrierE->data[r][c-1]) : FlowBarrierW->Drc;
                fb_x2 = bc2 ? qMax(FlowBarrierE->Drc, FlowBarrierE->data[r][c+1]) : FlowBarrierE->Drc;
                fb_y1 = br1 ? qMax(FlowBarrierN->Drc, FlowBarrierS->data[r-1][c]) : FlowBarrierN->Drc;
                fb_y2 = br2 ? qMax(FlowBarrierS->Drc, FlowBarrierN->data[r+1][c]) : FlowBarrierS->Drc;
            }

            // non-muscl solution, cell centres for boundaries in x and y directions
            hx1r = h_x1; hxl = H; hxr = H; hx2l = h_x2;  // x-1 r; x l; x r; x+1 l
            ux1r = u_x1; uxl = U; uxr = U; ux2l = u_x2;
            vx1r = v_x1; vxl = V; vxr = V; vx2l = v_x2;

            hy1d = h_y1; hyu = H; hyd = H; hy2u = h_y2;
            uy1d = u_y1; uyu = U; uyd = U; uy2u = u_y2;
            vy1d = v_y1; vyu = V; vyd = V; vy2u = v_y2;

            //======== MUSCL: on the 4 boundaties of a gridcell interpolate from the center values
            // called "reconstruction" in SWOF code
            if (SwitchMUSCL) {
                bool b2c1 ,b2c2 ,b2r1 ,b2r2;
                double h_xx1, h_xx2, u_xx1, u_xx2, v_xx1, v_xx2;
                double h_yy1, h_yy2, u_yy1, u_yy2, v_yy1, v_yy2;
                double dh, du, dv, dz_h;
                double delta_h1, delta_h2, delta_h3, delta_h4;
                double delta_u1, delta_u2, delta_u3, delta_u4;
                double delta_v1, delta_v2, delta_v3, delta_v4;

                b2c1 = c > 1 && !MV(r,c-2)         ;
                b2c2 = c < _nrCols-2 && !MV(r,c+2);
                b2r1 = r > 1 && !MV(r-2,c)        ;
                b2r2 = r < _nrRows-2 && !MV(r+2,c);

                if(b2c1) {
                    h_xx1 = h->data[r][c-2];
                    u_xx1 = u->data[r][c-2];
                    v_xx1 = v->data[r][c-2];
                } else {
                    h_xx1 = h_x1;
                    u_xx1 = u_x1;
                    v_xx1 = v_x1;
                }
                if(b2c2) {
                    h_xx2 = h->data[r][c+2];
                    u_xx2 = u->data[r][c+2];
                    v_xx2 = v->data[r][c+2];
                } else {
                    h_xx2 = h_x2;
                    u_xx2 = u_x2;
                    v_xx2 = v_x2;
                }
                if(b2r1) {
                    h_yy1 = h->data[r-2][c];
                    u_yy1 = u->data[r-2][c];
                    v_yy1 = v->data[r-2][c];
                } else {
                    h_yy1 = h_y1;
                    u_yy1 = u_y1;
                    v_yy1 = v_y1;
                }
                if(b2r2) {
                    h_yy2 = h->data[r+2][c];
                    u_yy2 = u->data[r+2][c];
                    v_yy2 = v->data[r+2][c];
                }else {
                    h_yy2 = h_y2;
                    u_yy2 = u_y2;
                    v_yy2 = v_y2;
                }

                //horizontal direction, leftn to right
                // x-1-x-2   x-x-1  x+1-x   x+2-x+1        always right minus left
                delta_h1 = h_x1 - h_xx1; delta_h2 = H-h_x1; delta_h3 = h_x2-H; delta_h4 = h_xx2-h_x2;
                delta_u1 = u_x1 - u_xx1; delta_u2 = U-u_x1; delta_u3 = u_x2-U; delta_u4 = u_xx2-u_x2;
                delta_v1 = v_x1 - v_xx1; delta_v2 = V-v_x1; delta_v3 = v_x2-V; delta_v4 = v_xx2-v_x2;

                // center cell, all boundaries
                dh = limiter(delta_h2, delta_h3);
                du = limiter(delta_u2, delta_u3);
                dv = limiter(delta_v2, delta_v3);
                hxl = H - 0.5*dh;
                hxr = H + 0.5*dh;
                if (H > he_ca) {
                    uxl = U - 0.5*du*hxl/H;
                    uxr = U + 0.5*du*hxr/H;
                    vxl = V - 0.5*dv*hxl/H;
                    vxr = V + 0.5*dv*hxr/H;
                } else {
                    uxl = U - 0.5*du;
                    uxr = U + 0.5*du;
                    vxl = V - 0.5*dv;
                    vxr = V + 0.5*dv;
                }

                dz_h = limiter(delta_h2 + (Z-dz_x1), delta_h3 + (dz_x2-Z));
                delzcx = (Z+0.5*(dz_h-dh))-(Z+0.5*(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!2*(dz_h - dh); //

                // left hand cell, right boundary
                dh = limiter(delta_h1, delta_h2);
                du = limiter(delta_u1, delta_u2);
                dv = limiter(delta_v1, delta_v2);
                hx1r = h_x1 + 0.5*dh;
                if (H > he_ca) {
                    ux1r = u_x1 + 0.5*du*hxl/H;
                    vx1r = v_x1 + 0.5*dv*hxl/H;
                } else {
                    ux1r = u_x1 + 0.5*du;
                    vx1r = v_x1 + 0.5*dv;
                }

                // right hand cell, left boundary
                dh = limiter(delta_h3, delta_h4);
                du = limiter(delta_u3, delta_u4);
                dv = limiter(delta_v3, delta_v4);
                hx2l = h_x2 - 0.5*dh;
                if (H > he_ca) {
                    ux2l = u_x2 - 0.5*du*hxr/H;
                    vx2l = v_x2 - 0.5*dv*hxr/H;
                } else {
                    ux2l = u_x2 - 0.5*du;
                    vx2l = v_x2 - 0.5*dv;
                }

                // vertical, direction from up to down
                // y-1 - y-2   y-y-1  y+1-y   y+2-y+1        always down minus up
                delta_h1 = h_y1 - h_yy1; delta_h2 = H-h_y1; delta_h3 = h_y2-H; delta_h4 = h_yy2-h_y2;
                delta_u1 = u_y1 - u_yy1; delta_u2 = U-u_y1; delta_u3 = u_y2-U; delta_u4 = u_yy2-u_y2;
                delta_v1 = v_y1 - v_yy1; delta_v2 = V-v_y1; delta_v3 = v_y2-V; delta_v4 = v_yy2-v_y2;

                // center cell, all boundaries
                dh = limiter(delta_h2, delta_h3);
                du = limiter(delta_u2, delta_u3);
                dv = limiter(delta_v2, delta_v3);
                hyu = H - 0.5*dh;
                hyd = H + 0.5*dh;
                if (H > he_ca) {
                    uyu = U - 0.5*du*hyu/H;
                    uyd = U + 0.5*du*hyd/H;
                    vyu = V - 0.5*dv*hyu/H;
                    vyd = V + 0.5*dv*hyd/H;
                } else {
                    uyu = U - 0.5*du;
                    uyd = U + 0.5*du;
                    vyu = V - 0.5*dv;
                    vyd = V + 0.5*dv;
                }

                dz_h = limiter(delta_h1 + (Z-dz_y1), delta_h2 + (dz_y2-Z));
                delzcy = (Z+0.5*(dz_h-dh))-(Z+0.5*(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!2*(dz_h-dh);

                // upper cell, lower boundary
                dh = limiter(delta_h1, delta_h2);
                du = limiter(delta_u1, delta_u2);
                dv = limiter(delta_v1, delta_v2);
                hy1d = h_y1 + 0.5*dh;
                if (H > he_ca) {
                    uy1d = u_y1 + 0.5*du*hyu/H;
                    vy1d = v_y1 + 0.5*dv*hyu/H;
                } else {
                    uy1d = u_y1 + 0.5*du;
                    vy1d = v_y1 + 0.5*dv;
                }

                // lower cell, up boundary
                dh = limiter(delta_h3, delta_h4);
                du = limiter(delta_u3, delta_u4);
                dv = limiter(delta_v3, delta_v4);
                hy2u = h_y2 - 0.5*dh;
                if (H > he_ca) {
                    uy2u = u_y2 - 0.5*du*hyd/H;
                    vy2u = v_y2 - 0.5*dv*hyd/H;
                } else {
                    uy2u = u_y2 - 0.5*du;
                    vy2u = v_y2 - 0.5*dv;
                }

            } //MUSCL

            //########### calculate Riemann valaues for all four boundaries of a cell ############

            // if muscl H and h_x1 etc become Hx1l and hx1r
            // z is blocking to prevent flow when water is flat and Z is not flat, described in article SWOF
            // barrier is ourown additiona, to vcreate flood walls.
            //result Riemann
            //  1st component [0]: Mass flux per meter ( dus (m3/s)/(m) = m2/s, unit discharge
            //  2e component [1]: Momentum flux direction of flow ( m4/s2)/(m) = m3/s2 = h*u*u)
            //  3d component [3]: Momentum flux perpendicular to flow ( (m4/s2)/(m) = m3/s2 = h*u*v)
            //  4th component[3]: celerity (time)


            //left and right hand side of c and c-1 (x and x1)
            if (bc1) {
                h_x1r = qMax(0.0, hx1r - qMax(0.0,  dz_x1 + fb_x1)); //rechts van c-1
                h_xl  = qMax(0.0, hxl  - qMax(0.0, -dz_x1 + fb_x1)); //links van het midden
                //fb1 is barrier height (m) between c and c-1 cell boundary
                // if h_x1r or h_xl < z+barrier then make it zero, no pressure on that boundary
                // dz_x1 = (Z - z_x1);
            } else {
                h_x1r = 0.0;
            }
            if (h_x1r == 0) {
                ux1r = 0;
                vx1r = 0;
            }
            if (h_xl == 0) {
                uxl = 0;
                vxl = 0;
            }
            hll_x1 = F_Riemann(h_x1r,ux1r,vx1r, h_xl,uxl,vxl); // c-1 (x1 right) and c (x1 left)

            //right and left hand side of c and c+1 (x and x2)
            if (bc2) {
                h_xr  = qMax(0.0, hxr  - qMax(0.0,  dz_x2 + fb_x2));
                h_x2l = qMax(0.0, hx2l - qMax(0.0, -dz_x2 + fb_x2));
            } else {
                h_x2l = 0.0;
            }
            if (h_xr == 0) {
                vxr = 0;
                uxr = 0;
            }
            if (h_x2l == 0) {
                vx2l = 0;
                ux2l = 0;
            }
            hll_x2 = F_Riemann(h_xr,uxr,vxr, h_x2l,ux2l,vx2l); // c and c+1

            if (br1) {
                h_y1d = qMax(0.0, hy1d - qMax(0.0,  dz_y1 + fb_y1));
                h_yu  = qMax(0.0, hyu  - qMax(0.0, -dz_y1 + fb_y1));
            } else {
                h_y1d = 0.0;
            }
            if (h_yu == 0) {
                uyu = 0;
                vyu = 0;
            }
            if (h_y1d == 0) {
                uy1d = 0;
                vy1d = 0;
            }
            hll_y1 = F_Riemann(h_y1d,vy1d,uy1d, h_yu,vyu,uyu); // r-1 (y1 down) and r (y up)
            // v and u chnaged places for y comnpared to x ? why? is also in swof code

            if (br2) {
                h_yd  = qMax(0.0, hyd  - qMax(0.0,  dz_y2 + fb_y2));// lower side of upper cell
                h_y2u = qMax(0.0, hy2u - qMax(0.0, -dz_y2 + fb_y2));// upper side of lower cell
            } else {
                h_y2u = 0.0;
            }
            if (h_yd == 0) {
                uyd = 0;
                vyd = 0;
            }
            if (h_y2u == 0) {
                uy2u = 0;
                vy2u = 0;
            }
            hll_y2 = F_Riemann(h_yd,vyd,uyd, h_y2u,vy2u,uy2u); // r and r+1

            // determine smallest dt in x and y for each cell
            double dtx = courant_factor*dx/qMax(hll_x1.v[3],hll_x2.v[3]);
            double dty = courant_factor*dy/qMax(hll_y1.v[3],hll_y2.v[3]);
            FloodDT->Drc = qMin(dtx, dty);

            // save the Riemann results in maps, needed for Saint-Venant
            // noite hxl, hxr, hyl, hyr are all equal to H when not using MUSCL, else they have a value based on the minmod limiter
            // so h_xl-hxl is the difference in height between the boundary of the cell and the mid of the cell
            gflowx->Drc = GRAV*0.5*( (h_xl-hxl)*(h_xl+hxl) + (hxr-h_xr)*(hxr+h_xr) + delzcx*(hxl+hxr) ); // delzcx = 0 if not muscl
            gflowy->Drc = GRAV*0.5*( (h_yu-hyu)*(h_yu+hyu) + (hyd-h_yd)*(hyd+h_yd) + delzcy*(hyu+hyd) );
            hllx12_0->Drc = hll_x1.v[0] - hll_x2.v[0];
            hlly12_0->Drc = hll_y1.v[0] - hll_y2.v[0];
            hllx21_1->Drc = hll_x2.v[1] - hll_x1.v[1];
            hllx21_2->Drc = hll_x2.v[2] - hll_x1.v[2];
            hlly21_1->Drc = hll_y2.v[1] - hll_y1.v[1];
            hlly21_2->Drc = hll_y2.v[2] - hll_y1.v[2];
        }
    }} // all cells done

    //find smallest dt in domain
    double dt_req_min = dt;
    #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        dt_req_min = qMin(dt_req_min, FloodDT->Drc);
    }}
    dt_req_min = qMax(TimestepfloodMin, qMin(dt, qMin(dt_req_min, _dt-timesum)));

    return dt_req_min;
}
//-----------------------------------------------------------------------------------------------------------
void TWorld::doSWOFStV(double dt, cTMap *h, cTMap *u, cTMap *v)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double dx = _dx;
        double dy = _dx;
        double Un = 0;
        double Vn = 0;
        double tx = dt/dx;
        double ty = dt/dy;

        double hn = qMax(0.0, h->Drc + tx*(hllx12_0->Drc) + ty*(hlly12_0->Drc));
        // mass balance, hll_....v[0] is the height

        // momentum balance for cells with water

        if(hn > he_ca) { // && qAbs(hn-h->Drc) > 1e-6
            // SWOF solution, delzc1 = 0 when not MUSCL
            double qxn = h->Drc*u->Drc - tx*(hllx21_1->Drc + gflowx->Drc) - ty*hlly21_2->Drc;
            double qyn = h->Drc*v->Drc - tx*hllx21_2->Drc - ty*(hlly21_1->Drc + gflowy->Drc);

            if (SwitchTimeavgV) {
                double nsq1 = (N->Drc)*(N->Drc)*GRAV/qMax(0.0001,std::pow(hn,4.0/3.0));
                double nsq = nsq1 * sqrt(u->Drc*u->Drc + v->Drc*v->Drc) * dt;

                Un = (qxn/(1.0+nsq))/qMax(0.0001,hn);
                Vn = (qyn/(1.0+nsq))/qMax(0.0001,hn);

                double fac = 0.5 + 0.5*qMin(1.0,4*hn)*qMin(1.0,4*hn); // if hn > 1 fac = 1
                fac = fac * exp(- qMax(1.0,dt) / nsq1);
                Un = fac * u->Drc + (1.0-fac) *Un;
                Vn = fac * v->Drc + (1.0-fac) *Vn;
            } else {
                double nsq1 = (N->Drc)*(N->Drc)*GRAV/std::pow(hn,4.0/3.0);
                double nsq = nsq1*sqrt(u->Drc*u->Drc + v->Drc*v->Drc)*dt;
                Un = (qxn/(1.0+nsq))/hn;
                Vn = (qyn/(1.0+nsq))/hn;
            }
        } else {
            // hn < ha
            hn = h->Drc; // if no fluxes then also no change in h
            Un = u->Drc;
            Vn = v->Drc;
        }

        // komt niet meer voor
        if (std::isnan(Un) || std::isnan(Vn)) {
            Un = 0;
            Vn = 0;
        }

        if (fabs(Vn) <= he_ca)
            Vn = 0;
        if (fabs(Un) <= he_ca)
            Un = 0;
        if (Vn == 0 && Un == 0)
            hn = h->Drc;

        h->Drc = hn;
        u->Drc = Un;
        v->Drc = Vn;
    }}
}
//-----------------------------------------------------------------------------------------------------------

