/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2022  Victor Jetten, Bastian van de Bout
**  contact: v.g.jetten@utwente.nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/
//This particular code uses much of the principles of the FullSWOF model
// https://arxiv.org/abs/1204.3210
//https://sourcesup.renater.fr/frs/?group_id=895&release_id=3901#fullswof_2d-_1.10.00-title-content
// the scheme is made suited for parallel processing

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"


//----------------------------------------------------------------------------------------
double TWorld::fullSWOF2openMUSCL(cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.75);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;
    int step = 0;

    sumh = getMass(h, 0);
    //        if (SwitchErosion)
    //            sumS = getMassSed(SSFlood, 0);

    do {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            FloodDT->Drc = dt_max;
            tmd->Drc = 0;

            tma->Drc = h->Drc;
            tmb->Drc = u->Drc;
            tmc->Drc = v->Drc;
        }}

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (h->Drc > F_minWH) {
                tmd->Drc = 1;
                if (c > 0 && !MV(r,c-1)        ) tmd->data[r][c-1] = 1;
                if (c < _nrCols-1 && !MV(r,c+1)) tmd->data[r][c+1] = 1;
                if (r > 0 && !MV(r-1,c)        ) tmd->data[r-1][c] = 1;
                if (r < _nrRows-1 && !MV(r+1,c)) tmd->data[r+1][c] = 1;

                if (c > 1 && !MV(r,c-2)        ) tmd->data[r][c-2] = 1;
                if (c < _nrCols-2 && !MV(r,c+2)) tmd->data[r][c+2] = 1;
                if (r > 1 && !MV(r-2,c)        ) tmd->data[r-2][c] = 1;
                if (r < _nrRows-2 && !MV(r+2,c)) tmd->data[r+2][c] = 1;
            }
        }}

        doSWOFLoop(step, dt_req_min, dt_max, tmd, h, u, v, z);

        // find smallest domain dt
        #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
        FOR_ROW_COL_MV_L {
            dt_req_min = std::min(dt_req_min, FloodDT->Drc);
        }}
        dt_req_min = std::min(dt_req_min, _dt-timesum);

        step += 1; // now we have a good dt min, do the real calculations

SwitchHeun = true;
        // 2nd order, sort of
        if (SwitchHeun) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tma->Drc = h->Drc;
                tmb->Drc = u->Drc;
                tmc->Drc = v->Drc;
            }}
            double d2 = dt_req_min;
            doSWOFLoop(step, dt_req_min, dt_max, tmd, h, u, v, z);

            #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                dt_req_min = std::min(dt_req_min, FloodDT->Drc);
            }}
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            if (d2 > dt_req_min) {
               // qDebug() << d2 << dt_req_min;
                //Heun, see SWOF doc
                #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    double tmp = 0.5*(h->Drc+tma->Drc);
                    if (tmp>=he_ca){
                      double q1 = 0.5*(tma->Drc*tmb->Drc + h->Drc*u->Drc);
                      u->Drc = q1/tmp;
                      double q2 = 0.5*(tma->Drc*tmc->Drc + h->Drc*v->Drc);
                      v->Drc = q2/tmp;
                      h->Drc = tmp;
                    }
                }}
            }
        }

        if (SwitchErosion) {
            SWOFSediment(dt_req_min, h,u,v);
        }

        if (Switch2DDiagonalFlow) {
            SWOFDiagonalFlowNew(dt_req_min, h, u, v);
        }

        timesum += dt_req_min;
        count++; // nr loops

        stop = timesum > _dt-0.001;
        if(count > F_MaxIter)
        stop = true;

    } while (!stop);

    correctMassBalance(sumh, h, 0);

    //qDebug() << _dt/count << count << dt_req_min;
    iter_n = std::max(1,count);
    return(count > 0 ? _dt/count : _dt);

}
//------------------------------------------------------------------------------------------------------
void TWorld::doSWOFLoop(int step, double dt, double dt_max, cTMap *activeCells, cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    //do all flow and state calculations
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (activeCells->Drc > 0) {
            double dx = ChannelAdj->Drc;
            double dy = DX->Drc;
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

            H = h->Drc;
            Z = z->Drc;
            U = u->Drc;
            V = v->Drc;

            bc1 = c > 0 && !MV(r,c-1)        ;
            bc2 = c < _nrCols-1 && !MV(r,c+1);
            br1 = r > 0 && !MV(r-1,c)        ;
            br2 = r < _nrRows-1 && !MV(r+1,c);

            z_x1, z_x2, z_y1, z_y2;
            h_x1, h_x2, h_y1, h_y2;
            u_x1, u_x2, u_y1, u_y2;
            v_x1, v_x2, v_y1, v_y2;
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

            dz_x1 = (Z - z_x1);
            dz_x2 = (z_x2 - Z);
            dz_y1 = (Z - z_y1);
            dz_y2 = (z_y2 - Z);

            if (SwitchFlowBarriers) {
                fb_x1 = bc1 ? std::max(FlowBarrierW->Drc, FlowBarrierE->data[r][c-1]) : FlowBarrierW->Drc;
                fb_x2 = bc2 ? std::max(FlowBarrierE->Drc, FlowBarrierE->data[r][c+1]) : FlowBarrierE->Drc;
                fb_y1 = br1 ? std::max(FlowBarrierN->Drc, FlowBarrierS->data[r-1][c]) : FlowBarrierN->Drc;
                fb_y2 = br2 ? std::max(FlowBarrierS->Drc, FlowBarrierN->data[r+1][c]) : FlowBarrierS->Drc;
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

            //#include "lisMUSCL.in"
            bool b2c1 ,b2c2 ,b2r1 ,b2r2;
            double h_xx1, h_xx2, u_xx1, u_xx2, v_xx1, v_xx2;
            double h_yy1, h_yy2, u_yy1, u_yy2, v_yy1, v_yy2;
            double dh, du, dv, dz_h;
            double delta_h1, delta_h2, delta_h3, delta_h4;
            double delta_u1, delta_u2, delta_u3, delta_u4;
            double delta_v1, delta_v2, delta_v3, delta_v4;

            b2c1 = c > 1 && !MV(r,c-2)        ;
            b2c2 = c < _nrCols-2 && !MV(r,c+2);
            b2r1 = r > 1 && !MV(r-2,c)        ;
            b2r2 = r < _nrRows-2 && !MV(r+2,c);

            if(b2c1) {
                h_xx1 = h->data[r][c-2];
                u_xx1 = u->data[r][c-2];
                v_xx1 = v->data[r][c-2];
            }
            if(b2c2) {
                h_xx2 = h->data[r][c+2];
                u_xx2 = u->data[r][c+2];
                v_xx2 = v->data[r][c+2];
            }
            if(b2r1) {
                h_yy1 = h->data[r-2][c];
                u_yy1 = u->data[r-2][c];
                v_yy1 = v->data[r-2][c];
            }
            if(b2r2) {
                h_yy2 = h->data[r+2][c];
                u_yy2 = u->data[r+2][c];
                v_yy2 = v->data[r+2][c];
            }

            if(b2c1 && b2c2) {

                // x-1-x-2   x-x-1  x+1-x   x+2-x+1        always right minus left
                delta_h1 = h_x1 - h_xx1; delta_h2 = H-h_x1; delta_h3 = h_x2-H; delta_h4 = h_xx2-h_x2;
                delta_u1 = u_x1 - u_xx1; delta_u2 = U-u_x1; delta_u3 = u_x2-U; delta_u4 = u_xx2-u_x2;
                delta_v1 = v_x1 - v_xx1; delta_v2 = V-v_x1; delta_v3 = v_x2-V; delta_v4 = v_xx2-v_x2;

                // center cell
                dh = limiter(delta_h2, delta_h3);
                du = limiter(delta_u2, delta_u3);
                dv = limiter(delta_v2, delta_v3);
                hxl = H - 0.5*dh;
                hxr = H + 0.5*dh;
                uxl = U - 0.5*du*hxl/H;
                uxr = U + 0.5*du*hxr/H;
                vxl = V - 0.5*dv*hxl/H;
                vxr = V + 0.5*dv*hxr/H;

                dz_h = limiter(delta_h2 + (Z-dz_x1), delta_h3 + (dz_x2-Z));
                delzcx = (Z+0.5*(dz_h-dh))-(Z+0.5*(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!2*(dz_h - dh); //

                // left hand cell, right boundary
                dh = limiter(delta_h1, delta_h2);
                du = limiter(delta_u1, delta_u2);
                dv = limiter(delta_v1, delta_v2);
                hx1r = h_x1 + 0.5*dh;
                ux1r = u_x1 + 0.5*du*hxl/H;
                vx1r = v_x1 + 0.5*dv*hxl/H;

                // right hand cell, left boundary
                dh = limiter(delta_h3, delta_h4);
                du = limiter(delta_u3, delta_u4);
                dv = limiter(delta_v3, delta_v4);
                hx2l = h_x2 - 0.5*dh;
                ux2l = u_x2 - 0.5*du*hxr/H;
                vx2l = v_x2 - 0.5*dv*hxr/H;
            }

            if (b2r1 && b2r2) {
                // vertical, direction from up to down
                // y-1 - y-2   y-y-1  y+1-y   y+2-y+1        always down minus up
                delta_h1 = h_y1 - h_yy1; delta_h2 = H-h_y1; delta_h3 = h_y2-H; delta_h4 = h_yy2-h_y2;
                delta_u1 = u_y1 - u_yy1; delta_u2 = U-u_y1; delta_u3 = u_y2-U; delta_u4 = u_yy2-u_y2;
                delta_v1 = v_y1 - v_yy1; delta_v2 = V-v_y1; delta_v3 = v_y2-V; delta_v4 = v_yy2-v_y2;

                // center cell
                dh = limiter(delta_h2, delta_h3);
                du = limiter(delta_u2, delta_u3);
                dv = limiter(delta_v2, delta_v3);
                hyu = H - 0.5*dh;
                hyd = H + 0.5*dh;
                uyu = U - 0.5*du*hyu/H;
                uyd = U + 0.5*du*hyd/H;
                vyu = V - 0.5*dv*hyu/H;
                vyd = V + 0.5*dv*hyd/H;

                dz_h = limiter(delta_h1 + (Z-dz_y1), delta_h2 + (dz_y2-Z));
                delzcy = (Z+0.5*(dz_h-dh))-(Z+0.5*(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!2*(dz_h-dh);

                // upper cell, down boundary
                dh = limiter(delta_h1, delta_h2);
                du = limiter(delta_u1, delta_u2);
                dv = limiter(delta_v1, delta_v2);
                hy1d = h_y1 + 0.5*dh;
                uy1d = u_y1 + 0.5*du*hyu/H;
                vy1d = v_y1 + 0.5*dv*hyu/H;

                // lower cell, up boundary
                dh = limiter(delta_h3, delta_h4);
                du = limiter(delta_u3, delta_u4);
                dv = limiter(delta_v3, delta_v4);
                hy2u = h_y2 - 0.5*dh;
                uy2u = u_y2 - 0.5*du*hyd/H;
                vy2u = v_y2 - 0.5*dv*hyd/H;
            }
            //########### calculate Riemann valaues for all four boundaries of a cell ############

            // if muscl H and h_x1 etc become Hx1l and hx1r
            // z is blocking to prevent flow when water is flat and Z is not flat, described in article SWOF
            // barrier is ourown additiona, to vcreate flood walls.

            //left and right hand side of c and c-1 (x and x1)
            if (bc1) {
                h_x1r = std::max(0.0, hx1r - std::max(0.0,  dz_x1 + fb_x1)); //rechts van c-1
                h_xl  = std::max(0.0, hxl  - std::max(0.0, -dz_x1 + fb_x1)); //links van het midden
            } else {
                h_x1r=ux1r=vx1r=0.0;
            } // if !inside = boundary
            hll_x1 = F_Riemann(h_x1r,ux1r,vx1r, h_xl,uxl,vxl); // c-1 (x1 right) and c (x1 left)

            //right and left hand side of c and c+1 (x and x2)
            if (bc2) {
                h_xr  = std::max(0.0, hxr  - std::max(0.0,  dz_x2 + fb_x2));
                h_x2l = std::max(0.0, hx2l - std::max(0.0, -dz_x2 + fb_x2));
            } else {
                h_x2l=ux2l=vx2l=0.0;
            }
            hll_x2 = F_Riemann(h_xr,uxr,vxr, h_x2l,ux2l,vx2l); // c and c+1

            if (br1) {
                h_y1d = std::max(0.0, hy1d - std::max(0.0,  dz_y1 + fb_y1));
                h_yu  = std::max(0.0, hyu  - std::max(0.0, -dz_y1 + fb_y1));
            } else {
                h_y1d=vy1d=uy1d=0.0;
            }
            hll_y1 = F_Riemann(h_y1d,vy1d,uy1d, h_yu,vyu,uyu); // r-1 (y1 down) and r (y up)
            // v and u chnaged places for y comnpared to x ? why? is also in swof code

            if (br2) {
                h_yd  = std::max(0.0, hyd  - std::max(0.0,  dz_y2 + fb_y2));
                h_y2u = std::max(0.0, hy2u - std::max(0.0, -dz_y2 + fb_y2));
            } else {
                h_y2u=vy2u=uy2u=0.0;
            }
            hll_y2 = F_Riemann(h_yd,vyd,uyd, h_y2u,vy2u,uy2u); // r and r+1

            // determine smallest dt in x and y for each cell
            double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
            double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]); // v[3] is max U and V in x and y
            FloodDT->Drc = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));
            // if step = 0 do not calculate new fluxes and states yet because the first dt is always dt_max
            // find a smallest dt of the flow domain first

            //########### after finding the smallest dt, do saint venant eq)
            if (step > 0) {
                double Un, Vn;
                double tx = dt/dx;
                double ty = dt/dy;

                double hn = std::max(0.0, H + tx*(hll_x1.v[0]-hll_x2.v[0]) + ty*(hll_y1.v[0]-hll_y2.v[0]));
                // mass balance, hll_....v[0] is the height

                // momentum balance for cells with water
                if(hn > he_ca) {
                    // SWOF solution, delzc1 = 0 when not MUSCL
                    double gflow_x = GRAV*0.5*( (h_xl-hxl)*(h_xl+hxl) + (hxr-h_xr)*(hxr+h_xr) + delzcx*(hxl+hxr)); // delzcx = 0 is not muscl
                    double gflow_y = GRAV*0.5*( (h_yu-hyu)*(h_yu+hyu) + (hyd-h_yd)*(hyd+h_yd) + delzcy*(hyu+hyd));

                    double qxn = H * U - tx*(hll_x2.v[1] - hll_x1.v[1] + gflow_x) - ty*(hll_y2.v[2] - hll_y1.v[2]);
                    double qyn = H * V - tx*(hll_x2.v[2] - hll_x1.v[2]) - ty*(hll_y2.v[1] - hll_y1.v[1] + gflow_y);

                    double vsq = sqrt(U*U + V*V);
                    double nsq1 = (N->Drc)*(N->Drc)*GRAV/pow(hn,4.0/3.0);//std::max(0.0001,pow(hn,4.0/3.0)); //
                    double nsq = nsq1*vsq*dt;

                    //Un = (qxn/(1.0+nsq))/std::max(0.0001,hn);
                    //Vn = (qyn/(1.0+nsq))/std::max(0.0001,hn);
                    Un = (qxn/(1.0+nsq))/hn;
                    Vn = (qyn/(1.0+nsq))/hn;

                    if (SwitchTimeavgV) {
                        double fac = 0.5 + 0.5*std::min(1.0,4*hn)*std::min(1.0,4*hn);
                        fac = fac * exp(- std::max(1.0,dt) / nsq1);
                        Un = fac * U + (1.0-fac) *Un;
                        Vn = fac * V + (1.0-fac) *Vn;
                    }

                } else { // hn < ha
                    hn = H; // if no fluxes then also no change in h
                    Un = 0;
                    Vn = 0;
                }

                // dan maar even met geweld!
                if (std::isnan(Un) || std::isnan(Vn)  )
                {
                    Un = 0;
                    Vn = 0;
                }
                if (FlowBoundaryType == 0 || (FlowBoundaryType == 2 && FlowBoundary->Drc == 0)) {
                    if (DomainEdge->Drc == 4 && Un < 0) Un = 0;
                    if (DomainEdge->Drc == 6 && Un > 0) Un = 0;
                    if (DomainEdge->Drc == 2 && Vn > 0) Vn = 0;
                    if (DomainEdge->Drc == 8 && Vn < 0) Vn = 0;
                }
                if (Vn == 0 && Un == 0)
                    hn = H;

                h->Drc = hn;
                u->Drc = Un;
                v->Drc = Vn;

            } // step > 0

        } // tmd > 0, active cells
    }} // all cells done
}

