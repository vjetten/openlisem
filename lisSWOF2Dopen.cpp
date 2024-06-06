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


// uses functions in auxillary : lisSWOFaux.cpp

//-------------------------------------------------------------------------------------------------
double TWorld::fullSWOF2open(cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.75);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;
    int step = 0;

//    Qout.clear();
//    FOR_ROW_COL_LDD5 {
//       Qout << 0.0;
//    }}

    if (startFlood)
    {

        sumh = getMass(h, 0);
//        if (SwitchErosion)
//            sumS = getMassSed(SSFlood, 0);

        do {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                hs->Drc = h->Drc;
                FloodDT->Drc = dt_max;
                //FloodT->Drc = 0;
                tmd->Drc = 0;
            }}

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (hs->Drc > F_minWH) {
                    tmd->Drc = 1;
                    if (c > 0 && !MV(r,c-1)        ) tmd->data[r][c-1] = 1;
                    if (c < _nrCols-1 && !MV(r,c+1)) tmd->data[r][c+1] = 1;
                    if (r > 0 && !MV(r-1,c)        ) tmd->data[r-1][c] = 1;
                    if (r < _nrRows-1 && !MV(r+1,c)) tmd->data[r+1][c] = 1;

                    if (c > 0 && r > 0 && !MV(r-1,c-1)                ) tmd->data[r-1][c-1]=1;
                    if (c < _nrCols-1 && r < _nrRows-1 && !MV(r+1,c+1)) tmd->data[r+1][c+1]=1;
                    if (r > 0 && c < _nrCols-1 && !MV(r-1,c+1)        ) tmd->data[r-1][c+1]=1;
                    if (c > 0 && r < _nrRows-1 && !MV(r+1,c-1)        ) tmd->data[r+1][c-1]=1;
                    if (SwitchMUSCL) {
                        if (c > 1 && !MV(r,c-2)        ) tmd->data[r][c-2] = 1;
                        if (c < _nrCols-2 && !MV(r,c+2)) tmd->data[r][c+2] = 1;
                        if (r > 1 && !MV(r-2,c)        ) tmd->data[r-2][c] = 1;
                        if (r < _nrRows-2 && !MV(r+2,c)) tmd->data[r+2][c] = 1;
                    }
                }
            }}

            //do all flow and state calculations
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {

                if (tmd->Drc > 0) {
                        //double dt = FloodDT->Drc; //dt_req_min;
                    double dt = dt_req_min;
                    double Un, Vn;

                    //FloodT->Drc += FloodDT->Drc;

                    vec4 hll_x1;
                    vec4 hll_x2;
                    vec4 hll_y1;
                    vec4 hll_y2;

                    double dx = _dx;//ChannelAdj->Drc;
                    double dy = _dx;//DX->Drc;

                    double H = hs->Drc;
                    double n = N->Drc;
                    double Z = z->Drc;
                    double U = u->Drc;
                    double V = v->Drc;

                    bool bc1 = c > 0 && !MV(r,c-1)        ;
                    bool bc2 = c < _nrCols-1 && !MV(r,c+1);
                    bool br1 = r > 0 && !MV(r-1,c)        ;
                    bool br2 = r < _nrRows-1 && !MV(r+1,c);

                    double z_x1 =  bc1 ? z->data[r][c-1] : Z;
                    double z_x2 =  bc2 ? z->data[r][c+1] : Z;
                    double z_y1 =  br1 ? z->data[r-1][c] : Z;
                    double z_y2 =  br2 ? z->data[r+1][c] : Z;

                    double h_x1 =  bc1 ? hs->data[r][c-1] : H;  //??? hs???
                    double h_x2 =  bc2 ? hs->data[r][c+1] : H;
                    double h_y1 =  br1 ? hs->data[r-1][c] : H;
                    double h_y2 =  br2 ? hs->data[r+1][c] : H;

                    double u_x1 = bc1 ? u->data[r][c-1] : U;
                    double u_x2 = bc2 ? u->data[r][c+1] : U;
                    double u_y1 = br1 ? u->data[r-1][c] : U;
                    double u_y2 = br2 ? u->data[r+1][c] : U;

                    double v_x1 = bc1 ? v->data[r][c-1] : V;
                    double v_x2 = bc2 ? v->data[r][c+1] : V;
                    double v_y1 = br1 ? v->data[r-1][c] : V;
                    double v_y2 = br2 ? v->data[r+1][c] : V;

                    double h_xx1, h_xx2, u_xx1, u_xx2, v_xx1, v_xx2;
                    double h_yy1, h_yy2, u_yy1, u_yy2, v_yy1, v_yy2;

                    if (SwitchMUSCL) {
                        h_xx1=h_x1; h_xx2=h_x2; u_xx1=u_x1; u_xx2=u_x2; v_xx1=v_x1; v_xx2=v_x2;
                        h_yy1=h_y1; h_yy2=h_y2; u_yy1=u_y1; u_yy2=u_y2; v_yy1=v_y1; v_yy2=v_y2;

                        if(c > 1 && !MV(r,c-2)) {
                            h_xx1 = hs->data[r][c-2];
                            u_xx1 =  u->data[r][c-2];
                            v_xx1 =  v->data[r][c-2];
                        }
                        if(c < _nrCols-2 && !MV(r,c+2)) {
                            h_xx2 = hs->data[r][c+2];
                            u_xx2 =  u->data[r][c+2];
                            v_xx2 =  v->data[r][c+2];
                        }
                        if(r > 1 && !MV(r-2,c)) {
                            h_yy1 = hs->data[r-2][c];
                            u_yy1 =  u->data[r-2][c];
                            v_yy1 =  v->data[r-2][c];
                        }
                        if(r < _nrRows-2 && !MV(r+2,c)) {
                            h_yy2 = hs->data[r+2][c];
                            u_yy2 =  u->data[r+2][c];
                            v_yy2 =  v->data[r+2][c];
                        }
                    }

                    double fb_x1=0,fb_x2=0,fb_y1=0,fb_y2=0;
                    if (SwitchFlowBarriers) {
                        fb_x1 = bc1 ? std::max(FlowBarrierW->Drc, FlowBarrierE->data[r][c-1]) : FlowBarrierW->Drc;
                        fb_x2 = bc2 ? std::max(FlowBarrierE->Drc, FlowBarrierE->data[r][c+1]) : FlowBarrierE->Drc;
                        fb_y1 = br1 ? std::max(FlowBarrierN->Drc, FlowBarrierS->data[r-1][c]) : FlowBarrierN->Drc;
                        fb_y2 = br2 ? std::max(FlowBarrierS->Drc, FlowBarrierN->data[r+1][c]) : FlowBarrierS->Drc;
                    }

                    double dz_x1 = (Z - z_x1);
                    double dz_x2 = (z_x2 - Z);
                    double dz_y1 = (Z - z_y1);
                    double dz_y2 = (z_y2 - Z);

                    // muscl
                    double delzcx = 0;
                    double delzcy = 0;

                    // all boundaries are done in all direction (so horizontal in x and y and vertical in x and why
                    // not only horizontal in x and vertical in y, because these are vectors so for each boundary there is an x and y direction

                    double hx1r, hxl, hxr, hx2l;  // c-1, c, c+1
                    double hy1d, hyu, hyd, hy2u;  // r-1, r, r+1

                    double ux1r, uxl, uxr, ux2l;
                    double uy1d, uyu, uyd, uy2u;

                    double vx1r, vxl, vxr, vx2l;
                    double vy1d, vyu, vyd, vy2u;

                    // all boundaries are done in all direction (so horizontal in x and y and vertical in x and why
                    // not only horizontal in x and vertical in y, because these are vectors so for each boundary there is an x and y direction

                    // non muscl solution, cell centres for boundaries in x and y directions
                    hx1r = h_x1; hxl = H; hxr = H; hx2l = h_x2;  // x-1 r; x l; x r; x+1 l
                    ux1r = u_x1; uxl = U; uxr = U; ux2l = u_x2;
                    vx1r = v_x1; vxl = V; vxr = V; vx2l = v_x2;

                    hy1d = h_y1; hyu = H; hyd = H; hy2u = h_y2;
                    uy1d = u_y1; uyu = U; uyd = U; uy2u = u_y2;
                    vy1d = v_y1; vyu = V; vyd = V; vy2u = v_y2;

                    //  MUSCL: on the 4 boundaties of a gridcell interpolate from the center values
                    if (SwitchMUSCL) {
                        double dh, du, dv, dz_h;
                        double delta_h1, delta_h2, delta_h3, delta_h4;
                        double delta_u1, delta_u2, delta_u3, delta_u4;
                        double delta_v1, delta_v2, delta_v3, delta_v4;
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
                        delzcx = 2*(dz_h - dh); // Z+(dz_h-dh)- (Z+(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!

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
                        delzcy = 2*(dz_h-dh);//Z+(dz_h-dh)- (Z+(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!

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
                    double h_x1r = std::max(0.0, hx1r - std::max(0.0,  dz_x1 + fb_x1)); //rechts van c-1
                    double h_xl  = std::max(0.0, hxl  - std::max(0.0, -dz_x1 + fb_x1)); //links van het midden
                    if(!bc1) { h_x1r=ux1r=vx1r=0.0; } // if !inside = boundary
                    hll_x1 = F_Riemann(h_x1r,ux1r,vx1r, h_xl,uxl,vxl); // c-1 (x1 right) and c (x1 left)

                    //right and left hand side of c and c+1 (x and x2)
                    double h_xr  = std::max(0.0, hxr  - std::max(0.0,  dz_x2 + fb_x2));
                    double h_x2l = std::max(0.0, hx2l - std::max(0.0, -dz_x2 + fb_x2));
                    if(!bc2) { h_x2l=ux2l=vx2l=0.0;}
                    hll_x2 = F_Riemann(h_xr,uxr,vxr, h_x2l,ux2l,vx2l); // c and c+1

                    double h_y1d = std::max(0.0, hy1d - std::max(0.0,  dz_y1 + fb_y1));
                    double h_yu  = std::max(0.0, hyu  - std::max(0.0, -dz_y1 + fb_y1));
                    if (!br1) {h_y1d=vy1d=uy1d=0.0;}
                    hll_y1 = F_Riemann(h_y1d,vy1d,uy1d, h_yu,vyu,uyu); // r-1 (y1 down) and r (y up)
                    // v and u chnaged places for y comnpared to x ? why? is also in swof code

                    double h_yd  = std::max(0.0, hyd  - std::max(0.0,  dz_y2 + fb_y2));
                    double h_y2u = std::max(0.0, hy2u - std::max(0.0, -dz_y2 + fb_y2));
                    if(!br2) { h_y2u=vy2u=uy2u=0.0; }
                    hll_y2 = F_Riemann(h_yd,vyd,uyd, h_y2u,vy2u,uy2u); // r and r+1

                    // determine smallest dt in x and y for each cell
                    double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                    double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]); // v[3] is max U and V in x and y

                    double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));
                    FloodDT->Drc = dt_req;
                    // dt does not need to be a map, left over from earlier code
                    // if step = 0 do not calculate new fluxes and states yet because the first dt is always dt_max
                    // find a smallest dt of the flow domain first

                    //########### after finding the smallest dt, do the st venant eq)
                    if (step > 0) {

                        double tx = dt/dx;
                        double ty = dt/dy;

                        double hn = std::max(0.0, H + dt/_dx*(hll_x1.v[0]-hll_x2.v[0] + hll_y1.v[0]-hll_y2.v[0]));
                        // mass balance, hll_....v[0] is the height

                        // momentum balance for cells with water
                        if(hn > he_ca) {
                            // SWOF solution, delzc1 = 0 when not MUSCL
                            double gflow_x = GRAV*0.5*( (h_xl-hxl)*(h_xl+hxl) + (hxr-h_xr)*(hxr+h_xr) + delzcx*(hxl+hxr));
                            double gflow_y = GRAV*0.5*( (h_yu-hyu)*(h_yu+hyu) + (hyd-h_yd)*(hyd+h_yd) + delzcy*(hyu+hyd));

                            double qxn = H * U - tx*(hll_x2.v[1] - hll_x1.v[1] + gflow_x) - ty*(hll_y2.v[2] - hll_y1.v[2]);
                            double qyn = H * V - tx*(hll_x2.v[2] - hll_x1.v[2]) - ty*(hll_y2.v[1] - hll_y1.v[1] + gflow_y);

                            double vsq = sqrt(U*U + V*V);
                            double nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.0001,pow(hn,4.0/3.0)); //pow(hn,4.0/3.0);//
                            double nsq = nsq1*vsq*dt;

                            Un = (qxn/(1.0+nsq))/std::max(0.0001,hn);
                            Vn = (qyn/(1.0+nsq))/std::max(0.0001,hn);

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

                            if (DomainEdge->Drc == 4 && Un < 0) {
                                Un = 0;
                            }
                            if (DomainEdge->Drc == 6 && Un > 0) {
                                Un = 0;
                            }
                            if (DomainEdge->Drc == 2 && Vn > 0) {
                                Vn = 0;
                            }
                            if (DomainEdge->Drc == 8 && Vn < 0) {
                                Vn = 0;
                            }

                        }
                        if (Vn == 0 && Un == 0)
                            hn = H;

                        h->Drc = hn;
                        u->Drc = Un;
                        v->Drc = Vn;

                    } // step > 0
                } // tmd > 0, active cells + 1
            }}

            // find smallest domain dt
            #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                dt_req_min = std::min(dt_req_min, FloodDT->Drc);
            }}

            dt_req_min = std::min(dt_req_min, _dt-timesum);

            if (step > 0) {

               if (SwitchErosion) {
                    SWOFSediment(dt_req_min, h,u,v);
                }

                if (Switch2DDiagonalFlow) {
                    if (Switch2DDiagonalFlowNew)
                        SWOFDiagonalFlowNew(dt_req_min, h, u, v);
                    else
                        SWOFDiagonalFlow(dt_req_min, h, u, v); //old, not used
                }


                timesum += dt_req_min;
                count++; // nr loops
            }

            step += 1; // now we have a good dt min, do the real calculations

            stop = timesum > _dt-0.001;
            if(count > F_MaxIter)
            stop = true;

        } while (!stop);

        correctMassBalance(sumh, h, 0);

    } // if floodstart

    //qDebug() << _dt/count << count << dt_req_min;
    iter_n = std::max(1,count);
    return(count > 0 ? _dt/count : _dt);
}
