

/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten, Bastian van de Bout
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
//https://sourcesup.renater.fr/frs/?group_id=895&release_id=3901#fullswof_2d-_1.10.00-title-content
#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

//-------------------------------------------------------------------------------------------------



// force flow when a diagonal solution exists and a DEM blockage is present
void TWorld::SWOFDiagonalFlowNew(double dt_req_min, cTMap *h, cTMap *vx, cTMap *vy)
{
	
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tma->Drc = 0;
        tmb->Drc = 0;
        tmc->Drc = 0;
    }}

    bool doit = false;

    #pragma omp parallel for num_threads(userCores)
    for(long i_= 0; i_ < dcr_.size(); i_++) {

        int r = dcr_[i_].r;
        int c = dcr_[i_].c;

        if (h->Drc > F_pitValue) {
            int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1,  0,  1};
            int dy[10] = {0,  1, 1, 1,  0, 0, 0, -1, -1, -1};
            doit = true;

            vec4 rec;
            int ldd = dcr_[i_].ldd;
            int rr = r+dy[ldd];
            int cr = c+dx[ldd];

            // h downstream cannot be updated inside parallel loop!
            // save these values and add later
            if (h->Drcr < h->Drc) {
                // 1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
                rec = F_Riemann(h->Drc, vx->Drc, vy->Drc, h->Drcr, vx->Drcr, vy->Drcr);
                double flux = std::abs(rec.v[0]);
                double dH = std::min(h->Drc *0.9, flux*dt_req_min/_dx);

                h->Drc -= dH;
                tmc->Drcr += dH;

                if (SwitchErosion) {
                    double dS = std::min(0.9*SSFlood->Drc, dH*CHAdjDX->Drc*SSCFlood->Drc);
                    SSFlood->Drc -= dS;
                    tma->Drcr += dS;
                    if (SwitchUse2Phase) {
                        double dBL = std::min(0.9*BLFlood->Drc, dH*CHAdjDX->Drc*BLCFlood->Drc);
                        BLFlood->Drc -= dBL;
                        tmb->Drcr += dBL;
                    }
                }
            }
        }
    }

    if (doit) {

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            h->Drc += tmc->Drc;
        }}

        if (SwitchErosion) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                SSFlood->Drc += tma->Drc;
                if (SwitchUse2Phase)
                    BLFlood->Drc += tmb->Drc;
            }}
        }
    }
}

//-------------------------------------------------------------------------------------------------
// force flow when a diagonal solution exists and a DEM blockage is present
void TWorld::SWOFDiagonalFlow(double dt_req_min, cTMap *h, cTMap *vx, cTMap *vy)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tma->Drc = 0;
        tmb->Drc = 0;
        tmc->Drc = 0;
    }}

    bool doit = false;

    #pragma omp parallel for num_threads(userCores)
    for(long i_= 0; i_ < dcr_.size(); i_++) {

        int r = dcr_[i_].r;
        int c = dcr_[i_].c;

        if (h->Drc > F_pitValue) {
            int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1,  0,  1};
            int dy[10] = {0,  1, 1, 1,  0, 0, 0, -1, -1, -1};
            doit = true;

            vec4 rec;
            int ldd = dcr_[i_].ldd;
            int rr = r+dy[ldd];
            int cr = c+dx[ldd];

            // h downstream cannot be updated inside parallel loop!
            // save these values and add later
            if (h->Drcr < h->Drc) {
                // 1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
                rec = F_Riemann(h->Drc, vx->Drc, vy->Drc, h->Drcr, vx->Drcr, vy->Drcr);
                double flux = std::abs(rec.v[0]);
                double dH = std::min(h->Drc *0.9, flux*dt_req_min/_dx);

                h->Drc -= dH;
                //h->Drcr += dH;
                tmc->Drcr += dH;
                //Qdiag->Drc = flux;

                if (SwitchErosion) {
                    double dS = std::min(0.9*SSFlood->Drc, dH*CHAdjDX->Drc*SSCFlood->Drc);
                    SSFlood->Drc -= dS;
                    //SSFlood->Drcr += dS;
                    tma->Drcr += dS;
                    if (SwitchUse2Phase) {
                        double dBL = std::min(0.9*BLFlood->Drc, dH*CHAdjDX->Drc*BLCFlood->Drc);
                        BLFlood->Drc -= dBL;
                        tmb->Drcr += dBL;
                    }
                }
            }
        }
    }

    if (doit) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (SwitchErosion) {
                SSFlood->Drc += tma->Drc;
                if (SwitchUse2Phase)
                    BLFlood->Drc += tmb->Drc;
            }
            h->Drc += tmc->Drc;
        }}
    }

}
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
                    hx1r = h_x1; hxl = H; hxr = H; hx2l = h_x2;
                    ux1r = u_x1; uxl = U; uxr = U; ux2l = u_x2;
                    vx1r = v_x1; vxl = V; vxr = V; vx2l = v_x2;

                    hy1d = h_y1; hyu = H; hyd = H; hy2u = h_y2;
                    uy1d = u_y1; uyu = U; uyd = U; uy2u = u_y2;
                    vy1d = v_y1; vyu = V; vyd = V; vy2u = v_y2;

                    //  MUSCL: on the 4 boundaties of a gridcell interpolate from the center values
                    if (SwitchMUSCL) {
                        double minh = 0.0001;
                        double dh, du, dv, dz_h;
                        double tmph, tmpu, tmpv;
                        double delta_h1, delta_h2;
                        double delta_u1, delta_u2;
                        double delta_v1, delta_v2;

                        // center cell
                        delta_h1 = H-h_x1;
                        delta_u1 = U-u_x1;
                        delta_v1 = V-v_x1;
                        delta_h2 = h_x2-H;
                        delta_u2 = u_x2-U;
                        delta_v2 = v_x2-V;
                        tmph = delta_h2;
                        tmpu = delta_u2;
                        tmpv = delta_v2;
                        dh = limiter(delta_h1, delta_h2);
                        du = limiter(delta_u1, delta_u2);
                        dv = limiter(delta_v1, delta_v2);
                        hxl = H - 0.5*dh;
                        hxr = H + 0.5*dh;
                        uxl = U - 0.5*du*hxl/H;
                        uxr = U + 0.5*du*hxr/H;
                        vxl = V - 0.5*dv*hxl/H;
                        vxr = V + 0.5*dv*hxr/H;
                        // if (H > minh) {
                        // uxl = U - 0.5*du*hxr/H;
                        // uxr = U + 0.5*du*hxl/H;
                        // vxl = V - 0.5*dv*hxr/H;
                        // vxr = V + 0.5*dv*hxl/H;
                        // } else {
                        //     uxl = U - 0.5*hxr*du;
                        //     uxr = U + 0.5*hxl*du;
                        //     vxl = V - 0.5*hxr*dv;
                        //     vxr = V + 0.5*hxl*dv;
                        // }

                        dz_h = limiter(delta_h1 + (Z-dz_x1), delta_h2 + (dz_x2-Z));
                        delzcx = Z+(dz_h-dh)- (Z+(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!

                        // left hand cell, right boundary
                        if(c > 1 && !MV(r,c-2)) {
                            delta_h2 = delta_h1;
                            delta_u2 = delta_u1;
                            delta_v2 = delta_v1;
                            delta_h1 = h_x1 - h_xx1;
                            delta_u1 = u_x1 - u_xx1;
                            delta_v1 = u_x1 - v_xx1;
                            dh = limiter(delta_h1, delta_h2);
                            du = limiter(delta_u1, delta_u2);
                            dv = limiter(delta_v1, delta_v2);
                            hx1r = h_x1 + 0.5*dh;
                            ux1r = u_x1 + 0.5*du*hxl/H;
                            vx1r = v_x1 + 0.5*dv*hxl/H;
                            // double hx1l = h_x1 - 0.5*dh;
                            // if (H > minh) {
                            //     ux1r = u_x1 + 0.5*du*hx1r/H;
                            //     vx1r = v_x1 + 0.5*dv*hx1r/H;
                            // } else {
                            //     ux1r = u_x1 + hx1l*du*0.5;
                            //     vx1r = v_x1 + hx1l*dv*0.5;
                            // }
                        }

                        // right hand cell, left boundary
                        if(c < _nrCols-2 && !MV(r,c+2)) {
                            delta_h2 = h_xx2 - h_x2;
                            delta_u2 = u_xx2 - u_x2;
                            delta_v2 = v_xx2 - v_x2;
                            delta_h1 = tmph;
                            delta_u1 = tmpu;
                            delta_v1 = tmpv;
                            dh = limiter(delta_h1, delta_h2);
                            du = limiter(delta_u1, delta_u2);
                            dv = limiter(delta_v1, delta_v2);
                            hx2l = h_x2 - 0.5*dh;
                            ux2l = u_x2 - 0.5*du*hxr/H;
                            vx2l = v_x2 - 0.5*dv*hxr/H;
                            // double hx2r = h_x2 - 0.5*dh;
                            // if (H > minh) {
                            //     ux2l = u_x2 - 0.5*du*hx2l/H;
                            //     vx2l = v_x2 - 0.5*dv*hx2l/H;
                            // } else {
                            //     ux2l = u_x2 - 0.5*hx2r*du;
                            //     vx2l = v_x2 - 0.5*hx2r*dv;
                            // }
                        }

                        // vertical, direction from up to down
                        // center cell
                        delta_h1 = H-h_y1;
                        delta_u1 = U-u_y1;
                        delta_v1 = V-v_y1;
                        delta_h2 = h_y2-H;
                        delta_u2 = u_y2-U;
                        delta_v2 = v_y2-V;
                        tmph = delta_h2;
                        tmpu = delta_u2;
                        tmpv = delta_v2;
                        dh = limiter(delta_h1, delta_h2);
                        du = limiter(delta_u1, delta_u2);
                        dv = limiter(delta_v1, delta_v2);
                        hyu = H - 0.5*dh;
                        hyd = H + 0.5*dh;
                        uyu = U - 0.5*du*hyu/H;
                        uyd = U + 0.5*du*hyd/H;
                        vyu = V - 0.5*dv*hyu/H;
                        vyd = V + 0.5*dv*hyd/H;
                        // if (H > minh) {
                        //     uyu = U - 0.5*du*hyd/H;
                        //     uyd = U + 0.5*du*hyu/H;
                        //     vyu = V - 0.5*dv*hyd/H;
                        //     vyd = V + 0.5*dv*hyu/H;
                        // } else {
                        //     uyu = U - 0.5*hyd*du;
                        //     uyd = U + 0.5*hyu*du;
                        //     vyu = V - 0.5*hyd*dv;
                        //     vyd = V + 0.5*hyu*dv;
                        // }
                        dz_h = limiter(delta_h1 + (Z-dz_y1), delta_h2 + (dz_y2-Z));
                        delzcy = Z+(dz_h-dh)- (Z+(dh-dz_h));// = (dz_h-dh)-(dh-dz_h) = 2*dz_h-2*dh; //!!!!

                        // upper cell, down boundary
                        if(r > 1 && !MV(r-2,c)) {
                            delta_h2 = delta_h1;
                            delta_u2 = delta_u1;
                            delta_v2 = delta_v1;
                            delta_h1 = h_y1 - h_yy1;
                            delta_u1 = u_y1 - u_yy1;
                            delta_v1 = u_y1 - v_yy1;
                            dh = limiter(delta_h1, delta_h2);
                            du = limiter(delta_u1, delta_u2);
                            dv = limiter(delta_v1, delta_v2);
                            hy1d = h_y1 + 0.5*dh;
                            uy1d = u_y1 + 0.5*du*hyu/H;
                            vy1d = v_y1 + 0.5*dv*hyu/H;
                            // double hy1u = h_y1 - 0.5*dh;
                            // if (H > minh) {
                            //     uy1d = u_y1 + 0.5*du*hy1d/H;
                            //     vy1d = v_y1 + 0.5*dv*hy1d/H;
                            // } else {
                            //     uy1d = u_y1 + hy1u*du*0.5;
                            //     vy1d = v_y1 + hy1u*dv*0.5;
                            // }
                        }

                        // lower cell, up boundary
                        if(r < _nrRows-2 && !MV(r+2,c)) {
                            delta_h2 = h_yy2 - h_y2;
                            delta_u2 = u_yy2 - u_y2;
                            delta_v2 = v_yy2 - v_y2;
                            delta_h1 = tmph;
                            delta_u1 = tmpu;
                            delta_v1 = tmpv;
                            dh = limiter(delta_h1, delta_h2);
                            du = limiter(delta_u1, delta_u2);
                            dv = limiter(delta_v1, delta_v2);
                            hy2u = h_y2 - 0.5*dh;
                            uy2u = u_y2 - 0.5*du*hyd/H;
                            vy2u = v_y2 - 0.5*dv*hyd/H;
                            // double hy2d = h_y2 + 0.5*dh;
                            // if (H > minh) {
                            //     uy2u = u_y2 - 0.5*du*hy2u/H;
                            //     vy2u = v_y2 - 0.5*dv*hy2u/H;
                            // } else {
                            //     uy2u = u_y2 - hy2d*du*0.5;
                            //     vy2u = v_y2 - hy2d*dv*0.5;
                            // }
                        }
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
                    hll_y2 = F_Riemann(h_yd, vyd, uyd, h_y2u,vy2u,uy2u); // r and r+1

                    // determine smallest dt in x and y for each cell
                    double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                    double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]); // v[3] is max U and V in x and y

                    double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));
                    FloodDT->Drc = dt_req; // dt does not need to be a map, left over from earlier code
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
                            //double gflow_x = GRAV*0.5*( (H_l-H)*(H_l+H)+(H-H_r)*(H+H_r));// + delzcx*(H_l+H_r) ); delzcx = 0 when no muscl
                            //double gflow_y = GRAV*0.5*( (H_u-H)*(H_u+H)+(H-H_d)*(H+H_d));// + delzcy*(H_u+H_d) );
                            double gflow_x = GRAV*0.5*( (h_xl-hxl)*(h_xl+hxl) + (hxr-h_xr)*(hxr+h_xr) + delzcx*(hxl+hxr)); //delzcx = 0 when no muscl
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

double TWorld::fullSWOF2openWS(int nr_, cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
{
    /*
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    //double sumS = 0;
    bool stop;
    double dt_req_min = dt_max;
    int step = 0;

    if (startFlood)
    {

        sumh = getMassWS(nr_, h, 0);

        do {
            // bool SwitchLimitSWOFVelocity = true;
            //double vmax = 100000;
            // if (SwitchLimitSWOFVelocity)
            //      vmax = std::min(courant_factor, 0.2) * _dx/dt_req_min;
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_LWS(nr_) {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
                FloodDT->Drc = dt_max;
                tmb->Drc = 0;
            }}

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_LWS(nr_) {
                if (hs->Drc > F_minWH) {
                    tmb->Drc = 1;
                    if (c > 0 && !MV(r,c-1)        ) tmb->data[r][c-1] = 1;
                    if (c < _nrCols-1 && !MV(r,c+1)) tmb->data[r][c+1] = 1;
                    if (r > 0 && !MV(r-1,c)        ) tmb->data[r-1][c] = 1;
                    if (r < _nrRows-1 && !MV(r+1,c)) tmb->data[r+1][c] = 1;

                    if (c > 0 && r > 0 && !MV(r-1,c-1)                ) tmb->data[r-1][c-1]=1;
                    if (c < _nrCols-1 && r < _nrRows-1 && !MV(r+1,c+1)) tmb->data[r+1][c+1]=1;
                    if (r > 0 && c < _nrCols-1 && !MV(r-1,c+1)        ) tmb->data[r-1][c+1]=1;
                    if (c > 0 && r < _nrRows-1 && !MV(r+1,c-1)        ) tmb->data[r+1][c-1]=1;
                }
            }}

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (tmb->Drc == 1 && WaterSheds->Drc != (double) nr_)
                    tmb = 0;
            }}

            //do all flow and state calculations
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_LWS(nr_) {
            if (tmb->Drc > 0) {

                double dt = dt_req_min;
                double Un, Vn;
                //  double vmax = std::min(courant_factor, 0.2) * _dx/dt_req_min;

                vec4 hll_x1;
                vec4 hll_x2;
                vec4 hll_y1;
                vec4 hll_y2;

                double dx = _dx;//ChannelAdj->Drc;
                double dy = _dx;//DX->Drc;

                double H = hs->Drc;
                double n = N->Drc;
                double Z = z->Drc;
                double Vx = vxs->Drc;
                double Vy = vys->Drc;

                bool bc1 = c > 0 && !MV(r,c-1)        ;
                bool bc2 = c < _nrCols-1 && !MV(r,c+1);
                bool br1 = r > 0 && !MV(r-1,c)        ;
                bool br2 = r < _nrRows-1 && !MV(r+1,c);

                double z_x1 =  bc1 ? z->data[r][c-1] : Z;
                double z_x2 =  bc2 ? z->data[r][c+1] : Z;
                double z_y1 =  br1 ? z->data[r-1][c] : Z;
                double z_y2 =  br2 ? z->data[r+1][c] : Z;

                double h_x1 =  bc1 ? hs->data[r][c-1] : H;
                double h_x2 =  bc2 ? hs->data[r][c+1] : H;
                double h_y1 =  br1 ? hs->data[r-1][c] : H;
                double h_y2 =  br2 ? hs->data[r+1][c] : H;

                double u_x1 = bc1 ? vxs->data[r][c-1] : Vx;
                double u_x2 = bc2 ? vxs->data[r][c+1] : Vx;
                double u_y1 = br1 ? vxs->data[r-1][c] : Vx;
                double u_y2 = br2 ? vxs->data[r+1][c] : Vx;

                double v_x1 = bc1 ? vys->data[r][c-1] : Vy;
                double v_x2 = bc2 ? vys->data[r][c+1] : Vy;
                double v_y1 = br1 ? vys->data[r-1][c] : Vy;
                double v_y2 = br2 ? vys->data[r+1][c] : Vy;

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

                // z is blocking to prevent flow when water is flat and Z is not flat, described in article SWOF
                double h_x1r = std::max(0.0, h_x1 - std::max(0.0,  dz_x1 + fb_x1));
                double H_l   = std::max(0.0, H    - std::max(0.0, -dz_x1 + fb_x1));
                if(bc1)
                    hll_x1 = F_Riemann(h_x1r,u_x1,v_x1, H_l,Vx,Vy); // c-1 and c  //
                else
                    hll_x1 = F_Riemann(0,0,0, H_l,Vx,Vy);

                double H_r   = std::max(0.0, H    - std::max(0.0,  dz_x2 + fb_x2));
                double h_x2l = std::max(0.0, h_x2 - std::max(0.0, -dz_x2 + fb_x2));
                if(bc2)
                    hll_x2 = F_Riemann(H_r,Vx,Vy, h_x2l,u_x2,v_x2); // c and c+1
                else
                    hll_x2 = F_Riemann(H_r,Vx,Vy, 0,0,0);

                double h_y1d = std::max(0.0, h_y1 - std::max(0.0,  dz_y1 + fb_y1));
                double H_u   = std::max(0.0, H    - std::max(0.0, -dz_y1 + fb_y1));
                if (br1)
                    hll_y1 = F_Riemann(h_y1d,v_y1,u_y1, H_u,Vy,Vx); // r-1 and r
                else
                    hll_y1 = F_Riemann(0,0,0, H_u,Vy,Vx);

                double H_d   = std::max(0.0, H    - std::max(0.0,  dz_y2 + fb_y2));
                double h_y2u = std::max(0.0, h_y2 - std::max(0.0, -dz_y2 + fb_y2));
                if(br2)
                    hll_y2 = F_Riemann(H_d,Vy,Vx, h_y2u,v_y2,u_y2); // r and r+1
                else
                    hll_y2 = F_Riemann(H_d,Vy,Vx, 0,0,0);

                // determine smallest dt in x and y for each cell
                double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]);
                double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));
                FloodDT->Drc = dt_req;

                // if step = 0 do not calculate new fluxes and states yet because the first dt is always dt_max
                // find a smallest dt of the flow domain first
                if (step > 0) {

                    double tx = dt/dx;
                    double ty = dt/dy;

                    double flux_x1 = +hll_x1.v[0]/_dx;
                    double flux_x2 = -hll_x2.v[0]/_dx;
                    double flux_y1 = +hll_y1.v[0]/_dx;
                    double flux_y2 = -hll_y2.v[0]/_dx;

                    // if cell drops < 0 then adjust timestep
                    double tot = dt*(flux_x1 + flux_x2 + flux_y1 + flux_y2);
                    if (H+tot < 0) {
                        dt = H/-tot*dt;
                        // qDebug() << "oei" << H-tot;
                    }

                    double hn = std::max(0.0, H + dt*(flux_x1 + flux_x2 + flux_y1 + flux_y2));
                    // mass balance

                    // momentum balance for cells with water
                    if(hn > he_ca) {
                        // SWOF solution, delzc1 = 0 when not MUSCL
                        //  GRAV*0.5*((h1g_-h1l_)*(h1g_+h1l_) + (h1r_-h1d_)*(h1r_+h1d_) + (h1l_+h1r_)*delzc1->Drc));
                        double gflow_x = GRAV*0.5*( (H_l-H)*(H_l+H)+(H-H_r)*(H+H_r) );
                        double gflow_y = GRAV*0.5*( (H_u-H)*(H_u+H)+(H-H_d)*(H+H_d) );

                        double qxn = H * Vx - tx*(hll_x2.v[1] - hll_x1.v[1] + gflow_x) - ty*(hll_y2.v[2] - hll_y1.v[2]);
                        double qyn = H * Vy - tx*(hll_x2.v[2] - hll_x1.v[2]) - ty*(hll_y2.v[1] - hll_y1.v[1] + gflow_y);

                        double vsq = sqrt(Vx * Vx + Vy * Vy);
                        double nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01,pow(hn,4.0/3.0));
                        double nsq = nsq1*vsq*dt;

                        Un = (qxn/(1.0+nsq))/std::max(0.01,hn);
                        Vn = (qyn/(1.0+nsq))/std::max(0.01,hn);

                        if (SwitchTimeavgV) {
                            double fac = 0.5+0.5*std::min(1.0,4*hn)*std::min(1.0,4*hn);
                            fac = fac *exp(- std::max(1.0,dt) / nsq1);
                            Un = fac * Vx + (1.0-fac) *Un;
                            Vn = fac * Vy + (1.0-fac) *Vn;
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


                    h->Drc = hn;
                    vx->Drc = Un;
                    vy->Drc = Vn;
                } // step > 0
            } // tmb > 0, active cells + 1
            }}

            // find smallest domain dt
            #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
            FOR_ROW_COL_MV_LWS(nr_) {
                dt_req_min = std::min(dt_req_min, FloodDT->Drc);
            }}
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            if (step > 0) {

//                if (SwitchErosion) {
//                    SWOFSediment(dt_req_min, h,vx,vy);
//                }

                if (Switch2DDiagonalFlow) {
                    SWOFDiagonalFlow(dt_req_min, h, vx, vy);
                }


                timesum += dt_req_min;
                count++; // nr loops
            }

            step += 1; // now we have a good dt min, do the real calculations

            stop = timesum > _dt-0.001;
            if(count > F_MaxIter)
            stop = true;

        } while (!stop);

       correctMassBalanceWS(nr_, sumh, h, 0);

    } // if floodstart

    //qDebug() << _dt/count << count << dt_req_min;
    iter_n = std::max(1,count);
    return(count > 0 ? _dt/count : _dt);
    */
    return 0;
}


