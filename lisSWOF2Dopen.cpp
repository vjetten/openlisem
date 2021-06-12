

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
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-10
#define ve_ca 1e-10

#define GRAV 9.8067


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
//                        BLFlood->Drcr += dBL;
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
double TWorld::fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    double sumS = 0;
    bool stop;
    double dt_req_min = dt_max;
    int step = 0;

    if (startFlood)
    {

        sumh = getMass(h, 0);
//        if (SwitchErosion)
//            sumS = getMassSed(SSFlood, 0);


//        #pragma omp parallel for num_threads(userCores)
//        FOR_ROW_COL_MV_L {
//           FloodDT->Drc = dt_max;
//        }}
fill(*tmc,0);
        do {
            // bool SwitchLimitSWOFVelocity = true;
            //double vmax = 100000;
            // if (SwitchLimitSWOFVelocity)
            //      vmax = std::min(courant_factor, 0.2) * _dx/dt_req_min;
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
                FloodDT->Drc = dt_max;
                tmc->Drc = h->Drc*sqrt(vx->Drc*vx->Drc+vy->Drc*vy->Drc);
                tmb->Drc = 0;
            }}

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (hs->Drc > 0) {
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

            //do all flow and state calculations
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
            if (tmb->Drc > 0) {

                double dt = dt_req_min;
                double vxn, vyn;
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

                double vx_x1 = bc1 ? vxs->data[r][c-1] : Vx;
                double vx_x2 = bc2 ? vxs->data[r][c+1] : Vx;
                double vx_y1 = br1 ? vxs->data[r-1][c] : Vx;
                double vx_y2 = br2 ? vxs->data[r+1][c] : Vx;

                double vy_x1 = bc1 ? vys->data[r][c-1] : Vy;
                double vy_x2 = bc2 ? vys->data[r][c+1] : Vy;
                double vy_y1 = br1 ? vys->data[r-1][c] : Vy;
                double vy_y2 = br2 ? vys->data[r+1][c] : Vy;

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

                // calculate Riemann valaues for all four boundaries of a cell

                //coding left right and up/down boundary h
                //h_x1r|H_l  H  H_r|h_x2l
                //_____|___________|_____

                // |h_y1d
                // |-----
                // |H_u
                // |
                // |H
                // |
                // |H_d
                // |-----
                // |h_y2u

                // z is blocking to prevent flow when water is flat and Z is not flat, described in article SWOF
                double h_x1r = std::max(0.0, h_x1 - std::max(0.0,  dz_x1 + fb_x1));
                double H_l   = std::max(0.0, H    - std::max(0.0, -dz_x1 + fb_x1));
                if(bc1)
                    hll_x1 = F_Riemann(h_x1r,vx_x1,vy_x1, H_l,Vx,Vy); // c-1 and c  //
                else
                    hll_x1 = F_Riemann(0,0,0, H_l,Vx,Vy);

                double H_r   = std::max(0.0, H    - std::max(0.0,  dz_x2 + fb_x2));
                double h_x2l = std::max(0.0, h_x2 - std::max(0.0, -dz_x2 + fb_x2));
                if(bc2)
                    hll_x2 = F_Riemann(H_r,Vx,Vy, h_x2l,vx_x2,vy_x2); // c and c+1
                else
                    hll_x2 = F_Riemann(H_r,Vx,Vy, 0,0,0);

                double h_y1d = std::max(0.0, h_y1 - std::max(0.0,  dz_y1 + fb_y1));
                double H_u   = std::max(0.0, H    - std::max(0.0, -dz_y1 + fb_y1));
                if (br1)
                    hll_y1 = F_Riemann(h_y1d,vy_y1,vx_y1, H_u,Vy,Vx); // r-1 and r
                else
                    hll_y1 = F_Riemann(0,0,0, H_u,Vy,Vx);

                double H_d   = std::max(0.0, H    - std::max(0.0,  dz_y2 + fb_y2));
                double h_y2u = std::max(0.0, h_y2 - std::max(0.0, -dz_y2 + fb_y2));
                if(br2)
                    hll_y2 = F_Riemann(H_d,Vy,Vx, h_y2u,vy_y2,vx_y2); // r and r+1
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

                    // limiting flux determines everything NOT NECESSARY?
                    //                        double C = 1.0;
                    //                        flux_x1 = std::max(-H * C,std::min(flux_x1,h_x1 * C));
                    //                        flux_x2 = std::max(-H * C,std::min(flux_x2,h_x2 * C));
                    //                        flux_y1 = std::max(-H * C,std::min(flux_y1,h_y1 * C));
                    //                        flux_y2 = std::max(-H * C,std::min(flux_y2,h_y2 * C));

                    //                        double factor_flowx1f = 1.0-std::min(1.0,std::max(0.0, dz_x1)/std::max(1e-6,h_x1));
                    //                        double factor_flowy1f = 1.0-std::min(1.0,std::max(0.0, dz_y1)/std::max(1e-6,h_y1));
                    //                        double factor_flowx2f = 1.0-std::min(1.0,std::max(0.0,-dz_x2)/std::max(1e-6,h_x2));
                    //                        double factor_flowy2f = 1.0-std::min(1.0,std::max(0.0,-dz_y2)/std::max(1e-6,h_y2));

                    //                        double factor_flowx1t = 1.0-std::min(1.0,std::max(0.0,-dz_x1)/std::max(1e-6,H));
                    //                        double factor_flowy1t = 1.0-std::min(1.0,std::max(0.0,-dz_y1)/std::max(1e-6,H));
                    //                        double factor_flowx2t = 1.0-std::min(1.0,std::max(0.0, dz_x2)/std::max(1e-6,H));
                    //                        double factor_flowy2t = 1.0-std::min(1.0,std::max(0.0, dz_y2)/std::max(1e-6,H));

                    //                        flux_x1 = std::max(-H * factor_flowx1t * C, std::min(flux_x1, h_x1 * factor_flowx1f * C));
                    //                        flux_x2 = std::max(-H * factor_flowx2t * C, std::min(flux_x2, h_x2 * factor_flowy1f * C));
                    //                        flux_y1 = std::max(-H * factor_flowy1t * C, std::min(flux_y1, h_y1 * factor_flowx2f * C));
                    //                        flux_y2 = std::max(-H * factor_flowy2t * C, std::min(flux_y2, h_y2 * factor_flowy2f * C));

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

                        vxn = (qxn/(1.0+nsq))/std::max(0.01,hn);
                        vyn = (qyn/(1.0+nsq))/std::max(0.01,hn);

                        if (SwitchTimeavgV) {
                            double fac = 0.5+0.5*std::min(1.0,4*hn)*std::min(1.0,4*hn);
                            fac = fac *exp(- std::max(1.0,dt) / nsq1);
                            vxn = fac * Vx + (1.0-fac) *vxn;
                            vyn = fac * Vy + (1.0-fac) *vyn;
                        }

                    } else { // hn < ha
                        hn = H; // if no fluxes then also no change in h
                        vxn = 0;
                        vyn = 0;
                    }

                    // dan maar even met geweld!
                    if (std::isnan(vxn) || std::isnan(vyn)  )
                    {
                        vxn = 0;
                        vyn = 0;
                    }

                    if (fabs(vxn) <= ve_ca)
                        vxn = 0;
                    if (fabs(vyn) <= ve_ca)
                        vyn = 0;

                    if (fabs(vxn) > 20 && hn < 0.01)
                        vxn = (vxn < 0 ? -1.0 : 1.0) * pow(hn,2.0/3.0)*sqrtGrad->Drc/N->Drc;
                    if (fabs(vyn) > 20 && hn < 0.01)
                        vyn = (vyn < 0 ? -1.0 : 1.0) * pow(hn,2.0/3.0)*sqrtGrad->Drc/N->Drc;


                    h->Drc = hn;
                    vx->Drc = vxn;
                    vy->Drc = vyn;
                } // step > 0
            } // tmb > 0, active cells + 1
            }}

            // find smallest domain dt
            #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                dt_req_min = std::min(dt_req_min, FloodDT->Drc);
            }}
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            if (step > 0) {

                if (SwitchErosion) {
                    SWOFSediment(dt_req_min, h,vx,vy);
                }

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

        correctMassBalance(sumh, h, 0);

//        if (SwitchErosion)
  //          correctMassBalanceSed(sumS, SSFlood, 0);

        //            double sumh1 = getMass(h, 0);
        //            qDebug() << sumh << sumh1 << (sumh-sumh1)/sumh;
//        if (SwitchErosion && !SwitchErosionInsideLoop) {
//            SWOFSediment(_dt, h,vx,vy);
//        }
    } // if floodstart

    //qDebug() << _dt/count << count << dt_req_min;
    iter_n = std::max(1,count);
    return(count > 0 ? _dt/count : _dt);
}

