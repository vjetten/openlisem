

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


double TWorld::fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;
    int step = 0;

    if (startFlood)
    {
        sumh = getMass(h, 0);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            FloodDT->Drc = dt_max;
            //FloodT->Drc = 0;
          //  tma->Drc = 0;
        }}

        do {

           // bool SwitchLimitSWOFVelocity = true;
            double vmax = 100000;
           // if (SwitchLimitSWOFVelocity)
          //      vmax = std::min(courant_factor, 0.2) * _dx/dt_req_min;

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                hs->Drc = h->Drc;
                vxs->Drc = std::max(-vmax, std::min(vmax,vx->Drc));
                vys->Drc = std::max(-vmax, std::min(vmax,vy->Drc));
                //limit V here, than not necessary later
                tmb->Drc = 0;
            }}
            // tmb is used as flag for cells that need processing

            // set tmb for all cells with surface water plus 1 surrounding cell
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
               if (tmb->Drc > 0)
               {
                    double dt = dt_req_min;
                    double vxn, vyn;
                    //  double vmax = std::min(courant_factor, 0.2) * _dx/dt_req_min;

                    //typedef struct vec4 { double v[4]; } vec4;
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

                        // limiting flux determines everything
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

                            // bastian

//                            double threshold = 0.01 * dx;
//                            if(hn < threshold)
//                            {
//                                double B = 0.5; //1.0 is theoretical max else faster than gravity
//                                double sx_zh_x1 = std::min(B, std::max(-B, (Z + H - z_x1 - h_x1)/dx));
//                                double sx_zh_x2 = std::min(B, std::max(-B, (z_x2 + h_x2 - Z - H)/dx));
//                                double sy_zh_y1 = std::min(B, std::max(-B, (Z + H - z_y1 - h_y1)/dy));
//                                double sy_zh_y2 = std::min(B, std::max(-B, (z_y2 + h_y2 - Z - H)/dy));

//                                // if B = 0.5 this can never be >1?
//                                double sx_zh = std::min(1.0,std::max(-1.0,limiter(sx_zh_x1, sx_zh_x2)));
//                                double sy_zh = std::min(1.0,std::max(-1.0,limiter(sy_zh_y1, sy_zh_y2)));

//                                double kinfac = std::max(0.0,(threshold - hn) / (0.025 * dx));
//                                double acc_eff = (vxn - Vx)/std::max(0.0001,dt);

//                                double v_kin = (sx_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sx_zh>0?sx_zh:-sx_zh))/(0.001+n);

//                                vxn = kinfac * v_kin + vxn*(1.0-kinfac);

//                                kinfac = std::max(0.0,(threshold - hn) / (0.025 * dx));
//                                acc_eff = (vyn - Vy)/std::max(0.0001,dt);

//                                v_kin = (sy_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sy_zh>0?sy_zh:-sy_zh))/(0.001+n);

//                                vyn = kinfac * v_kin + vyn*(1.0f-kinfac);

//                            }
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

                        h->Drc = hn;
                        vx->Drc = vxn;
                        vy->Drc = vyn;
                    } // step > 0
                } // tmb > 0, active cells + 1
            }}

            // find smallest domain dt
            #pragma omp parallel for reduction(min:dt_req_min) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                double res = FloodDT->Drc;
                dt_req_min = std::min(dt_req_min, res);
            }}
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            // put it back in FloodDT
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                FloodDT->Drc = dt_req_min;
            }}

            if (step > 0) {
                SWOFDiagonalFlow(dt_req_min, h, vx, vy);

                if (SwitchErosion && SwitchErosionInsideLoop) {
                    SWOFSediment(dt_req_min, FloodDT,h,vx,vy);
                }
                timesum += dt_req_min;
                count++; // nr loops
            }

            step = 1; // now we have a good dt min, do the real calculations

            stop = timesum > _dt-0.001;
            if(count > F_MaxIter)
                stop = true;

        } while (!stop);

        correctMassBalance(sumh, h, 0);
//            double sumh1 = getMass(h, 0);
//            qDebug() << sumh << sumh1 << (sumh-sumh1)/sumh;
        if (SwitchErosion && !SwitchErosionInsideLoop) {
            SWOFSediment(_dt, FloodDT,h,vx,vy);
        }
    } // if floodstart

    //qDebug() << _dt/count << count << dt_req_min;
    iter_n = std::max(1,count);
    return(count > 0 ? _dt/count : _dt);
}

//-----------------------------------------------------------------------------------------------------
void TWorld::makeChannelList()
{
    /*
    if(!SwitchIncludeChannel)
        return;

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    fill(*tma, -1);

    for (int rr = 0; rr < _nrRows; rr++)
        for (int cr = 0; cr < _nrCols; cr++) {
            if(LDDChannel->Drcr == 5) {
                //LDD_LINKEDLIST *chlist = nullptr;
                LDD_LINKEDLIST *temp = nullptr;
                chlist = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                chlist->prev = nullptr;
                chlist->rowNr = rr;
                chlist->colNr = cr;

                while (chlist != nullptr)
                {
                    int i = 0;
                    bool  subCatchDone = true;

                    int rowNr = chlist->rowNr;
                    int colNr = chlist->colNr;


                    for (i=1; i<=9; i++)
                    {
                        int r, c;
                        int ldd = 0;

                        // this is the current cell
                        if (i==5)
                            continue;

                        r = rowNr+dy[i];
                        c = colNr+dx[i];

                        if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                            ldd = (int) LDDChannel->Drc;

                        // check if there are more cells upstream, if not subCatchDone remains true
                        if (tma->Drc < 0 && ldd > 0 && FLOWS_TO(ldd, r, c, rowNr, colNr)) {
                            temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                            temp->prev = chlist;
                            chlist = temp;
                            chlist->rowNr = r;
                            chlist->colNr = c;
                            subCatchDone = false;
                        }
                    }

                    if (subCatchDone)
                    {
                        tma->data[rowNr][colNr] = 1; // flag done

                        temp=chlist;
                        chlist=chlist->prev;
                        //free(temp);
                    }
                }
            }
        }
    }
    */
}

void TWorld::ChannelSWOFopen()
{
    if(!SwitchIncludeChannel)
        return;

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    bool CorrectMassBalance = true;
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    bool stop;
    double dt = dt_max;
    double C = 0.25;//std::min(0.25, courant_factor);
    double B = 1.0;
    double dt_req = dt_max;

    double sum1 = 0;
    double qout = 0;
    if (CorrectMassBalance) {
        FOR_ROW_COL_MV_CH {
            if(ChannelWH->Drc > 0)
                sum1 += ChannelWH->Drc;//*ChannelWidth->Drc*ChannelDX->Drc;
        }
    }


    fill(*ChannelQn, 0);
    fill(*tmb, 0);

    do {
        stop = false;

        // do the whole channel
        fill(*tma, -1);
        dt_req = dt_max;

        for (int rr = 0; rr < _nrRows; rr++)
            for (int cr = 0; cr < _nrCols; cr++) {
                if(LDDChannel->Drcr == 5) {

                    LDD_LINKEDLIST *chlist = nullptr;
                    LDD_LINKEDLIST *temp = nullptr;
                    chlist = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                    chlist->prev = nullptr;
                    chlist->rowNr = rr;
                    chlist->colNr = cr;

                    while (chlist != nullptr)
                    {
                        int i = 0;
                        bool  subCatchDone = true;
                        int rowNr = chlist->rowNr;
                        int colNr = chlist->colNr;

                        for (i=1; i<=9; i++)
                        {
                            int r, c;
                            int ldd = 0;

                            // this is the current cell
                            if (i==5)
                                continue;

                            r = rowNr+dy[i];
                            c = colNr+dx[i];

                            if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                                ldd = (int) LDDChannel->Drc;

                            // check if there are more cells upstream, if not subCatchDone remains true
                            if (tma->Drc < 0 && ldd > 0 && FLOWS_TO(ldd, r, c, rowNr, colNr)) {
                                temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                                temp->prev = chlist;
                                chlist = temp;
                                chlist->rowNr = r;
                                chlist->colNr = c;
                                subCatchDone = false;
                            }
                        }

                        if (subCatchDone)
                        {
                            double Dx = ChannelDX->data[rowNr][colNr];
                            double tx = dt/Dx;

                            // ===== the current cell =====
                            double H = ChannelWH->data[rowNr][colNr];
                            double V = ChannelU->data[rowNr][colNr];
                            double Z = DEM->data[rowNr][colNr];
                            double W = ChannelWidth->data[rowNr][colNr];
                            double N = ChannelN->data[rowNr][colNr];
                          //  double Vol = W*Dx*H;

                            // new H and V
                            double Hn = H;
                            double Vn = V;

                            // ===== outflow to downstream cell =====
                            double ch_vadd = 0; //gravity component
                            double flux_out = 0; // flux to downstream cell as volume

                            vec4 hll_out = {0,0,0,0};
                            int ldd = (int)LDDChannel->data[rowNr][colNr];

                            // if ldd == 5 use these values
                            double Vo = V;//pow(H, 2.0/3.0)/N*sqrt(H/_dx+0.001);
                            double Ho = H;//std::max(0.0, H*(1-Vo/Dx*dt));
                            double Zo = Z;//*0.99;//std::max(0.99, (1.0-ChannelGrad->data[rowNr][colNr])); // at least a 1% slope
                            double Wo = W;


//                            float ch_slope = (ch_h)/dx;
//                            float ch_q = max(0.0f,dt * ch_h * ch_width * sqrt((float)(0.001+ch_h))/(dx * (0.05f)));
//                            ch_q = max(0.0f,min(0.25f * ch_vol,ch_q));
//                            chhn = chhn - 0.5 * ch_q/(ch_width * dx);
//                            flux_chx2 = flux_chx2 + 0.5 * ch_q;
//                            ch_vadd = ch_vadd + dt * 0.5 * GRAV * max(-1.0f,min(1.0f,(float)(ch_slope)));
//                            qfout = qfout + 0.5 * ch_q;

                            // if not pit get the downstream values
                            if (ldd != 5) {
                                int r = rowNr+dy[ldd];
                                int c = colNr+dx[ldd];
                                Ho = ChannelWH->Drc;
                                Vo = ChannelU->Drc;
                                Zo = DEM->Drc;
                                Wo = ChannelWidth->Drc;
                            }

//                            float3 hll_x1 = F_HLL2(ch_h,ch_v,0,chn_h,chn_v,0);
//                            float ch_q = (dt/dx)*(max(ch_width,chn_width)/dx)*((dx * 0.5*(chn_width +ch_width)) *hll_x1.x);
//                            ch_q = min(0.25f * ch_vol,ch_q);
//                            ch_q = max(-0.25f * chn_vol,ch_q);
//                            ch_q = ch_q * 0.5;
//                            float ch_slope = (z + ch_h - chn_z - chn_h)/dx;
//                            ch_vadd = ch_vadd + dt * 0.5 * GRAV * max(-1.0f,min(1.0f,(float)(ch_slope)));
//                            if(ch_q < 0)
//                            {
//                                    float new_ch_vol = chhn*(ch_width*dx);
//                                    chvn = (chvn * new_ch_vol - chn_v *(ch_q))/max(0.01f,new_ch_vol - ch_q);
//                            }
//                            chhn = chhn - ch_q/(ch_width * dx);
//                            flux_chx2 = flux_chx2 + ch_q;

                            hll_out = F_Riemann(H,V,0, Ho,Vo,0);
                            // 1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
                            double Q = tx * hll_out.v[0] * (std::max(W,Wo)/Dx * (Dx*0.5*(W+Wo)));
                            // s/m * m2/s * relatieve breedte (waarom max) * channel oppervlakte = volume
                            double Volo = Wo*Dx*Ho;
                            double Vol = W*Dx*H;
                            Q = std::max(-C*Volo, std::min(Q, C*Vol));
                            if(Q < 0) {
                                Vn = (Vn*Vol - Vo*Q)/std::max(0.01,Vol - Q);
                            }
                            double s_zh_out = std::min(B, std::max(-B, (H + Z - Zo - Ho)/Dx));
                            ch_vadd = ch_vadd + dt * 0.5 * GRAV * s_zh_out;
                            // gravity pressure deel


                            Hn = Hn - Q/(W*Dx);
                            flux_out = Q; //m3? must be m3/s

                            // ===== weighed sum inflow from upstream cells =====
                            double flux_in = 0;
                            double ch_vaddw = 0.5;
                            for (i = 1; i <= 9; i++)
                            {
                                int r, c, ldd = 0;

                                if (i==5)
                                    continue;

                                r = rowNr+dy[i];
                                c = colNr+dx[i];

                                if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                                    ldd = (int) LDDChannel->Drc;
                                else
                                    continue;

                                // if the cells flows into the c
                                if (ldd > 0 && ldd != 5 && FLOWS_TO(ldd, r,c,rowNr,colNr)){
                                    double Hi = ChannelWH->Drc;
                                    double Vi = ChannelU->Drc;
                                    double Zi = DEM->Drc;
                                    double Wi = ChannelWidth->Drc;

                                    vec4 hll_in = F_Riemann(Hi,Vi,0, H,V,0);
                                    double Q = tx * hll_in.v[0] * (std::max(W,Wi)/Dx * (Dx*0.5*(W+Wi)) );
                                    double Voli = Wi*Dx*Hi;
                                    double Vol = H*W*Dx;
                                    Q = std::max(-C*Vol,std::min(C*Voli,Q));
                                    if(Q > 0) {
                                        double Voln = Hn*W*Dx;
                                        Vn = (Vn*Voln + Vi*Q)/std::max(0.01,Voln + Q);
                                    }

                                    // gravity + pressure part
                                    double s_zh_in = std::min(B, std::max(-B, (Hi + Zi - Z - H)/Dx));
                                    ch_vadd = ch_vadd + dt * 0.5 * GRAV * s_zh_in;
                                    ch_vaddw = ch_vaddw + 0.5 * Wi/W;

                                    Hn = Hn + Q/(W * Dx);

                                    flux_in = flux_in + Q;
                                }
                            }

                            if(ch_vaddw > 1) {
                                ch_vadd = ch_vadd/ch_vaddw;
                            }

                            Hn = std::max(Hn, 0.0);
                            if (Hn > he_ca) {
                                Vn = Vn + ch_vadd;
                                double qv = sqrt(Vn*Vn);
                                double chnsq1 = (0.001+N)*(0.001+N)*GRAV/pow(Hn,4.0/3.0);
                                double chnsq = chnsq1*qv*dt;
                                Vn = (qv/(1.0+chnsq));
                                Vn = std::min(25.0,std::max(-25.0,Vn));

                                //                                if (SwitchTimeavgV) {
                                //                                    double fac = 0.5+0.5*std::min(1.0,4*Hn)*std::min(1.0,4*Hn);
                                //                                    fac = fac *exp(- std::max(1.0,dt) / chnsq1);
                                //                                    Vn = fac * V + (1.0-fac) *Vn;
                                //                                }
                            } else {
                                Vn = 0;
                                Hn = 0;
                            }

                            if (fabs(Vn) <= ve_ca)
                                Vn = 0;

                            dt_req = std::min(dt_req,courant_factor *Dx/( std::min(dt_max,std::max(0.01,fabs(Vn)))));
                            //std::max(TimestepfloodMin,

                            // gebruik riemann solver cfl
                            //   double dtx = Dx/hll_out.v[3];
                            //   dt_req = std::max(TimestepfloodMin, std::min(dt_req, courant_factor*dtx));


                            ChannelU->data[rowNr][colNr] = Vn;
                            ChannelWH->data[rowNr][colNr] = Hn;
                         //   ChannelQn->data[rowNr][colNr] = flux_out;
                            tmb->data[rowNr][colNr] += flux_out;

                            tma->data[rowNr][colNr] = 1; // flag done

                            temp=chlist;
                            chlist=chlist->prev;
                            free(temp);
                        }
                    }
                }
            }

        //  qDebug() << count << dt_req;
        dt = std::min(dt_req, _dt-timesum);
        //  qDebug() << timesum << dt;

        timesum = timesum + dt;
        if (timesum > _dt-1e-6)
            stop = true;
        count++;
        if (count > 200)
            stop = true;

    } while (!stop);

    //qDebug() << count;

    if (CorrectMassBalance) {
        double sum2 = 0;
        double n = 0;

        FOR_ROW_COL_MV_CH {
            if(ChannelWH->Drc > 0) {
                sum2 += ChannelWH->Drc;//*ChannelWidth->Drc*ChannelDX->Drc;
                if(ChannelWH->Drc > 0)
                    n += 1.0;
            }
        }
        double dhtot = sum2 > 0 ? ((sum1-qout) - sum2)/sum2 : 0;
        double sum3 = 0;
        qDebug() << sum1-sum2 << sum1-sum2-qout;

        FOR_ROW_COL_MV_CH {
            if(ChannelWH->Drc > 0) {
                ChannelWH->Drc = ChannelWH->Drc*(1.0 + dhtot);
                ChannelWH->Drc = std::max(ChannelWH->Drc , 0.0);
            }
            sum3 = sum3 + ChannelWH->Drc;
        }
        //    qDebug() << sum1-sum2 << sum1-sum3;
    }

    //#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_CH {
        ChannelV->Drc = fabs(ChannelU->Drc);
        ChannelQn->Drc = tmb->Drc/_dt;//ChannelV->Drc*ChannelWH->Drc*ChannelWidth->Drc;
        ChannelQ->Drc = ChannelQn->Drc;
        ChannelWaterVol->Drc = ChannelWH->Drc*ChannelWidth->Drc*ChannelDX->Drc;
        ChannelAlpha->Drc = ChannelWH->Drc*ChannelWidth->Drc/std::pow(ChannelQn->Drc, 0.6);
    }
    //report(*ChannelU,"chu");
}




