
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-10
#define ve_ca 1e-10

#define GRAV 9.8067
#define EPSILON 1e-10

// 2nd order without iteration dt1, dt2!
double TWorld::fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;
    int cnt;

    if (startFlood)
    {
        sumh = getMass(h);

#pragma omp parallel for collapse(2) num_threads(userCores)
        FOR_ROW_COL_MV_L {
            FloodDT->Drc = dt_max;
            FloodT->Drc = 0;
            tma->Drc = 0;
        }

        do {
            // make a copy
            double vmax = std::min(courant_factor, 0.2) * _dx/dt_req_min;
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                hs->Drc = h->Drc;
                vxs->Drc = std::max(-vmax, std::min(vmax,vx->Drc));
                vys->Drc = std::max(-vmax, std::min(vmax,vy->Drc));
                //limit V here, than not necessary later
            }

            // tmb is used as flag for cells that need processing
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tmb->Drc = 0;
            }

#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (hs->Drc > 0) {
                    tmb->Drc = 1;
                    if (c > 0 && !MV(r,c-1)        ) tmb->data[r][c-1] = 1;
                    if (c < _nrCols-1 && !MV(r,c+1)) tmb->data[r][c+1] = 1;
                    if (r > 0 && !MV(r-1,c)        ) tmb->data[r-1][c] = 1;
                    if (r < _nrRows-1 && !MV(r+1,c)) tmb->data[r+1][c] = 1;
                }
            }

            if (SwitchVariableTimestep) {
                double dt = dt_max;
                //#pragma omp parallel for reduction(min:dt) collapse(2) num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    dt = FloodDT->Drc;
                    if (c > 0 && !MV(r,c-1)        ) dt = std::min(dt, FloodDT->data[r][c-1]);
                    if (c < _nrCols-1 && !MV(r,c+1)) dt = std::min(dt, FloodDT->data[r][c+1]);
                    if (r > 0 && !MV(r-1,c)        ) dt = std::min(dt, FloodDT->data[r-1][c]);
                    if (r < _nrRows-1 && !MV(r+1,c)) dt = std::min(dt, FloodDT->data[r+1][c]);
                    double fc = 2.0;
                    if (c > 0 && r > 0 && !MV(r-1,c-1)                ) dt = std::min(dt, fc*FloodDT->data[r-1][c-1]);
                    if (c < _nrCols-1 && r < _nrRows-1 && !MV(r+1,c+1)) dt = std::min(dt, fc*FloodDT->data[r+1][c+1]);
                    if (r > 0 && c < _nrCols-1 && !MV(r-1,c+1)        ) dt = std::min(dt, fc*FloodDT->data[r-1][c+1]);
                    if (c > 0 && r < _nrRows-1 && !MV(r+1,c-1)        ) dt = std::min(dt, fc*FloodDT->data[r+1][c-1]);
                    fc = 1.5;
                    if (c > 1 && !MV(r,c-2)        ) dt = std::min(dt, fc*FloodDT->data[r][c-2]);
                    if (c < _nrCols-2 && !MV(r,c+2)) dt = std::min(dt, fc*FloodDT->data[r][c+2]);
                    if (r > 1 && !MV(r-2,c)        ) dt = std::min(dt, fc*FloodDT->data[r-2][c]);
                    if (r < _nrRows-2 && !MV(r+2,c)) dt = std::min(dt, fc*FloodDT->data[r+2][c]);

                    FloodDT->Drc = dt;

                    if (FloodT->Drc > _dt - 0.001)
                        tmb->Drc = 0;
                }
            }

            //flow for cells which have h and not done yet (FloodT < _dt)
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (tmb->Drc > 0) {
                    double dt = SwitchVariableTimestep ? FloodDT->Drc : dt_req_min;
                    double vxn, vyn;
                    //   double vmax = std::min(courant_factor, 0.2) * _dx/dt_req_min;

                    tma->Drc += 1.0; // nr times a cell is processed

                    //typedef struct vec4 { double v[4]; } vec4;
                    vec4 hll_x1;
                    vec4 hll_x2;
                    vec4 hll_y1;
                    vec4 hll_y2;

                    double dx = ChannelAdj->Drc;
                    double dy = DX->Drc;

                    double H = hs->Drc;
                    double n = N->Drc;
                    double Z = z->Drc;
                    double Vx = vxs->Drc;//std::max(-vmax, std::min(vmax, vxs->Drc));
                    double Vy = vys->Drc;//std::max(-vmax, std::min(vmax, vys->Drc));

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

                    double fac = DEMdz->Drc; // if Z is in a pit > 10m from the surrounding cells, reduce the effect of the DEM
                    double dz_x1 = fac*(Z - z_x1);
                    double dz_x2 = fac*(z_x2 - Z);
                    double dz_y1 = fac*(Z - z_y1);
                    double dz_y2 = fac*(z_y2 - Z);

                    double h_x1l = std::max(0.0, h_x1 - std::max(0.0,  dz_x1 + fb_x1));
                    double h_x1r = std::max(0.0, H    - std::max(0.0, -dz_x1 + fb_x1));
                    if(bc1)
                        hll_x1 = F_Riemann(h_x1l,vx_x1,vy_x1,h_x1r,Vx,Vy); // c-1 and c  //
                    else
                        hll_x1 = F_Riemann(0,0,0,h_x1r,Vx,Vy);

                    double h_x2l = std::max(0.0, H    - std::max(0.0,  dz_x2 + fb_x2));
                    double h_x2r = std::max(0.0, h_x2 - std::max(0.0, -dz_x2 + fb_x2));
                    if(bc2)
                        hll_x2 = F_Riemann(h_x2l,Vx,Vy,h_x2r,vx_x2,vy_x2); // c and c+1
                    else
                        hll_x2 = F_Riemann(h_x2l,Vx,Vy,0,0,0);

                    double h_y1u = std::max(0.0, h_y1 - std::max(0.0,  dz_y1 + fb_y1));
                    double h_y1d = std::max(0.0, H    - std::max(0.0, -dz_y1 + fb_y1));
                    if (br1)
                        hll_y1 = F_Riemann(h_y1u,vy_y1,vx_y1,h_y1d,Vy,Vx); // r-1 and r
                    else
                        hll_y1 = F_Riemann(0,0,0,h_y1d,Vy,Vx);

                    double h_y2u = std::max(0.0, H    - std::max(0.0,  dz_y2 + fb_y2));
                    double h_y2d = std::max(0.0, h_y2 - std::max(0.0, -dz_y2 + fb_y2));
                    if(br2)
                        hll_y2 = F_Riemann(h_y2u,Vy,Vx,h_y2d,vy_y2,vx_y2); // r and r+1
                    else
                        hll_y2 = F_Riemann(h_y2u,Vy,Vx,0,0,0);

                    double B = 0.5; //1.0 is theoretical max else faster than gravity
                    double sx_zh_x2 = std::min(B, std::max(-B, (z_x2 + h_x2 - Z - H)/dx));
                    double sy_zh_y1 = std::min(B, std::max(-B, (Z + H - z_y1 - h_y1)/dy));
                    double sx_zh_x1 = std::min(B, std::max(-B, (Z + H - z_x1 - h_x1)/dx));
                    double sy_zh_y2 = std::min(B, std::max(-B, (z_y2 + h_y2 - Z - H)/dy));

                    // if B = 0.5 this can never be >1?
                    double sx_zh = std::min(1.0,std::max(-1.0,limiter(sx_zh_x1, sx_zh_x2)));
                    double sy_zh = std::min(1.0,std::max(-1.0,limiter(sy_zh_y1, sy_zh_y2)));

                    double tx = dt/dx;
                    double ty = dt/dy;

                    double C = std::min(0.25, courant_factor);
                    double flux_x1 = std::max(-H * C,std::min(+tx*hll_x1.v[0],h_x1 * C));
                    double flux_x2 = std::max(-H * C,std::min(-tx*hll_x2.v[0],h_x2 * C));
                    double flux_y1 = std::max(-H * C,std::min(+ty*hll_y1.v[0],h_y1 * C));
                    double flux_y2 = std::max(-H * C,std::min(-ty*hll_y2.v[0],h_y2 * C));

                    double hn = std::max(0.0, H + flux_x1 + flux_x2 + flux_y1 + flux_y2);

                    if(hn > he_ca) {

                        double qxn = H * Vx - tx*(hll_x2.v[1] - hll_x1.v[1]) - ty*(hll_y2.v[2] - hll_y1.v[2])- 0.5 * GRAV *hn*sx_zh * dt;
                        double qyn = H * Vy - tx*(hll_x2.v[2] - hll_x1.v[2]) - ty*(hll_y2.v[1] - hll_y1.v[1])- 0.5 * GRAV *hn*sy_zh * dt;

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

                        double threshold = 0.01 * _dx; // was 0.01
                        if(hn < threshold)
                        {
                            double h23 = pow(hn, 2.0/3.0);//hn * sqrt(hn)
                            double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
                            double v_kin = (sx_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sx_zh > 0 ? sx_zh : -sx_zh))/(0.001+n);
                            vxn = kinfac * v_kin + vxn*(1.0-kinfac);
                            v_kin = (sy_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sy_zh > 0 ? sy_zh : -sy_zh))/(0.001+n);
                            vyn = kinfac * v_kin + vyn*(1.0-kinfac);
                        }

                    } else {
                        vxn = 0;
                        vyn = 0;
                        hn = 0;
                    }
                    // dan maar even met geweld!
                    if (std::isnan(vxn) || std::isnan(vyn)  )
                    {
                        vxn = 0;
                        vyn = 0;
                        hn= 0;
                    }

                    if (fabs(vxn) <= ve_ca)
                        vxn = 0;
                    if (fabs(vyn) <= ve_ca)
                        vyn = 0;
                    if (hn < he_ca) {
                        vxn = 0;
                        vyn = 0;
                        hn = 0;
                    }

                    double dt_req1 = courant_factor *_dx/( std::min(dt_max,std::max(0.01,sqrt(vxn*vxn + vyn*vyn))));
                    // gebruik riemann solver cfl
                    double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                    double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]);
                    double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));

                    FloodDT->Drc = std::min(dt_req1, dt_req);
                    // taking the smalklest works best for instabiliies!
                    h->Drc = hn;
                    vx->Drc = vxn;
                    vy->Drc = vyn;
                }
            }

#pragma omp parallel for reduction(min:dt_req_min) collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                double res = FloodDT->Drc;
                dt_req_min = std::min(dt_req_min, res);
            }
            dt_req_min = std::min(dt_req_min, _dt-timesum);
            timesum += dt_req_min;

#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (!SwitchVariableTimestep)
                    FloodDT->Drc = dt_req_min;
                else
                    FloodDT->Drc = FloodDT->Drc = std::max(TimestepfloodMin, std::min(FloodDT->Drc, _dt-FloodT->Drc));

                FloodT->Drc += FloodDT->Drc;
                if (FloodT->Drc > _dt)
                    FloodT->Drc = _dt;
            }

            if (SwitchErosion)
                SWOFSediment(FloodDT,hs,vxs,vys);

            if (SwitchVariableTimestep) {
                cnt = 0;
                // nr cells that need processing
                //#pragma omp parallel for reduction(+:cnt) collapse(2) num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    if (FloodT->Drc < _dt-0.001)
                        cnt++;
                }
                //qDebug() << cnt;
                stop = cnt < 1;
            } else {
                stop = timesum > _dt-0.001;
            }

            count++; // nr loops

            if(count > F_MaxIter)
                stop = true;
        } while (!stop);

        correctMassBalance(sumh, h);
    } // if floodstart

    // for screen output
    double avgdt = 0;
    double nc= 0;
    //#pragma omp parallel for reduction(+:avgdt) collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        FloodDT->Drc = tma->Drc > 0 ? _dt/tma->Drc : 0;
        FloodT->Drc = tma->Drc;
        if (tma->Drc > 0)
            nc += 1.0;
        avgdt = avgdt + FloodDT->Drc;
    }
    avgdt = avgdt/nc;
    iter_n = count;

    return(avgdt);//_dt/count
}


void TWorld::ChannelSWOFopen()
{
    if(!SwitchIncludeChannel)
        return;

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};


    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    bool stop;
    double dt = dt_max;
    //double C = std::min(0.25, courant_factor);
    double B = 0.5;
    double dt_req = dt_max;

    double sum1 = 0;
    #pragma omp parallel for reduction(+:sum1) collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_CH {
        if(ChannelWH->Drc > 0)
            sum1 += ChannelWH->Drc;
    }

    fill(*ChannelQn, 0);

    do {
        stop = false;

        // do the whole channel
        fill(*tma, -1);
        dt_req = dt_max;

        for (int rr = 0; rr < _nrRows; rr++)
            for (int cr = 0; cr < _nrCols; cr++) {
                if(LDDChannel->Drcr == 5) {
                    LDD_LINKEDLIST *list = nullptr;
                    LDD_LINKEDLIST *temp = nullptr;
                    list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                    list->prev = nullptr;
                    list->rowNr = rr;
                    list->colNr = cr;

                    while (list != nullptr)
                    {
                        int i = 0;
                        bool  subCatchDone = true;
                        int rowNr = list->rowNr;
                        int colNr = list->colNr;

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
                                temp->prev = list;
                                list = temp;
                                list->rowNr = r;
                                list->colNr = c;
                                subCatchDone = false;
                            }
                        }

                        if (subCatchDone)
                        {
                            double Dx = ChannelDX->data[rowNr][colNr];
                            double tx = dt/Dx;

                            // the current cell
                            double H = ChannelWH->data[rowNr][colNr];
                            double V = ChannelU->data[rowNr][colNr];
                            double Z = DEM->data[rowNr][colNr];
                            double W = ChannelWidth->data[rowNr][colNr];
                            double N = ChannelN->data[rowNr][colNr];

                            double Hn = 0;
                            double Vn = 0;
                            double flux_out = 0;

                            double ch_vadd = 0;

                            //downstream cell
                            int ldd = (int)LDDChannel->data[rowNr][colNr];
                            double DX2 = ldd % 2 == 0 ? Dx : Dx*sqrt(2);
                            int r = rowNr+dy[ldd];
                            int c = colNr+dx[ldd];
                            double Ho = ChannelWH->Drc;
                            double Vo = ChannelU->Drc;
                            double Zo = DEM->Drc;
                            double Wo = ChannelWidth->Drc;

                            vec4 hll_out = F_Riemann(H,V,0, Ho, Vo, 0);
                            double ch_q = tx * 0.5*(H+Ho)*std::min(W,Wo) * hll_out.v[0];
                            double Vol = W*Dx*H;
                            double Volo = Wo*Dx*H;
                            ch_q = std::max(-0.25*Volo, std::min(ch_q, 0.25*Vol));

                            double s_zh_out = std::min(B, std::max(-B, (Ho + Zo - Z - H)/DX2));
                            ch_vadd = ch_vadd + dt * 0.5 * GRAV * s_zh_out;
                            if(ch_q < 0) {
                                Vn = (V*Vol - Vo*ch_q)/std::max(0.01,Vol - ch_q);
                            }
                            //Hn = H - ch_q/(W * Dx);
                            flux_out = ch_q;

                            // upstream cells
                            double flux_in = 0;
                            double ch_vaddw = 1.0;
                            double s_zh_inavg = 0;
                            double nn = 0;
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

                                if (ldd > 0 && FLOWS_TO(ldd, r,c,rowNr,colNr)){
                                    double Hi = ChannelWH->Drc;
                                    double Vi = ChannelU->Drc;
                                    double Zi = DEM->Drc;
                                    double Wi = ChannelWidth->Drc;

                                    vec4 hll_in = F_Riemann(Hi,Vi,0, H,V,0);
                                    double ch_qi = tx * 0.5*(H+Hi)*std::min(W,Wi) * hll_in.v[0];
                                    double Voli = Wi*Dx*Hi;
                                    ch_qi = std::max(-0.25*Vol,std::min(0.25*Voli,ch_qi));

                                    double DX2 = ldd % 2 == 0 ? Dx : Dx*sqrt(2);
                                    double s_zh_in = std::min(B, std::max(-B, (H + Z - Zi - Hi)/DX2));
                                    s_zh_inavg += s_zh_in;//std::min(B, std::max(-B, (H + Z - Zi - Hi)/DX2));

                                    ch_vadd = ch_vadd + dt * 0.5 * GRAV * s_zh_in;
                                    if(ch_qi > 0) {
                                     //   double new_ch_vol =  Hn*(W*Dx); //H * W * Dx; //
                                        Vn = (V*Vol + Vi*ch_qi)/std::max(0.01,Vol + ch_qi);
                                    }
                                    flux_in = flux_in + ch_qi;
                                    nn += 1.0;
                                }
                            }

                            if(ch_vaddw > 1)
                                ch_vadd = ch_vadd/ch_vaddw;
                            if(nn > 0)
                                s_zh_inavg = s_zh_inavg/nn;

                            Hn = H + (flux_out + flux_in)/(W*Dx);
                            Hn = std::max(0.0,Hn);
                            // add incoming fluxes

                            if (Hn > he_ca) {
                                double s_zh = std::min(1.0,std::max(-1.0,limiter(s_zh_inavg, s_zh_out)));
                               // ch_vadd = dt * 0.5 * GRAV * s_zh;
   // qDebug() << shin << s_zh_out << ch_vadd;

                                Vn = Vn - ch_vadd;

                                double chnsq1 = (0.001+N)*(0.001+N)*GRAV/std::max(0.01,pow(Hn,4.0/3.0));
                                double  chnsq = chnsq1*sqrt(Vn*Vn)*dt;
                                Vn = std::min(25.0,std::max(-25.0,Vn/(1.0+chnsq)));

//                                if (SwitchTimeavgV) {
//                                    double fac = 0.5+0.5*std::min(1.0,4*Hn)*std::min(1.0,4*Hn);
//                                    fac = fac *exp(- std::max(1.0,dt) / chnsq1);
//                                    Vn = fac * V + (1.0-fac) *Vn;
//                                }

//                                double threshold = 0.01 * Dx;
//                                if(Hn < threshold)
//                                {
//                                    double h23 = pow(Hn, 2.0/3.0);
//                                    double kinfac = std::max(0.0,(threshold - Hn) / (0.025 * W));
//                                    double v_kin = (s_zh>0?1:-1) * h23 * std::max(0.001, sqrt(s_zh > 0 ? s_zh : -s_zh))/(0.001+N);
//                                    Vn = kinfac * v_kin + Vn*(1.0-kinfac);
//                                }

                            }

                            if (std::isnan(Vn)) {
                                Vn = 0;
                                Hn = 0;
                                qDebug() << "oei";
                            }
                            if (fabs(Vn) <= ve_ca)
                                Vn = 0;
                            if (Hn < he_ca) {
                                Vn = 0;
                                Hn = 0;
                            }
                            dt_req = std::min(dt_req,courant_factor *Dx/( std::min(dt_max,std::max(0.01,fabs(Vn)))));

                            ChannelU->data[rowNr][colNr] = Vn;
                            ChannelWH->data[rowNr][colNr] = Hn;

                            tma->data[rowNr][colNr] = 1; // flag done
                            temp=list;
                            list=list->prev;
                            free(temp);
                        }
                    }
                }
            }

        dt = std::min(dt_req, _dt-timesum);
      //  qDebug() << timesum << dt;

        timesum = timesum + dt;
        if (timesum > _dt-1e-6)
            stop = true;
        count++;
        if (count > 200)
            stop = true;

    } while (!stop);

    qDebug() << count;


//    double sum2 = 0;
//    double n = 0;
//#pragma omp parallel for reduction(+:sum2) collapse(2) num_threads(userCores)
//    FOR_ROW_COL_MV_CH {
//        if(ChannelWH->Drc > 0) {
//            sum2 += ChannelWH->Drc;
//            if(ChannelWH->Drc > 0)
//                n += 1.0;
//        }
//    }
//    double dhtot = sum2 > 0 ? (sum1 - sum2)/sum2 : 0;

//#pragma omp parallel for collapse(2) num_threads(userCores)
//    FOR_ROW_COL_MV_CH {
//        if(ChannelWH->Drc > 0) {
//            ChannelWH->Drc = ChannelWH->Drc*(1.0 + dhtot);
//            ChannelWH->Drc = std::max(ChannelWH->Drc , 0.0);
//        }
//    }


#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_CH {
        ChannelV->Drc = fabs(ChannelU->Drc);
        ChannelQn->Drc = ChannelV->Drc*ChannelWH->Drc*ChannelWidth->Drc;
        ChannelWaterVol->Drc = ChannelWH->Drc*ChannelWidth->Drc*ChannelDX->Drc;
    }
}

