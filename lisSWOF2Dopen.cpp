
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

    if (!startFlood)
        TimestepfloodLast = dt_max;

    if (startFlood)
    {
        sumh = getMass(h);

#pragma omp parallel for collapse(2)
        FOR_ROW_COL_MV_L {
            FloodDT->Drc = dt_max;
            //FloodT->Drc = 0;
        }

        do {

            // make a copy
            FOR_ROW_COL_MV_L {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
            }
            //flow
            #pragma omp parallel for collapse(2)
            FOR_ROW_COL_MV_L {
                double dt = TimestepfloodLast;//FloodDT->Drc;//
                double vxn, vyn;

                //typedef struct vec4 { double v[4]; } vec4;
                vec4 hll_x1;
                vec4 hll_x2;
                vec4 hll_y1;
                vec4 hll_y2;

                double dx = ChannelAdj->Drc;
                double dy = DX->Drc;
                double vmax = 0.5 * dx/dt;  // courant?

                double H = hs->Drc;
                double n = N->Drc;
                double Z = z->Drc;
                double Vx = std::max(-vmax, std::min(vmax, vx->Drc));
                double Vy = std::max(-vmax, std::min(vmax, vy->Drc));

                double z_x1 =  c > 0 && !MV(r,c-1)         ? z->data[r][c-1] : Z;
                double z_x2 =  c < _nrCols-1 && !MV(r,c+1) ? z->data[r][c+1] : Z;
                double z_y1 =  r > 0 && !MV(r-1,c)         ? z->data[r-1][c] : Z;
                double z_y2 =  r < _nrRows-1 && !MV(r+1,c) ? z->data[r+1][c] : Z;

                double h_x1 =  c > 0 && !MV(r,c-1)         ? hs->data[r][c-1] : hs->Drc;
                double h_x2 =  c < _nrCols-1 && !MV(r,c+1) ? hs->data[r][c+1] : hs->Drc;
                double h_y1 =  r > 0 && !MV(r-1,c)         ? hs->data[r-1][c] : hs->Drc;
                double h_y2 =  r < _nrRows-1 && !MV(r+1,c) ? hs->data[r+1][c] : hs->Drc;

                double vx_x1 = c > 0 && !MV(r,c-1)         ? vxs->data[r][c-1] : vxs->Drc;
                double vx_x2 = c < _nrCols-1 && !MV(r,c+1) ? vxs->data[r][c+1] : vxs->Drc;
                double vx_y1 = r > 0 && !MV(r-1,c)         ? vxs->data[r-1][c] : vxs->Drc;
                double vx_y2 = r < _nrRows-1 && !MV(r+1,c) ? vxs->data[r+1][c] : vxs->Drc;

                double vy_x1 = c > 0 && !MV(r,c-1)         ? vys->data[r][c-1] : vys->Drc;
                double vy_x2 = c < _nrCols-1 && !MV(r,c+1) ? vys->data[r][c+1] : vys->Drc;
                double vy_y1 = r > 0 && !MV(r-1,c)         ? vys->data[r-1][c] : vys->Drc;
                double vy_y2 = r < _nrRows-1 && !MV(r+1,c) ? vys->data[r+1][c] : vys->Drc;

                double fb_x1=0,fb_x2=0,fb_y1=0,fb_y2=0;
                if (SwitchFlowBarriers) {
                    fb_x1 = std::max(FlowBarrierW->Drc, FlowBarrierE->data[r][c-1]);
                    fb_x2 = std::max(FlowBarrierE->Drc, FlowBarrierE->data[r][c+1]);
                    fb_y1 = std::max(FlowBarrierN->Drc, FlowBarrierS->data[r-1][c]);
                    fb_y2 = std::max(FlowBarrierS->Drc, FlowBarrierN->data[r+1][c]);
                }

                vx_x1 = std::max(-vmax, std::min(vmax, vx_x1));
                vx_x2 = std::max(-vmax, std::min(vmax, vx_x2));
                vx_y1 = std::max(-vmax, std::min(vmax, vx_y1));
                vx_y2 = std::max(-vmax, std::min(vmax, vx_y2));

                vy_x1 = std::max(-vmax, std::min(vmax, vy_x1)); //left
                vy_x2 = std::max(-vmax, std::min(vmax, vy_x2)); //right
                vy_y1 = std::max(-vmax, std::min(vmax, vy_y1)); //up
                vy_y2 = std::max(-vmax, std::min(vmax, vy_y2)); //down

                // No effect of terrain: use for lakes?
//                hll_x1 = F_Riemann(h_x1,vx_x1,vy_x1,H,Vx,Vy); // c-1 and c
//                hll_x2 = F_Riemann(H,Vx,Vy,h_x2,vx_x2,vy_x2); // c and c+1
//                hll_y1 = F_Riemann(h_y1,vy_y1,vx_y1,H,Vy,Vx); // r-1 and r
//                hll_y2 = F_Riemann(H,Vy,Vx,h_y2,vy_y2,vx_y2); // r and r+1

                double fac = DEMdz->Drc;
                double dz_x1 = fac*(Z - z_x1);
                double dz_x2 = fac*(z_x2 - Z);
                double dz_y1 = fac*(Z - z_y1);
                double dz_y2 = fac*(z_y2 - Z);

                double h_x1l = std::max(0.0, h_x1 - std::max(0.0,  dz_x1 + fb_x1));
                double h_x1r = std::max(0.0, H    - std::max(0.0, -dz_x1 + fb_x1));
                if(c > 0 && !MV(r,c-1))
                    hll_x1 = F_Riemann(h_x1l,vx_x1,vy_x1,h_x1r,Vx,Vy); // c-1 and c
                else
                    hll_x1 = F_Riemann(0,0,0,h_x1r,Vx,Vy); // c-1 and c

                double h_x2l = std::max(0.0, H    - std::max(0.0,  dz_x2 + fb_x2));
                double h_x2r = std::max(0.0, h_x2 - std::max(0.0, -dz_x2 + fb_x2));
                if(c < _nrCols-1 && !MV(r,c+1))
                    hll_x2 = F_Riemann(h_x2l,Vx,Vy,h_x2r,vx_x2,vy_x2); // c and c+1
                else
                    hll_x2 = F_Riemann(h_x2l,Vx,Vy,0,0,0); // c and c+1

                double h_y1u = std::max(0.0, h_y1 - std::max(0.0,  dz_y1 + fb_y1));
                double h_y1d = std::max(0.0, H    - std::max(0.0, -dz_y1 + fb_y1));
                if (r > 0 && !MV(r-1,c))
                    hll_y1 = F_Riemann(h_y1u,vy_y1,vx_y1,h_y1d,Vy,Vx); // r-1 and r
                else
                    hll_y1 = F_Riemann(0,0,0,h_y1d,Vy,Vx); // r-1 and r

                double h_y2u = std::max(0.0, H    - std::max(0.0,  dz_y2 + fb_y2));
                double h_y2d = std::max(0.0, h_y2 - std::max(0.0, -dz_y2 + fb_y2));
                if( r < _nrRows-1 && !MV(r+1,c))
                    hll_y2 = F_Riemann(h_y2u,Vy,Vx,h_y2d,vy_y2,vx_y2); // r and r+1
                else
                    hll_y2 = F_Riemann(h_y2u,Vy,Vx,0,0,0); // r and r+1

                double B = 0.5; // dh+dz hydraulisch verschil mag max 0.5 zijn?
                //1.0 is theoretical max else falling faster than gravity
                double sx_zh_x2 = std::min(B,std::max(-B,(z_x2 + h_x2 - Z - H)/dx));
                double sy_zh_y1 = std::min(B,std::max(-B,(Z + H - z_y1 - h_y1)/dy));
                double sx_zh_x1 = std::min(B,std::max(-B,(Z + H - z_x1 - h_x1)/dx));
                double sy_zh_y2 = std::min(B,std::max(-B,(z_y2 + h_y2 - Z - H)/dy));

                // if B = 0.5 this can never be >1?
                double sx_zh = std::min(1.0,std::max(-1.0,limiter(sx_zh_x1, sx_zh_x2)));
                double sy_zh = std::min(1.0,std::max(-1.0,limiter(sy_zh_y1, sy_zh_y2)));
                double C = std::min(0.5, courant_factor);
                double tx = dt/dx;
                double ty = dt/dy;

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
                        double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
                        double v_kin = (sx_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001, sqrt(sx_zh > 0 ? sx_zh : -sx_zh))/(0.001+n);
                        vxn = kinfac * v_kin + vxn*(1.0-kinfac);
                        v_kin = (sy_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001, sqrt(sy_zh > 0 ? sy_zh : -sy_zh))/(0.001+n);
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

                // werkt allebij!
                double dt_req = courant_factor *_dx/( std::min(dt_max,std::max(0.01,sqrt(vxn*vxn + vyn*vyn))));

                double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]);
                dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));

                FloodDT->Drc = dt_req;
                h->Drc = hn;
                vx->Drc = vxn;
                vy->Drc = vyn;
            }


#pragma omp parallel for reduction(min:dt_req_min) collapse(2)
            FOR_ROW_COL_MV {
                double res = FloodDT->Drc;
                dt_req_min = std::min(dt_req_min, res);
            }

            dt_req_min = std::max(TimestepfloodMin, dt_req_min);
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            TimestepfloodLast = dt_req_min;
            timesum += dt_req_min;

            FOR_ROW_COL_MV_L {
                FloodDT->Drc = dt_req_min;
            }

            if (SwitchErosion)
                SWOFSediment(FloodDT,hs,vxs,vys);

            stop = timesum > _dt-0.001;
            count++;

            //            qDebug() << timesum << count;

            if(count > 1000) stop = true;
        } while (!stop);

    } // if floodstart

    correctMassBalance(sumh, h);
    if (count == 0) count =1;
    iter_n = count;

    return(_dt/count);
}


//              //  if (SwitchMUSCL) {
//                    double dh1   = 0.5*limiter(H-h_x1, h_x2-H);
//                    double dvx1  = 0.5*limiter(Vx-vx_x1, vx_x2-Vx);
//                    double dvy1  = 0.5*limiter(Vy-vy_x1, vy_x2-Vy);

//                    double hlh = H > he_ca ? (H+dh1)/H : 1.0;
//                    vx_x2 = Vx + hlh*dvx1; //xright
//                    vx_x1 = Vx - hlh*dvx1; //xleft
//                    vy_x2 = Vy + hlh*dvy1;
//                    vy_x1 = Vy - hlh*dvy1;

//                    // row -1 and +1
//                    double dh2   = 0.5*limiter(H-h_y1, h_y2-H);
//                    double dvx2  = 0.5*limiter(Vx-vx_y1, vx_y2-Vx);
//                    double dvy2  = 0.5*limiter(Vy-vy_y1, vy_y2-Vy);
//                    double hlh2 = H > he_ca ? (H+dh2)/H : 1.0;
//                    vy_x2 = Vx + hlh2*dvx2;
//                    vy_x1 = Vx - hlh2*dvx2;
//                    vy_y2 = Vy + hlh2*dvy2;
//                    vy_y1 = Vy - hlh2*dvy2;
//              //  }


//                u1r->Drc = _u->Drc + hlh * du;
//                u1l->Drc = _u->Drc - hrh * du;
//                v1r->Drc = _v->Drc + hlh * dv;
//                v1l->Drc = _v->Drc - hrh * dv;
//                h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
//                h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
//                rec = F_Riemann(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);


