
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-10
#define ve_ca 1e-10

#define GRAV 9.8067
#define EPSILON 1e-10

//double TWorld::minmod(double a, double b)
//{   double rec = 0.;
//    if (a >= 0 && b >= 0)
//        rec = std::min(a, b);
//    else
//        if (a <= 0 && b <= 0)
//            rec = std::max(a, b);
//    return rec;
//}

// 2nd order without iteration dt1, dt2!
double TWorld::fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;

    if (SwitchErosion) {
        FOR_ROW_COL_MV_L {
            SSFlood->Drc += DETSplash->Drc;
            SSCFlood->Drc = MaxConcentration(ChannelAdj->Drc * DX->Drc * h->Drc, &SSFlood->Drc, &DepFlood->Drc);
            // recalc concentration
        }
    }

    if (!startFlood)
        TimestepfloodLast = dt_max;

    if (startFlood)
    {
        sumh = getMass(h);
#pragma omp parallel for collapse(2)
        FOR_ROW_COL_MV_L {
            FloodDT->Drc = dt_max;
            FloodT->Drc = 0;
        }
qDebug() << "hier";
        do {

            // make a copy
            FOR_ROW_COL_MV_L {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
                FloodHMaskDer->Drc = 1;
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
                double vmax = 0.5 * dx/dt;

                double H = hs->Drc;
                double n = N->Drc;
                double Z = z->Drc;
                double Vx = std::max(-vmax, std::min(vmax, vx->Drc));
                double Vy = std::max(-vmax, std::min(vmax, vy->Drc));

                double z_x1 =  c > 0 && !MV(r,c-1)         ? z->data[r][c-1] : z->Drc;
                double z_x2 =  c < _nrCols-1 && !MV(r,c+1) ? z->data[r][c+1] : z->Drc;
                double z_y1 =  r > 0 && !MV(r-1,c)         ? z->data[r-1][c] : z->Drc;
                double z_y2 =  r < _nrRows-1 && !MV(r+1,c) ? z->data[r+1][c] : z->Drc;

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

                vx_x1 = std::max(-vmax, std::min(vmax, vx_x1));
                vx_x2 = std::max(-vmax, std::min(vmax, vx_x2));
                vx_y1 = std::max(-vmax, std::min(vmax, vx_y1));
                vx_y2 = std::max(-vmax, std::min(vmax, vx_y2));

                vy_x1 = std::max(-vmax, std::min(vmax, vy_x1));
                vy_x2 = std::max(-vmax, std::min(vmax, vy_x2));
                vy_y1 = std::max(-vmax, std::min(vmax, vy_y1));
                vy_y2 = std::max(-vmax, std::min(vmax, vy_y2));

                double B = 1.0;//0.5; // dh+dz hydraulisch verschil mag max 0.5 zijn?
                double sx_zh_x2 = std::min(B,std::max(-B,(z_x2 + h_x2 - Z - H)/dx));
                double sy_zh_y1 = std::min(B,std::max(-B,(Z + H - z_y1 - h_y1)/dy));
                double sx_zh_x1 = std::min(B,std::max(-B,(Z + H - z_x1 - h_x1)/dx));
                double sy_zh_y2 = std::min(B,std::max(-B,(z_y2 + h_y2 - Z - H)/dy));

                // if B = 0.5 this can never be >1?
                double sx_zh = std::min(1.0,std::max(-1.0,limiter(sx_zh_x1, sx_zh_x2)));
                double sy_zh = std::min(1.0,std::max(-1.0,limiter(sy_zh_y1, sy_zh_y2)));

                hll_x1 = F_Riemann(h_x1,vx_x1,vy_x1,H,Vx,Vy); // c-1 and c
                hll_x2 = F_Riemann(H,Vx,Vy,h_x2,vx_x2,vy_x2); // c and c+1
                hll_y1 = F_Riemann(h_y1,vy_y1,vx_y1,H,Vy,Vx); // r-1 and r
                hll_y2 = F_Riemann(H,Vy,Vx,h_y2,vy_y2,vx_y2); // r and r+1

                double C = std::max(0.5, courant_factor);;
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

                    double threshold = 0.001 * _dx; // was 0.01
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
                FloodHMaskDer->Drc = 1.0; // needed for sed
            }

            if (SwitchErosion)
                SWOFSediment(0, FloodDT,hs,vxs,vys);

            stop = timesum > _dt-0.001;
            count++;

            //            qDebug() << timesum << count;

            if(count > 1000) stop = true;
        } while (!stop);

    } // if floodstart

    correctMassBalance(sumh, h);

    iter_n = count;

    return(_dt/count);
}


/*
// 2nd order without iteration dt1, dt2!
double TWorld::fullSWOF2open2(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
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

        do {

            // make a copy
#pragma omp parallel for
            for(long i = 0; i < _nrCells; i++)
            {
                int c = xcl[i];
                int r = ycl[i];
                FloodDT->Drc = dt_max;
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
                FloodHMaskDer->Drc = 0;
            }

#pragma omp parallel for
            for(long i = 0; i < _nrCells; i++)
            {
                int c = xcl[i];
                int r = ycl[i];
                if(h->data[r][c] > HMIN)
                {
                    FloodHMaskDer->Drc = 1.0;

                    if(!OUTORMV(r+1,c))
                        FloodHMaskDer->data[r+1][c] = 1.0;
                    if(!OUTORMV(r-1,c))
                        FloodHMaskDer->data[r-1][c] = 1.0;
                    if(!OUTORMV(r,c+1))
                        FloodHMaskDer->data[r][c+1] = 1.0;
                    if(!OUTORMV(r,c-1))
                        FloodHMaskDer->data[r][c-1] = 1.0;
                }
            }
            //flow
#pragma omp parallel for
            for(long i = 0; i < _nrCells; i++)
            {
                int c = xcl[i];
                int r = ycl[i];
                if (FloodHMaskDer->Drc > 0) {

                    double dt = 0.5*TimestepfloodLast;
                    double vxn, vyn;

                    //typedef struct vec4 { double v[4]; } vec4;
                    vec4 hll_x1;
                    vec4 hll_x2;
                    vec4 hll_y1;
                    vec4 hll_y2;

                    double dx = ChannelAdj->Drc;
                    double dy = DX->Drc;
                    double vmax = 0.5 * dx/dt;

                    double H = hs->Drc;
                    double n = N->Drc;
                    double Z = z->Drc;
                    double Vx = std::max(-vmax, std::min(vmax, vx->Drc));
                    double Vy = std::max(-vmax, std::min(vmax, vy->Drc));

                    double z_x1 =  c > 0 && !MV(r,c-1)         ? z->data[r][c-1] : z->Drc;
                    double z_x2 =  c < _nrCols-1 && !MV(r,c+1) ? z->data[r][c+1] : z->Drc;
                    double z_y1 =  r > 0 && !MV(r-1,c)         ? z->data[r-1][c] : z->Drc;
                    double z_y2 =  r < _nrRows-1 && !MV(r+1,c) ? z->data[r+1][c] : z->Drc;

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

                    vx_x1 = std::max(-vmax, std::min(vmax, vx_x1));
                    vx_x2 = std::max(-vmax, std::min(vmax, vx_x2));
                    vx_y1 = std::max(-vmax, std::min(vmax, vx_y1));
                    vx_y2 = std::max(-vmax, std::min(vmax, vx_y2));

                    vy_x1 = std::max(-vmax, std::min(vmax, vy_x1));
                    vy_x2 = std::max(-vmax, std::min(vmax, vy_x2));
                    vy_y1 = std::max(-vmax, std::min(vmax, vy_y1));
                    vy_y2 = std::max(-vmax, std::min(vmax, vy_y2));


                    // col -1 and +1
                    //                double dh1   = 0.5*minmod(H-h_x1, h_x2-H);
                    //                double dvx1  = 0.5*minmod(Vx-vx_x1, vx_x2-Vx);
                    //                double dvy1  = 0.5*minmod(Vy-vy_x1, vy_x2-Vy);
                    //                double h1r = H+dh1;
                    //                double h1l = H-dh1;
                    //                double hlh1 = H > he_ca ? (H+dh1)/H : 1.0;
                    //                double vx1r = Vx + hlh1*dvx1;
                    //                double vx1l = Vx - hlh1*dvx1;
                    //                double vy1r = Vy + hlh1*dvy1;
                    //                double vy1l = Vy - hlh1*dvy1;
                    //                // row -1 and +1
                    //                double dh2   = 0.5*minmod(H-h_y1, h_y2-H);
                    //                double dvx2  = 0.5*minmod(Vx-vx_y1, vx_y2-Vx);
                    //                double dvy2  = 0.5*minmod(Vy-vy_y1, vy_y2-Vy);
                    //                double h2r = H+dh2;
                    //                double h2l = H-dh2;
                    //                double hlh2 = H > he_ca ? (H+dh2)/H : 1.0;
                    //                double vx2r = Vx + hlh2*dvx2;
                    //                double vx2l = Vx - hlh2*dvx2;
                    //                double vy2r = Vy + hlh2*dvy2;
                    //                double vy2l = Vy - hlh2*dvy2;

                    //                hll_x1 = F_HLL3(h_x1,vx_x1,vy_x1, h1l,vx1l,vy1l); // c-1 and c
                    //                hll_x2 = F_HLL3(h1r,vx1r,vy1r,    h_x2,vx_x2,vy_x2); // c and c+1
                    //                hll_y1 = F_HLL3(h_y1,vy_y1,vx_y1, h2l,vx2l,vy2l); // r-1 and r
                    //                hll_y2 = F_HLL3(h2r,vx2r,vy2r,    h_y2,vy_y2,vx_y2); // r and r+1



                    double B = 0.5; // dh+dz hydraulisch verschil mag max 0.5 zijn?
                    double sx_zh_x2 = std::min(B,std::max(-B,(z_x2 + h_x2 - Z - H)/dx));
                    double sy_zh_y1 = std::min(B,std::max(-B,(Z + H - z_y1 - h_y1)/dy));
                    double sx_zh_x1 = std::min(B,std::max(-B,(Z + H - z_x1 - h_x1)/dx));
                    double sy_zh_y2 = std::min(B,std::max(-B,(z_y2 + h_y2 - Z - H)/dy));

                    // if B = 0.5 this can never be >1?
                    double sx_zh = std::min(1.0,std::max(-1.0,limiter(sx_zh_x1, sx_zh_x2)));
                    double sy_zh = std::min(1.0,std::max(-1.0,limiter(sy_zh_y1, sy_zh_y2)));

                    hll_x1 = F_Riemann(h_x1,vx_x1,vy_x1,H,Vx,Vy); // c-1 and c
                    hll_x2 = F_Riemann(H,Vx,Vy,h_x2,vx_x2,vy_x2); // c and c+1
                    hll_y1 = F_Riemann(h_y1,vy_y1,vx_y1,H,Vy,Vx); // r-1 and r
                    hll_y2 = F_Riemann(H,Vy,Vx,h_y2,vy_y2,vx_y2); // r and r+1

                    double C = courant_factor;//0.2;
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

                        double threshold = 0.001 * _dx; // was 0.01
                        if(hn < threshold)
                        {
                            double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));

                            double v_kin = (sx_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001, sqrt(sx_zh > 0 ? sx_zh : -sx_zh))/(0.001+n);

                            vxn = kinfac * v_kin + vxn*(1.0-kinfac);
                        }

                        if(hn < threshold)
                        {
                            double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));

                            double v_kin = (sy_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001, sqrt(sy_zh > 0 ? sy_zh : -sy_zh))/(0.001+n);

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
                    //    double dt_req = courant_factor *_dx/( std::min(dt_max,std::max(0.01,sqrt(vxn*vxn + vyn*vyn))));

                    double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                    double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]);
                    double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));


                    //                if (SwitchErosion)
                    //                    SWOFSediment(0, FloodDT, h, Vx, Vy);  //TODO why not hs, us, vs

                    FloodDT->Drc = dt_req;
                    h->Drc = hn;
                    vx->Drc = vxn;
                    vy->Drc = vyn;
                }
            }

            dt_req_min = dt_max;

#pragma omp parallel for reduction(min:dt_req_min)
            for(long i = 0; i < _nrCells; i++)
            {
                int c = xcl[i];
                int r = ycl[i];
                double res = FloodDT->Drc;
                dt_req_min = std::min(dt_req_min, res);
            }

            dt_req_min = std::max(TimestepfloodMin, dt_req_min);
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            TimestepfloodLast = dt_req_min;
            timesum += dt_req_min;

            stop = timesum > _dt-0.001;
            count++;

            //            qDebug() << timesum << count;

            if(count > 1000) stop = true;
        } while (!stop);

    } // if floodstart

    correctMassBalance(sumh, h);

    iter_n = count;

    return(_dt/count);
}

*/

