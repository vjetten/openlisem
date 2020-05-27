
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-6
#define ve_ca 1e-6

#define GRAV 9.8067
#define EPSILON 1e-6


double TWorld::minmod(double a, double b)
{   double rec = 0.;
    if (a >= 0 && b >= 0)
        rec = std::min(a, b);
    else
        if (a <= 0 && b <= 0)
            rec = std::max(a, b);
    return rec;
}

// 2nd order without iteration dt1, dt2!
double TWorld::fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z)
{
    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    double sumh = 0;
    bool stop;
    double dt_req_min = dt_max;
    double vxn, vyn;

    QList <double> hll_x1;
    QList <double> hll_x2;
    QList <double> hll_y1;
    QList <double> hll_y2;

    if (startFlood)
    {
        sumh = getMass(h);

        do {

            double dt = dt_max;
            // make a copy
            FOR_ROW_COL_MV {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
            }


            //flow 1
            FOR_ROW_COL_MV {

                double dx = ChannelAdj->Drc;
                double dy = DX->Drc;
                double tx = dt/dx;
                double ty = dt/dy;
                double vmax = 0.1 * dx/dt;

                double H = hs->Drc;
                double n = N->Drc;
                double Z = z->Drc;
                double Vx = std::max(-vmax, std::min(vmax, vx->Drc));
                double Vy = std::max(-vmax, std::min(vmax, vy->Drc));
                  double qfout = Qn->Drc;

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


//                double dh   = 0.5*minmod(H-h_x1, h_x2-H);
//                double dz_h = 0.5*minmod(H-h_x1 + Z-z_x1, h_x2-H + z_x2-Z);
//                double dvx  = 0.5*minmod(Vx-vx_x1, vx_x2-Vx);
//                double dvy  = 0.5*minmod(Vy-vy_x1, vy_x2-Vy);

//                double h1r = H+dh;
//                double h1l = H-dh;
//                double z1r = Z+(dz_h-dh);
//                double z1l = Z+(dh-dz_h);
//                double hlh = H > he_ca ? (H+dh)/H : 1.0;
//                double vx1r = Vx + hlh*dvx;
//                double vx1l = Vx - hlh*dvx;
//                double vy1r = Vy + hlh*dvy;
//                double vy1l = Vy - hlh*dvy;
//                double h1d = std::max(0.0, h1r - std::max(0.0, z1l-z1r));
//                double h1g = std::max(0.0, h1l - std::max(0.0, z1r-z1l));
//                F_HLL2(h1d, vx1r, vy1r, h1g, vx1l, vy1l);

                double B = 0.5; // dh+dz hydraulisch verschil mag max 0.5 zijn?
                double sx_zh_x2 = std::min(B,std::max(-B,(z_x2 + h_x2 - Z - H)/dx));
                double sy_zh_y1 = std::min(B,std::max(-B,(Z + H - z_y1 - h_y1)/dy));
                double sx_zh_x1 = std::min(B,std::max(-B,(Z + H - z_x1 - h_x1)/dx));
                double sy_zh_y2 = std::min(B,std::max(-B,(z_y2 + h_y2 - Z - H)/dy));

                // if B = 0.5 this can never be >1?
                double sx_zh = std::min(1.0,std::max(-1.0,minmod(sx_zh_x1, sx_zh_x2)));
                double sy_zh = std::min(1.0,std::max(-1.0,minmod(sy_zh_y1, sy_zh_y2)));

                hll_x1.clear();
                hll_x2.clear();
                hll_y1.clear();
                hll_y2.clear();
                F_HLL2(h_x1,vx_x1,vy_x1,H,Vx,Vy);
                hll_x1 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                F_HLL2(H,Vx,Vy,h_x2,vx_x2,vy_x2);
                hll_x2 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                F_HLL2(h_y1,vy_y1,vx_y1,H,Vy,Vx);
                hll_y1 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                F_HLL2(H,Vy,Vx,h_y2,vy_y2,vx_y2);
                hll_y2 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                double C = 0.1;
                //void TWorld::F_HLL2(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)


                double flux_x1 = std::max(-H * C,std::min(+tx*hll_x1.at(0),h_x1 * C));
                double flux_x2 = std::max(-H * C,std::min(-tx*hll_x2.at(0),h_x2 * C));
                double flux_y1 = std::max(-H * C,std::min(+ty*hll_y1.at(0),h_y1 * C));
                double flux_y2 = std::max(-H * C,std::min(-ty*hll_y2.at(0),h_y2 * C));

//                double flux_x1 = std::min(+tx*hll_x1.at(0),h_x1 * C);
//                double flux_x2 = std::min(-tx*hll_x2.at(0),h_x2 * C);
//                double flux_y1 = std::min(+ty*hll_y1.at(0),h_y1 * C);
//                double flux_y2 = std::min(-ty*hll_y2.at(0),h_y2 * C);


                qfout = qfout + flux_x1;
                qfout = qfout + flux_x2;
                qfout = qfout + flux_y1;
                qfout = qfout + flux_y2;

                double hn = std::max(0.0, H + flux_x1 + flux_x2 + flux_y1 + flux_y2);
              //   hn = std::max(0.0, H + tx * (hll_x1.at(0)-hll_x2.at(0))+ty*(hll_y1.at(0)-hll_y2.at(0)) );

                if(hn > he_ca) {

                    double qxn = H * Vx - tx*(hll_x2.at(1) - hll_x1.at(1)) - ty*(hll_y2.at(2) - hll_y1.at(2))- 0.5 * GRAV *hn*sx_zh * dt;
                    double qyn = H * Vy - tx*(hll_x2.at(2) - hll_x1.at(2)) - ty*(hll_y2.at(1) - hll_y1.at(1))- 0.5 * GRAV *hn*sy_zh * dt;

                    double vsq = sqrt(Vx * Vx + Vy * Vy);
                    double nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01,pow(hn,4.0/3.0));
                    double nsq = nsq1*vsq*dt;

                    vxn = (qxn/(1.0+nsq))/std::max(0.01,hn);
                    vyn = (qyn/(1.0+nsq))/std::max(0.01,hn);

                    double fac = 0;
                    if (SwitchTimeavgV) {
                        fac = 0.5+0.5*std::min(1.0,4*hn)*std::min(1.0,4*hn);
                        fac = fac *exp(- std::max(1.0,dt) / nsq1);
                    }
                    vxn = fac * Vx + (1.0-fac) *vxn;
                    vyn = fac * Vy + (1.0-fac) *vyn;

                    double threshold = 0.01 * _dx;
                    if(hn < threshold)
                    {
                        double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
                        // float acc_eff = (vxn - vx)/std::max(0.0001,dt);

                        double v_kin = (sx_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sx_zh>0?sx_zh:-sx_zh))/(0.001+n);

                        vxn = kinfac * v_kin + vxn*(1.0-kinfac);
                    }

                    if(hn < threshold)
                    {
                        double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
                        // float acc_eff = (vyn -vy)/std::max(0.0001,dt);

                        double v_kin = (sy_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sy_zh>0?sy_zh:-sy_zh))/(0.001+n);

                        vyn = kinfac * v_kin + vyn*(1.0-kinfac);

                    }

                } else {
                    vxn = 0;
                    vyn = 0;
                  //  qfout = 0;
                }
                // dan maar even met geweld!
                if (std::isnan(vxn) || std::isnan(vyn)  )
                {
                    vxn = 0;
                    vyn = 0;
                    hn= 0;
                }

                double dt_req = courant_factor *_dx/( std::min(dt_max,std::max(0.01,sqrt(vxn*vxn + vyn*vyn))));

//                double dtx = std::min(dt_max, courant_factor*dx/(0.5*(hll_x1.at(3)+hll_x2.at(3))));
//                double dty = std::min(dt_max, courant_factor*dy/(0.5*(hll_y1.at(3)+hll_y2.at(3))));
//                dt_req = std::min(dtx, dty);

                FloodDT->Drc = dt_req;
                h->Drc = hn;
                vx->Drc = vxn;
                vy->Drc = vyn;
               // Qn->Drc = fabs(qfout);
            }
            setZeroOF(h, vx, vy);

            dt_req_min = dt_max;
            FOR_ROW_COL_MV {
                dt_req_min = std::min(dt_req_min,FloodDT->Drc);
            }
          //  dt_req_min = std::max(TimestepfloodMin, dt_req_min);
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            dt = dt_req_min;

            timesum += dt_req_min;
            stop = timesum > _dt-0.001;
            count++;
            if(count > 1000) stop = true;


        } while (!stop);
    } // if floodstart

    correctMassBalance(sumh, h);

    iter_n = count;

    return(count > 0? _dt/count : dt_max);
}


