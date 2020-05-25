
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
double TWorld::fullSWOF2open(cTMap *h, cTMap *u, cTMap *v, cTMap *z)
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
                vs->Drc = v->Drc;
                us->Drc = u->Drc;
            }


            //flow 1
            FOR_ROW_COL_MV {
                double tx = dt/ChannelAdj->Drc;
                double ty = dt/DX->Drc;
                double dx = ChannelAdj->Drc;
                double dy = DX->Drc;

                double H = hs->Drc;
              //  double qfout = Qn->Drc;
                double n = N->Drc;
                double Z = z->Drc;

                double zc_x1 = c > 0 && !pcr::isMV(LDD->data[r][c-1])         ? z->data[r][c-1] : z->Drc;
                double zc_x2 = c < _nrCols-1 && !pcr::isMV(LDD->data[r][c+1]) ? z->data[r][c+1] : z->Drc;
                double zc_y1 = r > 0 && !pcr::isMV(LDD->data[r-1][c])         ? z->data[r-1][c] : z->Drc;
                double zc_y2 = r < _nrRows-1 && !pcr::isMV(LDD->data[r+1][c]) ? z->data[r+1][c] : z->Drc;

                double h_x1 =  c > 0 && !pcr::isMV(LDD->data[r][c-1])          ? hs->data[r][c-1] : hs->Drc;
                double h_x2 =  c < _nrCols-1 && !pcr::isMV(LDD->data[r][c+1])  ? hs->data[r][c+1] : hs->Drc;
                double h_y1 =  r > 0 && !pcr::isMV(LDD->data[r-1][c])          ? hs->data[r-1][c] : hs->Drc;
                double h_y2 =  r < _nrRows-1 && !pcr::isMV(LDD->data[r+1][c])  ? hs->data[r+1][c] : hs->Drc;

                double vmax = 1000;//1.0 * _dx/_dt;    // !!!!! dt1?

                double vx = std::min(vmax,std::max(-vmax, u->Drc));
                double vy = std::min(vmax,std::max(-vmax, v->Drc));

                double vx_x1 = std::min(vmax,std::max(-vmax,c > 0 && !pcr::isMV(LDD->data[r][c-1])          ? us->data[r][c-1] : us->Drc ));
                double vy_x1 = std::min(vmax,std::max(-vmax,c < _nrCols-1 && !pcr::isMV(LDD->data[r][c+1])  ? us->data[r][c+1] : us->Drc ));
                double vx_x2 = std::min(vmax,std::max(-vmax,r > 0 && !pcr::isMV(LDD->data[r-1][c])          ? us->data[r-1][c] : us->Drc ));
                double vy_x2 = std::min(vmax,std::max(-vmax,r < _nrRows-1 && !pcr::isMV(LDD->data[r+1][c])  ? us->data[r+1][c] : us->Drc ));

                double vx_y1 = std::min(vmax,std::max(-vmax,c > 0 && !pcr::isMV(LDD->data[r][c-1])          ? vs->data[r][c-1] : vs->Drc ));
                double vy_y1 = std::min(vmax,std::max(-vmax,c < _nrCols-1 && !pcr::isMV(LDD->data[r][c+1])  ? vs->data[r][c+1] : vs->Drc ));
                double vx_y2 = std::min(vmax,std::max(-vmax,r > 0 && !pcr::isMV(LDD->data[r-1][c])          ? vs->data[r-1][c] : vs->Drc ));
                double vy_y2 = std::min(vmax,std::max(-vmax,r < _nrRows-1 && !pcr::isMV(LDD->data[r+1][c])  ? vs->data[r+1][c] : vs->Drc ));

                double sx_zh_x2 = std::min(0.5,std::max(-0.5,(zc_x2 + h_x2 - Z - H)/dx));
                double sy_zh_y1 = std::min(0.5,std::max(-0.5,(Z + H-zc_y1  -  h_y1)/dy));
                double sx_zh_x1 = std::min(0.5,std::max(-0.5,(Z + H-zc_x1  -  h_x1)/dx));
                double sy_zh_y2 = std::min(0.5,std::max(-0.5,(zc_y2 + h_y2 - Z - H)/dy));

                double sx_zh = std::min(1.0,std::max(-1.0,minmod(sx_zh_x1,sx_zh_x2)));
                double sy_zh = std::min(1.0,std::max(-1.0,minmod(sy_zh_y1,sy_zh_y2)));

                hll_x1.clear();
                hll_x2.clear();
                hll_y1.clear();
                hll_y2.clear();
                F_HLL2(h_x1,vx_x1,vy_x1,H,vx,vy);
                hll_x1 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                F_HLL2(H,vx,vy,h_x2,vx_x2,vy_x2);
                hll_x2 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                F_HLL2(h_y1,vy_y1,vx_y1,H,vy,vx);
                hll_y1 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                F_HLL2(H,vy,vx,h_y2,vy_y2,vx_y2);
                hll_y2 << HLL2_f1 << HLL2_f2 << HLL2_f3<< HLL2_cfl;
                double C = 0.1;

//                double flux_x1 = std::max(-H * C,std::min(+tx*hll_x1.at(0),h_x1 * C));
           //     double flux_x2 = std::max(-H * C,std::min(-tx*hll_x2.at(0),h_x2 * C));
           //     double flux_y1 = std::max(-H * C,std::min(+ty*hll_y1.at(0),h_y1 * C));
           //     double flux_y2 = std::max(-H * C,std::min(-ty*hll_y2.at(0),h_y2 * C));

                double flux_x1 = std::min(+tx*hll_x1.at(0),h_x1 * C);
                double flux_x2 = std::min(-tx*hll_x2.at(0),h_x2 * C);
                double flux_y1 = std::min(+ty*hll_y1.at(0),h_y1 * C);
                double flux_y2 = std::min(-ty*hll_y2.at(0),h_y2 * C);


//                qfout = qfout + flux_x1;
//                qfout = qfout + flux_x2;
//                qfout = qfout + flux_y1;
//                qfout = qfout + flux_y2;

                double hn = std::max(0.0, H + flux_x1 + flux_x2 + flux_y1 + flux_y2);
              //   hn = std::max(0.0, H + tx * (hll_x1.at(0)-hll_x2.at(0))+ty*(hll_y1.at(0)-hll_y2.at(0)) );

                if(hn > he_ca) {

                    double qxn = H * vx - tx*(hll_x2.at(1) - hll_x1.at(1)) - ty*(hll_y2.at(2) - hll_y1.at(2))- 0.5 * GRAV *hn*sx_zh * dt;
                    double qyn = H * vy - tx*(hll_x2.at(2) - hll_x1.at(2)) - ty*(hll_y2.at(1) - hll_y1.at(1))- 0.5 * GRAV *hn*sy_zh * dt;

                    double vsq = sqrt(vx * vx + vy * vy);
                    double nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01,pow(hn,4.0/3.0));
                    double nsq = nsq1*vsq*dt;

                    vxn = (qxn/(1.0+nsq))/std::max(0.01,hn);
                    vyn = (qyn/(1.0+nsq))/std::max(0.01,hn);

                    double fac = 0;
                    if (SwitchTimeavgV) {
                        fac = 0.5+0.5*std::min(1.0,4*hn)*std::min(1.0,4*hn);
                        fac = fac *exp(- std::max(1.0,dt) / nsq1);
                    }
                    vxn = fac * vx + (1.0-fac) *vxn;
                    vyn = fac * vy + (1.0-fac) *vyn;

                    //                double threshold = 0.01 * _dx;
                    //                if(hn < threshold)
                    //                {
                    //                    double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
                    //                   // float acc_eff = (vxn - vx)/std::max(0.0001,dt);

                    //                    double v_kin = (sx_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sx_zh>0?sx_zh:-sx_zh))/(0.001+n);

                    //                    vxn = kinfac * v_kin + vxn*(1.0-kinfac);
                    //                }

                    //                if(hn < threshold)
                    //                {
                    //                    double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
                    //                   // float acc_eff = (vyn -vy)/std::max(0.0001,dt);

                    //                    double v_kin = (sy_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sy_zh>0?sy_zh:-sy_zh))/(0.001+n);

                    //                    vyn = kinfac * v_kin + vyn*(1.0-kinfac);

                    //                }

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

                double dt_req = courant_factor *_dx/( std::min(100.0,std::max(0.01,(sqrt(vxn*vxn + vyn * vyn)))));

               // double dtx = std::min(dt_max, courant_factor*dx/(0.5*(hll_x1.at(3)+hll_x2.at(3))));
               // double dty = std::min(dt_max, courant_factor*dy/(0.5*(hll_y1.at(3)+hll_y2.at(3))));
               // dt_req = std::min(dtx, dty);

                FloodDT->Drc = dt_req;
                h->Drc = hn;
                u->Drc = vxn;
                v->Drc = vyn;
               // Qn->Drc = qfout;
            }

            dt_req_min = dt_max;
            FOR_ROW_COL_MV {
                dt_req_min = std::min(dt_req_min,FloodDT->Drc);
            }

            dt_req_min = std::min(dt_req_min, _dt-timesum);

            dt = dt_req_min;

            timesum += dt_req_min;
            stop = timesum > _dt-0.001;
            count++;
            if(count > 200) stop = true;


        } while (!stop);
    } // if floodstart

    correctMassBalance(sumh, h);

//    iter_n = n;
//    dt1 = n > 0? _dt/n : dt1;

    return(dt_req_min);
}


