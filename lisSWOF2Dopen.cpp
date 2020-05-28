
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-6
#define ve_ca 1e-6

#define GRAV 9.8067
#define EPSILON 1e-6

vec4 TWorld::F_HLL4(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl;
    double c;
    if (h_L<=0. && h_R<=0.){
        c = 0.;
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        c = std::max(fabs(u_L)+sqrt(GRAV*h_L),fabs(u_R)+sqrt(GRAV*h_R));
        double cd = c*0.5;
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        f1 = (q_L+q_R)*0.5-cd*(h_R-h_L);
        f2 = ((u_L*q_L)+(GRAV*0.5*h_L*h_L)+(u_R*q_R)+(GRAV*0.5*h_R*h_R))*0.5-cd*(q_R-q_L);
        f3 = (q_L*v_L+q_R*v_R)*0.5-cd*(h_R*v_R-h_L*v_L);
        cfl = c;//*tx;
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}
//---------------------------------------------------------------------------
vec4 TWorld::F_HLL3(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl, tmp = 0;
    if (h_L<=0. && h_R<=0.){
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }
    else
    {
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double sqrt_grav_h_L = sqrt(grav_h_L);  // wave velocity
        double sqrt_grav_h_R = sqrt(grav_h_R);
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;

        double c1 = std::min(u_L - sqrt_grav_h_L,u_R - sqrt_grav_h_R); //we already have u_L - sqrt_grav_h_L<u_L + sqrt_grav_h_L and u_R - sqrt_grav_h_R<u_R + sqrt_grav_h_R
        double c2 = std::max(u_L + sqrt_grav_h_L,u_R + sqrt_grav_h_R); //so we do not need all the eigenvalues to get c1 and c2
        tmp = 1./(c2-c1);
        double t1 = (std::min(c2,0.) - std::min(c1,0.))*tmp;
        double t2 = 1. - t1;
        double t3 = (c2*fabs(c1) - c1*fabs(c2))*0.5*tmp;

        f1 = t1*q_R + t2*q_L - t3*(h_R - h_L);
        f2 = t1*(q_R*u_R + grav_h_R*h_R*0.5) + t2*(q_L*u_L + grav_h_L*h_L*0.5) - t3*(q_R - q_L);
        f3 = t1*q_R*v_R + t2*q_L*v_L - t3*(h_R*v_R - h_L*v_L);
        cfl = std::max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition std::max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;

}
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

    if (!startFlood)
        TimestepfloodLast = dt_max;

    if (startFlood)
    {
        sumh = getMass(h);

        do {

            // make a copy
//#pragma omp parallel for collapse(2)
            for(int r = 0; r < _nrRows; r++)
                for (int c = 0; c < _nrCols; c++)
                if(!pcr::isMV(LDD->data[r][c])) {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
            }
//#pragma omp barrier

            //flow 1
#pragma omp parallel for collapse(2)
            for(int r = 0; r < _nrRows; r++)
                for (int c = 0; c < _nrCols; c++)
                if(!pcr::isMV(LDD->data[r][c])) {
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

                double B = 0.5; // dh+dz hydraulisch verschil mag max 0.5 zijn?
                double sx_zh_x2 = std::min(B,std::max(-B,(z_x2 + h_x2 - Z - H)/dx));
                double sy_zh_y1 = std::min(B,std::max(-B,(Z + H - z_y1 - h_y1)/dy));
                double sx_zh_x1 = std::min(B,std::max(-B,(Z + H - z_x1 - h_x1)/dx));
                double sy_zh_y2 = std::min(B,std::max(-B,(z_y2 + h_y2 - Z - H)/dy));

                // if B = 0.5 this can never be >1?
                double sx_zh = std::min(1.0,std::max(-1.0,minmod(sx_zh_x1, sx_zh_x2)));
                double sy_zh = std::min(1.0,std::max(-1.0,minmod(sy_zh_y1, sy_zh_y2)));

                hll_x1 = F_HLL3(h_x1,vx_x1,vy_x1,H,Vx,Vy);
                hll_x2 = F_HLL3(H,Vx,Vy,h_x2,vx_x2,vy_x2);
                hll_y1 = F_HLL3(h_y1,vy_y1,vx_y1,H,Vy,Vx);
                hll_y2 = F_HLL3(H,Vy,Vx,h_y2,vy_y2,vx_y2);

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
             //   double dt_req = courant_factor *_dx/( std::min(dt_max,std::max(0.01,sqrt(vxn*vxn + vyn*vyn))));

                double dtx = dx/std::max(hll_x1.v[3],hll_x2.v[3]);
                double dty = dy/std::max(hll_y1.v[3],hll_y2.v[3]);
                double dt_req = std::max(TimestepfloodMin, std::min(dt_max, courant_factor*std::min(dtx, dty)));

                FloodDT->Drc = dt_req;
                h->Drc = hn;
                vx->Drc = vxn;
                vy->Drc = vyn;
               // Qn->Drc = fabs(qfout);
            }
#pragma omp barrier
            dt_req_min = dt_max;
            FOR_ROW_COL_MV {
                dt_req_min = std::min(dt_req_min,FloodDT->Drc);
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
            FOR_ROW_COL_MV {
                hs->Drc = h->Drc;
                vxs->Drc = vx->Drc;
                vys->Drc = vy->Drc;
            }

            // do the tango
            FOR_ROW_COL_MV {
                double dt = 0.5*TimestepfloodLast;
                double vxn, vyn;

                //typedef struct vec4 { double v[4]; } vec4;
                vec4 hll_x1;
                vec4 hll_x2;
                vec4 hll_y1;
                vec4 hll_y2;

                double dh,dz_h,dvx,dvy,hlh;
                double h1r,h1l,h2r,h2l,z1r,z1l,z2r,z2l,vx1r,vx1l,vy1r,vy1l,delzc1, delzc2, delz1, delz2;

                double dx = ChannelAdj->Drc;
                double dy = DX->Drc;

                double H = hs->Drc;
                double n = N->Drc;
                double Z = z->Drc;
                double Vx = vx->Drc;
                double Vy = vy->Drc;

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

                //col -1 and +1
                dh   = 0.5*minmod(H-h_x1, h_x2-H);
                dz_h = 0.5*minmod(H-h_x1 + Z-z_x1, h_x2-H + z_x2-Z);
                dvx  = 0.5*minmod(Vx-vx_x1, vx_x2-Vx);
                dvy  = 0.5*minmod(Vy-vy_x1, vy_x2-Vy);

                h1r = H+dh;
                h1l = H-dh;
                z1r = Z+(dz_h-dh);
                z1l = Z+(dh-dz_h);
                delzc1 = z1r-z1l;
                delz1 = z1l - z_x1;

                hlh = H > he_ca ? (H+dh)/H : 1.0;
                vx1r = Vx + hlh*dvx;
                vx1l = Vx - hlh*dvx;
                vy1r = Vy + hlh*dvy;
                vy1l = Vy - hlh*dvy;

                hll_x1 = F_HLL3(h1l,vx1l,vy1l,H,Vx,Vy);
                hll_x2 = F_HLL3(H,Vx,Vy,h1r,vx1r,vy1r);

                // row -1 and +1
                dh   = 0.5*minmod(H-h_y1, h_y2-H);
                dz_h = 0.5*minmod(H-h_y1 + Z-z_y1, h_y2-H + z_y2-Z);
                dvx  = 0.5*minmod(Vx-vx_y1, vx_y2-Vx);
                dvy  = 0.5*minmod(Vy-vy_y1, vy_y2-Vy);

                h2r = H+dh;
                h2l = H-dh;
                z2r = Z+(dz_h-dh);
                z2l = Z+(dh-dz_h);
                delzc2 = z2r-z2l;
                delz2 = z2l - z_y1;

                hlh = H > he_ca ? (H+dh)/H : 1.0;
                vx1r = Vx + hlh*dvx;
                vx1l = Vx - hlh*dvx;
                vy1r = Vy + hlh*dvy;
                vy1l = Vy - hlh*dvy;

                hll_y1 = F_HLL3(h2l,vx1l,vy1l,H,Vx,Vy);
                hll_y2 = F_HLL3(H,Vx,Vy,h2r,vx1r,vy1r);

//                double sx_zh_x2 = (z_x2 + h_x2 - Z - H)/dx;
//                double sy_zh_y1 = (Z + H - z_y1 - h_y1)/dy;
//                double sx_zh_x1 = (Z + H - z_x1 - h_x1)/dx;
//                double sy_zh_y2 = (z_y2 + h_y2 - Z - H)/dy;
//                double sx_zh = minmod(sx_zh_x1, sx_zh_x2);
//                double sy_zh = minmod(sy_zh_y1, sy_zh_y2);

                double h1d = std::max(0.0, h1r - std::max(0.0,  delz1));
                double h1g = std::max(0.0, h1l - std::max(0.0, -delz1));
                double h2d = std::max(0.0, h2r - std::max(0.0,  delz2));
                double h2g = std::max(0.0, h2l - std::max(0.0, -delz2));

//                h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]));
//                h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]));

                double C = courant_factor;//0.2;
                double tx = dt/dx;
                double ty = dt/dy;

                double flux_x1 = std::max(-H * C,std::min(+tx*hll_x1.v[0],h_x1 * C));
                double flux_x2 = std::max(-H * C,std::min(-tx*hll_x2.v[0],h_x2 * C));
                double flux_y1 = std::max(-H * C,std::min(+ty*hll_y1.v[0],h_y1 * C));
                double flux_y2 = std::max(-H * C,std::min(-ty*hll_y2.v[0],h_y2 * C));

                double hn = std::max(0.0, H + flux_x1 + flux_x2 + flux_y1 + flux_y2);

                if(hn > he_ca) {
                    double qxn = H * Vx -ty*(hll_y2.v[2] - hll_y1.v[2]) - tx*(hll_x2.v[1] - hll_x1.v[1] +
                            GRAV*0.5*((h1g-h1l)*(h1g+h1l) + (h1r-h1d)*(h1r+h1d) + (h1l+h1r)*delzc1)) ;
                            // 0.5 * GRAV *((h1l+h1r)*(h1l-h1r)+ (h1l+h1r)*delzc1));
                    double qyn = H * Vy - tx*(hll_x2.v[2] - hll_x1.v[2]) - ty*(hll_y2.v[1] - hll_y1.v[1] +
                             // 0.5 * GRAV *((h2l+h2r)*(h1l-h2r)+(h2l+h2r)*delzc2));
                            GRAV*0.5*((h2g-h2l)*(h2g+h2l) + (h2r-h2d)*(h2r+h2d) + (h2l+h2r)*delzc2)) ;


//                    qes1 = he->Drc*ve1->Drc -
//                            ty*(_g2 - g2->Drc) -
//                            tx*(_f2 - f2->Drc +
//                                GRAV*0.5*((h1g->Drc-h1l->Drc)*(h1g->Drc+h1l->Drc) +
//                                          (h1r->Drc-h1d->Drc)*(h1r->Drc+h1d->Drc)
//                                          + (h1l->Drc+h1r->Drc)*delzc1->Drc)) ;

//                    qes2 = he->Drc*ve2->Drc -
//                            tx*(_f3 - f3->Drc) -
//                            ty*(_g3 - g3->Drc +
//                                GRAV*0.5*((h2g->Drc-h2l->Drc)*(h2g->Drc+h2l->Drc) +
//                                          (h2r->Drc-h2d->Drc)*(h2r->Drc+h2d->Drc)
//                                          + (h2l->Drc+h2r->Drc)*delzc2->Drc));


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

//                    double threshold = 0.01 * _dx;
//                    if(hn < threshold)
//                    {
//                        double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
//                        // float acc_eff = (vxn - vx)/std::max(0.0001,dt);

//                        double v_kin = (sx_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sx_zh>0?sx_zh:-sx_zh))/(0.001+n);

//                        vxn = kinfac * v_kin + vxn*(1.0-kinfac);
//                    }

//                    if(hn < threshold)
//                    {
//                        double kinfac = std::max(0.0,(threshold - hn) / (0.025 * _dx));
//                        // float acc_eff = (vyn -vy)/std::max(0.0001,dt);

//                        double v_kin = (sy_zh>0?1:-1) * hn * sqrt(hn) * std::max(0.001,sqrt(sy_zh>0?sy_zh:-sy_zh))/(0.001+n);

//                        vyn = kinfac * v_kin + vyn*(1.0-kinfac);

//                    }

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

                if (hn <= he_ca)
                {
                    hn = 0;
                    vxn = 0;
                    vyn = 0;
                }

                if (fabs(vxn) <= ve_ca)
                    vxn = 0;
                if (fabs(vyn) <= ve_ca)
                    vyn = 0;

                double dt_req = courant_factor *_dx/( std::min(dt_max,std::max(TimestepfloodMin,sqrt(vxn*vxn + vyn*vyn))));

//                double dtx = std::min(dt_max, courant_factor*dx/(0.5*(hll_x1.at(3)+hll_x2.at(3))));
//                double dty = std::min(dt_max, courant_factor*dy/(0.5*(hll_y1.at(3)+hll_y2.at(3))));
//                dt_req = std::min(dtx, dty);

                FloodDT->Drc = dt_req;
                h->Drc = hn;
                vx->Drc = vxn;
                vy->Drc = vyn;
               // Qn->Drc = fabs(qfout);
            }

            dt_req_min = dt_max;
            FOR_ROW_COL_MV {
                dt_req_min = std::min(dt_req_min,FloodDT->Drc);
            }
          //  dt_req_min = std::max(TimestepfloodMin, dt_req_min);
            dt_req_min = std::min(dt_req_min, _dt-timesum);

            TimestepfloodLast = dt_req_min;

            timesum += dt_req_min;
            stop = timesum > _dt-0.001;
            count++;
            if(count > 1000) stop = true;
        } while (!stop);

    } // if floodstart

    correctMassBalance(sumh, h);

    iter_n = count;

    return(_dt/count);
}



