

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


// NOT USED !!!

//-----------------------------------------------------------------------------------------------------
void TWorld::ChannelSWOFopen()
{
    if(!SwitchIncludeChannel)
        return;

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    bool CorrectMassBalance = false;
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


#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ChannelQn = 0;
            tmb->Drc = 0;
}}

//    fill(*ChannelQn, 0);
//    fill(*tmb, 0);

    do {
        stop = false;

        // do the whole channel
        //fill(*tma, -1);
        dt_req = dt_max;

#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            double Dx = ChannelDX->Drc;
            double tx = dt/Dx;

            // ===== the current cell =====
            double H = ChannelWH->Drc;
            double V = ChannelU->Drc;
            double Z = DEM->Drc;
            double W = ChannelWidth->Drc;
            double N = ChannelN->Drc;
            double ch_height = ChannelDepth->Drc;
            double Hmx = hmx->Drc;

            // new H and V
            double Hn = H;
            double Vn = V;

            // ===== outflow to downstream cell =====
            double ch_vadd = 0; //gravity component
            double flux_out = 0; // flux to downstream cell as volume

            vec4 hll_out = {0,0,0,0};
            int ldd = (int)LDDChannel->Drc;

            double Vo = V;//pow(H, 2.0/3.0)/N*sqrt(H/_dx+0.001);
            double Ho = H;//std::max(0.0, H*(1-Vo/Dx*dt));
            double Zo = Z;//*0.99;//std::max(0.99, (1.0-ChannelGrad->data[rowNr][colNr])); // at least a 1% slope
            double Wo = W;

            if (ldd == 5) {

                Vo = pow(H, 2.0/3.0)/N*sqrt(ChannelGrad->Drc);
                Zo = Z;
                Wo = W;

                double Q = W*H*Vo*dt;
                Q = std::min(C* W*Dx*H,Q);
                Hn = Hn - Q/(W*Dx);
                ch_vadd = ch_vadd + dt * 0.5 * GRAV * std::max(-B,std::min(B,H/Dx+0.001));
                flux_out = flux_out + Q;

            } else {

                int rr = r+dy[ldd];
                int cr = c+dx[ldd];
                Ho = ChannelWH->Drcr;
                Vo = ChannelU->Drcr;
                Zo = DEM->Drcr;
                Wo = ChannelWidth->Drcr;

                // float3 hll_x1 = F_HLL2(ch_h,ch_v,0,chn_h,chn_v,0);
                // float ch_q = (dt/dx)*(max(ch_width,chn_width)/dx)*((dx * 0.5*(chn_width +ch_width)) *hll_x1.x);
                // ch_q = min(0.25f * ch_vol,ch_q);
                // ch_q = max(-0.25f * chn_vol,ch_q);
                // ch_q = ch_q * 0.5;
                // float ch_slope = (z + ch_h - chn_z - chn_h)/dx;
                // ch_vadd = ch_vadd + dt * 0.5 * GRAV * max(-1.0f,min(1.0f,(float)(ch_slope)));
                // if(ch_q < 0)
                // {
                //         float new_ch_vol = chhn*(ch_width*dx);
                //         chvn = (chvn * new_ch_vol - chn_v *(ch_q))/max(0.01f,new_ch_vol - ch_q);
                // }
                // chhn = chhn - ch_q/(ch_width * dx);
                // flux_chx2 = flux_chx2 + ch_q;

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
                flux_out = flux_out + Q; //note m3? must be m3/s for use as Qn

            }

            // ===== weighed sum inflow from upstream cells =====
            double flux_in = 0;
            double ch_vaddw = 0.5;
            for (int i = 1; i <= 9; i++)
            {
                int rr, cr, ldd = 0;

                if (i==5)
                    continue;

                rr = r+dy[i];
                cr = c+dx[i];

                if (INSIDE(rr, cr) && !pcr::isMV(LDDChannel->Drcr))
                    ldd = (int) LDDChannel->Drcr;
                else
                    continue;

                // if the cells flows into the c
                if (ldd > 0 && ldd != 5 && FLOWS_TO(ldd, rr,cr,r,c)) {
                    double Hi = ChannelWH->Drcr;
                    double Vi = ChannelU->Drcr;
                    double Zi = DEM->Drcr;
                    double Wi = ChannelWidth->Drcr;

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
               // Vn = std::min(25.0,std::max(-25.0,Vn));

                // if (SwitchTimeavgV) {
                //     double fac = 0.5+0.5*std::min(1.0,4*Hn)*std::min(1.0,4*Hn);
                //     fac = fac *exp(- std::max(1.0,dt) / chnsq1);
                //     Vn = fac * V + (1.0-fac) *Vn;
                // }
            } else {
                Vn = 0;
                Hn = 0;
            }

            if (fabs(Vn) <= ve_ca)
                Vn = 0;

            if(Hn < ch_height && Hmx > 0)
            {
                double vol_room = (ch_height - Hn)*W;
                double v_chin = Vn;//Hmx * sqrt(Hmx) * sqrt(Hmx/_dx)/(N*N);
                double flux_in = std::min(vol_room, Hmx*(_dx*Dx) * dt*v_chin/(0.5*std::max(1.0,_dx - W)));

                Hn = Hn + flux_in/(Dx * W);
                Hmx = Hmx - flux_in/(_dx * Dx);
            } else
                if(Hn > ch_height)
                {
                    double vol_tomuch = (Hn - ch_height)*W;
                    Hn = ch_height;
                    Hmx = Hmx + vol_tomuch/(_dx * Dx);
                }
            dt_req = std::min(dt_req,courant_factor *Dx/( std::min(dt_max,std::max(0.01,fabs(Vn)))));
            //std::max(TimestepfloodMin,

            // gebruik riemann solver cfl
            //   double dtx = Dx/hll_out.v[3];
            //   dt_req = std::max(TimestepfloodMin, std::min(dt_req, courant_factor*dtx));

            ChannelU->Drc = Vn;
            ChannelWH->Drc = Hn;
            hmx->Drc = Hmx;

            tmb->Drc += flux_out;

        }}
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
                sum2 += ChannelWH->Drc; // *ChannelWidth->Drc*ChannelDX->Drc;
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

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        ChannelV->Drc = fabs(ChannelU->Drc);
        ChannelQn->Drc = std::max(0.0,ChannelU->Drc)*ChannelWH->Drc*ChannelWidth->Drc;//tmb->Drc/_dt;//
        ChannelQ->Drc = ChannelQn->Drc;
        ChannelWaterVol->Drc = ChannelWH->Drc*ChannelWidth->Drc*ChannelDX->Drc;
        ChannelAlpha->Drc = ChannelWH->Drc*ChannelWidth->Drc/std::pow(ChannelQn->Drc, 0.6);
    }}
    //report(*ChannelU,"chu");



    //channel inflow and outflow



}


void TWorld::KinematicSWOFopen(cTMap *_h, cTMap *_V)
{

    double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int count = 0;
    bool stop;
    double dt = dt_max;
    double dt_req = dt_max;
 //   double qout = 0;

    do {
        stop = false;

        dt_req = dt_max;

#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            double C = 0.25;//std::min(0.25, courant_factor);
            int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
            int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
            double Dx = DX->Drc;


            // ===== the current cell =====
            double H = _h->Drc;
            double V = _V->Drc;
            double W = ChannelAdj->Drc;//FlowWidth->Drc;
            double n = N->Drc;

            if(V == 0)
                V = pow(H, 2.0/3.0)/n*sqrt(Grad->Drc);

            // new H and V
            double Hn = H;
            double Vn = V;

            // ===== outflow to downstream cell =====
            vec4 hll_out = {0,0,0,0};
            int ldd = (int)LDD->Drc;

            double Vo = V;
            double Ho = H;
            double Wo = W;

            if (ldd == 5) {
                // if outlet use manning
                Vo = pow(H, 2.0/3.0)/n*sqrt(Grad->Drc);
                double Q = W*H*Vo*dt;
                Q = std::min(C* W*Dx*H,Q);
                Hn = Hn - Q/(W*Dx);
            } else {
                //downstream cell
                int rr = r+dy[ldd];
                int cr = c+dx[ldd];
                Ho = _h->Drcr;
                Vo = _V->Drcr;
                Wo = ChannelAdj->Drcr;//FlowWidth->Drcr;

                hll_out = F_Riemann(H,V,0, Ho,Vo,0);
                // 1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
                //double Q = tx * hll_out.v[0] * (std::max(W,Wo)/Dx * (Dx*0.5*(W+Wo))); // s/m * m2/s * m/m * m2
                double Q = dt * hll_out.v[0] * std::max(W,Wo);
                double Volo = Wo*Dx*Ho;
                double Vol = W*Dx*H;
                Q = std::max(-C*Volo, std::min(Q, C*Vol));
                if(Q < 0) {
                    Vn = (Vn*Vol - Vo*Q)/std::max(0.01,Vol - Q);
                }
                Hn = Hn - Q/(W*Dx);
            }

            // ===== weighed sum inflow from upstream cells =====
            for (int i = 1; i <= 9; i++)
            {
                int rr, cr, ldd = 0;

                if (i==5)
                    continue;

                rr = r+dy[i];
                cr = c+dx[i];

                if (INSIDE(rr, cr) && !pcr::isMV(LDD->Drcr))
                    ldd = (int) LDD->Drcr;
                else
                    continue;

                // if the cells flows into the c
                if (ldd > 0 && ldd != 5 && FLOWS_TO(ldd, rr,cr,r,c)) {
                    double Hi = _h->Drcr;
                    double Vi = _V->Drcr;
                    double Wi = ChannelAdj->Drcr;//FlowWidth->Drcr;
                    vec4 hll_in = F_Riemann(Hi,Vi,0, H,V,0);
                    double Q = dt * hll_in.v[0] * std::max(W,Wi);
                    double Voli = Wi*Dx*Hi;
                    double Vol = H*W*Dx;
                    Q = std::max(-C*Vol,std::min(C*Voli,Q));
                   // if (Q > 0)
                        Hn = Hn + Q/(W * Dx); // add all Qin to Hn
                }
            }

            Hn = std::max(Hn, 0.0);
            if (Hn > he_ca) {
                double qv = sqrt(Vn*Vn);
                double nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01,pow(Hn,4.0/3.0));
                double nsq = nsq1*qv*dt;
                Vn = (H*V)/(1.0+nsq)/std::max(0.01,Hn);
               // Vn = std::min(25.0,std::max(-25.0,Vn));

                 if (SwitchTimeavgV) {
                     double fac = 0.5+0.5*std::min(1.0,4*Hn)*std::min(1.0,4*Hn);
                     fac = fac *exp(- std::max(1.0,dt) / nsq1);
                     Vn = fac * V + (1.0-fac) *Vn;
                 }
            } else {
                Vn = 0;
                Hn = 0;
            }

            if (fabs(Vn) <= ve_ca)
                Vn = 0;
            _V->Drc = Vn;
            _h->Drc = Hn;
            dt_req = std::min(dt_req,courant_factor *Dx/( std::min(dt_max,std::max(0.01,fabs(Vn)))));

        }}

        //dt_req = 0.25f *min(dt_req,dx/( min(100.0f,max(0.01f,(sqrt(chvn*chvn))))));
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
qDebug() << count;




}


