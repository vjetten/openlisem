/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2022  Victor Jetten, Bastian van de Bout
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
//This particular code uses much of the principles of the FullSWOF model
// https://arxiv.org/abs/1204.3210
//https://sourcesup.renater.fr/frs/?group_id=895&release_id=3901#fullswof_2d-_1.10.00-title-content
// the scheme is made suited for parallel processing

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

//---------------------------------------------------------------------------

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
//---------------------------------------------------------------------------
/**
 * @brief TWorld::limiter: Flux limiters are used in high resolution schemes to avoid occilations
 * @param a slope on one side
 * @param b slope on oposite side
 * @return rec
 *
 * ONLY used in MUSCL or ENO
 * WIKI: Flux limiters are used in high resolution schemes, such as the MUSCL scheme, to avoid
 * the spurious oscillations (wiggles) that would otherwise occur with high order spatial
 * discretisation schemes due to shocks, discontinuities or sharp changes in the solution domain.
 * Use of flux limiters, together with an appropriate high resolution scheme, make the solutions
 * total variation diminishing (TVD).
 */
double TWorld::limiter(double a, double b)
{
    double eps = 1.e-15;
    double rec = 0.;
    // F_fluxLimiter=1;
    if (F_fluxLimiter == (int)MINMOD)
    {
        if (a >= 0 && b >= 0)
            rec = std::min(a, b);
        else
            if (a <= 0 && b <= 0)
                rec = std::max(a, b);
    }
    else
    {
        double ab = a*b;

        if (F_fluxLimiter == (int)VANLEER)
        {
            if (ab > 0)
                return (2*ab/(a+b));
        }
        else
            if (F_fluxLimiter == (int)VANALBEDA)
            {
                double aa = a*a;
                double bb = b*b;
                if (ab > 0)
                    rec=(a*(bb+eps)+b*(aa+eps))/(aa+bb+2*eps);
            }
    }
    return(rec);
}


//---------------------------------------------------------------------------

//  1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
//  2e component: Momentum flux in gelijke richting per meter per tijdseenheid  ( dus (m4/s2)/(m) = m3/s2 = h*u*u)
//  3d component: Momentum flux in loodrechte richting per meter per tijdseenheid  ( dus (m4/s2)/(m) = m3/s2 = h*u*v)

//f_hllc2.cpp in swof
vec4 TWorld::F_HLL3(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl;
    double c;
    if (h_L < he_ca && h_R < he_ca){
        c = 0.;
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;

        double sqrt_grav_h_L = sqrt(grav_h_L);  // wave velocity
        double sqrt_grav_h_R = sqrt(grav_h_R);
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        double c1;
        double c2;
        if(h_L < he_ca) {
            c1 = u_R - 2*sqrt_grav_h_R;//sqrt(GRAV*h_R);
        }else{
            c1 = std::min(u_L-sqrt_grav_h_L, u_R-sqrt_grav_h_R); // as u-sqrt(grav_h) <= u+sqrt(grav_h)
        }
        if(h_R < he_ca) {
            c2 = u_L + 2*sqrt_grav_h_L;//sqrt(GRAV*h_L);
        }else{
            c2 = std::max(u_L+sqrt_grav_h_L, u_R+sqrt_grav_h_R); // as u+sqrt(grav_h) >= u-sqrt(grav_h)
        }
        double tmp = 1./(c2-c1);
        double t1 = (std::min(c2,0.) - std::min(c1,0.))*tmp;
        double t2 = 1. - t1;
        double t3 = (c2*fabs(c1) - c1*fabs(c2))*0.5*tmp;
        double c_star = (c1*h_R *(u_R - c2) - c2*h_L *(u_L - c1))/(h_R *(u_R - c2) - h_L *(u_L - c1)) ;

        f1 = t1*q_R+t2*q_L-t3*(h_R-h_L);
        f2 = t1*(q_R*u_R+grav_h_R*h_R*0.5)+t2*(q_L*u_L+grav_h_L*h_L*0.5)-t3*(q_R-q_L);
        if(c_star > EPSILON) {
            f3=f1*v_L;
        }else{
            f3=f1*v_R;
        }
        cfl = std::max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}

vec4 TWorld::F_HLL2(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl;
    if (h_L < he_ca && h_R < he_ca){
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    } else {
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double sqrt_grav_h_L = sqrt(grav_h_L);
        double sqrt_grav_h_R = sqrt(grav_h_R);
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        double c1 = std::min(u_L-sqrt_grav_h_L,u_R-sqrt_grav_h_R);   // as u-sqrt(grav_h) <= u+sqrt(grav_h)
        double c2 = std::max(u_L+sqrt_grav_h_L,u_R+sqrt_grav_h_R);   // as u+sqrt(grav_h) >= u-sqrt(grav_h)
        double tmp = 1./(c2-c1);
        double t1 = (std::min(c2,0.)-std::min(c1,0.))*tmp;
        double t2 = 1.-t1;
        double t3 = (c2*fabs(c1)-c1*fabs(c2))*0.5*tmp;

        f1 = t1*q_R+t2*q_L-t3*(h_R-h_L);
        f2 = t1*(q_R*u_R+grav_h_R*h_R*0.5)+t2*(q_L*u_L+grav_h_L*h_L*0.5)-t3*(q_R-q_L);
        f3 = t1*q_R*v_R+t2*q_L*v_L-t3*(h_R*v_R-h_L*v_L);
        cfl = std::max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}

vec4 TWorld::F_HLL(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl;
    if (h_L < he_ca && h_R < he_ca){

        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double halfL = GRAV*h_L*h_L*0.5;
        double halfR = GRAV*h_R*h_R*0.5;
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        double c1 = std::min(u_L-sqrt(grav_h_L),u_R-sqrt(grav_h_R));
        double c2 = std::max(u_L+sqrt(grav_h_L),u_R+sqrt(grav_h_R));

        //cfl is the velocity to calculate the real cfl=std::max(fabs(c1),fabs(c2))*tx with tx=dt/dx
        if (fabs(c1)<EPSILON && fabs(c2)<EPSILON){              //dry state
            f1=0.;
            f2=0.;
            f3=0.;
            cfl=0.; //std::max(fabs(c1),fabs(c2))=0
        }else if (c1>=EPSILON){ //supercritical flow, from left to right : we have std::max(abs(c1),abs(c2))=c2>0
            f1=q_L;   //flux
            f2=q_L*u_L+halfL;  //flux*velocity + 0.5*(wave velocity squared)
            f3=q_L*v_L; //flux *velocity
            cfl=c2; //std::max(fabs(c1),fabs(c2))=c2>0
        }else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have std::max(abs(c1),abs(c2))=-c1>0
            f1=q_R;
            f2=q_R*u_R+halfR;
            f3=q_R*v_R;
            cfl=fabs(c1); //std::max(fabs(c1),fabs(c2))=fabs(c1)
        }else{ //subcritical flow
            double tmp = 1./(c2-c1);
            f1=(c2*q_L-c1*q_R)*tmp + c1*c2*(h_R-h_L)*tmp;
            f2=(c2*(q_L*u_L+halfL) - c1*(q_R*u_R+halfR))*tmp + c1*c2*(q_R-q_L)*tmp;
            f3=(c2*(q_L*v_L)-c1*(q_R*v_R))*tmp + c1*c2*(h_R*v_R-h_L*v_L)*tmp;
            cfl=std::max(fabs(c1),fabs(c2));
        }
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}
//  1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
//  2e component: Momentum flux in gelijke richting per meter per tijdseenheid  ( dus (m4/s2)/(m) = m3/s2 = h*u*u)
//  3d component: Momentum flux in loodrechte richting per meter per tijdseenheid  ( dus (m4/s2)/(m) = m3/s2 = h*u*v)

vec4 TWorld::F_Rusanov(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl, tmp = 0;
    if (h_L < he_ca && h_R < he_ca){

        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        cfl = std::max(fabs(u_L)+sqrt(GRAV*h_L),fabs(u_R)+sqrt(GRAV*h_R));
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        f1 = ((q_L+q_R) - cfl*(h_R-h_L))*0.5;
        f2 = ((u_L*q_L) + (GRAV_DEM*h_L*h_L) + (u_R*q_R) + (GRAV_DEM*h_R*h_R) - cfl*(q_R-q_L))*0.5;
        f3 = ((q_L*v_L+q_R*v_R) - cfl*(h_R*v_R-h_L*v_L))*0.5;
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}

vec4 TWorld::F_Riemann(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 rec;// = {0,0,0,0};
    //    if (F_scheme == 6)
    //    rec = F_ROE(h_L, u_L, v_L, h_R, u_R, v_R);
    //    else
    //    if (F_scheme == 5)
    //        rec = F_HLL4(h_L, u_L, v_L, h_R, u_R, v_R);
    //    else
    if (F_scheme == 4)
        rec = F_HLL3(h_L, u_L, v_L, h_R, u_R, v_R);
    else
        if (F_scheme == 3)
            rec = F_HLL2(h_L, u_L, v_L, h_R, u_R, v_R);
        else
            if (F_scheme == 2)
                rec = F_HLL(h_L, u_L, v_L, h_R, u_R, v_R);
            else
                if (F_scheme == 1)
                    rec = F_Rusanov( h_L, u_L, v_L, h_R, u_R, v_R);
    return (rec);
}

//--------------------------------------------------------------------------------------------
// correct mass balance
double TWorld::getMass(cTMap *M, double th)
{
    double sum2 = 0;
#pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
            sum2 += M->Drc*CHAdjDX->Drc;
    }}
return sum2;
}
//---------------------------------------------------------------------------
double TWorld::getMassSed(cTMap *M, double th)
{
    double sum2 = 0;
    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
            sum2 += M->Drc;
    }}
    return sum2;
}
//---------------------------------------------------------------------------
// correct mass balance
void TWorld::correctMassBalance(double sum1, cTMap *M, double th)
{
    double sum2 = 0;
    double n = 0;

    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
        {
            sum2 += M->Drc*CHAdjDX->Drc;
            n += 1;
        }
    }}
    sum2 = std::max(0.0, sum2);
    // total and cells active for M
    double dhtot = fabs(sum2) > 0 ? (sum1 - sum2)/sum2 : 0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
        {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            M->Drc = std::max(M->Drc , 0.0);
        }
    }}
}

void TWorld::correctMassBalanceSed(double sum1, cTMap *M, double th)
{
    double sum2 = 0;
    double n = 0;

    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
        {
            sum2 += M->Drc;
            n += 1;
        }
    }}
    // total and cells active for M
    double dhtot = fabs(sum2) > 0 ? (sum1 - sum2)/sum2 : 0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
        {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            M->Drc = std::max(M->Drc , 0.0);
        }
    }}
}
