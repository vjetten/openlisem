/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

// This particular code uses much of the principles of the FullSWOF model
// https://arxiv.org/abs/1204.3210
// https://sourcesup.renater.fr/frs/?group_id=895&release_id=3901#fullswof_2d-_1.10.00-title-content
// the scheme is made suited for parallel processing
// LICENCE: http://cecill.info/licences/Licence_CeCILL_V2-en.html

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

//---------------------------------------------------------------------------

// force flow when a diagonal solution exists and a DEM blockage is present
// runs from inside swof loop because of min dt
void TWorld::SWOFDiagonalFlowLDD(double dt_req_min, cTMap *z, cTMap *h, cTMap *vx, cTMap *vy)
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

            vec4 rec;
            int ldd = dcr_[i_].ldd;
            int rr = r+dy[ldd];
            int cr = c+dx[ldd];

            if (z->Drcr+h->Drcr < z->Drc+h->Drc) {
                // 1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
                rec = F_Riemann(h->Drc, vx->Drc, vy->Drc, h->Drcr, vx->Drcr, vy->Drcr);
                double flux = std::abs(rec.v[0]);
                double dH = qMin(h->Drc*0.5, flux*dt_req_min/_dx);
                double Hldd = z->Drcr+h->Drcr;
                double H = z->Drc+h->Drc;
                int cnt = 0;
                // if movning water causes an imbalance
                if (Hldd+dH > H-dH) {
                    while (Hldd+dH > H-dH && cnt < 100) {
                        dH -= 0.01;
                        cnt++;
                    }
                }

                h->Drc -= dH;
                h->Drc = qMax(0.0,h->Drc);
                tmc->Drcr += dH;

                doit = true;

                if (SwitchErosion) {
                    double dS = qMin(0.5*SSFlood->Drc, dH*CHAdjDX->Drc*SSCFlood->Drc);
                    SSFlood->Drc -= dS;
                    tma->Drcr += dS;
                    if (SwitchUse2Phase) {
                        double dBL = qMin(0.5*BLFlood->Drc, dH*CHAdjDX->Drc*BLCFlood->Drc);
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
    // check in diagonal direction, not with ldd
void TWorld::SWOFDiagonalFlow(double dt_req_min, cTMap *z, cTMap *h, cTMap *vx, cTMap *vy)
{

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tma->Drc = 0;
        tmb->Drc = 0;
        tmc->Drc = 0;
        tmd->Drc = 0;
    }}

    bool doit = false;

    #pragma omp parallel for num_threads(userCores)
    for(long i_= 0; i_ < dcr_.size(); i_++) {

        int r = dcr_[i_].r;
        int c = dcr_[i_].c;

        if (h->Drc > F_pitValue) {
            int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1,  0,  1};
            int dy[10] = {0,  1, 1, 1,  0, 0, 0, -1, -1, -1};

            double H = z->Drc+h->Drc;
            int rr, cr, k, j;

            // hydraulic potential in X and Y directions
            j = 2;
            rr = r+dy[j];
            cr = c+dx[j];
            double H2 = z->Drcr+h->Drcr;
            j = 4;
            rr = r+dy[j];
            cr = c+dx[j];
            double H4 = z->Drcr+h->Drcr;
            j = 6;
            rr = r+dy[j];
            cr = c+dx[j];
            double H6 = z->Drcr+h->Drcr;
            j = 8;
            rr = r+dy[j];
            cr = c+dx[j];
            double H8 = z->Drcr+h->Drcr;

            // if a pit then check the diagonals
            if (H < H2 && H < H4 && H < H6 && H < H8) {
                j = 1;
                rr = r+dy[j];
                cr = c+dx[j];
                double H1 = z->Drcr+h->Drcr;
                j = 3;
                rr = r+dy[j];
                cr = c+dx[j];
                double H3 = z->Drcr+h->Drcr;
                j = 7;
                rr = r+dy[j];
                cr = c+dx[j];
                double H7 = z->Drcr+h->Drcr;
                j = 7;
                rr = r+dy[j];
                cr = c+dx[j];
                double H9 = z->Drcr+h->Drcr;

                double dH1 = H-H1;
                double dH3 = H-H3;
                double dH7 = H-H7;
                double dH9 = H-H9;

                int k = 0;
                int dHfin = 0;
                if (dH1 > 0) {
                    dHfin = dH1;
                    k = 1;
                }
                if (dH3 > dHfin) {
                    dHfin = dH3;
                    k = 3;
                }
                if (dH7 > dHfin) {
                    dHfin = dH7;
                    k = 7;
                }
                if (dH9 > dHfin) {
                    dHfin = dH9;
                    k = 9;
                }
            }
            // a diagonal solution is found
            if (k > 0) {
                doit = true;

                vec4 rec;
                int rr = r+dy[k];
                int cr = c+dx[k];
                rec = F_Riemann(h->Drc, vx->Drc, vy->Drc, h->Drcr, vx->Drcr, vy->Drcr);
                double flux = std::abs(rec.v[0]);
                double dH = qMin(h->Drc*0.5, flux*dt_req_min/_dx);
                double Hdown = z->Drcr+h->Drcr;
                //double H = z->Drc+h->Drc;
                int cnt = 0;
                // if movning water causes an imbalance
                if (Hdown+dH > H-dH) {
                    while (Hdown+dH > H-dH && cnt < 100) {
                        dH -= 0.01;
                        cnt++;
                    }
                }

//                h->Drc -= dH;
//                h->Drc = qMax(0.0,h->Drc);
                tma->Drc = -dH;
                tmb->Drcr = dH;

                if (SwitchErosion) {
                    // just do suspended
                    double dS = qMin(0.5*SSFlood->Drc, dH*CHAdjDX->Drc*SSCFlood->Drc);
                    //SSFlood->Drc -= dS;
                    tmc->Drc = -dS;
                    tmd->Drcr = dS;
                    // if (SwitchUse2Phase) {
                    //     double dBL = qMin(0.5*BLFlood->Drc, dH*CHAdjDX->Drc*BLCFlood->Drc);
                    //     BLFlood->Drc -= dBL;
                    //     tmb->Drcr += dBL;
                    // }
                }
            } // found
        } // pit value
    } // LOOP

    if (doit) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            h->Drc += tma->Drc;
            h->Drc += tmb->Drc;
            h->Drc = qMax(0.0, h->Drc);
        }}

        if (SwitchErosion) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                SSFlood->Drc += tmc->Drc;
                SSFlood->Drc += tmd->Drc;
                // if (SwitchUse2Phase)
                //     BLFlood->Drc += tmb->Drc;
            }}
        }
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

// note, using something else then minmod slows down the simulation very much. no idea why
double TWorld::limiter(double a, double b)
{
    double eps = 1.e-12;

    if (F_fluxLimiter == (int)MINMOD) {
        if (a > 0 && b > 0)
            return qMin(a, b);
        else
            if (a < 0 && b < 0)
                return qMax(a, b);
            else
                return 0.;
    }
    else
        if (F_fluxLimiter == (int)VANLEER) {
            if ((a > 0 && b > 0) || (a < 0 && b < 0))
                return (2*a*b/(a+b));
            else
                return 0.;
        }
        else
            if (F_fluxLimiter == (int)VANALBEDA) {
                if (a*b < 0.)
                    return 0.;
                else
                    return  (a*(b*b+eps)+b*(a*a+eps))/(a*a+b*b+2.*eps);
            }

    return 0;
}


//---------------------------------------------------------------------------

//  1e component: Massa flux per meter ( dus (m3/s)/(m) = m2/s, wat dezelfde berekening is als momentum = h*u)
//  2e component: Momentum flux in gelijke richting per meter per tijdseenheid  ( dus (m4/s2)/(m) = m3/s2 = h*u*u)
//  3d component: Momentum flux in loodrechte richting per meter per tijdseenheid  ( dus (m4/s2)/(m) = m3/s2 = h*u*v)

//f_hllc.cpp in swof
vec4 TWorld::F_HLL4(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
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
    } else {
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double grav_2h_L = GRAV * h_L * h_L*0.5;
        double grav_2h_R = GRAV * h_R * h_R*0.5;
        double sqrt_grav_h_L = sqrt(grav_h_L);  // wave velocity
        double sqrt_grav_h_R = sqrt(grav_h_R);
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        double c1;
        double c2;
        if(h_L < he_ca) {
            c1 = u_R - 2*sqrt_grav_h_R;
        } else {
            c1 = qMin(u_L-sqrt_grav_h_L, u_R-sqrt_grav_h_R); // as u-sqrt(grav_h) <= u+sqrt(grav_h)
        }
        if(h_R < he_ca) {
            c2 = u_L + 2*sqrt_grav_h_L;//sqrt(GRAV*h_L);
        } else {
            c2 = qMax(u_L+sqrt_grav_h_L, u_R+sqrt_grav_h_R); // as u+sqrt(grav_h) >= u-sqrt(grav_h)
        }

        //cfl is the velocity to calculate the real cfl=max(qFabs(c1),qFabs(c2))*tx with tx=dt/dx
        if (qFabs(c1) < EPSILON && qFabs(c2) < EPSILON) {
            //dry state
            f1 = 0.;
            f2 = 0.;
            f3 = 0.;
            cfl = 0.;
        } else
            if (c1 >= EPSILON) {
                //supercritical flow, from left to right : we have max(abs(c1),abs(c2))=c2>0
                f1 = q_L;
                f2 = q_L*u_L+grav_2h_L;
                f3 = q_L*v_L;
                cfl = c2; //max(qFabs(c1),qFabs(c2))=c2>0
            }
            else
                if (c2 <= -EPSILON) {
                    //supercritical flow, from right to left : we have max(abs(c1),abs(c2))=-c1>0
                    f1 = q_R;
                    f2 = q_R*u_R+grav_2h_R;
                    f3 = q_R*v_R;
                    cfl = qFabs(c1); //max(qFabs(c1),qFabs(c2))=qFabs(c1)
                } else {
                    //subcritical flow
                    double c_star = (c1*h_R *(u_R - c2) - c2*h_L *(u_L - c1))/(h_R *(u_R - c2) - h_L *(u_L - c1));
                    double tmp = 1./(c2-c1);
                    f1 = (c2*q_L-c1*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;
                    f2 = (c2*(q_L*u_L+grav_2h_L)-c1*(q_R*u_R+grav_2h_R))*tmp+c1*c2*(q_R-q_L)*tmp;
                    if (c_star > EPSILON) {
                        f3 = f1*v_L;
                    } else {
                        f3 = f1*v_R;
                    }
                    cfl = qMax(qFabs(c1), qFabs(c2));
                }
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}

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
        } else {
            c1 = qMin(u_L-sqrt_grav_h_L, u_R-sqrt_grav_h_R); // as u-sqrt(grav_h) <= u+sqrt(grav_h)
        }
        if(h_R < he_ca) {
            c2 = u_L + 2*sqrt_grav_h_L;//sqrt(GRAV*h_L);
        } else {
            c2 = qMax(u_L+sqrt_grav_h_L, u_R+sqrt_grav_h_R); // as u+sqrt(grav_h) >= u-sqrt(grav_h)
        }
        double tmp = 1./(c2-c1);
        double t1 = (qMin(c2,0.) - qMin(c1,0.))*tmp;
        double t2 = 1. - t1;
        double t3 = (c2*qFabs(c1) - c1*qFabs(c2))*0.5*tmp;
        double c_star = (c1*h_R *(u_R - c2) - c2*h_L *(u_L - c1))/(h_R *(u_R - c2) - h_L *(u_L - c1)) ;

        f1 = t1*q_R+t2*q_L-t3*(h_R-h_L);
        f2 = t1*(q_R*u_R+grav_h_R*h_R*0.5)+t2*(q_L*u_L+grav_h_L*h_L*0.5)-t3*(q_R-q_L);
        if(c_star > EPSILON) {
            f3=f1*v_L;
        }else{
            f3=f1*v_R;
        }
        cfl = qMax(qFabs(c1),qFabs(c2)); //cfl is the velocity to compute the cfl condition max(qFabs(c1),qFabs(c2))*tx with tx=dt/dx
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}


//F_HLL2.cpp in fullswof
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
        double c1 = qMin(u_L-sqrt_grav_h_L,u_R-sqrt_grav_h_R);   // as u-sqrt(grav_h) <= u+sqrt(grav_h)
        double c2 = qMax(u_L+sqrt_grav_h_L,u_R+sqrt_grav_h_R);   // as u+sqrt(grav_h) >= u-sqrt(grav_h)
        double tmp = 1./(c2-c1);
        double t1 = (qMin(c2,0.)-qMin(c1,0.))*tmp;
        double t2 = 1.-t1;
        double t3 = (c2*qFabs(c1)-c1*qFabs(c2))*0.5*tmp;

        f1 = t1*q_R+t2*q_L-t3*(h_R-h_L);
        f2 = t1*(q_R*u_R+grav_h_R*h_R*0.5)+t2*(q_L*u_L+grav_h_L*h_L*0.5)-t3*(q_R-q_L);
        f3 = t1*q_R*v_R+t2*q_L*v_L-t3*(h_R*v_R-h_L*v_L);
        cfl = qMax(qFabs(c1),qFabs(c2)); //cfl is the velocity to compute the cfl condition max(qFabs(c1),qFabs(c2))*tx with tx=dt/dx
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}

// F_HLL.cpp in fullswof
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
        double c1 = qMin(u_L-sqrt(grav_h_L),u_R-sqrt(grav_h_R));
        double c2 = qMax(u_L+sqrt(grav_h_L),u_R+sqrt(grav_h_R));

        //cfl is the velocity to calculate the real cfl=qMax(qFabs(c1),qFabs(c2))*tx with tx=dt/dx
        if (qFabs(c1)<EPSILON && qFabs(c2)<EPSILON){              //dry state
            f1=0.;
            f2=0.;
            f3=0.;
            cfl=0.; //qMax(qFabs(c1),qFabs(c2))=0
        }else if (c1>=EPSILON){ //supercritical flow, from left to right : we have qMax(abs(c1),abs(c2))=c2>0
            f1=q_L;   //flux
            f2=q_L*u_L+halfL;  //flux*velocity + 0.5*(wave velocity squared)
            f3=q_L*v_L; //flux *velocity
            cfl=c2; //qMax(qFabs(c1),qFabs(c2))=c2>0
        }else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have qMax(abs(c1),abs(c2))=-c1>0
            f1=q_R;
            f2=q_R*u_R+halfR;
            f3=q_R*v_R;
            cfl=qFabs(c1); //qMax(qFabs(c1),qFabs(c2))=qFabs(c1)
        }else{ //subcritical flow
            double tmp = 1./(c2-c1);
            f1=(c2*q_L-c1*q_R)*tmp + c1*c2*(h_R-h_L)*tmp;
            f2=(c2*(q_L*u_L+halfL) - c1*(q_R*u_R+halfR))*tmp + c1*c2*(q_R-q_L)*tmp;
            f3=(c2*(q_L*v_L)-c1*(q_R*v_R))*tmp + c1*c2*(h_R*v_R-h_L*v_L)*tmp;
            cfl=qMax(qFabs(c1),qFabs(c2));
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
        cfl = qMax(qFabs(u_L)+sqrt(GRAV*h_L), qFabs(u_R)+sqrt(GRAV*h_R));
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
    vec4 rec;

    if (F_scheme == 5)
        rec = F_HLL4(h_L, u_L, v_L, h_R, u_R, v_R);
    else
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

vec3 TWorld::F_VFRoe(double h_L,double u_L,double h_R,double u_R)
{

    double cL=sqrt(GRAV*h_L);
    double cR=sqrt(GRAV*h_R);
    double umean=(u_L+u_R)/2.;
    double cmean=(cL+cR)/2.;
    double lamb1=umean-cmean;
    double lamb2=umean+cmean;
    double lamb1L=u_L-cL;
    double lamb2L=u_L+cL;
    double lamb1R=u_R-cR;
    double lamb2R=u_R+cR;
    double f1, f2, cfl, c, tx;
    vec3 res;

    if ( ((lamb1L < 0.0) && (lamb1R > 0.0)) ||
         ((lamb2L < 0.0) && (lamb2R > 0.0)) ) {
        //entropy correction with the Rusanov flux
        c = qMax(qFabs(u_L)+cL,qFabs(u_R)+cR);
        f1 = (h_L*u_L+h_R*u_R)*0.5-c*(h_R-h_L)*0.5;
        f2 = (u_L*u_L*h_L + (GRAV_DEM*h_L*h_L) + u_R*u_R*h_R + (GRAV_DEM*h_R*h_R) )*0.5 - c*(h_R*u_R-h_L*u_L)*0.5;
        cfl = c*tx;
    }
    else
        if (lamb1 >= 0.0){
            //supercritical flow from the left to the right
            f1 = h_L*u_L;
            f2 = h_L*u_L*u_L + GRAV_DEM*h_L*h_L;
            cfl = qMax(qFabs(u_L)+cL,qFabs(u_R)+cR)*tx;
        }
        else
            if (lamb2 <= 0.0){
                //supercritical flow from the right to the left
                f1 = h_R*u_R;
                f2 = h_R*u_R*u_R + GRAV_DEM*h_R*h_R;
                cfl = qMax(qFabs(u_L)+cL,qFabs(u_R)+cR)*tx;
            } else {
                //subcritical flow
                double lambmax=0.;
                double ustar=0.;
                double hstar=0.;

                lambmax = qMax(qFabs(lamb1),qFabs(lamb2));
                ustar = (u_L+u_R)/2.0-(cR-cL);
                double tmp = (cR+cL)/2.0-(u_R-u_L)/4.0;
                hstar = tmp*tmp/GRAV;
                f1 = hstar*ustar;
                f2 = hstar*ustar*ustar + GRAV_DEM*hstar*hstar;
                cfl = qMax(lambmax,qMax(qFabs(u_L)+cL,qFabs(u_R)+cR))*tx;
            }
    res.v[0] = f1;
    res.v[1] = f2;
    res.v[2] = cfl;
    return (res);
}


//--------------------------------------------------------------------------------------------
// correct mass balance
double TWorld::getMass(cTMap *M)
{
    double sum2 = 0;
    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > 0)
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
void TWorld::correctMassBalance(double sum1, cTMap *M)
{
    double sum2 = 0;

    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > 0)
            sum2 += M->Drc*CHAdjDX->Drc;
    }}

    double Mcorr = sum2 > 0 ? (1.0+(sum1 - sum2)/sum2) : 1.0;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > 0) {
            M->Drc = M->Drc*Mcorr;            // <- distribution weighted to h
            M->Drc = qMax(M->Drc , 0.0);
        }
    }}
}

void TWorld::correctMassBalanceSed(double sum1, cTMap *M, double th)
{
    double sum2 = 0;

    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
            sum2 += M->Drc;
    }}
    // total and cells active for M
    double Mcorr = qFabs(sum2) > 0 ? (1.0+(sum1 - sum2)/sum2) : 1.0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th) {
            M->Drc = M->Drc*Mcorr;
            M->Drc = qMax(M->Drc , 0.0);
        }
    }}
}
//---------------------------------------------------------------------------
#include <QQueue>
#include <QPoint>

void TWorld::floodFill(cTMap *raster, cTMap* labels, int row, int col, int currentlabel)
{
    int dr[4] = {-1,0,1,0};
    int dc[4] = {0,-1,0,1};
    QQueue<QPoint> q;
    q.enqueue(QPoint(row, col));
    //labels->data[row][col] = currentlabel;

    while (!q.isEmpty()) {
        QPoint pt = q.dequeue();
        int r = pt.x();
        int c = pt.y();

        // Directions: up, down, left, right
        for (int i= 0; i < 4; i++) {
            int nr = r + dr[i];
            int nc = c + dc[i];
            if (nr >= 0 && nr < _nrRows && nc >= 0 && nc < _nrCols) {
                if (raster->data[nr][nc] != 0 && labels->data[nr][nc] == 0) {
                    labels->data[nr][nc] = currentlabel;
                    q.enqueue(QPoint(nr, nc));
                }
            }
        }
    }
}

void TWorld::floodCount(cTMap *h)
{
    Fill(*tmb,0);
    int currentlabel = 1;
    FOR_ROW_COL_MV_L {
        if (h->Drc == 0 and tmb->Drc == 0)
            floodFill(h,tmb, r,c,currentlabel++);
    }}
}
