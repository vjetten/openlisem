
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-10
#define ve_ca 1e-10

#define GRAV 9.8067
#define EPSILON 1e-10
//--------------------------------------------------------------------------------------------
// correct mass balance
double TWorld::getMass(cTMap *M, double th)
{
    double sum2 = 0;
#pragma omp parallel for reduction(+:sum2) collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV
    {
        if(M->Drc > th)
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
    }
    return sum2;
}
//---------------------------------------------------------------------------
// correct mass balance
void TWorld::correctMassBalance(double sum1, cTMap *M, double th)
{
    double sum2 = 0;
    double n = 0;

#pragma omp parallel for reduction(+:sum2) collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(M->Drc > th)
        {
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
            n += 1;
        }
    }
    // total and cells active for M
    double dhtot = fabs(sum2) > 0 ? (sum1 - sum2)/sum2 : 0;

#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        if(M->Drc > th)
        {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            M->Drc = std::max(M->Drc , 0.0);
        }
    }
}
//---------------------------------------------------------------------------
// used in datainit, done once
void TWorld::prepareFloodZ(cTMap *z)
{
    prepareFlood = false;

    fill(*delz1,0);
    fill(*delz2,0);
    // diff between z cell and adjacent
    for (int r = 0; r < _nrRows; r++)
        for (int c = 1; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) && !pcr::isMV(LDD->data[r][c-1]))
            {
                delz1->data[r][c-1] = (z->Drc) - z->data[r][c-1];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }
    for (int r = 1; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) && !pcr::isMV(LDD->data[r-1][c]))
            {
                delz2->data[r-1][c] = z->Drc - z->data[r-1][c];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }

    //    FOR_ROW_COL_MV_L {
    //        if (!MV(r,c+1) && !MV(r,c-1) && c > 0 && c < _nrCols-1) {
    //            delzc1->Drc = limiter((z->data[r][c+1] - z->Drc),(z->Drc - z->data[r][c-1]));
    //        }
    //        if (!MV(r+1,c) && !MV(r-1,c) && r > 0 && r < _nrRows-1) {
    //            delzc2->Drc = limiter((z->data[r+1][c] - z->Drc),(z->Drc - z->data[r-1][c]));
    //        }
    //    }
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

/// Numerical flux calculation on which the new velocity is based
/// U_n+1 = U_n + dt/dx* [flux]  when flux is calculated by HLL, HLL2, Rusanov
/// HLL = Harten, Lax, van Leer numerical solution

vec4 TWorld::F_HLL3(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl,tmp;
    double c;
    if (h_L<=0. && h_R<=0.){
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
        double c1;// = std::min(u_L - sqrt_grav_h_L,u_R - sqrt_grav_h_R); //we already have u_L - sqrt_grav_h_L<u_L + sqrt_grav_h_L and u_R - sqrt_grav_h_R<u_R + sqrt_grav_h_R
        double c2;// = std::max(u_L + sqrt_grav_h_L,u_R + sqrt_grav_h_R); //so we do not need all the eigenvalues to get c1 and c2
        if(h_L < he_ca) {
            c1 = u_R - 2*sqrt(GRAV*h_R);
        }else{
            c1 = std::min(u_L-sqrt_grav_h_L,u_R-sqrt_grav_h_R); // as u-sqrt(grav_h) <= u+sqrt(grav_h)
        }
        if(h_R < he_ca) {
            c2 = u_L + 2*sqrt(GRAV*h_L);
        }else{
            c2 = std::max(u_L+sqrt_grav_h_L,u_R+sqrt_grav_h_R); // as u+sqrt(grav_h) >= u-sqrt(grav_h)
        }
        tmp = 1./std::max(0.1,c2-c1);
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
        tmp = 1./std::max(0.1,c2-c1);
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

vec4 TWorld::F_HLL(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    vec4 hll;
    double f1, f2, f3, cfl, tmp = 0;
    if (h_L<=0. && h_R<=0.){
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
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
            f2=q_L*u_L+GRAV*h_L*h_L*0.5;  //flux*velocity + 0.5*(wave velocity squared)
            f3=q_L*v_L; //flux *velocity
            cfl=c2; //std::max(fabs(c1),fabs(c2))=c2>0
        }else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have std::max(abs(c1),abs(c2))=-c1>0
            f1=q_R;
            f2=q_R*u_R+GRAV*h_R*h_R*0.5;
            f3=q_R*v_R;
            cfl=fabs(c1); //std::max(fabs(c1),fabs(c2))=fabs(c1)
        }else{ //subcritical flow
            tmp = 1./(c2-c1);
            f1=(c2*q_L-c1*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;
            f2=(c2*(q_L*u_L+GRAV*h_L*h_L*0.5)-c1*(q_R*u_R+GRAV*h_R*h_R*0.5))*tmp+c1*c2*(q_R-q_L)*tmp;
            f3=(c2*(q_L*v_L)-c1*(q_R*v_R))*tmp+c1*c2*(h_R*v_R-h_L*v_L)*tmp;
            cfl=std::max(fabs(c1),fabs(c2));
        }
    }
    hll.v[0] = f1;
    hll.v[1] = f2;
    hll.v[2] = f3;
    hll.v[3] = cfl;
    return hll;
}

vec4 TWorld::F_Rusanov(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
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
        f1 = (q_L+q_R)*0.5-cd*(h_R-h_L); //m*m/s
        f2 = ((u_L*q_L)+(GRAV*0.5*h_L*h_L)+(u_R*q_R)+(GRAV*0.5*h_R*h_R))*0.5-cd*(q_R-q_L); //m/s*m2/s
        f3 = (q_L*v_L+q_R*v_R)*0.5-cd*(h_R*v_R-h_L*v_L);
        cfl = c;//*tx;
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
    if (F_scheme == 1)
        rec = F_Rusanov( h_L, u_L, v_L, h_R, u_R, v_R);
    else
        if (F_scheme == 2)
            rec = F_HLL(h_L, u_L, v_L, h_R, u_R, v_R);
        else
            if (F_scheme == 3)
                rec = F_HLL2(h_L, u_L, v_L, h_R, u_R, v_R);
            else
                if (F_scheme == 4)
                    rec = F_HLL3(h_L, u_L, v_L, h_R, u_R, v_R);
    return (rec);
}
//---------------------------------------------------------------------------
//OBSOLETE
void TWorld::correctSpuriousVelocities(int r, int c, cTMap *hes, cTMap *ves1, cTMap *ves2) //, double thv, double dv, double dt)
{
    double sign1 = ves1->Drc > 0 ? 1.0 : -1.0;
    double sign2 = ves2->Drc > 0 ? 1.0 : -1.0;
    double s1 = Grad->Drc, s2 = Grad->Drc;
    double G = sqrt(2*GRAV*hes->Drc); // bernouilly pressure velocity

    if (sign1 < 0) {
        if (c > 0 && !MV(r, c-1))
            s1 = sin(atan(fabs(hes->Drc-hes->data[r][c-1])));
    } else {
        if (c < _nrCols-1 && !MV(r, c+1))
            s1 = sin(atan(fabs(hes->Drc-hes->data[r][c+1])));
    }
    if (sign2 < 0) {
        if (r > 0 && !MV(r-1, c))
            s2 = sin(atan(fabs(hes->Drc-hes->data[r-1][c])));
    } else {
        if (r < _nrRows-1 && !MV(r+1, c))
            s2 = sin(atan(fabs(hes->Drc-hes->data[r+1][c])));
    }

    double U1 =  (pow(hes->Drc,2.0/3.0)*sqrt(s1+Grad->Drc)/N->Drc + G);
    double V1 =  (pow(hes->Drc,2.0/3.0)*sqrt(s2+Grad->Drc)/N->Drc + G);

    ves1->Drc = sign1 * std::min(fabs(ves1->Drc), U1);
    ves2->Drc = sign2 * std::min(fabs(ves2->Drc), V1);
    // when V is much larger than kinematic wave V + pressure flow, limit it to that

}
//---------------------------------------------------------------------------

void TWorld::simpleSchemeOF(cTMap *_h,cTMap *_u,cTMap *_v)
{
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L{
        h1r->Drc = _h->Drc;
        u1r->Drc = _u->Drc;
        v1r->Drc = _v->Drc;
        h1l->Drc = _h->Drc;
        u1l->Drc = _u->Drc;
        v1l->Drc = _v->Drc;

        h2r->Drc = _h->Drc;
        u2r->Drc = _u->Drc;
        v2r->Drc = _v->Drc;
        h2l->Drc = _h->Drc;
        u2l->Drc = _u->Drc;
        v2l->Drc = _v->Drc;
    }
}

void TWorld::setZeroOF(cTMap *_h, cTMap *_u, cTMap *_v)
{
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L  {
        if (_h->Drc <= he_ca)
        {
            _h->Drc = 0;
            _u->Drc = 0;
            _v->Drc = 0;
        }

        if (fabs(_u->Drc) <= ve_ca)
            _u->Drc = 0;
        if (fabs(_v->Drc) <= ve_ca)
            _v->Drc = 0;
    }
}
//---------------------------------------------------------------------------

void TWorld::MUSCLOF(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z)
{

#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (_h->Drc > he_ca)
        {
            double delta_h1 = 0;
            double delta_u1 = 0;
            double delta_v1 = 0;
            double delta_h2 = 0;
            double delta_u2 = 0;
            double delta_v2 = 0;
            double dz1 = 0;
            double dz2 = 0;

            double hx = _h->Drc;
            double ux = _u->Drc;
            double vx = _v->Drc;
            double zx = _z->Drc;

            if(c > 0 && !MV(r,c-1)) {
                double hx1 = _h->data[r][c-1];
                double ux1 = _u->data[r][c-1];
                double vx1 = _v->data[r][c-1];
                double zx1 = _z->data[r][c-1];

                delta_h1 = hx - hx1;
                delta_u1 = ux - ux1;
                delta_v1 = vx - vx1;
                dz1 = zx - zx1;
            }
            if(c < _nrCols-1 && !MV(r,c+1)) {

                double hx2 = _h->data[r][c+1];
                double ux2 = _u->data[r][c+1];
                double vx2 = _v->data[r][c+1];
                double zx2 = _z->data[r][c+1];

                delta_h2 = hx2 - hx;
                delta_u2 = ux2 - ux;
                delta_v2 = vx2 - vx;
                dz2 = zx2 - zx;
            }

            double dh   = 0.5*limiter(delta_h1, delta_h2);
            double dz_h = 0.5*limiter(delta_h1 + dz1, delta_h2 + dz2);

            double du   = 0.5*limiter(delta_u1, delta_u2);
            double dv   = 0.5*limiter(delta_v1, delta_v2);

            double _h1r = hx+dh;
            double _h1l = hx-dh;

            double _z1r = zx+(dz_h-dh);
            double _z1l = zx+(dh-dz_h);

            double _delzc1 = _z1r-_z1l;

            double hlh = _h1l/hx;
            double hrh = _h1r/hx;

            double _u1r = ux + hlh * du;
            double _u1l = ux - hrh * du;
            double _v1r = vx + hlh * dv;
            double _v1l = vx - hrh * dv;

            h1r->Drc = _h1r;
            h1l->Drc = _h1l;
            z1r->Drc = _z1r;
            z1l->Drc = _z1l;
            delzc1->Drc = _delzc1;
            u1r->Drc = _u1r;
            u1l->Drc = _u1l;
            v1r->Drc = _v1r;
            v1l->Drc = _v1l;
        }
    }

//#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (_h->Drc > he_ca) {
            if (c > 0 && !MV(r,c-1))
                delz1->data[r][c-1] = z1l->Drc - z1r->data[r][c-1];
        }
    }

#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (_h->Drc > he_ca)
        {
            double delta_h1 = 0;
            double delta_u1 = 0;
            double delta_v1 = 0;
            double delta_h2 = 0;
            double delta_u2 = 0;
            double delta_v2 = 0;
            double dz1 = 0;
            double dz2 = 0;
            double hx = _h->Drc;
            double ux = _u->Drc;
            double vx = _v->Drc;
            double zx = _z->Drc;

            if(r > 0 && !MV(r-1,c)) {
                double hx1 = _h->data[r-1][c];
                double ux1 = _u->data[r-1][c];
                double vx1 = _v->data[r-1][c];
                double zx1 = _z->data[r-1][c];
                delta_h1 = hx - hx1;
                delta_u1 = ux - ux1;
                delta_v1 = vx - vx1;
                dz1 = zx - zx1;
            }
            if(r < _nrRows-1 && !MV(r+1, c)) {
                double hx2 = _h->data[r+1][c];
                double ux2 = _u->data[r+1][c];
                double vx2 = _v->data[r+1][c];
                double zx2 = _z->data[r+1][c];
                delta_h2 = hx2 - hx;
                delta_u2 = ux2 - ux;
                delta_v2 = vx2 - vx;
                dz2 = zx2 - zx;
            }
            double dh   = 0.5*limiter(delta_h1, delta_h2);
            double dz_h = 0.5*limiter(delta_h1+dz2,delta_h2+dz2);
            double du   = 0.5*limiter(delta_u1, delta_u2);
            double dv   = 0.5*limiter(delta_v1, delta_v2);

            double _h2r = hx+dh;
            double _h2l = hx-dh;

            double _z2r = zx+(dz_h-dh);
            double _z2l = zx+(dh-dz_h);

            double _delzc2 = _z2r-_z2l;

            double hlh = _h2l/hx;
            double hrh = _h2r/hx;

            double _u2r = ux + hlh * du;
            double _u2l = ux - hrh * du;
            double _v2r = vx + hlh * dv;
            double _v2l = vx - hrh * dv;

            h2r->Drc = _h2r;
            h2l->Drc = _h2l;
            z2r->Drc = _z2r;
            z2l->Drc = _z2l;
            delzc2->Drc = _delzc2;
            u2r->Drc = _u2r;
            u2l->Drc = _u2l;
            v2r->Drc = _v2r;
            v2l->Drc = _v2l;
        }
    }

//#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (_h->Drc > he_ca) {
            if(r > 0 && !MV(r-1,c))
                delz2->data[r-1][c] = z2l->Drc - z2r->data[r-1][c];
        }
    }
}
//---------------------------------------------------------------------------

double TWorld::maincalcfluxOF(cTMap *_h,double dt, double dt_max)
{
    double dt_tmp, dtx, dty;
    cTMap *fbw = FlowBarrierW;
    cTMap *fbe = FlowBarrierE;
    cTMap *fbn = FlowBarrierN;
    cTMap *fbs = FlowBarrierS;

    //#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double _f1=0;
        double _f2=0;
        double _f3=0;
        double _cflx=0;
        double _f1o=0;
        double _f2o=0;
        double _f3o=0;
        double _h1d = 0;
        double _h1g = 0;
        vec4 rec;
        double _h1r = h1r->Drc;
        double _h1l = h1l->Drc;
        double _fbe = fbe->Drc;
        double _fbw = fbw->Drc;
        double _u1r = u1r->Drc;
        double _u1l = u1l->Drc;
        double _v1r = v1r->Drc;
        double _v1l = v1l->Drc;

        if(_h->Drc > he_ca)
        {
            if(c > 0 && !MV(r,c-1)) {
                _h1r = h1r->data[r][c-1];
                double _delz1 = delz1->data[r][c-1];
                _fbe = fbe->data[r][c-1];
                _u1r = u1r->data[r][c-1];
                _v1r = v1r->data[r][c-1];
                _h1d = std::max(0.0, _h1r - std::max(0.0,  _delz1  + std::max(_fbw,_fbe)));
                _h1g = std::max(0.0, _h1l - std::max(0.0, -_delz1  + std::max(_fbw,_fbe)));

                rec = F_Riemann(_h1d, _u1r, _v1r,_h1g, _u1l, _v1l);
                _f1 =   rec.v[0];
                _f2 =   rec.v[1];
                _f3 =   rec.v[2];
                _cflx = rec.v[3];
            }

            // left hand side boundary
            if (c == 0 || MV(r,c-1))
            {
//                _h1l = h1l->Drc;
//                _fbe = fbe->Drc;
//                _u1l = u1l->Drc;
//                _v1l = v1l->Drc;
                _h1g = std::max(0.0, _h1l - _fbe);

                rec = F_Riemann(0,0,0, _h1g, _u1l, _v1l);
                _f1 = rec.v[0];
                _f2 = rec.v[1];
                _f3 = rec.v[2];
                _cflx = rec.v[3];
            }

            // right hand side boundary
            if(c == _nrCols-1 || MV(r, c+1)) {
//                _h1r = h1r->Drc;
//                _fbw = fbw->Drc;
//                _u1r = u1r->Drc;
//                _v1r = v1r->Drc;
                _h1d = std::max(0.0, _h1r - _fbw);

                rec = F_Riemann(_h1d,_u1r,_v1r,0.,0.,0.);
                _f1o = rec.v[0];
                _f2o = rec.v[1];
                _f3o = rec.v[2];
            }
        }
        f1->Drc = _f1;
        f2->Drc = _f2;
        f3->Drc = _f3;
        cflx->Drc = _cflx;
        f1o->Drc = _f1o;
        f2o->Drc = _f2o;
        f3o->Drc = _f3o;

        if(c > 0 && !MV(r,c-1))
            h1d->data[r][c-1] = _h1d;
        h1g->Drc = _h1g;
    }

    //#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double _g1=0;
        double _g2=0;
        double _g3=0;
        double _cfly=0;
        double _g1o=0;
        double _g2o=0;
        double _g3o=0;
        vec4 rec;
        double _h2d = 0;
        double _h2g = 0;
        double _h2r = h2r->Drc;
        double _h2l = h2l->Drc;
        double _fbn = fbn->Drc;
        double _fbs = fbs->Drc;
        double _u2r = u2r->Drc;
        double _u2l = u2l->Drc;
        double _v2r = v2r->Drc;
        double _v2l = v2l->Drc;


        if(_h->Drc > he_ca)
        {
            if(r > 0 && !MV(r-1,c)) {
                _h2r = h2r->data[r-1][c];
                 double _delz2 = delz2->data[r-1][c];
                _fbn = fbn->data[r-1][c];
                _u2r = u2r->data[r-1][c];
                _v2r = v2r->data[r-1][c];
                _h2d = std::max(0.0, _h2r - std::max(0.0,  _delz2  + std::max(_fbs,_fbn)));
                _h2g = std::max(0.0, _h2l - std::max(0.0, -_delz2  + std::max(_fbs,_fbn)));

                rec = F_Riemann(_h2d,_v2r,_u2r, _h2g,_v2l,_u2l);
                _g1 = rec.v[0];
                _g2 = rec.v[2]; //!!!!!!!
                _g3 = rec.v[1]; //!!!!!!!
                _cfly = rec.v[3];
            }

            // upper side boundary
            if (r == 0 || MV(r-1,c))
            {
                _h2g = std::max(0.0, _h2l - _fbn);

                rec = F_Riemann(0,0,0,_h2g,_v2l,_u2l);
                _g1 = rec.v[0];
                _g2 = rec.v[2];
                _g3 = rec.v[1];
                _cfly = rec.v[3];
            }
            // lower side boundary
            if (r == _nrRows-1 || MV(r+1, c)) {
                _h2d = std::max(0.0, _h2d - _fbs);

                rec = F_Riemann(_h2d,_v2l,_u2l,0.,0.,0.);
                _g1o = rec.v[0];
                _g2o = rec.v[2];
                _g3o = rec.v[1];
            }
        }

        g1->Drc = _g1;
        g2->Drc = _g2;
        g3->Drc = _g3;
        cfly->Drc = _cfly;
        g1o->Drc = _g1o;
        g2o->Drc = _g2o;
        g3o->Drc = _g3o;
        if(r > 0 && !MV(r-1,c))
            h2d->data[r-1][c] = _h2d;
        h2g->Drc = _h2g;

    }

    dtx = dt_max;
    dty = dt_max;
//    #pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV
            if (_h->Drc > he_ca)
    {
        double dx = _dx;//ChannelAdj->Drc;//FlowWidth->Drc;//
        if (qFabs(cflx->Drc*dt/dx) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dx/cflx->Drc;
        dtx = std::min(std::min(dt, dt_tmp), dtx);
//        dtx = std::min(dt_tmp, dtx);
    }

//#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV
            if (_h->Drc > he_ca)
    {
        double dy = _dx;//DX->Drc;
        if (qFabs(cfly->Drc*dt/dy) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dy/cfly->Drc;
        dty = std::min(std::min(dt, dt_tmp), dty);
//        dty = std::min(dt_tmp, dty);
    }

    return(std::max(TimestepfloodMin, std::min(dtx,dty)));
}

void TWorld::maincalcschemeOF(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,cTMap *hes, cTMap *ves1, cTMap *ves2)
{
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Hes, Ves1, Ves2;
        double tx = dt/ChannelAdj->Drc;
        double ty = dt/DX->Drc;
        double _f1=0, _f2=0, _f3=0, _g1=0, _g2=0, _g3=0;

        //choose left hand boundary and normal (f1), or right hand boundary values (f1o)
        if (c < _nrCols-1 && !MV(r, c+1)){
            _f1 = f1->data[r][c+1];
            _f2 = f2->data[r][c+1];
            _f3 = f3->data[r][c+1];
        }
        if (c == _nrCols-1 || MV(r, c+1)){
            _f1 = f1o->Drc;
            _f2 = f2o->Drc;
            _f3 = f3o->Drc;
        }
        if (r < _nrRows-1 && !MV(r+1, c)) {
            _g1 = g1->data[r+1][c];
            _g2 = g2->data[r+1][c];
            _g3 = g3->data[r+1][c];
        }
        if (r == _nrRows-1 || MV(r+1, c)) {
            _g1 = g1o->Drc;
            _g2 = g2o->Drc;
            _g3 = g3o->Drc;
        }

//        double C = 0.25;
//        double h_x1 =  c > 0 && !MV(r,c-1)         ? he->data[r][c-1] : he->Drc;
//        double h_x2 =  c < _nrCols-1 && !MV(r,c+1) ? he->data[r][c+1] : he->Drc;
//        double h_y1 =  r > 0 && !MV(r-1,c)         ? he->data[r-1][c] : he->Drc;
//        double h_y2 =  r < _nrRows-1 && !MV(r+1,c) ? he->data[r+1][c] : he->Drc;
//        double flux_x1 = std::max(C*-he->Drc, std::min(-tx*_f1,      C*h_x2));
//        double flux_x2 = std::max(C*-he->Drc, std::min(+tx* f1->Drc, C*h_x1));
//        double flux_y1 = std::max(C*-he->Drc, std::min(-ty*_g1,      C*h_y2));
//        double flux_y2 = std::max(C*-he->Drc, std::min(+ty* g1->Drc, C*h_y1));

        //Hes = std::max(0.0, he->Drc + flux_x1 + flux_x2 + flux_y1 + flux_y2);
        Hes = std::max(0.0, he->Drc - tx*(_f1 - f1->Drc) - ty*(_g1 - g1->Drc));

        if (Hes > he_ca)
        {
            //Solution of the equation of momentum (Second and third equation of Saint-venant)
            double qes1;
            double qes2;

            qes1 = he->Drc*ve1->Drc -
                    ty*(_g2 - g2->Drc) -
                    tx*(_f2 - f2->Drc +
                        GRAV*0.5*((h1g->Drc-h1l->Drc)*(h1g->Drc+h1l->Drc) +
                                  (h1r->Drc-h1d->Drc)*(h1r->Drc+h1d->Drc)
                                  + (h1l->Drc+h1r->Drc)*delzc1->Drc)) ;

            qes2 = he->Drc*ve2->Drc -
                    tx*(_f3 - f3->Drc) -
                    ty*(_g3 - g3->Drc +
                        GRAV*0.5*((h2g->Drc-h2l->Drc)*(h2g->Drc+h2l->Drc) +
                                  (h2r->Drc-h2d->Drc)*(h2r->Drc+h2d->Drc)
                                  + (h2l->Drc+h2r->Drc)*delzc2->Drc));


            double sqUV = qSqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc);
            double nsq1 = (0.001+N->Drc)*(0.001+N->Drc)*GRAV/std::max(0.01, qPow(Hes,4.0/3.0));
            double nsq = nsq1*sqUV*dt;

            Ves1 = (qes1/(1.0+nsq))/std::max(0.01, Hes);
            Ves2 = (qes2/(1.0+nsq))/std::max(0.01, Hes);

            if (SwitchTimeavgV) {
                double fac = 0;
                fac = 0.5+0.5*std::min(1.0,4*Hes)*std::min(1.0,4*Hes);
                fac = fac *exp(- std::max(1.0,dt) / nsq1);
                Ves1 = fac * ve1->Drc + (1.0-fac) *Ves1;
                Ves2 = fac * ve2->Drc + (1.0-fac) *Ves2;
            }

            //       correctSpuriousVelocities(r, c, hes, ves1, ves2);
            // gives instability!

//            double threshold = 0.001 * _dx;
//            if(Hes < threshold) {
//                double h23 = pow(Hes, 2.0/3.0);//Hes*sqrt(Hes);
//                double kinfac = std::max(0.0,(threshold - Hes) / (0.025 * _dx));
//                double sx_zh = delz1->Drc;
//                double sy_zh = delz2->Drc;
//                double v_kin = (sx_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sx_zh > 0 ? sx_zh : -sx_zh))/(0.001+N->Drc);
//                Ves1 = kinfac * v_kin + Ves1*(1.0-kinfac);
//                v_kin = (sy_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sy_zh > 0 ? sy_zh : -sy_zh))/(0.001+N->Drc);
//                Ves2 = kinfac * v_kin + Ves2*(1.0-kinfac);
//            }

            double vmax = 0.25*_dx/dt;
            Ves1 = std::max(-vmax, std::min(vmax, Ves1));
            Ves2 = std::max(-vmax, std::min(vmax, Ves2));
        }
        else
        {
            Hes = 0;
            Ves1 = 0;
            Ves2 = 0;
        }

        // dan maar even met geweld!
        if (std::isnan(Ves1) || std::isnan(Ves2))
        {
            //   qDebug() << "Ves nan" << Ves1 << Ves2;
            Ves1 = 0;
            Ves2 = 0;
            Hes = 0;
        }

        hes->Drc = Hes;
        ves1->Drc = Ves1;
        ves2->Drc = Ves2;
    }
}
//---------------------------------------------------------------------------
double TWorld::fullSWOF2RO(cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    double dt_max = std::min(_dt, _dx*0.5);
    double dt1 = dt_max, timesum = 0;
    int n = 0;
    double sumh = 0;
    bool stop;


    if (startFlood)
    {
        sumh = getMass(h, 0);

        do {

         //   dt1 = dt_max;

            setZeroOF(h, u, v);
            simpleSchemeOF(h,u,v);
            // MUSCL: build left and right pressures anbd velocities with flux limiters
            if (SwitchMUSCL)
                MUSCLOF(h,u,v,z);

            // non openmp version
            //if (SwitchMUSCL)
            //MUSCL(h,u,v,z);
            //dt1 = maincalcflux(h, dt1, dt_max);
            //dt1 = std::min(dt1, _dt-timesum);
            //maincalcscheme(dt1, h,u,v, hs,us,vs);

            // riemann solvers
            dt1 = maincalcfluxOF(h, dt1, dt_max);
            dt1 = std::min(dt1, _dt-timesum);
            // st venant equations
            maincalcschemeOF(dt1, h,u,v, hs,us,vs);

            // for erosion
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                FloodDT->Drc = dt1;
            }

            if (SwitchErosion)
                SWOFSediment(FloodDT,h,u,v);

            setZeroOF(hs, us, vs);
#pragma omp parallel for collapse(2) num_threads(userCores)
            FOR_ROW_COL_MV_L {
                h->Drc = hs->Drc;
                u->Drc = us->Drc;
                v->Drc = vs->Drc;
            }


            timesum = timesum + dt1;
            stop = timesum  > _dt-1e-6;
            n++;

            if (n > F_MaxIter)
                break;

        } while (!stop);
    } // if floodstart

    correctMassBalance(sumh, h, 0);

    iter_n = n;
    dt1 = n > 0? _dt/n : dt1;

    return(dt1);
}


