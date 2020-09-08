
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-6
#define ve_ca 1e-6

#define GRAV 9.8067
#define EPSILON 1e-6
//--------------------------------------------------------------------------------------------
// correct mass balance
double TWorld::getMass(cTMap *M)
{
    double sum2 = 0;
    #pragma omp parallel for reduction(+:sum2) collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
    }
    return sum2;
}
//---------------------------------------------------------------------------
// correct mass balance
void TWorld::correctMassBalance(double sum1, cTMap *M)
{
    double sum2 = 0;
    double n = 0;

    #pragma omp parallel for reduction(+:sum2) collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        if(M->Drc > 0)
        {
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
            if(M->Drc > 0)
                n += 1;
        }
    }
    // total and cells active for M

    //double dh = (n > 0 ? (sum1 - sum2)/n : 0);
    double dhtot = sum2 > 0 ? (sum1 - sum2)/sum2 : 0;
   // qDebug() << 1+dhtot;
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        if(M->Drc > 0)
        {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            //M->Drc += dh/(DX->Drc*ChannelAdj->Drc); // <- equal distribution error
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

    double U1 =  4.0*(pow(hes->Drc,2.0/3.0)*sqrt(s1)/N->Drc + G);
    double V1 =  4.0*(pow(hes->Drc,2.0/3.0)*sqrt(s2)/N->Drc + G);

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
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h, hlh, hrh, dz1, dz2;
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (_h->Drc > he_ca) {
            if(c > 0 && !MV(r,c-1) && c < _nrCols-1 && !MV(r,c+1)) {
                delta_h1 = _h->Drc - _h->data[r][c-1];
                delta_u1 = _u->Drc - _u->data[r][c-1];
                delta_v1 = _v->Drc - _v->data[r][c-1];
                dz1 = _z->Drc - _z->data[r][c-1];
                delta_h2 = _h->data[r][c+1] - _h->Drc;
                delta_u2 = _u->data[r][c+1] - _u->Drc;
                delta_v2 = _v->data[r][c+1] - _v->Drc;
                dz2 = _z->data[r][c+1] - _z->Drc;
            } else {
                delta_h1 = 0;
                delta_u1 = 0;
                delta_v1 = 0;
                delta_h2 = 0;
                delta_u2 = 0;
                delta_v2 = 0;
                dz1 = 0;
                dz2 = 0;
            }
            dh   = 0.5*limiter(delta_h1, delta_h2);
            dz_h = 0.5*limiter(delta_h1 + dz1, delta_h2 + dz2);

            du   = 0.5*limiter(delta_u1, delta_u2);
            dv   = 0.5*limiter(delta_v1, delta_v2);

            h1r->Drc = _h->Drc+dh;
            h1l->Drc = _h->Drc-dh;

            z1r->Drc = _z->Drc+(dz_h-dh);
            z1l->Drc = _z->Drc+(dh-dz_h);

            delzc1->Drc = z1r->Drc-z1l->Drc;
            if (c > 0 && !MV(r,c-1))
                delz1->data[r][c-1] = z1l->Drc - z1r->data[r][c-1];

            if (_h->Drc > he_ca) {
                hlh = h1l->Drc/_h->Drc;
                hrh = h1r->Drc/_h->Drc;
            }
            else {
                hlh = 1.0;
                hrh = 1.0;
            }

            u1r->Drc = _u->Drc + hlh * du;
            u1l->Drc = _u->Drc - hrh * du;
            v1r->Drc = _v->Drc + hlh * dv;
            v1l->Drc = _v->Drc - hrh * dv;
        }
    }
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (_h->Drc > he_ca) {
            if(r > 0 && MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c)) {
                delta_h1 = _h->Drc - _h->data[r-1][c];
                delta_u1 = _u->Drc - _u->data[r-1][c];
                delta_v1 = _v->Drc - _v->data[r-1][c];
                dz1 = _z->Drc - _z->data[r-1][c];
                delta_h2 = _h->data[r+1][c] - _h->Drc;
                delta_u2 = _u->data[r+1][c] - _u->Drc;
                delta_v2 = _v->data[r+1][c] - _v->Drc;
                dz2 = _z->data[r+1][c] - _z->Drc;
            } else {
                delta_h1 = 0;
                delta_u1 = 0;
                delta_v1 = 0;
                delta_h2 = 0;
                delta_u2 = 0;
                delta_v2 = 0;
                dz1 = 0;
                dz2 = 0;
            }
            dh   = 0.5*limiter(delta_h1, delta_h2);
            dz_h = 0.5*limiter(delta_h1+dz2,delta_h2+dz2);
            du   = 0.5*limiter(delta_u1, delta_u2);
            dv   = 0.5*limiter(delta_v1, delta_v2);

            h2r->Drc = _h->Drc+dh;
            h2l->Drc = _h->Drc-dh;

            z2r->Drc = _z->Drc+(dz_h-dh);
            z2l->Drc = _z->Drc+(dh-dz_h);

            delzc2->Drc = z2r->Drc - z2l->Drc;
            if(r > 0 && MV(r-1,c))
                delz2->data[r-1][c] = z2l->Drc - z2r->data[r-1][c];

            if (_h->Drc > he_ca) {
                hlh = h2l->Drc/_h->Drc;
                hrh = h2r->Drc/_h->Drc;
            }
            else {
                hlh = 1.0;
                hrh = 1.0;
            }

            u2r->Drc = _u->Drc + hlh * du;
            u2l->Drc = _u->Drc - hrh * du;
            v2r->Drc = _v->Drc + hlh * dv;
            v2l->Drc = _v->Drc - hrh * dv;
        }
    }
}
//---------------------------------------------------------------------------

double TWorld::maincalcfluxOF(cTMap *_h,double dt, double dt_max)
{
    vec4 rec;
    double dt_tmp, dtx, dty;
    cTMap *fbw = FlowBarrierW;
    cTMap *fbe = FlowBarrierE;
    cTMap *fbn = FlowBarrierN;
    cTMap *fbs = FlowBarrierS;

//#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(_h->Drc > he_ca){
            f1->Drc = 0;
            f2->Drc = 0;
            f3->Drc = 0;
            f1o->Drc = 0;
            f2o->Drc = 0;
            f3o->Drc = 0;

            if(c > 0 && !MV(r,c-1)) {
                h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                rec = F_Riemann(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
                f1->Drc =   rec.v[0];
                f2->Drc =   rec.v[1];
                f3->Drc =   rec.v[2];
                cflx->Drc = rec.v[3];
            } else {
                double _h1g = std::max(0.0, h1l->Drc - fbe->Drc);
                rec = F_Riemann(0,0,0, _h1g, u1l->Drc, v1l->Drc);
                f1->Drc = rec.v[0];
                f2->Drc = rec.v[1];
                f3->Drc = rec.v[2];
                cflx->Drc = rec.v[3];
            }

            // right hand side boundary
            if(c == _nrCols-1 || MV(r, c+1)){
                double _h1d = std::max(0.0, h1r->Drc - fbw->Drc);
                rec = F_Riemann(_h1d,u1r->Drc,v1r->Drc,0.,0.,0.);
                f1o->Drc = rec.v[0];
                f2o->Drc = rec.v[1];
                f3o->Drc = rec.v[2];
            }
        }
    }

//#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(_h->Drc > he_ca){
            g1->Drc = 0;
            g2->Drc = 0;
            g3->Drc = 0;
            g1o->Drc = 0;
            g2o->Drc = 0;
            g3o->Drc = 0;

            if(r > 0 && !MV(r-1,c)) {
                h2d->data[r-1][c] = std::max(0.0, h2r->data[r-1][c] - std::max(0.0,  delz2->data[r-1][c]  + std::max(fbs->Drc,fbn->data[r-1][c])));
                h2g->Drc          = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[r-1][c]  + std::max(fbs->Drc,fbn->data[r-1][c])));

                rec = F_Riemann(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c], h2g->Drc,v2l->Drc,u2l->Drc);

                g1->Drc = rec.v[0];
                g2->Drc = rec.v[2];
                g3->Drc = rec.v[1];
                cfly->Drc = rec.v[3];
            } else {
                double _h2g = std::max(0.0, h2l->Drc - fbn->Drc);
                rec = F_Riemann(0,0,0,_h2g,v2l->Drc,u2l->Drc);
                g1->Drc = rec.v[0];
                g2->Drc = rec.v[2];
                g3->Drc = rec.v[1];
                cfly->Drc = rec.v[3];
            }
            // left hand side boundary
            if (r == _nrRows-1 || MV(r+1, c)) {
                double _h2d = std::max(0.0, h2d->Drc - fbs->Drc);
                rec = F_Riemann(_h2d,v2l->Drc,u2l->Drc,0.,0.,0.);
                g1o->Drc = rec.v[0];
                g2o->Drc = rec.v[2];
                g3o->Drc = rec.v[1];
            }
        }
    }

    dtx = dt_max;
    dty = dt_max;
    FOR_ROW_COL_MV
            if (_h->Drc > he_ca)
    {
        double dx = _dx;//ChannelAdj->Drc;//FlowWidth->Drc;//
        if (qFabs(cflx->Drc*dt/dx) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dx/cflx->Drc;
        dtx = std::min(std::min(dt, dt_tmp), dtx);
    }

    FOR_ROW_COL_MV
            if (_h->Drc > he_ca)
    {
        double dy = _dx;//DX->Drc;
        if (qFabs(cfly->Drc*dt/dy) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dy/cfly->Drc;
        dty = std::min(std::min(dt, dt_tmp), dty);
    }

    return(std::max(TimestepfloodMin, std::min(dtx,dty)));
}

void TWorld::maincalcschemeOF(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,cTMap *hes, cTMap *ves1, cTMap *ves2)
{
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Hes, Ves1, Ves2;
        double tx = dt/_dx;//ChannelAdj->Drc; //FlowWidth->Drc;
        double ty = dt/_dx;//DX->Drc;
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
        else if (r == _nrRows-1 || MV(r+1, c)) {
            _g1 = g1o->Drc;
            _g2 = g2o->Drc;
            _g3 = g3o->Drc;
        }

        Hes = std::max(0.0, he->Drc - tx*(_f1 - f1->Drc) - ty*(_g1 - g1->Drc));
        //hes->Drc = std::max(0.0, he->Drc - tx*(_f1 - f1->Drc) - ty*(_g1 - g1->Drc));

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

            Ves1 = (qes1/(1.0+nsq))/Hes;
            Ves2 = (qes2/(1.0+nsq))/Hes;

            //NOTE ves1 is c-1, c+1    ves2 is r-1, r+1
            if (SwitchTimeavgV) {
                double fac = 0.5+0.5*std::min(1.0,4*Hes)*std::min(1.0,4*Hes);
                fac = fac *exp(- std::max(1.0,dt) / nsq1);
                Ves1 = fac * ve1->Drc + (1.0-fac) *Ves1;
                Ves2 = fac * ve2->Drc + (1.0-fac) *Ves2;
            }

            //correctSpuriousVelocities(r, c, hes, ves1, ves2);


            double threshold = 0.01 * _dx; // was 0.01
            if(Hes < threshold) {
            //    correctSpuriousVelocities(r, c, hes, ves1, ves2);
                double h23 = Hes*sqrt(Hes);//pow(Hes, 2.0/3.0);//hn * sqrt(hn)
                double kinfac = std::max(0.0,(threshold - Hes) / (0.025 * _dx));
                double sx_zh = delz1->Drc;
                double sy_zh = delz2->Drc;
                double v_kin = (sx_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sx_zh > 0 ? sx_zh : -sx_zh))/(0.001+N->Drc);
                Ves1 = kinfac * v_kin + Ves1*(1.0-kinfac);
                v_kin = (sy_zh>0?1:-1) * h23 * std::max(0.001, sqrt(sy_zh > 0 ? sy_zh : -sy_zh))/(0.001+N->Drc);
                Ves2 = kinfac * v_kin + Ves2*(1.0-kinfac);
            }

            double vmax = 0.5*_dx/dt;
            Ves1 = std::max(-vmax, std::min(vmax, Ves1));
            Ves2 = std::max(-vmax, std::min(vmax, Ves2));
        }
        else
        {
            // Case of height of water < ha.
            Ves1 = 0;
            Ves2 = 0;
        }
        // dan maar even met geweld!
        if (std::isnan(Ves1) || std::isnan(Ves2)  )
        {
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
    double dt1 = 0, timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int n = 0;
    double sumh = 0;
    bool stop;

    if (startFlood)
    {
        sumh = getMass(h);

        do {

            dt1 = dt_max;

            setZeroOF(h, u, v);
            simpleSchemeOF(h,u,v);
            if (SwitchMUSCL)
                MUSCLOF(h,u,v,z);

            dt1 = maincalcfluxOF(h, dt1, dt_max);
            dt1 = std::min(dt1, _dt-timesum);

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

    correctMassBalance(sumh, h);

    iter_n = n;
    dt1 = n > 0? _dt/n : dt1;

    return(dt1);
}


