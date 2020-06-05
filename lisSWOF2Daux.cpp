
#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-6
#define ve_ca 1e-6

#define GRAV 9.8067
#define EPSILON 1e-6


void TWorld::simpleSchemeOF(cTMap *_h,cTMap *_u,cTMap *_v)
{
    FOR_ROW_COL_MV{
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
    FOR_ROW_COL_MV  {
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

void TWorld::MUSCLOF(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z)
{
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h, hlh, hrh, dz1, dz2;

    FOR_ROW_COL_MV {
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

    FOR_ROW_COL_MV {
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

double TWorld::maincalcfluxOF(cTMap *_h,double dt, double dt_max)
{
    double dt_tmp, dtx, dty;
    cTMap *fbw = FlowBarrierW;
    cTMap *fbe = FlowBarrierE;
    cTMap *fbn = FlowBarrierN;
    cTMap *fbs = FlowBarrierS;

    FOR_ROW_COL_MV {
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
                F_HLL2(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
                f1->Drc = HLL2_f1;
                f2->Drc = HLL2_f2;
                f3->Drc = HLL2_f3;
                cflx->Drc = HLL2_cfl;
            } else {
                double _h1g = std::max(0.0, h1l->Drc - FlowBarrierE->Drc);
                F_HLL2(0,0,0, _h1g, u1l->Drc, v1l->Drc);
                f1->Drc = HLL2_f1;
                f2->Drc = HLL2_f2;
                f3->Drc = HLL2_f3;
                cflx->Drc = HLL2_cfl;
            }

            // right hand side boundary
            if(c == _nrCols-1 || MV(r, c+1)){
                double _h1d = std::max(0.0, h1r->Drc - fbw->Drc);
                F_HLL2(_h1d,u1r->Drc,v1r->Drc,0.,0.,0.);
                f1o->Drc = HLL2_f1;
                f2o->Drc = HLL2_f2;
                f3o->Drc = HLL2_f3;
            }
       }
    }

    FOR_ROW_COL_MV {
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
                F_HLL2(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c], h2g->Drc,v2l->Drc,u2l->Drc);
                g1->Drc = HLL2_f1;
                g2->Drc = HLL2_f3;
                g3->Drc = HLL2_f2;
                cfly->Drc = HLL2_cfl;
            } else {
                double _h2g = std::max(0.0, h2l->Drc - fbn->Drc);
                F_HLL2(0,0,0,_h2g,v2l->Drc,u2l->Drc);
                g1->Drc = HLL2_f1;
                g2->Drc = HLL2_f2;
                g3->Drc = HLL2_f3;
                cflx->Drc = HLL2_cfl;
            }
            // left hand side boundary
            if (r == _nrRows-1 || MV(r+1, c)) {
                double _h2d = std::max(0.0, h2d->Drc - fbs->Drc);
                F_HLL2(_h2d,v2l->Drc,u2l->Drc,0.,0.,0.);
                g1o->Drc = HLL2_f1;
                g2o->Drc = HLL2_f3;
                g3o->Drc = HLL2_f2;
            }
        }
    }

    dtx = dt_max;
    dty = dt_max;
    FOR_ROW_COL_MV
            if (_h->Drc > he_ca)
    {
        double dx = ChannelAdj->Drc;//FlowWidth->Drc;//
        if (qFabs(cflx->Drc*dt/dx) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dx/cflx->Drc;
        dtx = std::min(std::min(dt, dt_tmp), dtx);
    }

    FOR_ROW_COL_MV
            if (_h->Drc > he_ca)
    {
        double dy = DX->Drc;
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
    FOR_ROW_COL_MV {
        double tx = dt/ChannelAdj->Drc; //FlowWidth->Drc;
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
        else if (r == _nrRows-1 || MV(r+1, c)) {
            _g1 = g1o->Drc;
            _g2 = g2o->Drc;
            _g3 = g3o->Drc;
        }
        hes->Drc = std::max(0.0, he->Drc - tx*(_f1 - f1->Drc) - ty*(_g1 - g1->Drc));

        if (hes->Drc > he_ca)
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

            //

            qes2 = he->Drc*ve2->Drc -
                    tx*(_f3 - f3->Drc) -
                    ty*(_g3 - g3->Drc +
                        GRAV*0.5*((h2g->Drc-h2l->Drc)*(h2g->Drc+h2l->Drc) +
                                  (h2r->Drc-h2d->Drc)*(h2r->Drc+h2d->Drc)
                                  + (h2l->Drc+h2r->Drc)*delzc2->Drc));

            // double sqQ = qSqrt(qes1*qes1+qes2*qes2);

            double sqUV = qSqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc);
            double nsq1 = (0.001+N->Drc)*(0.001+N->Drc)*GRAV/qPow(hes->Drc,4.0/3.0);
            double nsq = nsq1*sqUV*dt;

            ves1->Drc = (qes1/(1.0+nsq))/hes->Drc;
            ves2->Drc = (qes2/(1.0+nsq))/hes->Drc;

            //NOTE ves1 is c-1, c+1    ves2 is r-1, r+1

            correctSpuriousVelocities(r, c, hes, ves1, ves2);

            double fac = 0;
            if (SwitchTimeavgV) {
                fac = 0.5+0.5*std::min(1.0,4*hes->Drc)*std::min(1.0,4*hes->Drc);
                fac = fac *exp(- std::max(1.0,dt) / nsq1);
            }
            ves1->Drc = fac * ve1->Drc + (1.0-fac) *ves1->Drc;
            ves2->Drc = fac * ve2->Drc + (1.0-fac) *ves2->Drc;

        }
        else
        {
            // Case of height of water < ha.
            ves1->Drc = 0;
            ves2->Drc = 0;
        }
        // dan maar even met geweld!
        if (std::isnan(ves1->Drc) || std::isnan(ves2->Drc)  )
        {
            ves1->Drc = 0;
            ves2->Drc = 0;
            hes->Drc = 0;
        }
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

    //F_MaxIter = (int) (_dt/TimestepfloodMin);

    //SwitchHeun = false;

    if (startFlood)
    {
//        if (SwitchErosion) {
//            FOR_ROW_COL_MV {
//                SSFlood->Drc += DETSplash->Drc;
//                SSCFlood->Drc = MaxConcentration(ChannelAdj->Drc * DX->Drc * h->Drc, &SSFlood->Drc, &DepFlood->Drc);
//            }

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


//            if (SwitchErosion)
//                SWOFSediment(dt1,h,u,v);

            setZeroOF(hs, us, vs);
            FOR_ROW_COL_MV {
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


