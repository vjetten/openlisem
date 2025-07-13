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

//#include <algorithm>
#include "model.h"


//--------------------------------------------------------------------------------
/// new matrix head: Gaussian elimination and backsubstitution
/**
tri diagonalmatrix to solve differential equations to get the new matrix potential
and the new moisture content theta
method is by gaussian elimination and backsubstitution, 2 times instead of iteration
matrix shape:
\code
|b0 c0 .  .| |h0| |F0|
|a1 .. cn-1|*|  |=|  |
|   an bn  | |hn| |Fn|
\endcode
F (thomf) is new moisture content theta
thoma ?
thomc has unit cm
thomb has unit dtheta/dh (like C)
dz and disZ are negative distances
*/
/*!
  \file swatstep.cpp
  \brief SWATRE: calculation of new h for every pixel and infiltration/ponding

  functions:\n
- double TWorld::NewTimeStep(double prevDt,const double *hLast,const double *h,int nrNodes, double dtMin)\n
- void TWorld::ComputeForPixel(PIXEL_INFO *pixel, SOIL_MODEL *s, double drainfraction)\n
- void TWorld::calcSinktermSWATRE(PIXEL_INFO *pixel, double *h, double *S)\n
*/


// units in SWATRE are cm and K cm/sec

//--------------------------------------------------------------------------------
//void TWorld::calcSinktermSWATRE(PIXEL_INFO *pixel, QVector<double> h, QVector<double> S)
void TWorld::calcSinktermSWATRE(PIXEL_INFO *pixel, double *h, double *S)
{
    int r = pixel->r;
    int c = pixel->c;

    if (Rain->Drc*3600000.0/_dt > rainfallETa_threshold) {
        ETa->Drc = 0;
        ETp->Drc = 0;
        for (int j = 0; j < pixel->profile->zone->nrNodes; j++)
            S[j] = 0;
    }

    // ETafactor is calculated at model level, before hydrology
    if (ETp->Drc*ETafactor > 0) {

        double ETp_ = ETp->Drc * ETafactor;// * 100; // potential ETp in meter/day to cm/day!
        double tot = 0;
        double etanet = ETp_;
        const ZONE *zone = pixel->profile->zone;

        //transpiration under Cover from rootzone
        etanet = ETp_*(Cover->Drc)*(1-fractionImperm->Drc);
        for (int j = 0; j < zone->nrNodes; j++) {
            S[j] = 0;
            if (zone->rootz[j] > 0) {
                // van genuchten H50 = -3.5 m
                double f = 1.0/(1.0+pow(h[j]/-350,1.5));
                if (h[j] > -10) f = 0; // near saturation
                if (h[j] < -16000) f = 0; // wilting point -16000 cm
                S[j] =  etanet * f * zone->rootz[j];
            }
        }

        // add surface evaporation (1-Cover) to top node if no ponding
        etanet = ETp_*(1-fractionImperm->Drc);
        if (h[0] < -1) {
            double the = FindValue(h[0], pixel->profile->horizon[0], H_COL, THETA_COL);
            double theS = FindValue(0, pixel->profile->horizon[0], H_COL, THETA_COL);
            S[0] += etanet * the/theS;
        }

        for (int j = 0; j < zone->nrNodes; j++) {
            tot += S[j];
        }
//if(r == 200 && c == 200) qDebug() << "swatre" << tot << S[0] << S[1] << S[2];

        ETa->Drc = tot;
        ETaCum->Drc += tot;
    }
}
//--------------------------------------------------------------------------------
double TWorld::NewTimeStep(double prevDt,const double *hLast,const double *h,int nrNodes, double dtMin, double precParam)
//double TWorld::NewTimeStep(double prevDt,QVector <double> hLast,QVector <double> h,int nrNodes, double dtMin, double precParam)//
{
   // double precParam = SwatrePrecision;
    // note "5" is a precision factor determining next timestep, set to 5 in old lisem
    // higher gives better results!
   // precParam = 5;
    double dt = _dt;
    double accur1 = qMax(0.0, 0.3 - 0.02 * precParam);
    double accur2 = 0.1*accur1;//0.03 - 0.002 * precParam; //SwatrePrecision;//

    for(int i=0; i < nrNodes; i++) {
        double mdih = accur1 + accur2 * qMax(1.0, fabs(h[i]));
        double dih  = fabs(h[i] - hLast[i]);
        // if difference is small
        // dih = e.g. 10 and h = -200 then mdih = 200*0.01 + 0.1 = 2.1
        // mdih/dih = 2.1/10 =0.21

        if (dih > 0.10)
            dt = qMin(dt, prevDt*mdih/dih);
    }
    return (qMax(dt, dtMin));
}
//--------------------------------------------------------------------------------
// Units are:
// Z and H in cm; table units K in cm/day converted to cm/sec, lisem time in seconds
// NOTE: dz is negative, disZ is negative!

// #define dz(j) p->zone->dz[j]
// #define disZ(j) p->zone->disnod[j]

void TWorld::ComputeForPixel(long i_, SOIL_MODEL *s)//, NODES l)
{
    PIXEL_INFO *pixel = &s->pixel[i_];
    const PROFILE *p = pixel->profile;
    int r = pixel->r;
    int c = pixel->c;
    int nN = p->zone->nrNodes;

    double dt = swatreDT;//_dt/2;//SwatrePrecision; // start dt, precision set to 6 like in old code.
                                     // A higher precision gives less infiltration
    double WH = pixel->wh*100; // m to cm
    int tnode = pixel->tilenode;
    double impfrac = fractionImperm->Drc;

    double elapsedTime = 0;
    double drainout = 0;
    double percolation = 0;

    //  qDebug() << i_ << r << c << p->profileId;

    // DO NOT USE STATIC ARRAYS, they are shared bgetween threads in omp, so race conditions occuur, striping
    //NODE_ARRAY kavg, k, C, theta, thetaPrev, h, hPrev, dz, disZ, S;
    //NODE_ARRAY thoma, thomb, thomc, thomf, beta;

    // slower:
    // QVector <double> theta(MAX_NODES+3, 0.0);
    // QVector <double> kavg(MAX_NODES+3, 0.0);
    // QVector <double> k(MAX_NODES+3, 0.0);
    // QVector <double> C(MAX_NODES+3, 0.0);
    // QVector <double> thetaPrev(MAX_NODES+3, 0.0);
    // QVector <double> h(MAX_NODES+3, 0.0);
    // QVector <double> hPrev(MAX_NODES+3, 0.0);
    // QVector <double> dz(MAX_NODES+3, 0.0);
    // QVector <double> disZ(MAX_NODES+3, 0.0);
    // QVector <double> S(MAX_NODES+3, 0.0);
    // QVector <double> thoma(MAX_NODES+3, 0.0);
    // QVector <double> thomb(MAX_NODES+3, 0.0);
    // QVector <double> thomc(MAX_NODES+3, 0.0);
    // QVector <double> thomf(MAX_NODES+3, 0.0);
    // QVector <double> beta(MAX_NODES+3, 0.0);


    double* theta = new double[MAX_NODES+3]();
    double* kavg  = new double[MAX_NODES+3]();
    double* k     = new double[MAX_NODES+3]();
    double* C     = new double[MAX_NODES+3]();
    double* thetaPrev = new double[MAX_NODES+3]();
    double* h     = new double[MAX_NODES+3]();
    double* hPrev = new double[MAX_NODES+3]();
    double* dz    = new double[MAX_NODES+3]();
    double* disZ  = new double[MAX_NODES+3]();
    double* S     = new double[MAX_NODES+3]();
    double* thoma = new double[MAX_NODES+3]();
    double* thomb = new double[MAX_NODES+3]();
    double* thomc = new double[MAX_NODES+3]();
    double* thomf = new double[MAX_NODES+3]();
    double* beta  = new double[MAX_NODES+3]();

    for (int j = 0; j < nN; j++) {
      h[j] = pixel->h[j];
      dz[j] = p->zone->dz[j];
      disZ[j] = p->zone->disnod[j];
    }
    // memcpy(h, pixel->h.data(), nN * sizeof(double));
    // memcpy(dz, p->zone->dz.data(), nN * sizeof(double));
    // memcpy(disZ, p->zone->disnod.data(), nN * sizeof(double));

    if (SwitchIncludeET && WH <= 0) {
        calcSinktermSWATRE(pixel, h, S);
    }
    // get sinkterm S

    while (elapsedTime < _dt) {

        bool isPonded, fltsat;    // flag if ponded or if profile fully saturated
        double qmax, qtop, qbot, ThetaSat;  // fluxes at top and bottom, max theta
        double qdrain; // tile drainage

        // get nodal values of theta, K, dif moist cap
        for (int j = 0; j < nN; j++) {
            k[j] = FindValue(h[j], p->horizon[j], H_COL, K_COL);
            // K in cm/sec from h, ksatcal filled with values for ksat1,2,3
            C[j] = FindValue(h[j], p->horizon[j], DMCH_COL, DMCC_COL);
                    //DmcNode(h[j], p->horizon[j],  true); // true is more detailed method, false is DMCH directly from H
            // differential moisture capacity d(theta)/d(h), tangent moisture retention curve
            theta[j] = FindValue(h[j], p->horizon[j], H_COL, THETA_COL);
            // moisture content from H
        }

        // per pixel correction of Ks and Pore for org mat and density
        // near saturated so for h > -1 cm, and only for topsoil, assumed to be 30 cm
        if (SwitchOMCorrection) {
            for (int j = 0; j < nN && p->zone->endComp[j] <= 30 && h[j] > -10; j++) {
                k[j] = pixel->corrKsOA*k[j] + pixel->corrKsOB;
                theta[j] = pixel->corrPOA*theta[j] + pixel->corrPOB;
               // theta gives mass balance error because this decouples Theta from H?
            }
        }

        if (SwitchDensCorrection) {
            for (int j = 0; j < nN  && p->zone->endComp[j] <= 30 && h[j] > -10.0; j++) {
                k[j] = pixel->corrKsDA*k[j] + pixel->corrKsDB;
                theta[j] = pixel->corrPDA*theta[j] + pixel->corrPDB;
            }
        }

        // do calibration after dens and OM calculations
        for (int j = 0; j < nN; j++) {
             k[j] *= p->KsatCal[j];
        }

        // average K for 1st to n-1 node, top node is done below
        // original swatre artithmetric mean, Vauclin nin Belmans says geometric mean!
        // for(int j = 1; j < nN; j++) {
        //     kavg[j] = (k[j]+k[j-1])/2.0;
        //     //kavg[j] = sqrt(k[j]*k[j-1]);
        // }
        switch (KavgType) {
            case 0: for(int j = 1; j < nN; j++) { kavg[j] = Aavg(k[j],k[j-1]);} break;
            case 1: for(int j = 1; j < nN; j++) { kavg[j] = Savg(k[j],k[j-1]);} break;
            case 2: for(int j = 1; j < nN; j++) { kavg[j] = Havg(k[j],k[j-1],dz[j],dz[j-1]); }break;
            case 3: for(int j = 1; j < nN; j++) { kavg[j] = Mavg(k[j],k[j-1]);} break;
        }

        //--- boundary conditions ---

        //----- TOP -----
        // 1st check flux against max flux

        // max possible flux with Ksat
        double Ksat = FindValue(0, p->horizon[0], H_COL, K_COL)*p->KsatCal[0]*(1.0-impfrac);
        if (SwitchOMCorrection)
            Ksat = pixel->corrKsOA*Ksat + pixel->corrKsOB;
        if (SwitchDensCorrection)
            Ksat = pixel->corrKsDA*Ksat + pixel->corrKsDB;

        kavg[0] = sqrt(Ksat * k[0]);
        kavg[0] *= (1.0-impfrac);

        // adjust kavg[0] for roads and houses, impermeable fraction
        // max possible always geometric mean
        // geometric avg of ksat and k[0] => is used for max possible

        qmax = kavg[0]*(WH-h[0])/disZ[0] - kavg[0];
        // Darcy: = -k(dh/dz+1) = -kdh/dz-k
        // disZ is negative !!!

        // check if ponded: 1st compare fluxes, 2nd compare store
        qtop = -WH/dt;
        // top flux is water/timestep (cm/sec), negative downward
        // only for non impermeable surfaces. if more than 0.99 impermeable, swatstep is not done in infiltration()!
        isPonded = (qtop < qmax);
        // if more flux then max possible flag ponded is true. both are negative

        //2nd check: isPonded layer depth against storage
        if (!isPonded) {
            // calculate available space in profile in cm: (pore-theta)*dz
            double space = 0;
            for(int i = 0; i < nN && space < WH; i++) {
                ThetaSat = FindValue(0, p->horizon[i], H_COL, THETA_COL);
                // if (SwitchDensCorrection && p->zone->endComp[i] <= 30)
                //     ThetaSat = pixel->corrPDA*ThetaSat + pixel->corrPDB;
                space += (ThetaSat - theta[i]) * -dz[i];
            }
            isPonded = WH > space;
        }

        // check if profile is completely saturated (flstsat)
        fltsat = true;
        for (int i = nN-1; i >= 0; i--) {
            if (h[i] < 0) {
                fltsat = false;
                break;
            }
        }
        if (fltsat && (qtop <= qbot))
            fltsat = false;

        //----- BOTTOM -----
        // bottom is 0 or copy of flux of last 2 layers
        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = kavg[nN-1]*(h[nN-1]-h[nN-2])/disZ[nN-1] - kavg[nN-1];

        for (int j = 0; j < nN; j++) {
          hPrev[j] = h[j];
          thetaPrev[j] = theta[j];
        }
        //std::memcpy(hPrev, h, nN * sizeof(double));
        //std::memcpy(thetaPrev, theta, nN * sizeof(double));

        //HeadCalc(p, h, &isPonded, fltsat, thetaPrev, hPrev, kavg, C, dt, WH, qtop, qbot);

        // First node : 0 (include boundary cond. qtop or pond)
        if (isPonded || fltsat) {
            // h at soil surface prescribed, ponding
            thomc[0] = -dt * kavg[1]/dz[0]/disZ[1];
            thomb[0] = -thomc[0] + C[0] + dt*kavg[0]/disZ[0]/dz[0];
            thomf[0] = C[0]*h[0] + dt/(-dz[0]) * (kavg[0] - kavg[1]) + dt*kavg[0]*WH/disZ[0]/dz[0];
        } else {
            //  q at soil surface prescribed, qtop = rainfall
            isPonded = false;
            thomc[0] = -dt * kavg[1] / (dz[0]*disZ[1]);
            thomb[0] = -thomc[0] + C[0];
            thomf[0] = C[0]*h[0] + dt/(-dz[0]) * (-qtop - kavg[1]) - dt*S[0];
        }

        // Intermediate nodes: i = 1 to n-2
        for (int i = 1; i < nN-1; i++) {
            thoma[i] = -dt*kavg[i]/dz[i]/disZ[i];
            thomc[i] = -dt*kavg[i+1]/dz[i]/disZ[i+1];
            thomb[i] = -thoma[i] - thomc[i] + C[i];
            thomf[i] = C[i]*h[i] + dt/-dz[i]*(kavg[i]-kavg[i+1]) - dt*S[i];  //add sinkterm according to Belmans
            // Belmans: E = h + (dt/C*dz)K+1/2 + (dt/C*dz)K-1/2 - (dt/C)*S;
            // F = C*E = Ch + dt*dz*K+1/2 +dt*dz*K-1/2  -dt*S
            //dh/dt = 1/C* 1/dz etc -dt*S/C eq 6 page 275
            // S is ET flux also in cm/day
        }

        // last node : nN-1 (include boundary cond. qbot)
        thoma[nN-1] = -dt*kavg[nN-1]/dz[nN-1]/disZ[nN-1];
        thomb[nN-1] = -thoma[nN-1] + C[nN-1];
        thomf[nN-1] = C[nN-1]*h[nN-1] + dt/(-dz[nN-1])*(kavg[nN-1]+qbot) - dt*S[nN-1];

        // Gaussian elimination and backsubstitution h - first time
        double alpha = thomb[0];
        h[0] = thomf[0] / alpha;
        for (int i = 1; i < nN; i++) {
            beta[i] = thomc[i-1] / alpha;
            alpha = thomb[i] - thoma[i] * beta[i];
            h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
        }
        for (int i = (nN-2); i >= 0; i--)
            h[i] -= beta[i+1] * h[i+1];

        // correct tridiagonal matrix
        for (int i = 0; i < nN; i++) {
            double thetaNew = FindValue(h[i], p->horizon[i], H_COL, THETA_COL);

            // if (SwitchDensCorrection && p->zone->endComp[i] <= 30 && h[i] > -10.0)
            //     thetaNew = pixel->corrPDA*thetaNew + pixel->corrPDB;
            // if (SwitchOMCorrection && p->zone->endComp[i] <= 30 && h[i] > -10.0)
            //     thetaNew = pixel->corrPOA*thetaNew + pixel->corrPOB;

            double CNew = FindValue(h[i], p->horizon[i], DMCH_COL, DMCC_COL);
            thomb[i] = thomb[i] - C[i] + CNew;
            thomf[i] = thomf[i] - C[i]*hPrev[i] + CNew*h[i]
                    - thetaNew + thetaPrev[i];
        }

        // Gaussian elimination and backsubstitution h - second time
        alpha = thomb[0];
        h[0] = thomf[0] / alpha;
        for (int i = 1; i < nN; i++) {
            beta[i] = thomc[i-1] / alpha;
            alpha = thomb[i] - thoma[i] * beta[i];
            h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
        }

        for (int i = (nN-2); i >= 0; i--)
            h[i] -= beta[i+1] * h[i+1];

        // we don't need this unless for output
        // for (int j = 0; j < nN; j++)
        //     theta[j] = FindValue(h[j], p->horizon[j], H_COL, THETA_COL);
        //

        // determine new boundary fluxes

        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = -kavg[nN-1]*(h[nN-1]-h[nN-2])/disZ[nN-1] - kavg[nN-1];
        //qbot = kavg[n-1]*(h[n-1]-h[n-2])/disZ[n-1] - kavg[n-1];
        percolation += qbot*dt;

        if (isPonded || fltsat)
             qtop = -kavg[0] * ((h[0] - WH)/disZ[0] + 1) * (1.0-impfrac);
        // else qtop is WH/dt !

        WH += qtop*dt;       // decrease pond with top flux
        WH = qMax(WH, 0.0);

        //influx += qmax*dt;
        // add max infil to influx (negative), to get potential infil
        // not used

        //--- calculate tile drain ---//
        //TODO: CHECK THIS
        if (SwitchIncludeTile && tnode > 0) {
            if (h[tnode] >= TileEntrySuction) {
                double vollayer = -disZ[tnode]*0.01 * CHAdjDX->Drc; // m3
                qdrain =  0.01*k[tnode]*dt*TileDiameter->Drc*DX->Drc; // m3
                double water = theta[tnode] * vollayer; // m3
                // total amonut of water available to drain in this node (m3)
                // note: distnode has a negative value (in cm so 0.01)
                qdrain = qMin(qdrain, water);
                // cannot have more drainage than water available
                water -= qdrain;
                theta[tnode] = water/vollayer; //m3/m3
                h[tnode] = FindValue(theta[tnode], p->horizon[tnode], THETA_COL, H_COL );
                hPrev[tnode] = h[tnode];
                // new h from theta

                drainout += qdrain;
                // add for all swatre timestps, in m3
            }
        }

        // estimate new dt within lisemtimestep
        dt = NewTimeStep(dt, hPrev, h, nN, swatreDT, SwatrePrecision);

        if (elapsedTime+dt >= _dt - TIME_EPS)
            dt = _dt - elapsedTime;

        elapsedTime += dt;

    } // elapsedTime < lisemTimeStep

    double sumth = 0;
    double n = 0;

    for (int i = 0; i < nN; i++) {
        if (p->zone->rootz[i] > 0){
            sumth += FindValue(h[i], p->horizon[i], H_COL, THETA_COL);
            n += 1.0;
        }
    }
    pixel->thetaroot = sumth/n;

    //put new h back into h
    //memcpy(pixel->h.data(), h, nN * sizeof(double));
    for (int j = 0; j < nN; j++) {
        pixel->h[j] = h[j];
    }
    // these variables can all be direcvtly saved to the maps, inflated pixel structure
    pixel->wh = WH*0.01; //convert cm to m
    pixel->tiledrain = drainout;
    pixel->percolation = -percolation*0.01; // cm to m, this is not a flux?

    // theta.clear();
    // kavg.clear();
    // k.clear();
    // C.clear();
    // thetaPrev.clear();
    // h.clear();
    // hPrev.clear();
    // dz.clear();
    // disZ.clear();
    // S.clear();
    // thoma.clear();
    // thomb.clear();
    // thomc.clear();
    // thomf.clear();
    // beta.clear();

   delete[] theta;
   delete[] kavg;
   delete[] k;
   delete[] C;
   delete[] thetaPrev;
   delete[] h;
   delete[] hPrev;
   delete[] dz;
   delete[] disZ;
   delete[] S;
   delete[] thoma;
   delete[] thomb;
   delete[] thomc;
   delete[] thomf;
   delete[] beta ;

}
//--------------------------------------------------------------------------------
