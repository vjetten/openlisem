/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
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

#include <algorithm>
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
thomc has unit cm
thomb has unit dtheta/dh (like dimoca)

*/
//OBSOLETE
void TWorld::HeadCalc(double *h, const PROFILE *p , bool *isPonded, bool fltsat,
              const double *thetaPrev, const double *hPrev, const double *kavg, const double *dimoca,
              double dt, double pond, double qtop, double qbot)
{
    int nN = NrNodes(p);
    //const double *dz = Dz(p), *disnod = DistNode(p);
    // dz is layer thickness, distnode is distance between centre of layers
    int i;
    NODE_ARRAY thoma, thomb, thomc, thomf, beta;
    double alpha;

    // First node : 0 (include boundary cond. qtop or pond)
    if ( *isPonded || fltsat ) {
        // h at soil surface prescribed, ponding
        thomc[0] = -dt * kavg[1] / Dz(p)[0] / DistNode(p)[1];
        thomb[0] = -thomc[0] + dimoca[0] + dt*kavg[0]/DistNode(p)[0]/Dz(p)[0];
        thomf[0] = dimoca[0]*h[0] + dt/(-Dz(p)[0]) * (kavg[0] - kavg[1]) +
                dt*kavg[0]*pond/DistNode(p)[0]/Dz(p)[0];
    } else {
        //  q at soil surface prescribed, qtop = rainfall
        *isPonded = false;
        thomc[0] = -dt * kavg[1] / Dz(p)[0] / DistNode(p)[1];
        thomb[0] = -thomc[0] + dimoca[0];
        thomf[0] = dimoca[0]*h[0] + dt/(-Dz(p)[0]) * (-qtop - kavg[1]); //(- qtop - kavg[1]);
    }

    // Intermediate nodes: i = 1 to n-2
    for (i = 1; i < nN-1; i++)
    {
        thoma[i] = -dt*kavg[i]/Dz(p)[i]/DistNode(p)[i];
        thomc[i] = -dt*kavg[i+1]/Dz(p)[i]/DistNode(p)[i+1];
        thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
        thomf[i] = dimoca[i]*h[i] + dt/(-Dz(p)[i])*(kavg[i]-kavg[i+1]);
    }

    // last node : nN-1 (include boundary cond. qbot)
    thoma[nN-1] = -dt*kavg[nN-1]/Dz(p)[nN-1]/DistNode(p)[nN-1];
    thomb[nN-1] = -thoma[nN-1] + dimoca[nN-1];
    thomf[nN-1] = dimoca[nN-1]*h[nN-1] + dt/(-Dz(p)[nN-1])*(kavg[nN-1]+qbot);

    // Gaussian elimination and backsubstitution h - first time
    alpha = thomb[0];
    h[0] = thomf[0] / alpha;
    for (i = 1; i < nN; i++) {
        beta[i] = thomc[i-1] / alpha;
        alpha = thomb[i] - thoma[i] * beta[i];
        h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
    }
    for (i = (nN-2); i >= 0; i--)
        h[i] -= beta[i+1] * h[i+1];

    // correct tridiagonal matrix
    for (i = 0; i < nN; i++) {
        double theta = TheNode(h[i], Horizon(p,i));
        double dimocaNew = DmcNode(h[i], Horizon(p,i));
        thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
        thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i]
                - theta + thetaPrev[i];
    }

    // Gaussian elimination and backsubstitution h - second time
    alpha = thomb[0];
    h[0] = thomf[0] / alpha;
    for (i = 1; i < nN; i++) {
        beta[i] = thomc[i-1] / alpha;
        alpha = thomb[i] - thoma[i] * beta[i];
        h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
    }

    for (i = (nN-2); i >= 0; i--)
        h[i] -= beta[i+1] * h[i+1];

}
//--------------------------------------------------------------------------------
/**
 * @brief
 *
 * @param prevDt
 * @param hLast
 * @param h
 * @param nrNodes
 * @param precParam
 * @param dtMin
 * @param dtMax
 * @return double
 */
double  TWorld::NewTimeStep(
        double 		 prevDt,
        const double *hLast,
        const double *h,
        int    		 nrNodes,
        double 		 precParam,
        double 		 dtMin,
        double 		 dtMax)
{
    int i;
    double dt = _dt;//dtMax;
    double accur1 = 0.3 - 0.02 * precParam;
    double accur2 = 0.03 - 0.002 * precParam;

    for(i=0; i < nrNodes; i++)
    {
        double mdih = accur1 + accur2 * std::max(1.0, fabs(h[i]));
        double dih  = fabs(h[i] - hLast[i]);
        // if difference is small
        // dih = e.g. 10 and h = -200 then mdih = 200*0.01 + 0.1 = 2.1
        // mdih/dih = 2.1/10 =0.21

        if (dih > 0.10)
            dt = std::min(dt, prevDt*mdih/dih);
    }
    return (std::max(dt, dtMin));
}
//--------------------------------------------------------------------------------
// Units are:
// Z and H in cm; table units K in cm/day converted to cm/sec, lisem time in seconds

void TWorld::ComputeForPixel(PIXEL_INFO *pixel, SOIL_MODEL *s, double drainfraction)
{
    NODE_ARRAY theta, thetaPrev, hPrev, dimoca, kavg, k;
    const PROFILE *p = pixel->profile;

    int nN = NrNodes(p);
    double dt = _dt/5;
    double pond = pixel->wh;
    double elapsedTime = 0;
    double influx = 0;
    double drainout = 0;
    double percolation = 0;
    double Theta = 0;
    int tnode = pixel->tilenode;

    if (SHOWDEBUG)
        qDebug() << "compute for pixel" << nN << pond << p->profileId;

    //double *h = new double [n];
    NODE_ARRAY h;
    for (int i = 0; i < nN; i++) {
        h[i] = pixel->h[i];
    }

    while (elapsedTime < _dt)
    {

        bool isPonded, fltsat;    // flag if ponded or if profile fully saturated
        double qmax, qtop, qbot, ThetaSat;  // fluxes at top and bottom, max theta
        double qdrain; // tile drainage

        // get nodal values of theta, K, dif moist cap
        for (int j = 0; j < nN; j++) {
            k[j] = HcoNode(h[j], p->horizon[j], ksatCalibration);
            // k[j] = FindNode(h[j], p->horizon[j], K_COL);
            // input tables in cm/day function returns in cm/sec !!
            dimoca[j] = DmcNode(h[j], p->horizon[j]);
            // dimoca[j] = FindNode(h[j], p->horizon[j], DMCC_COL);
            // differential moisture capacity d(theta)/d(h), tangent moisture retention curve
            theta[j] = TheNode(h[j], p->horizon[j]);
            // moisture content
            // if(SHOWDEBUG)
            // qDebug() << h[j] << k[j] << theta[j] << dimoca[j];
        }

        Theta = (theta[0]+theta[1])/2;

        // average K for 1st to n-1 node, top node is done below
        for(int j = 1; j < nN; j++) {
            kavg[j] = (k[j-1]+k[j])/2.0;
            //kavg[j] = std::sqrt(k[j-1]*k[j]);
        }

        //--- boundary conditions ---

        //----- TOP -----
        // check if ponded: 1st compare fluxes, 2nd compare store
        // top flux is water/timestep (cm/sec)
        qtop = -pond/dt;

        //----- BOTTOM -----
        // bottom is 0 or copy of flux of last 2 layers
        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = kavg[nN-1]*(h[nN-1]-h[nN-2])/DistNode(p)[nN-1] - kavg[nN-1];

        // 1st check flux against max flux
        double Ksat = HcoNode(0, p->horizon[0], ksatCalibration);
        //FindNode(0, p->horizon[0], K_COL)*ksatCalibration;
        //kavg[0]= sqrt( Ksat * k[0]);
        kavg[0]= ( Ksat + k[0])/2.0;
        // geometric avg of ksat and k[0] => is used for max possible

        qmax = -kavg[0]*(pond-h[0]) / DistNode(p)[0] - kavg[0];
        qmax = std::min(qmax, 0.0);
        // maximum possible flux, compare to real top flux available
        isPonded = (qtop < qmax);
        // if more flux then max possible flag ponded is true

        //2nd check: isPonded layer depth against storage
        if (!isPonded) {
            // calculate available space in profile in cm: (pore-theta)*dz
            double space = 0;
            for(int i = 0; i < nN && h[i] < 0 && space > pond; i++) {
                ThetaSat = TheNode(0, p->horizon[i]);        //FindNode(0, Horizon(p, i), THETA_COL);
                space += (ThetaSat - theta[i]) * (-Dz(p)[i]);
            }
            isPonded = pond > space;
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

        // save last h and theta, used in headcalc
        // for (int i = 0; i < n; i++)  {
        //    hPrev[i] = h[i];
        //    thetaPrev[i] = theta[i];
        // }

        std::memcpy(hPrev, h, nN * sizeof(double));
        std::memcpy(thetaPrev, theta, nN * sizeof(double));

        //HeadCalc(h, p , &isPonded, fltsat,thetaPrev, hPrev, kavg, dimoca,dt, pond, qtop, qbot);

        NODE_ARRAY thoma, thomb, thomc, thomf, beta;
        double alpha;

        /* First node : 0 (include boundary cond. qtop or pond) */
        if ( isPonded || fltsat) {
            // h at soil surface prescribed, ponding
            thomc[0] = -dt * kavg[1] / Dz(p)[0] / DistNode(p)[1];
            thomb[0] = -thomc[0] + dimoca[0] + dt*kavg[0]/DistNode(p)[0]/Dz(p)[0];
            thomf[0] = dimoca[0]*h[0] + dt/(-Dz(p)[0]) * (kavg[0] - kavg[1]) +
                    dt*kavg[0]*pond/DistNode(p)[0]/Dz(p)[0];
        } else {
            //  q at soil surface prescribed, qtop = rainfall
            isPonded = false;
            thomc[0] = -dt * kavg[1] / Dz(p)[0] / DistNode(p)[1];
            thomb[0] = -thomc[0] + dimoca[0];
            thomf[0] = dimoca[0]*h[0] + dt/(-Dz(p)[0]) * (-qtop - kavg[1]); //(- qtop - kavg[1]);
        }


        /* Intermediate nodes: i = 1 to n-2 */
        for (int i = 1; i < nN-1; i++)
        {
            thoma[i] = -dt*kavg[i]/Dz(p)[i]/DistNode(p)[i];
            thomc[i] = -dt*kavg[i+1]/Dz(p)[i]/DistNode(p)[i+1];
            thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
            thomf[i] = dimoca[i]*h[i] + dt/(-Dz(p)[i])*(kavg[i]-kavg[i+1]);
        }

        // last node : nN-1 (include boundary cond. qbot)
        thoma[nN-1] = -dt*kavg[nN-1]/Dz(p)[nN-1]/DistNode(p)[nN-1];
        thomb[nN-1] = -thoma[nN-1] + dimoca[nN-1];
        thomf[nN-1] = dimoca[nN-1]*h[nN-1] + dt/(-Dz(p)[nN-1])*(kavg[nN-1]+qbot);

        // Gaussian elimination and backsubstitution h - first time
        alpha = thomb[0];
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
            double theta = TheNode(h[i], Horizon(p,i));
            double dimocaNew = DmcNode(h[i], Horizon(p,i));
            thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
            thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i]
                    - theta + thetaPrev[i];
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


        for (int j = 0; j < nN; j++) {
            theta[j] = TheNode(h[j], p->horizon[j]);
        }

        // determine new boundary fluxes

        if (SwitchImpermeable)
            qbot = 0;
        else
            qbot = -kavg[nN-1]*(h[nN-1]-h[nN-2])/DistNode(p)[nN-1] - kavg[nN-1];
        //qbot = kavg[n-1]*(h[n-1]-h[n-2])/DistNode(p)[n-1] - kavg[n-1];
        // new qbot is actually not use but may come in handy later
        percolation += qbot*dt;

        if ( isPonded || fltsat)
            // qtop = -kavg[0] * ((h[0] - pond)/DistNode(p)[0] + 1);
            qtop = -kavg[0]*(pond - h[0]) / DistNode(p)[0] - kavg[0];
        qtop = std::min(0.0,qtop);

        // adjust top flux
        if (SHOWDEBUG) {
            QList <double> s;
            for (int j=0; j < nN; j++)
                s << h[j];

            qDebug()<< dt << qtop << "h" << s;
            //qDebug()<< qtop << "th" << theta[0] << theta[2] << theta[3] << theta[4] << theta[5] << theta[6] << theta[7] << theta[8] << theta[9];
        }

        pond += qtop*dt;       // decrease pond with top flux
        if (pond < POND_EPS)
            pond = 0;
        influx += qmax*dt;
        // add max infil to influx (negative), to get potential infil  WHY???

        //--- calculate tile drain ---//
        if (SwitchIncludeTile && tnode > 0) //VJ 110825 tnode = -1 if cell has no drainage
        {
            //options:
            qdrain =  k[tnode];
            //qdrain =  HcoNode(0, Horizon(p, tnode), 1.0);
            // drainage in node is ksat of node, no calibration!
            // drainage is cond of the node in cm/sec
            double water = theta[tnode] * -DistNode(p)[tnode] * drainfraction;
            // total amonut of water available to drain in this node (cm)
            // note: distnode has a negative value
            qdrain = std::min(qdrain, water/dt);
            // cannot have more drainage than water available

            theta[tnode] = std::max(0.001, theta[tnode] - (qdrain*dt)/DistNode(p)[tnode]* drainfraction);
            // adjust theta with drainage removed

            h[tnode] = HNode(theta[tnode], Horizon(p, tnode) );
            hPrev[tnode] = h[tnode];
            // new h from theta

            drainout += qdrain*dt;
            // add for all swatre timestps, in cm
        }

        elapsedTime += dt;
        /* estimate new dt within lisemtimestep */
        dt = NewTimeStep(dt, hPrev, h, nN, precision, s->minDt, _dt);

        if (elapsedTime+dt+TIME_EPS >= _dt)
            dt = _dt - elapsedTime;
        // elapsedTime = _dt+TIME_EPS;

    } // elapsedTime < lisemTimeStep

    /*
    if (pixel->dumpH>0)
       printf("Cell %4d : wh after infil. %8.5f (mm) infil. %8.5f (mm)\n"\
       ,pixel->dumpH,pond*10,-influx*10);
   //   if (pixel->dumpHid == 1)
   //      qDebug() << pond << influx << h[0] << h[1] << h[2] << h[3] << h[4] << h[5] << h[6] << h[7] << h[8] << h[9];
   */

    pixel->wh = pond;
    pixel->tiledrain = drainout;
    pixel->infil = influx;
    pixel->percolation = -percolation;
    pixel->theta = Theta;

    pixel->h.clear();
    for (int i = 0; i < nN; i++) {
        pixel->h.append(h[i]);
    }

}
//--------------------------------------------------------------------------------
// units in SWATRE are cm and cm/day
void TWorld::SwatreStep(long i_, int r, int c, SOIL_MODEL *s, cTMap *_WH, cTMap *_drain, cTMap *_theta)
{   
    double drainfraction = 0;
    QString dig;

    showr = r;
    showc = c;

    s->pixel[i_].wh = _WH->Drc*100;    // WH is in m, convert to cm
    s->pixel[i_].infil = 0;
    s->pixel[i_].tiledrain = 0;

    if (SwitchIncludeTile)
        drainfraction = TileWidth->Drc/_dx;

    ComputeForPixel(&s->pixel[i_], s, drainfraction);
    // estimate new h and theta at the end of dt

    //SwitchDumpH = true;
    if(SwitchDumpH || SwitchDumpTheta || SwitchDumpK) {
        if(s->pixel[i_].dumpHid > 0) {
            for (int i = 0; i < zone->nrNodes; i++) {
                QString name = QString("SwH%1").arg(runstep,2, 10, QLatin1Char('0'));
                dig = QString("%1").arg(i, 12-name.length(), 10, QLatin1Char('0'));
                name=name+dig;
                name.insert(8, ".");
                //qDebug() << name << dig;
            }
        }
    }

    _WH->Drc = s->pixel[i_].wh*0.01;
    _theta->Drc = s->pixel[i_].theta;

    // _fpot->Drc = std::max(0.0, std::abs(infil)*0.01);
    // infil is negative (downward flux * dt, in cm)
    //fpot is positive like in other infil  methods (in m)

    if (SwitchIncludeTile)
        _drain->Drc = s->pixel[i_].tiledrain*0.01;  // in m
    // drained water from the soil, already accounts for drainwidth versus i_l width
}
//--------------------------------------------------------------------------------



