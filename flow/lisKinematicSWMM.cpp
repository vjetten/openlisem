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

#include <math.h>
#include "model.h"

//#include "findroot.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXIT 60



// based on https://github.com/USEPA/Stormwater-Management-Model.git
//
//  Input:   j = link index
//           qinflow = inflow at current time (cfs)
//           tStep = time step (sec)
//  Output:  qoutflow = outflow at current time (cfs),
//           returns number of iterations used
//  Purpose: finds outflow over time step tStep given flow entering a
//           conduit using Kinematic Wave flow routing.
//
//
//                               ^ q3
//  t                            |
//  |          qin, ain |-------------------| qout, aout
//  |                   |  Flow --->        |
//  |----> x     q1, a1 |-------------------| q2, a2
//
//

/*
else ain = xsect_getAofS(pXsect, qin/Beta1) / Afull;

xsect_getAofS(TXsect* xsect, double s)
    Ψ is known as the section factor (Chow, 1959) and is a function of the flow area and conduit geometry. F
    psi = A*R^2/3 = area*hydr radius ^2/3 so that Q = beta*psi(A)  pagge 64 eq 4-6
    double psi = s / xsect->sFull;  = S/(A*R^2/3)
    if ( s <= 0.0 ) return 0.0;
    if ( s > xsect->sMax ) s = xsect->sMax;
    circ_getAofS(xsect, s);

        double circ_getAofS(TXsect* xsect, double s)
        {
            double psi = s / xsect->sFull;
            if (psi == 0.0) return 0.0;
            if (psi >= 1.0) return xsect->aFull;

            // --- use special function for small s/sFull
            if (psi <= 0.015) return xsect->aFull * getAcircular(psi);

            // --- otherwise use table
            else return xsect->aFull * invLookup(psi, S_Circ, N_S_Circ);
           // table C-1 page 153 hydraulic manual
        }
//s = qin/beta1 = qin/(beta/Qfull) = qin/(beta/((A/P)^5/3 * beta) = qin/(1/(A/P^5/3)) ???
*/

// Function to compute psi(theta)
double psi_func(double r, double theta) {
    double A = 0.5 * r * r * (theta - sin(theta));
    double P = r * theta;
    qDebug() << "ap" << r << theta << P << A;
    double R = A / P;
    return A * pow(R, 2.0 / 3.0);
    //(theta - sin(theta)) / (2 * PI) - area_ratio;//
}

// Derivative of psi with respect to theta (numerical)
double psi_derivative(double r, double theta, double h = 1e-6) {
    return (psi_func(r, theta + h) - psi_func(r, theta - h)) / (2 * h);
}

// f = (theta - sin(theta)) / (2 * pi) - area_ratio
// df = (1 - cos(theta)) / (2 * pi)
// theta_next = theta - f / df
// if abs(theta_next - theta) < tol:
//     return theta_next
// theta = theta_next

// Newton-Raphson to solve for theta
double solve_theta(double r, double psi_target, double tol = 1e-6, int max_iter = 100) {
    double theta = PI; // initial guess
    for (int i = 0; i < max_iter; ++i) {
        double f = psi_func(r, theta) - psi_target;
        double df = psi_derivative(r, theta);
        qDebug()<< psi_target << f << df << theta;
        double delta = f / df;
        theta -= delta;
        if (fabs(delta) < tol) {
            return theta;
        }
    }
    return theta;
}

double TWorld::getAfromS(DRAIN_PROP *dr, double s)
{
    double psi = s / dr->sFull;
    if (psi == 0.0) return 0.0;
    if (psi >= 1.0) return dr->Afull;

    double r =  TileDiameter->data[dr->r][dr->c]/2.0;
    double theta = solve_theta(r, psi);
    double A = 0.5 * r * r * (theta - sin(theta));
    return A;
}


void TWorld::PipeFlowSWMM()
{

    // walk down the network!!!

    downstream(crlinkedlddtile_, TileA, tma);
    downstream(crlinkedlddtile_, TileQ, tmb);
    Fill(*Qn,0);
    DRAIN_PROP *drain = new DRAIN_PROP;

    for(long i_ =  0; i_ < crlinkedlddtile_.size(); i_++)
    {
        int r = crlinkedlddtile_[i_].r;
        int c = crlinkedlddtile_[i_].c;
        double Qin = 0;
        // get inflow
        if (crlinkedlddtile_[i_].nr > 0) {
            for(int j = 0; j < crlinkedlddtile_[i_].nr; j++) {
                int rr = crlinkedlddtile_[i_].inn[j].r;
                int cr = crlinkedlddtile_[i_].inn[j].c;
                Qin += TileQn->Drcr;
            }
        }
        TileQin->Drc = Qin;

        if (Qin < 1e-12 && TileQ->Drc < 1e-12) {
            TileQn->Drc = 0;
            TileWaterVol->Drc = 0;
            return;
        }
/* see page 82 hydraulic
 *
 * y = flow depth, yFull is perimeter
 *     case CIRCULAR:
        xsect->yFull = p[0]/ucf;
        xsect->wMax  = xsect->yFull; // diameter
        xsect->aFull = PI / 4.0 * xsect->yFull * xsect->yFull;  (pi r^2 = pi*d/2*d/2)
        xsect->rFull = 0.2500 * xsect->yFull;
        xsect->sFull = xsect->aFull * pow(xsect->rFull, 2./3.);
        xsect->sMax  = 1.08 * xsect->sFull;
        xsect->ywMax = 0.5 * xsect->yFull;
        break;
    Link[j].qFull = Link[j].xsect.sFull * Conduit[k].beta;
    Conduit[k].qMax = Link[j].xsect.sMax * Conduit[k].beta;

import math

def solve_theta(area_ratio, tol=1e-6, max_iter=100):
    theta = pi  # Good starting guess for half full
    for _ in range(max_iter):
        f = (theta - sin(theta)) / (2 * pi) - area_ratio
        df = (1 - cos(theta)) / (2 * pi)
        theta_next = theta - f / df
        if abs(theta_next - theta) < tol:
            return theta_next
        theta = theta_next
    raise RuntimeError("Did not converge")

   double        yFull;           // depth when full (ft)
   double        wMax;            // width at widest point (ft)
   double        ywMax;           // depth at widest point (ft)
   double        aFull;           // area when full (ft2)
   double        rFull;           // hyd. radius when full (ft)
   double        sFull;           // section factor when full (ft^4/3)
   double        sMax;            // section factor at max. flow (ft^4/3)

*/
        drain->c = c;
        drain->r = r;
        drain->beta = sqrt(TileGrad->Drc)/TileN->Drc;  // s = qin/beta
        drain->Afull = TileArea->Drc;
        drain->sFull = drain->Afull * std::pow(0.25*TileDiameter->Drc,2.0/3.0);  // 0.5r=0.25D is hydrasulic radius when full
        drain->Qfull = drain->sFull * drain->beta;
        drain->Beta1 = drain->beta / drain->Qfull; // = 1/sFull =>qin/beta1 = qin/(beta/Qfull)
        drain->dxdt = _dx/_dt * drain->Afull / drain->Qfull;
        drain->sMax = 1.08 * drain->sFull;  // circular

        // s = (qin/beta1)
        //double psi = s / xsect->sFull;  = s/(A*R^2/3)

        int    result = 1;
        double dq;
        double WT = 0.6;
        double WX = 0.6;

        // --- normalize previous flows, averrage with downstream for now
        drain->q1 = TileQ->Drc / drain->Qfull;
        drain->q2 = ((TileQ->Drc + tmb->Drc)*0.5)/ drain->Qfull;
        // --- normalize inflow
        drain->qin = std::min(drain->Qfull, TileQin->Drc)/drain->Qfull;
        // in SWMM code the inflow is maximized to the possible incflow

        // --- compute evaporation and infiltration loss rate
       // double q3 = 0;//link_getLossRate(j, KW, qin*Qfull, tStep) / Qfull;

        // --- normalize previous areas, averrage with downstream
        drain->a1 = TileA->Drc / drain->Afull;
        drain->a2 = ((TileA->Drc + tma->Drc)*0.5)/ drain->Afull;

        // --- use full area when inlet flow >= full flow
        if ( drain->qin >= 1.0 )
            drain->ain = 1.0;
        else
            drain->ain = getAfromS(drain, drain->qin/drain->Beta1)/drain->Afull;
            // --- get normalized inlet area corresponding to inlet flow
         //   drain->ain = (drain->qin/drain->Beta1) / drain->Afull;
qDebug() << "ain" << drain->ain << drain->qin << drain->qin/drain->Beta1;
        // beta1 depends on shape

        // --- check for no flow
        if ( drain->qin < 1e-12 && drain->q2 < 1e-12 ) {
            drain->qout = 0.0;
            drain->aout = 0.0;
        } else {

            dq   = drain->q2 - drain->q1;
            drain->C1   = drain->dxdt * WT / WX;
            drain->C2   = (1.0 - WT) * (drain->ain - drain->a1);
            drain->C2   = drain->C2 - WT * drain->a2;
            drain->C2   = drain->C2 * drain->dxdt / WX;
            drain->C2   = drain->C2 + (1.0 - WX) / WX * dq - drain->qin;
            //drain->C2   = C2 + q3 / WX;

            // --- starting guess for aout is value from previous time step
            drain->aout = drain->a2;

            // --- solve continuity equation for aout
            result = solveContinuity(drain);

            // --- report error if continuity eqn. not solved
            if ( result == -1 )
            {
                //report_writeErrorMsg(ERR_KINWAVE, Link[j].ID);
                Error("Kinwave SWMM error solvecontinuity");
                //return;
            }
            if ( result <= 0 )
                result = 1;

            // --- compute normalized outlet flow from outlet area
            drain->qout = drain->Beta1 * drain->aout*drain->Afull;
            // if ( drain->qin > 1.0 )
            //     drain->qin = 1.0;
            //for the next in line!!!
        }
        TileQn->Drc = drain->qout;
        TileWaterVol->Drc = drain->aout * DX->Drc;
    }
    delete drain;
}

int TWorld::solveContinuity(DRAIN_PROP *dr)
//
//  Input:   qin = upstream normalized flow
//           ain = upstream normalized area
//           aout = downstream normalized area
//  Output:  new value for aout; returns an error code
//  Purpose: solves continuity equation f(a) = Beta1*S(a) + C1*a + C2 = 0
//           for 'a' using the Newton-Raphson root finder function.
//           Return code has the following meanings:
//           >= 0 number of function evaluations used
//           -1   Newton function failed
//           -2   flow always above max. flow
//           -3   flow always below zero
//
//     Note: pXsect (pointer to conduit's cross-section), and constants Beta1,
//           C1, and C2 are module-level shared variables assigned values
//           in kinwave_execute().
//
{
    int    n;                          // # evaluations or error code
    double aLo, aHi, aTmp;             // lower/upper bounds on a
    double fLo, fHi;                   // lower/upper bounds on f

    // --- first determine bounds on 'a' so that f(a) passes through 0.

    // --- set upper bound to area at full flow
    aHi = 1.0;
    fHi = 1.0 + dr->C1 + dr->C2;

    // --- try setting lower bound to area where section factor is maximum
    aLo = std::min(dr->a1,dr->a2)/dr->Afull;//1;//xsect_getAmax(pXsect) / Afull;
    if ( aLo < aHi )
    {
        fLo = ( dr->Beta1 * dr->sMax ) + (dr->C1 * aLo) + dr->C2;
    }
    else fLo = fHi;

    // --- if fLo and fHi have same sign then set lower bound to 0
    if ( fHi*fLo > 0.0 )
    {
        aHi = aLo;
        fHi = fLo;
        aLo = 0.0;
        fLo = dr->C2;
    }

    // --- proceed with search for root if fLo and fHi have different signs
    if ( fHi*fLo <= 0.0 )
    {
        // --- start search at midpoint of lower/upper bounds
        //     if initial value outside of these bounds
        if ( dr->aout < aLo || dr->aout > aHi ) dr->aout = 0.5*(aLo + aHi);

        // --- if fLo > fHi then switch aLo and aHi
        if ( fLo > fHi )
        {
            aTmp = aLo;
            aLo  = aHi;
            aHi  = aTmp;
        }

        // --- call the Newton root finder method passing it the
        //     evalContinuity function to evaluate the function
        //     and its derivatives
        n = findroot_Newton(dr, aLo, aHi); //,NULL);

        // --- check if root finder succeeded
        if ( n <= 0 ) n = -1;
    }

    // --- if lower/upper bound functions both negative then use full flow
    else if ( fLo < 0.0 )
    {
        if ( dr->qin > 1.0 ) dr->aout = dr->ain;
        else dr->aout = 1.0;
        n = -2;
    }

    // --- if lower/upper bound functions both positive then use no flow
    else if ( fLo > 0 )
    {
        dr->aout = 0.0;
        n = -3;
    }
    else n = -1;
    return n;
}


int TWorld::findroot_Newton(DRAIN_PROP *dr, double x1, double x2)
//
//  Using a combination of Newton-Raphson and bisection, find the root of a
//  function func bracketed between x1 and x2. The root, returned in rts,
//  will be refined until its accuracy is known within +/-xacc. func is a
//  user-supplied routine, that returns both the function value and the first
//  derivative of the function. p is a pointer to any auxilary data structure
//  that func may require. It can be NULL if not needed. The function returns
//  the number of function evaluations used or 0 if the maximum allowed
//  iterations were exceeded.
//
// NOTES:
// 1. The calling program must insure that the signs of func(x1) and func(x2)
//    are not the same, otherwise x1 and x2 do not bracket the root.
// 2. If func(x1) > func(x2) then the order of x1 and x2 should be
//    switched in the call to Newton.
//
{
    int j, n = 0;
    double df, dx, dxold, f, x;
    double temp, xhi, xlo;

    // Initialize the "stepsize before last" and the last step.
    x = dr->aout;
    xlo = x1;
    xhi = x2;
    dxold = fabs(x2-x1);
    dx = dxold;

    n++;

    // Loop over allowed iterations.
    for (j=1; j<=MAXIT; j++)
    {
        // Bisect if Newton out of range or not decreasing fast enough.
        if ( ( ( (x-xhi)*df-f)*((x-xlo)*df-f) >= 0.0 || (fabs(2.0*f) > fabs(dxold*df) ) ) ) {
            dxold = dx;
            dx = 0.5*(xhi-xlo);
            x = xlo + dx;
            if ( xlo == x ) break;
        } else {
            // Newton step acceptable. Take it.
            dxold = dx;
            dx = f/df;
            temp = x;
            x -= dx;
            if ( temp == x )
                break;
        }

        // Convergence criterion.
        if ( fabs(dx) < EPSILON )
            break;

        f = dr->Beta1 * (xlo*dr->Afull) + dr->C1*xlo + dr->C2;
        df = dr->Beta1*dr->Afull* (xlo*dr->Afull) + dr->C1;
        //xlo = a
        // *f  = (Beta1 * xsect_getSofA(pXsect, a*Afull)) + (C1 * a) + C2;
        //  Input:   xsect = ptr. to a cross section data structure
        //           a = area (ft2)
        //  Output:  returns section factor (ft^(8/3))

        // *df = (Beta1 * Afull * xsect_getdSdA(pXsect, a*Afull)) + C1;
        //  Purpose: computes xsection's section factor at a given area.
        //  Input:   xsect = ptr. to a cross section data structure
        //           a = area (ft2)
        //  Output:  returns derivative of section factor w.r.t. area (ft^2/3)
        //  Purpose: computes xsection's derivative of its section factor with
        //           respect to area at a given area.

        // this is like a newton raphson on height instead of discharge
        n++;
        if ( f < 0.0 )
            xlo = x;
        else
            xhi = x;
    }
    dr->aout = x;
    if (n <= MAXIT)
        return n;
    else
        return 0;
}
