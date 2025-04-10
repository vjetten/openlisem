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
void TWorld::PipeFlowSWMM()
{

    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {
        upstream(crlinkedlddtile_, TileQ, TileQin);
        upstream(crlinkedlddtile_, TileA, TileAin);
        double Afull = TileArea->Drc;
        double Qfull = std::pow(TileArea->Drc/TileDiameter,5.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;        
        double dxdt = _dx/_dt * Afull / Qfull;
        double Beta1 = /*Conduit[k].beta */0.6 / Qfull;
        int    k;
        int    result = 1;
        double dq;
        double ain, aout;
        double qin, qout;
        double a1, a2, q1, q2, q3;
        double WX = 0.6;
        double WT = 0.6;

        // --- normalize previous flows
        q1 = std::min(TileQin->Drc/Qfull) / Qfull;
        q2 = TileQ->Drc / Qfull;

        // --- normalize inflow
        //qin = (*qinflow) / Qfull;
        qin = TileQin->Drc/Qfull;

        // --- compute evaporation and infiltration loss rate
        q3 = 0;//link_getLossRate(j, KW, qin*Qfull, tStep) / Qfull;

        // --- normalize previous areas
        a1 = std::min(Afull, TileAin->Drc) / Afull;
        a2 = (TileA->Drc) / Afull;

        // --- use full area when inlet flow >= full flow
        if ( qin >= 1.0 ) ain = 1.0;

        // --- get normalized inlet area corresponding to inlet flow
        else
            ain = (qin/Beta1) / Afull;

        // --- check for no flow
        if ( qin <= 1e-12 && q2 <= 1e-12 )
        {
            qout = 0.0;
            aout = 0.0;
        }

        dq   = q2 - q1;
        C1   = dxdt * WT / WX;  // WT = 0.6; WX = 0.6
        C2   = (1.0 - WT) * (ain - a1);
        C2   = C2 - WT * a2;
        C2   = C2 * dxdt / WX;
        C2   = C2 + (1.0 - WX) / WX * dq - qin;
        C2   = C2 + q3 / WX;

        // --- starting guess for aout is value from previous time step
        aout = a2;

        // --- solve continuity equation for aout
        result = solveContinuity(r, c, C1, C2, Beta1, TileQin, TileAin, &aout);
    }}

}

int TWorld::solveContinuity(int r, int c, double C1, double C2, double Beta1, double qin, double ain, double* aout)
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
    double tol = 1e-10;//EPSIL;                // absolute convergence tol.

    // --- first determine bounds on 'a' so that f(a) passes through 0.

    // --- set upper bound to area at full flow
    aHi = 1.0;
    fHi = 1.0 + C1 + C2;

    // --- try setting lower bound to area where section factor is maximum
    aLo = TileArea->Drc / Afull;
    if ( aLo < aHi )
    {
        fLo = ( Beta1 * pXsect->sMax ) + (C1 * aLo) + C2; // Smax= m^4/3
    }
    else fLo = fHi;

    // --- if fLo and fHi have same sign then set lower bound to 0
    if ( fHi*fLo > 0.0 )
    {
        aHi = aLo;
        fHi = fLo;
        aLo = 0.0;
        fLo = C2;
    }

    // --- proceed with search for root if fLo and fHi have different signs
    if ( fHi*fLo <= 0.0 )
    {
        // --- start search at midpoint of lower/upper bounds
        //     if initial value outside of these bounds
        if ( *aout < aLo || *aout > aHi ) *aout = 0.5*(aLo + aHi);

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
        n = findroot_Newton(aLo, aHi, aout, tol, evalContinuity, nullptr);

        // --- check if root finder succeeded
        if ( n <= 0 ) n = -1;
    }

    // --- if lower/upper bound functions both negative then use full flow
    else if ( fLo < 0.0 )
    {
        if ( qin > 1.0 ) *aout = ain;
        else *aout = 1.0;
        n = -2;
    }

    // --- if lower/upper bound functions both positive then use no flow
    else if ( fLo > 0 )
    {
        *aout = 0.0;
        n = -3;
    }
    else n = -1;
    return n;
}

void TWorld::evalContinuity(double a, double* f, double* df, void* p)
//
//  Input:   a = outlet normalized area
//  Output:  f = value of continuity eqn.
//           df = derivative of continuity eqn.
//  Purpose: computes value of continuity equation (f) and its derivative (df)
//           w.r.t. normalized area for link with normalized outlet area 'a'.
//
{
    *f  = (Beta1 * xsect_getSofA(pXsect, a*Afull)) + (C1 * a) + C2;
    *df = (Beta1 * Afull * xsect_getdSdA(pXsect, a*Afull)) + C1;
}


int TWorld::findroot_Newton(double x1, double x2, double* rts, double xacc,
                    void (*func) (double x, double* f, double* df, void* p),void* p)
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
    x = *rts;
    xlo = x1;
    xhi = x2;
    dxold = fabs(x2-x1);
    dx = dxold;
    func(x, &f, &df, p);
    n++;

    // Loop over allowed iterations.
    for (j=1; j<=MAXIT; j++)
    {
        // Bisect if Newton out of range or not decreasing fast enough.
        if ( ( ( (x-xhi)*df-f)*((x-xlo)*df-f) >= 0.0
        || (fabs(2.0*f) > fabs(dxold*df) ) ) )
        {
            dxold = dx;
            dx = 0.5*(xhi-xlo);
            x = xlo + dx;
            if ( xlo == x ) break;
        }

        // Newton step acceptable. Take it.
        else
        {
            dxold = dx;
            dx = f/df;
            temp = x;
            x -= dx;
            if ( temp == x ) break;
        }

        // Convergence criterion.
        if ( fabs(dx) < xacc ) break;

        // Evaluate function. Maintain bracket on the root.
        func(x, &f, &df, p);
        n++;
        if ( f < 0.0 ) xlo = x;
        else           xhi = x;
    }
    *rts = x;
    if ( n <= MAXIT) return n;
    else return 0;
};


/*

int kinwave_execute(int j, double* qinflow, double* qoutflow, double tStep)
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
{
    int    k;
    int    result = 1;
    double dxdt, dq;
    double ain, aout;
    double qin, qout;
    double a1, a2, q1, q2, q3;

    // --- no routing for non-conduit link
    (*qoutflow) = (*qinflow);
    if ( Link[j].type != CONDUIT ) return result;

    // --- no routing for dummy xsection
    if ( Link[j].xsect.type == DUMMY ) return result;

    // --- assign module-level variables
    pXsect = &Link[j].xsect;
    Qfull = Link[j].qFull;
    Afull = Link[j].xsect.aFull;
    k = Link[j].subIndex;
    Beta1 = Conduit[k].beta / Qfull;

    // --- normalize previous flows
    q1 = Conduit[k].q1 / Qfull;
    q2 = Conduit[k].q2 / Qfull;

    // --- normalize inflow
    qin = (*qinflow) / Conduit[k].barrels / Qfull;

    // --- compute evaporation and infiltration loss rate
    q3 = link_getLossRate(j, KW, qin*Qfull, tStep) / Qfull;

    // --- normalize previous areas
    a1 = Conduit[k].a1 / Afull;
    a2 = Conduit[k].a2 / Afull;

    // --- use full area when inlet flow >= full flow
    if ( qin >= 1.0 ) ain = 1.0;

    // --- get normalized inlet area corresponding to inlet flow
    else ain = xsect_getAofS(pXsect, qin/Beta1) / Afull;

    // --- check for no flow
    if ( qin <= TINY && q2 <= TINY )
    {
        qout = 0.0;
        aout = 0.0;
    }

    // --- otherwise solve finite difference form of continuity eqn.
    else
    {
        // --- compute constant factors
        dxdt = link_getLength(j) / tStep * Afull / Qfull;
        dq   = q2 - q1;
        C1   = dxdt * WT / WX;  // WT = 0.6; WX = 0.6
        C2   = (1.0 - WT) * (ain - a1);
        C2   = C2 - WT * a2;
        C2   = C2 * dxdt / WX;
        C2   = C2 + (1.0 - WX) / WX * dq - qin;
        C2   = C2 + q3 / WX;

        // --- starting guess for aout is value from previous time step
        aout = a2;

        // --- solve continuity equation for aout
        result = solveContinuity(qin, ain, &aout);

        // --- report error if continuity eqn. not solved
        if ( result == -1 )
        {
            report_writeErrorMsg(ERR_KINWAVE, Link[j].ID);
            return 1;
        }
        if ( result <= 0 ) result = 1;

        // --- compute normalized outlet flow from outlet area
        qout = Beta1 * xsect_getSofA(pXsect, aout*Afull);
        if ( qin > 1.0 ) qin = 1.0;
    }

    // --- save new flows and areas
    Conduit[k].q1 = qin * Qfull;
    Conduit[k].a1 = ain * Afull;
    Conduit[k].q2 = qout * Qfull;
    Conduit[k].a2 = aout * Afull;
    Conduit[k].fullState =
        link_getFullState(Conduit[k].a1, Conduit[k].a2, Afull);
    (*qinflow)  = Conduit[k].q1 * Conduit[k].barrels;
    (*qoutflow) = Conduit[k].q2 * Conduit[k].barrels;
    return result;
}

//=============================================================================

int solveContinuity(double qin, double ain, double* aout)
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
    double tol = EPSIL;                // absolute convergence tol.

    // --- first determine bounds on 'a' so that f(a) passes through 0.

    // --- set upper bound to area at full flow
    aHi = 1.0;
    fHi = 1.0 + C1 + C2;

    // --- try setting lower bound to area where section factor is maximum
    aLo = xsect_getAmax(pXsect) / Afull;
    if ( aLo < aHi )
    {
        fLo = ( Beta1 * pXsect->sMax ) + (C1 * aLo) + C2;
    }
    else fLo = fHi;

    // --- if fLo and fHi have same sign then set lower bound to 0
    if ( fHi*fLo > 0.0 )
    {
        aHi = aLo;
        fHi = fLo;
        aLo = 0.0;
        fLo = C2;
    }

    // --- proceed with search for root if fLo and fHi have different signs
    if ( fHi*fLo <= 0.0 )
    {
        // --- start search at midpoint of lower/upper bounds
        //     if initial value outside of these bounds
        if ( *aout < aLo || *aout > aHi ) *aout = 0.5*(aLo + aHi);

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
        n = findroot_Newton(aLo, aHi, aout, tol, evalContinuity, NULL);

        // --- check if root finder succeeded
        if ( n <= 0 ) n = -1;
    }

    // --- if lower/upper bound functions both negative then use full flow
    else if ( fLo < 0.0 )
    {
        if ( qin > 1.0 ) *aout = ain;
        else *aout = 1.0;
        n = -2;
    }

    // --- if lower/upper bound functions both positive then use no flow
    else if ( fLo > 0 )
    {
        *aout = 0.0;
        n = -3;
    }
    else n = -1;
    return n;
}

//=============================================================================

void evalContinuity(double a, double* f, double* df, void* p)
//
//  Input:   a = outlet normalized area
//  Output:  f = value of continuity eqn.
//           df = derivative of continuity eqn.
//  Purpose: computes value of continuity equation (f) and its derivative (df)
//           w.r.t. normalized area for link with normalized outlet area 'a'.
//
{
    *f  = (Beta1 * xsect_getSofA(pXsect, a*Afull)) + (C1 * a) + C2;
    *df = (Beta1 * Afull * xsect_getdSdA(pXsect, a*Afull)) + C1;
}

//=============================================================================


*/
