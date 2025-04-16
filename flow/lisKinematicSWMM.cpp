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

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXIT 60

// see page 82 hydraulic SWMM manual part 2,
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


//---------------------------------------------------------------------------
void TWorld::TileFlowSWMM(void)
{
  if (!SwitchIncludeTile && !SwitchIncludeStormDrains)
    return;

  // get water from surface
  if (SwitchIncludeStormDrains) {
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {
      TileWaterVol->Drc += RunoffVolinToTile->Drc;
      // add water from the surface
    }}
  }

  // get water from soil
  if (SwitchIncludeTile) {
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {
      TileWaterVol->Drc += TileDrainSoil->Drc * TileDiameter->Drc * DX->Drc;
      // asume water can come from all sides!
      // add inflow to Tile in m3, tiledrainsoil is in m per timestep

      TileWaterVolSoil->Drc += TileDrainSoil->Drc * TileDiameter->Drc  * DX->Drc;
      // soil only used for MB correction

    }}
  }

  #pragma omp parallel for num_threads(userCores)
  FOR_ROW_COL_MV_TILEL {
    double Area = TileWaterVol->Drc / DX->Drc;
    TileA->Drc = Area;
    double a = Area/TileArea->Drc;
    double theta_next;
    double theta = PI;
    double tol = 1e-6;
    for (int j = 0; j < MAXIT; j++ ) {
       double f = (theta - sin(theta)) / (2 * PI) - a;
       double df = (1 - cos(theta)) / (2 * PI);
       theta_next = theta - f / df;
       if (abs(theta_next - theta) < tol)
           break;
       theta = theta_next;
    }

    double perim = TileDiameter->Drc/2.0*theta_next; // P = r*theta; A =
    if (perim < 1e-6)
        TileQ->Drc = 0;
    else
        TileQ->Drc = std::pow(Area/perim, 5.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;
  }}

  PipeFlowSWMM();

}
//---------------------------------------------------------------------------
//  calculate actual radius r from psi relative (psi = s / dr->sFull;)
double TWorld::psi_rel(double r, double theta) {
    // Partial flow
    double A = 0.5 * r * r * (theta - sin(theta));
    double P = r * theta;
    double R = A / P;

    // Full flow
    double A_full = PI * r * r;
    double R_full = r/2;//  A_full / P_full;

    return (A / A_full) * pow(R / R_full, 2.0 / 3.0);
}
//---------------------------------------------------------------------------
// Newton-Raphson to solve for theta from psi relative (psi = A*R^2/3 unit m^8/3
double TWorld::solve_theta(double r, double psi_target) {
    double theta = PI;
    double tol = 1e-6;
    for (int i = 0; i < MAXIT; ++i) {
        double f = psi_rel(r, theta) - psi_target;
        double df = (psi_rel(r, theta+tol) - psi_rel(r, theta-tol)) / (2*tol);
        double delta = f / df;
        theta -= delta;
        if (fabs(delta) < tol) {
            return theta;
        }
    }
    return theta;
}
//---------------------------------------------------------------------------
//get area A from sectiopn factor s, s = A*R^2/3
double TWorld::getAfromS(DRAIN_PROP *dr, double s)
{
    double psi = s / dr->sFull;
    if (psi == 0.0) return 0.0;
    if (psi >= 1.0) return dr->Afull;

    double r =  dr->diam/2.0;
    double theta = solve_theta(r, psi);

    return (0.5*r*r*(theta - sin(theta)));
}
//---------------------------------------------------------------------------
// do pipe flow according to confined kin wave in SWMM
void TWorld::PipeFlowSWMM()
{
    downstream(crlinkedlddtile_, TileA, tma);
    downstream(crlinkedlddtile_, TileQ, tmb);
    Fill(*Qn,0);

    DRAIN_PROP *drain = new DRAIN_PROP;

    double Qnout=0, Anout=0;

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

        drain->c = c;
        drain->r = r;
        drain->diam = TileDiameter->Drc;
        drain->beta = sqrt(TileGrad->Drc)/TileN->Drc;  // s = qin/beta
        drain->Afull = TileArea->Drc;
        drain->sFull = drain->Afull * std::pow(0.25*TileDiameter->Drc,2.0/3.0);  // 0.5r=0.25D is hydrasulic radius when full
        // section factor = A*R^2/3 units m^2*m^2/3 = m^6/3*m^2/3 = m^8/3
        drain->Qfull = drain->sFull * drain->beta;
        drain->Beta1 = drain->beta / drain->Qfull; // = 1/sFull =>qin/beta1 = qin/(beta/Qfull)
        drain->dxdt = _dx/_dt * drain->Afull / drain->Qfull;
        drain->sMax = 1.08 * drain->sFull;  // circular
        //drain->qin/drain->Beta1

        int    result = 1;
        double dq;
        double WT = 0.6;
        double WX = 0.6;

        // --- normalize previous flows, averrage with downstream for now
        drain->q1 = TileQ->Drc / drain->Qfull;
        drain->q2 = TileQ->Drc / drain->Qfull;
      //  drain->q2 = ((TileQ->Drc + tmb->Drc)*0.5)/ drain->Qfull;
        // --- normalize inflow
        drain->qin = std::min(drain->Qfull, TileQin->Drc)/drain->Qfull;
        // in SWMM code the inflow is maximized to the possible inflow

        // --- compute evaporation and infiltration loss rate
       // double q3 = 0;//link_getLossRate(j, KW, qin*Qfull, tStep) / Qfull;

        // --- normalize previous areas, averrage with downstream
        drain->a1 = TileA->Drc/drain->Afull;
        drain->a2 = TileA->Drc/drain->Afull;
      //  drain->a2 = ((TileA->Drc + tma->Drc)*0.5)/ drain->Afull;

        // --- use full area when inlet flow >= full flow
        if ( drain->qin >= 1.0 )
            drain->ain = 1.0;
        else
            drain->ain = getAfromS(drain, drain->qin/drain->Beta1)/drain->Afull;
        // --- get normalized inlet area corresponding to inlet flow
        //drain->qin/drain->Beta1 = qin/qfull * AR^2/3 / Afull

        // --- check for no flow
        if ( drain->qin < 1e-12 && drain->q2 < 1e-12 ) {
            drain->qout = 0.0;
            drain->aout = 0.0;
        } else {
            dq = drain->q2 - drain->q1;
            drain->C1 = drain->dxdt*WT/WX;
            drain->C2 = (1.0 - WT)*(drain->ain - drain->a1);
            drain->C2 = drain->C2 - WT*drain->a2;
            drain->C2 = drain->C2 * drain->dxdt/WX;
            drain->C2 = drain->C2 + ((1.0 - WX)/WX)*dq - drain->qin;
            //drain->C2   = C2 + q3/WX;

            // --- starting guess for aout is value from previous time step
            drain->aout = drain->a2;

            // --- solve continuity equation for aout
            result = solveContinuity(drain);

            // --- report error if continuity eqn. not solved
            if (result == -1) {
                Error("Kinwave SWMM error solvecontinuity");
            }
            if (result <= 0)
                result = 1;

            // --- compute normalized outlet flow from outlet area
            // polynomial approximation from table for circular pipe
            double a = drain->aout;
            double sfroma = -1.222*a*a*a + 1.9904*a*a + 0.312*a + 0.0039;
            if (a > 0.98)
              sfroma = 1.07662;
            if (a > 0.99)
              sfroma = 1.0;
            drain->qout = drain->Beta1 * sfroma;
            //qout = Beta1 * xsect_getSofA(pXsect, aout*Afull);
            //xsect->sFull * lookup(alpha, S_Circ, N_S_Circ);

           // drain->qout /= drain->Qfull;
            if (drain->qin > 1.0)
                drain->qin = 1.0;

        }

        TileQn->Drc = drain->qout*drain->Qfull;
        TileQn->Drc =  std::min(TileWaterVol->Drc/_dt + drain->qin*drain->Qfull, TileQn->Drc);

        TileWaterVol->Drc = TileWaterVol->Drc + _dt*(drain->qin*drain->Qfull - TileQn->Drc);
        TileWaterVol->Drc = std::max(0.0, TileWaterVol->Drc);
        TileWaterVol->Drc = std::min(TileWaterVol->Drc, TileArea->Drc * DX->Drc);

        Qnout = drain->qout/drain->Qfull;
        Anout = drain->aout;

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
    aLo = std::max(dr->a1,dr->a2)/dr->Afull;
    if (aLo < aHi)
        fLo = ( dr->Beta1 * dr->sMax ) + (dr->C1 * aLo) + dr->C2;
    else
        fLo = fHi;

    // --- if fLo and fHi have same sign then set lower bound to 0
    if (fHi*fLo > 0.0)  {
        aHi = aLo;
        fHi = fLo;
        aLo = 0.0;
        fLo = dr->C2;
    }

    // --- proceed with search for root if fLo and fHi have different signs
    if (fHi*fLo <= 0.0) {
        // --- start search at midpoint of lower/upper bounds
        //     if initial value outside of these bounds
        if (dr->aout < aLo || dr->aout > aHi)
            dr->aout = 0.5*(aLo + aHi);

        // --- if fLo > fHi then switch aLo and aHi
        if (fLo > fHi) {
            aTmp = aLo;
            aLo  = aHi;
            aHi  = aTmp;
        }

        n = findroot_Newton(dr, aLo, aHi);

        // --- check if root finder succeeded
        if ( n <= 0 )
            n = -1;
    } else
        // --- if lower/upper bound functions both negative then use full flow
        if ( fLo < 0.0 ) {
            if ( dr->qin > 1.0 )
                dr->aout = dr->ain;
            else
                dr->aout = 1.0;
            n = -2;
        } else
            // --- if lower/upper bound functions both positive then use no flow
            if ( fLo > 0 ){
                dr->aout = 0.0;
                n = -3;
            } else
                n = -1;
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
    int n = 0;
    double df, dx, dxold, f, x;
    double temp, xhi, xlo;

    // Initialize the "stepsize before last" and the last step.
    x = dr->a2; // first guess
    xlo = x1;
    xhi = x2;
    dxold = fabs(x2-x1);
    dx = dxold;

    double a = x;// /dr->Afull; //relative area
    double s = -1.1148*a*a*a + 1.6721*a*a + 0.4553*a - 0.0063;
    // polynomial fit of alpha and s
    f = dr->Beta1*s + dr->C1*xlo + dr->C2;
    double ds = -0.9519*a*a*a + 1.6474*a*a + 0.4234*a - 0.003;
    df = dr->Beta1*dr->Afull*ds + dr->C1;

    n++;

    // Loop over allowed iterations.
    for (int j=0; j < MAXIT; j++)
    {
        // Bisect if Newton out of range or not decreasing fast enough.
        if (((x-xhi)*df-f)*((x-xlo)*df-f) >= 0.0 || (fabs(2.0*f) > fabs(dxold*df))) {
            dxold = dx;
            dx = 0.5*(xhi-xlo);
            x = xlo + dx;
            if (xlo == x)
                break;
        } else {
            // Newton step acceptable. Take it.
            dxold = dx;
            dx = f/df;
            temp = x;
            x -= dx;
            if (temp == x)
                break;
        }

        // Convergence criterion.
        if (fabs(dx) < EPSILON)
            break;

        // x = relative area
        // in the m anual is a tabular approach with 50 steps but a perfect fit can be made with a polynomial for circular pipes
        double s = -1.1148*x*x*x + 1.6721*x*x + 0.4553*x - 0.0063; // xsect_getSofA
        // polynomial fit of alpha and s
        f = dr->Beta1*s + dr->C1*xlo + dr->C2;
        double ds = -0.9519*x*x*x + 1.6474*x*x + 0.4234*x - 0.003; // xsect_getdSdA
        df = dr->Beta1*dr->Afull*ds + dr->C1;

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

        n++;
        if (f < 0.0)
            xlo = x;
        else
            xhi = x;
    }
    dr->aout = x;
    if (n < MAXIT)
        return n;
    else
        return 0;
}

