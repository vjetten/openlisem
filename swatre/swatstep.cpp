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


/*!
  \file swatstep.cpp
  \brief SWATRE iteration and soil matrix solution

 Note that there is no real iteration: all calculations are done twice
 and if the solution is not good enough, the NEXT SWATRE internal timestep is
 decreased.
 Order of calculations:\
 LisInfiltration, swatre -> SwatreStep --> for each pixxel do: ComputeForPixel HeadCalc + NewTimestep
 functions:
- void TWorld::HeadCalc(double *h,bool *ponded,const PROFILE *p,const double  *thetaPrev, \n
      const double *hPrev,const double  *kavg,const double  *dimoca, \n
      bool fltsat, double dt, double pond, double qtop, double qbot) \n
- double  TWorld::NewTimeStep(double prevDt, const double *hLast, const double *h, \n
  int nrNodes, double precParam, double dtMin, double dtMax) \n
- void TWorld::ComputeForPixel(PIXEL_INFO *pixel, double *waterHeightIO, double *infil,
                             double *drain, SOIL_MODEL *s) \n
- void TWorld::SwatreStep(SOIL_MODEL *s, cTMap *_WH, cTMap *_fpot, cTMap *_drain, cTMap *where) \n
 */
#include <algorithm>
#include "model.h"
#include "swatreLookup.h"




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
/**
 * @brief
 *
 * @param h
 * @param ponded
 * @param p
 * @param thetaPrev
 * @param hPrev
 * @param kavg
 * @param dimoca
 * @param fltsat
 * @param dt
 * @param pond
 * @param qtop
 * @param qbot
 */
void TWorld::HeadCalc(double *h, bool *ponded,const PROFILE *p,const double *thetaPrev,
                      const double *hPrev, const double *kavg, const double *dimoca,
                      bool fltsat, double dt, double pond, double qtop, double qbot)
{
   int nN = NrNodes(p);
   //const double *dz = Dz(p), *disnod = DistNode(p);
   // dz is layer thickness, distnode is distance between centre of layers
   int i;
   NODE_ARRAY thoma, thomb, thomc, thomf, beta;
   double alpha;

   /* First node : 0 (include boundary cond. qtop or pond) */
    if ( (*ponded) || (fltsat && (qtop <= qbot)) ) {
		 /* h at soil surface prescribed, ponding */
         thomc[0] = -dt * kavg[1] / p->zone->dz[0] / p->zone->disnod[1];
         thomb[0] = -thomc[0] + dimoca[0] + dt*kavg[0]/p->zone->disnod[0]/p->zone->dz[0];
         thomf[0] = dimoca[0]*h[0] + dt/(-p->zone->dz[0]) * (kavg[0] - kavg[1]) +
                 dt*kavg[0]*pond/p->zone->disnod[0]/p->zone->dz[0];
    } else {
		 /*  q at soil surface prescribed, qtop = rainfall  */
         (*ponded) = false;
         thomc[0] = -dt * kavg[1] / p->zone->dz[0] / p->zone->disnod[1];
		 thomb[0] = -thomc[0] + dimoca[0];
         thomf[0] = dimoca[0]*h[0] + dt/(-p->zone->dz[0]) * (- qtop - kavg[1]);
	}


	 /* Intermediate nodes: i = 1 to n-2 */
	for (i = 1; i < nN-1; i++)
	{
         thoma[i] = -dt*kavg[i]/p->zone->dz[i]/p->zone->disnod[i];
         thomc[i] = -dt*kavg[i+1]/p->zone->dz[i]/p->zone->disnod[i+1];
		 thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
         thomf[i] = dimoca[i]*h[i] + dt/(-p->zone->dz[i])*(kavg[i]-kavg[i+1]);
	}

    // last node : nN-1 (include boundary cond. qbot)
    thoma[nN-1] = -dt*kavg[nN-1]/p->zone->dz[nN-1]/p->zone->disnod[nN-1];
    thomb[nN-1] = -thoma[nN-1] + dimoca[nN-1];
    thomf[nN-1] = dimoca[nN-1]*h[nN-1] + dt/(-p->zone->dz[nN-1])*(kavg[nN-1]+qbot);

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
   double dt = dtMax;
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
/**
 * @brief
 *
 * @param pixel
 * @param waterHeightIO
 * @param infil
 * @param drain
 * @param drainfraction
 * @param repel
 * @param Theta
 * @param s
 */


void TWorld::ComputeForPixel(PIXEL_INFO *pixel, double *waterHeightIO, double *infil,
                             double *drain, double drainfraction,
                             double *Theta, SOIL_MODEL *s)
{
   NODE_ARRAY theta, thetaPrev, hPrev, dimoca, kavg, k;
   const PROFILE *p = pixel->profile;
   int n = pixel->nrNodes;//p->zone->nrNodes;
   double dt = _dt/5;// std::max(_dt/5, pixel->currDt);
   double pond = *waterHeightIO;
   double elapsedTime = 0;
   double influx = 0;
   double drainout = 0;
   int tnode = pixel->tilenode;

   if (SHOWDEBUG)
       qDebug() << "compute for pixel" << n << pixel->MV << p->profileId;

   double *h = new double [n];
   for (int i = 0; i < n; i++) {
       h[i] = pixel->h[i];
   }


   while (elapsedTime < _dt)
   {

      bool ponded, fltsat;    // flag if ponded or if profile fully saturated
      double _max, qtop, qbot, ThetaSat;  // fluxes at top and bottom, max theta
      double qdrain; // tile drainage

      // get nodal values of theta, K, dif moist cap
      for (int j = 0; j < pixel->nrNodes; j++) {
         k[j] = HcoNode(h[j], p->horizon[j], ksatCalibration);
         // input tables in cm/day function returns in cm/sec !!
         dimoca[j] = DmcNode(h[j], p->horizon[j]);
         // differential moisture capacity d(theta)/d(h), tangent moisture retention curve
         theta[j] = TheNode(h[j], p->horizon[j]);
         // moisture content
        // qDebug() << h[j] << k[j] << theta[j] << dimoca[j];
      }     
      if (SHOWDEBUG) {
         qDebug() << elapsedTime <<dt;
          qDebug() << "h1" << h[0] << h[2] << h[3] << h[4] << h[5] << h[6] << h[7] << h[8] << h[9];
      }
      *Theta = (theta[0]+theta[1])/2;
      // avg water content of first two nodes, choice ...
      //TODO: WHY? FOR WHAT, pesticides, repelency?

      // average K for 1st to n-1 node, top node is done below
      //choice arithmetric average K, geometric in org. SWATRE

      for(int j = 1; j < nNodes; j++) {
        kavg[j] = std::sqrt(k[j-1]*k[j]);
        kavg[0] = kavg[1];

          // switch (KavgType) {
          // case 0: kavg[j] = Aavg(k[j-1],k[j]); break;
          // case 1: kavg[j] = Savg(k[j-1],k[j]); break;
          // case 2: kavg[j] = Havg(k[j-1],k[j],p->zone->disnod[j-1],p->zone->disnod[j]); break;
          // case 3: kavg[j] = Mavg(k[j-1],k[j]); break;
          // }
      }

      //--- boundary conditions ---

      //----- TOP -----
      // check if ponded: 1st compare fluxes, 2nd compare store
      // top flux is ponded layer / timestep, available water/rainfall, cm/sec
      qtop = -pond/dt;

      //----- BOTTOM -----
      // bottom is 0 or copy of flux of last 2 layers
      if (SwitchImpermeable)
         qbot = 0;
      else
         qbot = kavg[n-1]*(h[n-1]-h[n-2])/DistNode(p)[n-1] - kavg[n-1];

      // 1st check flux against max flux
      ThetaSat = TheNode(0.0, Horizon(p, 0));
      //kavg[0]= sqrt( (*repel) * HcoNode(0, Horizon(p, 0), ksatCalibration) * k[0]);
      kavg[0]= HcoNode(0, p->horizon[0], ksatCalibration);//sqrt( HcoNode(0, p->horizon[0], ksatCalibration) * k[0]);
      // geometric avg of ksat and k[0] => is used for max possible
     // _max = -kavg[0]*((h[0]-pond) / DistNode(p)[0] + 1);
      _max = kavg[0]*(pond-h[0]) / DistNode(p)[0] - kavg[0];
//      qmax = -kavg[0]*((h[0]-pond) / DistNode(p)[0] + 1);
      // maximum possible flux, compare to real top flux available
      ponded = (qtop < _max);
      // if more flux then max possible flag ponded is true
      // NOTE qtop and _max are both negative !
      if (SHOWDEBUG)
        qDebug() << "swatre" <<ponded << _max << qtop << DistNode(p)[0] << DistNode(p)[1];

      //2nd check: ponded layer depth against storage
      if (!ponded)
      {
         // calculate available space in profile in cm: (pore-theta)*dz
         double space = 0;
         for(int i = 0; i < n && h[i] < 0 /*&& space > pond*/; i++)
         {
            ThetaSat = TheNode(0, Horizon(p, i));
            space += (ThetaSat - theta[i]) * (-Dz(p)[i]);
         }
         ponded = pond > space;
      }

      // check if profile is completely saturated (flstsat)
      fltsat = true;
      for (int i = n-1; i >= 0; i--)
         if (h[i] < 0)
         {
            fltsat = false;
            break;
         }

      // save last h and theta, used in headcalc
      for (int i = 0; i < n; i++)
      {
         hPrev[i] = h[i];
         thetaPrev[i] = theta[i];
      }

      // calculate hew heads

      HeadCalc(h, &ponded, p, thetaPrev, hPrev, kavg, dimoca, fltsat, dt, pond, qtop, qbot);
      // calculate new h and theta with two times gaussian matrix
      // and back substitution abracadabra
 if (SHOWDEBUG) {
     qDebug() << "h" << h[0] << h[2] << h[3] << h[4] << h[5] << h[6] << h[7] << h[8] << h[9];
 }
      // determine new boundary fluxes

      if (SwitchImpermeable)
         qbot = 0;
      else
         qbot = kavg[n-1]*(h[n-1]-h[n-2])/DistNode(p)[n-1] - kavg[n-1];
      // new qbot is actually not use but may come in handy later

      if ( ponded || (fltsat && (qtop < qbot)) )
         qtop = -kavg[0] * ((h[0] - pond)/DistNode(p)[0] + 1);

      // adjust top flux

      pond += qtop*dt;
      // decrease pond with top flux
      if (pond < POND_EPS)  // 10-6 cm
         pond = 0;
      influx += _max*dt;
      // add max infil to influx (negative), to get potential infil

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
      dt = NewTimeStep(dt, hPrev, h, n, precision, s->minDt, _dt);
      pixel->currDt = dt;

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

   // save stuff for output in maps
  *waterHeightIO = pond; // waterlayer on surface
  *infil = influx; // total max influx in lisem timestep, fpot
  *drain = drainout; // tiledrain, is 0 when not activated

   for (int i = 0; i < n; i++) {
       pixel->h[i] = h[i];
   }
   delete[] h;

       if (SHOWDEBUG)
           qDebug() << "pond" << pond << influx;
}
//--------------------------------------------------------------------------------
// units in SWATRE are cm and cm/day
/**
 * @brief
 *
 * @param s
 * @param _WH
 * @param _fpot
 * @param _drain
 * @param _theta
 * @param where
 */
void TWorld::SwatreStep(int step, int r, int c, SOIL_MODEL *s, cTMap *_WH, cTMap *_fpot, cTMap *_drain, cTMap *_theta)
{   
    double wh, infil, drain, drainfraction = 0, Theta;
    QString dig;

    wh = _WH->Drc*100;
    // WH is in m, convert to cm
    infil = 0;
    drain = 0;
    showr = r;
    showc = c;

    long cel = r*_nrCols+c;

     // for (int i = 0; i < s->pixel[cel].nrNodes; i++) {
     //     qDebug() << i << s->pixel[cel].profile->horizon[i]->lut->hydro[H_COL];
     // }

    if (SwitchIncludeTile)
        drainfraction = TileWidth->Drc/_dx;

    ComputeForPixel(&s->pixel[cel], &wh, &infil, &drain, drainfraction,  &Theta, s);
    // estimate new h and theta at the end of dt

    //SwitchDumpH = true;
    if(SwitchDumpH || SwitchDumpTheta || SwitchDumpK) {
        if(s->pixel[cel].dumpHid > 0) {
            for (int i = 0; i < s->pixel[cel].nrNodes; i++) {
                QString name = QString("SwH%1").arg(step,2, 10, QLatin1Char('0'));
                dig = QString("%1").arg(i, 12-name.length(), 10, QLatin1Char('0'));
                name=name+dig;
                name.insert(8, ".");
                //qDebug() << name << dig;
            }
        }
    }


    _WH->Drc = wh*0.01;
    //back to m

    _fpot->Drc = std::max(0.0, -infil*0.01);
    // infil is negative (downward flux * dt, in cm)
    //fpot is positive like in other infil  methods (in m)

    if (SwitchIncludeTile)
        _drain->Drc = drain*0.01;  // in m
    // drained water from the soil, already accounts for drainwidth versus cell width

    //_theta->Drc = Theta;
    // save the average moisture content of the top two layers
    // used for repellency OBSOLETE
}
//--------------------------------------------------------------------------------



