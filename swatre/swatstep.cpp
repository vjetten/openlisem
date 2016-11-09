/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
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
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
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
void TWorld::HeadCalc(double *h, bool *ponded,const PROFILE *p,const double  *thetaPrev,
                      const double  *hPrev,const double  *kavg, const double  *dimoca,
                      bool fltsat, double dt, double pond, double qtop, double qbot)
{
   int nN = NrNodes(p);
   const double *dz = Dz(p), *disnod = DistNode(p);
   // dz is layer thickness, distnode is distance between centre of layers
   int i, last = nN-1; // nN nodes from 0 to nN-1 !
   NODE_ARRAY thoma, thomb, thomc, thomf, beta;
   double alpha;

   /* First node : 0 (include boundary cond. qtop or pond) */
	if ( (*ponded) || (fltsat && (qtop <= qbot)) )
	{
		 /* h at soil surface prescribed, ponding */
		 thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
		 thomb[0] = -thomc[0] + dimoca[0] +
					dt*kavg[0]/disnod[0]/dz[0];
		 thomf[0] = dimoca[0]*h[0] +
					dt/(-dz[0]) * (kavg[0] - kavg[1]) +
					dt*kavg[0]*pond/disnod[0]/dz[0];
	}
	else
	{
		 /*  q at soil surface prescribed, qtop = rainfall  */
		 (*ponded) = FALSE;
		 thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
		 thomb[0] = -thomc[0] + dimoca[0];
		 thomf[0] = dimoca[0]*h[0] +
					dt/(-dz[0]) * (- qtop - kavg[1]);
	}


	 /* Intermediate nodes: i = 1 to n-2 */
	for (i = 1; i < nN-1; i++)
	{
		 thoma[i] = -dt*kavg[i]/dz[i]/disnod[i];
		 thomc[i] = -dt*kavg[i+1]/dz[i]/disnod[i+1];
		 thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
		 thomf[i] = dimoca[i]*h[i] +
					dt/(-dz[i])*(kavg[i]-kavg[i+1]);
	}

	 /* last node : n-1 (include boundary cond. qbot)*/
	thoma[last] = -dt*kavg[last]/dz[last]/disnod[last];
	thomb[last] = -thoma[last] + dimoca[last];
	thomf[last] = dimoca[last]*h[last] +
				  dt/(-dz[last])*(kavg[last]+qbot);

	/* Gaussian elimination and backsubstitution h - first time */
	alpha = thomb[0];
	h[0] = thomf[0] / alpha;
	for (i = 1; i < nN; i++) {
		beta[i] = thomc[i-1] / alpha;
		alpha = thomb[i] - thoma[i] * beta[i];
		h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
	}
	for (i = (last-1); i >= 0; i--) /* CW (last-1) was last */
		h[i] -= beta[i+1] * h[i+1];

	/* correct tridiagonal matrix */
	for (i = 0; i < nN; i++) {
		double theta = TheNode(h[i], Horizon(p,i));
		double dimocaNew = DmcNode(h[i], Horizon(p,i));

		thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
		thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i]
				   - theta + thetaPrev[i];
	}

	/* Gaussian elimination and backsubstitution h - second time */
	alpha = thomb[0];
	h[0] = thomf[0] / alpha;
	for (i = 1; i < nN; i++) {
		beta[i] = thomc[i-1] / alpha;
		alpha = thomb[i] - thoma[i] * beta[i];
		h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
	}

	for (i = (last-1); i >= 0; i--) /* CW (last-1) was last */
		h[i] -= beta[i+1] * h[i+1];

   /*


   // First node : 0 (include boundary cond. qtop or pond)
   if ( (*ponded))// || (fltsat && (qtop <= qbot)) )
   {
      // h at soil surface prescribed, ponding
      thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
      thomb[0] = -thomc[0] + dimoca[0] + dt*kavg[0]/disnod[0]/dz[0];
      thomf[0] = dimoca[0]*h[0] + dt/(-dz[0])*(kavg[0] - kavg[1])
                                + dt*kavg[0]*pond/disnod[0]/dz[0];
   }
   else
   {
      //  q at soil surface prescribed, qtop = rainfall (or evap and then qtop positive)
      (*ponded) = false;
      thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
      thomb[0] = -thomc[0] + dimoca[0];
      thomf[0] = dimoca[0]*h[0] + dt/(-dz[0]) * (-qtop - kavg[1]);
   }

   // Intermediate nodes: i = 1 to n-2
   for (i = 1; i < nN-1; i++)
   {
      thoma[i] = -dt*kavg[i]/dz[i]/disnod[i];
      thomc[i] = -dt*kavg[i+1]/dz[i]/disnod[i+1];
      thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
      thomf[i] = dimoca[i]*h[i] + dt/(-dz[i])*(kavg[i]-kavg[i+1]);
      //+dt/(-dz[i])*(Sink[i])
      //add sink term to thomf as flux, Sink must be positive if evap, negative if e.g. drip irrigation
   }

   // last node : n-1 (include boundary cond. qbot)
   thoma[nN-1] = -dt*kavg[nN-1]/dz[nN-1]/disnod[nN-1];
   thomb[nN-1] = -thoma[nN-1] + dimoca[nN-1];
   thomf[nN-1] = dimoca[nN-1]*h[nN-1] + dt/(-dz[nN-1])*(kavg[nN-1]+qbot);

   // Gaussian elimination and backsubstitution h - first time
   alpha = thomb[0];
   h[0] = thomf[0] / alpha;
   for (i = 1; i < nN; i++)
   {
      beta[i] = thomc[i-1] / alpha;
      alpha = thomb[i] - thoma[i] * beta[i];
      h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
   }
   for (i = (nN-2); i >= 0; i--)
   {
      h[i] -= beta[i+1] * h[i+1];
   }

   // correct tridiagonal matrix
   for (i = 0; i < nN; i++)
   {
      double theta = TheNode(h[i], Horizon(p,i));
      double dimocaNew = DmcNode(h[i], Horizon(p,i));

      thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
      thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i] - theta + thetaPrev[i];
   }
   //thomf is new theta

   // Gaussian elimination and backsubstitution h - second time
   alpha = thomb[0];
   h[0] = thomf[0] / alpha;

   for (i = 1; i < nN; i++)
   {
      beta[i] = thomc[i-1] / alpha;
      alpha = thomb[i] - thoma[i] * beta[i];
      h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
      //alpha = tomb = new h/theta
      //thoma = alpha/beta is (new h/theta)/(old
   }

   for (i = nN-2; i >= 0; i--)
   {
      h[i] -= beta[i+1] * h[i+1];
   }

//         for (i = 1; i < nN; i++)
//         {
//            if (h[i] >= -0.001)
//            h[i] = -0.001;
//         }
   // niet doen want node 0 kan positief zijn door waterlaag
*/

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
                             double *drain, double drainfraction, double *repel,
                             double *Theta, SOIL_MODEL *s)
{
   NODE_ARRAY theta, thetaPrev, hPrev, dimoca, kavg, k;
   const PROFILE *p = pixel->profile;
   double *h = pixel->h;
   int    i, n = NrNodes(p);
   //double ThetaSat;
   double dt = pixel->currDt;
   double pond = *waterHeightIO;
   double elapsedTime = 0;
   double influx = 0;
   double drainout = 0;
   int tnode = pixel->tilenode;

   // while internal swatre timestep is not lisem timestep: (time is in seconds!!!)
   while (elapsedTime < _dt)
   {
      bool ponded, fltsat;    // flag idf ponded and if profile fully saturated
      double _max, qtop, qbot, ThetaSat;  // fluxes at top and bottom, max theta
      double qdrain; // tile drainage

      //--- get nodal values of theta, K, dif moist cap ---//
      for (i=0; i < n; i++)
      {
         k[i] = HcoNode(h[i], Horizon(p, i), ksatCalibration);
         // input tables in cm/day function returns in cm/sec !!
         dimoca[i] = DmcNode(h[i], Horizon(p, i));
         // differential moisture capacity dtheta/dh
         theta[i] = TheNode(h[i], Horizon(p, i));
         // moisture content
      }

      *Theta = (theta[0]+theta[1])/2;
      // avg water content of first two nodes, choice ...

      (*repel) = 1.0;
      if (SwitchWaterRepellency)
      {
         if (pixel->repellency == 1)
         {
            *repel = 1/(waterRep_d+pow(waterRep_a, 100*(*Theta-waterRep_b)));
            if (*Theta < waterRep_c)
               *repel = 1.0;
         }
         else
            *repel = 0;
         *repel = std::max(0.0,std::min(1-*repel, 1.0));

         k[0] = k[0] * (*repel);
      }

      // average K for 1st to n-1 node, top node is done below
      //choice arithmetric average K, geometric in org. SWATRE
      if (!SwitchGeometric)
      {
         for (i=1; i < n; i++)
            kavg[i] = 0.5*(k[i]+k[i-1]);
      }
      else
      {
         for(i=1; i < n; i++)
            kavg[i] = sqrt((k[i]*k[i-1]));
      }

      //--- boundary conditions ---//

      //----- TOP -----//
      // check if ponded: 1st compare fluxes, 2nd compare store
      qtop = -pond/dt;
      // top flux is ponded layer / timestep, available water, cm/sec

      //----- BOTTOM -----//
      // bottom is 0 or copy of flux of last 2 layers
      if (SwitchImpermeable)
         qbot = 0;
      else
         qbot = kavg[n-1]*(h[n-1]-h[n-2])/DistNode(p)[n-1] - kavg[n-1];
      //VJ 110122 why not last node? this gives nan, why? => n does not exist!



      // 1st check flux aginst max flux
      ThetaSat = TheNode(0.0, Horizon(p, 0));    
      kavg[0]= sqrt( (*repel) * HcoNode(0, Horizon(p, 0), ksatCalibration) * k[0]);
         
      //kavg[0]= sqrt( (*repel) * HcoNode(0.0, Horizon(p, 0), ksatCalibration) * k[0]);       
      
      // geometric avg of ksat and k[0]
      _max = kavg[0]*(pond-h[0]) / DistNode(p)[0] - kavg[0];
      // maximum possible flux, compare to real top flux available
      ponded = (qtop < _max);
      // if more flux then max possible flag ponded is true
      // NOTE qtop and _max are both negative !

      //2nd check: ponded layer depth against storage
      if (!ponded)
      {
         // calculate available space in profile in cm: (pore-theta)*dz
         double space = 0;
         for(i = 0; i < n && h[i] < 0 /*&& space > pond*/; i++)
         {
            ThetaSat = TheNode(0, Horizon(p, i));
                  //LUT_Highest(p->horizon[i]->lut, THETA_COL);
            space += (ThetaSat - theta[i]) * (-Dz(p)[i]);
         }
         //ponded = pond > space;
         ponded = ((-qtop) * dt) > space;
      }

      // check if profile is completely saturated (flstsat)
      fltsat = true;
      for (i = n-1; i >= 0; i--)
         if (h[i] < 0)
         {
            fltsat = false;
            break;
         }

//      for (i = 0; i < n && h[i] >= 0; i++)
         /* nothing */;
  //    fltsat = (i == n); /* TRUE if forall i h[i] >= 0 */
    //  if (fltsat)
      //   i = 1;

      // save last h and theta, used in headcalc
      for (i = 0; i < n; i++)
      {
         hPrev[i] = h[i];
         thetaPrev[i] = theta[i];
      }

      // calculate hew heads

      HeadCalc(h, &ponded, p, thetaPrev, hPrev, kavg, dimoca,
               fltsat, dt, pond, qtop, qbot);
      // calculate new h and theta with two times gaussian matrix
      // and back substitution abracadabra

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


   } /* elapsedTime < lisemTimeStep */

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
   //*drain = drainout; // tiledrain, is 0 when not activated
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
void TWorld::SwatreStep(int thread, SOIL_MODEL *s, cTMap *_WH, cTMap *_fpot, cTMap *_drain, cTMap *_theta, cTMap *where)
{   
   // map "where" is used as a flag here, it is the fraction of crust, compaction, grass
   // so that the additional calculations are not done everywhere
   // for normal soil surface where is always 1.
   // this prevents doing swatrestep for crusting for cells that are 0 for instance
   FOR_ROW_COL_2DMT
         if(where->Drc > 0) // flag to indicate if this pixel has to be done
         // for regular soil this is 1 so always done, for e.g. crusting only when larger than 0
   {
      double wh, infil, drain, drainfraction = 0, Theta, repellency;

      wh = _WH->Drc*100;
      // WH is in m, convert to cm
      infil = 0;
      drain = 0;


      ComputeForPixel(&s->pixel[r*_nrCols+c], &wh, &infil, &drain, drainfraction,
                      &repellency, &Theta, s);
      // estimate new h and theta at the end of dt

      _WH->Drc = wh*0.01;
      //back to m

      _fpot->Drc = std::max(0.0, -infil*0.01);
      // infil is negative (downward flux * dt, in cm)
      //fpot is positive like in other infil  methods (in m)

      _theta->Drc = Theta;
      // save the average moisture content of the top two layers

      if (SwitchWaterRepellency)
         RepellencyFraction->Drc = repellency;
   }}}
}
//--------------------------------------------------------------------------------



