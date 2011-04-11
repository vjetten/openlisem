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

 note that there is no real iteration: all calculations are done twice
 and if the solution is not good enough, the NEXT SWATRE internal timestep is
 decreased

 functions:
- void TWorld::HeadCalc(double *h,bool *ponded,const PROFILE *p,const double  *thetaPrev, \n
      const double *hPrev,const double  *kavg,const double  *dimoca, \n
      bool fltsat, double dt, double pond, double qtop, double qbot) \n
- double  TWorld::NewTimeStep(double prevDt, const double *hLast, const double *h, \n
  int nrNodes, double precParam, double dtMin, double dtMax) \n
- void TWorld::ComputeForPixel(PIXEL_INFO *pixel, double *waterHeightIO, double *infil,
                             double *drain, SOIL_MODEL *s) \n
- void TWorld::SwatreStep(SOIL_MODEL *s, TMMap *_WH, TMMap *_fpot, TMMap *_drain, TMMap *where) \n
 */
#include "csf.h"
#include "model.h"
#include "swatrelookup.h"



//--------------------------------------------------------------------------------
/*
matrix shape:
|b0 c0 .  .| |h0| |F0|
|a1 .. cn-1|*|  |=|  |
|   an bn  | |hn| |Fn|
*/
void TWorld::HeadCalc(double *h,bool *ponded,const PROFILE *p,const double  *thetaPrev,
				  const double  *hPrev,const double  *kavg,const double  *dimoca,
				  bool fltsat, double dt, double pond, double qtop, double qbot)
{
	int nN = NrNodes(p);
	const double *dz = Dz(p), *disnod = DistNode(p);
   int i; // nN nodes from 0 to nN-1 !
	NODE_ARRAY thoma, thomb, thomc, thomf, beta;
	double alpha;

	// First node : 0 (include boundary cond. qtop or pond)
	if ( (*ponded) || (fltsat && (qtop <= qbot)) )
	{
		// h at soil surface prescribed, ponding
		thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
		thomb[0] = -thomc[0] + dimoca[0] + dt*kavg[0]/disnod[0]/dz[0];
		thomf[0] = dimoca[0]*h[0] + dt/(-dz[0])*(kavg[0] - kavg[1]) + dt*kavg[0]*pond/disnod[0]/dz[0];
		////+ dt/(-dz[0])*Sink[0]
	}
	else
	{
		//  q at soil surface prescribed, qtop = rainfall
		(*ponded) = false;
		thomc[0] = -dt * kavg[1] / dz[0] / disnod[1];
		thomb[0] = -thomc[0] + dimoca[0];
		thomf[0] = dimoca[0]*h[0] + dt/(-dz[0]) * (-qtop - kavg[1]); //+ dt/(-dz[0])*Sink[0]
	}

	// Intermediate nodes: i = 1 to n-2
	for (i = 1; i < nN-1; i++)
	{
		thoma[i] = -dt*kavg[i]/dz[i]/disnod[i];
		thomc[i] = -dt*kavg[i+1]/dz[i]/disnod[i+1];
		thomb[i] = -thoma[i] - thomc[i] + dimoca[i];
		thomf[i] = dimoca[i]*h[i] + dt/(-dz[i])*(kavg[i]-kavg[i+1]);  //+dt/(-dz[i])*(Sink[i])
		//add sink term to thomf as flux, Sink must be negative
	}

	// last node : n-1 (include boundary cond. qbot)
   thoma[nN-1] = -dt*kavg[nN-1]/dz[nN-1]/disnod[nN-1];
   thomb[nN-1] = -thoma[nN-1] + dimoca[nN-1];
   thomf[nN-1] = dimoca[nN-1]*h[nN-1] + dt/(-dz[nN-1])*(kavg[nN-1]+qbot); //+dt/(-dz[nN-1])*(Sink[nN-1])

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
		h[i] -= beta[i+1] * h[i+1];

	// correct tridiagonal matrix
   if (SwitchBacksubstitution)
   {
      for (i = 0; i < nN; i++)
      {
         double theta = TheNode(h[i], Horizon(p,i));
         double dimocaNew = DmcNode(h[i], Horizon(p,i));

         thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
         thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i] - theta + thetaPrev[i];
      }

      // Gaussian elimination and backsubstitution h - second time
      alpha = thomb[0];
      h[0] = thomf[0] / alpha;
      for (i = 1; i < nN; i++)
      {
         beta[i] = thomc[i-1] / alpha;
         alpha = thomb[i] - thoma[i] * beta[i];
         h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
      }

      for (i = nN-2; i >= 0; i--)
         h[i] -= beta[i+1] * h[i+1];

   }

//   for (i = 2; i < nN; i++)
//      h[i] = min(0, h[i]);
   //?????????????????????????
}
//--------------------------------------------------------------------------------
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
		double mdih = accur1 + accur2 * max(1.0, fabs(h[i]));
      double dih  = fabs(h[i] - hLast[i]);
      // if difference is small
      // dih = e.g. 10 and h = -200 then mdih = 200*0.01 + 0.1 = 2.1
      // mdih/dih = 2.1/10 =0.21

		if (dih > 0.10)
			dt = min(dt, prevDt*mdih/dih);
	}
	return (max(dt, dtMin));
}
//--------------------------------------------------------------------------------
// Units are:
// Z and H in cm; table units K in cm/day converted to cm/sec, lisem time in seconds
void TWorld::ComputeForPixel(PIXEL_INFO *pixel, double *waterHeightIO, double *infil,
                             double *drain, SOIL_MODEL *s)
{
	NODE_ARRAY theta, thetaPrev, hPrev, dimoca, kavg, k;
	const PROFILE *p = pixel->profile;
	double *h = pixel->h;
	int    i, n = NrNodes(p);
	double dt = pixel->currDt;
	double pond = *waterHeightIO;
	double elapsedTime = 0;
	double influx = 0;
   double drainout = 0;

   tnode = 5;


	// time is in seconds!!!
   while (elapsedTime < _dt)
	{
		bool ponded, fltsat;
      double qmax, qtop, qbot, ThetaSat, qdrain;

      //--- get nodal values of theta, K, dif moist cap ---//
		for (i=0; i < n; i++)
		{
         k[i] = HcoNode(h[i], Horizon(p, i), ksatCalibration);
         // table in cm/day function returns in cm/sec
			dimoca[i] = DmcNode(h[i], Horizon(p, i));
			theta[i] = TheNode(h[i], Horizon(p, i));
		}

      //--- arithmetric average K, geometric in org. SWATRE ---//
		// average K for 1st to n-1 node
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

		//----- BOTTOM -----//
		// bottom is 0 or copy of flux of last 2 layers
      if (SwitchImpermeable)
			qbot = 0;
		else
         qbot = -kavg[n-1]*((h[n-2]-h[n-1])/DistNode(p)[n-1] + 1);
        // qbot = -kavg[n]*((h[n-1]-h[n])/DistNode(p)[n] + 1);
      //VJ 110122 why not last node? this gives nan, why? => n does not exist!

      // units now in cm/s

		//----- TOP -----//
		qtop = -pond/dt;
		// top flux is ponded layer / timestep, available water, cm/sec

		ThetaSat = TheNode(0.0, Horizon(p, 0));
      kavg[0]= sqrt( HcoNode(ThetaSat, Horizon(p, 0), ksatCalibration) * k[0]);
		// qmax of top node is always calculated with geometric average K
		qmax = -kavg[0]*((h[0]-pond) / DistNode(p)[0] + 1);
//      if (pond == 0)
//         qmax = -kavg[0]*((h[1]-h[0]) / DistNode(p)[0] + 1);
      //VJ 110122
     //actual infil rate Darcy
		//KLOPT eigenlijk niet als niet ponded is pond = 0,ipv een negatieve matrix potentiaal


		ponded = (qtop < qmax);
		// NOTE qtop and qmax are both negative !

		//correct: if more qtop than fits in profile then ponded
		if (!ponded)
		{
			double space = 0;
         for(i = 0; i < n && h[i] < 0 && space > pond; i++)
			{
            //ThetaSat = TheNode(0.0, Horizon(p, i));
            ThetaSat = LUT_Highest(p->horizon[i]->lut, THETA_COL);
				space += (ThetaSat - theta[i]) * (-Dz(p)[i]);
			}
			//ponded = ((-qtop) * dt) > space;
			ponded = pond > space;
		}

		/* check if profile is still completely saturated (flstsat) */
		fltsat = true;
//		for (i = 0; i < n; i++)
         for (i = n-1; i >= 0; i--)
         if (h[i] < 0)
      {
         fltsat = false;
         break;
      }

      //--- calculate tile drain ---//
      if (SwitchIncludeTile)
      {
         //TODO correct this for drainage width: width / _dx
//         qdrain = -kavg[tnode]*( (0-h[tnode])/DistNode(p)[tnode] + 1 );
//         qdrain = min(0, qdrain);
//         if (qdrain < 0)
//         {
//         //   theta[tnode] = theta[tnode] + (qdrain*dt)/DistNode(p)[tnode];
//         //   h[tnode] = HNode(theta[tnode], Horizon(p, tnode));
//         //   drainout += qdrain*dt; // a bit redundant
//         }
      }

		// save last h and theta, used in headcalc
		memcpy(hPrev, h, sizeof(double)*n);
		memcpy(thetaPrev, theta, sizeof(double)*n);

		// calculate hew heads

      HeadCalc(h, &ponded, p, thetaPrev, hPrev, kavg, dimoca,
               fltsat, dt, pond, qtop, qbot);

		// determine new boundary fluxes
      if (SwitchImpermeable)
			qbot = 0;
		else
         qbot = -kavg[n-1]*((h[n-2]-h[n-1])/DistNode(p)[n-1] + 1);

		if ( ponded || (fltsat && (qtop < qbot)) )
			qtop = -kavg[0] * ((h[0] - pond)/DistNode(p)[0] + 1);
		// adjust top flux
		pond += qtop*dt;
		// decrease pond with top flux
		if (pond < POND_EPS)
			pond = 0;
		influx += qmax*dt;
		// add max infil to influx (negative)

		elapsedTime += dt;
		/* estimate new dt within lisemtimestep */
      dt = NewTimeStep(dt, hPrev, h, n, precision, s->minDt, _dt);
		pixel->currDt = dt;

      if (elapsedTime+dt+TIME_EPS >= _dt)
         dt = _dt - elapsedTime;

	} /* elapsedTime < lisemTimeStep */

	/*
*  if (pixel->dumpH>0)
*   printf("Cell %4d : wh after infil. %8.5f (mm) infil. %8.5f (mm)\n"\
*   ,pixel->dumpH,pond*10,-influx*10);
*/
      if (pixel->dumpHid == 1)
         qDebug() << pond << influx << h[0] << h[1] << h[2] << h[3] << h[4] << h[5] << h[6] << h[7] << h[8] << h[9];

	*waterHeightIO = pond;
	*infil = influx; // total max influx in lisem timestep
   *drain = drainout;
}
//--------------------------------------------------------------------------------
// units in SWATRE are cm and cm/day
void TWorld::SwatreStep(SOIL_MODEL *s, TMMap *_WH, TMMap *_fpot, TMMap *_drain, TMMap *where)
{
      SwitchIncludeTile = true;
	FOR_ROW_COL_MV
   if(where->Drc > 0) // when there is crusting for instance
	{
      double wh, infil, drain;

      wh = _WH->Drc*100;
		// WH is in m, convert to cm
		infil = 0;
      drain = 0;

      ComputeForPixel(&s->pixel[r*_nrCols+c], &wh, &infil, &drain, s);
		//->minDt, s->precision, s->calibrationfactor, s->geometric);

      _WH->Drc = wh/100;
		//back to m

      _fpot->Drc = max(0, -infil/100);
		// infil is negative (downward flux * dt, in cm)
		//fpot is positive like in other infil  methods (in m)
//if (r == _nrRows/2 && c == _nrCols/2)
//      qDebug() << drain;
     if (SwitchIncludeTile)
         _drain->Drc = -drain;///_dt*0.01;  // in m/s

	}
   SwitchIncludeTile = false;
}
//--------------------------------------------------------------------------------
// calculates average soil moisture from surface to layernr
/* TODO: NOT CORRECT SHOULD BE RELATIVE TO THICKNESS OF LAYERS */
/*
void SwatreTheta(
		SOIL_MODEL *s,
		MEM_HANDLE *Theta,
		int layernr,
		int avg)
{
	int i = 0;
	layernr -= 1;  //????????
	avg -= 1; //????????

	for (int r = 0; r < Theta->nrRows; r++)
		for (int c = 0; c < Theta->nrCols; c++)
		{
		if (s->pixel[i].h != NULL)
		{
			PIXEL_INFO *p = (s->pixel)+i;
			const PROFILE *prof = p->profile;
			double *h = p->h;
			double theta = 0;
			int j;//, n = NrNodes(prof);  <= to check validity of layernr

			for (j = layernr; j < layernr+avg+1; j++)
				theta += (double)TheNode(h[j], Horizon(prof, j));
			theta /= (avg+1);

			Theta->Data[r][c] = theta;
		}
		i++;
	}
}
	*/
//--------------------------------------------------------------------------------



