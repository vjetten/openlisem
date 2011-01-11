/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

#include "csf.h"
#include "swatre_g.h"
#include "swatre_p.h"
#include "lookup.h"

#include "model.h"

//--------------------------------------------------------------------------------
/*
matrix shape:
|b0 c0 .  .| |h0| |F0|
|a1 .. cn-1|*|  |=|  |
|   an bn  | |hn| |Fn|
*/
void HeadCalc(double *h,bool *ponded,const PROFILE *p,const double  *thetaPrev,
				  const double  *hPrev,const double  *kavg,const double  *dimoca,
				  bool fltsat, double dt, double pond, double qtop, double qbot)
{
	int nN = NrNodes(p);
	const double *dz = Dz(p), *disnod = DistNode(p);
	int i, last = nN-1; // nN nodes from 0 to nN-1 !
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
	thoma[last] = -dt*kavg[last]/dz[last]/disnod[last];
	thomb[last] = -thoma[last] + dimoca[last];
	thomf[last] = dimoca[last]*h[last] + dt/(-dz[last])*(kavg[last]+qbot); //+dt/(-dz[last])*(Sink[last])

	// Gaussian elimination and backsubstitution h - first time
	alpha = thomb[0];
	h[0] = thomf[0] / alpha;
	for (i = 1; i < nN; i++) {
		beta[i] = thomc[i-1] / alpha;
		alpha = thomb[i] - thoma[i] * beta[i];
		h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
	}
	for (i = (last-1); i >= 0; i--)
		h[i] -= beta[i+1] * h[i+1];

	// correct tridiagonal matrix
	for (i = 0; i < nN; i++) {
		double theta = TheNode(h[i], Horizon(p,i));
		double dimocaNew = DmcNode(h[i], Horizon(p,i));

		thomb[i] = thomb[i] - dimoca[i] + dimocaNew;
		thomf[i] = thomf[i] - dimoca[i]*hPrev[i] + dimocaNew*h[i] - theta + thetaPrev[i];
	}

	// Gaussian elimination and backsubstitution h - second time
	alpha = thomb[0];
	h[0] = thomf[0] / alpha;
	for (i = 1; i < nN; i++) {
		beta[i] = thomc[i-1] / alpha;
		alpha = thomb[i] - thoma[i] * beta[i];
		h[i] = (thomf[i] - thoma[i] * h[i-1]) / alpha;
	}

	for (i = (last-1); i >= 0; i--)
		h[i] -= beta[i+1] * h[i+1];

}
//--------------------------------------------------------------------------------
static double  NewTimeStep(
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

		if (dih > 0.10)
			dt = min(dt, prevDt*mdih/dih);
	}
	return (max(dt, dtMin));
}
//--------------------------------------------------------------------------------
// Units are:
// Z and H in cm; table units K in cm/day converted to cm/sec, lisem time in seconds
void ComputeForPixel(PIXEL_INFO *pixel, double *waterHeightIO, double *infil, double *drain,
							double lisDT,SOIL_MODEL *s)
{
	NODE_ARRAY theta, thetaPrev, hPrev, dimoca, kavg, k;
	const PROFILE *p = pixel->profile;
	double *h = pixel->h;
	int    i, n = NrNodes(p);
	double dt = pixel->currDt;
	double pond = *waterHeightIO;
	double elapsedTime = 0;
	double influx = 0;

	// time is in seconds!!!
	while (elapsedTime < lisDT)
	{
		bool ponded, fltsat;
		double qmax, qtop, qbot, ThetaSat;

		//===== get nodal values of theta, K, cap =====//
		for (i=0; i < n; i++)
		{
			k[i] = HcoNode(h[i], Horizon(p, i), s->calibrationfactor, 86400);
			// table in cm/day now in cm/sec
			dimoca[i] = DmcNode(h[i], Horizon(p, i));
			theta[i] = TheNode(h[i], Horizon(p, i));
		}

		//===== arithmetric average K, geometric in org. SWATRE =====//
		// average K for 1st to n-1 node
		if (!s->geometric)
		{
			for (i=1; i < n; i++)
				kavg[i] = 0.5*(k[i]+k[i-1]);
		}
		else
		{
			for(i=1; i < n; i++)
				kavg[i] = sqrt((k[i]*k[i-1]));
		}

		//===== boundary conditions =====//

		//----- BOTTOM -----//
		// bottom is 0 or copy of flux of last 2 layers
		if (s->swatreBottomClosed)
			qbot = 0;
		else
			qbot = -kavg[n-1]*((h[n-2]-h[n-1])/DistNode(p)[n-1] + 1);
		// units now in cm/s

		//----- TOP -----//
		qtop = -pond/dt;
		// top flux is ponded layer / timestep, available water, cm/sec

		ThetaSat = TheNode(0.0, Horizon(p, 0));
		kavg[0]= sqrt( HcoNode(ThetaSat, Horizon(p, 0), s->calibrationfactor, 86400) * k[0]);
		// qmax of top node is always calculated with geometric average K
		qmax = -kavg[0]*((h[0]-pond) / DistNode(p)[0] + 1);
		//actual infil rate Darcy
		//KLOPT eigenlijk niet als niet ponded is pond = 0,ipv een negatieve matrix potentiaal
		// BIJV if (pond == 0) qmax = -kavg[0]*((h[0]-h[1]) / DistNode(p)[0] + 1);

		ponded = (qtop < qmax);
		// NOTE qtop and qmax are both negative !

		//correct: if more qtop than fits in profile then ponded
		if (!ponded)
		{
			double space = 0;
			for(i = 0; i < n && h[i] < 0; i++)
			{
				ThetaSat = TheNode(0.0, Horizon(p, i));
				space += (ThetaSat - theta[i]) * (-Dz(p)[i]);
			}
			//ponded = ((-qtop) * dt) > space;
			ponded = pond > space;
		}

		/* check if profile is still completely saturated (flstsat) */
		//                for (i = 0; i < n && h[i] >= 0; i++)
		//                        /* nothing */;
		//                fltsat = (i == n); /* TRUE if forall i h[i] >= 0 */
		fltsat = true;
		for (i = 0; i < n; i++)
			if (h[i] < 0)
		{
			fltsat = false;
			break;
		}

		// save last h and theta, used in headcalc
		memcpy(hPrev, h, sizeof(double)*n);
		memcpy(thetaPrev, theta, sizeof(double)*n);

		// calculate hew heads
		HeadCalc(h, &ponded, p, thetaPrev, hPrev, kavg, dimoca, /* Sink */
					fltsat, dt, pond, qtop, qbot);
		// determine new boundary fluxes
		if (s->swatreBottomClosed)
			qbot = 0;
		else
			qbot = -kavg[n-1]*((h[n-2]-h[n-1])/DistNode(p)[n-1] + 1);
		//         qbot = kavg[n]*(h[n]-h[n-1])
		//               / DistNode(p)[n] - kavg[n];

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
		dt = NewTimeStep(dt, hPrev, h, n, s->precision, s->minDt, lisDT);
		pixel->currDt = dt;

		if (elapsedTime+dt+TIME_EPS >= lisDT)
			dt = lisDT - elapsedTime;

	} /* elapsedTime < lisemTimeStep */

	/*
*  if (pixel->dumpH>0)
*   printf("Cell %4d : wh after infil. %8.5f (mm) infil. %8.5f (mm)\n"\
*   ,pixel->dumpH,pond*10,-influx*10);
*/
	*waterHeightIO = pond;
	*infil = influx; // total max influx in lisem timestep
}
//--------------------------------------------------------------------------------
// units in SWATRE are cm and cm/day
void TWorld::SwatreStep(SOIL_MODEL *s, TMMap *_WH, TMMap *_fpot, TMMap *_drain, TMMap *where)
{
	FOR_ROW_COL_MV
			if(where->Drc > 0) // when there is crusting for instance
	{
      double wh, infil, drain;

		wh = _WH->Data[r][c]*100;
		// WH is in m, convert to cm
		infil = 0;
      drain = 0;

      ComputeForPixel(&s->pixel[r*nrCols+c], &wh, &infil, &drain, _dt, s);
		//->minDt, s->precision, s->calibrationfactor, s->geometric);
		//DEBUGv(s->pixel[r*nrCols+c].currDt);
		_WH->Data[r][c] = wh/100;
		//back to m
		_fpot->Data[r][c] = max(0, -infil/100);
		// infil is negative (downward flux * dt, in cm)
		//fpot is positive like in other infil  methods (in m)
	}
}
//--------------------------------------------------------------------------------
// calculates average soil moisture from surface to layernr
/** TODO: NOT CORRECT SHOULD BE RELATIVE TO THICKNESS OF LAYERS */
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



