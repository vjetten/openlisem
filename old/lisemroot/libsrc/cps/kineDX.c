#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "lddutil.h"
#include "chkdata.h"

#include <math.h>      /*pow*/
//VJ 030827 changed from 1e-12 to 1e-6
//static    double    epsilon=1E-6; /* acceptable iteration error in m3/s */
#define MAX_ITERS 30
#define POW_BIG powl
#define MIN_FLUX 1e-20
//extern BOOL applyDiagonalInKinematic;
//BOOL applyDiagonalInKinematic;

#define Q0_ARG     "arg.nr. 4: Q-in"
#define q_ARG      "arg.nr. 5: q"
#define S0_ARG     "arg.nr. 7: S-in"
#define s_ARG      "arg.nr. 8: s"
#define alpha_ARG  "arg.nr. 10: alpha"
#define beta_ARG   "arg.nr. 11: beta"
#define dx_ARG     "arg.nr. 12: dx"

#define _beta 0.6000000000000000
#define _epsilon 1E-6

static double CalcS1(
    double Qj1i1,  /* Qj+1,i+1 : resultaat kin wave voor deze cell ;j = time, i = place */
    double Qj1i,   /* Qj+1,i   : som van alle bovenstroomse kinematic wave  */
    double Qji1,    /* Qj,i+1 Q voor de kinematic wave (t=j) in deze cell, heet Qin in LISEM */
    double q,
    double Sj1i,   /* Sj+1,i : som van alle bovenstroomse sediment */
    double Sji1,    /* Si,j+1 : voor de kinematic wave in de cell, heet Qsedin in LISEM */
    double s,
    double alpha,
    double beta,
    double dt,
    double dx)
{
    double Sj1i1, Cavg, Qavg, dQdt, aQb, abQb_1, A, B, C;
  //   if (Sj1i < 0)
 //      Sj1i = 0;
    // sum all sediment upstream

    if (Qj1i1 < MIN_FLUX)
       return (0);

    Qavg = 0.5*(Qji1+Qj1i);
    if (Qavg <= MIN_FLUX)
       return (0);

    Cavg = (Sj1i+Sji1)/(Qj1i+Qji1);
    aQb = alpha*powl(Qavg,beta);
    abQb_1 = alpha*beta*powl(Qavg,beta-1);

    A = dt*Sj1i;
    B = -dx*Cavg*abQb_1*(Qj1i1-Qji1);
    C = (Qji1 <= MIN_FLUX ? 0 : dx*aQb*Sji1/Qji1);
    Sj1i1 = 1/(dt+dx*aQb/Qj1i1)*(dx*dt*s+A+C+B);

    return max(0,Sj1i1);
}

static double IterateToQnew(
    double Qin, /* summed Q new in for all sub-cachments */
    double Qold,  /* current discharge */
    double q,
    double alpha,
    double beta,
    double deltaT,
    double deltaX)
{
  /* Using Newton-Raphson Method */
    double  ab_pQ, deltaTX, C;  //auxillary vars
    int   count;
    double Qkx; //iterated discharge, becomes Qnew
    double fQkx; //function
    double dfQkx;  //derivative

    /* if no input then output = 0 */
   if ((Qin + Qold) <= 0)
     return(0);

    /* common terms */
    ab_pQ = alpha*beta*POW_BIG(((Qold+Qin)/2),beta-1);
    // derivative of diagonal average (space-time)
    
    deltaTX = deltaT/deltaX;
    C = deltaTX*Qin + alpha*POW_BIG(Qold,beta) + deltaT*q;
    //dt/dx*Q = m3/s*s/m=m2; a*Q^b = A = m2; q*dt = s*m2/s = m2
    //C is unit volume of water

    /*  1. Initial guess Qk1.             */
    /*  2. Evaluate function f at Qkx.    */
    /*  3. Evaluate derivative df at Qkx. */
    /*  4. Check convergence.             */

    // first gues Qkx
    Qkx   = (deltaTX * Qin + Qold * ab_pQ + deltaT * q) / (deltaTX + ab_pQ);

  	// VJ 050704, 060830 infil so big all flux is gone
    if (Qkx < 0) 
      Qkx   = MIN_FLUX;
    
    Qkx   = MAX(Qkx, MIN_FLUX); 
    count = 0;
    do {
      fQkx  = deltaTX * Qkx + alpha * POW_BIG(Qkx, beta) - C;   /* Current k */
      dfQkx = deltaTX + alpha * beta * POW_BIG(Qkx, beta - 1);  /* Current k */
      Qkx   -= fQkx / dfQkx;                                /* Next k */
      Qkx   = MAX(Qkx, MIN_FLUX);
      count++;
    } while(fabs(fQkx) > _epsilon && count < MAX_ITERS);

    return Qkx;
}


static void CalculateKinematic(
      int pitRowNr,    /* r- Y-coordinate of OutflowPoint */
      int pitColNr,    /* r- X-coordinate of OutflowPoint */
      REAL4   **m_Q1, /* -w debiet at time j+1 */
CONST_REAL4_MAP m_Q0, /* debiet at time j   */
//CONST_REAL4_MAP m_q,  /* overland water that flows in channel*/
      REAL4   **m_q, /* -w debiet at time j+1 */
      REAL4   **m_S1, /* -w sediment at time j+1 */
CONST_REAL4_MAP m_S0, /* sediment at time j   */
CONST_REAL4_MAP m_s,  /* sediment that flows in channel*/
      REAL4   **m_StorVol,
      REAL4   **m_StorVolSed,

CONST_UINT1_MAP ldd,
      REAL4 **m_alpha, /* see kinematic wave formula */
CONST_REAL4_MAP m_beta,      /* see kinematic wave formula */
      double deltaT,    /*  interval time */
CONST_REAL4_MAP m_deltaX)    /*  deltaX cell width */
{
/***************************************************************************/
/* this funtion creates a list with a LDD path of all cells draining to the*/
/* cell colNr,rowNr, starting with the outflow point. For each cell in that*/
/* path the fluxes Q1 of the 8 neighbour cells is summed in Qin and a new  */
/* Q1 is estimated. The last cell is the outflow point.                    */
/***************************************************************************/
/*
       int count = 0;
        double h, h1 = 0.001;
        double Q = 0;
        double alpha;
        double _23 = 2.0/3.0;
        double _53 = 5.0/3.0;
        double nS;
        */
 double DeltaX;
 int nrRows, nrCols;
 Liststruct *list, *temp;
 RECMEM_HEAP *heap = NewRecMemHeap(sizeof(Liststruct), 200, NULL, NULL);
 list = (Liststruct *)NewRecord(heap);
 list->prev = NULL;
 list->rowNr = pitRowNr;
 list->colNr = pitColNr;

 nrRows = (int)RgiveNrRows();
 nrCols = (int)RgiveNrCols();

 while (list != NULL)
 {
	int i;
	BOOL  subCachDone = TRUE; /* are sub-catchment cells done ? */
	int rowNr = list->rowNr;
	int colNr = list->colNr;
       /* put all points that have to be calculated to
        * calculate the current point in the list,
	* before the current point
	*/
	for (i=1; i<=9; i++)
	{
		int r,c;

		if (i==5)  /* this is the current cell*/
			continue;
		r = ((int)rowNr)+LddData[i].deltaY;
		c = ((int)colNr)+LddData[i].deltaX;
		if (r>=0 && r<nrRows &&
		    c>=0 && c<nrCols &&
		    FLOWS_TO(ldd,r,c,rowNr, colNr) &&
		    IS_MV_REAL4(&m_Q1[r][c]) ) /* cell not computed */
		{
		    temp = NewRecord(heap);
		    temp->prev = list;
		    list = temp;
		    list->rowNr = r;
		    list->colNr = c;
		    subCachDone = FALSE;
		}
	}


	/* no other points are found that drain to point [rowNr,colNr] */

	if (subCachDone)
	{
		REAL8 Qin=0.0, Sin=0.0;
		REAL8 correctedDeltaX;
		REAL8 qVal=0.0;
		REAL8 sVal=0.0;
//      double factor = 0;
        	/* find for each point in list the neighbour
         	 * LDD points and sum Q1, starting with last point of prev
         	 */
		CheckReal4Cell(m_Q0,rowNr,colNr, Q0_ARG);
		CheckReal4Cell(m_S0,rowNr,colNr, S0_ARG);
		CheckReal4Cell(m_alpha,rowNr,colNr, alpha_ARG);
		CheckReal4Cell(m_beta,rowNr,colNr, beta_ARG);

		if (m_q != NULL)
		{
			CheckReal4Cell(m_q,rowNr,colNr, q_ARG);
			qVal = m_q[rowNr][colNr];
			CheckReal4Cell(m_s,rowNr,colNr, s_ARG);
			sVal = m_s[rowNr][colNr];
		}

		for (i=1;i<=9;i++) /* for all incoming cells */
		{
			int r = ((int)rowNr)+LddData[i].deltaY;
			int c = ((int)colNr)+LddData[i].deltaX;

			if (i==5)  /* Skip current cell */
				continue;

			if (r>=0 && r < nrRows &&
			    c>=0 && c < nrCols &&
			    FLOWS_TO(ldd,r,c,rowNr, colNr) &&
			    !IS_MV_REAL4(&m_Q0[r][c]))
			 {
				Qin += m_Q1[r][c];
				Sin += m_S1[r][c];
          }
		} /* eof all incoming cells */

		temp=list;
   	list=list->prev;
		FreeRecord(temp,heap);
/*
		correctedDeltaX = m_deltaX[rowNr][colNr];
 		if (applyDiagonalInKinematic && (ldd[rowNr][colNr] % 2)
		    && ldd[rowNr][colNr] != 5 && ldd[rowNr][colNr] != 0 )
            correctedDeltaX *= sqrt(2);
*/
		DeltaX = m_deltaX[rowNr][colNr];

      // fill up storage with incoming water
      if (m_StorVol[rowNr][colNr] > 0){
         m_StorVol[rowNr][colNr] -= Qin*deltaT;
         Qin = 0;
         // fill up with sediment
         m_StorVol[rowNr][colNr] -= Sin/2650.000*deltaT;

         m_StorVolSed[rowNr][colNr] += Sin*deltaT;
         Sin = 0;
      }else
      // storage is overflowing = <0
      if (m_StorVol[rowNr][colNr] < 0){
         Qin = -m_StorVol[rowNr][colNr]/deltaT;
         m_StorVol[rowNr][colNr] = 0;
      }
      //else if Qin is 0 then pass on Q

      // for one cell in the list iterate Qj + sum Qupstreamj to Qtimej+1
      m_Q1[rowNr][colNr] =
         (REAL4)IterateToQnew(Qin, m_Q0[rowNr][colNr], qVal,
            m_alpha[rowNr][colNr], _beta, deltaT, DeltaX);
      m_S1[rowNr][colNr] =
         (REAL4)CalcS1(m_Q1[rowNr][colNr], Qin, m_Q0[rowNr][colNr], qVal,
            Sin, m_S0[rowNr][colNr], sVal,
            m_alpha[rowNr][colNr], _beta, deltaT, DeltaX);

      m_q[rowNr][colNr] = (Qin+m_Q0[rowNr][colNr]-m_Q1[rowNr][colNr]);
      //VJ 050831 REPLACE infil with sum all fluxes, needed for infil calculation in main program




           /*
        //Newton iteration of discharge to get water height in meters from total Q
        Q = m_Q1[rowNr][colNr];//Qin + m_Q0[rowNr][colNr] + qVal;
        h1 = 0.001;//+ m_alpha[rowNr][colNr]*powl(Q, _beta)/DeltaX;
        nS = m_beta[rowNr][colNr];
        if (Q > 0) {
          do{
            h = h1;
            // function divided by derivative
            h1 = h - (1-Q/(nS *((powl(DeltaX*h,_53))/
                  (powl(DeltaX+2*h, _23)))))/
                  ((5*DeltaX+6*h)/(3*h*(DeltaX+2*h)));
            count++;
          }while(fabs(h1-h) > 1e-8 && count < 30);
          m_alpha[rowNr][colNr] = powl((nS * powl((DeltaX+2*h1), _23)), _beta);
       }
           */
/*
        factor = 0;
    	if (Qin + m_Q0[rowNr][colNr]+qVal > 0)
           factor = m_Q1[rowNr][colNr]/(Qin + m_Q0[rowNr][colNr]+qVal);
        m_S1[rowNr][colNr] = factor*max(0,Sin+m_S0[rowNr][colNr]+sVal);
*/
		/* cell rowNr, colNr is now done */
	}/* eof subcatchment done */
 } /* eowhile list != NULL */

 FreeAllRecords(heap);
}


static void KinematicRoutine(
	      REAL4   **Q1, /* -w discharge at time j+1 */
	CONST_REAL4_MAP Q0, /* discharge at time j   */
			REAL4 **q,
	      REAL4   **S1, /* -w sediment at time j+1 */
	CONST_REAL4_MAP S0, /* sediment at time j   */
	CONST_REAL4_MAP s,  /* overland sediment that flows in channel*/
	      REAL4   **StoreVol, /* -w sediment at time j+1 */
         REAL4   **StoreVolSed, /* -w sediment at time j+1 */
	CONST_UINT1_MAP ldd,
      	REAL4   **alpha, /* see kinematic wave formula */
	CONST_REAL4_MAP beta,    /* see kinematic wave formula */
	      double    deltaT,  /*  interval time */
	CONST_REAL4_MAP deltaX)  /*  deltaX cell width */
{
	size_t nrRows = RgiveNrRows();
	size_t nrCols = RgiveNrCols();
	size_t r, c; /* (r,c) becomes co-ordinate of outflowpoint */

        /* mark all cells as not done: */
	for(r = 0; r < nrRows; r++)
	   SetMemMV(Q1[r],nrCols,CR_REAL4);

	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
	  if ( ldd[r][c] == 5 ) {
	   CalculateKinematic(r,c, Q1, Q0, q, S1, S0, s, StoreVol,StoreVolSed,
	              (const UINT1 **) ldd, alpha, beta,
	                   deltaT, deltaX);

	  POSTCOND(!IS_MV_REAL4(Q1[r]+c));
	  POSTCOND(!IS_MV_REAL4(S1[r]+c));

	}
}



/* apply Kinematic wave to compute transport from time j to j+1
 * NOTE the lisem_model used Q2,Q1 where this module does use
 * Q1,Q0, matter of taste
 */
void KineDX(
		const char *srcFile,
		int    srcLineNr,
		MEM_HANDLE *Q1,   /* 3 -w (QNEW) discharge at time j+1 */
		MEM_HANDLE *Q0,   /* 4 r- (QOLD) discharge at time j   */
		MEM_HANDLE *q,    /* 5 r- overland water that flows in channel,
		                   * NULL_MAP if kinematic is not apllied to chanels
		                   */
		MEM_HANDLE *S1,   /* 6 -w sediment amount at time j+1 */
		MEM_HANDLE *S0,   /* 7 r- S0 (QOLD) sediment at time j   */
		MEM_HANDLE *s,    /* 8 r- s  sediment that flows in channel,
				   * NULL_MAP if kinematic is not apllied to chanels
				   */
      MEM_HANDLE *StoreVol,
      MEM_HANDLE *StoreVolSed,
		MEM_HANDLE *ldd,  /* 9  */
		MEM_HANDLE *alpha,/* 10 r- alpha: see kinematic wave formula */
		MEM_HANDLE *beta,      /*       r- beta:  see kinematic wave formula */
		double deltaT,    /* r- deltaT: interval time */
		MEM_HANDLE *deltaX)    /* r- deltaX: cell width */
/*
 * q has to be:         q(j+1,i+1) + q (j,i+1) / 2
 * alpha has to be:     (n * (b+2*h)^(2/3) / sqrt(S) )^ 0.6
 * beta has to be for manning's equation:  3/5
 */
{
	REAL4 **Q1Map, **S1Map, **StoreVolMap, **StoreVolSedMap, **alphaMap, **qMap;
	CONST_REAL4_MAP Q0Map,betaMap, S0Map,sMap, DXMap, SnMap;
	CONST_UINT1_MAP lddMap;

	SetCaller(srcFile,srcLineNr,"kineDX");

	/* both q and s have values or are NULL, NOT only one of two: */
	POSTCOND( (s == &NULL_MAP && q == &NULL_MAP)
		|| (s != &NULL_MAP && q != &NULL_MAP));

	Q1Map=(REAL4 **)CpsNormalHandle(Q1,CR_REAL4);
	S1Map=(REAL4 **)CpsNormalHandle(S1,CR_REAL4);
	StoreVolMap=(REAL4 **)CpsNormalHandle(StoreVol,CR_REAL4);
	StoreVolSedMap=(REAL4 **)CpsNormalHandle(StoreVolSed,CR_REAL4);
 	alphaMap=(REAL4 **)CpsNormalHandle(alpha,CR_REAL4);

	Q0Map=(CONST_REAL4_MAP)CpsNormalHandle(Q0,CR_REAL4);
	S0Map=(CONST_REAL4_MAP)CpsNormalHandle(S0,CR_REAL4);
	lddMap=(CONST_UINT1_MAP)CpsNormalHandle(ldd,CR_UINT1);
	betaMap=(CONST_REAL4_MAP)CpsNormalHandle(beta,CR_REAL4);
	DXMap=(CONST_REAL4_MAP)CpsNormalHandle(deltaX,CR_REAL4);

	if (q == &NULL_MAP){
		qMap = NULL;
      sMap = NULL;
   }
	else
	{
		qMap=(REAL4 **)CpsNormalHandle(q,CR_REAL4);//(CONST_REAL4_MAP)CpsNormalHandle(q,CR_REAL4);
		sMap=(CONST_REAL4_MAP)CpsNormalHandle(s,CR_REAL4);
	}

   KinematicRoutine(Q1Map,Q0Map,qMap,
	  	S1Map,S0Map,sMap,
      StoreVolMap,StoreVolSedMap,
      lddMap,alphaMap,betaMap,
      deltaT,DXMap);
}


