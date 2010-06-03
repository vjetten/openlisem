#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "lddutil.h"
#include "chkdata.h"

#include <math.h>      /*pow*/

static    double    epsilon=1E-12; /* iteration epsilon */
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



static double CalcS1(
    double Qj1i1,  /* Qj+1,i+1 : resultaat kin wave voor deze cell ;j = time, i = place */
    double Qj1i,   /* Qj+1,i   : som van alle bovenstroomse kinematic wave  */
    double Qji1,    /* Qj,i+1 Q voor de kinematic wave (t=j) in deze cell, heet Qin in LISEM */
    double q,
    double Sj1i,   /* Sj+1,i : som van alle bovenstroomse sediment */
    double Sji1,    /* Si,j+1 : voor de kinematic wave in de cell, heet Qsedin in LISEM */
    double s,
    double dep,
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
    double deltaX,
    double epsilon)
{
  /* Using Newton-Raphson Method */
    typedef long double REAL;
    REAL Qk1;      /* Q at loop k+1 for i+1, j+1 */
    REAL ab_pQ, deltaTX, C;
    int   count;

    REAL Qkx;
    REAL fQkx;
    REAL dfQkx;
    POSTCOND(sizeof(REAL) > 8);

    /* if no input then output = 0 */
    // MAG NIET q = in m2/s en de anderen in m3/s
//    if ((Qin+Qold+q) <= 0)  /* +q CW NEW! */
//        return(0);
   if((Qin+Qold) <= 0)
     return(0);

//    strcat(ErrorMessage, " KinWave power function");

    /* common terms */
    ab_pQ = alpha*beta*POW_BIG(((Qold+Qin)/2),beta-1);
    deltaTX = deltaT/deltaX;
    C = deltaTX*Qin + alpha*POW_BIG(Qold,beta) + deltaT*q;

    /*  1. Initial guess Qk1.             */
    /*  2. Evaluate function f at Qkx.    */
    /*  3. Evaluate derivative df at Qkx. */
    /*  4. Check convergence.             */

    /*
     * There's a problem with the first guess of Qkx. fQkx is only defined
     * for Qkx's > 0. Sometimes the first guess results in a Qkx+1 which is
     * negative or 0. In that case we change Qkx+1 to 1e-30. This keeps the
     * convergence loop healthy.
     */
    Qkx   = (deltaTX * Qin + Qold * ab_pQ + deltaT * q) / (deltaTX + ab_pQ);
    Qkx   = MAX(Qkx, 1e-20); /* added test-case calc::KinematicTest::iterate1 */
 /* fQkx  = deltaTX * Qkx + alpha * POW_BIG(Qkx, beta) - C;
    dfQkx = deltaTX + alpha * beta * POW_BIG(Qkx, beta - 1);
    Qkx   -= fQkx / dfQkx;
    Qkx   = MAX(Qkx, 1e-30);
  */
    count = 0;
    do {
      fQkx  = deltaTX * Qkx + alpha * POW_BIG(Qkx, beta) - C;   /* Current k */
      dfQkx = deltaTX + alpha * beta * POW_BIG(Qkx, beta - 1);  /* Current k */
      Qkx   -= fQkx / dfQkx;                                /* Next k */
      Qkx   = MAX(Qkx, 1e-20);
      count++;
    } while(fabs(fQkx) > epsilon && count < MAX_ITERS);

    Qk1 = Qkx;
    return(MAX(Qk1,0));
}


static double IterateToQ1(
	double Qin, /* summed Q new in for all sub-cachments */
	double Qold,  /* current discharge */
	double q,
	double alpha,
	double beta,
	double deltaT,
	double deltaX,
	double epsilon)
{
	REAL8 Qk1;      /* Q at loop k+1 for i+1, j+1 */
	REAL8 fQk1;	/* function f  applied to Qk1 so f(Q k+1) */
	REAL8 dfQk1;    /* derivative function f' applied to Qk */
	REAL8 ab_pQ, deltaTX, dk, C;
        int count = 0, maxcount = 10;


       /* if no input then output = 0 */
	if ((Qin+Qold+q) <= 0)
		return(0);

       /* common terms */
	ab_pQ = alpha*beta*pow(((Qold+Qin)/2),beta-1);
	deltaTX = deltaT/deltaX;
    	C = deltaTX*Qin + alpha*pow(Qold,beta) + deltaT*q;

	/* first guess for Qnew with eq 9.6.7 Ven te Chow p296*/
	Qk1 = (deltaTX*Qin + Qold*ab_pQ + deltaT*q)/(deltaTX + ab_pQ);

        fQk1 = deltaTX*Qk1 + alpha*pow(Qk1, beta) - C;
        dfQk1 = deltaTX + alpha*beta*pow(Qk1, beta-1);
        Qk1 -= fQk1/dfQk1;
        if (Qk1 < 1e-30)
           Qk1 = MAX(Qk1, 1e-30);


	do {
	    // function
            fQk1 = deltaTX*Qk1 + alpha*pow(Qk1,beta)-C;
	    // derivative
            dfQk1 = deltaTX + alpha*beta*pow(Qk1,beta-1);
            // newton-rapson
            dk = fQk1/dfQk1;
            Qk1 -= dk;
	    count++;
        } while ((fabs(dk) > Qk1*epsilon) && (count < maxcount) && (Qk1 > epsilon));


        if (Qk1 < epsilon)
  	     return(0);

/*
 *   if (epsilon > 1e-13)
 *   {
 *    printf("ITER Qold %e Qin %e Qnew %e\n",Qold,Qin,Qk1);
 *   }
 */
     // too many iterations
     return(MAX(Qk1,0));
}


static void CalculateKinematic(
      int pitRowNr,    /* r- Y-coordinate of OutflowPoint */
      int pitColNr,    /* r- X-coordinate of OutflowPoint */
      REAL4   **m_Q1, /* -w debiet at time j+1 */
CONST_REAL4_MAP m_Q0, /* debiet at time j   */
CONST_REAL4_MAP m_q,  /* overland water that flows in channel*/
      REAL4   **m_S1, /* -w sediment at time j+1 */
CONST_REAL4_MAP m_S0, /* sediment at time j   */
CONST_REAL4_MAP m_s,  /* sediment that flows in channel*/
CONST_REAL4_MAP m_dep,  /* sediment that flows in channel*/
CONST_UINT1_MAP ldd,
CONST_REAL4_MAP m_alpha, /* see kinematic wave formula */
CONST_REAL4_MAP m_beta,      /* see kinematic wave formula */
      double deltaT,    /*  interval time */
CONST_REAL4_MAP m_deltaX,    /*  deltaX cell width */
      double epsilon)   /* iteration epsilon */
{
/***************************************************************************/
/* this funtion creates a list with a LDD path of all cells draining to the*/
/* cell colNr,rowNr, starting with the outflow point. For each cell in that*/
/* path the fluxes Q1 of the 8 neighbour cells is summed in Qin and a new  */
/* Q1 is estimated. The last cell is the outflow point.                    */
/***************************************************************************/

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
      	        double factor = 0;
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

		correctedDeltaX = m_deltaX[rowNr][colNr];
 		if (applyDiagonalInKinematic && (ldd[rowNr][colNr] % 2)
		    && ldd[rowNr][colNr] != 5 && ldd[rowNr][colNr] != 0 )
            correctedDeltaX *= sqrt(2);

        // for one cell in the list iterate Qj + sum Qupstreamj to Qtimej+1
        m_Q1[rowNr][colNr] =
           (REAL4)IterateToQnew(Qin, m_Q0[rowNr][colNr], qVal,
              m_alpha[rowNr][colNr],m_beta[rowNr][colNr],deltaT,correctedDeltaX,epsilon);

        m_S1[rowNr][colNr] =
           (REAL4)CalcS3(m_Q1[rowNr][colNr], Qin, m_Q0[rowNr][colNr], qVal,
              Sin, m_S0[rowNr][colNr], sVal, m_dep[rowNr][colNr],m_alpha[rowNr][colNr],m_beta[rowNr][colNr],deltaT,correctedDeltaX);


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
          REAL4        *Qout,/* water discharge at outflow point */
          REAL4        *Sout,/* sediment discharge at outflow point */
	      REAL4   **Q1, /* -w discharge at time j+1 */
	CONST_REAL4_MAP Q0, /* discharge at time j   */
	CONST_REAL4_MAP q,  /* overland water that flows in channel*/
	      REAL4   **S1, /* -w sediment at time j+1 */
	CONST_REAL4_MAP S0, /* sediment at time j   */
	CONST_REAL4_MAP s,  /* overland sediment that flows in channel*/
	CONST_REAL4_MAP dep,  /* overland sediment that flows in channel*/
	CONST_UINT1_MAP ldd,
	CONST_REAL4_MAP alpha, /* see kinematic wave formula */
	CONST_REAL4_MAP beta,    /* see kinematic wave formula */
	      double    deltaT,  /*  interval time */
	CONST_REAL4_MAP deltaX,  /*  deltaX cell width */
	      double    epsilon) /* iteration epsilon */
{
	size_t nrRows   =(int) RgiveNrRows();
	size_t nrCols= RgiveNrCols();
	size_t r, c; /* (r,c) becomes co-ordinate of outflowpoint */

        /* mark all cells as not done: */
	for(r = 0; r < nrRows; r++)
	   SetMemMV(Q1[r],nrCols,CR_REAL4);

	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
	  if ( ldd[r][c] == 5 ) {
	   CalculateKinematic(r,c, Q1, Q0, q, S1, S0, s, dep,
	              (const UINT1 **) ldd, alpha, beta,
	                   deltaT, deltaX, epsilon);

	  POSTCOND(!IS_MV_REAL4(Q1[r]+c));
	  POSTCOND(!IS_MV_REAL4(S1[r]+c));

      *Qout = Q1[r][c];
	  *Sout = S1[r][c];
	}
}



/* apply Kinematic wave to compute transport from time j to j+1
 * NOTE the lisem_model used Q2,Q1 where this module does use
 * Q1,Q0, matter of taste
 */
void KineDXalt(
		const char *srcFile,
		int    srcLineNr,
		REAL4 *Qout,      /* 1 -w discharge amount outflow at pit */
		REAL4 *Sout,      /* 2 -w sediment amount outflow at pit */
		MEM_HANDLE *Q1,   /* 3 -w (QNEW) discharge at time j+1 */
		MEM_HANDLE *Q0,   /* 4 r- (QOLD) discharge at time j   */
		MEM_HANDLE *q,    /* 5 r- overland water that flows in channel,
		                   * NULL_MAP if kinematic is not apllied to chanels
		                   */
		MEM_HANDLE *S1,   /* 6 -w sediment amount at time j+1 */
		MEM_HANDLE *S0,   /* 7 r- S0 (QOLD) sediment at time j   */
		MEM_HANDLE *s,    /* 8 r- s  sediment that flows in channel,
				   * NULL_MAP if kinematic is not apllied to chanels */
       	MEM_HANDLE *dep,   /*  */

		MEM_HANDLE *ldd,  /* 9  */
		MEM_HANDLE *alpha,/* 10 r- alpha: see kinematic wave formula */
		MEM_HANDLE *beta,      /*       r- beta:  see kinematic wave formula */
		double deltaT,    /* r- deltaT: interval time */
		MEM_HANDLE *deltaX,    /* r- deltaX: cell width */
		REAL8 epsilon)
/*
 * q has to be:         q(j+1,i+1) + q (j,i+1) / 2
 * alpha has to be:     (n * (b+2*h)^(2/3) / sqrt(S) )^ beta
 * beta has to be for manning's equation:  3/5
 */
{
	REAL4 **Q1Map, **S1Map;
	CONST_REAL4_MAP Q0Map,qMap,alphaMap, betaMap, S0Map,sMap, DXMap, depMap;
	CONST_UINT1_MAP lddMap;

	SetCaller(srcFile,srcLineNr,"kineDXalt");

	/* both q and s have values or are NULL, NOT only one of two: */
	POSTCOND( (s == &NULL_MAP && q == &NULL_MAP)
		|| (s != &NULL_MAP && q != &NULL_MAP));

	Q1Map=(REAL4 **)CpsNormalHandle(Q1,CR_REAL4);
	S1Map=(REAL4 **)CpsNormalHandle(S1,CR_REAL4);

	Q0Map=(CONST_REAL4_MAP)CpsNormalHandle(Q0,CR_REAL4);
	S0Map=(CONST_REAL4_MAP)CpsNormalHandle(S0,CR_REAL4);
	lddMap=(CONST_UINT1_MAP)CpsNormalHandle(ldd,CR_UINT1);
	alphaMap=(CONST_REAL4_MAP)CpsNormalHandle(alpha,CR_REAL4);
	betaMap=(CONST_REAL4_MAP)CpsNormalHandle(beta,CR_REAL4);
	DXMap=(CONST_REAL4_MAP)CpsNormalHandle(deltaX,CR_REAL4);
	depMap=(CONST_REAL4_MAP)CpsNormalHandle(dep,CR_REAL4);

	if (q == &NULL_MAP)
		qMap = sMap = NULL;
	else
	{
		qMap=(CONST_REAL4_MAP)CpsNormalHandle(q,CR_REAL4);
		sMap=(CONST_REAL4_MAP)CpsNormalHandle(s,CR_REAL4);
	}

    KinematicRoutine(
	  	Qout, Sout,
	  	Q1Map,Q0Map,qMap,
	  	S1Map,S0Map,sMap,dep,
		lddMap,alphaMap,betaMap,
        deltaT,DXMap, epsilon);
}


