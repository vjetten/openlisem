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
    Qkx   = MAX(Qkx, MIN_FLUX); /* added test-case calc::KinematicTest::iterate1 */
 /* fQkx  = deltaTX * Qkx + alpha * POW_BIG(Qkx, beta) - C;
    dfQkx = deltaTX + alpha * beta * POW_BIG(Qkx, beta - 1);
    Qkx   -= fQkx / dfQkx;
    Qkx   = MAX(Qkx, MIN_FLUX);
  */
    count = 0;
    do {
      fQkx  = deltaTX * Qkx + alpha * POW_BIG(Qkx, beta) - C;   /* Current k */
      dfQkx = deltaTX + alpha * beta * POW_BIG(Qkx, beta - 1);  /* Current k */
      Qkx   -= fQkx / dfQkx;                                /* Next k */
      Qkx   = MAX(Qkx, MIN_FLUX);
      count++;
    } while(fabs(fQkx) > epsilon && count < MAX_ITERS);

    Qk1 = Qkx;
    return(MAX(Qk1,0));
}



static void CalculateKinematic(
          int pitRowNr,    /* r- Y-coordinate of OutflowPoint */
          int pitColNr,    /* r- X-coordinate of OutflowPoint */
          REAL4   **m_Qout, /* -w debiet at time j+1 */
    CONST_REAL4_MAP m_Qin, /* debiet at time j   */
    CONST_REAL4_MAP m_Qin1, /* debiet at time j   */
    CONST_REAL4_MAP m_q, /* debiet at time j   */
	      REAL4   **m_So0, /* -w sediment at time j+1 */
	      REAL4   **m_So1, /* -w sediment at time j+1 */
	      REAL4   **m_So2, /* -w sediment at time j+1 */
	CONST_REAL4_MAP m_Si0, /* sediment at time j   */
	CONST_REAL4_MAP m_Si1, /* sediment at time j   */
	CONST_REAL4_MAP m_Si2, /* sediment at time j   */

    CONST_UINT1_MAP ldd,
    CONST_REAL4_MAP m_alpha, /* see kinematic wave formula */
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

 //REAL8 diagonalDeltaX = deltaX*sqrt(2.0);
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
		    IS_MV_REAL4(&m_Qout[r][c]))
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
		REAL8 Qin=0.0, Sin0=0.0, Sin1=0.0,Sin2=0.0;
		REAL8 correctedDeltaX;
       	double factor = 0;
   		REAL8 qVal=0.0;
        double sVal = 0; /* temporary until kin wave is changed */

        	/* find for each point in list the neighbour
         	 * LDD points and sum Q1, starting with last point of prev
         	 */

		CheckReal4Cell(m_Qin,rowNr,colNr, "Qin");
		CheckReal4Cell(m_Si0,rowNr,colNr, "Phos");
		CheckReal4Cell(m_Si1,rowNr,colNr, "NH4");
		CheckReal4Cell(m_Si2,rowNr,colNr, "NO3");
		CheckReal4Cell(m_alpha,rowNr,colNr, "alpha");
		CheckReal4Cell(m_beta,rowNr,colNr, "beta");

		for (i=1;i<=9;i++) /* for all incoming cells */
		{
			int r = ((int)rowNr)+LddData[i].deltaY;
			int c = ((int)colNr)+LddData[i].deltaX;

			if (i==5)  /* Skip current cell */
				continue;

			if (r>=0 && r < nrRows &&
			    c>=0 && c < nrCols &&
			    FLOWS_TO(ldd, r,c, rowNr,colNr) &&
			    !IS_MV_REAL4(&m_Qin[r][c]))// &&
//                ldd[r][c] != MV_UINT1)
			 {
				Qin += m_Qout[r][c];
				Sin0 += m_So0[r][c];
				Sin1 += m_So1[r][c];
				Sin2 += m_So2[r][c];
             }
		} /* eof all incoming cells */

		temp=list;
   		list=list->prev;
		FreeRecord(temp,heap);

   		correctedDeltaX = m_deltaX[rowNr][colNr];
		if (applyDiagonalInKinematic && (ldd[rowNr][colNr] % 2)
		    && ldd[rowNr][colNr] != 5
            && ldd[rowNr][colNr] != 0 )
          correctedDeltaX *= sqrt(2);


         // for one cell in the list iterate Qj + sum Qupstreamj to Qtimej+1

        qVal = m_q[rowNr][colNr];
    	m_Qout[rowNr][colNr] = m_Qin1[rowNr][colNr];
/*
            (REAL4)IterateToQnew(Qin, m_Qin[rowNr][colNr], qVal,
				m_alpha[rowNr][colNr],m_beta[rowNr][colNr],deltaT,correctedDeltaX,epsilon);
*/
        m_So0[rowNr][colNr] =
            (REAL4)CalcS1(m_Qout[rowNr][colNr], Qin, m_Qin[rowNr][colNr], qVal,
                 Sin0, m_Si0[rowNr][colNr], sVal, m_alpha[rowNr][colNr],m_beta[rowNr][colNr],deltaT,correctedDeltaX);
        m_So1[rowNr][colNr] =
            (REAL4)CalcS1(m_Qout[rowNr][colNr], Qin, m_Qin[rowNr][colNr], qVal,
                 Sin1, m_Si1[rowNr][colNr], sVal, m_alpha[rowNr][colNr],m_beta[rowNr][colNr],deltaT,correctedDeltaX);
        m_So2[rowNr][colNr] =
            (REAL4)CalcS1(m_Qout[rowNr][colNr], Qin, m_Qin[rowNr][colNr], qVal,
                 Sin2, m_Si2[rowNr][colNr], sVal, m_alpha[rowNr][colNr],m_beta[rowNr][colNr],deltaT,correctedDeltaX);

   /*

        factor = 0;
    	if (Qin + m_Qin[rowNr][colNr] > 0)
           factor = m_Qout[rowNr][colNr]/(Qin + m_Qin[rowNr][colNr]);
        m_So0[rowNr][colNr] = factor*(Sin0+m_Si0[rowNr][colNr]);
        m_So1[rowNr][colNr] = factor*(Sin1+m_Si1[rowNr][colNr]);
        m_So2[rowNr][colNr] = factor*(Sin2+m_Si2[rowNr][colNr]);
  */
		/* cell rowNr, colNr is now done */
	}/* eof subcatchment done */
 } /* eowhile list != NULL */

 FreeAllRecords(heap);

}


static void KinematicRoutine(
	      REAL4   **Qout,/*  -w discharge at time j+1 */
	CONST_REAL4_MAP Qin, /* discharge at time j   */
	CONST_REAL4_MAP Qin1, /* discharge at time j   */
	CONST_REAL4_MAP q, /* discharge at time j   */
	      REAL4   **So0, /* -w sediment at time j+1 */
	      REAL4   **So1, /* -w sediment at time j+1 */
	      REAL4   **So2, /* -w sediment at time j+1 */
	CONST_REAL4_MAP Si0, /* sediment at time j   */
	CONST_REAL4_MAP Si1, /* sediment at time j   */
	CONST_REAL4_MAP Si2, /* sediment at time j   */
	CONST_UINT1_MAP ldd,
	CONST_REAL4_MAP alpha, /* see kinematic wave formula */
	CONST_REAL4_MAP beta,    /* see kinematic wave formula */
	      double    deltaT,  /*  interval time */
	CONST_REAL4_MAP deltaX)  /*  deltaX cell width */
{
	size_t nrRows   =(int) RgiveNrRows();
	size_t nrCols= RgiveNrCols();
	size_t r, c; /* (r,c) becomes co-ordinate of outflowpoint */

        /* mark all cells as not done: */
	for(r = 0; r < nrRows; r++)
	   SetMemMV(Qout[r],nrCols,CR_REAL4);

	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
	  if ( ldd[r][c] == 5 )
      {
	   CalculateKinematic(r,c, Qout, Qin, Qin1, q, 
                         So0,So1,So2,
                         Si0,Si1,Si2,
                         (const UINT1 **) ldd, alpha, beta,
	                     deltaT, deltaX);
      }
}



/* apply Kinematic wave to compute transport from time j to j+1
 * NOTE the lisem_model used Q2,Q1 where this module does use
 * Q1,Q0, matter of taste
 */

void Ninematic(
	const char *srcFile, int srcLineNr,
	MEM_HANDLE *Qout, MEM_HANDLE *Qin, MEM_HANDLE *Qin1, MEM_HANDLE *q,
    MEM_HANDLE *So0, MEM_HANDLE *So1, MEM_HANDLE *So2,
    MEM_HANDLE *Si0, MEM_HANDLE *Si1, MEM_HANDLE *Si2,
	MEM_HANDLE *ldd, MEM_HANDLE *alpha, MEM_HANDLE *beta,
        double deltaT, MEM_HANDLE *deltaX)
{
	REAL4 **So0Map, **So1Map, **So2Map, **QoutMap;
	CONST_REAL4_MAP DXMap, Qin1Map, QinMap, Si0Map, Si1Map, Si2Map, qMap, alphaMap, betaMap;
	CONST_UINT1_MAP lddMap;

	SetCaller(srcFile,srcLineNr,"ninematic");

	QoutMap=(REAL4 **)CpsNormalHandle(Qout,CR_REAL4);
	So0Map=(REAL4 **)CpsNormalHandle(So0,CR_REAL4);
	So1Map=(REAL4 **)CpsNormalHandle(So1,CR_REAL4);
	So2Map=(REAL4 **)CpsNormalHandle(So2,CR_REAL4);

	QinMap=(CONST_REAL4_MAP)CpsNormalHandle(Qin,CR_REAL4);
	Qin1Map=(CONST_REAL4_MAP)CpsNormalHandle(Qin1,CR_REAL4);
	qMap=(CONST_REAL4_MAP)CpsNormalHandle(q,CR_REAL4);
	Si0Map=(CONST_REAL4_MAP)CpsNormalHandle(Si0,CR_REAL4);
	Si1Map=(CONST_REAL4_MAP)CpsNormalHandle(Si1,CR_REAL4);
	Si2Map=(CONST_REAL4_MAP)CpsNormalHandle(Si2,CR_REAL4);

	lddMap=(CONST_UINT1_MAP)CpsNormalHandle(ldd,CR_UINT1);
	alphaMap=(CONST_REAL4_MAP)CpsNormalHandle(alpha,CR_REAL4);
	betaMap=(CONST_REAL4_MAP)CpsNormalHandle(beta,CR_REAL4);
	DXMap=(CONST_REAL4_MAP)CpsNormalHandle(deltaX,CR_REAL4);

    KinematicRoutine(
      QoutMap,QinMap,Qin1Map,q,
      So0Map, So1Map, So2Map,
      Si0Map, Si1Map, Si2Map,
      lddMap,alphaMap,betaMap,deltaT,DXMap);

}



