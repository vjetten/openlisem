#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "lddutil.h"
#include "chkdata.h"

#include <math.h>      /*pow*/
//VJ 030827 changed from 1e-12 to 1e-6
//static    double    epsilon=1E-6; /* acceptable iteration error in m3/s */
#define MAX_ITER 30
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

void SortLdd (MAP *ldd)
{
    REAL4 value  = 0;
    REAL4 value1 = 0;
    REAL4 value2 = 0;
    REAL4 value3 = 0;
    REAL4 value4 = 0;
    REAL4 value6 = 0;
    REAL4 value7 = 0;
    REAL4 value8 = 0;
    REAL4 value9 = 0;
    REAL4 flag;
    REAL4 checked = nodata;
    MAP *upstr;
    matrix *temp;
    size_t r,c;
    UINT4 nr = RgetNrRows(ldd);
    UINT4 nc = RgetNrCols(ldd);
    temp = new matrix(nr, nc); //init matrix
    sorted map2array; //instantiate the array to hold the sorted map
    sorted NotSortedArray; //array to hold the non-mv cells in "normal" order

    //------------------------------

    //calculates how many cells are in the ldd area. Needed for stop condition
    UINT4 counter = 0;
    for(r=1; r < nr-1; r++)
    {
     for(c=1; c < nc-1; c++)
     {
      value = temp->matrix[r][c];
      if (value >=0 && value <=9)
      {
          NotSortedArray.cell.push_back(cells(r, c, value));
          counter ++;
      }
     }
    }
    nonsorted = NotSortedArray;
    //-----------------------------------
      //-----------------------------------

    //this loop sorts the map filling in the array with the right order of the
    //cells and creates a coloured map with the incremental order of the cells.
    //It stops when the array has the same number of dll cells than the map
    UINT4 x = 0;
    UINT4 i = 0;
    do
    {
    x++;
      for(r=1; r < nr-1; r++)
      {
          for(c=1; c < nc-1; c++)
          {
          value = temp->matrix[r][c];
              if (value == checked)continue;
          value7 = temp->matrix[r-1][c-1];
          value8 = temp->matrix[r-1][c];
          value9 = temp->matrix[r-1][c+1];
          value4 = temp->matrix[r][c-1];
          value6 = temp->matrix[r][c+1];
          value1 = temp->matrix[r+1][c-1];
          value2 = temp->matrix[r+1][c];
          value3 = temp->matrix[r+1][c+1];
              if (value7 != 3 &&
                  value8 != 2 &&
                  value9 != 1 &&
                  value4 != 6 &&
                  value6 != 4 &&
                  value1 != 9 &&
                  value2 != 8 &&
                  value3 != 7 &&
                  value != nodata)
              {
               flag = x;
               //VJ 040616: commented out, it works
               //RputCell(upstr, r, c, &flag); //temporary stuff?. Remove when ready?
               map2array.cell.push_back(cells(r,c,value));
              }
           }
      }
      for(; i < map2array.cell.size(); i++)
      {
          r = map2array.cell[i].row;
          c = map2array.cell[i].col;
          temp->matrix[r][c] = (REAL4)checked;
      }
    }while(map2array.cell.size()< counter);
//VJ 040616: commented out, it works
//  Mclose(upstr);
    delete temp;
    return map2array;
}


static void CalculateKinematic(
      int pitRowNr,   
      int pitColNr,   
      REAL4   **m_Q1, 
CONST_REAL4_MAP m_Q0, 
      REAL4   **m_q, 
      REAL4   **m_S1,
CONST_REAL4_MAP m_S0,
CONST_REAL4_MAP m_s, 
      REAL4   **m_StorVol,
      REAL4   **m_StorVolSed,
CONST_UINT1_MAP ldd,
      REAL4 **m_alpha, 
CONST_REAL4_MAP m_beta,
      double dt,   
CONST_REAL4_MAP m_deltaX)
{
 double dx;
 int r, c, nrRows, nrCols;
  REAL8 Qji1, dtdx, abQ, avgQ, Qj1i, a, b, smallQ, dxc;
  REAL8 C, fQj1i1, dfQj1i1,Qk, Qk1;
  REAL8 Deltads = 0; //keeps the change in depression storage to add to small q
  REAL8 Infilt = 0; //keeps the exact amount of reinfiltration to update
                //infilt potential so we can send it to the infiltDeprStor function

 nrRows = (int)RgiveNrRows();
 nrCols = (int)RgiveNrCols();
 
 for (r = 0; r < nrRows; r++)
   for (c = 0; c < nrCols; c++)
   if(!IS_MV_REAL4(&m_Q0[r][c]))
   {
        int d = (int)ldd[r][c];
        int count = 0;

        
        switch (d) { //add de previously calculated discharge to the downstream cell
          case 1:   if (!IS_MV_REAL4(&m_Q1[r+1][c-1])) m_Q1[r+1][c-1]+= m_Q0[r][c];else continue; break;
          case 2:   if (!IS_MV_REAL4(&m_Q1[r+1][c]  )) m_Q1[r+1][c]  += m_Q0[r][c];else continue; break;
          case 3:   if (!IS_MV_REAL4(&m_Q1[r+1][c+1])) m_Q1[r+1][c+1]+= m_Q0[r][c];else continue; break;
          case 4:   if (!IS_MV_REAL4(&m_Q1[r][c-1]  )) m_Q1[r][c-1]  += m_Q0[r][c];else continue; break;
          case 5:   m_Q1[r][c] += 0;      break;//if it is the outlet simply send the value to hell
          case 6:   if (!IS_MV_REAL4(&m_Q1[r][c+1]  )) m_Q1[r][c+1]  += m_Q0[r][c];else continue; break;
          case 7:   if (!IS_MV_REAL4(&m_Q1[r-1][c-1])) m_Q1[r-1][c-1]+= m_Q0[r][c];else continue; break;
          case 8:   if (!IS_MV_REAL4(&m_Q1[r-1][c]  )) m_Q1[r-1][c]  += m_Q0[r][c];else continue; break;
          case 9:   if (!IS_MV_REAL4(&m_Q1[r-1][c+1])) m_Q1[r-1][c+1]+= m_Q0[r][c];else continue; break;
          default: return -1;
       }
       dx = m_deltaX[r][c];
       dtdx = dt/dx;
   
       Qj1i = m_Q1[r][c]; //discharge in m3/s j = timstep so j1i is total upstream
                                 // coming in, result at t+1
       Qji1 = m_Q0[r][c]; //discharge in m3/s ji1 = current at the cell itself, at t and x=i+1
       smallQ = m_q[r][c];  //discharge in m2/s
   
       b = 0.6000000;
   /*    
       if (!Is_ChFlow){
         //calculates alpha, depresion storage and flowwidth for overlandflow
         Infilt = 0;
         a = AlphaAndDeprStor(Qj1i, Qji1, smallQ, Deltads, Infilt, r, c);
         m_q[r][c] += Infilt;
         //VJ 040705: subtract "to channel flow" from smallQ
         smallQ = smallQ - _ChanOFdischarge->matrix[r][c]/dx - Deltads/(dt*dx);
       }else
         //calculates alpha, waterheight and depression storage for channel flow
         a = ChannelAlpha(Qj1i, Qji1, smallQ, r, c);
   */
   smallQ = m_q[r][c];
   a = m_alpha[r][c];

       avgQ = (Qji1+Qj1i)/2; //check behaviour if it is the outlet!!!!!

       if (avgQ==0) abQ=0;
       else abQ = a*b*powl(avgQ, b-1);

       Qk = ((dtdx*Qj1i)+(abQ*Qji1)+(dt*smallQ))/(dtdx+abQ);
       //so far was the linear solution which will be used as the initial guess
       //for the non-linear solution using newton-raphson iterations

       if (Qk <= 0) {
          m_Q1[r][c]=0;
          continue;
       } //account for infiltration and stuff

       C = dtdx*Qj1i+a*powl(Qji1, b)+dt*smallQ;
       //now we add the initial guess Qk
       Qk1 = Qk;
       do{
         Qk=Qk1;
         fQj1i1 = dtdx*Qk+a*powl(powl(Qk, 6.0),0.1)-C;
         dfQj1i1 = dtdx+a*b*powl(powl(Qk, 2.0),-1.0/5.0);
         Qk1 = Qk - (fQj1i1/dfQj1i1);
         //Qk1 = max(double(Qk1), MIN_FLUX); //avoids powl illegal operation
         fQj1i1 = dtdx*Qk+a*powl(powl(Qk, 6.0),0.1)-C;
         count++;
       }while(fabs(fQj1i1)>_epsilon && count < MAX_ITER);
   
       if (Qk1 <= MIN_FLUX) Qk1 = 0;
       m_Q1[r][c]= Qk1;
   }
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


