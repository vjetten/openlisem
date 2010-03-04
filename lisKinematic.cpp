#include "model.h"

#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
		( ldd != 0 && rFrom >= 0 && cFrom >= 0 &&\
		  rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

#define    D(r,c)  Data[r][c]
#define MAX_ITERS 20

/*
  local drain direction maps have values for directions as following:
    7  8  9
     \ | /
   4 - 5 - 6
     / | \
    1  2  3
*/

//---------------------------------------------------------------------------
static double CalcS2(
    double Qj1i1,  /* Qj+1,i+1 : resultaat kin wave voor deze cell ;j = time, i = place */
    double Qj1i,   /* Qj+1,i   : som van alle bovenstroomse kinematic wave  */
    double Sj1i,   /* Sj+1,i : som van alle bovenstroomse sediment */
    double dt,
    double vol,
    double sedvol)
{
    double Qsn = 0;
    double totsed = sedvol + Sj1i*dt;  // add upstream sed to sed present in cell
    double totwater = vol + Qj1i*dt;   // add upstream water to volume water in cell
    if (totwater <= 1e-10)
       return (Qsn);
    Qsn = min(totsed/dt, Qj1i1 * totsed/totwater);
    return (Qsn); // outflow is new concentration * new out flux

}
//---------------------------------------------------------------------------
static double CalcS1(
    double Qj1i1,  /* Qj+1,i+1 : resultaat kin wave voor deze cell ;j = time, i = place */
    double Qj1i,   /* Qj+1,i   : som van alle bovenstroomse kinematic wave  */
    double Qji1,   /* Qj,i+1 Q voor de kinematic wave (t=j) in deze cell, heet Qin in LISEM */
    double Sj1i,   /* Sj+1,i : som van alle bovenstroomse sediment */
    double Sji1,   /* Si,j+1 : voor de kinematic wave in de cell, heet Qsedin in LISEM */
    double alpha,
    double dt,
    double dx)
{
    double Sj1i1, Cavg, Qavg, aQb, abQb_1, A, B, C, beta = 0.6, s = 0;

    if (Qj1i1 < MIN_FLUX)
       return (0);

    Qavg = 0.5*(Qji1+Qj1i);
    if (Qavg <= MIN_FLUX)
       return (0);

    Cavg = (Sj1i+Sji1)/(Qj1i+Qji1);
    aQb = alpha*pow(Qavg,beta);
    abQb_1 = alpha*beta*pow(Qavg,beta-1);

    A = dt*Sj1i;
    B = -dx*Cavg*abQb_1*(Qj1i1-Qji1);
    C = (Qji1 <= MIN_FLUX ? 0 : dx*aQb*Sji1/Qji1);
    if (Qj1i1 > MIN_FLUX)
       Sj1i1 = (dx*dt*s+A+C+B)/(dt+dx*aQb/Qj1i1);
    else
       Sj1i1 = 0;

    return max(0,Sj1i1);
}
//---------------------------------------------------------------------------
static double IterateToQnew(
    double Qin, /* summed Q new in for all sub-cachments */
    double Qold,  /* current discharge */
    double q,
    double alpha,
    double deltaT,
    double deltaX)
{
  /* Using Newton-Raphson Method */
    double  ab_pQ, deltaTX, C;  //auxillary vars
    int   count;
    double Qkx; //iterated discharge, becomes Qnew
    double fQkx; //function
    double dfQkx;  //derivative
    double _epsilon = 1e-9;
    double beta = 0.6;

    /* if no input then output = 0 */
   if ((Qin + Qold) <= 0)
     return(0);

    /* common terms */
    ab_pQ = alpha*beta*pow(((Qold+Qin)/2),beta-1);
    // derivative of diagonal average (space-time)

    deltaTX = deltaT/deltaX;
    C = deltaTX*Qin + alpha*pow(Qold,beta) + deltaT*q;
    //dt/dx*Q = m3/s*s/m=m2; a*Q^b = A = m2; q*dt = s*m2/s = m2
    //C is unit volume of water
    // first gues Qkx
    Qkx   = (deltaTX * Qin + Qold * ab_pQ + deltaT * q) / (deltaTX + ab_pQ);

  	// VJ 050704, 060830 infil so big all flux is gone
    if (Qkx < 0)
      Qkx   = MIN_FLUX;

    Qkx   = max(Qkx, MIN_FLUX);
    count = 0;
    do {
      fQkx  = deltaTX * Qkx + alpha * pow(Qkx, beta) - C;   /* Current k */
      dfQkx = deltaTX + alpha * beta * pow(Qkx, beta - 1);  /* Current k */
      Qkx   -= fQkx / dfQkx;                                /* Next k */
      Qkx   = max(Qkx, MIN_FLUX);
      count++;
    } while(fabs(fQkx) > _epsilon && count < MAX_ITERS);

    return Qkx;
}
//---------------------------------------------------------------------------

void TWorld::Kinematic(int pitRowNr, int pitColNr,
TMMap *_LDD, TMMap *_Q, TMMap *_Qn, TMMap *_Qs, TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX
,TMMap *_Vol, TMMap*_SedVol
)
{
  int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
  int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

  Liststruct *list = NULL, *temp = NULL;
  list = (Liststruct *)malloc(sizeof(Liststruct));
  list->prev = NULL;
  list->rowNr = pitRowNr;
  list->colNr = pitColNr;

 while (list != NULL)
 {
      int i = 0;
      bool  subCachDone = true; /* are sub-catchment cells done ? */
      int rowNr = list->rowNr;
      int colNr = list->colNr;
   /* put all points that have to be calculated to
      calculate the current point in the list,
	   before the current point */
      for (i=1; i<=9; i++)
      {
        int r, c;
        int ldd = 0;

        if (i==5)  /* this is the current cell*/
            continue;

	r = rowNr+dy[i];
	c = colNr+dx[i];

        if (r>=0 && r<nrRows && c>=0 && c<nrCols &&
           !IS_MV_REAL4(&_LDD->Drc))
             ldd = (int) _LDD->Drc;
           else
             continue;
        if (r>=0 && r<nrRows &&
            c>=0 && c<nrCols &&
            FLOWS_TO(ldd, r, c, rowNr, colNr) &&
            IS_MV_REAL4(&_Qn->Drc) ) /* cell not computed */
        {
            temp = (Liststruct *)malloc(sizeof(Liststruct));
            temp->prev = list;
            list = temp;
            list->rowNr = r;
            list->colNr = c;
            subCachDone = false;
        }
      }

	if (subCachDone)
	{
	  double Qin=0.0, Sin=0.0;
          for (i=1;i<=9;i++) /* for all incoming cells */
          {
            int r, c;
            int ldd = 0;

            if (i==5)  /* Skip current cell */
              continue;

            r = rowNr+dy[i];
            c = colNr+dx[i];

            if (r>=0 && r<nrRows && c>=0 && c<nrCols &&
               !IS_MV_REAL4(&_LDD->Drc))
                 ldd = (int) _LDD->Drc;
               else
                 continue;

            if (r>=0 && r < nrRows &&
                c>=0 && c < nrCols &&
                FLOWS_TO(ldd, r,c,rowNr, colNr) &&
                 !IS_MV_REAL4(&_LDD->Drc) )
              {
                  Qin += _Qn->Drc;
                  if (SwitchErosion)
                    Sin += _Qsn->Drc;
              }
        } /* eof all incoming cells */

	temp=list;
    	list=list->prev;
        free(temp);

      _Qn->D(rowNr,colNr) = IterateToQnew(Qin, _Q->D(rowNr,colNr), _q->D(rowNr,colNr),
                               _Alpha->D(rowNr,colNr), _dt, _DX->D(rowNr,colNr));
       // Newton Rapson iteration for water of current cell
      if (SwitchErosion)
      {
        if (!SwitchSimpleSedKinWave)
          _Qsn->D(rowNr, colNr) = CalcS1(_Qn->D(rowNr,colNr), Qin, _Q->D(rowNr,colNr), Sin, _Qs->D(rowNr,colNr),
                                  _Alpha->D(rowNr,colNr), _dt, _DX->D(rowNr,colNr));
        else
          _Qsn->D(rowNr, colNr) = CalcS2(_Qn->D(rowNr,colNr), Qin, Sin, _dt,
                                  _Vol->D(rowNr,colNr), _SedVol->D(rowNr,colNr));
      }

      _q->D(rowNr,colNr) = Qin;
      //VJ 050831 REPLACE infil with sum of all incoming fluxes, needed for infil calculation
      // q is now in m3/s

      if (SwitchErosion)
      {
        _Qsn->D(rowNr, colNr) = min(_Qsn->D(rowNr, colNr),Sin+_SedVol->D(rowNr,colNr)/_dt);
      // no more sediment outflow than total sed in cell
        _SedVol->D(rowNr,colNr) = max(0, Sin*_dt + _SedVol->D(rowNr,colNr) - _Qsn->D(rowNr, colNr)*_dt);
      // new sed volume based on all fluxes and org sed present
      }

		/* cell rowNr, colNr is now done */
   }/* eof subcatchment done */
 } /* eowhile list != NULL */
}
//---------------------------------------------------------------------------
