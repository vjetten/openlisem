#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#define I(a,b,c) (((a)*(c))+b)

typedef struct ORD{ long i; REAL4 v; }ORD;
static ORD *ord;


int compared(const void *a, const void *b)
{
    long res = ((ORD *)b)->v - ((ORD *)a)->v;
    if (res < 0)
       return -1;
    else
    if (res > 0)
       return +1;
    else
       return 0;
}

void Frequency(MEM_HANDLE *OM, MEM_HANDLE *SM, double frequency, REAL4 subst)
{
	REAL4 **OMap, **SMap;
	size_t nrRows = RgiveNrRows();
	size_t nrCols = RgiveNrCols();
	size_t r, c; /* (r,c) becomes co-ordinate of outflowpoint */
        double cutof_value = 0;
        long i = 0;

        ord=(ORD *)malloc(sizeof(ORD)*nrRows*nrCols);

	OMap=(REAL4 **)CpsNormalHandle(OM,CR_REAL4);
	SMap=(REAL4 **)CpsNormalHandle(SM,CR_REAL4);


	for(r = 0; r < nrRows; r++)
	   SetMemMV(SMap[r],nrCols,CR_REAL4);

        for(r = 0; r < nrRows; r++)
         for(c = 0; c < nrCols; c++)
           if (!IS_MV_REAL4(&OMap[r][c]))
           {
              ord[i].v = OMap[r][c];
              ord[i].i = i;//I(r,c,nrCols);
              i++;
           }
        qsort(ord, i-1, sizeof(ORD), compared);

        cutof_value = ord[(int)frequency*(i-1)].v;

        for(r = 0; r < nrRows; r++)
         for(c = 0; c < nrCols; c++)
           if (!IS_MV_REAL4(&OMap[r][c]))
           {
             if (OMap[r][c] <= cutof_value)
                SMap[r][c] = OMap[r][c];
             else
                SMap[r][c] = subst;
           }
        free(ord);
}
