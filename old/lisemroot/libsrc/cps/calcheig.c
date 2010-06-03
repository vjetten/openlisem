#include "stddefx.h"

#include "cps.h"

#include <math.h> /*pow*/

static double IterateToNew_h(
	double Q,
	double b, /* Bw in Applied Hydrology */
	double h, /* previous h; y in Applied Hydrology */
	double z,
	double n,
	double s,
	double epsilon)
/* See Applied Hydrology p. 160 and further             */
/* modified by Victor Jetten -> check iteration on dk,  */
/* orig. checked on (Q-Qapprox), also h estimated here  */
/* and not in LISEM                                     */
{
	double Qapprox;

	/* common terms: */
	double _2div3 = 2.0/(double)3.0;
	double sqrtSn  = sqrt(s)/n;
	double sqrt1Z = sqrt(1+z*z);
    double zh, dk;
    int count = 0;


    /* first estimate of h as rectangular channel */
    if (b > 0)
       h = pow((Q*n)/(sqrt(s)*b),0.6);

    if (h < epsilon)
       return(0);

  	do {

        zh = b*h+z*h*h;
	    Qapprox = zh * pow(zh/(b+2*h*sqrt1Z),_2div3) * sqrtSn;

        dk = (1-(Q/Qapprox)) /
             (
	           ( (b+2*z*h)*(5*b+6*h*sqrt1Z)+(4*z*h*h*sqrt1Z) ) /
	           (  3*zh*(b+2*h*sqrt1Z) )
	         );

        h -= dk;

        count++;

	} while(fabs(dk) > epsilon*h && count < 20);

/*
    if (epsilon > 1e-13)
    {
     printf("ITER h %e dk %e \n",h,dk);
    }
*/

    if (h < epsilon)
       return(0);

	return(h);
}


void CalcHeight(MEM_HANDLE*Q,MEM_HANDLE*b, MEM_HANDLE*h, MEM_HANDLE*z, 
                MEM_HANDLE*n, MEM_HANDLE*s, double epsilon)
/*P Q r- is the Q wherefor h has to be appoximated*/
/*P b r- is the cell width,Bw in Applied Hydrology */
/*P h rw is the first approximated h y in Applied Hydrology */
/*P z rw  */
/*P n r- is manning's n*/
/*P s r- is the slope*/
/*P epsilon r- is the difference allowed between the appoxed Q and real Q*/
/* See Applied Hydrology p. 160 and further */
{
	size_t i,nrCells = GiveNrDefinedCells();
        REAL4 *Qmap=(REAL4 *)CpsCompressedHandle(Q,CR_REAL4);
        REAL4 *bmap=(REAL4 *)CpsCompressedHandle(b,CR_REAL4);
        REAL4 *hmap=(REAL4 *)CpsCompressedHandle(h,CR_REAL4);
        REAL4 *zmap=(REAL4 *)CpsCompressedHandle(z,CR_REAL4);
        REAL4 *nmap=(REAL4 *)CpsCompressedHandle(n,CR_REAL4);
        REAL4 *smap=(REAL4 *)CpsCompressedHandle(s,CR_REAL4);


	for (i=0;i<nrCells;nrCells++)
	  	if (!IS_MV_REAL4(Qmap+i))
		   hmap[i]=(REAL4)IterateToNew_h(
							Qmap[i],
							bmap[i],
							hmap[i],
							zmap[i],
							nmap[i],
							smap[i],
							epsilon);

}
