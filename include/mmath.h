/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

//---------------------------------------------------------------------------

#ifndef mmathH
#define mmathH
//---------------------------------------------------------------------------

//#include "mainall.h"
#include "csfmap.h"

#define ADD 0
#define SUB 1
#define MUL 2
#define DIV 3
#define POW 4



class TMMap : public cTMap
{
 //protected:
 public:
   cTMap *Mask;

   void fill(double value);
   void calcV(double v, int oper);
   void calc(cTMap *m, int oper);
   void calc2(cTMap *m1, cTMap *m2, int oper);
   void calc2V(cTMap *m1, double V, int oper);
   void copy(cTMap *m);
   void cover(double v);
   void setMV();
   double MapTotal();

   TMMap();
   ~TMMap();
};


#endif
