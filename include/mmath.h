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
  \file mmath.h
  \brief basic map mathematics on two maps or maps and variables
  */

#ifndef mmathH
#define mmathH

#include "csfmap.h"

#define ADD 0
#define SUB 1
#define MUL 2
#define DIV 3
#define POW 4
#define MIN 5  //VJ 041120 added this functionality
#define MAX 6
#define LARGER 7
#define SMALLER 8
#define LARGEREQUAL 9
#define SMALLEREQUAL 10


/// class defining some basic map algebra operations
class TMMap : public cTMap
{
 //protected:
 public:
   cTMap *Mask;

   void fill(double value);
   void calcValue(double v, int oper);
   void calcMap(cTMap *m, int oper);
   void calc2Maps(cTMap *m1, cTMap *m2, int oper);
   void calcMapValue(cTMap *m1, double V, int oper);
   void copy(cTMap *m);
   void cover(cTMap *m, double v);
   void setMV();
   void checkMap(int oper, double V, QString SS);
   double mapTotal();
   double mapAverage();
   double mapMinimum();
   double mapMaximum();
   void areaAverage(TMMap *area);

   TMMap();
   ~TMMap();
};


#endif
