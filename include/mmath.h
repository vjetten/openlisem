/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
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

/// #include "CsfMap.h"

/// #define ADD 0
/// #define SUB 1
/// #define MUL 2
/// #define DIV 3
/// #define POW 4
/// #define MIN 5  //VJ 041120 added this functionality
/// #define MAX 6
/// #define LARGER 7
/// #define SMALLER 8
/// #define LARGEREQUAL 9
/// #define SMALLEREQUAL 10
/// #define HIGHER 11
/// #define LOWER 12

/// #define CTMap cTMap

/*
/// class defining some basic map algebra operations
class CTMap : public cTMap
{
 //protected:
 public:
   //cTMap *Mask;

   void fill(double value);
   void calcValue(double v, int oper);
   void calcMap(cTMap *m, int oper);
   void calc2Maps(cTMap *m1, cTMap *m2, int oper);
   void calcMapValue(cTMap *m1, double V, int oper);
   void copy(cTMap *m);
 //  void shift(cTMap *M, int dr, int dc);
   void cover(cTMap *m, double v);
   void setMV();
   void checkMap(int oper, double V, QString SS);
   int countUnits();
   double mapTotal();
   double mapAverage();
   double mapMinimum();
   double mapMaximum();
 //  void areaAverage(CTMap *area);

   CTMap();
   ~CTMap();
};

*/
#endif
