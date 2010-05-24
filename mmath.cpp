/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/
/*
 * mmath provides a few convenient basic math function for maps
 * fill with a value, cover, add two maps etc
 */
#include "model.h"


//---------------------------------------------------------------------------
__fastcall TMMap::TMMap()
    : cTMap()
{

}
//---------------------------------------------------------------------------
__fastcall TMMap::~TMMap()
{

}
//---------------------------------------------------------------------------

void TMMap::fill(double value)
{
    int r, c;

     for (r = 0; r < nrRows; r++)
      for (c = 0; c < nrCols; c++)
      if (!IS_MV_REAL4(&Data[r][c]))
      {
          Data[r][c] = value;
      }
}
//---------------------------------------------------------------------------
void TMMap::cover(double value)
{
    int r, c;

     for (r = 0; r < nrRows; r++)
      for (c = 0; c < nrCols; c++)
      if (IS_MV_REAL4(&Data[r][c]))
      {
          Data[r][c] = value;
      }
 //     else
   //       SET_MV_REAL4(&Data[r][c]);
}
//---------------------------------------------------------------------------
void TMMap::calc(cTMap *M, int oper)
{
  for (int r = 0; r < M->nrRows; r++)
   for (int c = 0; c < M->nrCols; c++)
  if (!IS_MV_REAL4(&Data[r][c]))
  {
      if (!IS_MV_REAL4(&M->Data[r][c]))
      {
          switch (oper)
          {
          case ADD: Data[r][c] += M->Data[r][c]; break;
          case SUB: Data[r][c] -= M->Data[r][c]; break;
          case MUL: Data[r][c] *= M->Data[r][c]; break;
          case DIV: if (M->Data[r][c] > 0) Data[r][c] /= M->Data[r][c];
                       else SET_MV_REAL4(&Data[r][c]); break;
          case POW: Data[r][c] = powl(Data[r][c],M->Data[r][c]); break;
          }
      }
      else
          SET_MV_REAL4(&Data[r][c]);
  }
}
//---------------------------------------------------------------------------
void TMMap::calc2(cTMap *M1, cTMap *M2, int oper)
{
   for (int r = 0; r < nrRows; r++)
    for (int c = 0; c < nrCols; c++)
    if (!IS_MV_REAL4(&Data[r][c]))
    {
      if (!IS_MV_REAL4(&M1->Data[r][c]) && !IS_MV_REAL4(&M2->Data[r][c]))
      {
          switch (oper)
          {
          case ADD: Data[r][c] = M1->Data[r][c] + M2->Data[r][c]; break;
          case SUB: Data[r][c] = M1->Data[r][c] - M2->Data[r][c]; break;
          case MUL: Data[r][c] = M1->Data[r][c] * M2->Data[r][c]; break;
          case DIV: if (M2->Data[r][c] > 0) Data[r][c] = M1->Data[r][c] / M2->Data[r][c];
                       else SET_MV_REAL4(&Data[r][c]); break;
          case POW: Data[r][c] = pow(M1->Data[r][c],M2->Data[r][c]); break;
          }
      }
      else
          SET_MV_REAL4(&Data[r][c]);
    }
}
//---------------------------------------------------------------------------
void TMMap::calcV(double V, int oper)
{
     for (int r = 0; r < nrRows; r++)
      for (int c = 0; c < nrCols; c++)
      if (!IS_MV_REAL4(&Data[r][c]))
      {
          switch (oper)
          {
          case ADD: Data[r][c] += V; break;
          case SUB: Data[r][c] -= V; break;
          case MUL: Data[r][c] *= V; break;
          case DIV: if (V > 0) Data[r][c] /= V;
                       else SET_MV_REAL4(&Data[r][c]); break;
          case POW: Data[r][c] = pow(Data[r][c],V); break;
          }
      }
}
//---------------------------------------------------------------------------
double TMMap::MapTotal()
{
     double total = 0;
     for (int r = 0; r < nrRows; r++)
      for (int c = 0; c < nrCols; c++)
      if (!IS_MV_REAL4(&Data[r][c]))
      {
          total = total + Data[r][c];
      }
      return (total);

}
//---------------------------------------------------------------------------
void TMMap::copy(cTMap *M)
{
     for (int r = 0; r < M->nrRows; r++)
      for (int c = 0; c < M->nrCols; c++)
      if (!IS_MV_REAL4(&M->Data[r][c])&& !IS_MV_REAL4(&Data[r][c]))
      {
          Data[r][c] = M->Data[r][c];
      }
      else
          SET_MV_REAL4(&Data[r][c]);

}
//---------------------------------------------------------------------------
void TMMap::setMV()
{
   for(int r = 0; r < nrRows; r++)
	   SetMemMV(Data[r],nrCols,CR_REAL4);
}
//---------------------------------------------------------------------------

