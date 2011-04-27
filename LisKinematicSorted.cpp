
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
  \file lisKinematicSorted.cpp
  \brief A new kinematic wave routing based on a previously sorted LDD, might be faster.

functions: \n
   - void TWorld::KinematicSorted(LDD_POINT **_lddlist, long _lddlistnr,
                     TMMap *_Q, TMMap *_Qn, TMMap *_Qs,
                     TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap *SedVol,
                     TMMap *_StorVol, TMMap *_StorVolSed);\n
   - LDD_POINT** TWorld::makeSortedNetwork(TMMap *_LDD, long *lddlistnr); \n

 */

#include "model.h"

// check if cell From flows to To
#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
   ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

// check if cell is still inside the map boundaries
#define INSIDE(r, c) (r>=0 && r<_nrRows && c>=0 && c<_nrCols)

#define MAX_ITERS 10

/*
  local drain direction maps have values for directions as follows:
    7  8  9
     \ | /
   4 - 5 - 6
     / | \
    1  2  3
 */

//---------------------------------------------------------------------------
void TWorld::KinematicSorted(LDD_POINT **_lddlist, long _lddlistnr,
                        TMMap *_Q, TMMap *_Qn, TMMap *_Qs, TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX
                        ,TMMap *_Vol, TMMap*_Sed, TMMap *_StorVol, TMMap *_StorSed)
{
   for (long count = 0; count < _lddlistnr; count++)
   {
      double Qin=0.0, Sin=0.0;
      int rowNr = _lddlist[count][0].rowNr;
      int colNr = _lddlist[count][0].colNr;
      int i = 1;

      //for (i=1;i<=9;i++) // for all incoming cells of this cell
      while (_lddlist[count][i].rowNr > 0)
      {
         int r = _lddlist[count][i].rowNr;
         int c = _lddlist[count][i].colNr;
       //  if (r >= 0)
       //  {
            Qin += _Qn->Drc;
            if (SwitchErosion)
               Sin += _Qsn->Drc;
         //}
            i++;
      } // sum all incoming cells

      bool isBufferCellWater = false;
      bool isBufferCellSed = false;
      //double incoming = 0;

      //water in buffers
      if(SwitchBuffers)
      {
         //_StorVol is remaining space in buffers, not water in buffers.
         //_StorVol will go to 0
         if (BufferID->Data[rowNr][colNr] > 0 && _StorVol->Data[rowNr][colNr] > 0)
         {
            isBufferCellWater = true;
            //buffer still active
            _StorVol->Data[rowNr][colNr] -= Qin*_dt;
            // fill up storage with incoming water
            Qin = 0;
            // buffer is not full, no outflow
            if (_StorVol->Data[rowNr][colNr] < 0)  // store overflowing
            {
               Qin = -_StorVol->Data[rowNr][colNr]/_dt;
               // overflow part becomes flux again
               _StorVol->Data[rowNr][colNr] = 0;
               // remaining store = 0
               isBufferCellWater = false;
               //buffer is full, outflow
            }

            if (isBufferCellWater)
               _Qn->Data[rowNr][colNr] = 0;
         }
      }

      if(SwitchBuffers || SwitchSedtrap)
      {
         if (BufferID->Data[rowNr][colNr] > 0 && _StorSed->Data[rowNr][colNr] > 0)
         {
            isBufferCellSed = true;
            //buffer still active
            _StorSed->Data[rowNr][colNr] -= Sin*_dt;
            // add incoming to sed store, note: sed store calculated in datainit
            if (!SwitchSedtrap)
            {
               // TODO check this
               //incoming = Sin*_dt/2650;
               // fill water store up with sediment, decreasing volume
               // the bulkdensity does not matter, the volume taken up is related
               // to the particle desity dens, because the pores are filled
               // if we use bulk dens here we assume pores are empty!
               //		_StorVol->Data[rowNr][colNr] -= incoming;
               //	_StorVol->Data[rowNr][colNr] = max(0, _StorVol->Data[rowNr][colNr]);
               //	if (BufferVolInit->Data[rowNr][colNr] > 0)
               //		BufferVolInit->Data[rowNr][colNr] -= incoming;
               //	if (ChannelBufferVolInit->Data[rowNr][colNr] > 0)
               //	ChannelBufferVolInit->Data[rowNr][colNr] -= incoming;
               //adjust the total volume because it has decreased,
               //Note: the extra released water is not made avaliable
               // channel and slope are mutually exclusive, one or the other
            }
            Sin = 0;
            if (_StorSed->Data[rowNr][colNr] < 0)
            {
               Sin = -_StorSed->Data[rowNr][colNr]/_dt;
               _StorSed->Data[rowNr][colNr] = 0;
               isBufferCellSed = false;
               //buffer is full, outflow
            }

            if (isBufferCellSed)
            {
               _Qsn->Data[rowNr][colNr] = 0;
               _Sed->Data[rowNr][colNr] = max(0, Sin*_dt + _Sed->Data[rowNr][colNr] - _Qsn->Data[rowNr][colNr]*_dt);
            }

         }
      }

      if (!isBufferCellWater)
      {
         itercount = 0;
         _Qn->Data[rowNr][colNr] = IterateToQnew(Qin, _Q->Data[rowNr][colNr], _q->Data[rowNr][colNr],
                                             _Alpha->Data[rowNr][colNr], _dt, _DX->Data[rowNr][colNr]);
         // Newton Rapson iteration for water of current cell
         _q->Data[rowNr][colNr] = Qin;
         //VJ 050831 REPLACE infil with sum of all incoming fluxes, needed for infil calculation
         // q is now in m3/s
         tm->Data[rowNr][colNr] = itercount;
      }

      if (SwitchErosion && !isBufferCellSed)
      {
         if (!SwitchSimpleSedKinWave)
            _Qsn->Data[rowNr][colNr] = complexSedCalc(_Qn->Data[rowNr][colNr], Qin, _Q->Data[rowNr][colNr],
                                                   Sin, _Qs->Data[rowNr][colNr], _Alpha->Data[rowNr][colNr], _dt, _DX->Data[rowNr][colNr]);
         else
            _Qsn->Data[rowNr][colNr] = simpleSedCalc(_Qn->Data[rowNr][colNr], Qin, Sin, _dt,
                                                  _Vol->Data[rowNr][colNr], _Sed->Data[rowNr][colNr]);

         _Qsn->Data[rowNr][colNr] = min(_Qsn->Data[rowNr][colNr], Sin+_Sed->Data[rowNr][colNr]/_dt);
         // no more sediment outflow than total sed in cell
         _Sed->Data[rowNr][colNr] = max(0, Sin*_dt + _Sed->Data[rowNr][colNr] - _Qsn->Data[rowNr][colNr]*_dt);
         // new sed volume based on all fluxes and org sed present

      }
   }
  // tm->report("iter");
   // for debug analysis
}
//---------------------------------------------------------------------------
// this function makes an array of nrCells x 10 places for the LDD_POINT structure
LDD_POINT** TWorld::makeSortedNetwork(TMMap *_LDD, long *_lddlistnr)
{
   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
   LDD_POINT **_lddlist;
   long nr = 0;
   long count = 0;

//   LDD_POINT *lddp = new LDD_POINT;
//   listldd << lddp;

   tm->fill(-1);
   // flag none of the cells are processed

   for (int row = 0; row < _nrRows; row++)
      for (int col = 0; col < _nrCols; col++)
         if(!IS_MV_REAL8(&_LDD->Data[row][col]))
         {
            nr++; // count nr of non MV points in _LDD
         }

   *_lddlistnr = nr;
   _lddlist = new LDD_POINT*[nr];
   for(long _r=0; _r < nr; _r++)
   {
      _lddlist[_r] = new LDD_POINT[10];

      for (int i = 0; i < 10; i++)
      {
         _lddlist[_r][i].rowNr = -1;
         _lddlist[_r][i].colNr = -1;
      }
   }
   // make list structure and init to -1

   for (int row = 0; row < _nrRows; row++)
      for (int col = 0; col < _nrCols; col++)
         if(!IS_MV_REAL8(&_LDD->Data[row][col]))
         if (_LDD->Data[row][col] == 5) // is a pit
   {
      LDD_LINKEDLIST *list = NULL, *temp = NULL;

      list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
      list->prev = NULL;
      list->rowNr = row;
      list->colNr = col;
      // start the linked list with a pit

      // find the list from he outlet onward
      while (list != NULL)// && count < nr)
      {
         bool  subCachDone = true;
         int i = 0;
         int rowNr = list->rowNr;
         int colNr = list->colNr;

         // find all cells flowing to current point
         for (i=1; i<=9; i++)
         {
            int r, c;

            if (i==5)  /* this is the current cell*/
               continue;

            r = rowNr+dy[i];
            c = colNr+dx[i];

            // a cell that flows to the current cell is added to the linked list
            // untill no cells flow to the current cell (top of a branch)
            // then subcatchmentdone remains true
            if (INSIDE(r, c) && !IS_MV_REAL8(&_LDD->Drc) &&
                tm->Drc < 0 &&
                FLOWS_TO((int) _LDD->Drc, r, c, rowNr, colNr)
                )
            {
               temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
               temp->prev = list;
               list = temp;
               list->rowNr = r;
               list->colNr = c;
               subCachDone = false;
            }
         }

         // if no more upstream in this branch, walk the linked list back down
         if (subCachDone)
         {
            int j = 1;
            // put current cell in the list on column 0
            _lddlist[count][0].rowNr = list->rowNr;
            _lddlist[count][0].colNr = list->colNr;

            // look around this cell and save the r and c of inflowing cells
            for (i = 1; i <= 9; i++)
            {
               int r, c;

               if (i==5)  // Skip current cell itself
                  continue;

               r = rowNr+dy[i];
               c = colNr+dx[i];

               // add all cells flowing to current cell (between 0 and 8)
               // on columns 1 to 9 (5 is not used in fact but makes it easier)
               if (INSIDE(r, c) && !IS_MV_REAL8(&_LDD->Drc) &&
                   FLOWS_TO((int) _LDD->Drc, r, c, rowNr, colNr)
                   )
               {
                  _lddlist[count][j].rowNr = r;
                  _lddlist[count][j].colNr = c;
                  j++;
               }
            }

            tm->Data[list->rowNr][list->colNr] = 1;
            // flag this cell is done

            count++;

            temp=list;
            list=list->prev;
            free(temp);
            // go one back in the list and free the memory
            // it is made again for the nest catchment

         }
      }
   }
   return(_lddlist);
}
//---------------------------------------------------------------------------

