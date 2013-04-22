
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
  \file lisKinematic.cpp
  \brief kinematic wave routing functions and calculation of discharge and sed flux per cell.

  The routing functions use local variables because they are used for overland flow, channel flow, gully and tiledrain flow.

functions: \n
   - void TWorld::Kinematic(int pitRowNr, int pitColNr, TMMap *_LDD, TMMap *_Q, TMMap *_Qn, TMMap *_Qs,
                  TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap *SedVol,
                  TMMap *_StorVol, TMMap *_StorVolSed);\n
   - double TWorld::simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double dt, double vol, double sed);\n
   - double TWorld::complexSedCalc(double Qj1i1, double Qj1i, double Qji1, double Sj1i,
                         double Sji1, double alpha, double dt, double dx);\n
   - double TWorld::IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX);\n
 */

#include "model.h"

// check if cell From flows to To
#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

// check if cell is still inside the map boundaries
//moved to model.h
//#define INSIDE(r, c) (r>=0 && r<_nrRows && c>=0 && c<_nrCols)

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
/**
  Simple calculation of sediment outflux from a cell based on the sediment concentration multiplied by the new water flux,
  j = time and i = place: j1i1 is the new output, j1i is the new flux at the upstream 'entrance' flowing into the gridcell
   \param  Qj1i1 Qj+1,i+1 : result kin wave for this cell ;j = time, i = place
   \param  Qj1i  Qj+1,i   : sum of all upstreamwater from kin wave
   \param  Sj1i  Sj+1,i : sum of all upstream sediment
   \param  dt    timestep
   \param  vol   current volume of water in cell
   \param  sed    current mass of sediment in cell
 */
double TWorld::simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double dt, double vol, double sed)
{
    double Qsn = 0;
    double totsed = sed + Sj1i*dt;  // add upstream sed to sed present in cell
    double totwater = vol + Qj1i*dt;   // add upstream water to volume water in cell
    if (totwater <= 1e-10)
        return (Qsn);
    Qsn = min(totsed/dt, Qj1i1 * totsed/totwater);
    return (Qsn); // outflow is new concentration * new out flux

}
//---------------------------------------------------------------------------
/**
  Complex calculation of sediment outflux from a cell based on a explicit solution of the time/space matrix,
  j = time and i = place: j1i1 is the new output, j1i is the new flux at the upstream 'entrance' flowing into the gridcell
   \param Qj1i1   Qj+1,i+1 : result kin wave for this cell ;j = time, i = place
   \param Qj1i    Qj+1,i   : sum of all upstreamwater from kin wave
   \param Qji1    Qj,i+1 : incoming Q for kinematic wave (t=j) in this cell, map Qin in LISEM
   \param Sj1i    Sj+1,i : sum of all upstream sediment
   \param Sji1    Si,j+1 : incoming Sed for kinematic wave (t=j) in this cell, map Qsin in LISEM
   \param alpha   alpha calculated in LISEM from before kinematic wave
   \param dt      timestep
   \param dx      dx: length of the cell, corrected for slope (DX map in LISEM)
 */
double TWorld::complexSedCalc(double Qj1i1, double Qj1i, double Qji1,double Sj1i, double Sji1, double alpha, double dt,double dx)
{
    double Sj1i1, Cavg, Qavg, aQb, abQb_1, A, B, C, s = 0;
    const double beta = 0.6;


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
/** Newton Rapson iteration for new water flux in cell, based on Ven Te Chow 1987
\param Qin      summed Q new from upstream
\param Qold     current discharge in the cell  Qin in LISEM
\param q        infiltration surplus flux (in m2/s), has value <= 0
\param alpha    alpha calculated in LISEM from before kinematic wave
\param deltaT   timestep
\param deltaX   dx: length of the cell  corrected for slope (DX map in LISEM)
*/
double TWorld::IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX)
{
    /* Using Newton-Raphson Method */
    double  ab_pQ, deltaTX, C;  //auxillary vars
    int   count;
    double Qkx; //iterated discharge, becomes Qnew
    double fQkx; //function
    double dfQkx;  //derivative
    const double _epsilon = 1e-12;
    const double beta = 0.6;

    /* if no input then output = 0 */
    if ((Qin + Qold) <= q*deltaX)//0)
    {
        itercount = -1;
        return(0);
    }

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
    //VJ 110114 without this de iteration cannot be solved for very small values
    if (Qkx < MIN_FLUX)
    {
        itercount = -2;
        return(0);
    }

    Qkx   = max(Qkx, MIN_FLUX);

    count = 0;
    do {
        fQkx  = deltaTX * Qkx + alpha * pow(Qkx, beta) - C;   /* Current k */
        dfQkx = deltaTX + alpha * beta * pow(Qkx, beta - 1);  /* Current k */
        Qkx   -= fQkx / dfQkx;                                /* Next k */
        Qkx   = max(Qkx, MIN_FLUX);
        count++;
        //qDebug() << count << fQkx << Qkx;
    } while(fabs(fQkx) > _epsilon && count < MAX_ITERS);
    itercount = count;
    return Qkx;
}
//---------------------------------------------------------------------------
/**
Kinematic wave spatial part, used for slope, channel and tiledrain system: \n
Kinematic is called for each pit (i.e. a cell with value 5 in the LDD):\n
A linked list of all cells connected to the pit is made, after that it 'walks' through the list\n
calculating the fluxes from upstream to downstream.\n
\param   pitRowNr, pitColNr  row and col nr of the current pit (can be more than one, LISEM loops through all the pits, r_outlet and c_oulet is the main pit)
\param  	_LDD                LDD map, can be slope, channel or tiledrain system
\param  	_Q, _Qn             incoming at the start of the timestep and outgoing new discharge (m3/s)
\param  	_Qs, _Qsn	    incoming and outgoing sediment flux (kg/s)
\param  	_q                  infiltration surplus from infiltration functions, 0 or negative value
\param  	_Alpha              a in A=aQ^b
\param  	_DX                 dx corrected for slope
\param  	_Vol                water volume in cell (m3)
\param  	_Sed                sediment load in cell (kg)
\param  	_StorVol            water volume in buffers (m3)
\param  	_StorSed            sediment volume in bufers (kg)
*/
void TWorld::Kinematic(int pitRowNr, int pitColNr, TMMap *_LDD,
                       TMMap *_Q, TMMap *_Qn, TMMap *_Qs, TMMap *_Qsn,
                       TMMap *_q, TMMap *_Alpha, TMMap *_DX,
                       TMMap *_Vol, TMMap*_Sed,
                       TMMap *_StorVol, TMMap *_StorSed)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    /// Linked list of cells in order of LDD flow network, ordered from pit upwards
    LDD_LINKEDLIST *list = NULL, *temp = NULL;
    list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

    list->prev = NULL;
    /// start gridcell: outflow point of area
    list->rowNr = pitRowNr;
    list->colNr = pitColNr;

    if (SwitchErosion)
        _Qsn->fill(0);
    // set output sed flux to 0
    //   _Qn->setMV();
    // flag all Qn gridcell with MV

    while (list != NULL)
    {
        int i = 0;
        bool  subCachDone = true; // are sub-catchment cells done ?
        int rowNr = list->rowNr;
        int colNr = list->colNr;

        /** put all points that have to be calculated to calculate the current point in the list,
         before the current point */
        for (i=1; i<=9; i++)
        {
            int r, c;
            int ldd = 0;

            // this is the current cell
            if (i==5)
                continue;

            r = rowNr+dy[i];
            c = colNr+dx[i];

            if (INSIDE(r, c) && !IS_MV_REAL8(&_LDD->Drc))
                ldd = (int) _LDD->Drc;
            else
                continue;

            // check if there are more cells upstream, if not subCatchDone remains true
            if (IS_MV_REAL4(&_Qn->Drc) &&
                    FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                    INSIDE(r, c))
            {
                temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                temp->prev = list;
                list = temp;
                list->rowNr = r;
                list->colNr = c;
                subCachDone = false;
            }
        }

        // all cells above a cell are linked in a "sub-catchment or branch
        // continue with water and sed calculations
        // rowNr and colNr are the last upstreM cell linked
        if (subCachDone)
        {
            double Qin=0.0, Sin=0.0;

            // for all incoming cells of a cell, sum Q and Sed and put in Qin and Sin
            // note these are values of Qn and Qsn so the new flux
            for (i=1;i<=9;i++)
            {
                int r, c, ldd = 0;

                if (i==5)  // Skip current cell
                    continue;

                r = rowNr+dy[i];
                c = colNr+dx[i];

                if (INSIDE(r, c) && !IS_MV_REAL4(&_LDD->Drc))
                    ldd = (int) _LDD->Drc;
                else
                    continue;

                if (INSIDE(r, c) &&
                        FLOWS_TO(ldd, r,c,rowNr, colNr) &&
                        !IS_MV_REAL4(&_LDD->Drc) )
                {
                    Qin += _Qn->Drc; //_Qn->Drc;

                    if (SwitchErosion)
                        Sin += _Qsn->Drc; //_Qsn->Drc;

                    // ADD MC and NUTs HERE
                }
            }

            bool isBufferCellWater = false; // check if buffer is full, true means not full
            bool isBufferCellSed = false;

            //if buffers, add the incoming water Qin to the buffer
            if(SwitchBuffers)
            {
                //_StorVol is remaining space in buffers, not water in buffers!
                //_StorVol will go to 0
                if (_StorVol->Data[rowNr][colNr] > 0)
                {
                    isBufferCellWater = true;
                    //buffer still active, skip normal kin wave

                    _StorVol->Data[rowNr][colNr] -= Qin*_dt;
                    // fill up storage with incoming water

                    _q->Data[rowNr][colNr] = Qin;
                    Qin = 0;
                    // assume first buffer is not full, no outflow
                    // trick kin wave by saying that the buffer storage is a ink (_q) like infiltration, taken out of the flow
                    // and there is no incoming flux Qin to deal wth in the kin wave

                    // now check if the uffer was overflowing and there is a little bit of Qin!
                    if (_StorVol->Data[rowNr][colNr] < 0)  // negative means overflowing
                    {
                        Qin = -_StorVol->Data[rowNr][colNr]/_dt;
                        //_q->Data[rowNr][colNr] -= Qin;
                        // overflow part becomes flux again
                        _StorVol->Data[rowNr][colNr] = 0;
                        // remaining store = 0
                        isBufferCellWater = false;
                        //buffer is full, go to normal kin wave outflow
                    }

                    if (isBufferCellWater)
                        _Qn->Data[rowNr][colNr] = 0;
                    // if still volume left no outflow to downstream cell
                }
            }

            // sediment in buffers
            if(SwitchErosion && (SwitchBuffers || SwitchSedtrap))
            {
                // if there is water storage catch all the sediment
                if (_StorVol->Data[rowNr][colNr] > 0)
                {
                    isBufferCellSed = true;
                    //buffer still active

                    _StorSed->Data[rowNr][colNr] += Sin*_dt;
                    Sin = 0;

                    if (!SwitchSedtrap)
                    {
                        _StorVol->Data[rowNr][colNr] -= Sin/2600;
                        // decrease storvol with volume loss caused by incoming sediment
                        // the bulkdensity does not matter, the volume taken up is related
                        // to the particle desity dens, because the pores are filled with water
                        // if we use bulk dens here we assume pores are empty!

                        if (_StorVol->Data[rowNr][colNr] < 0)
                        {
                            Qin = -_StorVol->Data[rowNr][colNr]/_dt;
                            // overflow part becomes flux again
                            _StorVol->Data[rowNr][colNr] = 0;
                            // remaining store = 0
                            isBufferCellWater = false;
                            //buffer is full, outflow in kin wave
                            isBufferCellSed = false;
                            //buffer is full, sed outflow in kin wave
                        }
                    }
                }
            }

            // if cell is not a buffer cell or buffer is filled calc outflow with iteration
            if (!isBufferCellWater)
            {
                itercount = 0;
                _Qn->Data[rowNr][colNr] = IterateToQnew(Qin, _Q->Data[rowNr][colNr], _q->Data[rowNr][colNr],
                                                        _Alpha->Data[rowNr][colNr], _dt, _DX->Data[rowNr][colNr]);
                // Newton Rapson iteration for water of current cell

                _q->Data[rowNr][colNr] = Qin;
                //VJ 050831 REPLACE infil with sum of all incoming fluxes, needed for infil calculation below
                // q is now in m3/s
            }

            // if cell is not a buffer cell or buffer is filled calc SED outflow with iteration
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
            /* cell rowN, colNr is now done */

            temp=list;
            list=list->prev;
            free(temp);
            // go to the previous cell in the list

        }/* eof subcatchment done */
    } /* eowhile list != NULL */
}
//---------------------------------------------------------------------------
/**
Kinematic wave spatial part, used for slope, channel and tiledrain system: \n
Kinematic is called for each pit (i.e. a cell with value 5 in the LDD):\n
A linked list of all cells connected to the pit is made, after that it 'walks' through the list\n
calculating the fluxes from upstream to downstream.\n
\param   pitRowNr, pitColNr  row and col nr of the current pit (can be more than one, LISEM loops through all the pits, r_outlet and c_oulet is the main pit)
\param  	_LDD                LDD map, can be slope, channel or tiledrain system
\param  	_Q, _Qn             incoming at the start of the timestep and outgoing new discharge (m3/s)
\param  	_Qs, _Qsn	    incoming and outgoing sediment flux (kg/s)
\param  	_Alpha              a in A=aQ^b
\param  	_DX                 dx corrected for slope
\param  	_Sed                sediment load in cell (kg)
\param  	_StorVol            water volume in buffers (m3)
\param  	_StorSed            sediment volume in bufers (kg)
*/
void TWorld::routeSubstance(int pitRowNr, int pitColNr, TMMap *_LDD,
                            TMMap *_Q, TMMap *_Qn, TMMap *_Qs, TMMap *_Qsn,
                            TMMap *_Alpha, TMMap *_DX, TMMap*_Vol, TMMap*_Sed) //,TMMap *_StorVol, TMMap *_StorSed)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    /// Linked list of cells in order of LDD flow network, ordered from pit upwards
    LDD_LINKEDLIST *list = NULL, *temp = NULL;
    list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
    list->prev = NULL;
    /// start gridcell: outflow point of area
    list->rowNr = pitRowNr;
    list->colNr = pitColNr;

    _Qsn->fill(0);
    // VJ 12/12/12 !!!!!!!!!!!!!!!!! set output substance to zero

    while (list != NULL)
    {
        int i = 0;
        bool  subCachDone = true; // are sub-catchment cells done ?
        int rowNr = list->rowNr;
        int colNr = list->colNr;

        /** put all points that have to be calculated to calculate the current point in the list,
         before the current point */
        for (i=1; i<=9; i++)
        {
            int r, c;
            int ldd = 0;

            // this is the current cell
            if (i==5)
                continue;

            r = rowNr+dy[i];
            c = colNr+dx[i];

            if (INSIDE(r, c) && !IS_MV_REAL8(&_LDD->Drc))
                ldd = (int) _LDD->Drc;
            else
                continue;

            // check if there are more cells upstream, if not subCatchDone remains true
            if (IS_MV_REAL4(&_Qn->Drc) &&
                    FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                    INSIDE(r, c))
            {
                temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                temp->prev = list;
                list = temp;
                list->rowNr = r;
                list->colNr = c;
                subCachDone = false;
            }
        }

        // all cells above a cell are linked in a "sub-catchment or branch
        // continue with water and sed calculations
        // rowNr and colNr are the last upstreM cell linked
        if (subCachDone)
        {
            double Qin=0.0, Sin=0.0;

            // for all incoming cells of a cell, sum Q and Sed and put in Qin and Sin
            // note these are values of Qn and Qsn so the new flux
            for (i=1;i<=9;i++)
            {
                int r, c, ldd = 0;

                if (i==5)  // Skip current cell
                    continue;

                r = rowNr+dy[i];
                c = colNr+dx[i];

                if (INSIDE(r, c) && !IS_MV_REAL4(&_LDD->Drc))
                    ldd = (int) _LDD->Drc;
                else
                    continue;

                if (INSIDE(r, c) &&
                        FLOWS_TO(ldd, r,c,rowNr, colNr) &&
                        !IS_MV_REAL4(&_LDD->Drc) )
                {
                    Qin += _Qn->Drc;
                    Sin += _Qsn->Drc;
                }
            }

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

            /* cell rowN, colNr is now done */

            temp=list;
            list=list->prev;
            free(temp);
            // go to the previous cell in the list

        }/* eof subcatchment done */
    } /* eowhile list != NULL */
}

void TWorld::findFlood(int pitRowNr, int pitColNr, TMMap *_LDD)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    /// Linked list of cells in order of LDD flow network, ordered from pit upwards
    LDD_LINKEDLIST *list = NULL, *temp = NULL;
    list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

    list->prev = NULL;
    /// start gridcell: outflow point of area
    list->rowNr = pitRowNr;
    list->colNr = pitColNr;

    tm->setMV();

    while (list != NULL)
    {
        int i = 0;
        bool  subCachDone = true; // are sub-catchment cells done ?
        int rowNr = list->rowNr;
        int colNr = list->colNr;

        /** put all points that have to be calculated to calculate the current point in the list,
         before the current point */
        for (i=1; i<=9; i++)
        {
            int r, c;
            int ldd = 0;

            // this is the current cell
            if (i==5)
                continue;

            r = rowNr+dy[i];
            c = colNr+dx[i];

            if (INSIDE(r, c) && !IS_MV_REAL8(&_LDD->Drc))
                ldd = (int) _LDD->Drc;
            else
                continue;

            // check if there are more cells upstream, if not subCatchDone remains true
            if (IS_MV_REAL4(&tm->Drc) &&
                    FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                    INSIDE(r, c))
            {
                temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                temp->prev = list;
                list = temp;
                list->rowNr = r;
                list->colNr = c;
                subCachDone = false;
            }
        }

        if (subCachDone)
        {
            double edge = 1;

            for (i=1;i<=9;i++)
            {
                int r, c, ldd = 0;

                if (i==5)  // Skip current cell
                    continue;

                r = rowNr+dy[i];
                c = colNr+dx[i];

                if (INSIDE(r, c) && !IS_MV_REAL4(&_LDD->Drc))
                    ldd = (int) _LDD->Drc;
                else
                    continue;

                if (INSIDE(r, c) &&
                        FLOWS_TO(ldd, r,c,rowNr, colNr) &&
                        !IS_MV_REAL4(&_LDD->Drc) )
                {
                    if (FloodDomain->Drc == 0 &&
                            FloodDomain->Data[rowNr][colNr] > 0)
                        edge += 1;
                }
            }


            if (edge > 0)
                FloodDomain->Data[rowNr][colNr] = edge;
            tm->Data[rowNr][colNr] = 1;

            temp=list;
            list=list->prev;
            free(temp);
            // go to the previous cell in the list

        }/* eof subcatchment done */
    } /* eowhile list != NULL */
    FloodDomain->report("fd");
}
