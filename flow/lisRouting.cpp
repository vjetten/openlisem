/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file lisRouting.cpp
  \brief routing functions and calculation of discharge and sed flux per cell.


functions: \n
 */

#include "model.h"


//---------------------------------------------------------------------------
QVector <LDD_COORIN> TWorld::MakeLinkedList(cTMap *_LDD)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1,  0,  1};
    int dy[10] = {0,  1, 1, 1,  0, 0, 0, -1, -1, -1};

    QVector <LDD_COORIN> _crlinked_;
    _crlinked_.clear();

    Fill(*tma, -1); // flag

    FOR_ROW_COL_MV {
        if (_LDD->Drc == 5) {

            /// Linked list of cells in order of LDD flow network, ordered from pit upwards
            LDD_LINKEDLIST *list = nullptr, *temp = nullptr;
            list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

            list->prev = nullptr;
            /// start gridcell: outflow point of area
            list->rowNr = r;
            list->colNr = c;

            while (list != nullptr)
            {
                int i = 0;
                bool  subCachDone = true;
                int rowNr = list->rowNr;
                int colNr = list->colNr;

                for (i=1; i<=9; i++)
                {

                    // this is the current cell
                    if (i==5)
                        continue;

                    int ldd = 0;
                    int rr = rowNr+dy[i];
                    int cr = colNr+dx[i];

                    if (INSIDE(rr, cr) && !pcr::isMV(_LDD->Drcr))
                        ldd = static_cast<int> (_LDD->Drcr);
                    else
                        continue;

                    // check if there are more cells upstream, if not subCatchDone remains true
                    if (tma->Drcr == -1 && FLOWS_TO(ldd, rr, cr, rowNr, colNr))
                    {
                        temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                        temp->prev = list;
                        list = temp;
                        list->rowNr = rr;
                        list->colNr = cr;
                        subCachDone = false;
                    }
                }

                if (subCachDone)
                {
                    LDD_COORIN newcr;
                    newcr.r = rowNr;
                    newcr.c = colNr;
                    newcr.ldd = static_cast<int> (_LDD->data[rowNr][colNr]);/*
                    if (newcr.ldd < 1 || newcr.ldd > 9) {
                        Error("Invalid ldd found, outside range [1-9]");
                        throw 2;
                    }*/


                    newcr.nr = 0;
                    newcr.inn.clear();

                    int j = 0;
                    for (i=1;i<=9;i++)
                    {
                        if (i != 5) {

                            int rr = rowNr+dy[i];
                            int cr = colNr+dx[i];
                            int ldd = 0;
                            if (INSIDE(rr, cr)) {
                                if (!pcr::isMV(_LDD->Drcr)) {
                                    ldd = static_cast <int>(_LDD->Drcr);
                                    if (FLOWS_TO(ldd, rr,cr,rowNr,colNr))
                                    {
                                       LDD_COOR incr;
                                       incr.r = rr;
                                       incr.c = cr;
                                       newcr.inn << incr; // add the point that flows into the cell to inn
                                       j++;
                                    }
                                }
                            }
                        }
                    }
                    newcr.nr = j;

                    _crlinked_ << newcr;
                    tma->data[rowNr][colNr] = 0;

                    temp=list;
                    list=list->prev;
                    free(temp);
                    // go to the previous cell in the list

                }

            } /* eowhile list != nullptr */

        }
    }

    return(_crlinked_);
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::upstream(cTMap *_LDD, cTMap *_M, cTMap *out)
 * @brief Returns the sum of all values upstream
 *
 * Returns the sum of all values upstream using
 * the local drainage direction map (LDD)
 *
 * @param _M : Material map, can be any substance
 * @param out : Output map, sum of all upstream material
 *
 * @see LDD
 */
void TWorld::upstream(QVector <LDD_COORIN>_crlinked_, cTMap *_Q, cTMap *_Qn)
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _Qn->Drc = 0;
    }}

    for(long i_ =  0; i_ < _crlinked_.size(); i_++)
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;
        double Qin = 0;

        // get inflow
        if (_crlinked_[i_].nr > 0) {
            for(int j = 0; j < _crlinked_[i_].nr; j++) {
                int rr = _crlinked_[i_].inn[j].r;
                int cr = _crlinked_[i_].inn[j].c;
                Qin += _Q->Drcr;
            }
        }
        _Qn->Drc = Qin;
    }
}
//---------------------------------------------------------------------------
void TWorld::downstream(QVector <LDD_COORIN>_crlinked_, cTMap *_Q, cTMap *_Qn)
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _Qn->Drc = 0;
    }}

    for(long i_ =  0; i_ < _crlinked_.size(); i_++)
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;
        int cr, rr;
        double Qin = 0;
        int ldd = _crlinked_[i_].ldd;

        switch (ldd) {
            case 1: rr = r+1; cr = c-1; break;
            case 2: rr = r+1; cr = c  ; break;
            case 3: rr = r+1; cr = c+1; break;
            case 4: rr = r  ; cr = c-1; break;
            case 5: rr = r  ; cr = c  ; break;
            case 6: rr = r  ; cr = c+1; break;
            case 7: rr = r+1; cr = c-1; break;
            case 8: rr = r+1; cr = c  ; break;
            case 9: rr = r+1; cr = c+1; break;
        }
        if (!pcr::isMV(_Q->Drcr))
            _Qn->Drc = _Q->Drcr;
    }
}
//---------------------------------------------------------------------------
void TWorld::upstreamMax(QVector <LDD_COORIN>_crlinked_, cTMap *_MaxQ, cTMap *_Q, cTMap *_Qn)
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _Qn->Drc = 0;
    }}

    for(long i_ =  0; i_ < _crlinked_.size(); i_++)
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;
        double Qin = 0;

        // get inflow
        if (_crlinked_[i_].nr > 0) {
            for(int j = 0; j < _crlinked_[i_].nr; j++) {
                int rr = _crlinked_[i_].inn[j].r;
                int cr = _crlinked_[i_].inn[j].c;
                Qin += _Q->Drcr;
            }
        }
        _Qn->Drc = std::min(_MaxQ->Drc, Qin);
    }
}
//---------------------------------------------------------------------------
void TWorld::UpstreamAvg(QVector <LDD_COORIN>_crlinked_ , cTMap *_Q, cTMap *_Qn)
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _Qn->Drc = 0;
    }}

    for(long i_ =  0; i_ < _crlinked_.size(); i_++)
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;
        double Qin = 0;

        // get inflow
        if (_crlinked_[i_].nr > 0) {
            for(int j = 0; j < _crlinked_[i_].nr; j++) {
                int rr = _crlinked_[i_].inn[j].r;
                int cr = _crlinked_[i_].inn[j].c;
                Qin += _Q->Drcr;
            }
        }
        _Qn->Drc = Qin/ _crlinked_[i_].nr;
    }
}
//---------------------------------------------------------------------------
void TWorld::AccufluxGW(QVector <LDD_COORIN>_crlinked_ , cTMap *_Q, cTMap *_Qn, cTMap *_CW)
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _Qn->Drc = 0;
    }}

    for(long i_ =  0; i_ < _crlinked_.size(); i_++)
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;
        double Qin = 0;

        // get inflow
        if (_crlinked_[i_].nr >0) {
            for(int j = 0; j < _crlinked_[i_].nr; j++) {
                int rr = _crlinked_[i_].inn[j].r;
                int cr = _crlinked_[i_].inn[j].c;
                Qin += (_CW->Drcr > 0 ? 0.0 : _Qn->Drcr);
            }
        }
       _Qn->Drc = Qin + _Q->Drc;
    }
}
