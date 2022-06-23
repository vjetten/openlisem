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
  \file lisPesticide_MC.cpp
  \brief Transport and partitioning of pesticides with Euler forward method

functions: \n
- void TWorld::SimplePestConc() \n
-
*/

#include "model.h"
#include "operation.h"
#include <tuple>

// check if cell From flows to To
//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )


//---------------------------------------------------------------------------
/**
 * @fn double TWorld::simplePestCalc(double Crw_old, double Cmw_old, double Kfilm, double Qinf, double zm, double kr, double Kd,
                              double Crs_old, double Cms_old, double Ez, double Me, double A)
 * @brief Simple calculation of pesticide concentrations in a cell in both water and sediment
 * Simple calculation of pesticide concentrations in a cell in both water and sediment,
//update part below
 * @param  sed : current mass of sediment in cell
 * @return five concentrations for pesticides in sediment, water of runoff and mixing zone and infiltration
 *
 */
void TWorld::simplePestConc(double Crw_old, double Cmw_old, double Kfilm, double Qinf, double zm, double kr, double Kd,
                              double Crs_old, double Cms_old, double Ez, double Me, double A, double pore,
                            double Crw_n, double Crs_n, double Cmw_n, double Cms_n, double Cinf_n) // MC - this are the target values to be returned
{
    Crw_n = 0;
    Cinf_n = 0;
    Cmw_n = 0;
    Cms_n = 0;
    Crs_n = 0;
    double rho = 2650;

    Crw_n = Crw_old + _dt * (Kfilm * (Cmw_old - Crw_old));
    Cmw_n = Cmw_old + _dt * ((Kfilm * (Crw_old - Cmw_old) + Qinf * (Crw_old - Cmw_old) - zm * rho * kr * (Kd * Cmw_old - Cms_old))/ pore * zm);

    if (Ez > 0 ) {
        Cms_n = Cms_old + _dt * (kr * (Kd * Cmw_old - Cms_old));
        Crs_n = Crs_old + _dt * (((Crs_old * Me) + (Cms_old * Ez * rho * A))/(Me + (Ez * rho * A)));
    } else {
        Cms_n = Cms_old + _dt * (kr * (Kd * Cmw_old - Cms_old) + (((zm + Ez * Cms_old) - Crs_old * Ez)/ zm));
        Crs_n = Crs_old;
    }

    Cinf_n = 0.5 * (Crw_old + Crw_n);
}

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::MassPest(double WaterVolall, double Sed, double Qsn, double Qn, double Qinf)
 * @brief Calculation of pesticide mass at the end of each timestep.
//update part below
 * Simple calculation of pesticide concentrations in a cell in both water and sediment,
 * @param  sed : current mass of sediment in cell
 * @return sediment outflow in next timestep
 *
 */
double TWorld::MassPest(double WaterVolall, double Sed, double Qsn, double Qn, double Qinf, double PMtotI)
{
    double PMtot = 0;
    double PMerr = 0;

    PMtot = Pestinf + PestOutS + PestOutW + PMsoil + PMrw + PMrs + PMmw + PMms;
    PMerr = 1 - (PMtot/PMtotI);


}

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::MassPestInitial(double dx, cTMap *Cms, cTMap *Cmw, cTMap *zm, cTMap *zs, cTMap *ThetaI1)
 * @brief Calculation total pesticide mass initial in system.
 * @param dx: cell resolution [m]
 * @param Cms: map with initial pesticide concentration of soil in mixing zone
 * @param Cmw: map with initial pesticide concentration of water in mixing zone
 * @param zm: map with thickness of mixing layer [m]
 * @param zm: map with thickness of soil layer with pesticides [m]
 * @param ThetaI1: map with initial soil moisture of the first soil layer.
 * @return PMtotI
 *
 */
double TWorld::MassPestInitial(double dx, cTMap *PCms, cTMap *PCmw, cTMap *zm, cTMap *zs, cTMap *ThetaI1)
{
    double PMtotI = 0;
    double rho = 2650;
    // PMtotI = PMsoil + PMmw + PMms
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        PMms->Drc = PCms->Drc * dx * dx * rho * zm->Drc;
        PMmw->Drc = PCmw->Drc * ThetaI1->Drc * zm->Drc * dx * dx;
    }}

    PMtotI = MapTotal(*PMmw) + MapTotal(*PMms);
    return(PMtotI);
}


//---------------------------------------------------------------------------
/**
* @fn double TWorld::MassPestInitial(double dx, cTMap *Cms, cTMap *Cmw, cTMap *zm, cTMap *zs, cTMap *ThetaI1)
* @brief Calculation total pesticide mass initial in system.
* @param dx: cell resolution [m]
* @param Cms: map with initial pesticide concentration of soil in mixing zone
* @param Cmw: map with initial pesticide concentration of water in mixing zone
* @param zm: map with thickness of mixing layer [m]
* @param zm: map with thickness of soil layer with pesticides [m]
* @param ThetaI1: map with initial soil moisture of the first soil layer.
* @return PMtotI
*
*/
/*LDD_COOR *_crlinked_*/
void TWorld::KinematicPestMC(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD, cTMap *_Q, cTMap *_Qn, cTMap *_Qs, cTMap *_Qsn,
                             cTMap *_Alpha,cTMap *_DX, cTMap *_Sed, cTMap *_Qpn, cTMap *_Qpsn)
{
   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        SpinKW->Drc = 0;
      //  QpinKW->Drc = 0;
    }}

//#pragma omp parallel for reduction(+:Qin) num_threads(userCores)
    for(long i_ =  0; i_ < _crlinked_.size(); i_++) //_crlinked_.size()
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;

        double Qin = 0;
        double Sin = 0;
        double Qpin = 0;
        double Spin = 0;

        for (int i = 1; i <= 9; i++)
        {
            if (i != 5) {
                int ldd = 0;
                int rr = r+dy[i];
                int cr = c+dx[i];

                if (INSIDE(rr, cr) && !pcr::isMV(_LDD->Drcr)) {
                    ldd = (int) _LDD->Drcr;
                    // if the cells flow into
                    if (FLOWS_TO(ldd, rr,cr,r,c)) {
                        Qin += _Qn->Drcr;
                        Sin += _Qsn->Drcr;
                        Qpin += _Qpn->Drcr;
                        Spin += _Qpsn->Drcr;
                    }
                }
            }
        }


        SpinKW->Drc = Spin;
        QpinKW->Drc = Qpin;
// update the fomulas below (in simplePestCalc) to include the pesticide influx in the cells!!
        simplePestConc(Crw_old, Cmw_old, Kfilm, Qinf, zm, kr, Kd, Crs_old, Cms_old, Ez, Me, A, pore,
                     Crw_n, Crs_n, Cmw_n, Cms_n, Cinf_n);

        // The four new concentrations
        PCrw->Drc = Crw_n;
        PCrs->Drc = Crs_n;
        PCmw->Drc = Cmw_n;
        PCms->Drc = Cms_n;

        // The new fluxes
        // no more outflow than total pesticide in domain
        PQrs->Drc = std::min(PCrs->Drc * Qsn->Drc, SpinKW->Drc + PMrw->Drc/_dt);
        // for water first substract infiltration than runoff - TODO!!
        PQinf->Drc = Cinf_n * Qinf;
        PQrw->Drc = PCrw->Drc * Qn->Drc;

        // The new masses


        _Qsn->Drc = complexSedCalc(_Qn->Drc, Qin, _Q->Drc, Sin, _Qs->Drc, _Alpha->Drc, _DX->Drc);
        _Qsn->Drc = std::min(_Qsn->Drc, SinKW->Drc+_Sed->Drc/_dt);
        // no more sediment outflow than total sed in cell

        _Sed->Drc = std::max(0.0, SinKW->Drc*_dt + _Sed->Drc - _Qsn->Drc*_dt);
        // new sed volume based on all fluxes and org sed present

    }
}
