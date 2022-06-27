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
- double TWorld::MassPest(double WaterVolall, double Sed, double Qsn, double Qn, double Qinf, double PMtotI) \n
- double TWorld::MassPestInitial(cTMap *PCms, cTMap *PCmw, cTMap *zm, cTMap *zs, cTMap *ThetaI1) \n
- void TWorld::KinematicPestMC() \n
-
*/

#include "model.h"
#include "operation.h"

// check if cell From flows to To
//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::simplePestConc(double Crw_old, double Cmw_old, double Kfilm, double Qinf, double zm, double kr, double Kd,
                            double Crs_old, double Cms_old, double Ez, double Me, double A, double pore,
                            double Crw_in, double Crs_in,
                            double &Crw_n, double &Crs_n, double &Cmw_n, double &Cms_n, double &Cinf_n)
 * @brief Simple calculation of pesticide concentrations in a cell in both water and sediment
 * @return five concentrations for pesticides in sediment, water of runoff and mixing zone and infiltration
 *
 *
 * Layout of space (i) and time (j) towards next value
 *
 * Cj1i > > >  Cj1i1
 *  \            ^
 *  \            ^  ^time^
 *  \            ^
 * Cji -------- Cji1
 *     <-space->
 *
 * Cji = the concentration 1 cell upstream previous timestep
 * Cji1 = the concentration in current cell previous timestep
 * Cj1i = the concentration 1 cell upstream this timestep, we know it already because we calculate from
 *        upstream to downstream
 * Cj1i1 = the concentration in current cell this timestep - we want to calculate this.
 *
 * For flux calculations we use the average concentration of Cj1i and Cji1 to go to Cj1i1.
 */
void TWorld::simplePestConc(double Crw_old, double Cmw_old, double Kfilm, double Qinf, double zm, double kr, double Kd,
                            double Crs_old, double Cms_old, double Ez, double Me, double A, double pore,
                            double Crw_in, double Crs_in,
                            double &Crw_n, double &Crs_n, double &Cmw_n, double &Cms_n, double &Cinf_n) // MC - this are the target values to be returned
{
    Crw_n = 0;
    Cinf_n = 0;
    Cmw_n = 0;
    Cms_n = 0;
    Crs_n = 0;
    double rho = 2650;
    // calculate average between C of fluxes
    double Crw_avg = 0.5 * (Crw_old + Crw_in);
    double Crs_avg = 0.5 * (Crs_old + Crs_in);

    Crw_n = Crw_old + _dt * (Kfilm * (Cmw_old - Crw_avg));
    Cmw_n = Cmw_old + _dt * ((Kfilm * (Crw_avg - Cmw_old) + Qinf * (Crw_avg - Cmw_old) - zm * rho * kr * (Kd * Cmw_old - Cms_old))/ pore * zm);

    if (Ez > 0 ) {
        //erosion
        Cms_n = Cms_old + _dt * (kr * (Kd * Cmw_old - Cms_old));
        Crs_n = Crs_old + _dt * (((Crs_avg * Me) + (Cms_old * Ez * rho * A))/(Me + (Ez * rho * A)));
    } else {
        //deposition
        Cms_n = Cms_old + _dt * (kr * (Kd * Cmw_old - Cms_old) + (((zm + Ez * Cms_old) - Crs_old * Ez)/ zm));
        Crs_n = Crs_avg; // or Crs_old??
    }

    Cinf_n = 0.5 * (Crw_old + Crw_n); //Or Crw_avg??
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
void TWorld::MassPest(cTMap *PMrw, cTMap *PMrs, cTMap *PMms, cTMap *PMmw, cTMap *PMsoil, double PMtotI, double &PMerr, double &PMtot)
{

    // totals of outfluxes
    // at the moment only overland flow is accounted for. is channels etc must be included build that later.
    Pestinf += MapTotal(*PQinf) * _dt;
    FOR_ROW_COL_LDD5 {
        PQrw_dt += PQrw->Drc * _dt;
        PQrs_dt += PQrs->Drc * _dt;
    }}

    PestOutW += PQrw_dt;
    PestOutS += PQrs_dt;

    PMtot = Pestinf + PestOutS + PestOutW + MapTotal(*PMsoil) + MapTotal(*PMrw) + MapTotal(*PMrs) + MapTotal(*PMmw) + MapTotal(*PMms);

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
 * @param zs: map with thickness of soil layer with pesticides [m]
 * @param ThetaI1: map with initial soil moisture of the first soil layer.
 * @return PMtotI
 *
 */
double TWorld::MassPestInitial(cTMap *PCms, cTMap *PCmw, cTMap *zm, cTMap *zs, cTMap *ThetaI1)
{
    double PMtotI = 0;
    double rho = 2650;
    // PMtotI = PMsoil + PMmw + PMms
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        PMms->Drc = PCms->Drc * _dx * _dx * rho * zm->Drc;
        PMmw->Drc = PCmw->Drc * ThetaI1->Drc * zm->Drc * _dx * _dx;
        PMsoil->Drc = PCms->Drc * _dx * _dx * zs->Drc;
    }}

    PMtotI = MapTotal(*PMmw) + MapTotal(*PMms) + MapTotal(*PMsoil);
    return(PMtotI);
}


//---------------------------------------------------------------------------
/**
* @fn void TWorld::KinematicPestMC(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD, cTMap *_Qn, cTMap *_Qsn,
                             cTMap *_Qpn, cTMap *_Qpsn, cTMap *_PCmw, cTMap *_PCms, cTMap *_PCrw, cTMap *_PCrs,
                             cTMap *_Alpha,cTMap *_DX, cTMap *_Sed)
* @brief Adaptation of kinematic wave routing for pesticides.
* @return Concentrations, fluxes and new mass states of the pesticides in the different domains.
*
*/
/*LDD_COOR *_crlinked_*/
void TWorld::KinematicPestMC(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD, cTMap *_Qn, cTMap *_Qsn,
                             cTMap *_Qpn, cTMap *_Qpsn, cTMap *_PCmw, cTMap *_PCms, cTMap *_PCrw, cTMap *_PCrs,
                             cTMap *_Alpha,cTMap *_DX, cTMap *_Sed)
{
   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        SpinKW->Drc = 0;
        QpinKW->Drc = 0;
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
        double rho = 2650;

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
        // calculate concentrations for Crw_in and Crs_in
        double Crw_in = Qpin/Qin;
        double Crs_in = Spin/Sin;

        SpinKW->Drc = Spin;
        QpinKW->Drc = Qpin;
        // change of mass pesticide in soil under mixing layer
        double PMsoil_out = 0;
        if (Ez->Drc < 0) {
            //deposition
            PMsoil_out = Ez->Drc * PCms->Drc * _dx * _dx * rho;
        } else {
            // erosion
            PMsoil_out = -Ez->Drc * PCs->Drc * _dx * _dx;
        }
        // declare output vars for simplePestConc()
        double Kd = KdPestMC;
        double Kfilm = KfilmPestMC;
        double Kr = KrPestMC;
        double Crw_n, Crs_n, Cmw_n, Cms_n, Cinf_n = 0;
        simplePestConc(_PCrw->Drc, _PCmw->Drc, Kfilm, InfilVol->Drc, zm->Drc, Kr, Kd, _PCrs->Drc, _PCms->Drc, Ez->Drc, Sed->Drc, CHAdjDX->Drc, ThetaS1->Drc,
                       Crw_in, Crs_in,
                       Crw_n, Crs_n, Cmw_n, Cms_n, Cinf_n);
        // MC - do we use A = dx^2 or A = dx * flowwidth

        // The four new concentrations
        PCrw->Drc = Crw_n;
        PCrs->Drc = Crs_n;
        PCmw->Drc = Cmw_n;
        PCms->Drc = Cms_n;

        // The new fluxes
        // no more outflow than total pesticide in domain
        PQrs->Drc = std::min(PCrs->Drc * Qsn->Drc, SpinKW->Drc + PMrs->Drc/_dt);
        // for water first substract infiltration than runoff - TODO!!
        double Qinf = InfilVol->Drc;
        PQinf->Drc = std::min(Cinf_n * Qinf, QpinKW->Drc + PMrw->Drc/_dt);
        PMrw->Drc = std::max(0.0, PMrw->Drc - (PQinf->Drc * _dt) + (QpinKW->Drc * _dt));
        PQrw->Drc = std::min(PCrw->Drc * Qn->Drc, QpinKW->Drc + PMrw->Drc/_dt);

        // The new masses
        // new mass based on all fluxes and original pesticide present
        PMrw->Drc = std::max(0.0, PMrw->Drc - (PQrw->Drc * _dt) + (QpinKW->Drc * _dt));
        PMrs->Drc = std::max(0.0, PMrs->Drc + (SpinKW->Drc * _dt) - (PQrs->Drc * _dt));
        PMmw->Drc = PCmw->Drc * _dx * _dx * zm->Drc * ThetaS1->Drc;
        // assuming the mixing zone is always saturated
        PMms->Drc = PCms->Drc * _dx * _dx * zm->Drc * rho;
        PMsoil->Drc = PMsoil->Drc + PMsoil_out;


    }
}
