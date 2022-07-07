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
#include <tuple>

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
// NOT USED ANY MORE
template <typename... T>
void TWorld::simplePestConc(double Crw_old, double Cmw_old, double Kfilm, double Qinf, double zm, double kr, double Kd,
                            double Crs_old, double Cms_old, double Ez, double Me, double A, double pore,
                            double Crw_in, double Crs_in, std::tuple<T...> all_conc) // MC - this are the target values to be returned
{
    double rho = 2650;
        // calculate average between C of fluxes
        double Crw_avg = 0.5 * (Crw_old + Crw_in);
        double Crs_avg = 0.5 * (Crs_old + Crs_in);
        double Crw_n, Crs_n, Cmw_n, Cms_n, Cinf_n = 0;

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
        all_conc = std::make_tuple(Crw_n, Crs_n, Cmw_n, Cms_n, Cinf_n);
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
void TWorld::MassPest(double PMtotI, double &PMerr, double &PMtot)
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
double TWorld::MassPestInitial(void)
{
    double PMtotI = 0;
    double rho = 2650;
    // PMtotI = PMsoil + PMmw + PMms
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        PMms->Drc = PCms->Drc * SoilWidthDX->Drc * _dx * rho * zm->Drc;
        PMmw->Drc = PCmw->Drc * ThetaI1->Drc * zm->Drc * SoilWidthDX->Drc * _dx;
        PMsoil->Drc = PCms->Drc * SoilWidthDX->Drc * _dx * zs->Drc;
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
                             cTMap *_Qpn, cTMap *_Qpsn, cTMap *_PCmw, cTMap *_PCms, cTMap *_PCrw, cTMap *_PCrs, double _dx, cTMap *_Sed)
{
   //MC - change the function to use the map DX - which is _dx adjusted for slope.

   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        SpinKW->Drc = 0;
        QpinKW->Drc = 0;
        Ez->Drc = 0;
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
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L{
        SpinKW->Drc = Spin;
        QpinKW->Drc = Qpin;
        double Kd = KdPestMC;
        double Kfilm = KfilmPestMC;
        double Kr = KrPestMC;

//        std::tuple <double, double, double, double, double> all_conc;
//
//        simplePestConc(PCrw->Drc, PCmw->Drc, Kfilm, InfilVol->Drc, zm->Drc, Kr, Kd, PCrs->Drc, PCms->Drc, Ez->Drc, Sed->Drc, CHAdjDX->Drc, ThetaS1->Drc,
//                       Crw_in, Crs_in, all_conc);

        // MC - do we use A = dx^2 or A = dx * flowwidth or A = dx * SoilWidth
        // Dx = cellsize, SoilWidth = width of soil (mixing soil interaction and option for erosion),
        // FlowWidth = SoilWidth + Roads and hard surface, water flows over this area, but for hardsurface no interaction with mixing layer. Deposition on this area.

        // alternative to simplePestConc calculation. Instead of concentration calculate masses
        // per timestep there are 2 mass change processes: 1. exchange between domains and 2. influx and outflux.

        // 1. exchange between domains
        // calculate the correct C's based on the Qin and Qold etc.
        double Crw_avg, Crs_avg = 0;
        if (Qn->Drc + Qin > 0) {
        Crw_avg = (PMrw->Drc + (Qpin * _dt)) /((Qn->Drc + Qin) * _dt);
        }

        if (Qsn->Drc + Sin > 0) {
        Crs_avg = (PMrs->Drc + (Spin * _dt)) /((Qsn->Drc + Sin) * _dt);
        }

        // calculate the mass redistributions
        double Qinf = InfilVol->Drc;
        double Mmw_ex, Mms_ex, Mrw_ex, Mrs_ex = 0; // exchange mass (negative is reduction, positive is increase)
        // units : mg = mg - (m * sec-1 * (mg * m-3) * m * m * sec
        Mmw_ex -= ((zm->Drc * rho * Kr * (Kd * PCmw->Drc - PCms->Drc))/ (ThetaS1->Drc * zm->Drc)) *_dt * SoilWidthDX->Drc * _dx // exchange with soil in mixing (mg)
                 - (Kfilm*(PCmw->Drc - Crw_avg)* SoilWidthDX->Drc * _dx * _dt); // exchange between mixing water and runoff water (mg)

        // m * - * m * sec-1 * ( kg * m -3) / - * m
        // zm * rho * kr * (Kd * Cmw_old - Cms_old))/ pore * zm)
        Mrw_ex = (Kfilm*(PCmw->Drc - Crw_avg)* SoilWidthDX->Drc * _dx * _dt); // exchange with mixing layer water

        // calculate erosion depth
        Ez->Drc = (DEP->Drc + DETFlow->Drc + DETSplash->Drc);
        // change of mass pesticide in soil under mixing layer
        double PMsoil_out = 0;
        if (Ez->Drc < 0) {
            //deposition
            Ez->Drc = Ez->Drc / rho * _dx * FlowWidth->Drc;
            PMsoil_out = Ez->Drc * PCms->Drc * FlowWidth->Drc * _dx * rho;
            Mms_ex = ((zm->Drc * rho * Kr * (Kd * PCmw->Drc - PCms->Drc))/ (ThetaS1->Drc * zm->Drc)) *_dt * SoilWidthDX->Drc * _dx // exchange between soil water in mixing zone
                    - (Crs_avg * Ez->Drc * rho * _dx * FlowWidth->Drc); // added by deposition. problem with dep on roads!!!
            Mrs_ex = (Crs_avg * (Ez->Drc * rho * _dx * FlowWidth->Drc)); // loss by deposition
        } else {
            // erosion
            Ez->Drc = Ez->Drc / rho * _dx * SoilWidthDX->Drc;
            PMsoil_out = -Ez->Drc * PCs->Drc * SoilWidthDX->Drc * _dx;
            Mms_ex = ((zm->Drc * rho * Kr * (Kd * PCmw->Drc - PCms->Drc))/ (ThetaS1->Drc * zm->Drc)) *_dt * SoilWidthDX->Drc * _dx // mixing layer
                    - (PCms->Drc * Ez->Drc * rho * _dx * SoilWidthDX->Drc); // loss by erosion
            Mrs_ex = (PCms->Drc * Ez->Drc * rho * _dx * SoilWidthDX->Drc); // added by erosion
        }


        // 2. influx and outflux
        double Mrw_inf, Mmw_inf, Mrs_out, Mrw_out = 0;
        Mrw_inf -= Qinf * Crw_avg; // loss through infiltration from runoff
        Mmw_inf = Qinf * (Crw_avg - PCmw->Drc); // net loss through infiltration from mixing layer (mg)

        // no more outflow than total pesticide in domain
        PQrs->Drc = std::min(PCrs->Drc * Qsn->Drc, Spin + (PMrs->Drc - Mrs_ex)/_dt);
        PQrw->Drc = std::min(PCrw->Drc * Qn->Drc, Qpin + PMrw->Drc / _dt);
        Mrs_out = PQrs->Drc * _dt;
        Mrw_out = PQrw->Drc * _dt;

        // The new masses
        // new mass based on all fluxes and original pesticide present
        PMrw->Drc = std::max(0.0, PMrw->Drc + Mrw_ex - Mrw_out - Mmw_inf + (Qpin * _dt));
        PMrs->Drc = std::max(0.0, PMrs->Drc + Mrs_ex - Mrs_out + (Spin * _dt));
        PMmw->Drc = std::max(0.0, PMmw->Drc + Mmw_ex + Mmw_inf);
        // assuming the mixing zone is always saturated
        PMms->Drc = std::max(0.0, PMms->Drc + Mms_ex);
        PMsoil->Drc += PMsoil_out;

        // The four new concentrations
        double Volmw, Massms = 0;
        Volmw = zm->Drc * _dx * SoilWidthDX->Drc * ThetaS1->Drc;
        Massms = zm->Drc * _dx * SoilWidthDX->Drc * rho * (1-ThetaS1->Drc);
        PCrw->Drc = PMrw->Drc / WaterVolall->Drc;
        PCrs->Drc = PMrs->Drc / Sed->Drc;
        PCmw->Drc = PMmw->Drc / Volmw;
        PCms->Drc = PMms->Drc / Massms;
        }}
    }
}
