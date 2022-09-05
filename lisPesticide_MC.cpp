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
  \brief Transport and partitioning of pesticides with explicit method

functions: \n
- double TWorld::MassPest(double PMtotI, double &PMerr, double &PMtot) \n
- double TWorld::MassPestInitial(void) \n
- void TWorld::KinematicPestMC(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD,
                             cTMap *_Qn, cTMap *_Qsn, cTMap *_Qpwn, cTMap *_Qpsn,
                             cTMap *_DX, cTMap *_Alpha, cTMap *_Sed,
                             cTMap *_Q, cTMap *_Qs, cTMap *_Qpw, cTMap *_Qps) \n
- void TWorld::SimplePestConc() -- NOT USED ANYMORE--\n
*/

#include "model.h"
#include "operation.h"

// check if cell From flows to To
//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::MassPest(double WaterVolall, double Sed, double Qsn, double Qn, double Qinf)
 * @brief Calculation of pesticide mass at the end of each timestep.
 * Add all sinks and dynamic sources together to calculate total mass in the system.
 * @param  PMtotI: pesticide mass initially in the system - mg
 * @param  PMerr: the error of the mass balance - mg OR %
 * @param  PMtot: the total mass of pesticides in the system - this timestep - mg
 * @return updates PestOutW, PestOutS, PMtot and PMerr
 *
 */
void TWorld::MassPest(double PMtotI, double &PMerr, double &PMtot)
{
   // totals of outfluxes
   // at the moment only overland flow is accounted for. add channels etc later.
   FOR_ROW_COL_LDD5 {
       PQrw_dt = PQrw->Drc * _dt;
       if (SwitchErosion) {
            PQrs_dt = PQrs->Drc * _dt;
       }
    }}

    PestOutW += PQrw_dt; 
    Pestinf += mapTotal(*PMinf);
    PestPerc += mapTotal(*PMperc);
    double PMerosion {0.0};
    if (SwitchErosion) {
        PestOutS += PQrs_dt;
        PMerosion = mapTotal(*PMrs) + PestOutS;
    }

    PMtot = Pestinf + PestOutW + mapTotal(*PMsoil) + mapTotal(*PMrw)
            + mapTotal(*PMmw) + mapTotal(*PMms) + PestPerc + PMerosion;

  //  PMerr = 1 - (PMtot/PMtotI);
    PMerr = PMtot - PMtotI;
}

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::MassPestInitial(void)
 * @brief Calculation total pesticide mass initial in system.
 * @param DX: cell length adjusted for slope [m]
 * @param PCms: map with initial pesticide concentration of soil in mixing zone - mg kg-1
 * @param PCmw: map with initial pesticide concentration of water in mixing zone - mg L-1
 * @param zm: map with thickness of mixing layer [m]
 * @param zs: map with thickness of soil layer with pesticides [m]
 * @param ThetaI1: map with initial soil moisture of the first soil layer. - L L-1
 * @return PMtotI
 *
 */
double TWorld::MassPestInitial(void)
{
    double pmtot_i {0.0};
    double rho = rhoPestMC;
    // PMtotI = PMsoil + PMmw + PMms
    //#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        // mg = mg kg-1 * m * m * kg m-3 * m
        PMms->Drc = PCms->Drc * SoilWidthDX->Drc * DX->Drc * rho * zm->Drc;
        // mg = mg L-1 * m * m * m * 1000
        PMmw->Drc = PCmw->Drc * ThetaI1->Drc * zm->Drc * SoilWidthDX->Drc
                    * DX->Drc * 1000;
        // mg = mg kg-1 * m * m * m * kg m-3
        PMsoil->Drc = PCms->Drc * SoilWidthDX->Drc * DX->Drc
                      * (zs->Drc - zm->Drc) * rho;
    }}
    pmtot_i = mapTotal(*PMmw) + mapTotal(*PMms) + mapTotal(*PMsoil);
    return(pmtot_i);
}

//---------------------------------------------------------------------------
/**
* @fn void TWorld::KinematicPestMC(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD,
                             cTMap *_Qn, cTMap *_Qsn, cTMap *_Qpwn, cTMap *_Qpsn,
                             cTMap *_DX, cTMap *_Alpha, cTMap *_Sed,
                             cTMap *_Q, cTMap *_Qs, cTMap *_Qpw, cTMap *_Qps)
* @brief Do all the calculations for pesticide dynamics
* This includes:
*   - exchange of mass in the mixing layer
*   - infiltration and percolation
*   - discharge dissolved in water or adsorbed to sediment
* Kinamtic wave is solved with complexSedCalc()
* @return Concentrations, fluxes and new mass states of the
* pesticides in the different domains.
*
*/
/*LDD_COOR *_crlinked_*/
void TWorld::KinematicPestMC(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD,
                             cTMap *_Qn, cTMap *_Qsn, cTMap *_Qpwn, cTMap *_Qpsn,
                             cTMap *_DX, cTMap *_Alpha, cTMap *_Sed,
                             cTMap *_Q, cTMap *_Qs, cTMap *_Qpw, cTMap *_Qps)
// also add? : cTMap *_PMrw, cTMap *_PMrs, cTMap *_WaterVolall ???
{
   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

   FOR_ROW_COL_MV_L {
       _Qpwn->Drc = 0;
       _Qpsn->Drc = 0;
       QpinKW->Drc = 0;
       SpinKW->Drc = 0;
   }}

//#pragma omp parallel for reduction(+:Qin) num_threads(userCores)
    for(long i_ =  0; i_ < _crlinked_.size(); i_++) //_crlinked_.size()
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;

        double Qin {0}; //m3 sec-1
        double Qpin {0}; //mg sec-1
        double Sin {0}; //kg sec-1
        double Spin {0}; //mg sec-1

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
                        Qpin += _Qpwn->Drcr;
                         if (SwitchErosion) {
                            Sin += _Qsn->Drcr;
                            Spin += _Qpsn->Drcr;
                         }
                    }
                }
            }
        }
        //#pragma omp parallel for num_threads(userCores)
        QpinKW->Drc = Qpin;
        if (SwitchErosion) {
            Ez->Drc = 0;
            SpinKW->Drc = Spin;
        }
        double rho = rhoPestMC;     //kg m-3
        double Kd = KdPestMC;       // -
        double Kfilm = KfilmPestMC; // m sec-1
        double Kr = KrPestMC;       // sec
        double Mda_ex {0.0};        // mg - exchange mixing layer
        double Mwrm_ex {0.0};       // mg - exchange water runoff - mixing layer
        double mass_perc {0.0};     // mg - mass lost by percolation
        double Mrw_inf {0.0};       // mg - mass of infiltration from runoff
        double Mmw_inf {0.0};       // mg - mass of infiltration from mixing layer
        double Theta_mix {0.0};     // m3/m3 - soil moisture mixing layer
        double Msrm_ex {0.0};       // mg - exchange sediment runoff - mixing layer
        double PMsoil_ex {0.0};     // mass exchange deeper soil - mixing layer

        // MC - do we use A = dx^2 or A = dx * flowwidth or A = dx * SoilWidth
        // Dx = cellsize, SoilWidth = width of soil (mixing soil interaction and option for erosion),
        // FlowWidth = SoilWidth + Roads and hard surface, water flows over this area, b
        // but for hardsurface no interaction with mixing layer. Deposition on this area.
        // MC - 220815 - Used A = SoilWidthDX->Drc * _DX->Drc for now.

        // complexity is increased in situations (1 - 4)

// 1. no runoff, no erosion, no infiltration - only equilibrium in mixing zone
        // exchange between adsorbed and dissolved in mixing zone
        // unclear how exchange between sediment and water in mixing zone works for units.
        // for now the assumption is made that the resulting unit of
        // Kr * (Kd * PCmw->Drc - PCms->Drc) is mg * kg-1 * sec-1 !!!

        //Exchange in mixing layer
        // positive adds to absorbed, negative to dissolved.
        // mg = mg kg-1 sec-1 * kg m-3 * m * m * m * sec
        Mda_ex = (Kr * (Kd * PCmw->Drc - PCms->Drc) * rho * _dt
                  * SoilWidthDX->Drc * _DX->Drc * zm->Drc);

        // loss by percolation
        if (InfilVol->Drc < 1e-6) {
            if (Perc->Drc < 1e-7) {
                mass_perc = 0;
            } else {
            // calculate volume of percolated water from mixing layer
            double perc_vol {0.0}, perc_rat {0.0};
            perc_rat = Perc->Drc / (SoilDepth1->Drc - Lw->Drc);
            // L = (m * m * m * 1000)
            perc_vol = zm->Drc * _DX->Drc * SoilWidthDX->Drc * perc_rat * 1000;
            mass_perc = perc_vol * PCmw->Drc;
            }
            Theta_mix = Thetaeff->Drc; //percolation related theta mixing layer
        }

// 2. no runoff, no erosion, infiltration
        //no runoff - add leftover of mass in runoff water to mixing layer
        if (_Qn->Drc + Qin < 1e-6) {
            PCrw->Drc = 0.0;        //concentration = 0
            PMmw->Drc += PMrw->Drc; //add any leftover mass to mixing layer
            PMrw->Drc = 0.0;        // mass = 0
            PQrw->Drc = 0.0;        // discharge = 0
        }
        double Crw_avg {0.0};        // mg/L - no runoff; concentration = 0
        double Qinf = InfilVol->Drc; // m3 (per timestep)
        if (InfilVol->Drc > 1e-6) {
        // assume the mixing layer is saturated during infiltration or runoff.
        Theta_mix = ThetaS1->Drc;     
        // mg = m3 * 1000 * (mg L-1)
        Mmw_inf = Qinf * 1000 * PCmw->Drc; // infiltration mixing layer (mg)
        } // infiltration occurs

        // write infiltration and percolation to map
        PMinf->Drc = Mmw_inf;
        PMperc->Drc = mass_perc;
/* RUNOFF and SEDIMENT CALCULATIONS
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
 * Cj1i = the concentration 1 cell upstream this timestep, we know it already
 *        because we calculate from upstream to downstream
 * Cj1i1 = the concentration in current cell this timestep - we want to calculate this.
 *
 * For flux calculations we use the average concentration of Cj1i and Cji1 to go to Cj1i1.
 */

// 3. runoff, no erosion, infiltration
        if (_Qn->Drc + Qin > 1e-6) { // more than 1 ml - what is best definition of runoff?
            // calculate the correct C's based on the Qin and Qold etc.
            // mg L-1 = (mg + (mg sec-1 * sec)) / (m3 + (m3 sec-1 * sec)) * 1000(m3 -> L)
            Crw_avg = (PMrw->Drc + (Qpin * _dt))
                      / ((WaterVolall->Drc + (Qin * _dt)) * 1000);
            // only exchange with runoff water when there is a significant amount
            if (WH->Drc > 0.00001) { // more than 0.1 mm - slaat dit ergens op?
                // positive adds to runoff, negative to mixing layer.
                // mg = ((sec-1 * m * (mg L-1) / m) * m * m * m * sec * 1000 (m3 -> L)
                Mwrm_ex = (((Kfilm*(PCmw->Drc - Crw_avg))/zm->Drc)
                           * (SoilWidthDX->Drc * _DX->Drc * zm->Drc * _dt
                           * ThetaS1->Drc * 1000));
                // mg = m3 * 1000 (L->m3) * mg L-1
                Mrw_inf = Qinf * 1000 * Crw_avg; // loss through infiltration from runoff
            } // significant runoff

// 4. runoff, erosion and infiltration
        if (SwitchErosion) {
//                if (_Qsn->Drc + Sin < 1e-5) {
//                    PCrs->Drc = 0.0;        //concentration = 0
//                    PMms->Drc += PMrs->Drc; //add any leftover mass to mixing layer
//                    PMrs->Drc = 0.0;        // mass = 0
//                    PQrs->Drc = 0.0;        // discharge = 0
//                }
            double Crs_avg {0.0};
            if (_Sed->Drc > 1e-6) { // more than 0.01 gram _Qsn->Drc + Sin
            // mg kg-1 = (mg + (mg sec-1 * sec)) / (kg + kg sec-1 * sec)
                Crs_avg = (PMrs->Drc + (Spin * _dt))
                      / (_Sed->Drc + (Sin * _dt));
            // calculate erosion depth, no time component, per timestep
            // For now only use SoilWidth in formulas. Check what is doen with deposition on roads.
            // Can this be eroded after deposition or not?
            // option 1 - all deposition on roads add directly to sink
            // option 2 - deposition on roads can be eroded and added into the system...
                double eMass = (DEP->Drc + DETFlow->Drc + DETSplash->Drc); //kg/cell - sediment BulkDensity is part of runfile Conservation
                if (eMass < -1e-6) { //close to zero no calculations are done
                    //deposition
                    // m = kg / kg m-3 * m * m
                    Ez->Drc = eMass / (rho * _DX->Drc * FlowWidth->Drc); // also on road surface
                    PMsoil_ex = eMass * PCms->Drc; // what happens with pesticides on roads??
                    // mg = mg kg-1 * kg
                    Msrm_ex = Crs_avg * eMass; // loss by deposition
                } else if (eMass > 1e-6){
                    // erosion
                    Ez->Drc = eMass / (rho * _DX->Drc * SoilWidthDX->Drc); // only on soil surface
                    PMsoil_ex = eMass * PCs->Drc; // * SoilWidthDX->Drc * _DX->Drc * rho;
                    // mg = mg kg-1  kg
                    Msrm_ex = PCms->Drc * eMass; // * rho * _DX->Drc * SoilWidthDX->Drc); // added by erosion
                }
            // adjust mass for erosion or deposition
            // no more transport than mass in cell
                if (PMrs->Drc + Msrm_ex < 0) {
                    Msrm_ex = -PMrs->Drc;
                    PMrs->Drc = 0;
                    PCrs->Drc = 0;
                    _Qpsn->Drc = 0;
                } else   {
                    PMrs->Drc = std::max(0.0, PMrs->Drc + Msrm_ex);

                // concentration for outflux
                PCrs->Drc = PMrs->Drc / _Sed->Drc;
                //_Qps->Drc = _Qs->Drc * PCrs->Drc;
// OPTION 1 - simple extrapolation
                double totpests = PMrs->Drc + (Spin * _dt);
                double totsed = _Sed->Drc + (Sin * _dt);
                _Qpsn->Drc = std::min(totpests/_dt,
                                  _Qsn->Drc * (totpests / totsed));
                }
            // mg = mg sec-1 * sec
            PMrs->Drc = std::max(0.0, PMrs->Drc - (_Qpsn->Drc * _dt)
                                          + (Spin * _dt));
            PCrs->Drc = PMrs->Drc / _Sed->Drc;
            // adjust lower soil layer
            PMsoil->Drc -= PMsoil_ex;
            zs->Drc -= Ez->Drc;
            PCs->Drc = PMsoil->Drc
                       / (zs->Drc * rho * _DX->Drc * SoilWidthDX->Drc);
            } // erosion occurs
        } //switch erosion

            //substract infiltration and mixing layer exchange
            //mg = mg + mg - mg - mg + (mg sec-1 * sec)
            PMrw->Drc = std::max(0.0, PMrw->Drc + Mwrm_ex - Mrw_inf);
            // calculate concentration for new outflux
            PCrw->Drc = PMrw->Drc / (WaterVolall->Drc * 1000);
            _Qpw->Drc = _Q->Drc * 1000 * PCrw->Drc;
// OPTION 1 - simple extrapolation
//            // mg sec-1 = mg L-1 * m3 sec-1 * 1000, mg sec-1 + mg / sec
//            _Qpwn->Drc = std::min(PCrw->Drc * _Qn->Drc * 1000,
//                                  QpinKW->Drc + (PMrw->Drc / _dt));
// OPTION 2 - use complexSedCalc explicit method
            //use the Sediment explicit approach for concentration in water
            _Qpwn->Drc = complexSedCalc(_Qn->Drc, Qin, _Q->Drc, Qpin, _Qpw->Drc,
                                        _Alpha->Drc, _DX->Drc); //mg/sec
            _Qpwn->Drc = std::min(_Qpwn->Drc, QpinKW->Drc + PMrw->Drc / _dt);

            //substract all masses
            //mg = mg - (mg sec-1 * sec)
            PMrw->Drc = std::max(0.0, PMrw->Drc - (_Qpwn->Drc * _dt)
                                          + (QpinKW->Drc * _dt));
            // calculate new concentration
            PCrw->Drc = PMrw->Drc / (WaterVolall->Drc * 1000);
         } // runoff occurs

        // mg = mg - mg - mg
        PMmw->Drc = std::max(0.0, PMmw->Drc - Mda_ex - mass_perc
                                      - Mmw_inf + Mrw_inf - Mwrm_ex);
        // mg = mg + mg
        PMms->Drc = std::max(0.0, PMms->Drc + Mda_ex - Msrm_ex + PMsoil_ex);
        double Volmw {0.0}, Massms {0.0};
        // L = m * m * m * -- * 1000
        Volmw = zm->Drc * _DX->Drc * SoilWidthDX->Drc * Theta_mix * 1000;
        // kg = m * m * m * kg m_3 * --
        Massms = zm->Drc * _DX->Drc * SoilWidthDX->Drc * rho;
        //mg L-1 = mg / L
        PCmw->Drc = PMmw->Drc / Volmw;
        //mg kg-1 = mg / kg
        PCms->Drc = PMms->Drc / Massms;
    } // end loop over ldd
}


//---------------------------------------------------------------------------
// NOT USED ANY MORE
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
    double rho {1400};
    // calculate average between C of fluxes
    double Crw_avg = 0.5 * (Crw_old + Crw_in);
    double Crs_avg = 0.5 * (Crs_old + Crs_in);
    double Crw_n {0}, Crs_n {0} , Cmw_n {0}, Cms_n {0}, Cinf_n {0};

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
