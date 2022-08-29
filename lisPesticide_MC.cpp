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
#include "iostream"

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
   // at the moment only overland flow is accounted for. is channels etc must be included, build that later.
   FOR_ROW_COL_LDD5 {
        PQrw_dt += PQrw->Drc * _dt;
        PQrs_dt += PQrs->Drc * _dt;
    }}

    // Pestinf += MapTotal(*PMinf);
    PestOutW += PQrw_dt;
    PestOutS += PQrs_dt;
    //PestPerc += mapTotal(*perc_out);

    PMtot = Pestinf + PestOutS + PestOutW + mapTotal(*PMsoil) + mapTotal(*PMrw)
            + mapTotal(*PMrs) + mapTotal(*PMmw) + mapTotal(*PMms) + PestPerc;

  //  PMerr = 1 - (PMtot/PMtotI);
    PMerr = PMtot;

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
    double pmtot_i {0.0};
    double rho = rhoPestMC;
    // PMtotI = PMsoil + PMmw + PMms
    //#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        // mg = mg kg-1 * m * m * kg m-3 * m
        PMms->Drc = PCms->Drc * SoilWidthDX->Drc * _dx * rho * zm->Drc;
        // mg = mg L-1 * m * m * m * 1000
        PMmw->Drc = PCmw->Drc * ThetaI1->Drc * zm->Drc * SoilWidthDX->Drc * _dx * 1000;
        // mg = mg kg-1 * m * m * m * kg m-3
        PMsoil->Drc = PCms->Drc * SoilWidthDX->Drc * _dx * (zs->Drc - zm->Drc) * rho;      
    }}
    pmtot_i = mapTotal(*PMmw) + mapTotal(*PMms) + mapTotal(*PMsoil);
    return(pmtot_i);
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
                             cTMap *_Qpn, cTMap *_Qpsn, double _dx, cTMap *_Sed)
{
   //MC - change the function to use the map DX - which is _dx adjusted for slope.

   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

//#pragma omp parallel for reduction(+:Qin) num_threads(userCores)
    for(long i_ =  0; i_ < _crlinked_.size(); i_++) //_crlinked_.size()
    {
        int r = _crlinked_[i_].r;
        int c = _crlinked_[i_].c;

        double Qin {0}; //m3 * sec-1
        double Qpin {0}; //mg * sec-1
        if (SwitchErosion) {
            double Sin {0}; //kg * sec-1
            double Spin {0}; //mg * sec-1
        }
        double Sin {0}; //kg * sec-1
        double Spin {0}; //mg * sec-1
        double rho = rhoPestMC; //kg * m3-1

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
                        Qpin += _Qpn->Drcr;
                         if (SwitchErosion) {
                            Sin += _Qsn->Drcr;
                            Spin += _Qpsn->Drcr;
                         }
                    }
                }
            }
        }
        //#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L{
        QpinKW->Drc = 0;
        QpinKW->Drc = Qpin;
        if (SwitchErosion) {
            SpinKW->Drc = 0;
            Ez->Drc = 0;
            SpinKW->Drc = Spin;
        }
        double Kd = KdPestMC; // -
        double Kfilm = KfilmPestMC; // m sec-1
        double Kr = KrPestMC; // sec

//        std::tuple <double, double, double, double, double> all_conc;
//        simplePestConc(PCrw->Drc, PCmw->Drc, Kfilm, InfilVol->Drc, zm->Drc, Kr, Kd, PCrs->Drc, PCms->Drc, Ez->Drc, Sed->Drc, CHAdjDX->Drc, ThetaS1->Drc,
//                       Crw_in, Crs_in, all_conc);

        // MC - do we use A = dx^2 or A = dx * flowwidth or A = dx * SoilWidth
        // Dx = cellsize, SoilWidth = width of soil (mixing soil interaction and option for erosion),
        // FlowWidth = SoilWidth + Roads and hard surface, water flows over this area, but for hardsurface no interaction with mixing layer. Deposition on this area.
        // MC - 220815 - Used A = SoilWidthDX->Drc * _dx for now, can also include DX (adjusted for slope).

        // alternative to simplePestConc calculation. Instead of concentration calculate masses
        // per timestep there are 2 mass change processes: 1. exchange between domains and 2. influx and outflux.
        // complexity is increased in situations

        // 1. no runoff, no erosion, no infiltration - only equilibrium in mixing zone
        // exchange between adsorbed and dissolved in mixing zone
        // unclear how exchange between sediment and water in mixing zone works for units. for now the assumption is made that the resulting unit of
        // Kr * (Kd * PCmw->Drc - PCms->Drc) is mg * kg-1 * sec-1 !!!
        double Mda_ex {0.0}; // mg - mass dissolved - absorbed exchange
        double Mmw_ex {0.0}, Mms_ex {0.0}, Mrw_ex {0.0}, Mrs_ex {0.0}, a {0.0}, b {0.0}; // exchange mass (negative is reduction, positive is increase) - mg
        double Mrw_inf {0.0}, Mmw_inf {0.0}, Mrs_out {0.0}, Mrw_out {0.0}; // influx and outflux
        double Theta_mix {0.0}; // depending on the situation theta_mix is thetaeff or thetas
        double mass_perc {0.0};

        // positive adds to absorbed, negative to dissolved.
        // mg = mg kg-1 sec-1  * kg m-3 * m * m * m * sec
        Mda_ex = (Kr * (Kd * PCmw->Drc - PCms->Drc) * rho * _dt * SoilWidthDX->Drc * _dx * zm->Drc);

        // loss by percolation
        //perc_out->Drc = 0.0; // when there is infiltration percolation is neglected.
        if (InfilVol->Drc < 1e-6) {
            // calculate volume of percolated water from mixing layer
            double perc_vol {0.0}, perc_rat {0.0};
            perc_rat = Perc->Drc /(SoilDepth1->Drc - Lw->Drc);
            // L = (m * m * m * 1000)
            perc_vol = zm->Drc * _dx * SoilWidthDX->Drc * perc_rat * 1000;
            mass_perc = perc_vol * PCmw->Drc;

            Theta_mix = Thetaeff->Drc; //assign percolation related theta to mixing layer

        }
        // 2. no runoff, no erosion, infiltration
        if (InfilVol->Drc > 1e-6) {
        //infiltration
        double Qinf = InfilVol->Drc; // m3 (per timestep)
        // we assume that the mixing zone will be saturated if there is infiltration and/or runoff.
        Theta_mix = ThetaS1->Drc;
        double Crw_avg {0.0}; // mg/L - no runoff, so concentration will be 0
        // mg = m3 * 1000 (L->m3) * mg L-1
        Mrw_inf = Qinf * 1000 * Crw_avg; // loss through infiltration from runoff
        // mg = m3 * 1000 * (mg L-1)
        Mmw_inf = Qinf * 1000 * (PCmw->Drc - Crw_avg); // net loss through infiltration from mixing layer (mg)

//        // 3. runoff, no erosion, infiltration
//        // calculate the correct C's based on the Qin and Qold etc.
//        if (Qn->Drc + Qin > 1e-6) { // more than 1 ml - what is best definition of runoff?
//            // mg L-1 = (mg + (mg sec-1 * sec)) / m3 sec-1 * sec * 1000(m3 -> L)
//            Crw_avg = (PMrw->Drc + (Qpin * _dt)) /((Qn->Drc + Qin) * _dt * 1000);

//            // only exchange with runoff water when there is a significant amount
//            if (WH->Drc > 0.0001) { // more than 0.1 mm
//            // should we add zm in the next function?
//            // mg = ((sec-1 * m * (mg L-1) / m) * m * m * m * sec * 1000 (m3 -> L)
//            Mrw_ex = (((Kfilm*(PCmw->Drc - Crw_avg))/zm->Drc) * (SoilWidthDX->Drc * _dx * zm->Drc * _dt * ThetaS1->Drc * 1000)); // exchange with mixing layer water
//            // conc = mg * L-1 * sec-1 = ((m * sec-1 * (mg * L-1 - mg * L-1)) /  m * -
//            // mass = mg = (((m sec-1 * (mg L-1 - mg L-1)) / m * - ) * ( m * m * m * sec * 1000(m3 -> L)
//            b = - (((Kfilm *(PCmw->Drc - Crw_avg))/zm->Drc) * (SoilWidthDX->Drc * _dx * zm->Drc * _dt * ThetaS1->Drc * 1000)); // exchange between mixing water and runoff water (mg)
//            } // runoff occurs
//        // 4. runoff, erosion and infiltration
//        Crs_avg = 0;
//        if (SwitchErosion) {
//            if (Qsn->Drc + Sin > 0.0001) { // more than 0.1 gram
//            // mg kg-1 = (mg + (mg sec-1 * sec)) / kg sec-1 * sec
//            Crs_avg = (PMrs->Drc + (Spin * _dt)) /((Qsn->Drc + Sin) * _dt);

//            // calculate erosion depth, no time component in this formulas, this is already covered by Ez
//            Ez->Drc = (DEP->Drc + DETFlow->Drc + DETSplash->Drc);
//            // change of mass pesticide in soil under mixing layer
//            double PMsoil_out = 0;
//            if (Ez->Drc < -0.00001) { //close to zero no calculations are done
//                //deposition
//                Ez->Drc = Ez->Drc / rho * _dx * FlowWidth->Drc; // also on road surface
//                PMsoil_out = Ez->Drc * PCms->Drc * FlowWidth->Drc * _dx * rho;
//                // - (mg kg-1 * kg m-3 * -m * m * m)
//                a = - (Crs_avg * Ez->Drc * rho * _dx * FlowWidth->Drc) + PMsoil_out; // added by deposition. problem with dep on roads!!!
//                // mg = mg kg-1 * (-m * kg m-3 * m * m)
//                Mrs_ex = (Crs_avg * (Ez->Drc * rho * _dx * FlowWidth->Drc)); // loss by deposition
//            } else if (Ez->Drc > 0.00001){
//                // erosion
//                Ez->Drc = Ez->Drc / rho * _dx * SoilWidthDX->Drc; // only on soil surface
//                PMsoil_out = -Ez->Drc * PCs->Drc * SoilWidthDX->Drc * _dx * rho;
//                // (mg kg-1 * +m * kg m-3 * m * m)
//                a = - (PCms->Drc * Ez->Drc * rho * _dx * SoilWidthDX->Drc) + PMsoil_out; // loss by erosion
//                // mg = (mg kg-1 * +m + kg m-3 * m * m)
//                Mrs_ex = (PCms->Drc * Ez->Drc * rho * _dx * SoilWidthDX->Drc); // added by erosion
//            }
//            // mg = mg + mg - mg + (mg sec-1 * sec)
//            PMrs->Drc = std::max(0.0, PMrs->Drc + Mrs_ex - Mrs_out + (Spin * _dt));
//            PMsoil->Drc += PMsoil_out;
//            } // erosion occurs
//        } //switch erosion
//        //mg = mg + mg - mg - mg + (mg sec-1 * sec)
//        PMrw->Drc = std::max(0.0, PMrw->Drc + Mrw_ex - Mrw_out - Mrw_inf + (Qpin * _dt));
//        // no more outflow than total pesticide in domain
//        // mg sec-1 = mg kg-1 * kg sec-1, mg sec-1 + mg / sec
//        PQrs->Drc = std::min(PCrs->Drc * Qsn->Drc, Spin + (PMrs->Drc + Mrs_ex)/ _dt);
//        // mg sec-1 = mg L-1 * m3 sec-1 * 1000, mg sec-1 + mg / sec -- ??? already substract Mrw_ex here ???
//        PQrw->Drc = std::min(PCrw->Drc * Qn->Drc * 1000, Qpin + (PMrw->Drc + Mrw_ex)/ _dt);
//        // mg = mg sec-1 * sec
//        Mrs_out = PQrs->Drc * _dt;
//        Mrw_out = PQrw->Drc * _dt;
//        } // runoff occurs
        // The four new concentrations
//        PCrw->Drc = PMrw->Drc / WaterVolall->Drc;
//        PCrs->Drc = PMrs->Drc / Sed->Drc;

        } // infiltration occurs


        // write infiltration and percolation to map and sum for mass balance
        // MC 220826 -- the mass balance is correct now. Percolation and soilmosture are not exactly stable - for now not problematic.
        PMinf->Drc = Mmw_inf;
        Pestinf += Mmw_inf;
        PMperc->Drc = mass_perc; //map with percolation losses
        PestPerc += mass_perc;

        // mg = mg - mg - mg
        PMmw->Drc = std::max(0.0, PMmw->Drc - Mda_ex - mass_perc - Mmw_inf); // dit lijkt goed te werken (update PMmw to i)
        // mg = mg + mg
        PMms->Drc = std::max(0.0, PMms->Drc + Mda_ex);
        double Volmw {0.0}, Massms {0.0};
        // L = m * m * m * -- * 1000
        Volmw = zm->Drc * _dx * SoilWidthDX->Drc * Theta_mix * 1000;
        // kg = m * m * m * kg m_3 * --
        Massms = zm->Drc * _dx * SoilWidthDX->Drc * rho;
        //mg L-1 = mg / L
        PCmw->Drc = PMmw->Drc / Volmw;
        //mg kg-1 = mg / kg
        PCms->Drc = PMms->Drc / Massms;
       }}
    } // end loop over ldd
}
