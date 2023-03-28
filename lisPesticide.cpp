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
**  Author of the pesticide code: Meindert Commelin
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisPesticide.cpp
  \brief Transport and partitioning of pesticides

functions: \n
- void TWorld::PesticideDynamics(void)
- double TWorld::MassPest(double PMtotI, double &PMerr, double &PMtot) \n
- double TWorld::MassPestInitial(void) \n
- double TWorld::PesticidePercolation(double perc, double soildep, double lw,
                double zm, double dx, double swdx, double pcmw)
- void TWorld::KinematicPestDissolved(QVector <LDD_COORIN> _crlinked_,
               cTMap *_LDD, cTMap *_Qn, cTMap *_Qpwn, cTMap *_DX,
               cTMap *_Alpha, cTMap *_Q, cTMap *_Qpw, double _kfilm)
- void TWorld::KinematicPestAdsorbed(QVector <LDD_COORIN> _crlinked_,
                             cTMap *_LDD, cTMap *_Qsn, cTMap *_Qpsn, cTMap *_DX,
                             cTMap *_Alpha, cTMap *_Sed, cTMap *_Qs, cTMap *_Qps,
                                   double rho)
//OBSOLETE
- void TWorld::KinematicPest(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD,
                             cTMap *_Qn, cTMap *_Qsn, cTMap *_Qpwn, cTMap *_Qpsn,
                             cTMap *_DX, cTMap *_Alpha, cTMap *_Sed,
                             cTMap *_Q, cTMap *_Qs, cTMap *_Qpw, cTMap *_Qps) \n
*/

#include "model.h"
#include "operation.h"

// check if cell From flows to To
//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::MassPest(double PMtotI, double &PMerr, double &PMtot)
 * @brief Calculation of pesticide mass at the end of each timestep.
 * Add all sinks and dynamic sources together to calculate total mass in the system.
 * @param  PMtotI: pesticide mass initially in the system - mg
 * @param  PMerr: the error of the mass balance - mg OR %
 * @param  PMtot: the total mass of pesticides in the system - this timestep - mg
 * @return updates PestOutW, PestOutS, PMtot and PMerr
 *
 */
void TWorld::MassPest(double PMtotI, double &PMerr, double &PMtot, double &PMserr, double &PMwerr)
{
   // totals of outfluxes
   // at the moment only overland flow is accounted for. add channels etc later.
   FOR_ROW_COL_LDD5 {
       PQrw_dt = PQrw->Drc * _dt;
       if (SwitchErosion) {
            PQrs_dt = PQrs->Drc * _dt;
       }
    }}

    //parallel??
    //#pragma omp parallel for num_threads(userCores)
    PestOutW += PQrw_dt; 
    Pestinf += mapTotal(*PMinf);
    PestPerc += mapTotal(*PMperc);
    double PMerosion {0.0};
    if (SwitchErosion) {
        PestOutS += PQrs_dt;
        PMerosion = mapTotal(*PMrs) + PestOutS;
    }

    // all pesticide mass in system current timestep
    PMtot = Pestinf + PestOutW + mapTotal(*PMsoil) + mapTotal(*PMrw)
            + mapTotal(*PMmw) + mapTotal(*PMms) + PestPerc + PMerosion;
    PMerr = (PMtot - PMtotI) / PMtotI * 100;

    // mass balance for active adsorbed
    PMserr = 0;
    if (SwitchErosion) {
    double PMsdep {0.0};
    double PMsdet {0.0};
    PMsdep = mapTotal(*pmsdep);
    PMsdet = mapTotal(*pmsdet);
    PMserr = PMsdet > 0 ? (PMsdet + PMsdep - PMerosion) / PMsdet * 100 : 0;
    }

    // mass balance active dissolved
    double PMwactive {0.0};
    double PMwdep {0.0};
    double PMwdet {0.0};
    PMwactive = PestOutW + mapTotal(*PMrw);
    PMwdep = mapTotal(*pmwdep);
    PMwdet = mapTotal(*pmwdet);
    PMwerr = PMwdet > 0 ? (PMwdet + PMwdep - PMwactive) / PMwdet * 100 : 0;
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
    double rho = rhoPest;
    // PMtotI = PMsoil + PMmw + PMms
    //#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        // mg = mg kg-1 * m * m * kg m-3 * m
        PMms->Drc = PCms->Drc * SoilWidthDX->Drc * DX->Drc * rho * zm->Drc;
        // mg = mg L-1 * m * m * m * 1000
        PMmw->Drc = PCmw->Drc * ThetaI1->Drc * zm->Drc * SoilWidthDX->Drc
                    * DX->Drc * 1000;
        // mg = mg kg-1 * m * m * m * kg m-3
        PMsoil->Drc = PCs->Drc * SoilWidthDX->Drc * DX->Drc
                      * zs->Drc * rho;
    }}
    pmtot_i = mapTotal(*PMmw) + mapTotal(*PMms) + mapTotal(*PMsoil);
    return(pmtot_i);
}


//---------------------------------------------------------------------------
/**
* @fn void TWorld::PesticideDynamics(void)
* @brief Do all the calculations for pesticide dynamics
* This includes:
*   call functions for specific dynamics
*   update all masses etc.
*/

void TWorld::PesticideDynamics(void)
{
    double rho = rhoPest;     //kg m-3
    double Kd = KdPest;       // -
    double Kfilm = KfilmPest; // m sec-1
    double kr = KrPest;       // sec
    // this chunck can also be parallel
    //#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
       // no runoff, no erosion
       double mda_ex {0.0};        // mg - exchange mixing layer
       double mda_tot {0.0};       // mg - total mass in both phases
       double eql_ads {0.0};       // mg - adsorbed mass in equilibrium
       double eql_diss {0.0};      // mg - dissolved mass in equilibrium
       double m_diff {0.0};        // mg - current mass - eql mass
       double vol_w {0.0};         // l - volume water in mixing layer
       double mass_s {0.0};        // kg - mass sediment in mixing layer

       // exchange between adsorbed and dissolved in mixing zone
       // kg = m * m * m * kg L-1 * 1000
       vol_w = zm->Drc * DX->Drc * SoilWidthDX->Drc * Theta_mix->Drc * 1000;
       // kg = m * m * m * kg m-3
       mass_s = zm->Drc * DX->Drc * SoilWidthDX->Drc * rho;
       // for now the assumption is made that the resulting unit of
       // Kr * (Kd * PCmw->Drc - PCms->Drc) is mg * kg-1 * sec-1
       // this holds if 1L water = 1kg

       // positive adds to absorbed, negative to dissolved.
       // mg = mg kg-1 sec-1 *  sec * kg
       mda_ex = kr * (Kd * PCmw->Drc - PCms->Drc) * _dt
                * std::min(vol_w, mass_s);

       // calculate equilibrium mass division
       mda_tot = PMmw->Drc + PMms->Drc;
       eql_diss = mda_tot / (1 + (Kd / vol_w * mass_s));
       eql_ads = mda_tot - eql_diss;
       // mda_ex can not be larger than m_diff
       m_diff = eql_ads - PMms->Drc;
       mda_ex = std::abs(mda_ex) > std::abs(m_diff) ? m_diff : mda_ex;

       //percolation
       PMperc->Drc = 0.0; //does not need to be a map...
       if (InfilVol->Drc < 1e-6) {
           PMperc->Drc = PesticidePercolation(Perc->Drc, SoilDepth1->Drc,
                        Lw->Drc, zm->Drc, DX->Drc, SoilWidthDX->Drc, PCmw->Drc);
           Theta_mix->Drc = Thetaeff->Drc; //percolation related theta for mixing layer
       }

       //infiltration from mixing layer to deeper soil
       PMinf->Drc = 0.0; //does not need to be a map...
       if (InfilVol->Drc > tiny) {
           // assume the mixing layer is saturated during infiltration or runoff.
           Theta_mix->Drc = ThetaS1->Drc;
           // mg = m3 * 1000 * (mg L-1)
           PMinf->Drc = InfilVol->Drc * 1000 * PCmw->Drc; // infiltration mixing layer (mg)
       }

       // adjust masses if outflow is more than available mass.
       if (PMmw->Drc < PMinf->Drc + PMperc->Drc + mda_ex) {
           double tot = PMinf->Drc + PMperc->Drc + mda_ex;
           PMinf->Drc = (PMinf->Drc/tot) * PMmw->Drc;
           PMperc->Drc = (PMperc->Drc/tot) * PMmw->Drc;
           mda_ex = (mda_ex/tot) * PMmw->Drc;
           PCmw->Drc = 0;
       }
       // update mass after percolation and infiltration
       // mg = mg - mg - mg
       PMmw->Drc = std::max(0.0, PMmw->Drc - PMinf->Drc - PMperc->Drc - mda_ex);
       // mg = mg + mg
       PMms->Drc = std::max(0.0, PMms->Drc + mda_ex);
       // update PCmw before dissolved & PCms before adsorbed
       PCmw->Drc = PMmw->Drc / vol_w;
       PCms->Drc = PMms->Drc / mass_s;
    }}

    //runoff

    KinematicPestDissolved(crlinkedldd_, LDD, Qn, PQrw, DX, Alpha, Q, Qpw,
                        Kfilm);

    //erosion
    if(SwitchErosion){
        KinematicPestAdsorbed(crlinkedldd_, LDD, Qsn, PQrs, DX, Alpha, Sed,
                              Qs, Qps, rho);
    }
    // calculate new concentration
    // make parralell
    //#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
    double volmw {0.0};     //L - volume of water in mixing layer
    double massms {0.0};        // kg - mass of sediment in mixing layer
    if (WaterVolall->Drc > 0) {
    PCrw->Drc = PMrw->Drc / (WaterVolall->Drc * 1000);
    } else PCrw->Drc = 0.0;
    // L = m * m * m * -- * 1000
    volmw = zm->Drc * DX->Drc * SoilWidthDX->Drc * Theta_mix->Drc * 1000;
    PCmw->Drc = PMmw->Drc / volmw; //

    // kg = m * m * m * kg m_3 * --
    massms = zm->Drc * DX->Drc * SoilWidthDX->Drc * rho;
    //mg kg-1 = mg / kg
    PCms->Drc = PMms->Drc / massms;
    }}
}

//---------------------------------------------------------------------------
/**
* @fn double TWorld::PesticidePercolation(double perc, double soildep,
*               double lw, double zm, double dx, double swdx, double pcmw)
* @brief Calculate mass lost by percolation
*/

double TWorld::PesticidePercolation(double perc, double soildep, double lw,
                double zm, double dx, double swdx, double pcmw)
{
    double mass_perc {0.0}; // mg
    if (perc < tiny) {
        mass_perc = 0;
    } else {
        // calculate volume of percolated water from mixing layer
        double perc_vol {0.0}, perc_rat {0.0};
        perc_rat = perc / (soildep - lw);
        // L = (m * m * m * 1000)
        perc_vol = zm * dx * swdx * perc_rat * 1000;
        mass_perc = perc_vol * pcmw;
    }
    return(mass_perc);
}

//---------------------------------------------------------------------------
/**
* @fn double TWorld::KinematicPestDissolved(double perc, double soildep,
*               double lw, double zm, double dx, double swdx, double pcmw)
* @brief Calculate mass lost by percolation
*/

void TWorld::KinematicPestDissolved(QVector <LDD_COORIN> _crlinked_,
               cTMap *_LDD, cTMap *_Qn, cTMap *_Qpwn, cTMap *_DX,
               cTMap *_Alpha, cTMap *_Q, cTMap *_Qpw, double _kfilm)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_MV_L {
        _Qpwn->Drc = 0;
        QpinKW->Drc = 0;
    }}

//#pragma omp parallel for reduction(+:Qin) num_threads(userCores)
for(long i_ =  0; i_ < _crlinked_.size(); i_++) //_crlinked_.size()
{
    int r = _crlinked_[i_].r;
    int c = _crlinked_[i_].c;

    double Qin {0};         //m3 sec-1
    double Qpin {0};        //mg sec-1

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
                }
            }
        }
    }
    QpinKW->Drc = Qpin;
    //#pragma omp parallel for num_threads(userCores)
    double mwrm_ex {0.0};   //mg
    double mrw_inf {0.0};   //mg
    double Crw_avg {0.0};   // mg/L - no runoff; concentration = 0
    //no runoff - add leftover of mass in runoff water to mixing layer
    if (Qn->Drc + QinKW->Drc <= MIN_FLUX) {
        PCrw->Drc = 0.0;        //concentration = 0
        PMmw->Drc += PMrw->Drc; //add any leftover mass to mixing layer
        pmwdep->Drc -= PMrw->Drc;
        PMrw->Drc = 0.0;        // mass = 0
        PQrw->Drc = 0.0;        // discharge = 0

    }
    if (Qn->Drc + QinKW->Drc > MIN_FLUX) { // more than 1 ml - what is best definition of runoff?
        double vol_mw {0.0};    // volume of water in mixing layer [L]
        double vol_rw {0.0};    // volume of water in runoff [L]
        vol_mw = zm->Drc * Theta_mix->Drc * DX->Drc * SoilWidthDX->Drc * 1000;
        vol_rw = WaterVolin->Drc * 1000;

        // calculate the correct C's based on the Qin and Qold etc.
        // mg L-1 = (mg + (mg sec-1 * sec)) / (m3 + (m3 sec-1 * sec)) * 1000(m3 -> L)
        Crw_avg = (PMrw->Drc + (QpinKW->Drc * _dt))
                  / ((WaterVolin->Drc + (QinKW->Drc * _dt)) * 1000);

        // positive adds to runoff, negative to mixing layer.
        // mg = ((sec-1 * m-1 * (mg L-1) / m) * sec * L
        mwrm_ex = ((_kfilm * (PCmw->Drc - Crw_avg))/(zm->Drc + WHrunoff->Drc))
                   * _dt * std::min(vol_mw, vol_rw);

        double c_eql {0.0};
        double eql_mw {0.0};
        double m_diff {0.0};
        // equilibrium check
        // calculate equilibrium mass division
        c_eql = (PMmw->Drc + PMrw->Drc) / (vol_mw + vol_rw);
        eql_mw = c_eql * vol_mw; // mass in mixing layer at equilibrium
        // mwrm_ex can not be larger than m_diff
        m_diff = eql_mw - PMmw->Drc;
        mwrm_ex = std::abs(mwrm_ex) > std::abs(m_diff) ? m_diff : mwrm_ex;

        // mg = m3 * 1000 (L->m3) * mg L-1
        mrw_inf = InfilVol->Drc * 1000 * Crw_avg; // loss through infiltration from runoff

//        // adjust masses if outflow is more than available mass.
//        if (PMrw->Drc < mrw_inf - mwrm_ex) {
//            double tot = mrw_inf - mwrm_ex;
//            mrw_inf = (mrw_inf/tot) * PMrw->Drc;
//            mwrm_ex = (mwrm_ex/tot) * PMrw->Drc;
//        }
//        // adjust masses for PMMW
//        if (PMmw->Drc < mwrm_ex - mrw_inf) {
//            double tot = mwrm_ex - mrw_inf;
//            mrw_inf = (mrw_inf/tot) * PMmw->Drc;
//            mwrm_ex = (mwrm_ex/tot) * PMmw->Drc;
//        }


        //substract infiltration and mixing layer exchange
        double mrw_n {0.0}; // intermediate mass in runoff water
        //mg = mg + mg - mg - mg + (mg sec-1 * sec)
        mrw_n = std::max(0.0, PMrw->Drc + mwrm_ex - mrw_inf);
        // calculate concentration for new outflux
        PCrw->Drc = mrw_n / (WaterVolin->Drc * 1000);

        _Qpw->Drc = _Q->Drc * 1000 * PCrw->Drc;
        // use explicit backwards method from Chow
        _Qpwn->Drc = QpwSeparate(_Qn->Drc, QinKW->Drc, _Q->Drc, QpinKW->Drc, _Qpw->Drc,
                                    _Alpha->Drc, _DX->Drc, _dt); //mg/sec
        _Qpwn->Drc = std::min(_Qpwn->Drc, QpinKW->Drc + PMrw->Drc / _dt);

        // internal time loop

        //calculate courant number of standard timestep
        double Cr_rw {0.0};
        double Cr_mw {0.0};
        //double Cr_max = 0.9; // max Courant number
        //double dt_int_min = 0.1; // minimal timestep
        // runoff water
        Cr_rw = PMrw->Drc > 0 ? ((_Qpwn->Drc * _dt) - mwrm_ex + mrw_inf) / (PMrw->Drc + QpinKW->Drc * _dt) : 0.0;
        // mixing layer
        Cr_mw = PMmw->Drc > 0 ? (mwrm_ex - mrw_inf) / PMmw->Drc : 0.0;

        //start loop if one of the Cr's > Cr_max
        if (Cr_rw > Cr_max | Cr_mw > Cr_max) {
            // calculate steps and internal timestep
            double steps {0.0};
            double dt_int {0.0};
            steps = std::min(std::ceil(std::max(Cr_rw,Cr_mw)*(2/Cr_max)), _dt/dt_int_min);
            dt_int = _dt / steps;

            // fill intermediate concentrations and masses
            double int_Qpw, int_Qpwn, int_Cmw, int_Crw, int_Crw_avg; //
            double int_Mrw, int_Mmw, int_mwrm_ex, int_mrw_inf; //
            double sum_int_mwrm_ex {0.0}, sum_int_Qpwn {0.0}, sum_int_mrw_inf {0.0};

            int_Qpw = _Qpw->Drc;
            int_Cmw = PCmw->Drc;
            int_Mrw = PMrw->Drc;
            int_Mmw = PMmw->Drc;


            // make loop
            double count = 0;
            while (count < steps) {
                count++;

                // mg L-1 = (mg + (mg sec-1 * sec)) / (m3 + (m3 sec-1 * sec)) * 1000(m3 -> L)
                int_Crw_avg = (int_Mrw + (QpinKW->Drc * dt_int))
                          / ((WaterVolin->Drc + (QinKW->Drc * dt_int)) * 1000);
                // calculate mrwm_ex
                int_mwrm_ex = ((_kfilm * (int_Cmw - int_Crw_avg))/(zm->Drc + WHrunoff->Drc))
                              * dt_int * std::min(vol_mw, vol_rw);
                // equilibrium check
                // calculate equilibrium mass division
                c_eql = (int_Mmw + int_Mrw) / (vol_mw + vol_rw);
                eql_mw = c_eql * vol_mw; // mass in mixing layer at equilibrium
                // mwrm_ex can not be larger than m_diff
                m_diff = eql_mw - int_Mmw;
                int_mwrm_ex = std::abs(int_mwrm_ex) > std::abs(m_diff) ? m_diff : int_mwrm_ex;

                // calculate mrw_inf
                int_mrw_inf =  (InfilVol->Drc / steps) * 1000 * int_Crw_avg; // loss through infiltration

                // calculate concentration for new outflux
                int_Mrw = std::max(0.0, int_Mrw + int_mwrm_ex - int_mrw_inf);
                int_Crw = int_Mrw / (WaterVolin->Drc * 1000);

                // calculate Qpwn
                int_Qpwn = QpwSeparate(_Qn->Drc, QinKW->Drc, _Q->Drc, QpinKW->Drc, int_Qpw,
                                       _Alpha->Drc, _DX->Drc, dt_int);

                // add & substract all masses and update concentrations
                int_Qpwn = std::min(int_Qpwn, QpinKW->Drc + int_Mrw / dt_int);

                //substract infiltration and discharge
                //mg = mg + mg - mg - mg + (mg sec-1 * sec)
                int_Mrw = std::max(0.0, int_Mrw - int_Qpwn * dt_int + QpinKW->Drc * dt_int);
                int_Mmw = std::max(0.0, int_Mmw - int_mwrm_ex + int_mrw_inf);

                // for next internal loop Qpw = current Qpwn
                int_Qpw = int_Qpwn;
                // L = m * m * m * -- * 1000
                int_Cmw = int_Mmw / vol_mw;

                // mean Q and total exchange and infiltration
                sum_int_Qpwn += int_Qpwn;
                sum_int_mwrm_ex += int_mwrm_ex;
                sum_int_mrw_inf += int_mrw_inf;

            } // end internal time loop

            // calculate final concentrations and masses
            _Qpwn->Drc = sum_int_Qpwn / steps;
            mwrm_ex = sum_int_mwrm_ex;
            mrw_inf = sum_int_mrw_inf;

        } // end if Cr > Cr_max

        // mass balance
        mwrm_ex > 0 ? pmwdet->Drc += mwrm_ex : pmwdep->Drc += mwrm_ex;
        pmwdep->Drc -= mrw_inf;

        //substract discharge
        //mg = mg - (mg sec-1 * sec)
        PMrw->Drc = std::max(0.0, PMrw->Drc - (_Qpwn->Drc * _dt)
                                      + (QpinKW->Drc * _dt) + mwrm_ex - mrw_inf);
        PMmw->Drc = std::max(0.0, PMmw->Drc - mwrm_ex + mrw_inf);
       } //runoff occurs
    }//end ldd loop
}

//---------------------------------------------------------------------------
/**
* @fn double TWorld::KinematicPestAdsorbed(double perc, double soildep,
*               double lw, double zm, double dx, double swdx, double pcmw)
* @brief Calculate mass transported with runoff sediment
*/

void TWorld::KinematicPestAdsorbed(QVector <LDD_COORIN> _crlinked_,
                             cTMap *_LDD, cTMap *_Qsn, cTMap *_Qpsn, cTMap *_DX,
                             cTMap *_Alpha, cTMap *_Sed, cTMap *_Qs, cTMap *_Qps,
                                   double rho)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_MV_L {
        _Qpsn->Drc = 0;
        SpinKW->Drc = 0;
    }}

//#pragma omp parallel for reduction(+:Qin) num_threads(userCores)
for(long i_ =  0; i_ < _crlinked_.size(); i_++) //_crlinked_.size()
{
    int r = _crlinked_[i_].r;
    int c = _crlinked_[i_].c;

    double Sin {0}; //m3 sec-1
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
                    Sin += _Qsn->Drcr;
                    Spin += _Qpsn->Drcr;
                }
            }
        }
    }
    SpinKW->Drc = Spin;
    //} // end ldd loop
    //#pragma omp parallel for num_threads(userCores)
    //FOR_ROW_COL_MV_L {
    double Crs_avg {0.0};
    double msoil_ex {0.0};  // mass exchange between mixing layer and deeper soil
    double msrm_ex {0.0};   // mass exchange between mixing layer an suspended sediment

    //no erosion - add leftover of mass to mixing layer
    if (_Sed->Drc < tiny) {
        PCrs->Drc = 0.0;        //concentration = 0
        PMms->Drc += PMrs->Drc; //add any leftover mass to mixing layer
        pmsdep->Drc -= PMrs->Drc;
        PMrs->Drc = 0.0;        // mass = 0
        PQrs->Drc = 0.0;        // discharge = 0

    }

    if (_Sed->Drc > tiny) { // more than 0.01 gram _Qsn->Drc + Sin
        // positive = erosion, negative = deposition
        double eMass = (DEP->Drc + DETFlow->Drc + DETSplash->Drc); //kg/cell

        // mg kg-1 = (mg + (mg sec-1 * sec)) / (kg + kg sec-1 * sec)
        Crs_avg = (PMrs->Drc + (SpinKW->Drc * _dt))
                  / (_Sed->Drc - eMass + (SinKW->Drc * _dt));
        // substract eMass from _Sed for a more accurate concentration.
        // calculate erosion depth, no time component, per timestep
        // For now only use SoilWidth in formulas. Check what is done with deposition on roads.
        // Can this be eroded after deposition or not?
        // option 1 - all deposition on roads add directly to sink
        // option 2 - deposition on roads can be eroded and added into the system...

        if (eMass < 0) {
            //deposition
            // m = kg / kg m-3 * m * m
            Ez->Drc = eMass / (rho * _DX->Drc * SoilWidthDX->Drc); // also on road surface???
            msoil_ex = eMass * PCms->Drc; // what happens with pesticides on roads??
            // mg = mg kg-1 * kg
            msrm_ex = Crs_avg * eMass; // loss by deposition            
        } else if (eMass > 0){
            // erosion
            Ez->Drc = eMass / (rho * _DX->Drc * SoilWidthDX->Drc); // only on soil surface
            msoil_ex = eMass * PCs->Drc; //
            // mg = mg kg-1  kg
            msrm_ex = PCms->Drc * eMass; // added by erosion
        }
        // adjust mass for erosion or deposition
//        // no more transport than mass in cell domain
//        if (PMrs->Drc + msrm_ex < 0) {
//            msrm_ex = -PMrs->Drc;
//        }
//        if (PMms->Drc + msoil_ex < msrm_ex) {
//            msrm_ex = PMms->Drc + msoil_ex;
//        }

        // - simple extrapolation
        double totpests = std::max(0.0, PMrs->Drc + (SpinKW->Drc * _dt) + msrm_ex);
        double totsed = _Sed->Drc + (SinKW->Drc * _dt);
        _Qpsn->Drc = std::min(totpests/_dt,
                                  _Qsn->Drc * (totpests / totsed));


        // internal loop ----------------------------

        // calculate eMass rate : eMass/_dt
        //calculate courant number of standard timestep
        double Cr_rs {0.0};
        double Cr_ms {0.0};
        // runoff sediment
        Cr_rs = PMrs->Drc > 0 ? ((_Qpsn->Drc * _dt) - msrm_ex) / (PMrs->Drc + SpinKW->Drc * _dt) : 0.0;
        // mixing layer
        Cr_ms = PMmw->Drc > 0 ? (msrm_ex) / PMms->Drc : 0.0;

        //start loop if one of the Cr's > Cr_max
        if (Cr_rs > Cr_max | Cr_ms > Cr_max) {
            // calculate steps and internal timestep
            double steps {0.0};
            double dt_int {0.0};
            steps = std::min(std::ceil(std::max(Cr_rs,Cr_ms)*(2/Cr_max)), _dt/dt_int_min);
            dt_int = _dt / steps;

            // fill intermediate concentrations and masses
            double int_Qpsn {0.0}, int_Cms {0.0}, int_Crs {0.0}, int_Crs_avg {0.0}; //
            double int_Mrs {0.0}, int_Mms {0.0}, int_msrm_ex {0.0}, int_msoil_ex {0.0}; //
            double sum_int_msrm_ex {0.0}, sum_int_Qpsn {0.0}, sum_int_msoil_ex {0.0};
            double int_Ez {0.0}, int_eMass {0.0}, int_Cs {0.0}, int_zs {0.0}, int_Msoil {0.0};

            int_Cms = PCms->Drc;
            int_Mrs = PMrs->Drc;
            int_Mms = PMms->Drc;
            int_Cs = PCs->Drc;
            int_zs = zs->Drc;
            int_Msoil = PMsoil->Drc;

            int_eMass = (eMass / _dt) * dt_int;

            // make loop
            double count = 0;
            while (count < steps) {
                count++;

                // Crs
                // mg kg-1 = (mg + (mg sec-1 * sec)) / (kg + kg sec-1 * sec)
                int_Crs_avg = (int_Mrs + (SpinKW->Drc * dt_int))
                              / (_Sed->Drc - eMass + (SinKW->Drc * dt_int));
                // erosion or deposition
                if (int_eMass < 0) {
                    //deposition
                    // m = kg / kg m-3 * m * m
                    int_Ez = int_eMass / (rho * _DX->Drc * SoilWidthDX->Drc); // also on road surface???
                    int_msoil_ex = int_eMass * int_Cms; // what happens with pesticides on roads??
                    // mg = mg kg-1 * kg
                    int_msrm_ex = int_Crs_avg * int_eMass; // loss by deposition
                } else if (int_eMass > 0){
                    // erosion
                    int_Ez = int_eMass / (rho * _DX->Drc * SoilWidthDX->Drc); // only on soil surface
                    int_msoil_ex = int_eMass * int_Cs; //
                    // mg = mg kg-1  kg
                    int_msrm_ex = int_Cms * int_eMass; // added by erosion
                }

                // - simple extrapolation
                double totpests = std::max(0.0, int_Mrs + int_msrm_ex +
                                           (SpinKW->Drc * dt_int));
                double totsed = _Sed->Drc + (SinKW->Drc * dt_int);
                int_Qpsn = std::min(totpests/dt_int,
                                      _Qsn->Drc * (totpests / totsed));


               // update for next loop
               int_Mrs = std::max(0.0, int_Mrs - (int_Qpsn * dt_int)
                                           + (SpinKW->Drc * dt_int) + int_msrm_ex);
               int_Mms = std::max(0.0, int_Mms - int_msrm_ex + int_msoil_ex);
               int_Cms = int_Mms / (zm->Drc * rho * _DX->Drc * SoilWidthDX->Drc);
               int_zs -= int_Ez;
               int_Msoil = std::max(0.0, int_Msoil - int_msoil_ex);
               int_Cs  = int_Msoil / (int_zs * rho * _DX->Drc * SoilWidthDX->Drc);

               //sum for final values
               sum_int_Qpsn += int_Qpsn;
               sum_int_msoil_ex += int_msoil_ex;
               sum_int_msrm_ex += int_msrm_ex;


            } // end internal time loop

            //calculate final values
            _Qpsn->Drc = sum_int_Qpsn / steps;
            msoil_ex = sum_int_msoil_ex;
            msrm_ex = sum_int_msrm_ex;

        } // end if Cr > Cr_max

        // mass balance
        eMass < 0 ? pmsdep->Drc += msrm_ex: pmsdet->Drc += msrm_ex;
        // mg = mg sec-1 * sec
        PMrs->Drc = std::max(0.0, PMrs->Drc - (_Qpsn->Drc * _dt)
                                      + (SpinKW->Drc * _dt) + msrm_ex);
        PCrs->Drc = PMrs->Drc / _Sed->Drc;
        // adjust lower soil layer
        PMsoil->Drc = std::max(0.0, PMsoil->Drc - msoil_ex);
        zs->Drc -= Ez->Drc;
        PCs->Drc = PMsoil->Drc
                   / (zs->Drc * rho * _DX->Drc * SoilWidthDX->Drc);
        PMms->Drc = std::max(0.0, PMms->Drc - msrm_ex + msoil_ex);
        } // erosion occurs
    }// end ldd loop
}

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::QpwSeparate(double Qj1i1, double Qj1i, double Qji1,double Pj1i, double Pji1, double alpha, double dx)
 * @brief Explicit backward calculation of dissolved pesticide outflux from a cell
 *
 * Calculation of dissolved pesticide outflux from a cell based on a explicit solution of the time/space matrix,
 * j = time and i = place: j1i1 is the new output, j1i is the new flux at the upstream 'entrance' flowing into the gridcell
 *
 * @param Qj1i1 : result kin wave for this cell ( Qj+1,i+1 )
 * @param Qj1i : sum of all upstreamwater from kin wave ( Qj+1,i ),
 * @param Qji1 : incoming Q for kinematic wave (t=j) in this cell, map Q in LISEM (Qj,i+1)
 * @param Pj1i : sum of all upstream pesticide (Pj+1,i)
 * @param Pji1 : incoming dissolved pesticide for kinematic wave (t=j) in this cell, map Qpw in LISEM (Si,j+1)
 * @param alpha : alpha calculated in LISEM from before kinematic wave
 * @param dt : timestep
 * @param dx : length of the cell, corrected for slope (DX map in LISEM)
 * @return dissolved pesticide outflow in next timestep
 *
 */
double TWorld::QpwSeparate(double Qj1i1, double Qj1i, double Qji1,double Pj1i, double Pji1, double alpha, double dx, double dt)
{
    double Pj1i1, Cavg, Qavg, aQb, abQb_1, A, B, C;
    const double beta = 0.6;

    if (Qj1i1 < MIN_FLUX)
        return (0);

    Qavg = 0.5*(Qji1+Qj1i); //m3/sec
    if (Qavg <= MIN_FLUX)
        return (0);
    Cavg = (Pj1i+Pji1)/(Qj1i+Qji1); //mg/m3
    aQb = alpha*pow(Qavg,beta);
    abQb_1 = alpha*beta*pow(Qavg,beta-1);

    A = dt*Pj1i;
    B = -Cavg*abQb_1*(Qj1i1-Qji1)*dx;
    C = (Qji1 <= MIN_FLUX ? 0 : aQb*Pji1/Qji1)*dx;
    if (Qj1i1 > MIN_FLUX)
        Pj1i1 = (A+C+B)/(dt+aQb*dx/Qj1i1);
    else
        Pj1i1 = 0;
    return std::max(0.0 ,Pj1i1);
}
