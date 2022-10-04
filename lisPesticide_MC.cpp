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
  \brief Transport and partitioning of pesticides

functions: \n
- void TWorld::PesticideDynamicsMC(void)
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
- void TWorld::KinematicPestMC(QVector <LDD_COORIN> _crlinked_, cTMap *_LDD,
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
* @fn void TWorld::PesticideDynamicsMC(void)
* @brief Do all the calculations for pesticide dynamics
* This includes:
*   call functions for specific dynamics
*   update all masses etc.
*/

void TWorld::PesticideDynamicsMC(void)
{
    double rho = rhoPestMC;     //kg m-3
    double Kd = KdPestMC;       // -
    double Kfilm = KfilmPestMC; // m sec-1
    double kr = KrPestMC;       // sec

    FOR_ROW_COL_MV_L{
       // no runoff, no erosion
       double mda_ex {0.0};        // mg - exchange mixing layer

       // exchange between adsorbed and dissolved in mixing zone

       // for now the assumption is made that the resulting unit of
       // Kr * (Kd * PCmw->Drc - PCms->Drc) is mg * kg-1 * sec-1
       // this hold if 1L water = 1kg

       // positive adds to absorbed, negative to dissolved.
       // mg = mg kg-1 sec-1 * kg m-3 * m * m * m * sec
       mda_ex = (kr * (Kd * PCmw->Drc - PCms->Drc) * rho * _dt
                 * SoilWidthDX->Drc * DX->Drc * zm->Drc);

       //percolation
       PMperc->Drc = 0.0; //does not need to be a map...
       if (InfilVol->Drc < 1e-6) {
           PMperc->Drc = PesticidePercolation(Perc->Drc, SoilDepth1->Drc,
                        Lw->Drc, zm->Drc, DX->Drc, SoilWidthDX->Drc, PCmw->Drc);
           Theta_mix->Drc = Thetaeff->Drc; //percolation related theta for mixing layer
       }

       //infiltration from mixing layer to deeper soil
       PMinf->Drc = 0.0; //does not need to be a map...
       if (InfilVol->Drc > 1e-6) {
           // assume the mixing layer is saturated during infiltration or runoff.
           Theta_mix->Drc = ThetaS1->Drc;
           // mg = m3 * 1000 * (mg L-1)
           PMinf->Drc = InfilVol->Drc * 1000 * PCmw->Drc; // infiltration mixing layer (mg)
       }

       PMmw->Drc = std::max(0.0, PMmw->Drc - mda_ex);

       // adjust masses if outflow if more than available mass.
       if (PMmw->Drc < PMinf->Drc + PMperc->Drc) {
           double tot = PMinf->Drc + PMperc->Drc;
           PMinf->Drc = (PMinf->Drc/tot) * PMmw->Drc;
           PMperc->Drc = (PMperc->Drc/tot) * PMmw->Drc;
           PCmw->Drc = 0;
       }
       // update mass after percolation and infiltration
       // mg = mg - mg - mg
       PMmw->Drc = std::max(0.0, PMmw->Drc - PMinf->Drc - PMperc->Drc);
       // mg = mg + mg
       PMms->Drc = std::max(0.0, PMms->Drc + mda_ex);
       // update PCmw before dissolved
       PCmw->Drc = PMmw->Drc / (zm->Drc * DX->Drc * SoilWidthDX->Drc * Theta_mix->Drc * 1000);
    }}

    //runoff

    KinematicPestDissolvedCombined(crlinkedldd_, LDD, Qn, PQrw, DX, Alpha, Q, Qpw,
                        Kfilm);

    //erosion
    if(SwitchErosion){
        KinematicPestAdsorbed(crlinkedldd_, LDD, Qsn, PQrs, DX, Alpha, Sed,
                              Qs, Qps, rho);
    }
    // calculate new concentration
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
    if (perc < 1e-7) {
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
    if (Qn->Drc + QinKW->Drc <= 1e-6) {
        PCrw->Drc = 0.0;        //concentration = 0
        PMmw->Drc += PMrw->Drc; //add any leftover mass to mixing layer
        pmwdep->Drc -= PMrw->Drc;
        PMrw->Drc = 0.0;        // mass = 0
        PQrw->Drc = 0.0;        // discharge = 0

    }
    if (Qn->Drc + QinKW->Drc > 1e-6) { // more than 1 ml - what is best definition of runoff?
        // calculate the correct C's based on the Qin and Qold etc.
        // mg L-1 = (mg + (mg sec-1 * sec)) / (m3 + (m3 sec-1 * sec)) * 1000(m3 -> L)
        Crw_avg = (PMrw->Drc + (QpinKW->Drc * _dt))
                  / ((WaterVolin->Drc + (QinKW->Drc * _dt)) * 1000);
        // positive adds to runoff, negative to mixing layer.
        // mg = ((sec-1 * m-1 * (mg L-1) / m) * m * m * m * sec * 1000 (m3 -> L)
        mwrm_ex = (((_kfilm * (PCmw->Drc - Crw_avg))/zm->Drc)
                   * (SoilWidthDX->Drc * _DX->Drc * zm->Drc * _dt
                   * ThetaS1->Drc * 1000));
        // mg = m3 * 1000 (L->m3) * mg L-1
        mrw_inf = InfilVol->Drc * 1000 * Crw_avg; // loss through infiltration from runoff
        // adjust masses if outflow if more than available mass.
        if (PMrw->Drc < mrw_inf - mwrm_ex) {
            double tot = mrw_inf - mwrm_ex;
            mrw_inf = (mrw_inf/tot) * PMrw->Drc;
            mwrm_ex = (mwrm_ex/tot) * PMrw->Drc;
        }
        mwrm_ex > 0 ? pmwdet->Drc += mwrm_ex : pmwdep->Drc += mwrm_ex;
        pmwdep->Drc -= mrw_inf;
        //substract infiltration and mixing layer exchange
        //mg = mg + mg - mg - mg + (mg sec-1 * sec)
        PMrw->Drc = std::max(0.0, PMrw->Drc + mwrm_ex - mrw_inf);
        // calculate concentration for new outflux
        PCrw->Drc = PMrw->Drc / (WaterVolin->Drc * 1000);
        _Qpw->Drc = _Q->Drc * 1000 * PCrw->Drc;
        // use explicit backwards method from Chow
        _Qpwn->Drc = QpwSeparate(_Qn->Drc, QinKW->Drc, _Q->Drc, QpinKW->Drc, _Qpw->Drc,
                                    _Alpha->Drc, _DX->Drc); //mg/sec
        _Qpwn->Drc = std::min(_Qpwn->Drc, QpinKW->Drc + PMrw->Drc / _dt);

        //substract discharge
        //mg = mg - (mg sec-1 * sec)
        PMrw->Drc = std::max(0.0, PMrw->Drc - (_Qpwn->Drc * _dt)
                                      + (QpinKW->Drc * _dt));
        PMmw->Drc = std::max(0.0, PMmw->Drc - mwrm_ex + mrw_inf);
       } //runoff occurs
    }//end ldd loop
}

//---------------------------------------------------------------------------
/**
* @fn double TWorld::KinematicPestAdsorbed(double perc, double soildep,
*               double lw, double zm, double dx, double swdx, double pcmw)
* @brief Calculate mass lost by percolation
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
    double msoil_ex {0.0};
    double msrm_ex {0.0};

    //no erosion - add leftover of mass to mixing layer
    if (Sed->Drc < 1e-6) {
        PCrs->Drc = 0.0;        //concentration = 0
        PMms->Drc += PMrs->Drc; //add any leftover mass to mixing layer
        pmsdep->Drc -= PMrs->Drc;
        PMrs->Drc = 0.0;        // mass = 0
        PQrs->Drc = 0.0;        // discharge = 0

    }

    if (_Sed->Drc > 1e-6) { // more than 0.01 gram _Qsn->Drc + Sin
        // mg kg-1 = (mg + (mg sec-1 * sec)) / (kg + kg sec-1 * sec)
        Crs_avg = (PMrs->Drc + (SpinKW->Drc * _dt))
                  / (_Sed->Drc + (SinKW->Drc * _dt));
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
            msoil_ex = eMass * PCms->Drc; // what happens with pesticides on roads??
            // mg = mg kg-1 * kg
            msrm_ex = Crs_avg * eMass; // loss by deposition            
        } else if (eMass > 1e-6){
            // erosion
            Ez->Drc = eMass / (rho * _DX->Drc * SoilWidthDX->Drc); // only on soil surface
            msoil_ex = eMass * PCs->Drc; // * SoilWidthDX->Drc * _DX->Drc * rho;
            // mg = mg kg-1  kg
            msrm_ex = PCms->Drc * eMass; // * rho * _DX->Drc * SoilWidthDX->Drc); // added by erosion            
        }
        // adjust mass for erosion or deposition
        // no more transport than mass in cell
        if (PMrs->Drc + msrm_ex < 0) {
            msrm_ex = -PMrs->Drc;
        }
        PMrs->Drc = std::max(0.0, PMrs->Drc + msrm_ex);
        eMass < 0 ? pmsdep->Drc += msrm_ex: pmsdet->Drc += msrm_ex;

        // concentration for outflux
        PCrs->Drc = PMrs->Drc / _Sed->Drc;
        //_Qps->Drc = _Qs->Drc * PCrs->Drc;
        // - simple extrapolation
        double totpests = PMrs->Drc + (SpinKW->Drc * _dt);
        double totsed = _Sed->Drc + (SinKW->Drc * _dt);
        _Qpsn->Drc = std::min(totpests/_dt,
                                  _Qsn->Drc * (totpests / totsed));

        // mg = mg sec-1 * sec
        PMrs->Drc = std::max(0.0, PMrs->Drc - (_Qpsn->Drc * _dt)
                                      + (SpinKW->Drc * _dt));
        PCrs->Drc = PMrs->Drc / _Sed->Drc;
        // adjust lower soil layer
        PMsoil->Drc -= msoil_ex;
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
double TWorld::QpwSeparate(double Qj1i1, double Qj1i, double Qji1,double Pj1i, double Pji1, double alpha, double dx)
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

    A = _dt*Pj1i;
    B = -Cavg*abQb_1*(Qj1i1-Qji1)*dx;
    C = (Qji1 <= MIN_FLUX ? 0 : aQb*Pji1/Qji1)*dx;
    if (Qj1i1 > MIN_FLUX)
        Pj1i1 = (A+C+B)/(_dt+aQb*dx/Qj1i1);
    else
        Pj1i1 = 0;
    return std::max(0.0 ,Pj1i1);
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
double TWorld::QpwInfExCombined(double Qj1i1, double Qj1i, double Qji1,
                                double Pj1i, double Pji1, double alpha,
                                double dx, double zm, double kfilm, double qinf,
                                double cmw)
{
    double Pj1i1 {0}, Cp_avg {0}, Q_avg {0};
    double A {0}, B {0}, C {0}, D {0}; // aux vars
    double E {0}, F {0}, G {0}, H {0}; // aux vars
    const double beta = 0.6;
    cmw = cmw * 1000; // mg/L -> mg/m3

    Q_avg = 0.5*(Qji1+Qj1i); //m3/sec
    if (Q_avg <= MIN_FLUX)
        return (0);
    Cp_avg = (Pj1i+Pji1)/(Qj1i+Qji1); //mg/m3
    // calculate all parts of formula
    A = Pj1i * _dt;
    B = alpha * pow(Q_avg, beta) * dx * Pji1 / Qji1;
    C = alpha * beta * pow(Q_avg, beta-1) * (Qj1i1 - Qji1) * dx;
    D = kfilm * cmw * zm * dx * _dt;
    //
    E = _dt;
    F = alpha * pow(Q_avg, beta) * dx / Qj1i1;
    G = kfilm * zm * dx * _dt / Qj1i1;
    H = qinf * dx * _dt * zm / Qj1i1;

    // calculate new Qp - mg/sec
    Pj1i1 = (A+B-C-D) / (E+F-G-H);
    return std::max(0.0 ,Pj1i1);
}

//---------------------------------------------------------------------------
/**
* @fn double TWorld::KinematicPestDissolved(double perc, double soildep,
*               double lw, double zm, double dx, double swdx, double pcmw)
* @brief Calculate dissolved concentration of pesticides in runoff water
* with infiltration, exchange and runoff in same implicit solution.
*/

void TWorld::KinematicPestDissolvedCombined(QVector <LDD_COORIN> _crlinked_,
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
    double mwrm_ex {0.0};   //mg - mass exchange
    double mrw_inf {0.0};   //mg - mass infiltration
    double mrw_q {0.0};     //mg - mass new outflux

    //no runoff - add leftover of mass in runoff water to mixing layer
    if (Qn->Drc + QinKW->Drc <= 1e-6) {
        PCrw->Drc = 0.0;        //concentration = 0
        PMmw->Drc += PMrw->Drc; //add any leftover mass to mixing layer
        pmwdep->Drc -= PMrw->Drc;
        PMrw->Drc = 0.0;        // mass = 0
        PQrw->Drc = 0.0;        // discharge = 0

    }
    if (Qn->Drc + QinKW->Drc > 1e-6) { // more than 1 ml - what is best definition of runoff?
        Qpw->Drc = PCrw->Drc * 1000 * Q->Drc;
        // calculate new Qp
        _Qpwn->Drc = QpwInfExCombined(Qn->Drc, QinKW->Drc, Q->Drc,
                         QpinKW->Drc, Qpw->Drc, Alpha->Drc,
                         DX->Drc, zm->Drc, KfilmPestMC, fact->Drc,
                         PCmw->Drc);
        Crwn->Drc = Qn->Drc > 0 ? _Qpwn->Drc / (Qn->Drc * 1000) : 0; // mg/L

        //calculate new masses

        // positive adds to runoff, negative to mixing layer.
        // mg = ((sec-1 * m * (mg L-1) / m) * m * m * m * sec * 1000 (m3 -> L)
        mwrm_ex = (((_kfilm * (PCmw->Drc - Crwn->Drc))/zm->Drc)
                   * (SoilWidthDX->Drc * _DX->Drc * zm->Drc * _dt
                   * ThetaS1->Drc * 1000));
        // mg = m3 * 1000 (L->m3) * mg L-1
        mrw_inf = InfilVol->Drc * 1000 * Crwn->Drc; // loss through infiltration
        // mg = mg/L * L/sec * sec
        mrw_q = PQrw->Drc * _dt;

        // adjust masses if outflow is more than available mass.
        if (PMrw->Drc < mrw_inf - mwrm_ex) {
            double tot = mrw_inf - mwrm_ex;
            mrw_inf = (mrw_inf/tot) * PMrw->Drc;
            mwrm_ex = (mwrm_ex/tot) * PMrw->Drc;
        }
        mwrm_ex > 0 ? pmwdet->Drc += mwrm_ex : pmwdep->Drc += mwrm_ex;
        pmwdep->Drc -= mrw_inf;
//        if (PMrw->Drc + QpinKW->Drc * _dt < mrw_inf + mrw_q) {
//            double tot = mrw_inf + mrw_q;
//            mrw_inf = (mrw_inf/tot) * PMrw->Drc;
//            mrw_q = (mrw_q/tot) * PMrw->Drc;
//        }
//        PQrw->Drc = mrw_q / _dt;
        PMrw->Drc = std::max(0.0, PMrw->Drc - mrw_inf + mwrm_ex);
        PQrw->Drc = std::min(PQrw->Drc, QpinKW->Drc + PMrw->Drc / _dt);

        //substract infiltration and discharge
        //mg = mg + mg - mg - mg + (mg sec-1 * sec)
        PMrw->Drc = std::max(0.0, PMrw->Drc - PQrw->Drc * _dt + QpinKW->Drc * _dt);
        PMmw->Drc = std::max(0.0, PMmw->Drc - mwrm_ex + mrw_inf);
       } //runoff occurs
    }//end ldd loop
}
