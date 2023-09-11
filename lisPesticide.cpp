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
- double TWorld::MassPest(double PMtotI, double &PMerr, double &PMtot) \n
- double TWorld::MassPestInitial(void) \n
- void TWorld::PesticideCellDynamics(void) \n
- void TWorld::PesticideFlow1D(void) \n
- double TWorld::PesticidePercolation(double perc, double soildep, double lw,
                double zm, double dx, double swdx, double pcmw) \n
- void TWorld::KinematicPestDissolved(QVector <LDD_COORIN> _crlinked_,
               cTMap *_LDD, cTMap *_Qn, cTMap *_Qpwn, cTMap *_DX,
               cTMap *_Alpha, cTMap *_Q, cTMap *_Qpw, double _kfilm) \n
- void TWorld::KinematicPestAdsorbed(QVector <LDD_COORIN> _crlinked_,
                             cTMap *_LDD, cTMap *_Qsn, cTMap *_Qpsn, cTMap *_DX,
                             cTMap *_Alpha, cTMap *_Sed, cTMap *_Qs, cTMap *_Qps,
                                   double rho) \n
- double TWorld::QpwSeparate(double Qj1i1, double Qj1i, double Qji1,double Pj1i,
                             double Pji1, double alpha, double dx, double dt)
- void TWorld::PesticideSplashDetachment() \n
- void TWorld::PesticideFlowDetachment(double rho) \n
*/

#include "model.h"
#include "operation.h"

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
    //PestPerc += mapTotal(*PMperc);
    double PMerosion {0.0};
    if (SwitchErosion) {
        PestOutS += PQrs_dt;
        PMerosion = mapTotal(*PMrs) + PestOutS;
    }

    // all pesticide mass in system current timestep
    PMtot = Pestinf + PestOutW + mapTotal(*PMsoil) + mapTotal(*PMrw)
            + mapTotal(*PMmw) + mapTotal(*PMms) + PMerosion;
    PMerr = (PMtot - PMtotI) / PMtotI * 100;

    // mass balance for active sorbed
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

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        // mg = mg kg-1 * m * m * kg m-3 * m
        PMms->Drc = PCms->Drc * SoilWidthDX->Drc * DX->Drc * rho * zm->Drc;
        // mg = mg L-1 * m * m * m * 1000
        PMmw->Drc = PCmw->Drc * ThetaS1->Drc * zm->Drc * SoilWidthDX->Drc
                    * DX->Drc * 1000;
        // we use ThetaS because we assume saturation when the mixing zone is active.
        // mg = mg kg-1 * m * m * m * kg m-3
        PMsoil->Drc = PCs->Drc * SoilWidthDX->Drc * DX->Drc
                      * zs->Drc * rho;
    }}
    pmtot_i = mapTotal(*PMmw) + mapTotal(*PMms) + mapTotal(*PMsoil);
    return(pmtot_i);
}


//---------------------------------------------------------------------------
/**
* @fn void TWorld::PesticideCellDynamics(void)
* @brief Do all the calculations for dissolved pesticide dynamics in a cell
* This includes:
*   partitioning in mixing layer
*   losses by percolation and infiltration
*   update all masses etc.
*/

void TWorld::PesticideCellDynamics(void)
{
    double rho = rhoPest;     //kg m-3
    double Kd = KdPest;       // -
    double Kfilm = KfilmPest; // m sec-1
    double kr = KrPest;       // sec
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
       double mda_ex {0.0};        // mg - exchange mixing layer
       double mrw_inf {0.0};       // mg - mass from runoff to mixing layer by infiltration
       double mda_tot {0.0};       // mg - total mass in both phases
       double eql_ads {0.0};       // mg - sorbed mass in equilibrium
       double eql_diss {0.0};      // mg - dissolved mass in equilibrium
       double m_diff {0.0};        // mg - current mass - eql mass
       double vol_w {0.0};         // l - volume water in mixing layer
       double mass_s {0.0};        // kg - mass sediment in mixing layer

       // assume the mixing layer is saturated during infiltration or runoff.
       Theta_mix->Drc = ThetaS1->Drc;

       //infiltration from runoff through mixing layer to deeper soil
       PMinf->Drc = 0.0; //does not need to be a map...
       if (InfilVol->Drc > 0.0) {
           // mg = m3 * 1000 * (mg L-1)
           PMinf->Drc = InfilVol->Drc * 1000 * PCmw->Drc; // infiltration mixing layer (mg)

           // mg = m3 * 1000 (L->m3) * mg L-1
           mrw_inf = InfilVol->Drc * 1000 * PCrw->Drc; // loss through infiltration from runoff
       }

       // update mass after percolation and infiltration
       mrw_inf > PMrw->Drc ? mrw_inf = PMrw->Drc : mrw_inf;
       PMinf->Drc > PMmw->Drc ? PMinf->Drc = PMmw->Drc : PMinf->Drc;
       // mg = mg - mg - mg
       PMmw->Drc = std::max(0.0, PMmw->Drc - PMinf->Drc + mrw_inf);

       PMrw->Drc = std::max(0.0, PMrw->Drc - mrw_inf);

       pmwdep->Drc -= mrw_inf;
       // update PCmw before partitioning
       PCmw->Drc = PMmw->Drc / vol_w;

       // partitioning between sorbed and dissolved in mixing layer
       // L = m * m * m * [-] * 1000
       vol_w = zm->Drc * DX->Drc * SoilWidthDX->Drc * Theta_mix->Drc * 1000;
       // kg = m * m * m * kg m-3
       mass_s = zm->Drc * DX->Drc * SoilWidthDX->Drc * rho;
       // for now the assumption is made that the resulting unit of
       // Kr * (Kd * PCmw->Drc - PCms->Drc) is mg * kg-1 * sec-1
       // this holds if 1L water = 1kg

       // positive adds to sorbed, negative to dissolved.
       // mg = mg kg-1 sec-1 *  sec * kg
       // based on PhD lefrancq should always be mass_s
       mda_ex = kr * (Kd * PCmw->Drc - PCms->Drc) * _dt
                    * mass_s;

       // calculate equilibrium mass division
       mda_tot = PMmw->Drc + PMms->Drc;
       eql_diss = mda_tot / (1 + (Kd / vol_w * mass_s));
       eql_ads = mda_tot - eql_diss;
       // mda_ex can not be larger than m_diff
       m_diff = eql_ads - PMms->Drc;
       mda_ex = std::abs(mda_ex) > std::abs(m_diff) ? m_diff : mda_ex;

       //update masses
       PMmw->Drc = std::max(0.0, PMmw->Drc - mda_ex);
       // mg = mg + mg
       PMms->Drc = std::max(0.0, PMms->Drc + mda_ex);
       // update PCmw beforelateral transport
       PCmw->Drc = PMmw->Drc / vol_w;
       PCms->Drc = PMms->Drc / mass_s;

//----- OBSOLETE -------------------------------------------------
       // 2023-04-05 we exclude percolation transport since it does not improve
       // the model, and conceptually is not clear (MC)

       //percolation
//       PMperc->Drc = 0.0; //does not need to be a map...
//       if (InfilVol->Drc < 1e-6) {
//          // PMperc->Drc = PesticidePercolation(Perc->Drc, SoilDepth1->Drc,
//          //              Lw->Drc, zm->Drc, DX->Drc, SoilWidthDX->Drc, PCmw->Drc);
//           Theta_mix->Drc = Thetaeff->Drc; //percolation related theta for mixing layer
//       }
// --------------------------------------------------------------

       // mixing layer -- runoff water exchange
       double mwrm_ex {0.0};   //mg
       double A_mix {0.0};    // surface area of mixing transfer
       // if the water volume in a cell is too small, we cannot assume a film
       // over the full surface of the cell. This would overestimate mixing
       // mass transfer. When water height is smaller than 0.1 mm we assume the
       // surface area for mass transfer decreases.
       PCrw->Drc = PMrw->Drc / WaterVolall->Drc;
       if (WaterVolall->Drc > 0.0) {
           if (WH->Drc < 1e-4) {
               A_mix = WaterVolall->Drc / 1e-4;
           } else A_mix = DX->Drc * SoilWidthDX->Drc;
       // positive adds to runoff, negative to mixing layer.
       // mg = ((m sec-1 (mg m-3)) m2 * sec
       mwrm_ex = (Kfilm * (PCmw->Drc - PCrw->Drc) * 1000) * A_mix * _dt;

       double c_eql {0.0};
       double eql_mw {0.0};
       double vol_mw {0.0};
       // equilibrium check
       // calculate equilibrium mass division
       vol_mw = zm->Drc * Theta_mix->Drc * DX->Drc * SoilWidthDX->Drc * 1000;
       c_eql = (PMmw->Drc + PMrw->Drc) / (vol_mw + WaterVolall->Drc);
       eql_mw = c_eql * vol_mw; // mass in mixing layer at equilibrium
       // mwrm_ex can not be larger than m_diff
       m_diff = eql_mw - PMmw->Drc;
       mwrm_ex = std::abs(mwrm_ex) > std::abs(m_diff) ? m_diff : mwrm_ex;
       }
       // mass balance
       mwrm_ex > 0 ? pmwdet->Drc += mwrm_ex : pmwdep->Drc += mwrm_ex;

       PMmw->Drc = std::max(0.0, PMmw->Drc - mwrm_ex);
       PMrw->Drc = std::max(0.0, PMrw->Drc + mwrm_ex);
    }}
}

//---------------------------------------------------------------------------
/**
* @fn void TWorld::PesticideFlow1D(void)
* @brief Do all the calculations for pesticide dynamics between cells
* This includes:
*   call functions for adsorbed and dissolved dynamics
*   update all masses etc.
*/

void TWorld::PesticideFlow1D(void) {

    //double Kfilm = KfilmPest; // m sec-1
    double rho = rhoPest;     //kg m-3

    //runoff
    KinematicPestDissolved(crlinkedldd_, LDD, Qn, PQrw, DX, Alpha, Q, Qpw,
                        PMmw);

    //erosion
    if(SwitchErosion){
        KinematicPestAdsorbed(crlinkedldd_, LDD, Qsn, PQrs, DX, Alpha, SedMassIn,
                              Qs, Qps, PMms);
    }
    // calculate new concentration
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
    double volmw {0.0};         // L - volume of water in mixing layer
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
* @fn double TWorld::KinematicPestDissolved(double perc, double soildep,
*               double lw, double zm, double dx, double swdx, double pcmw)
* @brief explicit kinematic wave for dissolved pesticides
*/

void TWorld::KinematicPestDissolved(QVector <LDD_COORIN> _crlinked_,
               cTMap *_LDD, cTMap *_Qn, cTMap *_Qpwn, cTMap *_DX,
               cTMap *_Alpha, cTMap *_Q, cTMap *_Qpw, cTMap *_PMW)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _Qpwn->Drc = 0;
        QpinKW->Drc = 0;
    }}

// loop over ldd
for(long i_ =  0; i_ < _crlinked_.size(); i_++)
{
    int r = _crlinked_[i_].r;
    int c = _crlinked_[i_].c;

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
                    Qpin += _Qpwn->Drcr;
                }
            }
        }
    }
    QpinKW->Drc = Qpin;

    if (Qn->Drc + QinKW->Drc >= MIN_FLUX) { // more than 1 ml - what is best definition of runoff?
        // calculate concentration for new outflux
        PCrw->Drc = PMrw->Drc / (WaterVolall->Drc * 1000); // use watervollall and not watervolin for concentration

        _Qpw->Drc = _Q->Drc * 1000 * PCrw->Drc;
        // use explicit backwards method from Chow
        _Qpwn->Drc = ChowSubstance(_Qn->Drc, QinKW->Drc, _Q->Drc, QpinKW->Drc, _Qpw->Drc,
                                    _Alpha->Drc, _DX->Drc, _dt); //mg/sec
        _Qpwn->Drc = std::min(_Qpwn->Drc, QpinKW->Drc + _PMW->Drc / _dt);
       } //runoff occurs
    //substract discharge
    //mg = mg - (mg sec-1 * sec)
    PMrw->Drc = std::max(0.0, PMrw->Drc - (_Qpwn->Drc * _dt)
                                  + (QpinKW->Drc * _dt));
    }//end ldd loop
}

//---------------------------------------------------------------------------
/**
* @fn double TWorld::KinematicPestAdsorbed(double perc, double soildep,
*               double lw, double zm, double dx, double swdx, double pcmw)
* @brief Calculate adsorbed pesticide mass transported with runoff sediment
*/

void TWorld::KinematicPestAdsorbed(QVector <LDD_COORIN> _crlinked_,
                             cTMap *_LDD, cTMap *_Qn, cTMap *_Qpsn, cTMap *_DX,
                             cTMap *_Alpha, cTMap *_Sed, cTMap *_Q, cTMap *_Qps,
                                   cTMap *_PMS)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _Qpsn->Drc = 0;
        SpinKW->Drc = 0;
    }}

for(long i_ =  0; i_ < _crlinked_.size(); i_++)
{
    int r = _crlinked_[i_].r;
    int c = _crlinked_[i_].c;

    double Spin {0.0}; //mg sec-1

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
                    Spin += _Qpsn->Drcr;
                }
            }
        }
    }
    SpinKW->Drc = Spin;

    if (_Sed->Drc > 0 | SinKW->Drc > 0.0) { //
        if (Qn->Drc >= MIN_FLUX) {
//        // - simple extrapolation
//        double totpests = std::max(0.0, PMrs->Drc + (SpinKW->Drc * _dt));
//        double totsed = _Sed->Drc + (SinKW->Drc * _dt);
//        _Qpsn->Drc = std::min(totpests/_dt,
//                                  _Qsn->Drc * (totpests / totsed));

        // Chow applied for concentration of adsorbed in water!
        // Instead of the adsorbed concentration of pesticides in the sediment
        // mg/kg we use the concentration in suspended sediment multiplied by
        // the suspended sediment concentration - resulting in the adsorbed
        // pesticide concentration in the runoff water. This is suitable to be
        // solved with the explicit Chow equation. And takes flow speed into
        // acount when reditributing the adsorbed pesticide. The 'simple extrapolation
        // above does not do that and causes extreme concentration peaks at the
        // rising limb of the discharge.
        // mg sec-1 = m3 sec -1 * (mg m-3)
        _Qps->Drc = Q->Drc * (PMrs->Drc / WaterVolall->Drc);
        // use explicit backwards method from Chow
        _Qpsn->Drc = ChowSubstance(_Qn->Drc, QinKW->Drc, _Q->Drc, SpinKW->Drc, _Qps->Drc,
                                 _Alpha->Drc, _DX->Drc, _dt); //mg/sec
        _Qpsn->Drc = std::min(_Qpsn->Drc, SpinKW->Drc + _PMS->Drc / _dt);
        }
        } // erosion occurs
    // can move outside ldd loop to parralel section
    // mg = mg sec-1 * sec
    PMrs->Drc = std::max(0.0, PMrs->Drc - (_Qpsn->Drc * _dt)
                                  + (SpinKW->Drc * _dt));
    PCrs->Drc = Sed->Drc > 1e-6 ? PMrs->Drc / Sed->Drc : 0.0; // divide by Sed after kin wave
    // 0,001 g
    }// end ldd loop
}

//---------------------------------------------------------------------------
/**
* @fn double TWorld::PesticideSplashDetachment(;
* @brief Calculate adsorbed pesticide mass added to flow by splash erosion
*/

void TWorld::PesticideSplashDetachment() {
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
         double msoil_ex {0.0};  // mass exchange between mixing layer and deeper soil
        // add mass to pesticide in flow
        PMsplash->Drc = DETSplash->Drc * PCms->Drc;
        PMrs->Drc += PMsplash->Drc;
        PCrs->Drc = Sed->Drc > 1e-6 ? PMrs->Drc / Sed->Drc : 0.0; // if Sed > 0.001 gr

        // update mass and concentration in mixing layer
        msoil_ex = DETSplash->Drc * PCs->Drc;
        PMms->Drc = PMms->Drc - PMsplash->Drc + msoil_ex;
        PCms->Drc = PMms->Drc / (zm->Drc * DX->Drc * SoilWidthDX->Drc * rhoPest);
        // adjust lower soil layer
        PMsoil->Drc = std::max(0.0, PMsoil->Drc - msoil_ex);
        PCs->Drc = PMsoil->Drc / (zs->Drc * DX->Drc * SoilWidthDX->Drc * rhoPest);

        //mass balance
        pmsdet->Drc += PMsplash->Drc;
    }}
}




//---------------------------------------------------------------------------
/**
 * @fn double TWorld::ChowSubstance(double Qj1i1, double Qj1i, double Qji1,double Pj1i, double Pji1, double alpha, double dx)
 * @brief Explicit backward calculation of substance in water outflux from a cell
 *
 * Calculation of substance outflux from a cell based on a explicit solution of the time/space matrix,
 * j = time and i = place: j1i1 is the new output, j1i is the new flux at the upstream 'entrance' flowing into the gridcell
 *
 * @param Qj1i1 : result kin wave for this cell ( Qj+1,i+1 )
 * @param Qj1i : sum of all upstreamwater from kin wave ( Qj+1,i ),
 * @param Qji1 : incoming Q for kinematic wave (t=j) in this cell, map Q in LISEM (Qj,i+1)
 * @param Pj1i : sum of all upstream substance (Pj+1,i)
 * @param Pji1 : incoming substance for kinematic wave (t=j) in this cell, map Qpw in LISEM (Si,j+1)
 * @param alpha : alpha calculated in LISEM from before kinematic wave
 * @param dt : timestep
 * @param dx : length of the cell, corrected for slope (DX map in LISEM)
 * @return substance outflow in next timestep
 *
 */
double TWorld::ChowSubstance(double Qj1i1, double Qj1i, double Qji1,double Pj1i,
                             double Pji1, double alpha, double dx, double dt)
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

//---------------------------------------------------------------------------
/**
* @fn double TWorld::PesticideDetachment(double rho);
* @brief Calculate mass exchange by erosion and deposition with soil
*/

void TWorld::PesticideFlowDetachment(double rho) {

  // mass exchange between mixing layer an suspended sediment
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L{
        double msoil_ex {0.0};  // mass exchange between mixing layer and deeper soil
        double msrm_ex {0.0};
        // For now only use SoilWidth in formulas. Check what is done with deposition on roads.
        // Can this be eroded after deposition or not?
        // option 1 - all deposition on roads add directly to sink
        // option 2 - deposition on roads can be eroded and added into the system...

        PMdep->Drc = 0.0;
        PMflow->Drc = 0.0;

        if (DEP->Drc < 0) {
            //deposition
            msoil_ex = DEP->Drc * PCms->Drc; // what happens with pesticides on roads??
            // mg = mg kg-1 * kg
            msrm_ex = (DEP->Drc/SedAfterSplash->Drc) * PMrs->Drc; // loss by deposition
            // no more transport than mass in cell domain
            if (PMrs->Drc + msrm_ex < 0) {
                msrm_ex = -PMrs->Drc;
            }
            PMdep->Drc = msrm_ex;
        } else if (DETFlow->Drc > 0){
            // erosion
            msoil_ex = DETFlow->Drc * PCs->Drc; //
            // mg = mg kg-1  kg
            msrm_ex = PCms->Drc * DETFlow->Drc; // added by erosion
            // no more transport than mass in cell domain
            if (PMms->Drc + msoil_ex < msrm_ex) {
                msrm_ex = PMms->Drc + msoil_ex;
            }
            PMflow->Drc = msrm_ex;
        }

        // mass balance
        DEP->Drc < 0 ? pmsdep->Drc += msrm_ex: pmsdet->Drc += msrm_ex;

        // pesticides in suspended sediment
        PMrs->Drc = std::max(0.0, PMrs->Drc + msrm_ex);
        PCrs->Drc = SedMassIn->Drc > 1e-6 ? PMrs->Drc / SedMassIn->Drc : 0.0; //

        // adjust mass lower soil layer for mass balance
        PMsoil->Drc = std::max(0.0, PMsoil->Drc - msoil_ex);
        PCs->Drc = PMsoil->Drc / (zs->Drc * DX->Drc * SoilWidthDX->Drc * rhoPest);

        // pesticides in mixing layer
        PMms->Drc = std::max(0.0, PMms->Drc - msrm_ex + msoil_ex);
        PCms->Drc = PMms->Drc / (zm->Drc * DX->Drc * SoilWidthDX->Drc * rhoPest);

    }}
}
