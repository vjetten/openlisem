/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/
/*!
  \file lisInfiltration.cpp
  \brief Simplified infiltraton processes: Green and Ampt, Smith and Parlanage, both 1 and 2 layer. SWATRE has separate files.

functions: \n
- void TWorld::InfilEffectiveKsat(void)
- void TWorld::InfilSwatre(cTMap *_WH)
- void TWorld::InfilMorelSeytoux1(cTMap *_WH) Not working yet!
- void TWorld::InfilMethods(cTMap * _Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull)
- double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, double *L1p, double *L2p, double *FFullp)
- void TWorld::Infiltration(void)
- void TWorld::InfiltrationFloodNew(void)
- void TWorld::SoilWater()
 */

#include <algorithm>
#include "lisemqt.h"
#include "global.h"
#include "model.h"
#include "operation.h"



#define TINY 0.0001

//---------------------------------------------------------------------------
// Done outside timeloop, move inside when crusting is made dynamic!
void TWorld::InfilEffectiveKsat(void)
{
    if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Ksateff->Drc = Ksat1->Drc;
            Poreeff->Drc = ThetaS1->Drc;
            Thetaeff->Drc = std::max(0.025*Poreeff->Drc,ThetaI1->Drc);

            // affected surfaces
            if (SwitchInfilCompact) {
                Ksateff->Drc = Ksateff->Drc*(1-CompactFraction->Drc) + KsatCompact->Drc*CompactFraction->Drc;
                Poreeff->Drc = Poreeff->Drc*(1-CompactFraction->Drc) + PoreCompact->Drc*CompactFraction->Drc;
            }

            if (SwitchInfilCrust) {
                Ksateff->Drc = Ksateff->Drc*(1-CrustFraction->Drc) + KsatCrust->Drc*CrustFraction->Drc;
                Poreeff->Drc = ThetaS1->Drc*(1-CrustFraction->Drc) + PoreCrust->Drc*CrustFraction->Drc;
            }

            if (SwitchGrassStrip) {
                Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;
                Poreeff->Drc = ThetaS1->Drc*(1-GrassFraction->Drc) + PoreGrass->Drc*GrassFraction->Drc;
            }
            if (SwitchHouses) {
                Ksateff->Drc *= (1-HouseCover->Drc);
                Poreeff->Drc *= (1-HouseCover->Drc);
            }

            // impermeable surfaces
//            if (SwitchHardsurface) {
//                Ksateff->Drc *= (1-HardSurface->Drc);
//                Poreeff->Drc *= (1-HardSurface->Drc);
//            }

//            if (SwitchRoadsystem) {
//                Ksateff->Drc *= (1-RoadWidthDX->Drc/_dx);
//                Poreeff->Drc *= (1-RoadWidthDX->Drc/_dx);
//            }
//DO NOT INCLUDE ROADS AND HARDSURF HERE BECUAE SOILWIDTH ALREADY EXCLUDES THOSE (but not houses!)

            Ksateff->Drc = std::max(0.0, Ksateff->Drc);
            Poreeff->Drc = std::max(0.3, Poreeff->Drc);
            Thetaeff->Drc = std::min(1.0,Poreeff->Drc/ThetaS1->Drc) * ThetaI1->Drc;
//            bca1->Drc = 5.55*qPow(Ksateff->Drc,-0.114);
            Ksateff->Drc *= _dt/3600000;
            if(SwitchTwoLayer) {
//                bca2->Drc = 5.55*qPow(Ksat2->Drc,-0.114);
                Ksat2->Drc *= _dt/3600000;
            }
            // percolation coefficient
        }}
    }
}
//---------------------------------------------------------------------------
/*!
 \brief Main infiltration function, calls infiltration types (SWATRE, Green and Ampt,
  Smith and Parlange, Ksat subtraction. Calculates effective Ksat based on different
  surface types (crust, compaction).

  Main infiltration function that calculates\n
  - Use ksateff which accounts for different surface types: grass strips, compaction, crusting, roads, hard surface
  - do SWATRE or one of the other methods, SWATRE is a different set of functions
  - call one of the infiltration functions for the actual infiltration rate
  - calc infiltration surplus for the kinematic wave
  - increase of infiltration depth/wetting front, same function for each infiltration model: L1, L2, Fcum
  - decrease of surface water layer WH and calculate infiltration volume\n
  */
void TWorld::Infiltration()
{    
    //NOTE fact and fpot have a unit of m (not m/s)
    if (InfilMethod == INFIL_SWATRE) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            cell_InfilSwatre(r, c);
        }}
    }
    else
    if (InfilMethod != INFIL_NONE) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            cell_InfilMethods(r, c);
        }}
    }
}


void TWorld::cell_InfilSwatre(int r, int c)
{

    if (FloodDomain->Drc == 0)
        tm->Drc = WH->Drc;
    else
        tm->Drc = hmx->Drc;

    WHbef->Drc = tm->Drc;

    SwatreStep(op.runstep, r,c, SwatreSoilModel, tm, fpot, TileDrainSoil, thetaTop);
    // WH and fpot done in swatrestep, for normal surface swatre should be done in all cells

    fact->Drc = (WHbef->Drc - tm->Drc);
    // actual infil is dif between WH before and after

    if (FloodDomain->Drc == 0)
        WH->Drc = tm->Drc;
    else
        hmx->Drc = tm->Drc;

    if (CrustFraction->Drc > 0) {
        tm->Drc = WHbef->Drc;
        tma->Drc = 0;  // gpot for crusted
        tmb->Drc = 0;
        tmc->Drc = 0;  //thetatop
        tmd->Drc = 0;

        SwatreStep(op.runstep, r, c, SwatreSoilModelCrust, tm, tma, tmb, tmc);//, CrustFraction);
        // calculate crust SWATRE and get the soil moisture of the top node

        if (FloodDomain->Drc == 0)
            tmd->Drc = WH->Drc;
        else
            tmd->Drc = hmx->Drc;
        // water level on crusted areas

        tmd->Drc = tm->Drc*CrustFraction->Drc + tmd->Drc*(1-CrustFraction->Drc);
        // weighted average
        if (FloodDomain->Drc == 0)
            WH->Drc = tmd->Drc;
        else
            hmx->Drc = tmd->Drc;

        fact->Drc = (WHbef->Drc - tmc->Drc);
        fpot->Drc = tma->Drc*CrustFraction->Drc + fpot->Drc*(1-CrustFraction->Drc);
        thetaTop->Drc = tmc->Drc*CrustFraction->Drc + thetaTop->Drc*(1-CrustFraction->Drc);
    }

    if (SwitchInfilCompact)
    {
        tm->Drc = WHbef->Drc;
        tma->Drc = 0; // fpot
        tmb->Drc = 0; // tile drain
        tmc->Drc = 0; // theta top layer for repellency
        tmd->Drc = 0;

        SwatreStep(op.runstep, r, c, SwatreSoilModelCompact, tm, tma, tmb, tmc);//, CompactFraction);

        if (FloodDomain->Drc == 0)
            tmd->Drc = WH->Drc;
        else
            tmd->Drc = hmx->Drc;
        tmd->Drc = tm->Drc*CompactFraction->Drc + tmd->Drc*(1-CompactFraction->Drc);
        if (FloodDomain->Drc == 0)
            WH->Drc = tmd->Drc;
        else
            hmx->Drc = tmd->Drc;

        fact->Drc = (WHbef->Drc - tmd->Drc);
        fpot->Drc = tma->Drc*CompactFraction->Drc + fpot->Drc*(1-CompactFraction->Drc);
        thetaTop->Drc = tmc->Drc*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
    }

    if (SwitchGrassStrip)
    {
            tm->Drc = WHbef->Drc;//WHGrass->Drc;
            tma->Drc = 0;
            tmb->Drc = 0;
            tmc->Drc = 0;
            tmd->Drc = 0;

        SwatreStep(op.runstep, r,c, SwatreSoilModelGrass, tm, tma, tmb, tmc);//, GrassFraction);

        if (FloodDomain->Drc == 0)
            tmd->Drc = WH->Drc;
        else
            tmd->Drc = hmx->Drc;
        tmd->Drc = tm->Drc*GrassFraction->Drc + tmd->Drc*(1-GrassFraction->Drc);
        if (FloodDomain->Drc == 0)
            WH->Drc = tmd->Drc;
        else
            hmx->Drc = tmd->Drc;

        fact->Drc = (WHbef->Drc - tmd->Drc);
        fpot->Drc = tma->Drc*GrassFraction->Drc + fpot->Drc*(1-GrassFraction->Drc);
        thetaTop->Drc = tmc->Drc*GrassFraction->Drc + thetaTop->Drc*(1-GrassFraction->Drc);
    }

    if (SwitchWaterRepellency)
    {
        //      FOR_ROW_COL_MV
        //      {
        //         RepellencyFraction->Drc = 1 - 1/(waterRep_d+pow(waterRep_a, 100*(thetaTop->Drc-waterRep_b)));
        //         //         if (thetaTop->Drc < waterRep_c)
        //         //            RepellencyFraction->Drc = 0;//1.0;
        //      }
        //        thetaTop->report("thtop");
        //        RepellencyFraction->report("repelfr");
    }
}


//---------------------------------------------------------------------------
/// SWATRE infiltration, takes WH and calculateds new WH and infiltration surplus for kin wave
void TWorld::InfilSwatre()
{
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (FloodDomain->Drc == 0)
            tm->Drc = WH->Drc;
        else
            tm->Drc = hmx->Drc;

        WHbef->Drc = tm->Drc;

        SwatreStep(op.runstep, r,c, SwatreSoilModel, tm, fpot, TileDrainSoil, thetaTop);
        // WH and fpot done in swatrestep
        // for normal surface swatre should be done in all cells
        fact->Drc = (WHbef->Drc - tm->Drc);
        // actual; infil is dif between WH before and after
        if (FloodDomain->Drc == 0)
            WH->Drc = tm->Drc;
        else
            hmx->Drc = tm->Drc;

        if (CrustFraction->Drc > 0) {
            tm->Drc = WHbef->Drc;
            tma->Drc = 0;
            tmb->Drc = 0;
            tmc->Drc = 0;
            tmd->Drc = 0; // WH or hmx
            SwatreStep(op.runstep, r, c, SwatreSoilModelCrust, tm, tma, tmb, tmc);//, CrustFraction);
            // calculate crust SWATRE and get the soil moisture of the top node

            if (FloodDomain->Drc == 0)
                tmd->Drc = WH->Drc;
            else
                tmd->Drc = hmx->Drc;
            tmd->Drc = tm->Drc*CrustFraction->Drc + tmd->Drc*(1-CrustFraction->Drc);
            if (FloodDomain->Drc == 0)
                WH->Drc = tmd->Drc;
            else
                hmx->Drc = tmd->Drc;
            fact->Drc = (WHbef->Drc - tmc->Drc);
            fpot->Drc = tma->Drc*CrustFraction->Drc + fpot->Drc*(1-CrustFraction->Drc);
            thetaTop->Drc = tmc->Drc*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
        }
    }}

    //calculate a new crustfraction for water repellency
    // formula = f = 1/(1+1.2^(theta-30)), theta in %

    if (SwitchInfilCrust)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            tm->Drc = WHbef->Drc;
            tma->Drc = 0;
            tmb->Drc = 0;
            tmc->Drc = 0;
            tmd->Drc = 0; // WH or hmx
            if (CrustFraction->Drc > 0)
                SwatreStep(op.runstep, r, c, SwatreSoilModelCrust, tm, tma, tmb, tmc);//, CrustFraction);
            // calculate crust SWATRE and get the soil moisture of the top node
            // CrustFraction is cells > 0


        }}


        // calculate average cell values
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            if (FloodDomain->Drc == 0)
                tmd->Drc = WH->Drc;
            else
                tmd->Drc = hmx->Drc;
            tmd->Drc = tm->Drc*CrustFraction->Drc + tmd->Drc*(1-CrustFraction->Drc);
            if (FloodDomain->Drc == 0)
                WH->Drc = tmd->Drc;
            else
                hmx->Drc = tmd->Drc;
            fact->Drc = (WHbef->Drc - tmc->Drc);
            fpot->Drc = tma->Drc*CrustFraction->Drc + fpot->Drc*(1-CrustFraction->Drc);
            thetaTop->Drc = tmc->Drc*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
        }}
    }

    if (SwitchInfilCompact)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            tm->Drc = WHbef->Drc;
            tma->Drc = 0; // fpot
            tmb->Drc = 0; // tile drain
            tmc->Drc = 0; // theta top layer for repellency
            tmd->Drc = 0;
            if (CompactFraction->Drc > 0)
                SwatreStep(op.runstep, r, c, SwatreSoilModelCompact, tm, tma, tmb, tmc);//, CompactFraction);
        }}


        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            if (FloodDomain->Drc == 0)
                tmd->Drc = WH->Drc;
            else
                tmd->Drc = hmx->Drc;
            tmd->Drc = tm->Drc*CompactFraction->Drc + tmd->Drc*(1-CompactFraction->Drc);
            if (FloodDomain->Drc == 0)
                WH->Drc = tmd->Drc;
            else
                hmx->Drc = tmd->Drc;

            fact->Drc = (WHbef->Drc - tmd->Drc);
            fpot->Drc = tma->Drc*CompactFraction->Drc + fpot->Drc*(1-CompactFraction->Drc);
            thetaTop->Drc = tmc->Drc*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
        }}
    }

    if (SwitchGrassStrip)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            tm->Drc = WHbef->Drc;//WHGrass->Drc;
            tma->Drc = 0;
            tmb->Drc = 0;
            tmc->Drc = 0;
            tmd->Drc = 0;
            if (GrassFraction->Drc > 0)
                SwatreStep(op.runstep, r, c, SwatreSoilModelGrass, tm, tma, tmb, tmc);//, GrassFraction);
        }}

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            if (FloodDomain->Drc == 0)
                tmd->Drc = WH->Drc;
            else
                tmd->Drc = hmx->Drc;
            tmd->Drc = tm->Drc*GrassFraction->Drc + tmd->Drc*(1-GrassFraction->Drc);
            if (FloodDomain->Drc == 0)
                WH->Drc = tmd->Drc;
            else
                hmx->Drc = tmd->Drc;

            fact->Drc = (WHbef->Drc - tmd->Drc);
            fpot->Drc = tma->Drc*GrassFraction->Drc + fpot->Drc*(1-GrassFraction->Drc);
            thetaTop->Drc = tmc->Drc*GrassFraction->Drc + thetaTop->Drc*(1-GrassFraction->Drc);
        }}
    }

    if (SwitchWaterRepellency)
    {
        //      FOR_ROW_COL_MV
        //      {
        //         RepellencyFraction->Drc = 1 - 1/(waterRep_d+pow(waterRep_a, 100*(thetaTop->Drc-waterRep_b)));
        //         //         if (thetaTop->Drc < waterRep_c)
        //         //            RepellencyFraction->Drc = 0;//1.0;
        //      }
        //        thetaTop->report("thtop");
        //        RepellencyFraction->report("repelfr");

    }
}

//---------------------------------------------------------------------------
/*!
\brief function to increase wetting front and deal with 2nd layer and impermeable subsoil
 returns actual infiltration rate.

 this function is called form all infiltration functions except Swatre:\n
 - one layer or two layers
  - returns depth of the wetting front (Lw)\n
 - returns actual infiltration in mm, NOT rate in mm/h
*/

double TWorld::IncreaseInfiltrationDepthNew(double fact_in, int r, int c) //, double fact, double *L1p, double *L2p, double *FFullp)
{
    double dtheta1 = std::max(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double L = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double fact_out = 0;
    bool passing = false;

    if (SwitchTwoLayer) {

        double SoilDep2 = SoilDepth2->Drc;
        double dtheta2 = std::max(0.0,ThetaS2->Drc-ThetaI2->Drc);
        double dfact2 = 0;

        if (SwitchImpermeable && L > SoilDep2 - 0.001) {
            Lw->Drc = SoilDep2;
            return 0;
        }                
        if (SwitchImpermeable && dtheta2 < 0.001) {
            Lw->Drc = SoilDep2;
            return 0;
        }

        if (L < SoilDep1) {
            double space = (SoilDep1-L)*dtheta1;
            if(fact_in > (SoilDep1-L)*space) {
                passing = true;
                dfact2 = fact_in - space;
            } else {
                // still in SD1
                if (dtheta1 > 0.001)
                    L = L + fact_in/dtheta1;

                fact_out = fact_in;
            }
        } else {
            //L already in layer 2
            double space2 = (SoilDep2-L)*dtheta2;
            if (dtheta2 > 0.001)
                L = L + fact_in/dtheta2;
            else
                L = SoilDep2;

            if (L > SoilDep2-0.001) {
                fact_out = space2;
                L = SoilDep2;                
            } else {
                fact_out = fact_in;
                // everything fitted
            }

        }

        if (passing) {
            double space2 = (SoilDep2-SoilDep1)*dtheta2;
            dfact2 = std::min(dfact2, space2);
            if (dtheta2 > 0.001)
                L = SoilDep1 + dfact2/dtheta2; // increase L with remaining fact
            else
                L = SoilDep2;

            if (L > SoilDep2-0.001) {
                fact_out = space2;
                L = SoilDep2;
            } else
                fact_out = fact_in; // everything fitted
        }
        L = std::min(SoilDep2,std::max(0.0, L));
        Lw->Drc = L;
        return std::max(0.0,fact_out);

    } else {

        //===== single layer =====

        if (SwitchImpermeable && L > SoilDep1 - 0.001) {
            Lw->Drc = SoilDep1;
            return 0;
        }
        if (SwitchImpermeable && dtheta1 < 0.001) {
            Lw->Drc = SoilDep1;
            return 0;
        }

        if(L < SoilDep1-0.001) {
            // not full
            double space1 = (SoilDep1 - L)*dtheta1;
            if (dtheta1 > 0.001)
                L = L + fact_in/dtheta1; // increase wetting front
            else
                L = SoilDep1;

            if (L > SoilDep1-0.001) {
                fact_out = space1;
                L = SoilDep1;
            } else
                fact_out = fact_in;
        }

        L = std::min(SoilDep1,std::max(0.0, L));
        Lw->Drc = L;
        return std::max(0.0, fact_out);
    }

    return 0;
}

//---------------------------------------------------------------------------
// Infiltration by Green and Ampt,Smith and Parlange
// All the same except for calculation of the potential infiltration fpot
// 1 layer and 2 layers
/*!
\brief function to calculate potential and actula infiltration rate according to
Green and Mapt, Smith and Parlange or Ksat subtraction.

This function calculates the potential infiltration according to Ksat, G&A or S&P \n
then calls IncreaseInfiltrationDepth to increase the wetting front.
*/

void TWorld::cell_InfilMethods(int r, int c)
{
    // default vars are first layer vars
    double Ks = Ksateff->Drc;//*_dt/3600000.0;  //in m
    double Psi = Psi1->Drc/100; // in m
    double fwh = 0;
    double SW = 0;
    double fpot_ = 0;
    double fact_ = 0;// = fact->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    if (FloodDomain->Drc == 0) {
        fwh = WH->Drc; //runoff in kinwave or dyn wave
        SW = SoilWidthDX->Drc; //runoff in kinwave or dyn wave
    } else {
        fwh = hmx->Drc; // flood in kin wave
        SW = SoilWidthDX->Drc; // flood in kin wave
    }
    // select the appropriate domain water height for overpressure

    //calculate potential infiltration rate fpot
    if (SwitchTwoLayer ) {

        // if wetting front in second layer set those vars
        if (Lw->Drc > SoilDep1)
        {
            Ks = std::min(Ks, Ksat2->Drc);//*_dt/3600000.0);
            // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
            Psi = Psi2->Drc/100;
        }
    }

    if (InfilMethod == INFIL_GREENAMPT)
        fpot_ = Ks*(1.0+(Psi+fwh)/std::max(1e-4, Lw->Drc));
    else {
        double space = Poreeff->Drc-Thetaeff->Drc;
        if (Lw->Drc > SoilDepth1->Drc)
            space = ThetaS2->Drc-ThetaI2->Drc;
        double B = (fwh + Psi)*space;
        if (B > 0.01) {
            fpot_ = Ks*exp(Fcum->Drc/B)/(exp(Fcum->Drc/B)-1);
        } else
            fpot_ = Ks;
    }

    fact_ = std::min(fpot_, fwh);
    if (fact_ < 1e-10)
        fact_ = 0;
    // actual infil in m, cannot have more infil than water on the surface

    if (fact_ > 0)
        fact_ = IncreaseInfiltrationDepthNew(fact_, r, c);
    // adjust fact and increase Lw, for twolayer, impermeable etc

    if (fwh < fact_)
    {
        fact_ = fwh;
        fwh = 0;
    }
    else
        fwh -= fact_;

    // adjust the WH in the correct domain with new fact
    if(FloodDomain->Drc == 0)
        WH->Drc = fwh;
    else
        hmx->Drc = fwh;

    Fcum->Drc += fact_; // for Smith and Parlange
    fact->Drc = fact_;

    // increase cumulative infil in m
    InfilVol->Drc = fact_*SW*DX->Drc;
    // calc infiltrated volume for mass balance

    // calc surplus infiltration (negative in m) for kin wave
    if(SwitchKinematic2D != K2D_METHOD_DYN) {
        double space = 0;
        if (Lw->Drc < SoilDep1)
            space = (SoilDep1 - Lw->Drc)*(Poreeff->Drc-Thetaeff->Drc);
        if (SwitchTwoLayer) {
            if (Lw->Drc > SoilDep1 && Lw->Drc < SoilDepth2->Drc)
                space = (SoilDepth2->Drc - Lw->Drc)*(ThetaS2->Drc-ThetaI2->Drc);
        }

        FSurplus->Drc = -1.0 * std::min(space, fact_);//std::max(0.0, fpot_-fact_));
        // negative and smallest of space or fpot-fact ???
    }
}
//---------------------------------------------------------------------------
// NOT USED
void TWorld::infilInWave(cTMap *_h, double dt1)
{
    if (InfilMethod == INFIL_SWATRE || InfilMethod == INFIL_NONE)
        return;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(FFull->Drc==0) {
            double cdx = DX->Drc;
            double cdy = ChannelAdj->Drc;

            //calculate infiltration in time step
            double infil = -1.0*FSurplus->Drc*dt1/_dt;
            if (_h->Drc < infil)
                infil = _h->Drc;
            _h->Drc -= infil;
            _h->Drc = std::max(_h->Drc , 0.0);
            FSurplus->Drc += infil;//*SoilWidthDX->Drc/cdy;
            FSurplus->Drc = std::min(0.0, FSurplus->Drc);

            Fcum->Drc -= infil;//*SoilWidthDX->Drc/cdy; //VJ !!!

            //keep track of infiltration
            InfilVolKinWave->Drc = (infil*cdx*cdy);
        }
    }}
}
//---------------------------------------------------------------------------

