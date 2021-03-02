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
            Thetaeff->Drc = ThetaI1->Drc;

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

            // impermeable surfaces
//            if (SwitchHardsurface) {
//                Ksateff->Drc *= (1-HardSurface->Drc);
//                Poreeff->Drc *= (1-HardSurface->Drc);
//            }

            if (SwitchHouses) {
                Ksateff->Drc *= (1-HouseCover->Drc);
                Poreeff->Drc *= (1-HouseCover->Drc);
            }

//            if (SwitchRoadsystem) {
//                Ksateff->Drc *= (1-RoadWidthDX->Drc/_dx);
//                Poreeff->Drc *= (1-RoadWidthDX->Drc/_dx);
//            }
//DO NOT INCLUDE ROADS AND HARDSURF HERE BECUAE SOILWIDTH ALREADY EXCLUDES THOSE (but not houses!)

            Ksateff->Drc = std::max(0.0, Ksateff->Drc);
            Poreeff->Drc = std::max(0.35, Poreeff->Drc);
            Thetaeff->Drc = std::min(1.0,Poreeff->Drc/ThetaS1->Drc) * ThetaI1->Drc;
                    //std::min(Thetaeff->Drc, Poreeff->Drc);
            // prob unneccessary

            Ksateff->Drc *= ksatCalibration;
            // apply runfile/iface calibration factor

            if(SwitchTwoLayer)
                bca->Drc = 5.55*qPow(Ksat2->Drc,-0.114);
            else
                bca->Drc = 5.55*qPow(Ksateff->Drc,-0.114);
            // percolation coefficient
        }}
    }
report(*Ksateff,"ksateff1.map");
report(*Poreeff,"poreeff1.map");
report(*Thetaeff,"thetaeff1.map");
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
    double dtheta1 = Poreeff->Drc-Thetaeff->Drc; // space in the top layer
    double L = Lw->Drc;
    double fact_out = 0;
    bool passing = false;
    double SoilDep1 = SoilDepth1->Drc;

    if (SwitchTwoLayer && L < SoilDepth2->Drc) {

        double SoilDep2 = SoilDepth2->Drc;
        double dtheta2 = ThetaS2->Drc-ThetaI2->Drc;
        double dfact1 = 0, dfact2 = 0;

        if (SwitchImpermeable && L > SoilDep2 - 1e-6) {
            Lw->Drc = SoilDep2;
            return 0;
        }

        if (L < SoilDep1) {
            // still in first layer
            L = L + fact_in/std::max(0.01,dtheta1);
//if (r == 258 && c == 630)
//    qDebug() << "1st layer" << r << c << L;

            if (L > SoilDep1) {
                // moving into second layer
                dfact1 = (SoilDep1 - Lw->Drc) * dtheta1;
                dfact2 = fact_in - dfact1; // remaining going into layer 2
                passing = true;
            } else {
                fact_out = fact_in;
            }
        } else {
            // already in 2nd layer
            L = L + fact_in/std::max(0.01,dtheta2);
//if (r == 258 && c == 630)
//qDebug() << "2nd layer" << r << c << L;
            if (L > SoilDep2) {
                fact_out = (SoilDep2 - Lw->Drc) * dtheta2;
                L = SoilDep2;
            } else {
                fact_out = fact_in;
            }
        }       

        // moving from layer 1 to 2
        if (passing) {
//if (r == 258 && c == 630)
//qDebug() << "passing" << r << c << L;
            L = SoilDep1 + dfact2/std::max(0.01,dtheta2); // increase L with remaining fact

            if (L > SoilDep2) {
                fact_out = dfact1 + (SoilDep2 - Lw->Drc) * dtheta2;
                L = SoilDep2;
            }
        } else {
            fact_out = fact_in;
        }
//if (r == 258 && c == 630)
//qDebug() << "return" << r << c << L << fact_in << fact_out << SoilDep1 << SoilDep2;
        Lw->Drc = L;
        return fact_out;

    } else {
        //single layer
        if (SwitchImpermeable && L > SoilDep1 - 1e-6) {
            Lw->Drc = SoilDep1;
            return 0;
        }

        if(L < SoilDep1) {
            // not full
            L = L + fact_in/std::max(0.01,dtheta1); // increase wetting front
            if (L > SoilDep1) {
                fact_out = (SoilDep1 - Lw->Drc) * dtheta1;
                L = SoilDep1;
            } else {
                fact_out = fact_in;
            }
        }
        Lw->Drc = L;
        return fact_out;
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

//    #pragma omp parallel for num_threads(userCores)
//    FOR_ROW_COL_MV_L {
        // default vars are first layer vars
        double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
        double Psi = Psi1->Drc/100; // in m
        double fwh = 0;
        double SW = 0;
        double fpot = 0;
        double fact = 0;// = fact->Drc;
        double SoilDep1 = SoilDepth1->Drc;

        if (FloodDomain->Drc == 0) {
            fwh = WH->Drc; //runoff
            SW = SoilWidthDX->Drc; //runoff
            // hard surf are already in ksaeff, so infil is affected
        } else {
            fwh = hmx->Drc; // flood
            SW = ChannelAdj->Drc;
        }
        // select the appropriate domain water height for overpressure

        //calculate potential insiltration rate fpot
        if (SwitchTwoLayer ) {

            // if wetting front in second layer set those vars
            if (Lw->Drc > SoilDep1)
            {
                Ks = std::min(Ks, Ksat2->Drc*_dt/3600000.0);
                // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
                Psi = Psi2->Drc/100;
            }
        }

        if (InfilMethod == INFIL_GREENAMPT)
            fpot = Ks*(1.0+(Psi+fwh)/std::max(1e-4, Lw->Drc));
        else {
            double space = SwitchTwoLayer ? std::max(ThetaS2->Drc-ThetaI2->Drc, 0.0) :
                                            std::max(Poreeff->Drc-Thetaeff->Drc, 0.0);
            double B = (fwh + Psi)*space;
            if (B > 0.01) {
                fpot = Ks*exp(Fcum->Drc/B)/(exp(Fcum->Drc/B)-1);
            } else
                fpot = Ks;
        }

        fact = std::min(fpot, fwh);
        if (fact < 1e-10) fact = 0;
        // actual infil in m, cannot have more infil than water on the surface

        if (fact > 0)
            fact = IncreaseInfiltrationDepthNew(fact, r, c);
        // adjust fact and increase Lw, for twolayer, impermeable etc

        if (fwh < fact)
        {
            fact = fwh;
            fwh = 0;
        }
        else
            fwh -= fact;

        // adjust the WH in the correct domain with new fact
        if(FloodDomain->Drc == 0)
            WH->Drc = fwh;
        else
            hmx->Drc = fwh;

        Fcum->Drc += fact; // for Smith and Parlange

        // increase cumulative infil in m
        InfilVol->Drc = fact*SW*DX->Drc;
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

            FSurplus->Drc = -1.0 * std::min(space, fact);//std::max(0.0, fpot_-fact_));
            // negative and smallest of space or fpot-fact ???
        }
 //   }}
}
//---------------------------------------------------------------------------

/*!
 \brief Calculates changes in soilwater with percolation from the bottom of the profile.

  Calculates changes in soilwater with percolation from the bottom of the profile, \n
  resulting in the soil becoming dryer. Based on BrooksCorey type of percolation: \n
  percolation = ksat*(theta/pore)*bca, where bca = 5.55*qPow(Ksat2->Drc,-0.114); \n
  This is completely undocumented. The soil is either impermeable or has percolation. \n
*/

void TWorld::SoilWater()
{
    if (InfilMethod == INFIL_SWATRE || InfilMethod == INFIL_NONE)
        return;
    if (SwitchImpermeable)
        return;

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Percolation, dL, pore, theta, thetar, theta_E, Ks;

        Percolation = 0;

        if(SwitchTwoLayer) {
            thetar = 0.025 * ThetaS2->Drc;
            pore = ThetaS2->Drc;
            theta = ThetaI2->Drc;
            Ks = Ksat2->Drc*_dt/3600000.0;

            if(theta > thetar) {
                theta_E = (theta-thetar)/(pore-thetar);
                Percolation = Ks* pow(theta_E, bca->Drc);
                // percolation in m

                if (Lw->Drc > SoilDepth1->Drc)
                    dL = SoilDepth2->Drc - Lw->Drc;
                else
                    dL = SoilDepth2->Drc - SoilDepth1->Drc;
                // if Wet Fr still in first layer percolation only make 2nd drier

                double moisture = dL*(theta-thetar);

                if (moisture > Percolation) {
                    // decrease thetaeff because of percolation
                    moisture -= Percolation;
                    theta = moisture/dL+thetar;
                } else {
                    // wetting front = soildepth1, dL = 0, moisture = 0
                    // assume theta goes back to 0.7 pore and decrease the wetting fornt
                    theta = 0.7*(pore - thetar);
                    Lw->Drc -= std::max(0.0, Percolation/(pore - theta));
                }
                ThetaI2->Drc = theta;
            }
        } else {
            // one layer
            thetar = 0.025 * Poreeff->Drc;
            double pore = Poreeff->Drc;
            double theta = Thetaeff->Drc;

            if(theta > thetar) {
                theta_E = (theta-thetar)/(pore-thetar);
                Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca->Drc);
            }

            if (Percolation > 0) {
                dL = std::max(0.0, SoilDepth1->Drc - Lw->Drc);
                double moisture = dL*(theta-thetar);
                if (moisture > Percolation) {
                    // wetting front has not reached bottom, make soil drier
                    // decrease thetaeff because of percolation
                    moisture -= Percolation;
                    theta = moisture/dL+thetar;
                } else {
                    // wetting front = soildepth1, dL = 0, moisture = 0
                    // assume tehta goes back to half pore and decrease the wetting fornt
                    theta = 0.7*(pore - thetar);
                    Lw->Drc -= std::max(0.0, Percolation/(pore - theta));
                }
                Thetaeff->Drc = theta;
            }
        }

        if (Percolation > 0) {
            double moisture = dL*(theta-thetar);
            if (moisture > Percolation) {
                // wetting front has not reached bottom, make soil drier
                // decrease thetaeff because of percolation
                moisture -= Percolation;
                theta = moisture/dL+thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume tehta goes back to half pore and decrease the wetting fornt
                theta = 0.7*(pore - thetar);
                Lw->Drc -= std::max(0.0, Percolation/(pore - theta));
            }
        }
        Perc->Drc = Percolation;
    }}
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

