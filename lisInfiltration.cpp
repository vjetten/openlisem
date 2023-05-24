

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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/
/*!
  \file lisInfiltration.cpp
  \brief Simplified infiltraton processes: Green and Ampt, Smith and Parlanage, both 1 and 2 layer. SWATRE has separate files.

functions: \n
- void TWorld::InfilEffectiveKsat(void)
- void TWorld::InfilSwatre(cTMap *_WH)
- void TWorld::InfilMethods(cTMap * _Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull)
- double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, double *L1p, double *L2p, double *FFullp)
- void TWorld::Infiltration(void)
 */

#include <algorithm>
#include "lisemqt.h"
#include "global.h"
#include "model.h"
#include "operation.h"

//---------------------------------------------------------------------------
// Done outside timeloop, move inside when crusting is made dynamic!
void TWorld::InfilEffectiveKsat(bool first)
{
    // todo, move to datainit!
    if (first) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Ksat1->Drc *= _dt/3600000.0;
            if (SwitchTwoLayer)
                Ksat2->Drc *= _dt/3600000.0;
            if (SwitchInfilCrust)
                KsatCrust->Drc *= _dt/3600000.0;
            if (SwitchInfilCompact)
                KsatCompact->Drc *= _dt/3600000.0;
        }}
    }


    if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Ksateff->Drc = Ksat1->Drc;
            Poreeff->Drc = ThetaS1->Drc;

            // exponential crusting proces with cumulative rainfall
            if (SwitchInfilCrust) {
                //double KSc = Ksat1->Drc * (0.3+0.7*exp(-0.05*RainCum->Drc*1000));
                double ksatdiff = std::max(0.0,Ksat1->Drc - KsatCrust->Drc);
                double factor = RainCum->Drc > 0.01 ? exp(-0.05*(RainCum->Drc-0.01)*1000) : 1.0;
                double KSc = KsatCrust->Drc + ksatdiff * factor;
                // exponential decline until crust value, RainCum is in meters

                //Ksateff->Drc = (1-Cover->Drc) * KSc + Cover->Drc * Ksat1->Drc;
                Ksateff->Drc = KSc;
                // only on bare fraction of soil, depends on crop. We need basal cover! ???
                double porediff = std::max(0.0,ThetaS1->Drc - PoreCrust->Drc);
                Poreeff->Drc = PoreCrust->Drc + porediff * factor;
                        //ThetaS1->Drc*(1-CrustFraction->Drc) + PoreCrust->Drc*CrustFraction->Drc;
            }
            Thetaeff->Drc = std::max(0.025*Poreeff->Drc,ThetaI1->Drc);

            // affected surfaces
            if (SwitchInfilCompact) {
                Ksateff->Drc = Ksateff->Drc*(1-CompactFraction->Drc) + KsatCompact->Drc*CompactFraction->Drc;
                Poreeff->Drc = Poreeff->Drc*(1-CompactFraction->Drc) + PoreCompact->Drc*CompactFraction->Drc;
            }

            if (SwitchGrassStrip) {
                Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;
                Poreeff->Drc = ThetaS1->Drc*(1-GrassFraction->Drc) + PoreGrass->Drc*GrassFraction->Drc;
            }


            if (SwitchHouses) {
                Ksateff->Drc *= (1-HouseCover->Drc);
             //   Poreeff->Drc *= (1-HouseCover->Drc);
            }

            //these surfaces are excluded from infiltration so not necessary to adjust Ksat and Pore
//            // impermeable surfaces
//            if (SwitchHardsurface) {
//                Ksateff->Drc *= (1-HardSurface->Drc);
//             //   Poreeff->Drc *= (1-HardSurface->Drc);
//            }

//            if (SwitchRoadsystem) {
//                Ksateff->Drc *= (1-RoadWidthDX->Drc/_dx);
//             //   Poreeff->Drc *= (1-RoadWidthDX->Drc/_dx);
//            }

            Ksateff->Drc = std::max(0.0, Ksateff->Drc);
            Poreeff->Drc = std::max(0.3, Poreeff->Drc);
           // Thetaeff->Drc = std::min(1.0,Poreeff->Drc/ThetaS1->Drc) * ThetaI1->Drc;
           // tma->Drc =  Ksateff->Drc;
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

// this function is not used!
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

//---------------------------------------------------------------------------
// Infiltration by Green and Ampt,Smith and Parlange
// All the same except for calculation of the potential infiltration fpot
// 1 layer and 2 layers
/*!
\brief function to calculate potential and actula infiltration rate according to
Green and Ampt, or Smith and Parlange.

This function calculates the potential infiltration according to G&A or S&P \n
then calls IncreaseInfiltrationDepth to increase the wetting front.
*/
void TWorld::cell_InfilMethods(int r, int c)
{
    // default vars are first layer vars
    double Ks = Ksateff->Drc;  //in m
    double Psi = Psi1->Drc; // in m
    double fwh = 0;
    double fpot_ = 0;
    double fact_ = 0;
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SoilDepth2->Drc;

    if (FloodDomain->Drc == 0) {
        fwh = WH->Drc; //runoff in kinwave or dyn wave
    } else {
        fwh = hmx->Drc; // flood in kin wave
    }
    // select the appropriate domain water height for overpressure

    // only do infiltration on permeable soils
    if (SoilWidthDX->Drc > 0 && fwh > 0) {

        //calculate potential infiltration rate fpot
        if (SwitchTwoLayer || SwitchThreeLayer) {
            // if wetting front in second layer set those vars
            if (Lw->Drc > SoilDep1 && Lw->Drc < SoilDep2) {
                //weighed harmonic mean:
                //https://corporatefinanceinstitute.com/resources/data-science/harmonic-mean/
                // sum (weights) / sum (weight/variable)

                Ks = Lw->Drc/(SoilDep1/Ksateff->Drc+(Lw->Drc-SoilDep1)/Ksat2->Drc);
                // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
                Psi = Psi2->Drc; //in m
            }
        }

        if (InfilMethod == INFIL_GREENAMPT)
            fpot_ = Ks*(1.0+(Psi+fwh)/std::max(1e-3, Lw->Drc));
        else {
            // smith parlange, not really tested
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
        // actual infil in m, cannot have more infil than water on the surface, includes rainfall

        if (fact_ > 0) {
//            if (SwitchThreeLayer)
//                fact_ = IncreaseInfiltrationDepthNew3(fact_, r, c);
//            else
                if (SwitchTwoLayer)
                    fact_ = IncreaseInfiltrationDepthNew2(fact_, r, c);
                else
                    fact_ = IncreaseInfiltrationDepthNew1(fact_, r, c);
        }
        // adjust fact and increase Lw, for twolayer, impermeable etc

        if (fwh < fact_)
        {
            fact_ = fwh;
            fwh = 0;
        }
        else
            fwh -= fact_;

        if(FloodDomain->Drc == 0)
            WH->Drc = fwh;
        else
            hmx->Drc = fwh;
        // adjust the WH in the correct domain with new fact

        Fcum->Drc += fact_; // for Smith and Parlange
        // increase cumulative infil in m
       // fact->Drc = fact_;
        InfilVol->Drc = fact_* SoilWidthDX->Drc * DX->Drc;
        // calc infiltrated volume for mass balance
    } else {
       // fact->Drc = 0;
        InfilVol->Drc = 0;
    }

    // calc surplus infiltration (negative in m) for kin wave
    // no longer used
    FSurplus->Drc = 0;
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
double TWorld::IncreaseInfiltrationDepthNew1(double fact_in, int r, int c)
{
    double dtheta1 = std::max(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double L = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;

    // impermeable and L reached SD1, no more infil
    if (SwitchImpermeable && L > SoilDep1 - 0.001) {
        Lw->Drc = SoilDep1;
        return 0;
    }

    Lnew = L + fact_in/std::max(dtheta1,0.01);
    // increase wetting front
    space = (SoilDep1 - L)*dtheta1;
    if(Lnew > SoilDep1 || space < fact_in) {
        if (SwitchImpermeable)
            // if impermeable remaining space is infiltration
            fact_out = space;
        else
            fact_out = Perc->Drc;
        Lnew = SoilDep1;
    } else {
        fact_out = fact_in;
    }

    Lw->Drc = std::min(SoilDep1,std::max(0.0, Lnew));
    return std::max(0.0, fact_out);
}
//---------------------------------------------------------------------------
double TWorld::IncreaseInfiltrationDepthNew2(double fact_in, int r, int c)
{
    double dtheta1 = std::max(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double dtheta2 = std::max(0.0,ThetaS2->Drc-ThetaI2->Drc);
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;
    double L = Lw->Drc;
    double dfact2 = 0;
    bool passing = false;
    double space2 = 0;

    // profile is full
    if (SwitchImpermeable && L > SoilDep2 - 0.001) {
        Lw->Drc = SoilDep2;
        return 0;
    }

    if (SwitchGWflow) {
       if (L >= SoilDep1 && GWWH->Drc >= SoilDepth2init->Drc-HMIN) {
           Lw->Drc = SoilDep1;
           return 0;
       }
       // when GWWH fills osildep2 than soildep2 is 0 anyway
    }


    // L is in layer 1
    if (L <= SoilDep1) {
        Lnew = L + fact_in/std::max(0.01,dtheta1);
        space = (SoilDep1-L)*dtheta1;

        if(fact_in > space || Lnew > SoilDep1) {
            // water is moving into layer 2
            passing = true;
            dfact2 = fact_in - space;
            // remaining water for layer 2
        } else {
            // all remains SD1
            fact_out = fact_in;
        }
    }

    // L is in layer 2
    if (L > SoilDep1) {
        //L already in layer 2

        Lnew = L + fact_in/std::max(0.01,dtheta2);
        space2 = (SoilDep2-L)*dtheta2;

        if (Lnew > SoilDep2 || fact_in > space2) {
            if (SwitchImpermeable)
                fact_out = space2;
            else
                fact_out = Perc->Drc;

            Lnew = SoilDep2;
            // L at bottom
        } else {
            fact_out = fact_in;
            // everything fitted
        }
    }
    // Lnew is now soildep2 or the actual depth

    // L is moving from layer 1 into 2 in this timestep
    if (passing) {
        // second layer still at initial
        space2 = (SoilDep2-SoilDep1)*dtheta2;
        Lnew = SoilDep1 + dfact2/std::max(0.01,dtheta2);
        dfact2 = std::min(dfact2, space2);

        if (dtheta2 < 0.01 || Lnew > SoilDep2) {
            if (SwitchImpermeable)
                fact_out = space+space2;
            else
                fact_out = Perc->Drc;
            Lnew = SoilDep2;
        } else
            fact_out = fact_in; // everything fitted
    }

    Lw->Drc = std::min(SoilDep2,std::max(0.0, Lnew));
    return std::max(0.0,fact_out);
}
//---------------------------------------------------------------------------
// 3 layer infiltration! not used yet
double TWorld::IncreaseInfiltrationDepthNew3(double fact_in, int r, int c)
{
    double dtheta1 = std::max(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double dtheta2 = std::max(0.0,ThetaS2->Drc-ThetaI2->Drc);
    double dtheta3 = std::max(0.0,ThetaS3->Drc-ThetaI3->Drc);
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double SoilDep3 = SoilDepth3->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;
    double L = Lw->Drc;
    double dfact12 = 0;
    double dfact23 = 0;
    bool passing12 = false;
    bool passing23 = false;
    double space2 = 0;
    double space3 = 0;

    // profile is full
    if (SwitchImpermeable && L > SoilDep2 - 0.001) {
        Lw->Drc = SoilDep2;
        return 0;
    }

    // L is in layer 1
    if (L <= SoilDep1) {
        Lnew = L + fact_in/std::max(0.01,dtheta1);
        space = (SoilDep1-L)*dtheta1;

        if(fact_in > space || Lnew > SoilDep1) {
            // water is moving into layer 2
            passing12 = true;
            dfact12 = fact_in - space;
            // remaining water for layer 2
        } else {
            // all remains SD1
            fact_out = fact_in;
        }
    }

    // L is in layer 2
    if (L > SoilDep1 && L <= SoilDep2) {
        //L already in layer 2 but not in 3
        Lnew = L + fact_in/std::max(0.01,dtheta2);
        space2 = (SoilDep2-L)*dtheta2;

        if (fact_in > space2 || Lnew > SoilDep2) {
            passing23 = true;
            dfact23 = fact_in - space2;
        } else {
            // all remains SD2
            fact_out = fact_in;
        }
    }

    // L is in layer 3
    if (L > SoilDep2 && L <= SoilDep3) {
        //L already in layer 2 but not in 3
        Lnew = L + fact_in/std::max(0.01,dtheta3);
        space3 = (SoilDep3-L)*dtheta3;

        if (fact_in > space3 || Lnew > SoilDep3) {
            if (SwitchImpermeable)
                fact_out = space3;
            else
                fact_out = Perc->Drc;

            Lnew = SoilDep3;
            // L at bottom
        } else {
            // all remains SD3
            fact_out = fact_in;
        }
    }
    // Lnew is now soildep3 or the actual depth

    // L is moving from layer 1 into 2 in this timestep
    if (passing12) {
        // second layer still at initial
        space2 = (SoilDep2-SoilDep1)*dtheta2;
        Lnew = SoilDep1 + dfact12/std::max(0.01,dtheta2);
        dfact12 = std::min(dfact12, space2);

        if (dtheta2 < 0.01 || Lnew > SoilDep2) {
            passing23 = true;
            dfact23 = fact_in - space2;
            // also does not fit in SD2, passing to SD3
        } else {
            fact_out = fact_in;
            // everything fitted in SD2
        }
    }

    // L is moving from layer 2 into 3 in this timestep
    if (passing23) {
        // second layer still at initial
        space3 = (SoilDep3-SoilDep2)*dtheta3;
        Lnew = SoilDep2 + dfact23/std::max(0.01,dtheta3);
        dfact23 = std::min(dfact23, space3);

        if (dtheta3 < 0.01 || Lnew > SoilDep2) {
            if (SwitchImpermeable)
                fact_out = space+space3;
            else
                fact_out = Perc->Drc;
            Lnew = SoilDep3;
        } else
            fact_out = fact_in; // everything fitted
    }

    Lw->Drc = std::min(SoilDep3,std::max(0.0, Lnew));
    return std::max(0.0,fact_out);

}

//---------------------------------------------------------------------------
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

        // not done, experimental
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
