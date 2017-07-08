
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
  \file lisInfiltration.cpp
  \brief Simplified infiltraton processes: Green and Ampt, Smith and Parlanage, both 1 and 2 layer. SWATRE has separate files.

functions: \n
- void TWorld::InfilEffectiveKsat(void)
- void TWorld::InfilSwatre(cTMap *_WH)
- void TWorld::InfilMorelSeytoux1(cTMap *_WH) Not working yet!
- void TWorld::InfilMethods(cTMap * _Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull)
- double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFullp)
- void TWorld::Infiltration(void)
- void TWorld::InfiltrationFloodNew(void)
- void TWorld::SoilWater()
 */

#include <algorithm>
#include "model.h"
#include "operation.h"

//NOTE fact and fpot have a unit of m (not m/s)

#define tiny 0.0001

//1e-8

//---------------------------------------------------------------------------
void TWorld::InfilEffectiveKsat(void)
{
    if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
    {
        FOR_ROW_COL_MV
        {
            Ksateff->Drc = Ksat1->Drc;
            if (SwitchInfilCrust)
                Ksateff->Drc = Ksat1->Drc*(1-CrustFraction->Drc) + KsatCrust->Drc*CrustFraction->Drc;
            if (SwitchInfilCompact)
                Ksateff->Drc = Ksat1->Drc*(1-CompactFraction->Drc) + KsatCompact->Drc*CompactFraction->Drc;
            if (SwitchInfilCrust && SwitchInfilCompact)
                Ksateff->Drc = Ksat1->Drc*(1-CompactFraction->Drc-CrustFraction->Drc) +
                        KsatCrust->Drc*CrustFraction->Drc + KsatCompact->Drc*CompactFraction->Drc;

            if (SwitchHardsurface)
                Ksateff->Drc *= (1-HardSurface->Drc);

            if (SwitchHouses)
                Ksateff->Drc *= (1-HouseCover->Drc);

            if (RoadWidthDX->Drc > 0)
                Ksateff->Drc *= (1-RoadWidthDX->Drc/_dx);

            if (GrassFraction->Drc > 0)
                Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;

            Ksateff->Drc = std::max(0.0, Ksateff->Drc);
            // incase combining all fractions lead to less than zero

            Ksateff->Drc *= ksatCalibration;
            // apply runfile/iface calibration factor

            if (SwitchBuffers && !SwitchSedtrap && SwitchBuffersImpermeable)
                if(BufferID->Drc > 0)
                    Ksateff->Drc = 0;
        }
    }
}
//---------------------------------------------------------------------------
/// SWATRE infiltration, takes WH and calculateds new WH and infiltration surplus for kin wave
void TWorld::InfilSwatre(cTMap *_WH)
{
    copy(*WHbef, *_WH); // copy water height before infil

    fill(*tma, 1.0); // flag to indicate where the swatrestep model should be run

    //calculate a new crustfraction for water repellency
    // formula = f = 1/(1+1.2^(theta-30)), theta in %

    // for normal surface swatre should be done in all cells
    SwatreStep(SwatreSoilModel, _WH, fpot, TileDrainSoil, thetaTop, tma);
    // NOTE WH changes in SWATRE
    // tiledrainsoil is in m per timestep, if not switchtiles then contains 0

    // WH and fpot done in swatrestep
    FOR_ROW_COL_MV
            fact->Drc = (WHbef->Drc - _WH->Drc);
    // actual; infil is dif between WH before and after

    if (SwitchInfilCrust)
    {
        copy(*tm, *WHbef);
        fill(*tma, 0.0);
        fill(*tmb, 0.0);
        fill(*tmc, 0.0);

        SwatreStep(SwatreSoilModelCrust, tm, tma, tmb, thetaTop, CrustFraction);
        // calculate crust SWATRE and get the soil moisture of the top node
        // CrustFraction is cells > 0

        // calculate average cell values
        FOR_ROW_COL_MV
        {
            //tm = WH on crust and tma = fpot crust
            WH->Drc = tm->Drc*CrustFraction->Drc + WH->Drc*(1-CrustFraction->Drc);
            fact->Drc = (WHbef->Drc - WH->Drc);
            //fact->Drc = (WHbef->Drc - tm->Drc)*CrustFraction->Drc + fact->Drc*(1-CrustFraction->Drc);
            fpot->Drc = tma->Drc*CrustFraction->Drc + fpot->Drc*(1-CrustFraction->Drc);
        }
    }

    if (SwitchInfilCompact)
    {
        copy(*tm, *WHbef);
        fill(*tma, 0.0);
        fill(*tmb, 0.0);
        fill(*tmc, 0.0);

        SwatreStep(SwatreSoilModelCompact, tm, tma, tmb, tmc, CompactFraction);

        FOR_ROW_COL_MV
        {
            WH->Drc = tm->Drc*CompactFraction->Drc + WH->Drc*(1-CompactFraction->Drc);
            // tm is WH on compacted area
            fact->Drc = (WHbef->Drc - WH->Drc);
            //fact->Drc = (tm->Drc - WHbef->Drc)*CompactFraction->Drc + fact->Drc*(1-CompactFraction->Drc);
            fpot->Drc = tma->Drc*CompactFraction->Drc + fpot->Drc*(1-CompactFraction->Drc);

            thetaTop->Drc = tmc->Drc*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
        }
    }

    if (SwitchGrassStrip)
    {
        copy(*tm, *WHGrass);
        fill(*tmb, 0.0);
        fill(*tmc, 0.0);

        SwatreStep(SwatreSoilModelGrass, WHGrass, fpotgr, tmb, tmc, GrassFraction);

        FOR_ROW_COL_MV
        {
            factgr->Drc = (tm->Drc - WHGrass->Drc);
            thetaTop->Drc = tmc->Drc*GrassFraction->Drc + thetaTop->Drc*(1-GrassFraction->Drc);
        }
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
 - impermeable bottom or not, returns flag profile full
 - returns depth of the wetting front (L1 or L1+L2)\n
 - returns actual infiltration
*/
double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFullp)
{
    double dL1, dL2; // increase in wetting front layer 1 and 2 in m
    double L1, L2, FFull; // wetting front depth layer 1 and 2 in m
    double store1=(ThetaS1->Drc-ThetaI1->Drc);

    L1 = *L1p;
    L2 = *L2p;
    FFull = *FFullp;

    if (!SwitchTwoLayer)
    {
        if (L1 < SoilDepth1->Drc)
        {
            FFull = 0;
            // because drainage can reset moisture content
            dL1 = store1 > tiny ? fact/store1 : 0;
            // increase in depth (m) is actual infiltration/available porespace
            // do this always, correct if 1st layer is full
            L1 += dL1;
            // reaches bottom in this timestep, dL1 is remaning space (can be 0)
            // if impermeable fact = 0, if not impermeable fact is calculated with fixed L1
            if (L1 > SoilDepth1->Drc)
            {
                fact = (SoilDepth1->Drc-L1)*store1;
                L1 = SoilDepth1->Drc;
                if (SwitchImpermeable)
                    FFull = 1;
            }
        }
        else
        {
            L1 = SoilDepth1->Drc;
            if (SwitchImpermeable)
            {
                FFull = 1;
                fact = 0;
            }
        }
    }
    else  //twolayer
    {
        double store2=(ThetaS2->Drc-ThetaI2->Drc);
        //layer 1 available space

        // if not reached bottomof first layer
        if (L1 < SoilDepth1->Drc)
        {
            FFull = 0;
            dL1 = store1 > tiny ? fact/store1 : 0;
            L1 += dL1;
            if (L1 > SoilDepth1->Drc)
            {
                fact = (L1-SoilDepth1->Drc)*store1;
                L1 = SoilDepth1->Drc;
            }
        }
        else  // deeper than L1
        {
            if (L1+L2 < SoilDepth2->Drc)
            {
                FFull = 0;
                dL2 = store2 > tiny ? fact/store2 : 0;
                // increase in 2nd layer
                L2+=dL2;

                if ((L1+L2) > SoilDepth2->Drc)
                {
                    fact = ((L2+L1)-SoilDepth2->Drc)*store2;
                    L2 = SoilDepth2->Drc-L1;
                    if (SwitchImpermeable)
                        FFull = 1;
                }
            }
            else
            {
                L2 = SoilDepth2->Drc-L1;
                if (SwitchImpermeable)
                {
                    FFull = 1;
                    fact= 0;
                }
            }
        }
    }

    *L1p = (REAL8)L1;
    *L2p = (REAL8)L2;
    *FFullp = (REAL8)FFull;

    return fact;
    // return new fact

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
void TWorld::Infiltration(void)
{
    //Ksateff calculated before loop
    //NOTE: if crusting is calculated during event then move to the loop!

    switch (InfilMethod)
    {
    case INFIL_NONE :
        fill(*fact, 0.0);
        fill(*fpot, 0.0);
        return;
    case INFIL_SWATRE :
        fill(*tm, 0);
        FOR_ROW_COL_MV
        {
            if (FloodDomain->Drc == 0)
                tm->Drc = WH->Drc;
            else
                tm->Drc = hmx->Drc;
        }
        InfilSwatre(tm); // includes grasstrips, compaction etc, results in fact and fpot
        break;
    default:
        InfilMethodsNew();
    }

    // calc infiltrated volume for mass balance
    FOR_ROW_COL_MV
    {
        double cy = SoilWidthDX->Drc;
        if(FloodDomain->Drc > 0)
            cy = ChannelAdj->Drc;

        InfilVol->Drc = fact->Drc*cy*DX->Drc;
        // infil volume is WH before - water after
        // infil volume in flow equations is a separate variable
    }

}


//---------------------------------------------------------------------------
/*!
\brief function to calculate potential and actual infiltration rate according to
Green and Ampt, Smith and Parlange or Ksat subtraction.

This function calculates the actual potential infiltration according to Ksat, G&A or S&P \n
1 or 2 layers, then calls IncreaseInfiltrationDepth to increase the wetting front. \n
Then the appropriate domain water heights are adjusted and FSurplus is calculated.
*/
void TWorld::InfilMethodsNew()
{

    FOR_ROW_COL_MV
    {
        double fact1 = 0;
        double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
        double fwh = 0;
        double Psi = Psi1->Drc;
        double space = std::max(ThetaS1->Drc-ThetaI1->Drc, 0.0);

        if (space < tiny || Ks == 0 || SoilDepth1->Drc < tiny)
        {
            fpot->Drc = 0;
            fact->Drc = 0;
            FSurplus->Drc = 0;
            FFull->Drc = 1;
        }
        else
        {
            if (FloodDomain->Drc == 0)
                fwh = WH->Drc;
            else
                fwh = hmx->Drc;
            // select the appropriate domain water height for overpressure

            if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc)
            {
                Ks = std::min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
                // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
                Psi = Psi2->Drc;
                space = std::max(ThetaS2->Drc-ThetaI2->Drc, 0.0);
            }
            // two layers

            // calculate potential infiltration fpot in m, give a very large fpot in the beginning to start of the infil
            // process, the actual fact is then chosen anyway.
            switch (InfilMethod)
            {
            case INFIL_KSAT : fpot->Drc = Ks; break;
            case INFIL_GREENAMPT :
            case INFIL_GREENAMPT2 :
                fpot->Drc = L1->Drc+L2->Drc > tiny ? Ks*(1.0+(Psi+fwh)/(L1->Drc+L2->Drc)) : 1e10;
                break;
            case INFIL_SMITH :
            case INFIL_SMITH2 :
                double B = (fwh + Psi)*space;
                double Cdexp = B > 0.1 ? exp(Fcum->Drc/B) : 1.0;
                fpot->Drc = Cdexp > 1 ? Ks*Cdexp/(Cdexp-1): 1e10;
                break;
            }

            fact1 = std::min(fpot->Drc, fwh);
            // actual infil in m, cannot have more infil than water on the surface

            fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc, &FFull->Drc);
            // adjust fact and increase L1 and L2, for twolayer, impermeable etc

            // adjust the WH in the domains with new fact
            if(FloodDomain->Drc == 0)
            {
                if (WH->Drc < fact->Drc) // in case of rounding of errors, fact is equal to WH
                {
                    fact->Drc = WH->Drc;
                    WH->Drc = 0;
                }
                else
                    WH->Drc -= fact->Drc;
            }
            else
            {
                if (hmx->Drc < fact->Drc) // in case of rounding of errors, fact
                {
                    fact->Drc = hmx->Drc;
                    hmx->Drc = 0;
                }
                else
                    hmx->Drc -= fact->Drc;
            }

            Fcum->Drc += fact->Drc;
            // cumulative infil in m

            // calc surplus infiltration (negative in m) for kin wave
            if (FFull->Drc == 1)
                FSurplus->Drc = 0;
            else
            {
                FSurplus->Drc = std::min(0.0, fact->Drc - fpot->Drc);

                //limit surplus to available room in soil
                double room = 0;
                if (!SwitchTwoLayer || L1->Drc <  SoilDepth1->Drc)
                {
                    room = (SoilDepth1->Drc - L1->Drc)*(ThetaS1->Drc-ThetaI1->Drc);\
                    FSurplus->Drc = FSurplus->Drc < -room ? -room : FSurplus->Drc;
                }
                else
                    if (SwitchTwoLayer && L1->Drc == SoilDepth1->Drc)
                    {
                        room = (SoilDepth2->Drc - L2->Drc - L1->Drc)*(ThetaS1->Drc-ThetaI1->Drc);
                        FSurplus->Drc = FSurplus->Drc < -room ? -room : FSurplus->Drc;
                    }
            }
        }
    }
}
//---------------------------------------------------------------------------
/*!
 \brief Calculates changes in soilwater with percolation from the bottom of the profile.

  Calculates changes in soilwater with percolation from the bottom of the profile, \n
  resulting in the soil becoming dryer. Based on BrooksCorey type of percolation: \n
  percolation = ksat*(theta/pore)*bca, where bca = 5.55*qPow(Ksat2->Drc,-0.114); \n
  This is completely undocumented. Effects in ksateff do not influence percolation \n
  so ksat1 or ksat2 are used.
  */
void TWorld::SoilWater()
{
    if (!SwitchPercolation
            || InfilMethod == INFIL_SWATRE
            || SwitchImpermeable
            || InfilMethod == INFIL_NONE)
        return;

    FOR_ROW_COL_MV
    {
        double Percolation, bca;

        if (FFull->Drc == 1)
            continue;

        if (!SwitchTwoLayer)
        {
            bca = 5.55*qPow(Ksat1->Drc,-0.114);
            // Brooks corey value based on non lin regeression Saxton stuff
            double theta = 0.5*(ThetaSub->Drc + ThetaI1->Drc);
            Percolation = Ksat1->Drc * pow(theta/ThetaS1->Drc, bca)*_dt/3600000.0;
            // percolation = unsaturated K per timestep, Brooks Corey estimate

            if (L1->Drc > SoilDepth1->Drc-tiny)
            {
                L1->Drc = std::max(0.01, SoilDepth1->Drc-Percolation/(ThetaS1->Drc-ThetaI1->Drc+0.01));
                // cannot be less than 0.01= 1 cm to avoid misery
                // add 0.01 for safety to avoid division by zero
            }

            Soilwater->Drc = (SoilDepth1->Drc - L1->Drc)*ThetaI1->Drc;
            // max available water = unsat zone depth * thetai
            Percolation = std::min(Percolation, Soilwater->Drc);
            // cannot have more percolation than available water

            if (Soilwater->Drc-Percolation > 0.01*ThetaS1->Drc)
                Soilwater->Drc -= Percolation;
            // subtract percolation, cannot be less than residual theta is assumed 1 % of porosity

            ThetaI1->Drc = (SoilDepth1->Drc - L1->Drc > 0 ? Soilwater->Drc/(SoilDepth1->Drc - L1->Drc) : ThetaS1->Drc);
            ThetaI1->Drc = std::min(ThetaI1->Drc, ThetaS1->Drc);\
            // recalc thetai1 with new soilwater
        }
        else
        {
            bca = 5.55*qPow(Ksat2->Drc,-0.114);
            double theta = 0.5*(ThetaSub->Drc + ThetaI2->Drc);
            Percolation = Ksat2->Drc * pow(theta/ThetaS2->Drc, bca);

            if (L2->Drc > SoilDepth2->Drc-tiny)
            {
                L2->Drc = std::max(0.01, SoilDepth2->Drc-Percolation/(ThetaS2->Drc-ThetaI2->Drc+0.01));
                // cannot be less than 0.01= 1 cm to avoid misery
                // add 0.01 for safety to avoid division by zero
            }
            Soilwater->Drc = (SoilDepth2->Drc - L2->Drc)*ThetaI2->Drc;
            // max available water = unsat zone depth * thetai
            Percolation = std::min(Percolation, Soilwater->Drc);
            // cannot have more percolation than available water

            if (Soilwater->Drc-Percolation > 0.01*ThetaS2->Drc)
                Soilwater->Drc -= Percolation;
            // subtract percolation, cannot be less than residual theta is assumed 1 % of porosity

            ThetaI2->Drc = (SoilDepth2->Drc - L2->Drc > 0 ? Soilwater->Drc/(SoilDepth2->Drc - L2->Drc) : ThetaS2->Drc);
            ThetaI2->Drc = std::min(ThetaI2->Drc, ThetaS2->Drc);\
            // recalc thetai2 with new soilwater

        }

    }
}
