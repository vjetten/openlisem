
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
- double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, double *L1p, double *L2p, double *FFullp)
- void TWorld::Infiltration(void)
- void TWorld::InfiltrationFloodNew(void)
- void TWorld::SoilWater()
 */

#include <algorithm>
#include "lisemqt.h"
#include "global.h"
//#include "model.h"
#include "operation.h"

//NOTE fact and fpot have a unit of m (not m/s)

#define TINY 0.0001

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

            Poreeff->Drc = ThetaS1->Drc;
            Thetaeff->Drc = ThetaI1->Drc;

            if (SwitchInfilCompact) {
                Poreeff->Drc = ThetaS1->Drc*(1-CompactFraction->Drc) + PoreCompact->Drc*CompactFraction->Drc;
                Thetaeff->Drc = ThetaI1->Drc*(1-CompactFraction->Drc) +
                        ThetaI1->Drc*PoreCompact->Drc/ThetaS1->Drc *CompactFraction->Drc;
            }
            if (SwitchInfilCrust) {
                Poreeff->Drc = ThetaS1->Drc*(1-CrustFraction->Drc) + PoreCrust->Drc*CrustFraction->Drc;
                Thetaeff->Drc = ThetaI1->Drc*(1-CrustFraction->Drc) +
                        ThetaI1->Drc*PoreCrust->Drc/ThetaS1->Drc *CrustFraction->Drc;
            }
            if (SwitchGrassStrip) {
                Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;
                Poreeff->Drc = ThetaS1->Drc*(1-GrassFraction->Drc) + PoreGrass->Drc*GrassFraction->Drc;
                Thetaeff->Drc = ThetaI1->Drc*(1-GrassFraction->Drc) +
                        ThetaI1->Drc*PoreGrass->Drc/ThetaS1->Drc *GrassFraction->Drc;
            }

            if (SwitchHardsurface)
                Ksateff->Drc *= (1-HardSurface->Drc);

            if (SwitchHouses) {
                Ksateff->Drc *= (1-HouseCover->Drc);
            }

            Ksateff->Drc = std::max(0.0, Ksateff->Drc);
            Poreeff->Drc = std::max(0.0, Poreeff->Drc);
            Thetaeff->Drc = std::min(Thetaeff->Drc, Poreeff->Drc);
            // prob unneccessary

            Ksateff->Drc *= ksatCalibration;
            // apply runfile/iface calibration factor

        }
    }
//    report(*ThetaI2,"ti2.map");
//     report(*Thetaeff,"te.map");
//    report(*Ksateff, "kseff.map");
//        report(*Poreeff,"porec.map");
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
    SwatreStep(op.runstep, SwatreSoilModel, _WH, fpot, TileDrainSoil, thetaTop, tma);
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

        SwatreStep(op.runstep, SwatreSoilModelCrust, tm, tma, tmb, thetaTop, CrustFraction);
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

        SwatreStep(op.runstep, SwatreSoilModelCompact, tm, tma, tmb, tmc, CompactFraction);

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

        SwatreStep(op.runstep, SwatreSoilModelGrass, WHGrass, fpotgr, tmb, tmc, GrassFraction);

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
 - returns actual infiltration in mm, NOT rate in mm/h
*/

double TWorld::IncreaseInfiltrationDepthNew(int r, int c) //, double fact, double *L1p, double *L2p, double *FFullp)
{
    double store1 = (Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double store2 = 0;
    double L = L1->Drc + L2->Drc;
    double dfact = 0;
    double fact_out = 0;
    bool passing = false;

     if (SwitchTwoLayer)
         store2 = (ThetaS2->Drc-ThetaI2->Drc);

    if (FFull->Drc == 1 && SwitchImpermeable)
        return 0;  // act infil is zero, soil is full

    //L = L + fact->Drc/std::max(store, TINY);
    // add infil to wetting front
    FFull->Drc = 0;
    if (SwitchTwoLayer) {

        if (L < SoilDepth1->Drc) { // still in first layer
            L = L + fact->Drc/store1;

            if (L > SoilDepth1->Drc) { // moving into second layer
                dfact = fact->Drc - (SoilDepth1->Drc - L1->Drc) * store1;
                L1->Drc = SoilDepth1->Drc;
                passing = true;
            } else {
                fact_out = fact->Drc;
                L1->Drc = L;
            }
        } else {  // already in 2nd layer
             L = L + fact->Drc/store2;

             if (L > SoilDepth2->Drc) {
                 fact_out = fact->Drc - (SoilDepth2->Drc - L1->Drc + L2->Drc) * store2;
                 L2->Drc = SoilDepth2->Drc-SoilDepth1->Drc;
                 if (SwitchImpermeable)
                     FFull->Drc = 1;
             } else
                 fact_out = fact->Drc;
        }
        if (passing) {  // moving from depth 1 to 2
            L = L1->Drc + dfact/store2; // increase L with remaining fact
            if (L > SoilDepth2->Drc) {
                fact_out = fact->Drc - (SoilDepth2->Drc - L1->Drc + L2->Drc) * store2;
                L2->Drc = SoilDepth2->Drc-SoilDepth1->Drc;
                if (SwitchImpermeable)
                    FFull->Drc = 1;
            }
            else
                fact_out = fact->Drc;
        }
    } else { //single layer
        L = L + fact->Drc/store1;
        if (L > SoilDepth1->Drc) {
            fact_out = (SoilDepth1->Drc - L1->Drc) * store1;
            L1->Drc = SoilDepth1->Drc;
            if (SwitchImpermeable)
                FFull->Drc = 1;
        } else {
            L1->Drc = L;
            fact_out = fact->Drc;
        }
    }

    return fact_out;

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
        InfilSwatre(tm);

        FOR_ROW_COL_MV
        {
            if (FloodDomain->Drc == 0)
                WH->Drc = tm->Drc;
            else
                hmx->Drc = tm->Drc;
        }

        break;
    default:
        InfilMethodsNew();
        // this function results in an actual infiltration "fact" (in m) and
        // potential infiltration "fpot" (in m), according to G&A, S&P etc.
        // It deals with 1 or 2 layers and increase of water depth
    }

    // calc infiltrated volume for mass balance
    #pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        double cy = SoilWidthDX->Drc;
        if(FloodDomain->Drc > 0)
            cy = ChannelAdj->Drc;

        InfilVol->Drc = fact->Drc*cy*DX->Drc;
        // infil volume is WH before - water after
    }


}
//---------------------------------------------------------------------------
// Infiltration by Green and Ampt,Smith and Parlange or Ksat.
// All the same except for calculation of the potential infiltration fpot
// 1 layer and 2 layers
/*!
\brief function to calculate potential and actula infiltration rate according to
Green and Mapt, Smith and Parlange or Ksat subtraction.

This function calculates the potential infiltration according to Ksat, G&A or S&P \n
then calls IncreaseInfiltrationDepth to increase the wetting front.
*/

void TWorld::InfilMethodsNew()
{
    #pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        // default vars are first layer vars
        double fact1 = 0;
        double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
        double fwh = 0;
        double Psi = Psi1->Drc/100; // in m
        double space = std::max(Poreeff->Drc-Thetaeff->Drc, 0.0);
        double L = L1->Drc + L2->Drc;

        // get the correct water layer
        if (FloodDomain->Drc == 0)
            fwh = WH->Drc; //runoff
        else
            fwh = hmx->Drc; // flood
        // select the appropriate domain water height for overpressure

        //calculate potential insiltration rate fpot
        if (SwitchTwoLayer ) {

            // if wetting front in second layer set those vars
            if (L1->Drc > SoilDepth1->Drc)
            {
                Ks = std::min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
                // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
                Psi = Psi2->Drc/100;
                space = std::max(ThetaS2->Drc-ThetaI2->Drc, 0.0);
            }
        }

        switch (InfilMethod)
        {
        case INFIL_KSAT : fpot->Drc = Ks; break;
        case INFIL_GREENAMPT :
        case INFIL_GREENAMPT2 : fpot->Drc = Ks*(1.0+(Psi+fwh)/std::max(TINY, L)); break;
        case INFIL_SMITH :
        case INFIL_SMITH2 :
            double B = (fwh + Psi)*space;
            if (B > 0.01) {
                double Cdexp =exp(Fcum->Drc/B);
                fpot->Drc = Ks*exp(Fcum->Drc/B)/(exp(Fcum->Drc/B)-1);
            } else
                fpot->Drc = Ks;
            break;
        }

        fact->Drc = std::min(fpot->Drc, fwh);
        // actual infil in m, cannot have more infil than water on the surface

        fact->Drc = IncreaseInfiltrationDepthNew(r, c);
        // adjust fact and increase L1 and L2, for twolayer, impermeable etc

        // adjust the WH in the correct domain with new fact
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
        // increase cumulative infil in m

        // calc surplus infiltration (negative in m) for kin wave
        if (FFull->Drc == 1)
            FSurplus->Drc = 0;
        else
        {
            space = (SoilDepth1->Drc - L1->Drc)*(Poreeff->Drc-Thetaeff->Drc);
            if (SwitchTwoLayer && L2->Drc > 0)
                space = /* space +*/ (SoilDepth2->Drc - L2->Drc)*(ThetaS2->Drc-ThetaI2->Drc);
            // total spce left in the soil in m

            FSurplus->Drc = -1.0*std::min(space, std::max(0.0, fpot->Drc-fact->Drc));
            // negative and smallest of space or fpot-fact
        }
    }
//report(*fpot,"fpot");
//report(*L1,"la");
//report(*L2,"lb");
//report(*fact,"factb");
}
//---------------------------------------------------------------------------

/*!
 \brief Calculates changes in soilwater with percolation from the bottom of the profile.

  Calculates changes in soilwater with percolation from the bottom of the profile, \n
  resulting in the soil becoming dryer. Based on BrooksCorey type of percolation: \n
  percolation = ksat*(theta/pore)*bca, where bca = 5.55*qPow(Ksat2->Drc,-0.114); \n
  This is completely undocumented. The soil is either impermeable or has percolation. \n
  */
//
void TWorld::SoilWater()
{
    if (InfilMethod == INFIL_SWATRE || InfilMethod == INFIL_NONE)
        return;
    if (SwitchImpermeable)
        return;
#pragma omp parallel for collapse(2) num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        double Percolation, Ks, bca, dL, thetar, theta_E;

        Percolation = 0;

        if(SwitchTwoLayer) {
            thetar = 0.025 * ThetaS2->Drc;
            if(ThetaI2->Drc > thetar) {
                Ks = Ksat2->Drc*_dt/3600000.0;
                bca = 5.55*qPow(Ksat2->Drc,-0.114);
                theta_E = (ThetaI2->Drc-thetar)/(ThetaS2->Drc-thetar);
                Percolation = Ks * pow(theta_E, bca);
                // percolation in m

                dL = std::min(SoilDepth2->Drc - SoilDepth1->Drc, SoilDepth2->Drc - L1->Drc - L2->Drc);
                // dL is the second layer when L2 is 0 or the remainder of the second layer when L2 > 0
                // so we only have to deal with ThetaI2
                double moisture = dL*(ThetaI2->Drc-thetar);

                if (moisture > Percolation) {
                    // decrease thetaeff because of percolation
                    moisture -= Percolation;
                    ThetaI2->Drc = moisture/dL+thetar;
                } else {
                    // decrease thetaeff because of percolation and decrease L1 with the remaining bit
                    double dP = Percolation - moisture;
                    ThetaI2->Drc = thetar;
                    L2->Drc -= std::max(0.0, dP/(ThetaS2->Drc - thetar));
                }

                if (L1->Drc+L2->Drc < SoilDepth2->Drc)
                    FFull->Drc = 0;
            }
        } else {
            thetar = 0.025 * Poreeff->Drc;
            if(Thetaeff->Drc > thetar) {
                Ks = Ksateff->Drc*_dt/3600000.0;
                bca = 5.55*qPow(Ksateff->Drc,-0.114);
                theta_E = (Thetaeff->Drc-thetar)/(Poreeff->Drc-thetar);
                Percolation = Ks * pow(theta_E, bca);
            }

            dL = std::max(0.0, SoilDepth1->Drc - L1->Drc);
            double moisture = dL*(Thetaeff->Drc-thetar);
            if (moisture > Percolation) {
                // decrease thetaeff because of percolation
                moisture -= Percolation;
                Thetaeff->Drc = moisture/dL+thetar;
            } else {
                // decrease thetaeff because of percolation and decrease L1 with the remaining bit
                double dP = Percolation - moisture;
                Thetaeff->Drc = thetar;
                L1->Drc -= std::max(0.0, dP/(ThetaS1->Drc - thetar));
            }


            if (L1->Drc < SoilDepth1->Drc)
                FFull->Drc = 0;

        }
        Perc->Drc = Percolation;
    }
}


//---------------------------------------------------------------------------
void TWorld::infilInWave(cTMap *_h, double dt1)
{
    if (InfilMethod == INFIL_SWATRE || InfilMethod == INFIL_NONE)
        return;
#pragma omp parallel for collapse(2) num_threads(userCores)
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
    }
}
//---------------------------------------------------------------------------
double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFullp)
{
    double dL1, dL2; // increase in wetting front layer 1 and 2 in m
    double L1, L2, FFull; // wetting front depth layer 1 and 2 in m
    double store1=(Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double dfact = 0;

    L1 = *L1p;
    L2 = *L2p;
    FFull = *FFullp; // flag profile full

    // infiltration in top layer
    if (L1 < SoilDepth1->Drc) { // still room in layer 1
        FFull = 0;
        // because drainage can reset moisture content
        dL1 = store1 > TINY ? fact/store1 : 0;
        // increase in depth (m) is actual infiltration/available porespace
        double L1d = L1; // store old
        L1 += dL1; // add increase

        if (L1 > SoilDepth1->Drc-TINY) // fills up in this timestep
        {
            dfact = fact;
            fact = (SoilDepth1->Drc-L1d)*store1; //remaining
            dfact -= fact; // left over for infiltration in layer 2 is this exists
            L1 = SoilDepth1->Drc;
            FFull = 1;
        }
    } else {
        if (!SwitchTwoLayer){
            L1 = SoilDepth1->Drc;
            FFull = 1;
            fact = 0;
            dfact = 0;
            if (store1 < TINY)
                fact = Ksateff->Drc;
        }
    }


    if (SwitchTwoLayer)
    {
        double store2=(ThetaS2->Drc-ThetaI2->Drc);
        //layer 2 available space

        if (dfact > 0 || L1 > SoilDepth1->Drc-TINY) // infil beyond first layer
        {
            if (L1+L2 < SoilDepth2->Drc) {  // if there is room in layer 2
                FFull = 0;
                if (dfact > 0)
                    dL2 = store2 > TINY ? dfact/store2 : 0; // finish overflowing tiumestep
                else
                    dL2 = store2 > TINY ? fact/store2 : 0;  // L1 is full continue with L2
                // increase in 2nd layer
                double L2d = L2;
                L2+=dL2;

                if ((L1+L2) > SoilDepth2->Drc) // if beyond soildepth 2 in this timestep
                {
                    fact = (SoilDepth2->Drc-L1-L2d)*store2;
                    L2 = SoilDepth2->Drc-L1;
                    if (SwitchImpermeable)
                        FFull = 1;
                }
            }
            else {  // if L2 also full
                L2 = SoilDepth2->Drc-L1;
                if (SwitchImpermeable)
                {
                    FFull = 1;
                    fact = 0;
                }
            }
        }
    }

    *L1p = (REAL8)L1;
    *L2p = (REAL8)L2;
    *FFullp = (REAL8)FFull;

    return fact+dfact;
    // return new fact
}
