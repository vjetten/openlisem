
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
- void TWorld::InfilSwatre(TMMap *_WH)
- void TWorld::InfilMorelSeytoux1(TMMap *_WH) Not working yet!
- void TWorld::InfilMethods(TMMap * _Ksateff, TMMap *_WH, TMMap *_fpot, TMMap *_fact, TMMap *_L1, TMMap *_L2, TMMap *_FFull)
- double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFullp)
- void TWorld::Infiltration(void)
- void TWorld::InfiltrationFloodNew(void)
- void TWorld::SoilWater()
 */


#include "model.h"

//NOTE fact and fpot have a unit of m (not m/s)

#define tiny 1e-8

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

            Ksateff->Drc = max(0, Ksateff->Drc);
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
void TWorld::InfilSwatre(TMMap *_WH)
{
    WHbef->copy(_WH); // copy water height before infil

    tma->fill(1); // flag to indicate where the swatrestep model should be run

    //calculate a new crustfraction for water repellency
    // formula = f = 1/(1+1.2^(theta-30)), theta in %

    // for normal surface swatre should be done in all cells
    SwatreStep(SwatreSoilModel, _WH, fpot, TileDrainSoil, thetaTop, tma);
    // NOTE WH changes in SWATRE
    // tiledrainsoil is in m per timestep, if not switchtiles then contains 0
    // TileDrainSoil->report("drain");

    // WH and fpot done in swatrestep
    FOR_ROW_COL_MV
            fact->Drc = (WHbef->Drc - _WH->Drc);
    // actual; infil is dif between WH before and after

    if (SwitchInfilCrust)
    {
        tm->copy(WHbef);
        tma->fill(0);
        tmb->fill(0);
        tmc->fill(0);

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
        tm->copy(WHbef);
        tma->fill(0);
        tmb->fill(0);
        tmc->fill(0);

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
        tm->copy(WHGrass);
        tmb->fill(0);
        tmc->fill(0);

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
/// DOESN'T WORK YET
void TWorld::InfilMorelSeytoux1(TMMap *_WH)
{
    FOR_ROW_COL_MV
    {
        double fact1;
        double Ks = Ksateff->Drc/3600000.0;  //in m/s
        double fwh = _WH->Drc; // in m, in WH is old WH + net rainfall
        double rt = fwh/_dt; // m/s pseudo rainfall, all water/dt = rate
        double A = 0, B, tp;
        double tt=time-BeginTime;
        // double Psi = Psi1->Drc;

        if (SoilDepth1->Drc <= tiny)
            continue;
        // outcrops etc: no infil

        if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
        {
            Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
            // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
            // Psi = Psi2->Drc;
        }

        B = (fwh+Psi1->Drc)*(ThetaS1->Drc-ThetaI1->Drc); //m
        tp = max(0, Ks*B/(rt*rt - Ks*rt)); //sec

        if (Ks < rt )
        {
            if (tp < tt+_dt)
                A = pow(B+Fcum->Drc,2)/(2*Ks*B*pow(rt/Ks-1,2));
            fpot->Drc = 0.5*sqrt(2*Ks*pow(B+Fcum->Drc,2)/B)*(1/sqrt(tt-tp+A))+Ks;
        }
        else
            fpot->Drc = Ks;
        // potential infiltration in m
        // psi input map is in cm so multiply 0.01 for m

        fact1 = min(fpot->Drc, fwh);
        // actual infil in m, cannot have more infil than water on the surface

        fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc, &FFull->Drc);

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
void TWorld::InfilMethods(TMMap * _Ksateff, TMMap *_WH, TMMap *_fpot, TMMap *_fact, TMMap *_L1, TMMap *_L2, TMMap *_FFull)
{
    FOR_ROW_COL_MV
    {
        double fact1 = 0;
        double Ks = _Ksateff->Drc*_dt/3600000.0;  //in m
        double fwh = _WH->Drc; // in m, WH is old WH + net rainfall
        double Psi = Psi1->Drc;
        double space = max(ThetaS1->Drc-ThetaI1->Drc, 0);

        if (SoilDepth1->Drc <= tiny)
        {
            _fpot->Drc = 0;
            _fact->Drc = 0;
            continue;
        }
        // outcrops etc: no infil in this cell

        if (SwitchTwoLayer && _L1->Drc > SoilDepth1->Drc)
        {
            Ks = min(_Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
            // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
            Psi = Psi2->Drc;
            space = max(ThetaS2->Drc-ThetaI2->Drc, 0);
        }
        // two layers

        // calculate potential infiltration fpot in m
        switch (InfilMethod)
        {
        case INFIL_KSAT : _fpot->Drc = Ks; break;
        case INFIL_GREENAMPT :
        case INFIL_GREENAMPT2 :
            _fpot->Drc = Ks*(1.0+(Psi+fwh)/(_L1->Drc+_L2->Drc+tiny)); break;
        case INFIL_SMITH :
        case INFIL_SMITH2 :
            double B = (fwh + Psi)*space;
            double Cdexp = exp(Fcum->Drc/B);
            _fpot->Drc = Ks*Cdexp/(Cdexp-1);
            break;
        }

        fact1 = min(_fpot->Drc, fwh);
        // actual infil in m, cannot have more infil than water on the surface

        _fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &_L1->Drc, &_L2->Drc, &_FFull->Drc);
        // adjust fact and increase L1 and L2, for twolayer, impermeable etc
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
            dL1 = fact/max(tiny, store1);
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
            dL1 = fact/max(tiny, store1);
            // tiny avoids dvision by 0
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
                dL2 = fact/max(tiny, store2);
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

    if (InfilMethod == INFIL_NONE)
        return;

    //Ksateff calculated before loop
    //NOTE: if crusting is calculated during event then move to the loop!

    FOR_ROW_COL_MV
    {
        InfilVol->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
    }

    switch (InfilMethod)
    {
    case INFIL_NONE : fact->fill(0); fpot->fill(0); break;
    case INFIL_SWATRE : InfilSwatre(WH); break;   // includes grasstrips, compaction etc
    case INFIL_KSAT :
    case INFIL_GREENAMPT :
    case INFIL_GREENAMPT2 :
    case INFIL_SMITH :
    case INFIL_SMITH2 :
        InfilMethods(Ksateff, WH, fpot, fact, L1, L2, FFull);
        break;
    case INFIL_MOREL :
    case INFIL_HOLTAN : break;
    }
    // this function results in an actual infiltration "fact" (in m) and
    // potential infiltration "fpot" (in m), according to G&A, S&P etc.
    // It deals with 1 or 2 layers and increase of water depth

    FOR_ROW_COL_MV
    {
        if (InfilMethod != INFIL_SWATRE)
        {
            //WH->Drc = max(0, WH->Drc - fact->Drc);
            WH->Drc -= fact->Drc;
            if (WH->Drc < 0) // in case of rounding of errors
            {
                fact->Drc += WH->Drc;
                WH->Drc = 0;
            }
            // decrease WH with actual infil, note fact is already limited to WH in InfilMethods
            Fcum->Drc += fact->Drc;
            // cumulative infil in m

        }

        // calc surplus infiltration (negative in m) for kin wave
        if (FFull->Drc == 1 || FloodDomain->Drc > 0)
            FSurplus->Drc = 0;
        else
        {
            FSurplus->Drc = min(0, fact->Drc - fpot->Drc);

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
        // negative surplus of infiltration in m for kinematic wave in m
        //VJ 101216 if soil full and impermeable: no surplus and no extra infil in kin wave
        //VJ 131222 limit so that smaller than available room!

        InfilVol->Drc -= DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
        // infil volume is WH before - water after
    }

//    FFull->report("full.map");
//    L1->report("l1.map");
//    L2->report("l2.map");
//    //    SoilDepth1->report("sd1.map");
//    //     SoilDepth2->report("sd2.map");
//    Fcum->report("fc.map");
//   FSurplus->report("fs.map");
}

//---------------------------------------------------------------------------
// do infiltration for flooded area exclusing channel cells
void TWorld::InfiltrationFloodNew(void)
{
    if (!SwitchChannelFlood)
        return;

    if (InfilMethod == INFIL_NONE)
        return;


    FOR_ROW_COL_MV
    {
        if (FloodDomain->Drc > 0)
        {
            InfilVolFlood->Drc = hmx->Drc * _dx*DX->Drc;//ChannelAdj->Drc*DX->Drc;
            // potential water volume on surface before infil
            // is everywhere on cell except channel

            if (ChannelWidthUpDX->Drc > 0)
                Ksateff->Drc *= (1-ChannelWidthUpDX->Drc/_dx);
            // if width is whole cell than adjust ksateff for channel cells
        }
        else
        {
            InfilVolFlood->Drc = 0;
        }
    }

    switch (InfilMethod)
    {
    case INFIL_NONE : fact->fill(0); fpot->fill(0); break;
    case INFIL_SWATRE : InfilSwatre(hmx); break;   // includes grasstrips, compaction etc
    case INFIL_KSAT :
    case INFIL_GREENAMPT :
    case INFIL_GREENAMPT2 :
    case INFIL_SMITH :
    case INFIL_SMITH2 :
        InfilMethods(Ksateff, hmx, fpot, fact, L1, L2, FFull); break;
    }

    FOR_ROW_COL_MV
            if(FloodDomain->Drc > 0)
    {
        if (InfilMethod != INFIL_SWATRE)
        {
            hmx->Drc = max(0, hmx->Drc - fact->Drc);

            Fcum->Drc += fact->Drc;
            // cumulative infil in m used in G&A infil function
        }

        FSurplus->Drc = 0;
        // no surplus in floods, no kin wave!

        //        InfilVolFlood->Drc -= ChannelAdj->Drc*hmx->Drc*DX->Drc;
        InfilVolFlood->Drc -= hmx->Drc * _dx*DX->Drc;
        // infil volume is WH before - water after
        // used for water balance and screen display
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
            Percolation = Ksat1->Drc * pow(ThetaI1->Drc/ThetaS1->Drc, bca)*_dt/3600000.0;
            // percolation = unsaturated K per timestep, Brooks Corey estimate

            if (L1->Drc > SoilDepth1->Drc-tiny)
            {
                L1->Drc = max(0.01, SoilDepth1->Drc-Percolation/(ThetaS1->Drc-ThetaI1->Drc+0.01));
                // cannot be less than 0.01= 1 cm to avoid misery
                // add 0.01 for safety to avoid division by zero
            }

            Soilwater->Drc = (SoilDepth1->Drc - L1->Drc)*ThetaI1->Drc;
            // max available water = unsat zone depth * thetai
            Percolation = min(Percolation, Soilwater->Drc);
            // cannot have more percolation than available water

            if (Soilwater->Drc-Percolation > 0.01*ThetaS1->Drc)
                Soilwater->Drc -= Percolation;
            // subtract percolation, cannot be less than residual theta is assumed 1 % of porosity

            ThetaI1->Drc = (SoilDepth1->Drc - L1->Drc > 0 ? Soilwater->Drc/(SoilDepth1->Drc - L1->Drc) : ThetaS1->Drc);
            ThetaI1->Drc = min(ThetaI1->Drc, ThetaS1->Drc);\
            // recalc thetai1 with new soilwater
        }
        else
        {
            bca = 5.55*qPow(Ksat2->Drc,-0.114);
            Percolation = Ksat2->Drc * pow(ThetaI2->Drc/ThetaS2->Drc, bca);

            if (L2->Drc > SoilDepth2->Drc-tiny)
            {
                L2->Drc = max(0.01, SoilDepth2->Drc-Percolation/(ThetaS2->Drc-ThetaI2->Drc+0.01));
                // cannot be less than 0.01= 1 cm to avoid misery
                // add 0.01 for safety to avoid division by zero
            }
            Soilwater->Drc = (SoilDepth2->Drc - L2->Drc)*ThetaI2->Drc;
            // max available water = unsat zone depth * thetai
            Percolation = min(Percolation, Soilwater->Drc);
            // cannot have more percolation than available water

            if (Soilwater->Drc-Percolation > 0.01*ThetaS2->Drc)
                Soilwater->Drc -= Percolation;
            // subtract percolation, cannot be less than residual theta is assumed 1 % of porosity

            ThetaI2->Drc = (SoilDepth2->Drc - L2->Drc > 0 ? Soilwater->Drc/(SoilDepth2->Drc - L2->Drc) : ThetaS2->Drc);
            ThetaI2->Drc = min(ThetaI2->Drc, ThetaS2->Drc);\
            // recalc thetai2 with new soilwater

        }

    }
}
//---------------------------------------------------------------------------
/* OBSOLETE
/// Solution Eurosem v2 manual page 10, Morgan et al 1998
void TWorld::InfilSmithParlange1(TMMap *_WH)
{
    FOR_ROW_COL_MV
    {
        double fact1;
        double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
        double fwh = _WH->Drc; // in m, in WH is old WH + net rainfall
        double Psi = Psi1->Drc;
        //VJ 110118 added psi of layer 2, was a bug
        double Cdexp, B;

        B = (fwh + Psi)*max(ThetaS1->Drc-ThetaI1->Drc, tiny);
        // TODO what abbout 2 layers?

        if (SoilDepth1->Drc <= tiny)
            continue;
        // outcrops etc: no infil


        if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
        {
            Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
            // smallest of the two ksats in two laters blocks the flow
            Psi = Psi2->Drc;
            //VJ 110118 added psi of layer 2, was a bug
        }
        Cdexp = exp(Fcum->Drc/B);
        fpot->Drc = Ks*Cdexp/(Cdexp-1);
        // potential infiltration in m

        fact1 = min(fpot->Drc, fwh);
        // actual infil in m, cannot have more infil than water on the surface

        fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc, &FFull->Drc);
        // adjust fact->Drc for twolayer, impermeable etc

        if(SwitchGrassStrip)
        {
            fwh = WHGrass->Drc;
            Ks = KsatGrass->Drc*_dt/3600000.0;  //in m

            if (SwitchTwoLayer && L1gr->Drc > SoilDepth1->Drc - tiny)
                Ks = min(KsatGrass->Drc, Ksat2->Drc)*_dt/3600000.0;

            B = (fwh + Psi)*max(ThetaS1->Drc-ThetaI1->Drc, tiny);

            Cdexp = exp(Fcum->Drc/B);
            fpotgr->Drc = Ks*Cdexp/(Cdexp-1);
            // potential infiltration in m

            fact1 = min(fpotgr->Drc, fwh);

            factgr->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc, &FFull->Drc);
        }
    }
}
//---------------------------------------------------------------------------
/// Solution Kutilek and Nielsen 2004 pag 138
void TWorld::InfilGreenAmpt1(TMMap *_WH)
{
    FOR_ROW_COL_MV
    {
        double fact1;
        double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
        double fwh = _WH->Drc; // in m, WH is old WH + net rainfall
        double Psi = Psi1->Drc;

        if (SoilDepth1->Drc <= tiny)
        {
            fpot->Drc = 0;
            fact->Drc = 0;
            continue;
        }
        // outcrops etc: no infil

        if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
        {
            Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
            // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
            Psi = Psi2->Drc;
            //VJ 110118 added psi of layer 2, was a bug
        }

        fpot->Drc = Ks*(1.0+(Psi+fwh)/(L1->Drc+L2->Drc));
        // potential infiltration in m, Darcy : Q = K * (dh/dz + 1)
        // L1 initialized at 1e-10, psi in cm so in datainit multiplied by 0.01 for m

        fact1 = min(fpot->Drc, fwh);
        // actual infil in m, cannot have more infil than water on the surface

        fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc, &FFull->Drc);
        // adjust fact and increase L1 and L2, for twolayer, impermeable etc

        if(SwitchGrassStrip)
        {
            fwh = WHGrass->Drc; // in m, WH is old WH + net rainfall
            Ks = KsatGrass->Drc*_dt/3600000.0;  //in m

            if (SwitchTwoLayer && L1gr->Drc > SoilDepth1->Drc - tiny)
                Ks = min(KsatGrass->Drc, Ksat2->Drc)*_dt/3600000.0;

            fpotgr->Drc = Ks*(1.0+(Psi+fwh)/(L1gr->Drc+L2->Drc));

            fact1 = min(fpotgr->Drc, fwh);

            factgr->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc, &FFull->Drc);
        }
    }
}
//---------------------------------------------------------------------------
/// Direct subtraction of Ksat, added for testing purposes!
void TWorld::InfilKsat(TMMap *_WH)
{
    FOR_ROW_COL_MV
    {
        double fact1;
        double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
        double fwh = _WH->Drc; // in m

        if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
            Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;

        fpot->Drc = Ks;
        // potential infil equals Ksat during this timestep, unit is m (not a flux)

        fact1 = min(fpot->Drc, fwh);
        // actual infil in m, cannot have more infil than water on the surface
        fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc, &FFull->Drc);
        // adjust fact for twolayer, impermeable etc

        if(SwitchGrassStrip)
        {
            fwh = WHGrass->Drc; // in m
            Ks = KsatGrass->Drc*_dt/3600000.0;  //in m

            if (SwitchTwoLayer && L1gr->Drc > SoilDepth1->Drc - tiny)
                Ks = min(KsatGrass->Drc, Ksat2->Drc)*_dt/3600000.0;

            fpotgr->Drc = Ks;

            fact1 = min(fpotgr->Drc, fwh);

            factgr->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc, &FFull->Drc);
        }
    }
}
//---------------------------------------------------------------------------
// do infiltration for flooded area exclusing channel cells
void TWorld::InfiltrationFlood(void)
{
    if (!SwitchChannelFlood)
        return;

    if (InfilMethod == INFIL_NONE)
        return;

    InfilVolFlood->fill(0);

    FOR_ROW_COL_MV
            if (FloodDomain->Drc > 0)// && ChannelDepth->Drc == 0)
    {
        //  hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        // potentially wrong!! flood height can be larger than plant height so no interception.
        // separate house interception

        InfilVolFlood->Drc = hmx->Drc * ChannelAdj->Drc*DX->Drc;
        // potential water volume on surface before infil

        // calculate effective ksat for various situations
        if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
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

            Ksateff->Drc *= ksatCalibration;
            // apply runfile/iface calibration factor

            if (SwitchBuffers && !SwitchSedtrap && SwitchBuffersImpermeable)
                if(BufferID->Drc > 0)
                    Ksateff->Drc = 0;
        }
    } //row col

    switch (InfilMethod)
    {
    case INFIL_NONE : ffact->fill(0); ffpot->fill(0); break;
    case INFIL_SWATRE : InfilSwatre(); break;   // includes grasstrips, compaction etc
    case INFIL_KSAT :
    case INFIL_GREENAMPT :
    case INFIL_GREENAMPT2 :
    case INFIL_SMITH :
    case INFIL_SMITH2 :
        InfilMethods(Ksateff, hmx, ffpot, ffact, Lf1, Lf2, FfFull); break;

        // WHY SEPARATE VARS FOR FLOOD AREAS? DEPTH is DEPTH and FULL IS FULL
    }

    FOR_ROW_COL_MV
            if(FloodDomain->Drc > 0)// && ChannelDepth->Drc == 0)
    {
        if (InfilMethod != INFIL_SWATRE)
        {
            if (RoadWidthDX->Drc == _dx)
            {
                ffact->Drc = 0;
                ffpot->Drc = 0;
            }

            hmx->Drc -= ffact->Drc;
            // decrease wh with net infil
            if (hmx->Drc < 0)
            {
                ffact->Drc += hmx->Drc;
                hmx->Drc = 0;
            }

            Ffcum->Drc += ffact->Drc;
            // cumulative infil in m used in G&A infil function
        }

        FfSurplus->Drc = min(0, ffact->Drc - ffpot->Drc);
        // negative surplus of infiltration in m for kinematic wave in m
        if (FfFull->Drc == 1)
            FfSurplus->Drc = 0;

        InfilVolFlood->Drc -= ChannelAdj->Drc*hmx->Drc*DX->Drc;
        // infil volume is WH before - water after
        // used for water balance and screen display
    }
}*/
