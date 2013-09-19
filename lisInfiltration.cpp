
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
- void TWorld::InfilSwatre(void) \n
- void TWorld::InfilMorelSeytoux1(void) Not working yet! \n
- void TWorld::InfilSmithParlange1(void) \n
- void TWorld::InfilGreenAmpt1(void) \n
- void TWorld::InfilKsat(void) \n
- double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p)
      Called from all infilttration function, increases wetting front and calculates actual infiltration rate. \n
- void TWorld::Infiltration(void) Main infiltration function calling all infiltration methods,
      add net precipitation to WH and calculating new WH. \n
- void TWorld::SoilWater() Not mplemented yet! \n
 */


#include "model.h"

//NOTE fact and fpot have a unit of m (not m/s)

#define tiny 1e-8

//---------------------------------------------------------------------------
/// SWATRE infiltration, takes WH and calculateds new WH and infiltration surplus for kin wave
void TWorld::InfilSwatre(void)
{
    WHbef->copy(WH); // copy water height before infil

    tma->fill(1); // flag to indicate where the swatrestep model should be run

    //calculate a new crustfraction for water repellency
    // formula = f = 1/(1+1.2^(theta-30)), theta in %

    // for normal surface swatre should be done in all cells
    SwatreStep(SwatreSoilModel, WH, fpot, TileDrainSoil, thetaTop, tma);
    // NOTE WH changes in SWATRE
    // tiledrainsoil is in m per timestep, if not switchtiles then contains 0
    // TileDrainSoil->report("drain");

    // WH and fpot done in swatrestep
    FOR_ROW_COL_MV
            fact->Drc = (WHbef->Drc - WH->Drc);
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
        thetaTop->report("thtop");
        RepellencyFraction->report("repelfr");

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
        double Psi = Psi1->Drc;

        if (SoilDepth1->Drc <= tiny)
            continue;
        // outcrops etc: no infil

        if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
        {
            Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
            // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
            Psi = Psi2->Drc;
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

This function deals with 1 or 2 layer infiltration, calls advance of wetting front
and is called by normal soil infil and grass strip infil.
*/
void TWorld::InfilMethods(TMMap * _Ksateff, TMMap *_WH, TMMap *_fpot, TMMap *_fact, TMMap *_L1, TMMap *_L2, TMMap *_FFull)
{
    FOR_ROW_COL_MV
    {
        double fact1 = 0;
        double Ks = _Ksateff->Drc*_dt/3600000.0;  //in m
        double fwh = _WH->Drc; // in m, WH is old WH + net rainfall
        double Psi = Psi1->Drc;
        double space = max(ThetaS1->Drc-ThetaI1->Drc, tiny);

        if (SoilDepth1->Drc <= tiny)
        {
            _fpot->Drc = 0;
            _fact->Drc = 0;
            continue;
        }
        // outcrops etc: no infil in this cell

        if (SwitchTwoLayer && _L1->Drc > SoilDepth1->Drc - tiny)
        {
            Ks = min(_Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
            // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
            Psi = Psi2->Drc;
            space = max(ThetaS2->Drc-ThetaI2->Drc, tiny);
        }
        // two layers

        // calculate potential infiltration fpot in m
        switch (InfilMethod)
        {
        case INFIL_KSAT : _fpot->Drc = Ks; break;
        case INFIL_GREENAMPT :
        case INFIL_GREENAMPT2 :
            _fpot->Drc = Ks*(1.0+(Psi+fwh)/(_L1->Drc+_L2->Drc)); break;
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

 this function is called form all infiltration functions except Swatre
possible situations:
 one layer and not impermeable: L1+dL1
 one layer and impermeable: L1+dL1 untill L = soildepth1
 2 layer and not impermeable: L1+dL1 if L <= soildepth1 and L2+dL2 if L > soildepth1
 2 layer and impermeable: L1+dL1 if L <= soildepth1 and L2+dL2 untill L = soildepth2
*/
double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFullp)
{
    double dL1, dL2; // increase in wetting front layer 1 and 2 in m
    double L1, L2, FFull; // wetting front depth layer 1 and 2 in m
    L1 = *L1p;
    L2 = *L2p;
    FFull = *FFullp;

    if (FFull == 1)
        return 0;

    dL1 = fact/max(tiny, ThetaS1->Drc-ThetaI1->Drc);
    // increase in depth (m) is actual infiltration/available porespace
    // do this always, correct if 1st layer is full

    // correct dL1 in case of 1 layer and water reached the bottom
    if (!SwitchTwoLayer)
    {
        // reaches bottom in this timestep, dL1 is remaing space (can be 0)
        // if impermeable fact = 0, if not impermeable fact is calculated with fixed L1
        if (L1+dL1 > SoilDepth1->Drc)
        {
            dL1 = max(0, SoilDepth1->Drc - L1);
            if (SwitchImpermeable)
            {
                FFull = 1;
                fact = dL1 * (ThetaS1->Drc-ThetaI1->Drc);
            }
        }
    }

    L1 += dL1;
    // increase infiltration depth L1  = fact/avail pore space

    //infiltration can go on in the second layer & first layer is full
    if (SwitchTwoLayer)
    {
        // water enters into second layer
        if (L1 > SoilDepth1->Drc - tiny)
        {
            dL2 = fact/max(tiny, (ThetaS2->Drc-ThetaI2->Drc));
            // increase in 2nd layer

            // reaches bottom in this timestep, dL2 is remaing space (can be 0)
            // if impermeable fact = 0, if not impermeable fact is calculated with fixed L1+L2
            if (L2+dL2 > SoilDepth2->Drc && L2 < SoilDepth2->Drc)
            {
                dL2 = max(0, SoilDepth2->Drc - L2);
                if (SwitchImpermeable)
                {
                    FFull = 1;
                    fact = dL2 * (ThetaS2->Drc-ThetaI2->Drc);
                }
            }
            L2 += dL2;
            // increase infiltration depth L2  = fact/avail pore space
        }
    }

    *L1p = (REAL8)L1;
    *L2p = (REAL8)L2;

    return fact;
}
//---------------------------------------------------------------------------
/*!
 \brief Main infiltration function, calls infiltration types (SWATRE, Green and Ampt,
  Smith and Parlange, Ksat subtraction. Calculates effective Ksat based on different
  surface types (crust, compaction).

  Main infiltration function\n
  - add net rainfall and snowmelt to surface water WH\n
  - keep track of different surface types: grass strips, compaction, crusting, roads, hard surface
  - do SWATRE or one of the others, SWATRE is a different set of functions\n
  - for the other functions (Green and Ampt, Smith and Parlange, Ksat)
  - calculate effective Ksat\n
  - call one of the infiltration functions for the actual infiltration rate
    and the infiltration surplus for the kinematic wave\n
  - increase of infiltration depth/wetting front, same function for each infiltration model: L1, L2, Fcum
  - decrease of surface water layer WH and calculate infiltration volume
  */
void TWorld::Infiltration(void)
{
    InfilVol->fill(0);

    FOR_ROW_COL_MV
      //      if(FloodDomain->Drc == 0)
    {
        InfilVol->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

        // potential water volume on surface before infil

        // calculate effective ksat for various situations
        if (InfilMethod != INFIL_SWATRE  && InfilMethod != INFIL_NONE)
        {
            Ksateff->Drc = Ksat1->Drc;
//            if (SwitchInfilCrust)
//                Ksateff->Drc = Ksat1->Drc*(1-CrustFraction->Drc) + KsatCrust->Drc*CrustFraction->Drc;
//            if (SwitchInfilCompact)
//                Ksateff->Drc = Ksat1->Drc*(1-CompactFraction->Drc) + KsatCompact->Drc*CompactFraction->Drc;
//            if (SwitchInfilCrust && SwitchInfilCompact)
//                Ksateff->Drc = Ksat1->Drc*(1-CompactFraction->Drc-CrustFraction->Drc) +
//                        KsatCrust->Drc*CrustFraction->Drc + KsatCompact->Drc*CompactFraction->Drc;
            if (SwitchInfilCrust || SwitchInfilCompact)
                Ksateff->Drc = Ksat1->Drc*(1-CompactFraction->Drc-CrustFraction->Drc) +
                        KsatCrust->Drc*CrustFraction->Drc + KsatCompact->Drc*CompactFraction->Drc;
            // when not switched on fractions and ksat are 0
            // avg ksat of "normal" surface with crusting and compaction fraction
            // adjust effective infil for crusting and compaction
            //VJ 110106 adapted this calculation

            if (SwitchHardsurface && HardSurface->Drc > 0)
                Ksateff->Drc = (1-HardSurface->Drc)*Ksateff->Drc;// =  0;
            //VJ 110111 no infiltration on hard surfaces

            //houses
            if (SwitchHouses)
                Ksateff->Drc = Ksateff->Drc * (1-HouseCover->Drc);
            //VJ decrease ksat for celss with houses

            if (GrassFraction->Drc > 0)
                Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;

            Ksateff->Drc *= ksatCalibration;
            // apply runfile/iface calibration factor

            if (SwitchBuffers && !SwitchSedtrap)
                if(BufferID->Drc > 0)
                    Ksateff->Drc = 0;
            //VJ 1000608 no infil in buffers, , but sedtrap can have infil
        }
    } //row col

    switch (InfilMethod)
    {
    case INFIL_NONE : fact->fill(0); fpot->fill(0); break;
    case INFIL_SWATRE : InfilSwatre(); break;   // includes grasstrips, compaction etc
    case INFIL_KSAT :
    case INFIL_GREENAMPT :
    case INFIL_GREENAMPT2 :
    case INFIL_SMITH :
    case INFIL_SMITH2 :
        InfilMethods(Ksateff, WH, fpot, fact, L1, L2, FFull);
//        if(SwitchGrassStrip)
//            InfilMethods(KsatGrass, WHGrass, fpotgr, factgr, L1gr, L2gr, FFull);
        break;
    case INFIL_MOREL :
    case INFIL_HOLTAN : break;
    }
    // this function results in an actual infiltration "fact" (in m) and
    // potential infiltration "fpot" (in m), according to G&A, S&P etc.
    // It deals with 1 or 2 layers and increase of water depth
    // grass strips as a separate infiltration process (SWATRE grass is included)

    FOR_ROW_COL_MV
         //   if(FloodDomain->Drc == 0)
    {
        if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
        {
            if (RoadWidthDX->Drc == _dx)
            {
                fact->Drc = 0;
                fpot->Drc = 0;
            }
            // to make sure WH is not decreasd when all is road

            WH->Drc -= fact->Drc;
            if (WH->Drc < 0) // in case of rounding of errors
            {
                fact->Drc += WH->Drc;
                WH->Drc = 0;
            }
            // subtract fact->Drc from WH, cannot be more than WH

            Fcum->Drc += fact->Drc;
            // cumulative infil in m

//            if (GrassFraction->Drc > 0)
//            {
//                WHGrass->Drc -= factgr->Drc;
//                if (WHGrass->Drc < 0) // in case of rounding of errors
//                {
//                    factgr->Drc += WHGrass->Drc;
//                    WHGrass->Drc = 0;
//                }
//                Fcumgr->Drc += factgr->Drc;
//            }
            // calculate and correct water height on grass strips
        }

//        if (GrassFraction->Drc > 0)
//            WH->Drc = (1-GrassFraction->Drc) * WH->Drc + GrassFraction->Drc * WHGrass->Drc;
//        // average water height if grasstrip present

        FSurplus->Drc = min(0, fact->Drc - fpot->Drc);
        // negative surplus of infiltration in m for kinematic wave in m

//        if (GrassFraction->Drc > 0)
//            FSurplus->Drc = min(0, factgr->Drc - fpotgr->Drc);
//        // if grasstrip present use grasstrip surplus as entire surplus

        if (FFull->Drc == 1)
            FSurplus->Drc = 0;
        //VJ 101216 if soil full and impermeable: no surplus and no extra infil in kin wave

        InfilVol->Drc -= DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
        // infil volume is WH before - water after
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
    }

    FOR_ROW_COL_MV
            if(FloodDomain->Drc > 0)// && ChannelDepth->Drc == 0)
    {
        if (InfilMethod != INFIL_SWATRE)
        {
            if (RoadWidthDX->Drc == _dx)
            {
                fact->Drc = 0;
                fpot->Drc = 0;
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
}
//---------------------------------------------------------------------------
// calculates decrease in soil moisture under wetting front
// based on percolation = ksat*(theta/pore)*bca
// bca is related to ksat when looking at brooks corey
// effects in ksateff do not influence percolation so ksat1 or ksat2 are used
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
/*
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
*/
