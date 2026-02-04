/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/
/*!
  \file lisInfiltration.cpp
  \brief infiltraton processes: Green and Ampt, Smith and Parlanage,  1 and 2 layer. main SWATRE call

functions: \n
- void TWorld::InfilEffectiveKsat(void)
- void TWorld::InfilDynamicCrusting()
- void TWorld::cell_InfilMethods(int r, int c)
- double TWorld::IncreaseInfiltrationDepthNew1(double fact_in, int r, int c)
- double TWorld::IncreaseInfiltrationDepthNew2(double fact_in, int r, int c)
- void TWorld::InfilSwatre()
 */

#include <algorithm>
#include "lisemqt.h"
#include "global.h"
#include "model.h"
#include "operation.h"

//---------------------------------------------------------------------------
void TWorld::InfilEffectiveKsat()
{
    if (!SwitchInfiltration || InfilMethod == INFIL_SWATRE)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Ksat1->Drc *= _dt/3600000.0; // mm/h to m oper timestep
        if (SwitchTwoLayer)
            Ksat2->Drc *= _dt/3600000.0;
        if (SwitchThreeLayer)
            Ksat3->Drc *= _dt/3600000.0;
        if (SwitchInfilCrust)
            KsatCrust->Drc *= _dt/3600000.0;
        if (SwitchInfilCompact)
            KsatCompact->Drc *= _dt/3600000.0;
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Ksateff->Drc = Ksat1->Drc;
        Poreeff->Drc = ThetaS1->Drc;
        Thetaeff->Drc = qMax(ThetaR1->Drc,ThetaI1->Drc);  // this resets the thetaeff to thetai1 all the time which is false!
        // moved to datainit

        // static crusted surfaces
        if (SwitchInfilCrust) {
            Ksateff->Drc = Ksateff->Drc*(1-CrustFraction->Drc) + KsatCrust->Drc*CrustFraction->Drc;
            Poreeff->Drc = Poreeff->Drc*(1-CrustFraction->Drc) + PoreCrust->Drc*CrustFraction->Drc;
        }

        // compacted surfaces
        if (SwitchInfilCompact) {
            Ksateff->Drc = Ksateff->Drc*(1-CompactFraction->Drc) + KsatCompact->Drc*CompactFraction->Drc;
            Poreeff->Drc = Poreeff->Drc*(1-CompactFraction->Drc) + PoreCompact->Drc*CompactFraction->Drc;
        }

        // grass strips? old concept?
        if (SwitchGrassStrip) {
            Ksateff->Drc = Ksateff->Drc*(1-GrassFraction->Drc) + KsatGrass->Drc*GrassFraction->Drc;
            Poreeff->Drc = ThetaS1->Drc*(1-GrassFraction->Drc) + PoreGrass->Drc*GrassFraction->Drc;
        }

        Ksateff->Drc = qMax(0.0, Ksateff->Drc); // ???? waarom

        Ksateff->Drc *= 1.0-fractionImperm->Drc;
        //fractionImperm was made for SWATRE, total of houses, roads, hard surfaces

        // to avoid pore is less than thetaR else nan in redistribution
        if (Poreeff->Drc < ThetaR1->Drc)
            ThetaR1->Drc = 0.5*Poreeff->Drc;

        // may be a problem in for instance redistribution
        if (SwitchWaveUser) {
            // when incoming wave, no infil in that area
            if (WHboundarea->Drc > 0) {
                Ksateff->Drc = 0;
                Poreeff->Drc = 0;
                Ksat1->Drc = 0;
                Ksat2->Drc = 0;
                ThetaS1->Drc = 0;
                ThetaS2->Drc = 0;
                ThetaI1->Drc = 0;
                ThetaI2->Drc = 0;
            }
        }

    }}
    report(*Ksateff,"ksateff.map");
}
//---------------------------------------------------------------------------
// Calculate effective Ksat based on surface structure, impermeable etc.
void TWorld::InfilDynamicCrusting()
{
    if (!SwitchInfiltration || InfilMethod == INFIL_SWATRE)
        return;

    if (!SwitchInfilCrust || !SwitchDynamicCrusting)
        return;

    // recalc ksateff and poreeff
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        CrustFraction->Drc = qMin(1.0, CrustFraction0->Drc + (1.0-exp(-0.2*qMax(0.0, RainCumCrust->Drc*1000))));
        // cumulative rain larger than 5 mm/h
        // exponential crusting proces with cumulative rainfall
        // from no crusting to full crusting at ~ 30 mm,
        // old research Jean Boiffin, multiple rainfall events in a growing season, progressive crusting

        Ksateff->Drc = KsatCrust->Drc*CrustFraction->Drc + Ksat1->Drc*(1-CrustFraction->Drc);
        Poreeff->Drc = PoreCrust->Drc*CrustFraction->Drc + ThetaS1->Drc*(1-CrustFraction->Drc);
        // use crustfraction in line with SWATRE
    }}
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

    if (Ksateff->Drc == 0)
        return;

    if (FloodDomain->Drc == 0) {
        fwh = WH->Drc; //runoff in kinwave or dyn wave
    } else {
        fwh = hmx->Drc; // flood in kin wave
    }
    // select the appropriate domain water height for overpressure

    fwh += MBm->Drc; // mass balance correction
    fwh = qMax(0.0,fwh);

    // only do infiltration on permeable soils, is now incorporated in ksateff
    //if (SoilWidthDX->Drc > 0 && fwh > 0) {
    if (fwh > 0) {
        //calculate potential infiltration rate fpot
        if (SwitchTwoLayer || SwitchThreeLayer) {
             // if wetting front in second layer calculate average Ks
            if (Lw->Drc > SoilDep1 && Lw->Drc < SoilDepth2->Drc) {
                switch (KavgType) {
                    case 0: Ks = ARITHavg(Ksateff->Drc,Ksat2->Drc); break;
                    case 1: Ks = SQRTavg(Ksateff->Drc,Ksat2->Drc); break;
                    case 2: Ks = HARMavg(Ksateff->Drc,Ksat2->Drc,SoilDep1,Lw->Drc-SoilDep1); break;
                    case 3: Ks = MINavg(Ksateff->Drc,Ksat2->Drc); break;
                }
                Psi = Psi2->Drc;
            }
            // if wetting front in third layer calculate average Ks
            if (Lw->Drc > SoilDepth2->Drc && Lw->Drc < SoilDepth3->Drc) {
                switch (KavgType) {
                    case 0: Ks = ARITHavg(Ks,Ksat3->Drc); break;
                    case 1: Ks = SQRTavg(Ks,Ksat3->Drc); break;
                    case 2: Ks = HARMavg(Ks,Ksat3->Drc,SoilDep1+SoilDepth2->Drc,Lw->Drc-SoilDep1-SoilDepth2->Drc); break;
                    case 3: Ks = MINavg(Ks,Ksat3->Drc); break;
                }
                Psi = Psi3->Drc;
            }
        }

        if (InfilMethod == INFIL_GREENAMPT)
            fpot_ = Ks*(1.0+(Psi+fwh)/qMax(1e-3, Lw->Drc));
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

        fact_ = qMin(fpot_, fwh);
        if (fact_ < 1e-10)
            fact_ = 0;
        // actual infil in m, cannot have more infil than water on the surface, includes rainfall

        if (fact_ > 0) {
            if (SwitchThreeLayer)
                fact_ = IncreaseInfiltrationDepthNew3(fact_, r, c);
            else
                if (SwitchTwoLayer)
                    fact_ = IncreaseInfiltrationDepthNew2(fact_, r, c);
                else
                    fact_ = IncreaseInfiltrationDepthNew1(fact_, r, c);
        }
        // adjust fact and increase Lw, for twolayer, impermeable etc

        if (fwh < fact_) {
            fact_ = fwh;
            fwh = 0;
        }
        else
            fwh -= fact_;

        if(FloodDomain->Drc == 0)
            WH->Drc = fwh;
        else
            hmx->Drc = fwh;
        hmxWH->Drc = WH->Drc + hmx->Drc;
        // adjust the WH in the correct domain with new fact

        Fcum->Drc += fact_; // for Smith and Parlange
        // increase cumulative infil in m
       // fact->Drc = fact_;
        //InfilVol->Drc = fact_* SoilWidthDX->Drc * DX->Drc;
        InfilVol->Drc = fact_* FlowWidth->Drc * DX->Drc;
        // calc infiltrated volume for mass balance
        // use flowwidth because Ksateff included impermeable surfaces anyway
    } else {
       // fact->Drc = 0;
        InfilVol->Drc = 0;
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
double TWorld::IncreaseInfiltrationDepthNew1(double fact_in, int r, int c)
{
    double dtheta1 = qMax(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the top layer
    double L = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;

    // impermeable and L reached SD1, no more infil, should also catch dtheta1 = 0;
    if (SwitchImpermeable && L > SoilDep1 - 0.001) {
        // profle filled
        Lw->Drc = SoilDep1;
        return 0;
    }

    if (SwitchGWflow) {
        // profile filled with GW
        if (GWWH->Drc >= SoilDepth1init->Drc-HMIN) {
            Lw->Drc = 0; //?? check
            return 0;
        }
    }

    Lnew = L + fact_in/qMax(dtheta1,0.01);
    // increase wetting front
    space = (SoilDep1 - L)*dtheta1;
    if(Lnew > SoilDep1 || space < fact_in) {
        // if the new L fills up the profile
        if (SwitchImpermeable)
            // if impermeable remaining space is infiltration
            fact_out = space;
        else
            fact_out = space + Perc->Drc;
            // was only percolation but actual infiltration is filled up space plus percolation
        Lnew = SoilDep1;
    } else {
        fact_out = fact_in;
    }

    Lw->Drc = qBound(0.0, Lnew, SoilDep1); // should not be necessary!
    return qBound(0.0, fact_out, fact_in);
}
//---------------------------------------------------------------------------
// order of processes:
// check if profile is full
// check if wettingfront is in layer 1, calculate actual infiltration, flag if infiltration moves into SL2
// if wetting front is in SL2, calc infil from remaining space if impermeaable, or remaining space + percolation of prev timestep if not impermeable
// if wetting front passes from SL1 into SL2, calculate all and check if the profile fills up
// if the wetting front is SL2-0.001 (depth - 1 mm) then flag full
double TWorld::IncreaseInfiltrationDepthNew2(double fact_in, int r, int c)
{
    double dtheta1 = qMax(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the layers
    double dtheta2 = qMax(0.0,ThetaS2->Drc-ThetaI2->Drc);
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;
    double L = Lw->Drc;
    double dfact2 = 0;
    bool passing = false;
    double space2 = 0;
    double thmin = 0.001;

    // profile is full
    if (SwitchImpermeable && L > SoilDep2 - 0.001) {
        Lw->Drc = SoilDep2;
        return 0;
    }

    // GW has filled up profile
    if (SwitchGWflow) {
       if (/*L >= SoilDep1 && */GWWH->Drc >= SoilDepth2init->Drc-HMIN) {
           Lw->Drc = 0;//SoilDep1;
           return 0;
       }
       // when GWWH fills soildep2 then soildep2 is 0 anyway
    }

    // L is in layer 1
    if (L <= SoilDep1) {
        Lnew = L + fact_in/qMax(thmin,dtheta1);
        space = (SoilDep1-L)*dtheta1;

        if(fact_in > space || Lnew > SoilDep1) {
            passing = true;
            // water is moving into layer 2
            dfact2 = fact_in - space;
            // remaining water for layer 2
        } else {
            fact_out = fact_in;
            // all remains SD1
        }
    }

    // L is in layer 2
    if (L > SoilDep1) {
        Lnew = L + fact_in/qMax(thmin,dtheta2);
        space2 = (SoilDep2-L)*dtheta2;

        if (Lnew > SoilDep2 || fact_in > space2) {
            if (SwitchImpermeable)
                fact_out = space2;
            else
                fact_out = space2 + Perc->Drc;

            Lnew = SoilDep2;
            // L at bottom
        } else {
            fact_out = fact_in;
            // everything fitted
        }
    }

    // L is moving from layer 1 into 2 in this timestep, with dfact2
    if (passing) {
        // second layer still at initial
        space2 = (SoilDep2-SoilDep1)*dtheta2;
        Lnew = SoilDep1 + dfact2/qMax(thmin,dtheta2);

        // if it moves from SL1 all the way and fills SL2
        if (dtheta2 < thmin || Lnew > SoilDep2-0.001) {
            if (SwitchImpermeable)
                fact_out = space + space2;
            else
                fact_out = space + space2 + Perc->Drc;
            Lnew = SoilDep2;
        } else
            fact_out = fact_in; // everything fitted
    }

    Lw->Drc = qBound(0.0, Lnew, SoilDep2); // should not be necessary, may hide errors!
    return qBound(0.0, fact_out, fact_in);
}
//---------------------------------------------------------------------------
// 3 layer infiltration! not used yet
// check if the profile is full
// check if wetting front is in layer 1, or already in layer 2 or already in layer 3
// then check if the infiltration passes from layer 1 into layer 2, may result in water passing from 2 to 3
// then check if the infiltration passes from layer 2 into layer 3

//3 LAYER NEEDS TO BE CHECKED
double TWorld::IncreaseInfiltrationDepthNew3(double fact_in, int r, int c)
{
    double dtheta1 = qMax(0.0,Poreeff->Drc-Thetaeff->Drc); // space in the layers
    double dtheta2 = qMax(0.0,ThetaS2->Drc-ThetaI2->Drc);
    double dtheta3 = qMax(0.0,ThetaS3->Drc-ThetaI3->Drc);
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double SoilDep3 = SoilDepth3->Drc;
    double fact_out = 0;
    double Lnew = 0;
    double L = Lw->Drc;
    double dfact12 = 0;
    double dfact23 = 0;
    bool passing12 = false;
    bool passing23 = false;
    double space = 0;
    double space2 = 0;
    double space3 = 0;
    double thmin = 0.001;

    // profile is full
    if (SwitchImpermeable && L > SoilDep3 - 0.001) {
        Lw->Drc = SoilDep3;
        return 0;
    }

    // GW has filled up profile
    if (SwitchGWflow) {
       if (/*L >= SoilDep2 && */GWWH->Drc >= SoilDepth3init->Drc-0.001) {
           Lw->Drc = 0;//SoilDep2;
           return 0;
       }
    }

    // L is in layer 1
    if (L <= SoilDep1) {
        Lnew = L + fact_in/qMax(thmin,dtheta1);
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

    // if L is in layer 2
    if (L > SoilDep1 && L <= SoilDep2) {
        //L already in layer 2 but not in 3
        Lnew = L + fact_in/qMax(thmin,dtheta2);
        space2 = (SoilDep2-L)*dtheta2;

        // passing frm SL2 into SL3
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
        Lnew = L + fact_in/qMax(thmin,dtheta3);
        space3 = (SoilDep3-L)*dtheta3;

        if (fact_in > space3 || Lnew > SoilDep3) {
            if (SwitchImpermeable)
                fact_out = space3;
            else
                fact_out = space3 + Perc->Drc;

            Lnew = SoilDep3;
            // L at bottom
        } else {
            // all remains SD3
            fact_out = fact_in;
        }
    }

    // L is moving from layer 1 into 2 in this timestep
    if (passing12) {
        // second layer still at initial
        space2 = (SoilDep2-SoilDep1)*dtheta2;
        Lnew = SoilDep1 + dfact12/qMax(thmin,dtheta2);

        // moves all the way into soillayer 3
        // note that the use can make e.g. the second layer saturated!, hence check dtheta2 < 0.01
        if (dtheta2 < thmin || Lnew > SoilDep2) {
            passing23 = true;
            dfact23 = fact_in - space2 - space;
        } else {
            fact_out = fact_in;
            // everything fitted in SD2
        }
    }

    // L is moving from layer 2 into 3 in this timestep
    if (passing23) {
        // second layer still at initial
        space3 = (SoilDep3-SoilDep2)*dtheta3;
        Lnew = SoilDep2 + dfact23/qMax(thmin,dtheta3);

        if (dtheta3 < thmin || Lnew > SoilDep3-0.001) {
            if (SwitchImpermeable)
                fact_out = space + space2 + space3; // note that space can be 0
            else
                fact_out = space + space2 + space3 + Perc->Drc;
            Lnew = SoilDep3;
        } else
            fact_out = fact_in; // everything fitted
    }

    Lw->Drc = qBound(0.0, Lnew, SoilDep3); // bounding should not be necessary!
    return qBound(0.0, fact_out, fact_in);


}
