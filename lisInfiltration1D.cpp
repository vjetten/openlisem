

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

#include <algorithm>
#include "lisemqt.h"
#include "global.h"
#include "model.h"

void TWorld::InfilEffectiveKsat1D()
{

    if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            vKsateff[i_] = vKsat1[i_];
            vPoreeff[i_] = vThetaS1[i_];

            // exponential crusting proces with cumulative rainfall
            // crust fraction has no role now!!!
            if (SwitchInfilCrust) {
                //double KSc = Ksat1[i_] * (0.3+0.7*exp(-0.05*RainCum[i_]*1000));
                double ksatdiff = std::max(0.0,vKsat1[i_] - vKsatCrust[i_]);
                double factor = RainCum->Drc > 0.01 ? exp(-0.05*(RainCum->Drc-0.01)*1000) : 1.0;
                vKsateff[i_] = vKsatCrust[i_] + ksatdiff * factor;
                // exponential decline until crust value, RainCum is in meters
                //Ksateff[i_] = (1-Cover[i_]) * KSc + Cover[i_] * Ksat1[i_];

                // only on bare fraction of soil, depends on crop. We need basal cover! ???
                double porediff = std::max(0.0,vThetaS1[i_] - vPoreCrust[i_]);
                vPoreeff[i_] = vPoreCrust[i_] + porediff * factor;
            }

            // todo: dynamic and static crusting?
            //        if (SwitchInfilCrust) {
            //            vKsateff[i_] = vKsateff[i_]*(1-vCrustFraction[i_]) + vKsatCrust[i_]*vCrustFraction[i_];
            //            vPoreeff[i_] = vPoreeff[i_]*(1-vCrustFraction[i_]) + vPoreCrust[i_]*vCrustFraction[i_];
            //        }


            // affected surfaces, assumption that a compacted surface is not crusted
            if (SwitchInfilCompact) {
                vKsateff[i_] = vKsat1[i_]*(1-vCompactFraction[i_]) + vKsatCompact[i_]*vCompactFraction[i_];
                vPoreeff[i_] = vThetaS1[i_]*(1-vCompactFraction[i_]) + vPoreCompact[i_]*vCompactFraction[i_];
            }

            // a grass trip is not compacted or crusted so thetas1 and ksat1
            if (SwitchGrassStrip) {
                vKsateff[i_] = vKsat1[i_]*(1-vGrassFraction[i_]) + vKsatGrass[i_]*vGrassFraction[i_];
                vPoreeff[i_] = vThetaS1[i_]*(1-vGrassFraction[i_]) + vPoreGrass[i_]*vGrassFraction[i_];
            }


            if (SwitchHouses) {
                vKsateff[i_] *= (1-HouseCover->Drc);
            }

            //these surfaces are excluded from infiltration so not necessary to adjust Ksat and Pore
            // because they are not part of soilwidthDX!
            //            // impermeable surfaces
            //            if (SwitchHardsurface) {
            //                Ksateff[i_] *= (1-HardSurface[i_]);
            //             //   Poreeff[i_] *= (1-HardSurface[i_]);
            //            }

            //            if (SwitchRoadsystem) {
            //                Ksateff[i_] *= (1-RoadWidthDX[i_]/_dx);
            //             //   Poreeff[i_] *= (1-RoadWidthDX[i_]/_dx);
            //            }

            vKsateff[i_] = std::max(0.0, vKsateff[i_]);
            vPoreeff[i_] = std::max(0.3, vPoreeff[i_]);

        }}
    }

}


void TWorld::cell_InfilMethods1D(long i_, int r, int c)
{
    // default vars are first layer vars
    double Ks = vKsateff[i_];  //in m
    double Psi = vPsi1[i_]; // in m
    double fwh = 0;
    double fpot_ = 0;
    double fact_ = 0;
    double SoilDep1 = vSoilDepth1[i_];
    double SoilDep2 = 0;


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
            SoilDep2 = vSoilDepth2[i_];
            // if wetting front in second layer set those vars
            if (vLw[i_] > SoilDep1 && vLw[i_] < SoilDep2) {
                //weighed harmonic mean:
                //https://corporatefinanceinstitute.com/resources/data-science/harmonic-mean/
                // sum (weights) / sum (weight/variable)
                //               Ks = Havg(Ksateff->Drc,Ksat2->Drc,SoilDep1,Lw->Drc-SoilDep1);
                Ks = vLw[i_]/(SoilDep1/Ksateff->Drc+(vLw[i_]-SoilDep1)/vKsat2[i_]);
                // if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2
                Psi = vPsi2[i_]; //in m
            }
        }

        if (InfilMethod == INFIL_GREENAMPT)
            fpot_ = Ks*(1.0+(Psi+fwh)/std::max(1e-3, vLw[i_]));
        else {
            // smith parlange, not really tested
            double space = vPoreeff[i_]-vThetaeff[i_];
            if (vLw[i_] > vSoilDepth1[i_])
                space = vThetaS2[i_]-vThetaI2[i_];
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
            if (SwitchTwoLayer)
                fact_ = IncreaseInfiltrationDepthNew2_1D(fact_, i_,r,c);
            else
                fact_ = IncreaseInfiltrationDepthNew1_1D(fact_, i_,r,c);
        }
        // adjust fact and increase Lw, for twolayer, impermeable etc

        if (fwh < fact_)
        {
            fact_ = fwh;
            fwh = 0;
        }
        else
            fwh -= fact_;

        // copy back to maps
        if(FloodDomain->Drc == 0)
            WH->Drc = fwh;
        else
            hmx->Drc = fwh;
        // adjust the WH in the correct domain with new fact

        Fcum->Drc += fact_; // for Smith and Parlange
        // increase cumulative infil in m
        InfilVol->Drc = fact_* SoilWidthDX->Drc * DX->Drc;
        // calc infiltrated volume for mass balance
    } else {
        InfilVol->Drc = 0;
    }

    // calc surplus infiltration (negative in m) for kin wave
    // no longer used
    //FSurplus->Drc = 0;
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

//---------------------------------------------------------------------------
double TWorld::IncreaseInfiltrationDepthNew2_1D(double fact_in, long i_, int r, int c)
{
    double dtheta1 = std::max(0.0,vPoreeff[i_]-vThetaeff[i_]); // space in the top layer
    double dtheta2 = std::max(0.0,vThetaS2[i_]-vThetaI2[i_]);
    double SoilDep1 = vSoilDepth1[i_];
    double SoilDep2 = vSoilDepth2[i_];
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;
    double L = vLw[i_];
    double dfact2 = 0;
    bool passing = false;
    double space2 = 0;

    // profile is full
    if (SwitchImpermeable && L > SoilDep2 - 0.001) {
        vLw[i_] = SoilDep2;
        return 0;
    }

    if (SwitchGWflow) {
        if (L >= SoilDep1 && GWWH->Drc >= vSoilDepth2init[i_]-HMIN) {
            vLw[i_] = SoilDep1;
            return 0;
        }
        // when GWWH fills soildep2 than soildep2 is 0 anyway
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
                fact_out = vPerc[i_];

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
                fact_out = vPerc[i_];
            Lnew = SoilDep2;
        } else
            fact_out = fact_in; // everything fitted
    }

    vLw[i_] = std::min(SoilDep2,std::max(0.0, Lnew));
    return std::max(0.0,fact_out);
}
//---------------------------------------------------------------------------
// 3 layer infiltration! not used yet
double TWorld::IncreaseInfiltrationDepthNew3_1D(double fact_in, long i_, int r , int c)
{
    double dtheta1 = std::max(0.0,vPoreeff[i_]-vThetaeff[i_]); // space in the top layer
    double dtheta2 = std::max(0.0,vThetaS2[i_]-vThetaI2[i_]);
    double dtheta3 = std::max(0.0,vThetaS3[i_]-vThetaI3[i_]);
    double SoilDep1 = vSoilDepth1[i_];
    double SoilDep2 = vSoilDepth2[i_];
    double SoilDep3 = vSoilDepth3[i_];
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;
    double L = vLw[i_];
    double dfact12 = 0;
    double dfact23 = 0;
    bool passing12 = false;
    bool passing23 = false;
    double space2 = 0;
    double space3 = 0;

    // profile is full
    if (SwitchImpermeable && L > SoilDep2 - 0.001) {
        vLw[i_] = SoilDep2;
        return 0;
    }

    if (SwitchGWflow) {
        if (L >= SoilDep2 && GWWH->Drc >= vSoilDepth3init[i_]-HMIN) {
            vLw[i_] = SoilDep2;
            return 0;
        }
        // when GWWH fills soildep3 than soildep3 is 0 anyway
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
                fact_out = vPerc[i_];

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
                fact_out = vPerc[i_];
            Lnew = SoilDep3;
        } else
            fact_out = fact_in; // everything fitted
    }

    vLw[i_] = std::min(SoilDep3,std::max(0.0, Lnew));
    return std::max(0.0,fact_out);

}
//---------------------------------------------------------------------------
double TWorld::IncreaseInfiltrationDepthNew1_1D(double fact_in, long i_, int r, int c)
{
    double dtheta1 = std::max(0.0,vPoreeff[i_]-vThetaeff[i_]); // space in the top layer
    double L = vLw[i_];
    double SoilDep1 = vSoilDepth1[i_];
    double fact_out = 0;
    double space = 0;
    double Lnew = 0;

    // impermeable and L reached SD1, no more infil
    if (SwitchImpermeable && L > SoilDep1 - 0.001) {
        vLw[i_] = SoilDep1;
        return 0;
    }

    if (SwitchGWflow) {
        if (GWWH->Drc >= vSoilDepth1init[i_]-HMIN) {
            return 0;
        }
    }

    Lnew = L + fact_in/std::max(dtheta1,0.01);
    // increase wetting front
    space = (SoilDep1 - L)*dtheta1;
    if(Lnew > SoilDep1 || space < fact_in) {
        if (SwitchImpermeable)
            // if impermeable remaining space is infiltration
            fact_out = space;
        else
            fact_out = vPerc[i_];
        Lnew = SoilDep1;
    } else {
        fact_out = fact_in;
    }

    vLw[i_] = std::min(SoilDep1,std::max(0.0, Lnew));
    return std::max(0.0, fact_out);
}
//---------------------------------------------------------------------------

