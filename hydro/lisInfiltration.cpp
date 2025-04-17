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
        Thetaeff->Drc = std::max(ThetaR1->Drc,ThetaI1->Drc);  // this resets the thetaeff to thetai1 all the time which is false!
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

        // density factor and OM corrections directly in LISEM (instead of dbase creator)
        // because SWATRE also needs this
        // these correction come from calculations based on Saxton and Rawls
        // note ksat is in m/timestep, affects B of the regression eq for Ks, 0.001/3600.0*_dt
        if (SwitchOMCorrection) {
            double OM2 = OMcorr->Drc*OMcorr->Drc;
            double corrKsOA = 0.0026*OM2 + 0.0359*OMcorr->Drc + 1;
            double corrKsOB = 0.001/3600*_dt*(0.253*OM2 + 2.9368*OMcorr->Drc + 0.0007);
            double corrPOA  = -0.001*OM2 + 0.1014*OMcorr->Drc + 1.0;
            double corrPOB  = 0.0006*OM2 - 0.0282*OMcorr->Drc;
            Ksateff->Drc = corrKsOA*Ksateff->Drc + corrKsOB;
            Poreeff->Drc = corrPOA*Poreeff->Drc + corrPOB;
        }
        if (SwitchDensCorrection) {
            double D2 = DensFact->Drc*DensFact->Drc;
            double corrKsDA = 3.1429*D2 - 9.5657*DensFact->Drc + 7.4229;
            double corrKsDB = 0.001/3600.0*_dt*(135.4*D2 - 311.07*DensFact->Drc + 175.67);
            double corrPDA  = DensFact->Drc;
            double corrPDB   = -1.0 * DensFact->Drc + 1.0;
            Ksateff->Drc = corrKsDA*Ksateff->Drc + corrKsDB;
            Poreeff->Drc = corrPDA*Poreeff->Drc + corrPDB;
        }
        Ksateff->Drc = std::max(0.0, Ksateff->Drc); // ???? waarom

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

    if (!SwitchInfilCrust && !SwitchDynamicCrusting)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        //reset Ksateff and Poreeff
        Ksateff->Drc = Ksat1->Drc;
        Poreeff->Drc = ThetaS1->Drc;

        CrustFraction->Drc = 1.0-exp(-0.2*std::max(0.0, RainCumCrust->Drc*1000-5.0));  //
        // exponential crusting proces with cumulative rainfall
        // from no crusting to full crusting at ~ 30 mm,
        // old research Jean Boiffin, multiple rainfall events in a growing season, progressive crusting

        // double ksatdiff = std::max(0.0,Ksat1->Drc - KsatCrust->Drc);
        // Ksateff->Drc = KsatCrust->Drc + ksatdiff * factor;

        // double porediff = std::max(0.0,ThetaS1->Drc - PoreCrust->Drc);
        // Poreeff->Drc = PoreCrust->Drc + porediff * factor;
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
    double SoilDep2 = 0;

    if (Ksateff->Drc == 0)
        return;

    if (FloodDomain->Drc == 0) {
        fwh = WH->Drc; //runoff in kinwave or dyn wave
    } else {
        fwh = hmx->Drc; // flood in kin wave
    }
    // select the appropriate domain water height for overpressure

    fwh += MBm->Drc; // mass balance correction
    fwh = std::max(0.0,fwh);

    // only do infiltration on permeable soils, is now incorporated in ksateff
    //if (SoilWidthDX->Drc > 0 && fwh > 0) {
    if (fwh > 0) {
        //calculate potential infiltration rate fpot
        if (SwitchTwoLayer || SwitchThreeLayer) {
            SoilDep2 = SoilDepth2->Drc;
            // if wetting front in second layer set those vars
            if (Lw->Drc > SoilDep1 && Lw->Drc < SoilDep2) {
                //weighed harmonic mean:
                //https://corporatefinanceinstitute.com/resources/data-science/harmonic-mean/
                // sum (weights) / sum (weight/variable)
 //               Ks = Havg(Ksateff->Drc,Ksat2->Drc,SoilDep1,Lw->Drc-SoilDep1);
                Ks = Lw->Drc/(SoilDep1/Ksateff->Drc+(Lw->Drc-SoilDep1)/Ksat2->Drc);
                // if wetting front > layer 1 than ksat is determined weighted average (harmonic mean)
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

    if (SwitchGWflow) {
        if (GWWH->Drc >= SoilDepth1init->Drc-HMIN) {
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
            fact_out = Perc->Drc;
        Lnew = SoilDep1;
    } else {
        fact_out = fact_in;
    }

    Lnew = std::min(SoilDep1,std::max(0.0, Lnew));

    if (SwitchIncludeTile) {
      if (Lnew > TileDepth->Drc) {
          double vol = Dx->Drc*Ksat1->Drc*TileDiameter->Drc;
          // fraction of volume in layer 1, assuming full saturation so pore is draining
          TileWaterVolSoil->Drc = vol;



        }
      }


    Lw->Drc = Lnew;
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
void TWorld::InfilSwatre()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        // profile 0 is for impermeable surfaces
        if (ProfileID->Drc <= 0 || fractionImperm->Drc > 0.999) {
            InfilVol->Drc = 0;
            continue;
        }

        double WHorig;
        if (FloodDomain->Drc == 0)
            WHorig = WH->Drc;
        else
            WHorig = hmx->Drc;

        double drainfraction = 0;
        if (SwitchIncludeTile)
            drainfraction = TileWidth->Drc/_dx;
        SwatreSoilModel->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
        SwatreSoilModel->pixel[i_].tiledrain = 0;

        ComputeForPixel(i_, SwatreSoilModel, drainfraction);

        double WHN = SwatreSoilModel->pixel[i_].wh*0.01;
       //qDebug() << i_ << WHN;
        thetaTop->Drc = SwatreSoilModel->pixel[i_].theta; // not used!
        Perc->Drc= SwatreSoilModel->pixel[i_].percolation*0.01;
        if (SwitchIncludeTile)
            TileDrainSoil->Drc = SwatreSoilModel->pixel[i_].tiledrain*0.01;  // in m

        //TODO test infil swatre for crusts and compaction
        if (SwitchInfilCrust) {
            if (SwitchDynamicCrusting && ProfileIDCrust->Drc > 0) {
                CrustFraction->Drc = std::min(1.0, CrustFraction0->Drc + (1.0-exp(-0.2*std::max(0.0, RainCumCrust->Drc*1000-5.0))));
            }

            if (ProfileIDCrust->Drc > 0 && CrustFraction->Drc > 0) {
                SwatreSoilModelCrust->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
                SwatreSoilModelCrust->pixel[i_].tiledrain = 0;

                ComputeForPixel(i_, SwatreSoilModelCrust, 0.0);

                double WHcrust = SwatreSoilModel->pixel[i_].wh*0.01;

                double thetacrust = SwatreSoilModel->pixel[i_].theta;

                // weighed average
                WHN = WHcrust*CrustFraction->Drc + WHN*(1-CrustFraction->Drc);
                thetaTop->Drc = thetacrust*CrustFraction->Drc + thetaTop->Drc*(1-CrustFraction->Drc);
            }
        }

        if (SwitchInfilCompact) {
            if (ProfileIDCompact->Drc > 0 &&  CompactFraction->Drc > 0) {

                SwatreSoilModelCompact->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
                SwatreSoilModelCompact->pixel[i_].tiledrain = 0;

                ComputeForPixel(i_, SwatreSoilModelCompact, 0.0);

                double WHcompact = SwatreSoilModelCompact->pixel[i_].wh*0.01;
                double thetacompact = SwatreSoilModelCompact->pixel[i_].theta; // for pesticides ?

                // weighted average
                WHN = WHcompact*CompactFraction->Drc + WHN*(1-CompactFraction->Drc);
                thetaTop->Drc = thetacompact*CompactFraction->Drc + thetaTop->Drc*(1-CompactFraction->Drc);
            }
        }

        if (SwitchGrassStrip) {
            if (ProfileIDGrass->Drc > 0 &&  GrassFraction->Drc > 0) {
                SwatreSoilModelGrass->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
                SwatreSoilModelGrass->pixel[i_].tiledrain = 0;

                ComputeForPixel(i_, SwatreSoilModelGrass, 0.0);

                double WHgrass = SwatreSoilModelCompact->pixel[i_].wh*0.01;
                double thetagrass = SwatreSoilModelCompact->pixel[i_].theta; // for pesticides ?

                // weighted average
                WHN = WHgrass*GrassFraction->Drc + WHN*(1-GrassFraction->Drc);
                thetaTop->Drc = thetagrass*GrassFraction->Drc + thetaTop->Drc*(1-GrassFraction->Drc);
            }
        }

        if (FloodDomain->Drc == 0)
            WH->Drc = WHN;
        else
            hmx->Drc = WHN;

        InfilVol->Drc = (WHorig - WHN) * FlowWidth->Drc * DX->Drc;
        // use flowwidth because impermeable is done separately

    }}

    Copy(*thetaTop,*ThetaI1a);
    //for display

    //find depth wetting front, estimated at depth where h is initial value, very crude
    Fill(*Lwmm,0);
    for (int j = 0; j < SwatreSoilModel->pixel[0].profile->zone->nrNodes; j++) {
        cTMap *map = inith->at(j);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (j > 0 && SwatreSoilModel->pixel[i_].h[j] > map->Drc+1.0) {
                double l = SwatreSoilModel->pixel[i_].profile->zone->endComp[j-1]*10; // in mm
                double l1 = SwatreSoilModel->pixel[i_].profile->zone->endComp[j]*10; // in mm
                Lwmm->Drc = 0.5*(l+l1);
            }
        }}
    }

    // dump a map with h at every node
    if(SwitchDumphead) {
        for (int i = 0; i < SwatreSoilModel->pixel[0].profile->zone->nrNodes; i++) {

            QString dig = QString("%1").arg(i+1, 3, 10, QLatin1Char('0'));
            QString hname = QString("head0000.") + dig;
            QString tname = QString("theta000.") + dig;

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                hSwatre->Drc = SwatreSoilModel->pixel[i_].h[i];
                thetaSwatre->Drc = FindValue(hSwatre->Drc, SwatreSoilModel->pixel[i_].profile->horizon[i], H_COL, THETA_COL);
            }}
            report(*hSwatre, hname);
            report(*thetaSwatre, tname);
        }
    }
}
