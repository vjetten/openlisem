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

#include "lisemqt.h"
#include "model.h"

#define ss_space 0.001

void TWorld::MoistureContent()
{
    thetai1cur = 0;
    #pragma omp parallel for reduction(+:thetai2cur) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Lw_ = std::max(SoilDepth1->Drc - Lw->Drc,0.0);
        thetai1cur += Lw_*ThetaI1->Drc*CHAdjDX->Drc;
    }}

    if (SwitchTwoLayer) {
        thetai2cur = 0;
        #pragma omp parallel for reduction(+:thetai2cur) num_threads(userCores)
        FOR_ROW_COL_MV_L {
            double Lw_ = Lw->Drc < SoilDepth1->Drc ? SoilDepth1->Drc : Lw->Drc;
            thetai2cur += (SoilDepth2->Drc - Lw_)*ThetaI2->Drc*CHAdjDX->Drc;
        }}
    }
}

// redistribution of soilwater after infiltration
// output is new Lw and new Thetaeff and ThetaI2
void TWorld::cell_Redistribution(int r, int c)
{

    double Percolation, pore, theta, thetar, theta_E;
    double SoilDep2, pore2, theta2, thetar2, FC, FC2, DL2;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    pore = Poreeff->Drc;
    thetar = ThetaR1->Drc;
    theta = Thetaeff->Drc;
    FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;
    Percolation = 0;

    double factor = 1.0;//0.5;
    // max amount of redistribution allowed

    if (SwitchTwoLayer) {
        pore2 = ThetaS2->Drc;
        thetar2 = ThetaR2->Drc;
        theta2 = ThetaI2->Drc;
        SoilDep2 = SoilDepth2->Drc;
        FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore2;
        DL2 = SoilDep2-SoilDep1;
    }

    if (SwitchImpermeable) {
        if (SwitchTwoLayer && Lw_ > SoilDepth2->Drc-0.01)
            return;
        if (!SwitchTwoLayer && Lw_ > SoilDepth1->Drc-0.01)
            return;
    }

    if(SwitchTwoLayer) {
        // no redistrib of Lw when too little
        if (Lw_ > 0.05) {
            if (Lw_ < SoilDep1) {

                // percolation flux
                theta_E = (theta-thetar)/(pore-thetar);
                Percolation = Ksateff->Drc * pow(theta_E, bca1->Drc); // m/timestep
                //   Percolation = sqrt(Percolation * Ksateff->Drc);
                Percolation = 0.5*(Percolation + Ksateff->Drc);
                //flux across boundary is ks+ke/2, average after Swatre
                double moistw = Lw_ * (pore-thetar); //available sat moisture above Lw_
                double dm = (pore-FC)*Lw_ * factor;
                Percolation = std::min(dm, Percolation);

                moistw -= Percolation; // decrease moistw with percolation, can be 0

                Lw_ = moistw/(pore-thetar); // new Lw_

                double store = (SoilDep1 - Lw_) * (pore-theta); // space in SD1 under Lw_
                if (Percolation <= store) {
                    // if percolation fits in store layer 1 under the Lw
                    double m1 = (theta-thetar)*(SoilDep1-Lw_) + Percolation;
                    theta = m1/(SoilDep1-Lw_) + thetar;

                    if (theta >= pore) {
                        theta = pore;
                        Lw_= SoilDep1;
                    }

                } else {
                    // some spills over in layer 2, Lw_ is in layer 1
                    double m1 = (theta-thetar)*(SoilDep1-Lw_);
                    double m2 = (theta2-thetar2)*DL2;
                    double Perc1 = m1/(m1+m2)*Percolation;
                    double Perc2 = m2/(m1+m2)*Percolation;

                    m1 = (theta-thetar)*(SoilDep1-Lw_) + Perc1;
                    theta = m1/(SoilDep1-Lw_)+thetar;

                    if (theta >= pore) {
                        theta = pore;
                        Lw_ = SoilDep1;
                    }

                    m2 = (theta2-thetar2)*DL2 + Perc2;
                    theta2 = m2/DL2 + thetar2;

                    if (theta2 >= pore2) {
                        theta2 = pore2;
                        Lw_ = SoilDep2;
                    }
                }
            } else {
                //Lw_ > SoilDep1

                theta_E = (theta2-thetar2)/(pore2-thetar2);
                Percolation = Ksat2->Drc* pow(theta_E, bca2->Drc); // m/timestep
                //  Percolation = sqrt(Percolation * Ksat2->Drc);
                Percolation = 0.5*(Percolation + Ksat2->Drc);

                double moist1 = SoilDep1 * (pore - thetar); // moisture above wettingfront in layer 1
                double moist2 = (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2

                double dm1 = (pore - FC)*SoilDep1*factor; //available moisture to move
                double dm2 = (pore2 - FC2)*(Lw_-SoilDep1)*factor;

                Percolation = std::min(Percolation, dm1+dm2);
                //cannot have more flux than available water

                double perc2 = std::min(Percolation, dm2); // part taken from layer 2
                double perc1 = std::max(0.0, Percolation-perc2); // part taken from layer 1, can be 0

                moist2 -= perc2;
                double m2 = (SoilDep2-Lw_)*(theta2-thetar2) + perc2;
                Lw_ = moist2/(pore2-thetar2) + SoilDep1;
                theta2 = m2/(SoilDep2-Lw_) + thetar2;
                theta = pore;

                if (theta2 >= pore2) {
                    Lw_ = SoilDep2;
                    theta2 = pore2;
                }

                if (perc1 > 0) {
                    // percolation moves into layer 1
                    moist1 -= perc1;
                    // sat layer 1 now becomes free until FC
                    Lw_ = moist1/(pore-thetar);
                    theta = FC;
                }
            }
            Thetaeff->Drc = theta;
            ThetaI2->Drc = theta2;
            Lw->Drc = Lw_;
        } // Lw_ > 0.1

        // redistribute the unsat zone between SD1 and SD2
        if (Lw_ < SoilDep1) { // && theta > theta2) {
            theta_E = (theta-thetar)/(pore-thetar);
            double Perc1 = Ksateff->Drc * pow(theta_E, bca1->Drc); // m/timestep
            theta_E = (theta2-thetar2)/(pore2-thetar2);
            double Perc2 = Ksat2->Drc * pow(theta_E, bca2->Drc); // m/timestep
            Percolation = 0.5*(Perc1+Perc2);

            double m1 = (SoilDep1-Lw_)*(theta-thetar);  // max moist
            Percolation = std::min(Percolation, m1);
            double m2 = DL2*(pore2-theta2); // max fit
            Percolation = std::min(Percolation, m2);

            m1 -= Percolation;
            theta = thetar + m1/(SoilDep1-Lw_);

            m2 = DL2*(theta2-thetar2) + Percolation;
            theta2 = thetar2 + m2/DL2;

            Thetaeff->Drc = theta;
            ThetaI2->Drc = theta2;
        }

    } else {
        // not SwitchTwoLayer
        if (Lw_ > 0.05) {
            if (Lw_ < SoilDep1-0.001) {
                theta_E = (theta-thetar)/(pore-thetar);// MC - percolation should depend on the difference between the two zones??
                //theta_E = 1; // MC - this percolation is from saturated zone to unsaturated - so theta_E = 1 ???
                Percolation = Ksateff->Drc * pow(theta_E, bca1->Drc); // m/timestep
                //  Percolation = sqrt(Percolation * Ksateff->Drc);
                Percolation = 0.5*(Percolation + Ksateff->Drc);
                //flux across boundary is ks+ke/2, average after Swatre

                double moisture = Lw_ * (pore-thetar); //available sat moisture above Lw_
                double dm = (pore-FC)*Lw_ * factor;
                Percolation = std::min(dm, Percolation);

                moisture -= Percolation;
                Lw_ = moisture/(pore-thetar); // new Lw_

                double store = (SoilDep1 - Lw_) * (pore-theta); // space in SD1 under Lw_
                Percolation = std::min(Percolation, store);

                // if percolation fits in store layer 1 under the Lw
                double m1 = (theta-thetar)*(SoilDep1-Lw_) + Percolation;
                theta = m1/(SoilDep1-Lw_) + thetar;

                if (theta >= pore) {
                    theta = pore;
                    Lw_= SoilDep1;
                }
            } // Lw < SD1
        } //Lw_ > 0.1
        Thetaeff->Drc= theta;
        Lw->Drc= Lw_;
    }// 1 layer
}
//---------------------------------------------------------------------------
// percolation from the bottom of the soil profile
//factor is for use of GW recharge
double TWorld::cell_Percolation(int r, int c, double factor)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    if(SwitchTwoLayer) {

        pore = ThetaS2->Drc;
        thetar = ThetaR2->Drc;
        theta = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;
        double ksat = factor*Ksat2->Drc;
        double FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore;

        if(theta > thetar) {
            // percolation in m per timestep
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, bca2->Drc);

            if (Lw_ > SoilDep1)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // if Wet Fr still in first layer percolation only make 2nd drier

            if (Lw_ < SoilDep2-0.001) {
                // decrease thetaeff because of percolation
                double moisture = dL*(theta-thetar);//*0.5);
                Percolation = std::min(Percolation, moisture);
                moisture -= Percolation;
                theta = moisture/dL + thetar;
                //theta = std::max(thetar, theta);
            } else {
                // wetting front = soildepth2, dL = 0, moisture = 0
                // assume theta goes back to FC2 and decrease the wetting fornt
                Percolation = ksat;
                theta = FC2;
                double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - Percolation/(pore - theta));
                Percolation = (Lwo-Lw_)*(pore-theta);
            }
            ThetaI2->Drc = theta;
            Lw->Drc = Lw_;
            return(Percolation);
        }
    } else {
        // one layer
        pore = Poreeff->Drc;
        thetar = ThetaR1->Drc;
        theta = Thetaeff->Drc;
        double ksat = factor*Ksateff->Drc;

        if (theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, bca1->Drc);

            if (Lw_ < SoilDep1-0.001) {
                // wetting front has not reached bottom, make soil drier
                // decrease thetaeff because of percolation
                    double moisture = (SoilDep1 - Lw_)*(theta-thetar);
                    Percolation = std::min(Percolation, moisture);
                    moisture -= Percolation;
                    theta = moisture/(SoilDep1 - Lw_) + thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume theta goes back to FC and decrease the wetting fornt
                double FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;
                theta = FC;
                double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - ksat/(pore - theta));
                Percolation = (Lwo-Lw_)*(pore-theta);
             }
            Thetaeff->Drc = theta;
            Lw->Drc = Lw_;
            return(Percolation);
        }
    }
    return(0);
}

//---------------------------------------------------------------------------

/*!
 \brief Calculates changes in soilwater with percolation from the bottom of the profile.

  Calculates changes in soilwater with percolation from the bottom of the profile, \n
  resulting in the soil becoming dryer. Based on BrooksCorey type of percolation: \n
  percolation = ksat*(theta/pore)*bca, where bca = 5.55*qPow(Ksat2->Drc,-0.114); \n
  This is completely undocumented. The soil is either impermeable or has percolation. \n
*/
// this function is not used!
void TWorld::SoilWater()
{
    if (InfilMethod == INFIL_SWATRE || InfilMethod == INFIL_NONE)
        return;
    if (SwitchImpermeable)
        return;

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        cell_Redistribution(r, c);

        Perc->Drc = cell_Percolation(r, c, 1.0);
    }}
}
