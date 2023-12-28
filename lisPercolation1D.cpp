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

#include "lisemqt.h"
#include "model.h"

//---------------------------------------------------------------------------
// percolation from the bottom of the soil profile
//factor is for use of GW recharge
double TWorld::cell_Percolation1D(long i_, int r, int c, double factor)
{
    double Percolation, dL, pore, theta, thetar, theta_E;
    double Lw_ = vLw[i_];
    double SoilDep1 = vSoilDepth1[i_];

    if(SwitchTwoLayer) {

        if (SwitchGWflow) {
            if (GWWH->Drc > vSoilDepth2[i_]-HMIN)
                return 0;
        }
        // no percolation to second layer if it is full with GW

        pore = vThetaS2[i_];
        thetar = vThetaR2[i_];
        theta = vThetaI2[i_];
        double SoilDep2 = vSoilDepth2[i_];
        double ksat = factor*vKsat2[i_];

        if(theta > thetar) {
            // percolation in m per timestep
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, 3.0+2.0/vlambda2[i_]);

            if (Lw_ > SoilDep1)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // assumption: if Wet Fr still in first layer percolation only make 2nd drier

            if (Lw_ < SoilDep2-0.001) {
                // decrease thetaeff because of percolation
                double moisture = dL*(theta-thetar);
                Percolation = std::min(Percolation, moisture);
                moisture -= Percolation;
                theta = moisture/dL + thetar;
            } else {
                // wetting front = soildepth2, dL = 0, moisture = 0
                // assume theta goes back to FC2 and decrease the wetting fornt
                theta = vThetaFC2[i_];
                //double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - ksat/(pore - theta));
                Percolation = ksat; //(Lwo-Lw_)*(pore-theta);
            }
            vThetaI2[i_] = theta;
            vLw[i_] = Lw_;

            return(Percolation);
        }
    } else {
        // one layer
        pore = vPoreeff[i_];
        thetar = vThetaR1[i_];
        theta = vThetaeff[i_];
        double ksat = factor*vKsateff[i_];

        if (SwitchGWflow && GWWH->Drc > vSoilDepth1[i_]-HMIN)
            return 0;


        if (theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, 3.0+2.0/vlambda1[i_]);

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
                theta = vThetaFC1[i_];
                Lw_ = std::max(0.0, Lw_ - ksat/(pore - theta));
                Percolation = ksat;//(Lwo-Lw_)*(pore-theta);
            }
            vThetaeff[i_] = theta;
            //DO NOT RECALCULATE PSI
            // Psi1[i_] = 0.01 * 10.2 * Psia1[i_] * psiCalibration * std::max(1.0, pow((theta-thetar)/(pore-thetar), -1.0/lambda1[i_]));

            vLw[i_] = Lw_;
            return(Percolation);
        }
    }
    return(0);
}


//---------------------------------------------------------------------------

void TWorld::cell_Redistribution1_1D(long i_, int r, int c)
{
    double Lw_ = vLw[i_];
    double L_min = 0.05; // minimum L before percolation starts

    if (SwitchImpermeable) {
        if (Lw_ > vSoilDepth1[i_]-0.001)
            return;
    }
    if (Lw_ < L_min)
        return;

    double Percolation, theta_E;

    double pore = vThetaeff[i_];
    double thetar = vThetaR1[i_];
    double theta = vThetaeff[i_];
    double SoilDep1 = vSoilDepth1[i_];
    double FC = vThetaFC1[i_];

    if (Lw_ < SoilDep1-0.001) {
        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = vKsateff[i_] * pow(theta_E, 3.0+2.0/vlambda1[i_]); // m/timestep
        //  Percolation = sqrt(Percolation * Ksateff[i_]);
        Percolation = Aavg(Percolation, vKsateff[i_]);
        //flux across boundary is ks+ke/2, average after Swatre

        double moisture = Lw_ * (pore-thetar); //available sat moisture above Lw_
        double dm = (pore-FC)*Lw_;
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
    }

    vThetaeff[i_] = theta;
    vLw[i_] = Lw_;
}
//---------------------------------------------------------------------------

void TWorld::cell_Redistribution2_1D(long i_, int r, int c)
{
   double Lw_ = vLw[i_];

   if (SwitchImpermeable) {
        if (Lw_ > vSoilDepth2[i_]-0.001)
        return;
   }

   double Percolation, theta_E;

   double pore = vPoreeff[i_];
   double thetar = vThetaR1[i_];
   double theta = vThetaeff[i_];
   double SoilDep1 = vSoilDepth1[i_];
   double FC1 = vThetaFC1[i_];

   double pore2 = vThetaS2[i_];
   double thetar2 = vThetaR2[i_];
   double theta2 = vThetaI2[i_];
   double SoilDep2 = vSoilDepth2[i_];
   double FC2 = vThetaFC2[i_];
   double DL2 = SoilDep2-SoilDep1;

   // if Lw still in layer 1
   if (Lw_ < SoilDep1) {
        // unsaturated flow between layer 1 and 2
        // avg percolation flux between layers
        // if there is room in layer 2
        if (theta2 < pore2-0.01) {
            theta_E = (theta-thetar)/(pore-thetar);
            double Perc1 = vKsateff[i_] * pow(theta_E, 3.0+2.0/vlambda1[i_]); // m/timestep
            theta_E = (theta2-thetar2)/(pore2-thetar2);
            double Perc2 = vKsat2[i_] * pow(theta_E, 3.0+2.0/vlambda2[i_]); // m/timestep
            Percolation = Aavg(Perc1, Perc2);

            double m1 = (SoilDep1-Lw_)*(theta-thetar);  // max moist
            Percolation = std::min(Percolation, m1);
            double m2 = DL2*(pore2-theta2); // max fit
            Percolation = std::min(Percolation, m2);

            theta = theta - Percolation/(SoilDep1-Lw_);
            theta = std::max(theta,thetar);

            theta2 = theta2 + Percolation/DL2;
            theta2 = std::min(pore2,theta2);
        }

        // decrease L with flow into the unsat zone beneath
        // percolation flux, avg K from L into unsat SD1
        theta_E = (theta-thetar)/(pore-thetar);
        double Percolation = vKsateff[i_] * pow(theta_E, 3.0+2.0/vlambda1[i_]); // m/timestep
        Percolation = Aavg(Percolation, vKsateff[i_]);

        double moistw = Lw_ * (pore-thetar);
        //available sat moisture above Lw_
        double dm = (pore-FC1)*Lw_;
        // max that can move
        Percolation = std::min(dm, Percolation);
        // not more percolation than putting moisture content at field capacity

        Lw_ = std::max(0.0,moistw-Percolation)/(pore-thetar);
        // new Lw_

        double store = (SoilDep1 - Lw_) * (pore-theta); // space in SD1 under Lw_
        if (Percolation <= store) {
            // if percolation fits in store layer 1 under the Lw
            theta = theta + Percolation/(SoilDep1-Lw_);

            // cannot happen, you cannot have flow from L into the unsaturated zone that saturates the layer under Lw!
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

            theta = theta + Perc1/(SoilDep1-Lw_);
            theta2 = theta2 + Perc2/DL2;

            // cannot happen!
            if (theta >= pore) {
                theta = pore;
                Lw_ = SoilDep1;
            }
            if (theta2 >= pore2) {
                theta2 = pore2;
                Lw_ = SoilDep2;
            }
        }
   } else {
        //Lw_ > SoilDep1, water from wettng szone into unsat below wetting zone in layer 2

        theta_E = (theta2-thetar2)/(pore2-thetar2);
        Percolation = vKsat2[i_]* pow(theta_E, 3.0+2.0/vlambda2[i_]); // m/timestep
        Percolation = Aavg(Percolation, vKsat2[i_]);
        double moist1 = SoilDep1 * (pore - thetar); // moisture above wettingfront in layer 1
        double moist2 = (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2

        double dm1 = (pore - FC1)*SoilDep1;
        double dm2 = (pore2 - FC2)*(Lw_-SoilDep1);
        //max available moisture to move from wetting front in layer 1 and 2
        Percolation = std::min(Percolation, dm1+dm2);
        //cannot have more flux than available water

        double perc2 = std::min(Percolation, dm2); // part taken from layer 2
        double perc1 = std::max(0.0, Percolation-perc2); // part taken from layer 1, can be 0

        // if so much percolation that Lw goes back into layer 1
        if (perc1 > 0) {
            Lw_ = (moist1-perc1)/(pore-thetar);
            theta = FC1 + perc1/(SoilDep1-Lw_);
            theta2 = theta2 + perc2/DL2;
            // unsat zone below Lw becomes wetter
        } else {
            Lw_ = (moist2-perc2)/(pore2-thetar2) + SoilDep1;
            // new wetting front decreased with perc2
            theta2 = theta2 + perc2/(SoilDep2-Lw_);
            // unsat zone below Lw becomes wetter
        }
   }

   vThetaeff[i_] = theta;
   vThetaI2[i_] = theta2;

   vLw[i_] = Lw_;
}
//---------------------------------------------------------------------------
void TWorld::avgTheta1D()
{
#pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_L {
        double Lw_ = vLw[i_];
        double SoilDep1 = vSoilDepth1[i_];
        vThetaI1a[i_] = vThetaeff[i_];

        if (Lw_ > 0 && Lw_ < SoilDep1 - 1e-3) {
            double f = Lw_/SoilDep1;
            vThetaI1a[i_] = f * vPoreeff[i_] + (1-f) *vThetaeff[i_];
        }
        if (Lw_ > SoilDep1 - 1e-3)
            vThetaI1a[i_] = vPoreeff[i_];

        if (SwitchTwoLayer) {
            double SoilDep2 = vSoilDepth2[i_];
            vThetaI2a[i_] = vThetaI2[i_];
            if (Lw_ > SoilDep1 && Lw_ < SoilDep2 - 1e-3) {
                double f = (Lw_-SoilDep1)/(SoilDep2-SoilDep1);
                vThetaI2a[i_] = f * vThetaS2[i_] + (1-f) *vThetaI2[i_];
            }
            if (Lw_ > SoilDep2 - 1e-3)
                vThetaI2a[i_] = vThetaS2[i_];
        }
   }}
}
