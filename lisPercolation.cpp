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

// redistribution of soilwater after infiltration
// out[put is new Lw and new Thetaeff and ThetaI2
void TWorld::cell_Redistribution1(int r, int c)
{
    double Percolation, pore, theta, thetar, theta_E;
    double SoilDep2, pore2, theta2, thetar2, FC, FC2;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    pore = Poreeff->Drc;
    thetar = ThetaR1->Drc;
    theta = Thetaeff->Drc;
    FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;

    int ww = 0; //flag situation

    Percolation = 0;

    double factor = 0.1;
    // max amount of redistribution allowed

    if (SwitchTwoLayer) {
        pore2 = ThetaS2->Drc;
        thetar2 = ThetaR2->Drc;
        theta2 = ThetaI2->Drc;
        SoilDep2 = SoilDepth2->Drc;
        FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore2;
    }

    if (Lw_ < 0.1)
       return;
    // nothing to redistribute

    if (fact->Drc > 0)
        return;
    // no redistribution while there is infiltration

    if (SwitchImpermeable) {
        if (SwitchTwoLayer && Lw_ > SoilDepth2->Drc-0.01)
            return;
        if (!SwitchTwoLayer && Lw_ > SoilDepth1->Drc-0.01)
            return;
    }

    if(SwitchTwoLayer) {
        if (Lw_ < SoilDep1) {

            double store = (SoilDep1 - Lw_) * (pore-theta); // space in SD1 under Lw_

            if (store == 0) {
                Lw_ = SoilDep1;
                ww = -1;
            } else {
                theta_E = (theta-thetar)/(pore-thetar);
                Percolation = Ksateff->Drc*_dt/3600000.0 * pow(theta_E, bca1->Drc); // m/timestep
                Percolation = sqrt(Percolation * (Ksateff->Drc*_dt/3600000.0));
                //flux across boundary is ks+ke/2, average after Swatre
                double moistw = Lw_ * (pore-thetar); //available sat moisture above Lw_
                double dm = (pore-FC)*Lw_ * factor;
                Percolation = std::min(dm, Percolation);

                moistw -= Percolation; // decrease moistw with percolation, can be 0

                Lw_ = moistw/(pore-thetar); // new Lw_

                if (Percolation <= store) {
                    // if percolation fits in store layer 1 under the Lw
                    double m1 = (theta-thetar)*(SoilDep1-Lw_) + Percolation;
                    theta = m1/(SoilDep1-Lw_) + thetar;

                    if (theta >= pore) {
                        theta = pore;
                        Lw_= SoilDep1;
                    }
                    ww = 1;
                } else {

                    //double Perc1 = store;
                    //double Perc2 = std::max(0.0,Percolation - Perc1);
                    // some spills over in layer 2
                    // Lw_ is in layer 1
                    double m1 = (theta-thetar)*(SoilDep1-Lw_);
                    double m2 = (theta2-thetar2)*(SoilDep2-SoilDep1);
                    double Perc1 = m1/(m1+m2)*Percolation;
                    double Perc2 = m2/(m1+m2)*Percolation;

                    m1 = (theta-thetar)*(SoilDep1-Lw_) + Perc1;
                    theta = m1/(SoilDep1-Lw_)+thetar;

                    if (theta >= pore) {
                        theta = pore;
                        Lw_ = SoilDep1;
                    }

                    m2 = (theta2-thetar2)*(SoilDep2-SoilDep1) + Perc2;
                    theta2 = m2/(SoilDep2-SoilDep1) + thetar2;

                    if (theta2 >= pore2) {
                        theta2 = pore2;
                        Lw_ = SoilDep2;
                    }
                    ww = 2;
                }
            }
        } else {
            //Lw_ > SoilDep1

            // get the potential flux in m
            theta_E = (theta2-thetar2)/(pore2-thetar2);
            Percolation = Ksat2->Drc*_dt/3600000.0 * pow(theta_E, bca2->Drc); // m/timestep
            Percolation = sqrt(Percolation * (Ksat2->Drc*_dt/3600000.0));
            //flux across boundary is (ks+ke)/2, average after Swatre

            double moist1 = SoilDep1 * (pore - thetar); // moisture above wettingfront in layer 1
            double moist2 = (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2
            //available moisture
            double dm1 = (pore - FC)*SoilDep1*factor;
            double dm2 = (pore2 - FC2)*(Lw_-SoilDep1)*factor;
            //available moisture to move

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

            ww = 3;

            if (perc1 > 0) {
                // opercolation moves into layer 1
                moist1 -= perc1;
                // sat layer 1 now becomes free until FC
                Lw_ = moist1/(pore-thetar);
                theta = FC;// (FC)/(SoilDep1-Lw_);
                ww = 4;
            }
        }
        Thetaeff->Drc = theta;
        ThetaI2->Drc = theta2;
        Lw->Drc = Lw_;
    } // SwitchTwoLayer
}
//---------------------------------------------------------------------------
// redistribution of water from the wetting front to the layer below. Lw decreases and thetai increases
// percolation is done after this
void TWorld::cell_Redistribution(int r, int c)
{
    double Percolation;
    double Lw_ = Lw->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double SoilDep2 = (SwitchTwoLayer ? SoilDepth2->Drc : 0);

    Percolation = 0;

    double factor = 0.2;

    if (Lw_ < factor*SoilDep1)
       return;
    // nothing to redistribute

//    if (fact->Drc > 0)
  //      return;
    // no redistribution while there is infiltration

    if (SwitchImpermeable) {
        if (SwitchTwoLayer && Lw_ > SoilDep2-0.01)
            return;
        if (!SwitchTwoLayer && Lw_ > SoilDep1-0.01)
            return;
    }

    if (Lw_ < SoilDep1 && Lw_ > 0.1) {
        double pore = Poreeff->Drc;
        double thetar = ThetaR1->Drc;
        double ks = Ksateff->Drc*_dt/3600000.0;

        Percolation = ks * pow((Thetaeff->Drc-thetar)/(pore-thetar), bca1->Drc); // m/timestep
        Percolation = sqrt(Percolation * ks);
        //flux across boundary is sqrt(ks*ke), geom. average after Swatre

        double m = Lw_ * (pore-thetar); //available sat moisture above Lw_
        Percolation = std::min(m*factor, Percolation);

        double tn = (m+Percolation)/(SoilDep1-Lw_);
        if (tn >= thetar) {
            Thetaeff->Drc = tn;
            Lw->Drc = (m-Percolation)/(pore-thetar);
        }
        // subtract from moistw and recalc Lw
    }

    if (SwitchTwoLayer) {
        if (Lw_ > SoilDep2 + 0.1) {
            double pore2 = ThetaS2->Drc;
            double thetar2 = ThetaR2->Drc;
            double ks = Ksat2->Drc*_dt/3600000.0;

            // get the potential flux in m
            Percolation = ks * pow((ThetaI2->Drc-thetar2)/(pore2-thetar2), bca2->Drc); // m/timestep
            Percolation = sqrt(Percolation * ks);
            //flux across boundary is sqrt(ks*ke), geom. average after Swatre

            //get the available moisture
            double m =  (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2
            Percolation = std::min(factor * m, Percolation);

            double tn = (m+Percolation)/(SoilDep2-Lw_);
            if (tn >= thetar2) {
                ThetaI2->Drc = tn;
                Lw->Drc = (m-Percolation)/(pore2-thetar2) + SoilDep1;
            }
        }
    }
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
        double ksat = factor*Ksat2->Drc*_dt/3600000;

        if(theta > thetar) {
            // percolation in m per timestep
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, bca2->Drc);

            if (Lw->Drc > SoilDepth1->Drc)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // if Wet Fr still in first layer percolation only make 2nd drier

            double moisture = std::max(0.0, dL*(theta-thetar));
            if (Lw_ < SoilDep2-0.001) {
                // decrease thetaeff because of percolation
                Percolation = std::min(Percolation, moisture*0.5);
                moisture -= Percolation;
                theta = moisture/dL + thetar;
                theta = std::max(thetar, theta);
            } else {
                // wetting front = soildepth2, dL = 0, moisture = 0
                // assume theta goes back to 0.7 pore and decrease the wetting fornt
                Percolation = ksat;
                double FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore;
                theta = FC2;
                double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - Percolation/(pore - theta));
                Percolation = (Lwo-Lw_)*(pore-FC2);
            }            
            ThetaI2->Drc = theta;
            Lw->Drc = Lw_;
            return(Percolation);
        }
    } else {
        // one layer
        double pore = Poreeff->Drc;
        thetar = ThetaR1->Drc;
        double theta = Thetaeff->Drc;
        double ksat = factor*Ksateff->Drc*_dt/3600000;

        if (theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, bca1->Drc);
            if (Lw_ < SoilDep1-0.001) {
                // wetting front has not reached bottom, make soil drier
                // decrease thetaeff because of percolation
                dL = SoilDep1 - Lw_;
                double moisture = dL*(theta-thetar);
                Percolation = std::min(Percolation, moisture);
                moisture -= Percolation;
                theta = moisture/dL + thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume theta goes back to 0.7 pore and decrease the wetting fornt
                double FC = 0.7867*exp(-0.012*Ksateff->Drc)*pore;
                theta = FC;
                double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - ksat/(pore - theta));
                Percolation = (Lwo-Lw_)*(pore-FC);
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
