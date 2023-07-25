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

#define Aavg(a,b)  (0.5*(a+b))
#define Havg(a,b)  (2.0/(1.0/a+1.0/b))

double TWorld::SoilWaterMass()
{
    double totsatm3 = 0;
    double totunsatm3 = 0;

    FOR_ROW_COL_MV_L {
        double totsat = 0;
        double totunsat = 0;
        if (SwitchTwoLayer) {
            if (Lw->Drc <= SoilDepth1->Drc) {
                totsat = totsat + Lw->Drc * Poreeff->Drc;
                totunsat = totunsat + (SoilDepth1->Drc - Lw->Drc) * Thetaeff->Drc;
                totunsat = totunsat + (SoilDepth2->Drc - SoilDepth1->Drc) * ThetaI2->Drc;
            } else {
                totsat = totsat + SoilDepth1->Drc * Poreeff->Drc;
                totsat = totsat + (Lw->Drc-SoilDepth1->Drc) * ThetaS2->Drc;
                totunsat = totunsat + (SoilDepth1->Drc - SoilDepth2->Drc) * ThetaI2->Drc;
            }
        } else {
            totsat = totsat + Lw->Drc * Poreeff->Drc;
            totunsat = totunsat + (SoilDepth1->Drc - Lw->Drc) * Thetaeff->Drc;
        }
        totsatm3 += totsat * CHAdjDX->Drc;
        totunsatm3 += totsat * CHAdjDX->Drc;
    }}

    return totsatm3+totunsatm3;
}

//---------------------------------------------------------------------------
void TWorld::cell_Channelinfow1(int r, int c)
{
    ChannelQSide->Drc = 0.0;

//    if (ChannelWH->Drc > ChannelDepth->Drc - 0.05)
//        return;

    bool doUnsat = false;

    if (/* !doUnsat && */ Lw->Drc < 0.01)
        return;

   // double massbal = 0;
  //  double massbal2 = 0;

    double Lw_ = Lw->Drc;

    double pore = Poreeff->Drc;
    double thetar = ThetaR1->Drc;
    double theta = Thetaeff->Drc;
   // double SoilDep1 = SoilDepth1->Drc;
    double CHin1 = 0;
    double CHin2 = 0;
    double ChannelDep = ChannelDepth->Drc - ChannelWH->Drc - 0.05; // effective channel depth
    double K1 = Ksateff->Drc * pow((theta-thetar)/(pore-thetar), 3.0+2.0/lambda1->Drc); // m/timestep
    double DX_= DX->Drc;
    double dL = 0.5*ChannelAdj->Drc;

    CHin1 = Ksateff->Drc*2.0;
    CHin2 = K1 * 2.0;

    double h = std::min(ChannelDep,Lw_);

    double moist = Lw_*(pore-thetar);
    double dh = CHin1 * h*DX_/CHAdjDX->Drc * h/dL; // ks*cross section /cellsurface * Darcy pressure
    dh = std::min(dh, moist);
    moist -= dh;
    Lw_ = moist/(pore-thetar); // new Lw

//    double h2 = std::max(0.0, ChannelDep-Lw_);
//    if (doUnsat && theta > 0.95*pore && h2 > 0.01) {
//        moist = h2*(theta-thetar);
//        dh = CHin2*h2*DX_/CHAdjDX->Drc * h2/dL;
//        dh = std::min(dh, moist);
//        theta = thetar + moist/h2;
//    } else
//        CHin2 = 0;

    ChannelQSide->Drc = DX_*(CHin1*h*h/dL);// + CHin2*h2*h2/dL); // m3
}
//---------------------------------------------------------------------------
// Side inflow into channel from saturated part of the soil (Lw_), causes decrease of Lw_
// the assumption is thart the Darcy flow pressure difference dH/dL is 1.0
// afactor 2.0 is applied to Ksat because the flow is from both sides
void TWorld::cell_Channelinfow2(int r, int c)
{
    ChannelQSide->Drc = 0.0;

    if (ChannelWH->Drc > ChannelDepth->Drc - 0.05)
        return;

    if (Lw->Drc < 0.01)
        return;

    double Lw_ = Lw->Drc;
    double pore = Poreeff->Drc;
    double thetar = ThetaR1->Drc;
    double theta = Thetaeff->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double pore2 = ThetaS2->Drc;
    double thetar2 = ThetaR2->Drc;
    double theta2 = ThetaI2->Drc;
    double CHin1 = 0;
    double CHin2 = 0;
    double ChannelDep = ChannelDepth->Drc - ChannelWH->Drc; // effective channel depth
    double DX_= DX->Drc;
    double dL = 0.5* ChannelAdj->Drc;

    if (ChannelDep <= SoilDep1) {
        CHin1 = Ksateff->Drc*2.0;
        // sat layer 1
        double h = std::min(ChannelDep,Lw_);
        if (Lw_ > 0.01) {
            double moist = h*(pore-thetar);
            double frac = (h*DX_)/CHAdjDX->Drc; //= Lw_/ChannelAdj->Drc;
            double pressgrad = (h/dL);
            double dh = CHin1 * frac * pressgrad;
            dh = std::min(dh, moist);
            moist -= dh;
            Lw_ = moist/(pore-thetar); // new Lw
            CHin1 = dh/frac/pressgrad;
        } else
            CHin1 = 0;

        ChannelQSide->Drc = DX_*(CHin1*h*h/dL);

    } else {
        // chan > soildep1
        if (Lw_ <= SoilDep1) {

            CHin1 = Ksateff->Drc*2.0;

            // sat layer 1
            double h = Lw_;
            if (Lw_ > 0.01) {
                double moist = h*(pore-thetar);
                double frac = (h*DX_)/CHAdjDX->Drc; // Lw_/ChannelAdj->Drc;
                double pressgrad = (h/dL);
                double dh = CHin1 * frac * pressgrad;
                dh = std::min(dh, moist);
                moist -= dh;
                Lw_ = moist/(pore-thetar); // new Lw
                CHin1 = dh/frac/pressgrad;
            } else
                CHin1 = 0;

            ChannelQSide->Drc = DX_*(CHin1*h*h/dL);

        } else {
            // both chandep and Lw > soildep1

            CHin1 = Ksateff->Drc*2.0;
            CHin2 = Ksat2->Drc*2.0;
            double L = 0;
            double L2 =0;

            // layer 1 saturated
            double moist1 = SoilDep1*(pore-thetar);
            double frac = (SoilDep1*DX_)/CHAdjDX->Drc; // Lw_/ChannelAdj->Drc;
            double pressgrad = (SoilDep1/dL);
            double dh = CHin1 * frac * pressgrad;
            dh = std::min(dh, moist1);
            moist1 -= dh;
            L = moist1/(pore-thetar);
            CHin1 = dh/frac/pressgrad;

            // layer 2 saturated part, but not deeper than chandep
            double h2 = std::max(0.0,Lw_-SoilDep1);
            h2 = std::min(h2,ChannelDep-SoilDep1);
            if (h2 > 0.001) {
                double moist2 = h2*(pore2-thetar2);
                frac = (h2*DX_)/CHAdjDX->Drc; // Lw_/ChannelAdj->Drc;
                pressgrad = (h2/dL);
                dh = CHin2 * frac * pressgrad;
                dh = std::min(dh, moist2);
                moist2 -= dh;
                L2 = moist2/(pore-thetar);
                CHin2 = dh/frac/pressgrad;
            } else {
                h2 = 0;
                CHin2 = 0;
            }

            Lw_ = L + L2;

            ChannelQSide->Drc = DX->Drc*(CHin1*SoilDep1*SoilDep1/dL + CHin2*h2*h2/dL);
        }

    }

    if (!std::isnan(Lw_)) {
        Lw->Drc = Lw_;
        Thetaeff->Drc = theta;
        ThetaI2->Drc = theta2;
    }

    if (std::isnan(ChannelQSide->Drc)) {
        ChannelQSide->Drc = 0.0;
        //qDebug() << r << c << "nan" << CHin1 << CHin2 << CHin3 << Lw_ << i;
    }
    // update channel side inflow, sometimes nan occurs  in lw

}
//---------------------------------------------------------------------------
void TWorld::cell_Redistribution2psi(int r, int c)
{
    double Lw_ = Lw->Drc;

    if (SwitchImpermeable) {
        if (Lw_ > SoilDepth2->Drc-0.001)
            return;
    }

    double Percolation, theta_E;

    double pore1 = Poreeff->Drc;
    double thetar1 = ThetaR1->Drc;
    double theta1 = Thetaeff->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double FC1 = ThetaFC1->Drc;

    double pore2 = ThetaS2->Drc;
    double thetar2 = ThetaR2->Drc;
    double theta2 = ThetaI2->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double FC2 = ThetaFC2->Drc;
    double DL2 = SoilDep2-SoilDep1;


    double h0, h1, h2, h3;
    double dh1, dh2, dh3;
    double dz1, dz2, dz3;
    double K1, K2, K3;
    double Q1, Q2, Q3;
    double m0, m1, m2, m3;

    if (Lw_ < SoilDep1) {
        theta_E = (theta1-thetar1)/(pore1-thetar1);
        //se = (h / h b ) (-lambda)
        h1 = psi1ae->Drc*pow(theta_E, -1.0/lambda1->Drc);
        h0 = h1;
        K1 = Ksateff->Drc * pow(theta_E, 3.0+2.0/lambda1->Drc);
        K0 = K1;
        m1 = (SoilDep1-Lw_)*(theta1-thetar1);
        m0 = m1;
        if (Lw_ > 0) {
            h0 = 0;
            K0 = Ksateff->>Drc;
            m0 = Lw_*(pore1-thetar1);
        }

        theta_E = (theta2-thetar2)/(pore2-thetar2);
        h2 = psi2ae->Drc*pow(theta_E, -1.0/lambda2->Drc);
        h3 = h3;
        K2 = Ksat2->Drc * pow(theta_E, 3.0+2.0/lambda2->Drc);
        K3 = K2;
        m2 = (DL2-GWWH->Drc)*(theta2-thetar2);
        m3 = m2;
        if (GWWH->Drc > 0) {
            h3 = 0;
            K3 = Ksat2->Drc;
            m3 = GWWH->Drc*(pore2-thetar2);
        }

        dz1 = 0.5*SoilDep1;
        dz2 = 0.5*(DL2-GWWH->Drc) + (LW_ + 0.5*(SoilDep1 -Lw_));
        dz3 = 0.5*DL2;//0.5*GWWH-Drc + 0.5*(DL2-GWWH->Drc);
        dh1 = h0-h1;
        dh2 = h1-h2;
        dh3 = h2-h3;
        Q1 = -Aavg(K0, K1) * (dh1/dz1+1);
        Q2 = -Aavg(K1, K2) * (dh2/dz2+1);
        Q3 = -Aavg(K2, K3) * (dh3/dz3+1);


    }



    // if Lw still in layer 1
    if (Lw_ < SoilDep1) {
        // decrease L with flow into the unsat zone beneath

        // percolation flux, avg K from L into unsat SD1
        theta_E = (theta-thetar)/(pore-thetar);
        //se = (h / h b ) (-lambda)

        double psi1L = psi1ae->Drc*pow(theta_E, -1.0/lambda1->Drc);
        double Percolation = Ksateff->Drc * pow(theta_E, 3.0+2.0/lambda1->Drc); // m/timestep
        double dz = 0.5*SoilDep1; //Lw_ + 0.5*(SoilDep1-Lw_) - 0.5*Lw_;
        Percolation = -Percolation * (0-psi1L)/dz;
            //Lw_ + 0.5*SoilDep1 - 0.5*Lw_ - 0.5*Lw_ = 0.5*SD1?


        double moistw = Lw_ * (pore-thetar);
        //available sat moisture above Lw_
        double dm = (pore-FC1)*Lw_;
        Percolation = std::min(dm, Percolation);
        // not more percolation than putting moisture content at field capacity

        moistw -= Percolation;
        // decrease moistw with percolation, can be 0
        Lw_ = moistw/(pore-thetar);
        // new Lw_                

        double store = (SoilDep1 - Lw_) * (pore-theta); // space in SD1 under Lw_
        if (Percolation <= store) {
            // if percolation fits in store layer 1 under the Lw
            double m1 = (theta-thetar)*(SoilDep1-Lw_) + Percolation;
            theta = m1/(SoilDep1-Lw_) + thetar;

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

            m1 = (theta-thetar)*(SoilDep1-Lw_) + Perc1;
            theta = m1/(SoilDep1-Lw_)+thetar;

            // cannot happen!
            if (theta >= pore) {
                theta = pore;
                Lw_ = SoilDep1;
            }

            m2 = (theta2-thetar2)*DL2 + Perc2;
            theta2 = m2/DL2 + thetar2;

            // cannot happen!
            if (theta2 >= pore2) {
                theta2 = pore2;
                Lw_ = SoilDep2;
            }
        }

        // unsaturated flow between layer 1 and 2
        // avg percolation flux between layers
        theta_E = (theta-thetar)/(pore-thetar);
        double Perc1 = Ksateff->Drc * pow(theta_E, 3.0+2.0/lambda1->Drc); // m/timestep
        theta_E = (theta2-thetar2)/(pore2-thetar2);
        double Perc2 = Ksat2->Drc * pow(theta_E, 3.0+2.0/lambda2->Drc); // m/timestep
        Percolation = Aavg(Perc1, Perc2);

        double m1 = (SoilDep1-Lw_)*(theta-thetar);  // max moist
        Percolation = std::min(Percolation, m1);
        double m2 = DL2*(pore2-theta2); // max fit
        Percolation = std::min(Percolation, m2);

        m1 -= Percolation;
        theta = thetar + m1/(SoilDep1-Lw_);

        m2 = DL2*(theta2-thetar2) + Percolation;
        theta2 = thetar2 + m2/DL2;
        theta2 = std::min(pore2,theta2);

    } else {
        //Lw_ > SoilDep1

        theta_E = (theta2-thetar2)/(pore2-thetar2);
        Percolation = Ksat2->Drc* pow(theta_E, 3.0+2.0/lambda2->Drc); // m/timestep
        Percolation = Aavg(Percolation, Ksat2->Drc);
        double moist1 = SoilDep1 * (pore - thetar); // moisture above wettingfront in layer 1
        double moist2 = (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2

        double dm1 = (pore - FC1)*SoilDep1;
        double dm2 = (pore2 - FC2)*(Lw_-SoilDep1);
        //available moisture to move from layer 1 and 2
        Percolation = std::min(Percolation, dm1+dm2);
        //cannot have more flux than available water

        double perc2 = std::min(Percolation, dm2); // part taken from layer 2
        double perc1 = std::max(0.0, Percolation-perc2); // part taken from layer 1, can be 0

        moist2 -= perc2;
        // decrease WF moisture with outflow
        double m2 = (SoilDep2-Lw_)*(theta2-thetar2);// why add perc2??? + perc2;
        // m2 is moisture under Lw

        Lw_ = moist2/(pore2-thetar2) + SoilDep1;
        // new wetting front decreased with perc2 and CHin2
        theta2 = m2/(SoilDep2-Lw_) + thetar2;
        // new moisture under wetting front
        theta = pore;
        // moisture layer 1 is saturated

        if (theta2 >= pore2) {
            Lw_ = SoilDep2;
            theta2 = pore2;
        }

        // if so mucvh percolation that it goes back into layer 1
        if (perc1 > 0) {
            // percolation uses layer 1
            moist1 -= perc1;
            // sat layer 1 now becomes free until FC
            Lw_ = moist1/(pore-thetar);
            theta = FC2;
        }
    }

    Thetaeff->Drc = theta;
    ThetaI2->Drc = theta2;

    //DO NOT RECALCULATE PSI.
    // Psi1->Drc = 0.01 * 10.2 * Psia1->Drc * std::max(1.0, pow((theta-thetar)/(pore-thetar), -1.0/lambda1->Drc));
  //  Psi2->Drc = 0.01 * 10.2 * Psia2->Drc * std::max(1.0, pow((theta2-thetar2)/(pore2-thetar2), -1.0/lambda2->Drc));

    Lw->Drc = Lw_;


}
//---------------------------------------------------------------------------
void TWorld::cell_Redistribution1(int r, int c)
{
    double Lw_ = Lw->Drc;
    double L_min = 0.05; // minimum L before percolation starts

    if (SwitchImpermeable) {
        if (Lw_ > SoilDepth2->Drc-0.001)
            return;
    }
    if (Lw_ < L_min)
        return;

    double Percolation, theta_E;

    double pore = Thetaeff->Drc;
    double thetar = ThetaR1->Drc;
    double theta = Thetaeff->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double FC = ThetaFC1->Drc;

    if (Lw_ < SoilDep1-0.001) {
        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = Ksateff->Drc * pow(theta_E, 3.0+2.0/lambda1->Drc); // m/timestep
        //  Percolation = sqrt(Percolation * Ksateff->Drc);
        Percolation = Aavg(Percolation, Ksateff->Drc);
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

    Thetaeff->Drc = theta;
    //DO NOT RECALCULATE PSI
   // Psi1->Drc = 0.01 * 10.2 * Psia1->Drc * psiCalibration * std::max(1.0, pow((theta-thetar)/(pore-thetar), -1.0/lambda1->Drc));
    Lw->Drc = Lw_;
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

        if (SwitchGWflow) {
            if (GWWH->Drc > SoilDepth2->Drc-HMIN)
                return 0;
        }
        // no percolation to second layer if it is full with GW

        pore = ThetaS2->Drc;
        thetar = ThetaR2->Drc;
        theta = ThetaI2->Drc;
        double SoilDep2 = SoilDepth2->Drc;
        double ksat = factor*Ksat2->Drc;
        double FC2 = 0.7867*exp(-0.012*Ksat2->Drc)*pore;

        if(theta > thetar) {
            // percolation in m per timestep
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, 3.0+2.0/lambda2->Drc);

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
                theta = FC2;
                //double Lwo = Lw_;
                Lw_ = std::max(0.0, Lw_ - ksat/(pore - theta));
                Percolation = ksat; //(Lwo-Lw_)*(pore-theta);
            }
            ThetaI2->Drc = theta;
             //DO NOT RECALCULATE PSI
            //Psi2->Drc = 0.01 * 10.2 * Psia2->Drc * psiCalibration * std::max(1.0, pow((theta-thetar)/(pore-thetar), -1.0/lambda2->Drc));
            Lw->Drc = Lw_;
            return(Percolation);
        }
    } else {
        // one layer
        pore = Poreeff->Drc;
        thetar = ThetaR1->Drc;
        theta = Thetaeff->Drc;
        double ksat = factor*Ksateff->Drc;

        if (SwitchGWflow && GWWH->Drc > SoilDepth1->Drc-HMIN)
            return 0;


        if (theta > thetar) {
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, 3.0+2.0/lambda1->Drc);

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
                Lw_ = std::max(0.0, Lw_ - ksat/(pore - theta));
                Percolation = ksat;//(Lwo-Lw_)*(pore-theta);
            }
            Thetaeff->Drc = theta;
            //DO NOT RECALCULATE PSI
           // Psi1->Drc = 0.01 * 10.2 * Psia1->Drc * psiCalibration * std::max(1.0, pow((theta-thetar)/(pore-thetar), -1.0/lambda1->Drc));

            Lw->Drc = Lw_;
            return(Percolation);
        }
    }
    return(0);
}

//---------------------------------------------------------------------------
// aletrnative, compacter writing and possibility for three layer
// percolation is always from the lowest layer

// NOT USED OR TESTED YET

double TWorld::cell_PercolationMulti(int r, int c, double factor)
{
    double Percolation, dL, theta_E;
    double Lw_ = Lw->Drc;
    cTMap *SoilDepth = SoilDepth1;
    cTMap *pore = Thetaeff;
    cTMap *theta = ThetaI1;
    cTMap *thetar = ThetaR1;
    cTMap *ksat = Ksateff;
    cTMap *lambda = lambda1;
    cTMap *FC = ThetaFC1;
    cTMap *Psi = Psi1;
   // cTMap *Psia = Psia1;

    if(SwitchTwoLayer) {
        SoilDepth = SoilDepth2;
        pore = ThetaS2;
        theta = ThetaI2;
        thetar = ThetaR2;
        FC = ThetaFC2;
        ksat = Ksat2;
        lambda = lambda2;
        Psi = Psi2;
      //  Psia = Psia2;
    }

    if (SwitchGWflow) {
       if (GWWH->Drc > SoilDepth->Drc - HMIN)
           return 0;
    }

    if(theta->Drc > thetar->Drc) {
        // field capacity, after Saxton and Rawls 2006
        // correlation mad ein excel sheet
        double ksat_ = factor*ksat->Drc;

        // percolation in m per timestep
        theta_E = (theta-thetar)/(pore-thetar);
        Percolation = ksat_ * pow(theta_E, 3.0+2.0/lambda->Drc);

        if (SwitchThreeLayer)
            dL = SoilDepth3->Drc - std::max(SoilDepth2->Drc, Lw_);
        else
            if (SwitchTwoLayer)
                dL = SoilDepth2->Drc - std::max(SoilDepth1->Drc, Lw_);
            else
                dL = SoilDepth - Lw;
        // assumption: if Wet Fr still in first layer percolation only make 2nd drier

        if (Lw_ < SoilDepth->Drc - 0.001) {
            double moisture = dL*(theta->Drc - thetar->Drc);
            // available moisture in last layer
            Percolation = std::min(Percolation, moisture);
            moisture -= Percolation;
            theta->Drc = moisture/dL + thetar->Drc;
            // adjust theta of last layer
        } else {
            // wetting front = soildepth, dL = 0, moisture = 0
            // assume theta goes back to field capacity and decrease the wetting fornt
            theta->Drc = FC->Drc;

            Lw_ = std::max(0.0, Lw_ - ksat_/(pore->Drc - theta->Drc));
            Percolation = ksat_;
        }
        Lw->Drc = Lw_;
        return(Percolation);
    }

    return(0);
}

//---------------------------------------------------------------------------

/*!
 \brief Calculates changes in soilwater with percolation from the bottom of the profile.

  Calculates changes in soilwater with percolation from the bottom of the profile, \n
  resulting in the soil becoming dryer. Based on BrooksCorey type of percolation: \n
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
        if (SwitchTwoLayer)
        cell_Redistribution2(r, c);
        else
        cell_Redistribution1(r, c);

        Perc->Drc = cell_Percolation(r, c, 1.0);
    }}
}


void TWorld::cell_SlopeStability(int r, int c)
{

//    grad = slope(DEM)+0.005;
//    cosS = cos(atan(grad));
//    sinS = sin(atan(grad));

//    bulk_w = 9.8;
//    # bulk density water in kN/m3

//    Mu = GWDepth/1000;
//    # pore pressure in m
//    report D = soildepth/1000;
//    # soil depth in m
//    report S = (coh+(D*bulk - Mu*bulk_w)*(cosS**2)*TanPhi);
//    # shear strength
//    report T = D*bulk*sinS*cosS;


//    F = S/T;
//    #safety factor, strength/stress, F >=1 means stable
//    F = if(outcrop, 2, F);
//    # no instability on outcrops
//    report F = min(2,F);
//    #Safety Factor based on Coulomb, cut off at <= 2 for display
//    report FDays = FDays + if (F lt 1, 1, 0);
//    # cumulative days in year when unstable
//    report FdayTot = FDays;
//    # report the last timestep, cumulative unstable days

    double F = 0;
    if (CohesionSoil->Drc > 0) {
        double cosGrad_ = cosGrad->Drc;
   //  qDebug() << cosGrad_;
     double soilbulk = SoilDepth2->Drc*BulkDensity->Drc;
        double S = CohesionSoil->Drc + (soilbulk - GWWH->Drc * 1000.0)*(cosGrad_*cosGrad_)*AngleFriction->Drc; // shear strength kPa

           double T = soilbulk *Grad->Drc*cosGrad_;// shear stress kPa
      //  qDebug() << S << soilbulk << Grad->Drc << cosGrad_; //T;
        F = Grad->Drc > 0.01 ? S/T : 0.0;
    }

   FSlope->Drc = F;

}

//---------------------------------------------------------------------------
void TWorld::cell_Redistribution2(int r, int c)
{
   double Lw_ = Lw->Drc;

   if (SwitchImpermeable) {
        if (Lw_ > SoilDepth2->Drc-0.001)
        return;
   }

   double Percolation, theta_E;

   double pore = Poreeff->Drc;
   double thetar = ThetaR1->Drc;
   double theta = Thetaeff->Drc;
   double SoilDep1 = SoilDepth1->Drc;
   double FC1 = ThetaFC1->Drc;

   double pore2 = ThetaS2->Drc;
   double thetar2 = ThetaR2->Drc;
   double theta2 = ThetaI2->Drc;
   double SoilDep2 = SoilDepth2->Drc;
   double FC2 = ThetaFC2->Drc;
   double DL2 = SoilDep2-SoilDep1;

   // if Lw still in layer 1
   if (Lw_ < SoilDep1) {
        // decrease L with flow into the unsat zone beneath

        // percolation flux, avg K from L into unsat SD1
        theta_E = (theta-thetar)/(pore-thetar);
        double Percolation = Ksateff->Drc * pow(theta_E, 3.0+2.0/lambda1->Drc); // m/timestep
        Percolation = Aavg(Percolation, Ksateff->Drc);
        // harmonic mean = Havg

        double moistw = Lw_ * (pore-thetar);
        //available sat moisture above Lw_
        double dm = (pore-FC1)*Lw_;
        Percolation = std::min(dm, Percolation);
        // not more percolation than putting moisture content at field capacity

        moistw -= Percolation;
        // decrease moistw with percolation, can be 0
        Lw_ = moistw/(pore-thetar);
        // new Lw_

        double store = (SoilDep1 - Lw_) * (pore-theta); // space in SD1 under Lw_
        if (Percolation <= store) {
        // if percolation fits in store layer 1 under the Lw
        double m1 = (theta-thetar)*(SoilDep1-Lw_) + Percolation;
        theta = m1/(SoilDep1-Lw_) + thetar;

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

        m1 = (theta-thetar)*(SoilDep1-Lw_) + Perc1;
        theta = m1/(SoilDep1-Lw_)+thetar;

        // cannot happen!
        if (theta >= pore) {
                theta = pore;
                Lw_ = SoilDep1;
        }

        m2 = (theta2-thetar2)*DL2 + Perc2;
        theta2 = m2/DL2 + thetar2;

        // cannot happen!
        if (theta2 >= pore2) {
                theta2 = pore2;
                Lw_ = SoilDep2;
        }
        }

        // unsaturated flow between layer 1 and 2
        // avg percolation flux between layers
        theta_E = (theta-thetar)/(pore-thetar);
        double Perc1 = Ksateff->Drc * pow(theta_E, 3.0+2.0/lambda1->Drc); // m/timestep
        theta_E = (theta2-thetar2)/(pore2-thetar2);
        double Perc2 = Ksat2->Drc * pow(theta_E, 3.0+2.0/lambda2->Drc); // m/timestep
        Percolation = Aavg(Perc1, Perc2);

        double m1 = (SoilDep1-Lw_)*(theta-thetar);  // max moist
        Percolation = std::min(Percolation, m1);
        double m2 = DL2*(pore2-theta2); // max fit
        Percolation = std::min(Percolation, m2);

        m1 -= Percolation;
        theta = thetar + m1/(SoilDep1-Lw_);

        m2 = DL2*(theta2-thetar2) + Percolation;
        theta2 = thetar2 + m2/DL2;
        theta2 = std::min(pore2,theta2);

   } else {
        //Lw_ > SoilDep1

        theta_E = (theta2-thetar2)/(pore2-thetar2);
        Percolation = Ksat2->Drc* pow(theta_E, 3.0+2.0/lambda2->Drc); // m/timestep
        Percolation = Aavg(Percolation, Ksat2->Drc);
        double moist1 = SoilDep1 * (pore - thetar); // moisture above wettingfront in layer 1
        double moist2 = (Lw_-SoilDep1) * (pore2 - thetar2); // moisture above wettingfront in layer 2

        double dm1 = (pore - FC1)*SoilDep1;
        double dm2 = (pore2 - FC2)*(Lw_-SoilDep1);
        //available moisture to move from layer 1 and 2
        Percolation = std::min(Percolation, dm1+dm2);
        //cannot have more flux than available water

        double perc2 = std::min(Percolation, dm2); // part taken from layer 2
        double perc1 = std::max(0.0, Percolation-perc2); // part taken from layer 1, can be 0

        moist2 -= perc2;
        // decrease WF moisture with outflow
        double m2 = (SoilDep2-Lw_)*(theta2-thetar2);// why add perc2??? + perc2;
        // m2 is moisture under Lw

        Lw_ = moist2/(pore2-thetar2) + SoilDep1;
        // new wetting front decreased with perc2 and CHin2
        theta2 = m2/(SoilDep2-Lw_) + thetar2;
        // new moisture under wetting front
        theta = pore;
        // moisture layer 1 is saturated

        if (theta2 >= pore2) {
        Lw_ = SoilDep2;
        theta2 = pore2;
        }

        // if so mucvh percolation that it goes back into layer 1
        if (perc1 > 0) {
        // percolation uses layer 1
        moist1 -= perc1;
        // sat layer 1 now becomes free until FC
        Lw_ = moist1/(pore-thetar);
        theta = FC2;
        }
   }

   Thetaeff->Drc = theta;
   ThetaI2->Drc = theta2;

   //DO NOT RECALCULATE PSI.
   // Psi1->Drc = 0.01 * 10.2 * Psia1->Drc * psiCalibration * std::max(1.0, pow((theta-thetar)/(pore-thetar), -1.0/lambda1->Drc));
   //  Psi2->Drc = 0.01 * 10.2 * Psia2->Drc * psiCalibration * std::max(1.0, pow((theta2-thetar2)/(pore2-thetar2), -1.0/lambda2->Drc));

   Lw->Drc = Lw_;


}
