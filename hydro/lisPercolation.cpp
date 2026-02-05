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

#include "lisemqt.h"
#include "model.h"

// calc average soil moisture content for output to screen and folder
void TWorld::avgTheta()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Lw_ = Lw->Drc;
        double SoilDep1 = SoilDepth1->Drc;
        ThetaI1a->Drc = Thetaeff->Drc;

        if (Lw_ > 0 && Lw_ < SoilDep1 - 1e-3) {
            double f = Lw_/SoilDep1;
            ThetaI1a->Drc = f * Poreeff->Drc + (1-f) *Thetaeff->Drc;
        }
        if (Lw_ > SoilDep1 - 1e-3)
            ThetaI1a->Drc = Poreeff->Drc;

        if (SwitchTwoLayer) {
            double SoilDep2 = SoilDepth2->Drc;
            ThetaI2a->Drc = ThetaI2->Drc;
            if (Lw_ > SoilDep1 && Lw_ < SoilDep2 - 1e-3) {
                double f = (Lw_-SoilDep1)/(SoilDep2-SoilDep1);
                ThetaI2a->Drc = f * ThetaS2->Drc + (1-f) *ThetaI2->Drc;
            }
            if (Lw_ > SoilDep2 - 1e-3)
                ThetaI2a->Drc = ThetaS2->Drc;
        }
    }}
}
//---------------------------------------------------------------------------
//CURRENTLY NOT USED
double TWorld::SoilWaterMass()
{
    double totsatm3 = 0;
    double totunsatm3 = 0;

    FOR_ROW_COL_MV_L {
        double totsat = 0;
        double totunsat = 0;
        if (SwitchTwoLayer) {
            if (Lw->Drc <= SoilDepth1->Drc) {
                totsat = totsat + Lw->Drc * Poreeff->Drc; //layer 1
                totunsat = totunsat + (SoilDepth1->Drc - Lw->Drc) * Thetaeff->Drc; //layer 1
                totunsat = totunsat + (SoilDepth2->Drc - SoilDepth1->Drc) * ThetaI2->Drc; //layer 2
            } else {
                totsat = totsat + SoilDepth1->Drc * Poreeff->Drc; //layer 1
                totsat = totsat + (Lw->Drc-SoilDepth1->Drc) * ThetaS2->Drc; //layer 2
                totunsat = totunsat + (SoilDepth1->Drc - SoilDepth2->Drc) * ThetaI2->Drc; //layer 2
            }
        } else {
            totsat = totsat + Lw->Drc * Poreeff->Drc;
            totunsat = totunsat + (SoilDepth1->Drc - Lw->Drc) * Thetaeff->Drc;
        }
        totsatm3 += totsat * CHAdjDX->Drc;
        totunsatm3 += totunsat * CHAdjDX->Drc;
    }}

    return totsatm3+totunsatm3;
}
//---------------------------------------------------------------------------
// percolation is always from the lowest layer
// calc Percolation based on Ksat, compare with available moisture in the lowest layer
// adjust Lw and theta
// factor is from calibration factor from GW recharge, else if no GW it is 1.0
double TWorld::cell_PercolationMulti(int r, int c, double factor)
{
    double Lw_ = Lw->Drc;
    double Percolation;

    // 1 layer
    double pore = Poreeff->Drc;
    double thetar = ThetaR1->Drc;
    double theta = Thetaeff->Drc;
    double SoilDep = SoilDepth1->Drc;
    double SoilDepa = 0;
    double FC = ThetaFC1->Drc;
    double Ksat = Ksateff->Drc*factor;
    double lambda = lambda1->Drc;

    if (SwitchTwoLayer) {
        pore = ThetaS2->Drc;
        thetar = ThetaR2->Drc;
        theta = ThetaI2->Drc;
        SoilDep = SoilDepth2->Drc;
        SoilDepa = SoilDepth1->Drc;
        FC = ThetaFC2->Drc;
        Ksat = Ksat2->Drc*factor;
        lambda = lambda2->Drc;
    }

    if(SwitchThreeLayer) {
        pore = ThetaS3->Drc;
        thetar = ThetaR3->Drc;
        theta = ThetaI3->Drc;
        SoilDep = SoilDepth3->Drc;
        SoilDepa = SoilDepth2->Drc;
        FC = ThetaFC3->Drc;
        Ksat = Ksat3->Drc*factor;
        lambda = lambda3->Drc;
    }

    if (SwitchGWflow) {
        if (GWWH->Drc > SoilDep - HMIN)
            return 0;
        // soil is full with GW, no percolation
    }

    if(theta > thetar) {

        // percolation in m per timestepbaased on ksat lowest layer
        double theta_E = (theta-thetar)/(pore-thetar);
        Percolation = Ksat * pow(theta_E, 3.0+2.0/lambda);

        // calculate max amount of moisture available for percolation
        // moisture in last layer dL or less if Wettubg front is in last layer
        if (Lw_ < SoilDep - 0.001) {
            double dL = SoilDep - qMax(SoilDepa, Lw_);
            double moisture = dL*(theta - thetar);
            // available moisture in last layer
            Percolation = qMin(Percolation, moisture);
            moisture -= Percolation;
            theta = moisture/dL + thetar;
            // adjust theta of last layer
        } else {
            // wetting front = soildepth, dL = 0, moisture = 0
            // assume theta goes back to field capacity and decrease the wetting fornt
            // assume percolation is Ksat
            theta = FC;
            Lw_ = qMax(0.0, Lw_ - Ksat/(pore - theta));
            Percolation = Ksat;
        }

        if (SwitchThreeLayer)
            ThetaI3->Drc = theta;
        else
            if (SwitchTwoLayer)
                ThetaI2->Drc = theta;
            else
                Thetaeff->Drc = theta;

        Lw->Drc = Lw_;
        return(Percolation);
    }

    return(0);
}
//---------------------------------------------------------------------------
// calculates flux from wetting front in a layer to underlying unsat zone and adjusts theta and Lw
// function is only called when Lw is in the layer (1,2 or 3)
void TWorld::adjustLWTheta(int r, int c, double SoilDepAbove, cTMap *Ksat, cTMap *pore, cTMap *theta, cTMap *thetar, cTMap *FC, cTMap *SoilDep, cTMap *lambda)
{
    double Percolation = Ksat->Drc * pow((theta->Drc-thetar->Drc)/(pore->Drc-thetar->Drc), 3.0+2.0/lambda->Drc); // m/timestep
    Percolation = ARITHavg(Percolation, Ksat->Drc);

    // available sat moisture above Lw_
    double moist = (pore->Drc - FC->Drc)*(Lw->Drc-SoilDepAbove);
    // max that can move assuming the freed space goes to FC
    Percolation = qMin(moist, Percolation);
    // space in SD1 under Lw_
    double store = (SoilDep->Drc - Lw->Drc) * (pore->Drc-theta->Drc);
    // not more than fits into SoilDep1-Lw
    Percolation = qMin(store, Percolation);

    double moisture = qMax(0.0, Lw->Drc * (pore->Drc - thetar->Drc) - Percolation);
    //Lw_ = qMax(0.0,moisture-Percolation)/(pore-thetar);
    Lw->Drc = moisture/(pore->Drc-thetar->Drc);
    // new Lw_
    theta->Drc = qBound(thetar->Drc, theta->Drc + Percolation/(SoilDep->Drc-Lw->Drc), pore->Drc);
    // increase moisture under Lw
}
//---------------------------------------------------------------------------
// water flowing from wetting front into underlying zone, Lw decreases, theta increases
void TWorld::cell_Redistribution1(int r, int c)
{
    if (Lw->Drc == 0)
        return;
    // nothing to redistribute

    if (WH->Drc > he_ca)
        return;
    // no redistribution while infiltration because fluctuations

    if (SwitchImpermeable) {
        if (Lw->Drc > SoilDepth1->Drc-0.001)
            return;
    }
    // profile full

    if (Ksateff->Drc == 0 || Poreeff->Drc == 0)
        return;
    // impermeable

    double Lwmin = qMin(0.1,SoilDepth1->Drc/10);
    // only redistribute if the Lw is advanced a bit into the layer to avoid spurious fluctuations
    if (Lw_ > Lwmin ) {

        adjustLWTheta(r, c, 0.0, Ksateff, Poreeff, Thetaeff, ThetaR1, ThetaFC1, SoilDepth1, lambda1);
        // function calculates how much water flows from wetting zone to underlying unsat zone and adjusts Lw and Theta

    }

}
//---------------------------------------------------------------------------
//water moving from wetting front into underlying unsat zone, in layer 1 or in layer 2
void TWorld::cell_Redistribution2(int r, int c)
{
    if (Lw->Drc == 0)
        return;
    // nothing to redistribute

    if (WH->Drc > he_ca)
        return;
    // no redistribution while infiltration because fluctuations

    if (SwitchImpermeable) {
        if (Lw->Drc > SoilDepth2->Drc-0.001)
            return;
    }
    // profile full

    if (Ksateff->Drc == 0 || Poreeff->Drc == 0)
        return;
    // impermeable

    // if Lw still in layer 1
    if (Lw->Drc < SoilDepth1->Drc) {
        // flow from the saturated zone into the unsat layer below. but only the remainder of layer 1
        double Lwmin = qMin(0.1,SoilDepth1->Drc/10);
        // only redistribute if the Lw is advanced a bit into the layer to avoid spurious fluctuations
        if (Lw->Drc > Lwmin ) {

            adjustLWTheta(r,c, 0.0, Ksateff, Poreeff, Thetaeff, ThetaR1, ThetaFC1, SoilDepth1, lambda1);

        }
    } else {
        // [3] flow from wetting zone into unsat below wetting zone in layer 2
        // only if infil procvess stopped
        // consider ONLY layer 2
        double Lwmin = SoilDepth1->Drc + qMin(0.1,(SoilDepth2->Drc-SoilDepth1->Drc)/10);
        // only redistribute if the Lw is advanced a bit into the layer to avoid spurious fluctuations
        if (Lw->Drc > Lwmin ) {

            adjustLWTheta(r, c, SoilDepth1->Drc, Ksat2, ThetaS2, ThetaI2, ThetaR2, ThetaFC2, SoilDepth2, lambda2);

        }
    }

    // Thetaeff->Drc = theta;
    // ThetaI2->Drc = theta2;
    // Lw->Drc = Lw_;
}
//---------------------------------------------------------------------------
// redistribution of water in the wetting front to the unsat zone below
// only if infiltration process has stopped to avoid fluctuations and spurious behaviour
void TWorld::cell_Redistribution3(int r, int c)
{
    if (Lw->Drc == 0)
        return;
    // nothing to redistribute

    if (WH->Drc > he_ca)
        return;
    // no redistribution while infiltration because fluctuations

    if (SwitchImpermeable) {
        if (Lw->Drc > SoilDepth3->Drc-0.001)
            return;
    }
    // profile full

    if (Ksateff->Drc == 0 || Poreeff->Drc == 0)
        return;
    // impermeable

    // if Lw still in layer 1
    if (Lw->Drc < SoilDepth1->Drc) {
        double Lwmin = qMin(0.1,SoilDepth1->Drc/10);
        // only redistribute if the Lw is advanced a bit into the layer to avoid spurious fluctuations
        if (Lw->Drc > Lwmin ) {

            adjustLWTheta(r,c, 0.0, Ksateff, Poreeff, Thetaeff, ThetaR1, ThetaFC1, SoilDepth1, lambda1);

        }
    } else
    if (Lw->Drc < SoilDepth2->Drc) {
        // [3] flow from wetting zone into unsat below wetting zone in layer 2
        // consider ONLY layer 2
        double Lwmin = SoilDepth1->Drc + qMin(0.1,(SoilDepth2->Drc-SoilDepth1->Drc)/10);
        // only redistribute if the Lw is advanced a bit into the layer to avoid spurious fluctuations
        if (Lw->Drc > Lwmin ) {

            adjustLWTheta(r, c, SoilDepth1->Drc, Ksat2, ThetaS2, ThetaI2, ThetaR2, ThetaFC2, SoilDepth2, lambda2);

       }
    } else
        if (Lw->Drc < SoilDepth3->Drc) {
            // [4] flow from wetting zone into unsat below wetting zone in layer 3
            //consider ONLY layer 3
            double Lwmin = SoilDepth1->Drc + qMin(0.1,(SoilDepth2->Drc-SoilDepth1->Drc)/10);
            // only redistribute if the Lw is advanced a bit into the layer to avoid spurious fluctuations
            if (Lw->Drc > Lwmin ) {

                adjustLWTheta(r, c, SoilDepth2->Drc, Ksat3, ThetaS3, ThetaI3, ThetaR3, ThetaFC3, SoilDepth3, lambda3);

           }
        } else {
            // profile is completely saturated
            Lw->Drc = SoilDepth3->Drc;
        }

}
//---------------------------------------------------------------------------
//unsaturated flow between layers, 2 or 3 soil layers
void TWorld::cell_RedistributionUnsat(int r, int c)
{
    double Lw_ = Lw->Drc;
    double Percolation = 0;
    double pore = Poreeff->Drc;
    double thetar = ThetaR1->Drc;
    double theta = Thetaeff->Drc;
    double SoilDep1 = SoilDepth1->Drc;

    double pore2 = ThetaS2->Drc;
    double thetar2 = ThetaR2->Drc;
    double theta2 = ThetaI2->Drc;
    double SoilDep2 = SoilDepth2->Drc;
    double DL2 = SoilDep2-SoilDep1;

    // if Lw still in layer 1
    // [1] unsaturated flow between layer 1 and 2
    if (Lw_ < SoilDep1 && theta > thetar && theta2 < pore2-0.001) {
        // if there is room in layer 2 and layer 1 is not too dry
        // avg percolation flux between layers, theta1 decreases, theta2 increases
        double Perc1 = Ksateff->Drc * pow((theta-thetar)/(pore-thetar),   3.0+2.0/lambda1->Drc); // m/timestep
        double Perc2 = Ksat2->Drc * pow((theta2-thetar2)/(pore2-thetar2), 3.0+2.0/lambda2->Drc); // m/timestep
        Percolation = ARITHavg(Perc1, Perc2);

        double moist1 = (SoilDep1-Lw_)*(theta-thetar);  // max moist, if Lw_ = SoilDep1 than m1 = 0
        Percolation = qMin(Percolation, moist1);
        double moist2 = DL2*(pore2-theta2); // max fit
        Percolation = qMin(Percolation, moist2);
        if (Percolation > 0) {
            moist1 -= Percolation;
            moist2 += Percolation;
            theta = thetar + moist1/(SoilDep1-Lw_);
            theta2 = qMin(pore2, moist2/DL2);
        }
    }
    Thetaeff->Drc = theta;
    ThetaI2->Drc = theta2;

    if (SwitchThreeLayer) {
        double pore3 = ThetaS3->Drc;
        double thetar3 = ThetaR3->Drc;
        double theta3 = ThetaI3->Drc;
        double SoilDep3 = SoilDepth3->Drc;
        double DL3 = SoilDep3-SoilDep2;
        if(Lw_ < SoilDep2 && theta2 > thetar2 && theta3 < pore3-0.001) {
            // if there is room in layer 3 and layer 2 is not too dry
            // avg percolation flux between layers, theta2 decreases, theta3 increases
            double Perc2 = Ksateff->Drc * pow((theta2-thetar2)/(pore2-thetar2), 3.0+2.0/lambda2->Drc); // m/timestep
            double Perc3 = Ksat2->Drc * pow((theta3-thetar3)/(pore3-thetar3), 3.0+2.0/lambda3->Drc); // m/timestep
            Percolation = ARITHavg(Perc2, Perc3);

            double moist2 = DL2*(pore2-theta2); // max fit
            Percolation = qMin(Percolation, moist2);
            double moist3 = DL3*(pore3-theta3); // max fit
            Percolation = qMin(Percolation, moist3);
            if (Percolation > 0) {
                moist2 -= Percolation;
                moist3 += Percolation;
                theta2 = thetar2 + moist2/DL2;
                theta3 = qMin(pore3, moist3/DL3);
            }
        }
        ThetaI3->Drc = theta3;
    }
}
//---------------------------------------------------------------------------
void TWorld::cell_Tiledrain1(int r, int c)
{
    if (!SwitchIncludeTile)
        return;
    if (Lw->Drc < 0.05)
        return;

    double Lw_ = Lw->Drc;
    double pore = Poreeff->Drc;
    double thetar = ThetaR1->Drc;
    double theta = Thetaeff->Drc;
    double SoilDep1 = SoilDepth1->Drc;
    double FC = ThetaFC1->Drc;

    if (Lw_ > TileDepth->Drc) {
        double vol = DX->Drc*Ksateff->Drc*TileDiameter->Drc;
        // volume draining, assuming full saturation so Ksat is draining
        double volsoil = Lw_*(pore-FC)*CHAdjDX->Drc;
        // available volume, not drier than FC, gravity
        vol = qMin(volsoil, vol);
        double tiledm = vol/CHAdjDX->Drc; // removal in m

        double moisture = Lw_*(pore-thetar); //available sat moisture above Lw_
        moisture -= tiledm; // okay because removal limited to FC
        double newLw_ = moisture/(pore-thetar); // new Lw_

        theta = (Lw_-newLw_)*FC + (SoilDep1-Lw_)*theta;
        // new moisture content is weighed avg

        Thetaeff->Drc = theta;
        Lw->Drc = newLw_;
        TileWaterVolSoil->Drc = vol;
    }
}
//---------------------------------------------------------------------------
void TWorld::cell_Tiledrain2(int r, int c)
{
    if (!SwitchIncludeTile)
        return;
    if (Lw->Drc < 0.05)
        return;

    double Lw_ = Lw->Drc;

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

    if (Lw_ > TileDepth->Drc && TileDepth->Drc <= SoilDep1) {
        double vol = DX->Drc*Ksateff->Drc*TileDiameter->Drc;
        // volume draining, assuming full saturation so Ksat is draining
        double volsoil = Lw_*(pore-FC1)*CHAdjDX->Drc;
        // available volume, not drier than FC, gravity
        vol = qMin(volsoil, vol);
        double tiledm = vol/CHAdjDX->Drc; // removal in m

        double moisture = Lw_*(pore-thetar); //available sat moisture above Lw_
        moisture -= tiledm; // okay because removal limited to FC
        double newLw_ = moisture/(pore-thetar); // new Lw_

        theta = (Lw_-newLw_)*FC1 + (SoilDep1-Lw_)*theta;
        // new moisture content is weighed avg

        Thetaeff->Drc = theta;
        Lw->Drc = newLw_;
        TileWaterVolSoil->Drc = vol;
    }

    if (Lw_ > TileDepth->Drc && TileDepth->Drc > SoilDep1) {
        double vol = DX->Drc*Ksat2->Drc*TileDiameter->Drc;
        // volume draining, assuming full saturation so Ksat is draining
        double volsoil = (Lw_-SoilDep1)*(pore2-FC2)*CHAdjDX->Drc;
        // available volume, not drier than FC, gravity
        vol = qMin(volsoil, vol);
        double tiledm = vol/CHAdjDX->Drc; // removal in m

        double moisture = (Lw_-SoilDep1)*(pore2-thetar2); //available sat moisture above Lw_
        moisture -= tiledm; // okay because removal limited to FC
        double newLw_ = moisture/(pore2-thetar2); // new Lw_

        theta2 = (Lw_-newLw_)*FC2 + (SoilDep2-Lw_)*theta2;
        // new moisture content is weighed avg

        ThetaI2->Drc = theta2;
        Lw->Drc = newLw_;
        TileWaterVolSoil->Drc = vol;
    }
}

//---------------------------------------------------------------------------

void TWorld::cell_Channelinfow1(int r, int c)
{
    /*
   ChannelQSide->Drc = 0.0;

   //    if (ChannelWH->Drc > ChannelDepth->Drc - 0.05)
   //        return;

   bool doUnsat = false;

   if (Lw->Drc < 0.01)
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

   double h = qMin(ChannelDep,Lw_);

   double moist = Lw_*(pore-thetar);
   double dh = CHin1 * h*DX_/CHAdjDX->Drc * h/dL; // ks*cross section /cellsurface * Darcy pressure
   dh = qMin(dh, moist);
   moist -= dh;
   Lw_ = moist/(pore-thetar); // new Lw

   //    double h2 = qMax(0.0, ChannelDep-Lw_);
   //    if (doUnsat && theta > 0.95*pore && h2 > 0.01) {
   //        moist = h2*(theta-thetar);
   //        dh = CHin2*h2*DX_/CHAdjDX->Drc * h2/dL;
   //        dh = qMin(dh, moist);
   //        theta = thetar + moist/h2;
   //    } else
   //        CHin2 = 0;

   ChannelQSide->Drc = DX_*(CHin1*h*h/dL);// + CHin2*h2*h2/dL); // m3
   */
}
//---------------------------------------------------------------------------
// Side inflow into channel from saturated part of the soil (Lw_), causes decrease of Lw_
// the assumption is thart the Darcy flow pressure difference dH/dL is 1.0
// afactor 2.0 is applied to Ksat because the flow is from both sides
void TWorld::cell_Channelinfow2(int r, int c)
{
    /*
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
        double h = qMin(ChannelDep,Lw_);
        if (Lw_ > 0.01) {
            double moist = h*(pore-thetar);
            double frac = (h*DX_)/CHAdjDX->Drc; //= Lw_/ChannelAdj->Drc;
            double pressgrad = (h/dL);
            double dh = CHin1 * frac * pressgrad;
            dh = qMin(dh, moist);
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
                dh = qMin(dh, moist);
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
            dh = qMin(dh, moist1);
            moist1 -= dh;
            L = moist1/(pore-thetar);
            CHin1 = dh/frac/pressgrad;

            // layer 2 saturated part, but not deeper than chandep
            double h2 = qMax(0.0,Lw_-SoilDep1);
            h2 = qMin(h2,ChannelDep-SoilDep1);
            if (h2 > 0.001) {
                double moist2 = h2*(pore2-thetar2);
                frac = (h2*DX_)/CHAdjDX->Drc; // Lw_/ChannelAdj->Drc;
                pressgrad = (h2/dL);
                dh = CHin2 * frac * pressgrad;
                dh = qMin(dh, moist2);
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
*/
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

   // double F = 0;
   // if (CohesionSoil->Drc > 0) {
   //      double cosGrad_ = cosGrad->Drc;
   //      //  qDebug() << cosGrad_;
   //      double soilbulk = SoilDepth2->Drc*BulkDensity->Drc;
   //      double S = CohesionSoil->Drc + (soilbulk - GWWH->Drc * 1000.0)*(cosGrad_*cosGrad_)*AngleFriction->Drc; // shear strength kPa

   //      double T = soilbulk *Grad->Drc*cosGrad_;// shear stress kPa
   //      //  qDebug() << S << soilbulk << Grad->Drc << cosGrad_; //T;
   //      F = Grad->Drc > 0.01 ? S/T : 0.0;
   // }

   // FSlope->Drc = F;

}

//---------------------------------------------------------------------------
// percolation from the bottom of the soil profile
// factor is for use of GW recharge

// NO LONGER USED
double TWorld::cell_Percolation(int r, int c, double factor)
{
 /*
  *    double Percolation, dL, pore, theta, thetar, theta_E;
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
        double FC2 = ThetaFC2->Drc;//0.7867*exp(-0.012*Ksat2->Drc)*pore;

        if(theta > thetar) {
            // percolation in m per timestep, assume it equals kunsat Brooks Corey
            theta_E = (theta-thetar)/(pore-thetar);
            Percolation = ksat * pow(theta_E, 3.0+2.0/lambda2->Drc);

            if (Lw_ > SoilDep1)
                dL = SoilDep2 - Lw_;
            else
                dL = SoilDep2 - SoilDep1;
            // assumption: if Wet Fr still in first layer percolation only make 2nd drier

            if (Lw_ < SoilDep2-0.001) {
                // decrease theta because of percolation
                double moisture = dL*(theta-thetar); // unsat moisture
                Percolation = qMin(Percolation, moisture);
                moisture -= Percolation;
                theta = moisture/dL + thetar;
            } else {
                // wetting front = soildepth2, dL = 0, moisture = 0
                // assume theta goes back to FC2 and decrease the wetting fornt
                theta = FC2;
                //double Lwo = Lw_;
                Percolation = ksat;
                Lw_ = qMax(0.0, Lw_ - Percolation/(pore - theta));
            }
            ThetaI2->Drc = theta;
            if (std::isnan(ThetaI2->Drc)) {
                 qDebug() << "nan" << FC2 << thetar;
            }
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
                Percolation = qMin(Percolation, moisture);
                moisture -= Percolation;
                theta = moisture/(SoilDep1 - Lw_) + thetar;
            } else {
                // wetting front = soildepth1, dL = 0, moisture = 0
                // assume theta goes back to FC and decrease the wetting fornt
                theta = ThetaFC1->Drc;
                Percolation = ksat;
                Lw_ = qMax(0.0, Lw_ - Percolation/(pore - theta));
            }

            Thetaeff->Drc = theta;
            Lw->Drc = Lw_;
            return(Percolation);
        }
    }
    return(0);
    */
}
