
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

#include "model.h"


void TWorld::Seismic()
{
    if(!SwitchSeismic)
    {
        return;
    }

    FOR_ROW_COL_MV
    {
        if(!(PGATiming->Drc > (time/ 60.0)) && PGAInitiated->Drc == 0)
        {

            double loss1 = 0;
            double loss10 = 0;

            double loss2 = 0;
            double loss20 = 0;

            double finalloss1 = 0;
            double finalloss2 = 0;

            //get the estimated loss in cohesion and internal friction angle

            double h = DFSoilDepth->Drc;
            double s = DFSlope->Drc;
            double ifa = DFSoilInternalFrictionAngle->Drc;
            double coh = DFSoilCohesion->Drc + DFPlantCohesion->Drc;
            double theta = DFSoilDepth->Drc > 0? DFWaterHeight->Drc/DFSoilDepth->Drc : 0.0;
            double apeak = PGA->Drc;
            double psi = DFWaterSuction->Drc;
            double forcingdown = DFForcing->Drc;
            double forcingup = DFForcingUp->Drc;
            double sf_cal = DFSFCalibration->Drc;
            double gamma = DFSoilDensity->Drc;
            double gammaw = 1000.0;

            loss1 = GetSeismicStrengthLoss(h,s,ifa,coh,gamma,gammaw,theta,apeak,psi,forcingdown,forcingup,sf_cal);
            loss10 = GetSeismicStrengthLoss(h,s,ifa,coh,gamma,gammaw,theta,0,psi,forcingdown,forcingup,sf_cal);

            finalloss1 = (loss1 - loss10)/(1.0-loss10);

            if(SwitchBedrock)
            {
                double h = DFSoilDepth2->Drc;
                double s = DFSlope->Drc;
                double ifa = DFSoilInternalFrictionAngle2->Drc;
                double coh = 0.0;
                double apeak = PGA->Drc;
                double psi = 0.0;
                double forcingdown = DFForcing2->Drc;
                double forcingup = DFForcingUp2->Drc;
                double sf_cal = DFSFCalibration2->Drc;
                double gamma = DFSoilDensity2->Drc;
                double gammaw = 1000.0;

                loss2 = GetSeismicStrengthLoss(h,s,ifa,coh,gamma,gammaw,theta,apeak,psi,forcingdown,forcingup,sf_cal);
                loss20 = GetSeismicStrengthLoss(h,s,ifa,coh,gamma,gammaw,theta,0,psi,forcingdown,forcingup,sf_cal);

                finalloss2 = (loss2 - loss20)/(1.0-loss20);
            }

            //re-scale the actual soil properties

            StrengthLoss->Drc = finalloss1;
            StrengthLoss->Drc = finalloss2;

            //mark this cell as being done
            PGAInitiated->Drc = 1;
            PGACurrent->Drc = PGA->Drc;

        }else
        {

            PGACurrent->Drc = 0;
            //StrengthLoss->Drc = 0;
            //if(SwitchBedrock)
            //{
            //    StrengthLoss2->Drc = 0;
            //}


        }

    }

}

double TWorld::GetSeismicStrengthLoss(double H, double s, double ifa, double c, double gamma, double gammaw, double theta, double apeak, double psi, double forcingdown, double forcingup, double sf_cal)
{
    double a = forcingup + gamma*(1 - theta)*std::pow(std::cos(s),2.0) -
                 apeak*( gamma*(1 - theta) + gammaw*theta)*std::cos(s)*std::sin(s)
                +  gammaw*theta*std::pow(std::cos(s),2.0);
    double b = (forcingdown + apeak*(gamma + gammaw*theta) *  std::pow(std::cos(s),2.0) +
                std::cos(s)*std::sin(s) *(gamma - gammaw*theta));

    double loss = (-1/(H*2*(b - a*std::tan(ifa))))*
         std::exp(-((4.0*std::pow(std::atan(std::sqrt(1.0 + (a*a)/(b*b)) - a/b),2.0))/std::pow(ifa,2.0))) +
         (-1 + std::exp(4.0*std::pow(std::atan(std::sqrt(1.0 + (a*a)/(b*b)) - a/b),2.0))/std::pow(ifa,2.0))
         *(-2*b*H +c*std::sqrt(3.14)*std::erf((H*(b - a*std::tan(ifa)))/c) +2*a*H*std::tan(ifa));

    return std::max(std::min(loss,1.0),0.0);

}
