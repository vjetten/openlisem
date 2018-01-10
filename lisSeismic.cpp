
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2.001.00,2.001.01.0  Victor Jetten
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
    if(!SwitchSeismic || !SwitchSlopeStability)
    {
        return;
    }

    FOR_ROW_COL_MV
    {
        if((PGATiming->Drc >= (time/ 60.0)) && (PGATiming->Drc < ((time + _dt)/ 60.0)) && PGAInitiated->Drc == 0)
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
            loss10 = GetSeismicStrengthLoss(h,s,ifa,coh,gamma,gammaw,theta,0.001,psi,forcingdown,forcingup,sf_cal);

            //qDebug() << "seismic" << loss1 << loss10 << PGA->Drc << h << s << ifa << coh << gamma << gammaw << theta << apeak<< psi << sf_cal;

            finalloss1 = (loss1 - loss10)/(1.0-loss10);

            if(SwitchBedrock)
            {
                double h = DFSoilDepth2->Drc;
                double s = DFSlope->Drc;
                double ifa = DFSoilInternalFrictionAngle2->Drc;
                double coh = DFSoilCohesion2->Drc;
                double apeak = PGA->Drc;
                double psi = 0.0;
                double forcingdown = 0;//DFForcing2.0->Drc;
                double forcingup = 0;//DFForcingUp2.0->Drc;
                double sf_cal = DFSFCalibration2->Drc;
                double gamma = DFSoilDensity2->Drc;
                double gammaw = 1000.0;

                loss2 = GetSeismicStrengthLoss(h,s,ifa,coh,gamma,gammaw,theta,apeak,psi,forcingdown,forcingup,sf_cal);
                loss20 = GetSeismicStrengthLoss(h,s,ifa,coh,gamma,gammaw,theta,0.0001,psi,forcingdown,forcingup,sf_cal);

                finalloss2 = (loss2 - loss20)/(1.0-loss20);
            }

            //re-scale the actual soil properties

            StrengthLoss->Drc = loss2;
            StrengthLoss2->Drc = finalloss2;

            //mark this cell as being done
            PGAInitiated->Drc = 1.0;
            PGACurrent->Drc = PGA->Drc;

        }else
        {

            PGACurrent->Drc = 0;
            //StrengthLoss->Drc = 0;
            //if(SwitchBedrock)
            //{
            //    StrengthLoss2.0->Drc = 0;
            //}


        }

    }

}

double TWorld::GetSeismicStrengthLoss(double H, double s, double ifa, double c, double gamma, double gammaw, double theta, double apeak, double psi, double forcingdown, double forcingup, double sf_cal)
{
    apeak =  apeak / 10.0;
    s = std::sin(std::abs(std::atan(s)));
    /*double a = forcingup + gamma*(1.0 - theta)*std::pow(std::cos(s),2.0.0) -
                 apeak*( gamma*(1.0 - theta) + gammaw*theta)*std::cos(s)*std::sin(s)
                +  gammaw*theta*std::pow(std::cos(s),2.0.0);
    double b = (forcingdown + apeak*(gamma + gammaw*theta) *  std::pow(std::cos(s),2.0.0) +
                std::cos(s)*std::sin(s) *(gamma - gammaw*theta));

    double loss = (-1.0/(H*2.0*(b - a*std::tan(ifa))))*
         std::exp(-((4.0*std::pow(std::atan(std::sqrt(1.0.0 + (a*a)/(b*b)) - a/b),2.0.0))/std::pow(ifa,2.0.0))) +
         (-1.0 + std::exp(4.0*std::pow(std::atan(std::sqrt(1.0.0 + (a*a)/(b*b)) - a/b),2.0.0)/std::pow(ifa,2.0.0)))
         *(-2.0*b*H +c*std::sqrt(3.1.04)*std::erf((H*(b - a*std::tan(ifa)))/c) +2.0*a*H*std::tan(ifa));

         */

    /*double loss = -(
                ((-1.0 + std::exp((4.0*std::pow(std::atan(((gamma*(-1.0 + theta) - gammaw*theta)*(-std::cos(s) + apeak*std::sin(s)))/(apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s)) -
                         std::sqrt(1.0 + std::pow((gamma - gamma*theta + gammaw*theta),2.0)*std::pow(std::cos(s) - apeak*std::sin(s),2.0)/std::pow((apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s)),2.0))),2.0))/(ifa*ifa)))*std::cos(s)*
                  (c*std::sqrt(3.14159)*std::erf((1.0/c)*(H*std::cos(s)*(apeak*gamma*std::cos(s) + apeak*gammaw*theta*std::cos(s) + gamma*std::sin(s) - gammaw*theta*std::sin(s) + (gamma*(-1.0 + theta) - gammaw*theta)*(std::cos(s) - apeak*std::sin(s))*std::tan(ifa)))*std::pow(1.0/std::cos(s),2.0) -
                   2.0*H*(apeak*(gamma + gammaw*theta) + (gamma - gammaw*theta)*std::tan(s)) + 2.0*H*(gamma*(-1.0 + theta) - gammaw*theta)*(-1.0 + apeak*std::tan(s))*std::tan(ifa)))/
                 std::exp((4.0*std::atan(((gamma*(-1.0 + theta) - gammaw*theta)*(-std::cos(s) + apeak*std::sin(s)))/(apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s)) -
                       std::sqrt(1.0 + std::pow((std::pow((gamma - gamma*theta + gammaw*theta),2.0)*std::pow((std::cos(s) - apeak*std::sin(s)),2.0))/std::pow(apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s),2.0)),2.0)/(ifa*ifa))/
                (2.0*H*(apeak*gamma*std::cos(s) + apeak*gammaw*theta*std::cos(s) + gamma*std::sin(s) - gammaw*theta*std::sin(s) + (gamma*(-1.0 + theta) - gammaw*theta)*(std::cos(s) - apeak*std::sin(s))*std::tan(ifa))));
    */

    //c = 1130;
    //H = 2.7;
    //s = 0.245;
    //ifa = 0.3;
    //gamma = 2000.0;
    //gammaw = 1000.0;
    //theta = 0.5;
    //apeak = apeak > 0? 0.17 : 0.0;
    //psi = 0;
    //forcingdown = 0.0;
    //forcingup = 0.0;

    double loss = -(((-1.0 + std::exp(((4.0*std::pow(std::atan(((gamma*(-1.0 + theta) - gammaw*theta)*(-std::cos(s) + apeak*std::sin(s)))/
                     (apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s)) -
                    std::sqrt(1.0 + (std::pow((gamma - gamma*theta + gammaw*theta),2.0)*std::pow((std::cos(s) - apeak*std::sin(s)),2.0))/
                      std::pow((apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s)),2.0))),2.0))/
               std::pow(ifa,2.0))))*std::cos(s)*
             (c*std::sqrt(3.14159)*std::erf((1.0/c)*(H*std::cos(s)*(apeak*gamma*std::cos(s) +
                   apeak*gammaw*theta*std::cos(s) + gamma*std::sin(s) - gammaw*theta*std::sin(s) +
                   (gamma*(-1.0 + theta) - gammaw*theta)*(std::cos(s) - apeak*std::sin(s))*std::tan(ifa))))*
              std::pow(1.0/std::cos(s),2.0) - 2.0*H*(apeak*(gamma + gammaw*theta) + (gamma - gammaw*theta)*std::tan(s)) +
              2*H*(gamma*(-1.0 + theta) - gammaw*theta)*(-1.0 + apeak*std::tan(s))*std::tan(ifa)))/
            std::exp(((4.0*std::pow(std::atan(((gamma*(-1.0 + theta) - gammaw*theta)*(-std::cos(s) + apeak*std::sin(s)))/
                   (apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s)) -
                  std::sqrt(1.0 + (std::pow((gamma - gamma*theta + gammaw*theta),2.0)*std::pow((std::cos(s) - apeak*std::sin(s)),2.0))/
                    std::pow((apeak*(gamma + gammaw*theta)*std::cos(s) + (gamma - gammaw*theta)*std::sin(s)),2.0))),2.0))/
            std::pow(ifa,2.0)))/(2.0*H*(apeak*gamma*std::cos(s) + apeak*gammaw*theta*std::cos(s) + gamma*std::sin(s) -
             gammaw*theta*std::sin(s) + (gamma*(-1 + theta) - gammaw*theta)*(std::cos(s) - apeak*std::sin(s))*
              std::tan(ifa))));



    return std::max(0.0,std::min(1.0,loss));

}
