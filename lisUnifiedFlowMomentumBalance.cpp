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


/*!
 \file lisUnifiedFlowMomentum.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

double TWorld::UF_Friction(double dt,double velx,double vely, double NN, double h, double slope)
{

    double nsq = UF_MANNINGCOEFFICIENT * (0.1+NN)*(0.1+NN)*UF_Gravity*sqrt(velx*velx + vely*vely)*UF2D_MinimumDT/pow(std::max(0.0,h),4.0/3.0);

    double kinfac = pow(std::max(0.0,std::min(1.0,(h/0.25))),0.5);
    //double vkin = (slope  > 0? -1.0: 1.0) *sqrt(std::fabs(slope)) * pow(h,3.0/2.0)/NN;
    return velx/(1.0+ kinfac* nsq); //vkin * kinfac + (1-kinfac) *velx/(1.0+nsq);
}
/*double TWorld::UF_Friction2(double dt,double a, double velx,double vely, double NN, double h)
{
    //h = (h < 1.0)? std::min(1.0,pow(h,0.75)) : (h);
    double nsq = UF_MANNINGCOEFFICIENT * (0.1+NN)*(0.1+NN)*UF_Gravity/pow(std::max(UF_VERY_SMALL,h),4.0/3.0);
    double v = sqrt(a/nsq);
    return v;
}*/

double TWorld::UF2D_MomentumBalanceFluid(bool x, double _f,double _s,double fu, double fv, double su, double sv, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbf, double SlopeX, double SlopeY,
                                 double dhfdx,double dhfdy, double dh2pbdx, double dh2pbdy, double dsfdx, double dsfdy, double ddsfdxx, double ddsfdyy, double ddsfdxy,
                                 double dfudx, double dfudy, double dfvdx, double dfvdy, double ddfudxx, double ddfudyy, double ddfvdxy, double ddfvdxx, double ddfvdyy, double ddfudxy,
                                 double dsudx, double dsudy, double dsvdx,double dsvdy)
{
    double h = (_f + _s)/(_dx*_dx);
    if(x) {

        return (-UF_Gravity * sin(SlopeX) -//UF_Gravity *dhfdx-
                UF_Aspect *(
                    (dh2pbdx)/h
                   + (pbf * SlopeX)
                    -(1.0/(ff * Nr))*(
                        2.0*ddfudxx + ddfvdxy + ddfudyy - UF_Chi * fu/(UF_Aspect*UF_Aspect * h* h)
                        )
                    +(1.0/(ff * Nra))*(
                        2.0 *(dsfdx*(dfudx - dsudx) + ddsfdxx*(fu - su))
                        +(dsfdx*(dfvdy - dsvdy) + ddsfdxy*(fv - sv))
                        +(dsfdy*(dfudy - dsudy) + ddsfdyy*(fu - su))
                        + UF_Chi * fu/(UF_Aspect*UF_Aspect*ff*Nra*h*h)
                        )
                    -(UF_Ksi*sf*(fu - su)/(UF_Aspect*UF_Aspect*ff*Nra*h*h))
                    )
                );

    } else {

        return (-UF_Gravity * sin(SlopeY) - //UF_Gravity *dhfdy-
                UF_Aspect * (
                    (dh2pbdy)/h
                   + (pbf * SlopeY)
                    -(1.0/(ff * Nr))*(
                        2.0*ddfudyy + ddfvdxy + ddfudxx - UF_Chi * fv/(UF_Aspect*UF_Aspect * h* h)
                        )
                    +(1.0/(ff * Nra))*(
                        2.0 *(dsfdy*(dfvdy - dsvdy) + ddsfdyy*(fv - sv))
                        +(dsfdy*(dfudx - dsudx) + ddsfdxy*(fu - su))
                        +(dsfdx*(dfvdx - dsvdx) + ddsfdxx*(fv - sv))
                        + UF_Chi * fv/(UF_Aspect*UF_Aspect*ff*Nra*h*h)
                        )
                    -(UF_Ksi*sf*(fv - sv)/(UF_Aspect*UF_Aspect*ff*Nra*h*h))
                    )
                );
    }



}

double TWorld::UF2D_MomentumBalanceSolid(bool x, double _f,double _s,double fu, double fv, double su, double sv, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbs,double pbf, double SlopeX, double SlopeY,
                                 double dhsdx, double dhsdy, double dhdx, double dhdy, double dbdx, double dbdy)
{

    if(x)
    {
        return (UF_Gravity * sin(SlopeX) - (su < 0? 1.0 : -1.0)*std::tan(ifa)*pbs-UF_Aspect*pbs*SlopeX
                -UF_Aspect * gamma * pbf * ( dhdx +  SlopeX )
                 );

    }else
    {
        return (UF_Gravity * sin(SlopeY) - (sv < 0? 1.0 : -1.0)*std::tan(ifa)*pbs-UF_Aspect*pbs*SlopeY
                -UF_Aspect * gamma * pbf * ( dhdy +  SlopeY )

               );
    }

}

double TWorld::UF1D_MomentumBalanceFluid(double _f,double _s,double fu, double su, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbf, double Slope,
                                 double dhfdx, double dh2pbdx, double dsfdx, double ddsfdxx, double dfudx, double ddfudxx, double dsudx)
{
    double h = (_f + _s)/(_dx*_dx);
    return
        (UF_Gravity * sin(Slope) -
         UF_Aspect *(
             (dh2pbdx)/h
             +(pbf * Slope)
             -(1.0/(ff * Nr))*(
                 2.0*ddfudxx - UF_Chi * fu/(UF_Aspect*UF_Aspect * h* h)
                 )
             +(1.0/(ff * Nra))*(
                 2.0 *(dsfdx*(dfudx - dsudx) + ddsfdxx*(fu - su))
                 + UF_Chi * fu/(UF_Aspect*UF_Aspect*ff*Nra*h*h)
                 )
             -(UF_Ksi*sf*(fu - su)/(UF_Aspect*UF_Aspect*ff*Nra*h*h))
             )
          );


}

double TWorld::UF1D_MomentumBalanceSolid(double _f,double _s,double fu, double su, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbs,double pbf, double Slope,
                                 double dhsdx, double dhdx, double dbdx)
{

    return
        (UF_Gravity * sin(Slope) - (su > 0? 1.0 : -1.0)*std::tan(ifa)*pbs-UF_Aspect*pbs*Slope
        -UF_Aspect * gamma * pbf * ( dhdx +  Slope)
         );


}

