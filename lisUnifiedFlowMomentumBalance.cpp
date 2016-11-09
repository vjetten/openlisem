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

double TWorld::UF_Friction(double a,double dt,double velx,double vely, double NN, double h, double slope, bool solid, bool channel, double flowwidth)
{

    /*if(channel)
    {
            velx += a;
            double manning = solid?UF_MANNINGCOEFFICIENT_SOLID:UF_MANNINGCOEFFICIENT_FLUID;
            double nsq = (0.1+NN)*(0.1+NN)*UF_Gravity*sqrt(std::fabs(velx*velx))*pow(dt,1.0) /pow(std::max(UF_VERY_SMALL,h),4.0/3.0);

            if(channel)
            {
                if(flowwidth > 0)
                {
                    nsq = nsq * (flowwidth + 0.5 * h)/(flowwidth);
                }
            }
            double kinfac = 0.5 +  0.5 * pow(std::max(0.0,std::min(1.0,(h/0.25))),2.0);

            return (velx) /(1.0+ kinfac*nsq);
    }*/

    double velo = velx;
    double signa = a>0?1.0:-1.0;
    a = std::min(std::fabs(a)/dt,5.0 * h);

    double nsq = (0.1+NN)*(0.1+NN)*UF_Gravity/pow(std::max(UF_VERY_SMALL,h),4.0/3.0);

    if(channel)
    {
        if(flowwidth > 0)
        {
            //nsq = nsq * (flowwidth + h/4.0)/(flowwidth);
        }
    }
    double kinfac = std::max(0.05,(0.5 +  0.5 * pow(std::max(0.0,std::min(1.0,(h/0.25))),2.0)));
    velx = sqrt(a)/sqrt(kinfac *nsq); //-nsq * dt *a + sqrt(nsq)*sqrt(std::fabs(a))*sqrt(4 + a * dt*dt*nsq)/(2.0*nsq);
    velx = signa *velx ;

    double fac = exp(-dt / std::max(20.0,20.0 * h));
    velx = fac * velo + (1.0-fac) *velx;
    //qDebug() << nsq << a << velx << NN <<pow(std::max(0.1,h),4.0/3.0) ;
    return velx;


        //double vkin = (slope  > 0? -1.0: 1.0) *sqrt(std::fabs(slope)) * pow(h,3.0/2.0)/NN;
        /*double fac = std::min(1.0,std::max(nsq,0.0));
        if(std::fabs(a) > 0)
        {
            return (1-fac)*velo + (fac) *(10.0 * (a/dt)/(kinfac*nsq));
        }else
        {
            return (1-fac)* velo;
        }*/
}

double TWorld::UF_Friction2(double vel,double dt,double velx,double vely, double NN, double h, double slope, bool solid, bool channel, double flowwidth)
{
    return UF_Friction(dt,velx,vely,NN,h,slope,channel,flowwidth);
    /*double manning = solid?UF_MANNINGCOEFFICIENT_SOLID:UF_MANNINGCOEFFICIENT_FLUID;
    double nsq = manning * (0.1+NN)*(0.1+NN)*UF_Gravity*sqrt(vel * vel)*pow(dt/UF2D_MinimumDT,1.0/3.0)/pow(std::max(0.0,h),4.0/3.0);

    if(channel)
    {
        if(flowwidth > 0)
        {
            nsq = nsq * (flowwidth + 2.0 * h)/(flowwidth);
            if(h > 1)
            {
                 nsq = nsq * h;
            }
        }
    }
    double kinfac = 0.5 +  0.5 * pow(std::max(0.0,std::min(1.0,(h/0.25))),2.0);
    //double vkin = (slope  > 0? -1.0: 1.0) *sqrt(std::fabs(slope)) * pow(h,3.0/2.0)/NN;
    return vel/(1.0+  kinfac * nsq); //vkin * kinfac + (1-kinfac) *velx/(1.0+nsq);*/
}

double TWorld::UF2D_MomentumBalanceFluid(bool x, double _f,double _s,double fu, double fv, double su, double sv, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbf, double SlopeX, double SlopeY,
                                 double dhfdx,double dhfdy, double dh2pbdx, double dh2pbdy, double dsfdx, double dsfdy, double ddsfdxx, double ddsfdyy, double ddsfdxy,
                                 double dfudx, double dfudy, double dfvdx, double dfvdy, double ddfudxx, double ddfudyy, double ddfvdxy, double ddfvdxx, double ddfvdyy, double ddfudxy,
                                 double dsudx, double dsudy, double dsvdx,double dsvdy)
{

    double h = (_f + _s)/(_dx*_dx);
    if(h < UF_VERY_SMALL)
    {
        return 0;
    }
    if(x) {

        return (-UF_Gravity * sin(SlopeX + dhfdx) - UF_Gravity *dhfdx
                -UF_Aspect *(
                    //(dh2pbdx)/h
                   + (pbf * (SlopeX+ dhfdx))
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

        return (-UF_Gravity * sin(SlopeY + dhfdy) - UF_Gravity *dhfdy
                -UF_Aspect * (
                    //(dh2pbdy)/h
                   + (pbf * (SlopeY+ dhfdy))
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
    double vel = sqrt(su*su + sv*sv);
    if(_f < UF_VERY_SMALL)
    {
        return 0;
    }
    if(x)
    {
        return (-UF_Gravity * sin(SlopeX + dhdx) +(vel > 0? su/vel : 0.0)*std::tan(ifa)*pbs +UF_Aspect*pbs*(SlopeX + dhdx)
                -UF_Aspect * gamma * pbf * ( dhdx +  SlopeX )
                 );
    }else
    {
        return (-UF_Gravity * sin(SlopeY + dhdy) +(vel > 0? sv/vel : 0.0)*std::tan(ifa)*pbs+UF_Aspect*pbs*(SlopeY + dhdy)
                +UF_Aspect * gamma * pbf * ( dhdy +  SlopeY )
               );
    }

}

double TWorld::UF1D_MomentumBalanceFluid(double _f,double _s,double fu, double su, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbf, double Slope,
                                 double dhfdx, double dh2pbdx, double dsfdx, double ddsfdxx, double dfudx, double ddfudxx, double dsudx)
{
    double h = (_f + _s)/(_dx*_dx);
    if(h < UF_VERY_SMALL)
    {
        return 0;
    }
    return
        (-UF_Gravity * sin(Slope + dhfdx) - //h * UF_Gravity *dhfdx   -
        UF_Aspect *(
             +(dh2pbdx)/h
             +(pbf * (Slope + dhfdx))
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
        (-UF_Gravity * sin(Slope + dhdx) - (std::fabs(su) > UF_VERY_SMALL? (su > 0? 1.0 : -1.0) :0.0)*std::tan(ifa)*pbs+UF_Aspect*pbs*(Slope + dhdx)
        -UF_Aspect * gamma * pbf * ( dhdx +  Slope)
         );


}

