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
 \file lisUnifiedFlowEntrainmentDeposition.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"




void TWorld::UF_FlowEntrainment(int thread)
{

    if(SwitchErosion && UF_SOLIDPHASE && SwitchEntrainment)
    {
        cTMap*_dem = UF2D_DEM;
        cTMap*_ldd = UF1D_LDD;

        FOR_ROW_COL_UF2DMT
        {
            UF2D_ST->Drc = UF_FlowEntrainmentST(r,c,false);
        }}}

        if(UF_1DACTIVE)
        {
            FOR_ROW_COL_UF1DMT
            {
                UF1D_ST->Drc = UF_FlowEntrainmentST(r,c,true);
            }}}
        }

        UF_LateralEntrainment(thread);


        FOR_ROW_COL_UF2DMT
        {

            UF_FlowEntrainment(UF2D_DT->Drc,r,c,false);


        }}}
        if(UF_1DACTIVE)
        {
            FOR_ROW_COL_UF1DMT
            {
                UF_FlowEntrainment(UF1D_DT->Drc,r,c,true);
            }}}
        }
    }

}

double TWorld::UF_FlowEntrainmentST(int r, int c, bool channel)
{
    //get all parameters for entrainment either from 1D or 2D maps (channel or not)
    double f = channel? UF1D_f->Drc : UF2D_f->Drc;
    double s = channel? UF1D_s->Drc : UF2D_s->Drc;
    double sf = channel? (UF1D_ssm->Drc + UF1D_blm->Drc) : (UF2D_ssm->Drc + UF2D_blm->Drc);

    double sconc = (f+s) >0 ? s/(f+s) : 0.0;//UF_5CellAverage(UF2D_tsf,r,c);

    double width = channel? UF1D_LDDw->Drc : _dx;
    double area = width * DX->Drc;
    double velocity = channel? std::fabs(UF1D_fu->Drc) : sqrt(UF2D_fu->Drc*UF2D_fu->Drc + UF2D_fv->Drc*UF2D_fv->Drc);
    double velocitys = channel? std::fabs(UF1D_su->Drc) : sqrt(UF2D_su->Drc*UF2D_su->Drc + UF2D_sv->Drc*UF2D_sv->Drc);
    //velocitys = (s+sf)>0? (s*velocitys + sf * velocity)/(s+sf):0.0;
    double visc = channel? UF1D_visc->Drc : UF2D_visc->Drc;
    double density = channel? (std::max(1000.0,UF1D_d->Drc) * (f+s) + UF_DENSITY_SUSPENDED * (sf/UF_DENSITY_SUSPENDED))/(s+f+ (sf/UF_DENSITY_SUSPENDED) )
                            : (std::max(1000.0,UF2D_d->Drc) * (f+s) + UF_DENSITY_SUSPENDED * (sf/UF_DENSITY_SUSPENDED))/(s+f+ (sf/UF_DENSITY_SUSPENDED) );
    double rocksize = channel? UF1D_rocksize->Drc : UF2D_rocksize->Drc;
    double ifa = channel? UF1D_ifa->Drc : UF2D_ifa->Drc;
    double bed_density = SoilRockDensity->Drc;
    double bed_ifa= SoilRockIFA->Drc;
    double slope = channel? std::fabs(UF1D_Slope->Drc) : std::max(std::fabs(UF2D_SlopeX->Drc),std::fabs(UF2D_SlopeY->Drc));
    double bed_cohesion = Cohesion->Drc * st_scCalibration;


    return UnifiedFlowActiveEntrainmentST(slope,f,s,area,velocity,velocitys,sconc,visc,density,ifa,rocksize,bed_density, bed_ifa, bed_cohesion, RootCohesion->Drc,N->Drc, r, c);

}
double TWorld::UnifiedFlowActiveEntrainmentST(double slope, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, double coh_bed, double veg_coh, double manning, int r, int c)
{
    double h = (_f+_s)/area;
    return UF_Gravity * h * d * (_fv*_fv + _sv*_sv)*0.5*(manning*manning/(pow(h,4.0/3.0)) + _sc * ifa) * area;

}

void TWorld::UF_LateralEntrainment(int thread)
{
    cTMap * temp = ThreadPool->UF_t1.at(thread);
    FOR_ROW_COL_UF2DMT
    {
         temp->Drc = UF2D_ST->Drc;
         UF2D_STL->Drc = 0;
         UF2D_STLA->Drc = 0;
         UF2D_STLH->Drc = 0;
    }}}

    FOR_ROW_COL_UF2DMT
    {
        double h = (UF2D_f->Drc + UF2D_s->Drc)/(_dx * DX->Drc);
        double fac = std::min(3.0,std::max(1.0,sqrt(h)));
        double width_factor = 0.01;

        if(!OUTORMV(r-1,c))
        {
            UF2D_STL->Drc += width_factor *fac * temp->data[r-1][c] * (std::min(h *_dx,UF_LateralEntrainmentArea(r-1,c,1,0)));
            UF2D_STLA->Drc += std::min(h *_dx,UF_LateralEntrainmentArea(r-1,c,1,0));
            UF2D_STLH->Drc = std::max(UF2D_STLH->Drc,(UF2D_f->data[r-1][c] + UF2D_s->data[r-1][c])/(_dx * DX->data[r-1][c]) );
        }
        if(!OUTORMV(r,c-1))
        {
            UF2D_STL->Drc += width_factor *fac * temp->data[r][c-1] * (std::min(h*_dx,UF_LateralEntrainmentArea(r,c-1,0,1)));
            UF2D_STLA->Drc += std::min(h*_dx,UF_LateralEntrainmentArea(r,c-1,0,1));
            UF2D_STLH->Drc = std::max(UF2D_STLH->Drc,(UF2D_f->data[r][c-1] + UF2D_s->data[r][c-1])/(_dx * DX->data[r][c-1]) );
        }
        if(!OUTORMV(r+1,c))
        {
            UF2D_STL->Drc += width_factor *fac * temp->data[r+1][c] * (std::min(h*_dx,UF_LateralEntrainmentArea(r+1,c,-1,0)));
            UF2D_STLA->Drc += std::min(h*_dx,UF_LateralEntrainmentArea(r+1,c,-1,0));
            UF2D_STLH->Drc = std::max(UF2D_STLH->Drc,(UF2D_f->data[r+1][c] + UF2D_s->data[r+1][c])/(_dx * DX->data[r+1][c]) );
        }
        if(!OUTORMV(r,c+1))
        {
            UF2D_STL->Drc += width_factor *fac * temp->data[r][c+1] * (std::min(h*_dx,UF_LateralEntrainmentArea(r,c+1,0,-1)));
            UF2D_STLA->Drc += std::min(h*_dx,UF_LateralEntrainmentArea(r,c+1,0,-1));
            UF2D_STLH->Drc = std::max(UF2D_STLH->Drc,(UF2D_f->data[r][c+1] + UF2D_s->data[r][c+1])/(_dx * DX->data[r][c+1]) );
        }
    }}}
}

double TWorld::UF_LateralEntrainmentArea(int r, int c, int dr, int dc)
{
    return _dx * std::max(0.0,UF2D_DEM->data[r+dr][c+dc]-UF2D_DEM->Drc);
}

void TWorld::UF_FlowEntrainment(double dt, int r, int c, bool channel)
{

    //get all parameters for entrainment either from 1D or 2D maps (channel or not)
    double f = channel? UF1D_f->Drc : UF2D_f->Drc;
    double s = channel? UF1D_s->Drc : UF2D_s->Drc;
    double sf = 0;//channel? (UF1D_ssm->Drc + UF1D_blm->Drc) : (UF2D_ssm->Drc + UF2D_blm->Drc);
    double sconc = f> 0? (s + (sf/UF_DENSITY_SUSPENDED))/(f+s + (sf/UF_DENSITY_SUSPENDED)) : 0.0;
    double width = channel? UF1D_LDDw->Drc : _dx;
    double area = width * DX->Drc;
    double velocity = channel? std::fabs(UF1D_fu->Drc) : sqrt(UF2D_fu->Drc*UF2D_fu->Drc + UF2D_fv->Drc*UF2D_fv->Drc);
    double velocitys = channel? std::fabs(UF1D_su->Drc) : sqrt(UF2D_su->Drc*UF2D_su->Drc + UF2D_sv->Drc*UF2D_sv->Drc);
    //velocitys = (s+sf)>0? (s*velocitys + sf * velocity)/(s+sf):0.0;
    double visc = channel? UF1D_visc->Drc : UF2D_visc->Drc;
    double density = channel? (std::max(1000.0,UF1D_d->Drc) * (f+s) + UF_DENSITY_SUSPENDED * (sf/UF_DENSITY_SUSPENDED))/(s+f+ (sf/UF_DENSITY_SUSPENDED) )
                            : (std::max(1000.0,UF2D_d->Drc) * (f+s) + UF_DENSITY_SUSPENDED * (sf/UF_DENSITY_SUSPENDED))/(s+f+ (sf/UF_DENSITY_SUSPENDED) );
    double rocksize = channel? UF1D_rocksize->Drc : UF2D_rocksize->Drc;
    double ifa = channel? UF1D_ifa->Drc : UF2D_ifa->Drc;
    double bed_density = SoilRockDensity->Drc;
    double bed_ifa= SoilRockIFA->Drc;

    double slope = channel? std::fabs(UF1D_Slope->Drc) : std::max(std::fabs(UF2D_SlopeX->Drc),std::fabs(UF2D_SlopeY->Drc));
    double slope_lat = 10.0 * slope;

    double availabledepth = channel? 0:UnifiedFlowEntrainmentAvailableDepth(r,c,UF2D_fu->Drc,UF2D_fv->Drc);
    double vegetationcover = Cover->Drc;
    double entrained_depth = TotalEntrainmentDep->Drc/ (_dx * _dx);
    double vegetationcohesion = (channel? std::max(0.0,((UF_ENTRAINMENTROOTDEPTH-entrained_depth)/UF_ENTRAINMENTROOTDEPTH)) : 1.0) * RootCohesion->Drc;
    double bed_cohesion = Cohesion->Drc * st_scCalibration;

    double shearstress = channel? UF1D_ST->Drc : UF2D_ST->Drc;

    double entrainment = UnifiedFlowActiveEntrainment(dt,shearstress, slope,f,s,area,velocity,velocitys,sconc,visc,density,ifa,rocksize,bed_density, bed_ifa, bed_cohesion, vegetationcohesion,N->Drc, r, c);
    double entrainment_lat = 0.0;//channel? 0.0:UnifiedFlowActiveEntrainmentLat(dt,UF2D_STL->Drc, slope_lat,UF2D_STLH->Drc,f,s,area,velocity,velocitys,sconc,visc,density,ifa,rocksize,bed_density, bed_ifa, bed_cohesion, vegetationcohesion,N->Drc, r, c);
    double entrainment_sf = 0.0;//channel? 0.0:UF_EntrainmentSideSlopeFailure(dt,r,c);

    double deposition = 0;

    if(SwitchDeposition)
    {
        deposition = UnifiedFlowActiveDeposition(dt,slope,f,s,area,velocity,velocitys,sconc,visc,density,ifa,rocksize,bed_density, bed_ifa,r,c);
    }

    entrainment = std::min(entrainment + entrainment_lat + entrainment_sf,0.05 * area * availabledepth);



    if(!channel)
    {
        UF2D_EntrainmentSF->Drc  =entrainment;

    }else
    {

    }

    if(entrainment > 0 && SoilRockMaterial->Drc > 0)
    {

        UF_RockTake(r,c,entrainment,channel);
    }else
    {

        if(deposition > 0)
        {
            UF_RockAdd(r,c,deposition,channel);
        }
    }

}

double TWorld::UnifiedFlowEntrainmentAvailableDepth(int r,int c, double vx, double vy)
{
    int dr = vy > 0? 1.0:-1.0;
    int dc = vx > 0? 1.0:-1.0;

    double dem = UF2D_DEM->data[r][c];

    double depth1 = 0;
    double depth2 = 0;
    //double depth3 = 0;
    //double depth4 = 0;

    if(!OUTORMV(r + dr,c))
    {
        if(UF2D_DEM->data[r+dr][c] < dem)
        {
            depth1 = std::max(0.0,dem - UF2D_DEM->data[r+dr][c]);
        }
    }

    if(!OUTORMV(r,c+dc))
    {
        if(UF2D_DEM->data[r][c+dc] < dem)
        {
            depth2 = std::max(0.0,dem - UF2D_DEM->data[r][c+dc]);
        }
    }


    return std::max(depth1,depth2);
}

double TWorld::UnifiedFlowActiveEntrainmentLat(double dt,double st, double slope, double h, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, double coh_bed, double veg_coh, double manning, int r, int c)
{
    double entrainment = 0;

    double UF_SOILROCKPOROSITY = 0.65;

    //Hungr
    /*double dgamma = d/d_bed;
    double sf = _f > 0? (_s/_f) :1.0;
    entrainment =(sf > UF_MAXSOLIDCONCENTRATION)?0.0: std::max(0.0,std::min(h * (UF_MAXSOLIDCONCENTRATION - sf),dt * area * ( UF_ENTRAINMENTCONSTANT * (0.5 *_fv +sf * _sv))));*/

    //RAMMS
    /*double dgamma = d/d_bed;
    double sf = _f > 0? (_s/_f) :1.0;
    entrainment =(sf > UF_MAXSOLIDCONCENTRATION)?0.0: std::max(0.0,std::min( h * (UF_MAXSOLIDCONCENTRATION - sf),dt * area * h * ( UF_ENTRAINMENTCONSTANT * (0.5 *_fv +sf * _sv))));*/

    ////Egashira

    /*double densdiff = (d - 1000.0);
    double tanifa = tan(ifa_bed);
    double tanalpha = tan(slope);
    double tanalphae = tan( atan(densdiff/(densdiff + 1000.0))*tanifa);
    entrainment = UF_ENTRAINMENTCONSTANT * area  * _sv * UF_SOILROCKPOROSITY * (tanalphae - tanalpha);*/
    //Egashira can be negative, usefull to include?

    ////Pudasaini

    /*double entrainment_solid = UF_ENTRAINMENTCONSTANT *((_f + _s) > 0? (_s/(_f + _s)): 0.0) * 0.003 * area * std::sqrt(h * cos(slope) * UF_Gravity);
    double entrainment_fluid = UF_ENTRAINMENTCONSTANT *((_f + _s) > 0? (_f/(_f + _s)): 0.0) * 0.002 * area * _fv;

    entrainment = entrainment_solid + entrainment_fluid;*/

    ////Takahashi
    //first get maximum solids concentration that still allows entrainment

    double MaxCSF = std::max(0.0,std::min(0.8,(UF_ENTRAINMENTCCONSTANT)*slope > tan(ifa_bed)? 1.0:(1000.0 * (UF_ENTRAINMENTCCONSTANT)*slope)/((d - 1000)*(tan(ifa_bed)-(UF_ENTRAINMENTCCONSTANT)*slope))));

    //shear stress
    double gamma = std::min(1.0,d > UF_VERY_SMALL? 1000.0/d : 1.0);
    double pbs = (1-gamma)*(-UF_Gravity * h);
    double dc = UF_DragCoefficient(_f/(_f+_s),_sc,gamma ,visc,rocksize,d);

    double t = UF_Gravity * h * d * (_fv*_fv + _sv*_sv)*0.5*(manning*manning/(pow(h,4.0/3.0)) + _sc * ifa);

    double Coeff_Susp = 0.5;

    double coh = coh_bed + veg_coh;

    //critical shear stress
    double tc = (coh + (1-Coeff_Susp) *_sc * (d - 1000.0) * UF_Gravity * h * (cos(slope)*cos(slope) * tan(ifa_bed)));

    //get the actual scouring rate
    //double scourat = (UF_ENTRAINMENTCONSTANT * h * std::sqrt( _fv*_fv + _sv*_sv)*(MaxCSF - (_s/_f+_s)))/((UF_SOILROCKPOROSITY - MaxCSF)*rocksize);
    double scourat = std::max(0.0,UF_ENTRAINMENTCONSTANT * (t- UF_ENTRAINMENTTHRESHOLDCONSTANT *tc));
    //get entrainment in cubic meters
    entrainment = std::max(0.0,std::min(0.5 * (MaxCSF - _sc)*area * h,scourat *area*dt));



    if(area < UF_VERY_SMALL)
    {
        return 0;
    }
    if(!(h > UF_MINIMUMENTRAINMENTHEIGHT))
    {
        return 0;
    }

    //returns volume of entrainment
    return std::max(0.0,entrainment);
}

double TWorld::UnifiedFlowActiveEntrainment(double dt,double st, double slope, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, double coh_bed, double veg_coh, double manning, int r, int c)
{

    double entrainment = 0;
    double h = (_f + _s)/(area);


    double UF_SOILROCKPOROSITY = 0.65;

    //Hungr
    /*double dgamma = d/d_bed;
    double sf = _f > 0? (_s/_f) :1.0;
    entrainment =(sf > UF_MAXSOLIDCONCENTRATION)?0.0: std::max(0.0,std::min(h * (UF_MAXSOLIDCONCENTRATION - sf),dt * area * ( UF_ENTRAINMENTCONSTANT * (0.5 *_fv +sf * _sv))));*/

    //RAMMS
    /*double dgamma = d/d_bed;
    double sf = _f > 0? (_s/_f) :1.0;
    entrainment =(sf > UF_MAXSOLIDCONCENTRATION)?0.0: std::max(0.0,std::min( h * (UF_MAXSOLIDCONCENTRATION - sf),dt * area * h * ( UF_ENTRAINMENTCONSTANT * (0.5 *_fv +sf * _sv))));*/

    ////Egashira

    /*double densdiff = (d - 1000.0);
    double tanifa = tan(ifa_bed);
    double tanalpha = tan(slope);
    double tanalphae = tan( atan(densdiff/(densdiff + 1000.0))*tanifa);
    entrainment = UF_ENTRAINMENTCONSTANT * area  * _sv * UF_SOILROCKPOROSITY * (tanalphae - tanalpha);*/
    //Egashira can be negative, usefull to include?

    ////Pudasaini

    /*double entrainment_solid = UF_ENTRAINMENTCONSTANT *((_f + _s) > 0? (_s/(_f + _s)): 0.0) * 0.003 * area * std::sqrt(h * cos(slope) * UF_Gravity);
    double entrainment_fluid = UF_ENTRAINMENTCONSTANT *((_f + _s) > 0? (_f/(_f + _s)): 0.0) * 0.002 * area * _fv;

    entrainment = entrainment_solid + entrainment_fluid;*/


    ////Takahashi
    //first get maximum solids concentration that still allows entrainment

    double MaxCSF = std::max(0.1,std::min(0.8,(UF_ENTRAINMENTCCONSTANT)*slope > tan(ifa_bed)? 1.0:(1000.0 * (UF_ENTRAINMENTCCONSTANT)*slope)/((d - 1000)*(tan(ifa_bed)-(UF_ENTRAINMENTCCONSTANT)*slope))));

    //shear stress
    double gamma = std::min(1.0,d > UF_VERY_SMALL? 1000.0/d : 1.0);
    double pbs = (1-gamma)*(-UF_Gravity * h);
    double dc = UF_DragCoefficient(_f/(_f+_s),_sc,gamma ,visc,rocksize,d);

    double t = UF_Gravity * h * d * (_fv*_fv + _sv*_sv)*0.5*(manning*manning/(pow(h,4.0/3.0)) + _sc * ifa);

    double Coeff_Susp = 0.5;

    double coh = coh_bed + veg_coh;

    //critical shear stress
    double tc = (coh + (1-Coeff_Susp) *_sc * (d - 1000.0) * UF_Gravity * h * (cos(slope)*cos(slope) * tan(ifa_bed)));


    //get the actual scouring rate
    //double scourat = (UF_ENTRAINMENTCONSTANT * h * std::sqrt( _fv*_fv + _sv*_sv)*(MaxCSF - (_s/_f+_s)))/((UF_SOILROCKPOROSITY - MaxCSF)*rocksize);
    double scourat = std::max(0.0,UF_ENTRAINMENTCONSTANT * (t- UF_ENTRAINMENTTHRESHOLDCONSTANT *tc));
    //get entrainment in cubic meters
    entrainment = std::max(0.0,std::min(0.5 * (MaxCSF - _sc)*area * h,scourat *area*dt));



    Entrainmentshearstressc->Drc = MaxCSF;
    Entrainmentshearstress->Drc = _sc;

    if(area < UF_VERY_SMALL)
    {
        return 0;
    }
    if(!(h > UF_MINIMUMENTRAINMENTHEIGHT))
    {
        return 0;
    }

    //returns volume of entrainment
    return std::max(0.0,entrainment);
}

double TWorld::UnifiedFlowActiveDeposition(double dt,double slope, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, int r, int c)
{
    double h = (_f + _s)/(area);
    _sc = _s/(_f+_s);

    double d_s = (_sc > UF_VERY_SMALL)? (d - 1000 * (1-_sc))/_sc : 2000.0;
    double vc = _sc > UF_VERY_SMALL? ( 0.15 + (1/5.0) *sqrt(UF_Gravity ) * pow(h,0.5)) : 0.0; //sin(atan((_sc *(d_s - 1000.0) * tan(ifa)/(_sc*(d_s - 1000.0) + 1000.0))))* d /(0.02 * d_s)) *(pow(0.7/_sc,1/3.0)-1.0)

    double MaxCSF = std::max(0.25 * std::min(1.0,_sv * h),std::min(0.8,UF_ENTRAINMENTCCONSTANT*slope > tan(ifa_bed)? 1.0:(1000.0 *UF_ENTRAINMENTCCONSTANT* slope)/((d - 1000)*(tan(ifa_bed)-UF_ENTRAINMENTCCONSTANT*slope))));

    double deporat = UF_DEPOSITIONCONSTANT *std::max(_sc - (MaxCSF)*_f,std::max(0.0 ,(1.0-_fv/(UF_DEPOSITIONTHRESHOLDCONSTANT*vc)))*std::max(0.0,((_sc-MaxCSF)/0.7)*_fv));


    //returns volume of deposition
    return std::min(0.5*_s, deporat * area * dt );
}

double TWorld::UF_RockTake(int r, int c, double entrainment, bool channel)
{


    if(channel)
    {
        //convert to kg, and limimt to present material

        entrainment = std::min(RSoilRockMaterial->Drc, entrainment);
        double theta = SoilRockMaterial->Drc >0?RSoilRockWater->Drc/RSoilRockMaterial->Drc:0.0;
        RSoilRockMaterial->Drc -= entrainment;
        RSoilRockWater->Drc -= entrainment *theta;

        ChannelEntrainmentDet->Drc = entrainment;
        LDDChange->Drc -= entrainment/(_dx*UF1D_LDDw->Drc);


        //back to volume
        UF1D_d->Drc = (UF1D_s->Drc + entrainment) > UF_VERY_SMALL? (UF1D_s->Drc * UF1D_d->Drc + entrainment * RSoilRockDensity->Drc)/(UF1D_s->Drc + entrainment) : UF1D_d->Drc;
        UF1D_rocksize->Drc = (UF1D_s->Drc + entrainment) > UF_VERY_SMALL? (UF1D_s->Drc * UF1D_rocksize->Drc + entrainment * RSoilRockSize->Drc)/(UF1D_s->Drc + entrainment) : UF1D_rocksize->Drc;
        UF1D_ifa->Drc = (UF1D_s->Drc + entrainment) > UF_VERY_SMALL? (UF1D_s->Drc * UF1D_ifa->Drc + entrainment * RSoilRockIFA->Drc)/(UF1D_s->Drc + entrainment) : UF1D_ifa->Drc;
        UF1D_s->Drc += entrainment;
        UF1D_f->Drc += entrainment *theta;

        return entrainment;
    }else
    {
        //convert to kg, and limimt to present material

        entrainment = std::min(SoilRockMaterial->Drc, entrainment);
        double theta = SoilRockMaterial->Drc >0?SoilRockWater->Drc/SoilRockMaterial->Drc:0.0;
        SoilRockMaterial->Drc -= entrainment;
        SoilRockWater->Drc -= entrainment *theta;
        EntrainmentDet->Drc += entrainment;
        DEMChange->Drc -= entrainment/(_dx*_dx);


        if(entrainment > 0)
        {
            //back to volume
            UF2D_d->Drc = (UF2D_s->Drc + entrainment) > UF_VERY_SMALL? (UF2D_s->Drc * UF2D_d->Drc + entrainment * SoilRockDensity->Drc)/(UF2D_s->Drc + entrainment) : UF2D_d->Drc;
            UF2D_rocksize->Drc = (UF2D_s->Drc + entrainment) > UF_VERY_SMALL? (UF2D_s->Drc * UF2D_rocksize->Drc + entrainment * SoilRockSize->Drc)/(UF2D_s->Drc + entrainment) : UF2D_rocksize->Drc;
            UF2D_ifa->Drc = (UF2D_s->Drc + entrainment) > UF_VERY_SMALL? (UF2D_s->Drc * UF2D_ifa->Drc + entrainment * SoilRockIFA->Drc)/(UF2D_s->Drc + entrainment) : UF2D_ifa->Drc;
            UF2D_s->Drc += entrainment;
            UF2D_f->Drc += entrainment *theta;
        }
        return entrainment;
    }


}

double TWorld::UF_RockAdd(int r, int c, double deposition, bool channel)
{
    deposition = std::fabs(deposition);
    if(channel)
    {
        //convert to kg, and limimt to present material
        deposition = std::min(UF1D_s->Drc, deposition);
        double depositionw = std::min(UF1D_f->Drc, deposition * UF_FLOWCOMPACTION_DEPOSITIONPOROSITY);
        EntrainmentDep->Drc += deposition;
        UF1D_s->Drc -= deposition;
        UF1D_f->Drc -= depositionw;
        LDDChange->Drc += deposition/(_dx*UF1D_LDDw->Drc);
        RSoilRockDensity->Drc = (RSoilRockMaterial->Drc + deposition) > UF_VERY_SMALL? (deposition * UF1D_d->Drc + RSoilRockMaterial->Drc * RSoilRockDensity->Drc)/(RSoilRockMaterial->Drc + deposition) : RSoilRockDensity->Drc;
        RSoilRockSize->Drc = (RSoilRockMaterial->Drc + deposition) > UF_VERY_SMALL? (deposition * UF1D_rocksize->Drc + RSoilRockMaterial->Drc * RSoilRockSize->Drc)/(RSoilRockMaterial->Drc + deposition) : RSoilRockSize->Drc;
        RSoilRockIFA->Drc = (RSoilRockMaterial->Drc + deposition) > UF_VERY_SMALL? (deposition * UF1D_ifa->Drc + RSoilRockMaterial->Drc * RSoilRockIFA->Drc)/(RSoilRockMaterial->Drc + deposition) : RSoilRockIFA->Drc;
        RSoilRockMaterial->Drc += deposition;
    }else
    {
        //convert to kg, and limimt to present material
        deposition = std::min(UF2D_s->Drc, deposition);
        double depositionw = std::min(UF2D_f->Drc, deposition * UF_FLOWCOMPACTION_DEPOSITIONPOROSITY);
        EntrainmentDep->Drc += deposition;
        UF2D_s->Drc -= deposition;
        UF2D_f->Drc -= depositionw;
        DEMChange->Drc += deposition/(_dx*_dx);

        SoilRockDensity->Drc = (SoilRockMaterial->Drc + deposition) > UF_VERY_SMALL? (deposition * UF2D_d->Drc + SoilRockMaterial->Drc * SoilRockDensity->Drc)/(SoilRockMaterial->Drc + deposition) : SoilRockDensity->Drc;
        SoilRockSize->Drc = (SoilRockMaterial->Drc + deposition) > UF_VERY_SMALL? (deposition * UF2D_rocksize->Drc + SoilRockMaterial->Drc * SoilRockSize->Drc)/(SoilRockMaterial->Drc + deposition) : SoilRockSize->Drc;
        SoilRockIFA->Drc = (SoilRockMaterial->Drc + deposition) > UF_VERY_SMALL? (deposition * UF2D_ifa->Drc + SoilRockMaterial->Drc * SoilRockIFA->Drc)/(SoilRockMaterial->Drc + deposition) : SoilRockIFA->Drc;
        SoilRockMaterial->Drc += deposition;
    }



}

void TWorld::UF_FlowCompaction(int thread)
{
    if(SwitchCompaction)
    {

        if(SwitchErosion && UF_SOLIDPHASE )
        {
            cTMap*_dem = UF2D_DEM;
            cTMap*_ldd = UF1D_LDD;

            FOR_ROW_COL_UF2DMT
            {
                UF_FlowCompactionThreshold(UF2D_DT->Drc,r,c,false);

                UF_FlowCompaction(UF2D_DT->Drc,r,c,false);
            }}}

            if(UF_1DACTIVE)
            {
                FOR_ROW_COL_UF1DMT
                {
                    UF_FlowCompaction(UF1D_DT->Drc,r,c,true);
                }}}
            }
        }
    }

}

void TWorld::UF_FlowCompactionThreshold(double dt, int r, int c, bool channel)
{


    if(channel)
    {
        return;
    }

    double f = channel? UF1D_f->Drc : UF2D_f->Drc;
    double s = channel? UF1D_s->Drc : UF2D_s->Drc;
    double width = channel? UF1D_LDDw->Drc : _dx;
    double area = width * DX->Drc;

    double h = (s)/(area);
    double ff = f/(f+s);
    double sf = s/(f+s);
    double gamma = channel? UF1D_d->Drc : UF2D_d->Drc;
    gamma = gamma > UF_VERY_SMALL? 1000.0/gamma: 0.5;
    double dc = UF_DragCoefficient(ff,sf,gamma,channel?UF1D_visc->Drc:UF2D_visc->Drc, channel?UF1D_rocksize->Drc:UF2D_rocksize->Drc, channel?UF1D_d->Drc:UF2D_d->Drc);
    double pbf = -UF_Gravity;

    double pbs = (1-gamma)*pbf;
    double ifa = UF2D_ifa->Drc < 0.01? 0.3:UF2D_ifa->Drc;

    double xslope = UF2D_Derivative(UF2D_DEM,UF2D_DEM,r,c,UF_DIRECTION_X);
    double yslope = UF2D_Derivative(UF2D_DEM,UF2D_DEM,r,c,UF_DIRECTION_Y);

    double xh1 = !OUTORMV(r,c-1)? UF2D_s->data[r][c-1] / (_dx*DX->data[r][c-1]) : h;
    double xh2 = !OUTORMV(r,c+1)? UF2D_s->data[r][c+1] / (_dx*DX->data[r][c+1]) : h;
    double yh1 = !OUTORMV(r-1,c)? UF2D_s->data[r-1][c] / (_dx*DX->data[r-1][c]) : h;
    double yh2 = !OUTORMV(r+1,c)? UF2D_s->data[r+1][c] / (_dx*DX->data[r+1][c]) : h;

    double dhdx = ((h-xh1)+( xh2-h))/(_dx*2.0);
    double dhdy = ((h-yh1)+( yh2-h))/(_dx*2.0);

    double su = UF2D_su->Drc;
    double sv = UF2D_sv->Drc;

    double sux1 = UF2D_su->Drc;
    double svx1 = UF2D_sv->Drc;

    double sux2 = UF2D_su->Drc;
    double svx2 = UF2D_sv->Drc;

    double suy1 = UF2D_su->Drc;
    double svy1 = UF2D_sv->Drc;

    double suy2 = UF2D_su->Drc;
    double svy2 = UF2D_sv->Drc;

    double vel = sqrt(su*su + sv*sv);

    //x friction
    double friction_x = std::tan(ifa)* UF_Gravity;

    //y friction
    double friction_y = std::tan(ifa)* UF_Gravity;

    //x accaleration
    double acc_x = (-UF_Gravity * sin(xslope + dhdx) +UF_Aspect*pbs*(xslope + dhdx)
                -UF_Aspect * gamma * pbf * ( dhdx ));


    //y accaleration
    double acc_y = (-UF_Gravity * sin(yslope + dhdy)+UF_Aspect*pbs*(yslope + dhdy)
                +UF_Aspect * gamma * pbf * ( dhdy +yslope ));

    double frictioncontribution = std::min(0.05,std::max(0.0,(std::fabs(friction_x) + std::fabs(friction_y)) - (std::fabs(acc_x) + std::fabs(acc_y)))/(std::fabs(acc_x) + std::fabs(acc_y)));

    if(std::fabs(UF2D_sqx->Drc) + std::fabs(UF2D_sqy->Drc) < 0.5 && s > UF_VERY_SMALL && UF2D_EntrainmentSF->Drc == 0)
    {
        UF2D_Compaction->Drc = frictioncontribution;
    }else
    {
        UF2D_Compaction->Drc = 0;
    }

}


void TWorld::UF_FlowCompaction(double dt, int r, int c, bool channel)
{

    double s = channel? UF1D_s->Drc : UF2D_s->Drc;
    double width = channel? UF1D_LDDw->Drc : _dx;
    double area = width * DX->Drc;

    if(channel)
    {
        return;
    }

    double compaction = UF2D_Compaction->Drc;
    if(!OUTORMV(r+1,c))
    {
        compaction = std::max(compaction,0.5 * UF2D_Compaction->data[r+1][c]);
    }
    if(!OUTORMV(r-1,c))
    {
        compaction = std::max(compaction,0.5 * UF2D_Compaction->data[r-1][c]);
    }
    if(!OUTORMV(r,c+1))
    {
        compaction = std::max(compaction,0.5 * UF2D_Compaction->data[r][c+1]);
    }
    if(!OUTORMV(r,c-1))
    {
        compaction = std::max(compaction,0.5 * UF2D_Compaction->data[r][c-1]);
    }

    if(compaction > 0)
    {
        DepositionT->Drc += dt;

        if(DepositionT->Drc > 5.0)
        {
            {
                UF_RockAdd(r,c,s * std::min(0.1,0.01*dt),channel);
            }
        }
    }else
    {
        DepositionT->Drc = 0;
    }
}

double TWorld::UnifiedFlowDepositionAvailableDepth(int r, int c)
{
    return std::max(0.0,SolveStableDepthAt(r,c) - GetTotalSoilDepth(r,c));
}

double TWorld::UF_EntrainmentSideSlopeFailure(double dt, int r, int c)
{
    if(SwitchSlopeStability)
    {
        double volume = 0;
        double h = (UF2D_f->Drc + UF2D_s->Drc) /(_dx*DX->Drc);
        double sat = SoilRockWater->Drc/SoilRockMaterial->Drc;

        int dx[4] = {0,1,0,-1};
        int dy[4] = {1,0,-1,0};

        for(int i = 0; i < 4; i++)
        {
            int r2 = r + dy[i];
            int c2 = c + dx[i];

            if(OUTORMV(r2,c2))
            {
                continue;
            }

            double el_this = UF2D_DEM->Drc;
            double el_side = UF2D_DEM->data[r2][c2];

            double d_this = SoilRockMaterial->Drc / (_dx * DX->Drc);
            double d_side = SoilRockMaterial->data[r2][c2] / (_dx * DX->data[r2][c2]);

            double dd_this = DEMChange->Drc;
            double dd_side = DEMChange->data[r2][c2];

            if(dd_this < dd_side)
            {
                continue;
            }

            double dx = _dx/5.0;
            double slope =std::max(0.0,(dd_this-dd_side)/dx);

            double sf = 0.0;
            if(d_this < 0.1 || slope < 0.01)
            {
                sf = 1000.0;
            }else
            {
                sf = CalculateSafetyFactor(slope,d_this,h,10.0,SoilRockIFA->Drc,sat * d_this,0.0,SoilRockDensity->Drc,RootCohesion->Drc,0.0,1.0);
            }

            if(sf < 0.9)
            {
                double sd = SolveStableDepth(slope,d_this,h,10.0,SoilRockIFA->Drc,sat * d_this,0.0,SoilRockDensity->Drc,RootCohesion->Drc,0.0,1.0,1.0,1.0);
                volume += std::max(0.0,d_this - sd) * _dx * DX->Drc;
            }


        }

        return volume;
    }else
    {
        return 0.0;
    }

}
