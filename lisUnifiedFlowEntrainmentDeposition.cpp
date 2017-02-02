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

void TWorld::UF_FlowEntrainment(double dt, int r, int c, bool channel)
{

    //get all parameters for entrainment either from 1D or 2D maps (channel or not)
    double f = channel? UF1D_f->Drc : UF2D_f->Drc;
    double s = channel? UF1D_s->Drc : UF2D_s->Drc;
    double sf = channel? (UF1D_ssm->Drc + UF1D_blm->Drc) : (UF2D_ssm->Drc + UF2D_blm->Drc);
    double sconc = f> 0? (s + (sf/UF_DENSITY_SUSPENDED))/(f+s + (sf/UF_DENSITY_SUSPENDED)) : 0.0;
    double width = channel? UF1D_LDDw->Drc : _dx;
    double area = width * DX->Drc;
    double velocity = channel? std::fabs(UF1D_fu->Drc) : sqrt(UF2D_fu->Drc*UF2D_fu->Drc + UF2D_fv->Drc*UF2D_fv->Drc);
    double velocitys = channel? std::fabs(UF1D_su->Drc) : sqrt(UF2D_su->Drc*UF2D_su->Drc + UF2D_sv->Drc*UF2D_sv->Drc);
    velocitys = (s+sf)>0? (s*velocitys + sf * velocity)/(s+sf):0.0;
    double visc = channel? UF1D_visc->Drc : UF2D_visc->Drc;
    double density = channel? (std::max(1000.0,UF1D_d->Drc) * (f+s) + UF_DENSITY_SUSPENDED * (sf/UF_DENSITY_SUSPENDED))/(s+f+ (sf/UF_DENSITY_SUSPENDED) )
                            : (std::max(1000.0,UF2D_d->Drc) * (f+s) + UF_DENSITY_SUSPENDED * (sf/UF_DENSITY_SUSPENDED))/(s+f+ (sf/UF_DENSITY_SUSPENDED) );
    double rocksize = channel? UF1D_rocksize->Drc : UF2D_rocksize->Drc;
    double ifa = channel? UF1D_ifa->Drc : UF2D_ifa->Drc;
    double bed_density = SoilRockDensity->Drc;
    double bed_ifa= SoilRockIFA->Drc;
    double slope = channel? std::fabs(UF1D_Slope->Drc) : std::max(std::fabs(UF2D_SlopeX->Drc),std::fabs(UF2D_SlopeY->Drc));
    double availabledepth = channel? 0:UnifiedFlowEntrainmentAvailableDepth(r,c,UF2D_su->Drc,UF2D_sv->Drc);
    double vegetationcover = Cover->Drc;
    double vegetationcohesion = RootCohesion->Drc;
    double bed_cohesion = Cohesion->Drc * st_scCalibration;

    double entrainment = UnifiedFlowActiveEntrainment(dt,slope,f,s,area,velocity,velocitys,sconc,visc,density,ifa,rocksize,bed_density, bed_ifa, bed_cohesion, RootCohesion->Drc,N->Drc, r, c);
    double deposition = UnifiedFlowActiveDeposition(dt,slope,f,s,area,velocity,velocitys,sconc,visc,density,ifa,rocksize,bed_density, bed_ifa,r,c);

    entrainment = std::min(entrainment,area * availabledepth);

    if(entrainment > 0)
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

    double dem = DEM->data[r][c];

    double depth1 = 0;
    double depth2 = 0;
    if(!OUTORMV(r + dr,c))
    {
        if(DEM->data[r+dr][c] < dem)
        {
            depth1 = std::max(0.0,dem - DEM->data[r+dr][c]);
        }
    }

    if(!OUTORMV(r,c+dc))
    {
        if(DEM->data[r][c+dc] < dem)
        {
            depth2 = std::max(0.0,dem - DEM->data[r][c+dc]);
        }
    }

    if(OUTORMV(r + dr,c) && OUTORMV(r + dr,c))
    {
        if(!OUTORMV(r - dr,c))
        {
            if(DEM->data[r-dr][c] < dem)
            {
                depth1 = std::max(0.0,dem - DEM->data[r-dr][c]);
            }
        }

        if(!OUTORMV(r,c-dc))
        {
            if(DEM->data[r][c-dc] < dem)
            {
                depth2 = std::max(0.0,dem - DEM->data[r][c-dc]);
            }
        }
        if(OUTORMV(r - dr,c) && OUTORMV(r - dr,c))
        {
            return 0;
        }
    }

    return std::max(depth1,depth2);
}

double TWorld::UnifiedFlowActiveEntrainment(double dt,double slope, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, double coh_bed, double veg_coh, double manning, int r, int c)
{

    double entrainment = 0;
    double h = (_f + _s)/(area);


    double UF_SOILROCKPOROSITY = 0.65;

    //Hungr
    //double dgamma = d/d_bed;
    //double sf = _f > 0? (_s/_f) :1.0;
    //entrainment =(sf > UF_MAXSOLIDCONCENTRATION)?0.0: std::max(0.0,std::min(h * (UF_MAXSOLIDCONCENTRATION - sf),dt * area * ( UF_ENTRAINMENTCONSTANT * (0.5 *_fv +sf * _sv))));

    //Egashira

    /*double densdiff = (d - 1000.0);
    double tanifa = tan(ifa_bed);
    double tanalpha = tan(slope);
    double tanalphae = tan( atan(densdiff/(densdiff + 1000.0))*tanifa);
    entrainment = UF_ENTRAINMENTCONSTANT * area  * _sv * UF_SOILROCKPOROSITY * (tanalphae - tanalpha);*/

    //Egashira can be negative, usefull to include?
    //returns volume of entrainment


    ////Takahashi
    //first get maximum solids concentration that still allows entrainment

    double MaxCSF = std::max(0.0,std::min(0.8,slope > tan(ifa_bed)? 1.0:(1000.0 * slope)/((d - 1000)*(tan(ifa_bed)-slope))));

    //shear stress
    double gamma = std::min(1.0,d > UF_VERY_SMALL? 1000.0/d : 1.0);
    double pbs = (1-gamma)*(-UF_Gravity * h);
    double dc = UF_DragCoefficient(_f/(_f+_s),_sc,gamma ,visc,rocksize,d);


    double t = UF_Gravity * h * d * (_fv*_fv + _sv*_sv)*0.5*(manning*manning/(pow(h,4.0/3.0)) + _sc * ifa);

    double Coeff_Susp = 0.5;

    double coh = coh_bed + veg_coh;

    //critical shear stress
    double tc = coh + (1-Coeff_Susp) *_sc * (d - 1000.0) * UF_Gravity * h * (cos(slope)*cos(slope) * tan(ifa_bed));

    //get the actual scouring rate
    //double scourat = (UF_ENTRAINMENTCONSTANT * h * std::sqrt( _fv*_fv + _sv*_sv)*(MaxCSF - (_s/_f+_s)))/((UF_SOILROCKPOROSITY - MaxCSF)*rocksize);
    double scourat = std::max(0.0,UF_ENTRAINMENTCONSTANT * (t-tc));
    //get entrainment in cubic meters
    entrainment = std::max(0.0,std::min(0.5 * (MaxCSF - _sc)*area * h,scourat *area*dt));

    EntrainmentTC->Drc = MaxCSF;


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

    double MaxCSF = std::max(0.25 * std::min(1.0,_fv * h),std::min(0.8,slope > tan(ifa_bed)? 1.0:(1000.0 * slope)/((d - 1000)*(tan(ifa_bed)-slope))));

    double deporat = UF_DEPOSITIONCONSTANT *std::max(0.0,(1.0-_fv/(UF_DEPOSITIONTHRESHOLDCONSTANT*vc)))*std::max(0.0,((_sc-MaxCSF)/0.7)*_fv);



    Entrainmentshearstressc->Drc = vc;
    Entrainmentshearstress->Drc = (1.0-_fv/(UF_DEPOSITIONTHRESHOLDCONSTANT*vc));

    //returns volume of deposition
    return std::min(0.5*_s, deporat * area * dt );
}

double TWorld::UF_RockTake(int r, int c, double entrainment, bool channel)
{


    if(channel)
    {
        //convert to kg, and limimt to present material

        entrainment = std::min(RSoilRockMaterial->Drc, entrainment);
        double theta = SoilRockMaterial->Drc >0?SoilRockWater->Drc/SoilRockMaterial->Drc:0.0;
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
        EntrainmentDet->Drc = entrainment;
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
        EntrainmentDep->Drc = deposition;
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
        EntrainmentDep->Drc = deposition;
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

    /*if(SwitchErosion && UF_SOLIDPHASE)
    {
        cTMap*_dem = UF2D_DEM;
        cTMap*_ldd = UF1D_LDD;

        FOR_ROW_COL_UF2DMT
        {
            UF_FlowCompaction(UF2D_DT->Drc,r,c,false);
        }}}

        if(UF_1DACTIVE)
        {
            FOR_ROW_COL_UF1DMT
            {
                UF_FlowCompaction(UF1D_DT->Drc,r,c,true);
            }}}
        }
    }*/

}

void TWorld::UF_FlowCompaction(double dt, int r, int c, bool channel)
{

    double availabledepth = UnifiedFlowDepositionAvailableDepth(r,c);

    double f = channel? UF1D_f->Drc : UF2D_f->Drc;
    double s = channel? UF1D_s->Drc : UF2D_s->Drc;
    double sf = channel? (UF1D_ssm->Drc + UF1D_blm->Drc) : (UF2D_ssm->Drc + UF2D_blm->Drc);
    double sconc = f> 0? (s + (sf/UF_DENSITY_SUSPENDED))/(f+s + (sf/UF_DENSITY_SUSPENDED)) : 0.0;
    double width = channel? UF1D_LDDw->Drc : _dx;
    double area = width * DX->Drc;
    double velocity = channel? std::fabs(UF1D_fu->Drc) : sqrt(UF2D_fu->Drc*UF2D_fu->Drc + UF2D_fv->Drc*UF2D_fv->Drc);
    double velocitys = channel? std::fabs(UF1D_su->Drc) : sqrt(UF2D_su->Drc*UF2D_su->Drc + UF2D_sv->Drc*UF2D_sv->Drc);
    velocitys = (s+sf)>0? (s*velocitys + sf * velocity)/(s+sf):0.0;

    double sc = s/(f +s);
    if(sc > UF_FLOWCOMPACTION_MAXCONCENTRATION)
    {
        double desired_depth = 0;

        //solve (_s-d)/((_f-p*d)+(_s-d)) = UF_FLOWCOMPACTION_MAXCONCENTRATION;
        double p = 0.3;// (water taken with compacted soil)
        double desired_volume = std::max(0.0,(UF_FLOWCOMPACTION_MAXCONCENTRATION *(f+s) -s)/(-1 + (1+p)*UF_FLOWCOMPACTION_MAXCONCENTRATION));
        desired_depth = desired_volume/area;

        double v_c = UF_FLOWCOMPACTION_CRITICALVELOCITY;
        double compaction = std::max(0.0,(1-velocitys/v_c) *std::min(desired_depth,availabledepth) * dt/10.0);

        if(compaction > 0)
        {
            UF_RockAdd(r,c,compaction,channel);
        }
    }



}

double TWorld::UnifiedFlowDepositionAvailableDepth(int r, int c)
{
    return std::max(0.0,SolveStableDepthAt(r,c) - GetTotalSoilDepth(r,c));

}
