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
 \file lisUnifiedFlow.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

void TWorld::UF2D_FluidSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_f,cTMap *out_visc)
{
    FOR_ROW_COL_UF2D
    {
        out_f->Drc = _f->Drc;
        out_visc->Drc = _f->Drc;
    }

}

void TWorld::UF2D_SolidSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_s,cTMap * out_d,cTMap *  out_ifa,cTMap *  out_rocksize)
{
    FOR_ROW_COL_UF2D
    {
        out_s->Drc = _s->Drc;
        out_d->Drc = _d->Drc;
        out_ifa->Drc = _ifa->Drc;
        out_rocksize->Drc = _rocksize->Drc;
    }
}


void TWorld::UF1D_FluidSource(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_f,cTMap * out_visc)
{
    FOR_ROW_COL_UF1D
    {
        out_f->Drc = _f->Drc;
        out_visc->Drc = _f->Drc;
    }
}
void TWorld::UF1D_SolidSource(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_s,cTMap * out_d,cTMap *  out_ifa,cTMap *  out_rocksize)
{
    FOR_ROW_COL_UF1D
    {
        out_s->Drc = _s->Drc;
        out_d->Drc = _d->Drc;
        out_ifa->Drc = _ifa->Drc;
        out_rocksize->Drc = _rocksize->Drc;
    }
}

void TWorld::UF_FlowDetachment(double dt)
{
    cTMap*_dem = UF2D_DEM;
    cTMap*_ldd = UF1D_LDD;

    FOR_ROW_COL_UF2D
    {
        if(!SwitchUseGrainSizeDistribution)
        {
            UF_FlowDetachment(dt,r,c,-1,false);
        }else
        {
            FOR_GRAIN_CLASSES
            {
                UF_FlowDetachment(dt,r,c,d,false);
            }
        }
    }
    FOR_ROW_COL_UF1D
    {
        if(!SwitchUseGrainSizeDistribution)
        {
            UF_FlowDetachment(dt,r,c,-1,true);
        }else
        {
            FOR_GRAIN_CLASSES
            {
                UF_FlowDetachment(dt,r,c,d,true);
            }
        }
    }


}

void TWorld::UF_FlowDetachment(double dt, int r, int c,int d, bool channel)
{
    //set maps for this grain class
    cTMap * TBLTC;
    cTMap * TSSTC;
    cTMap * TBL;
    cTMap * TSS;
    cTMap * TW;
    double TSettlingVelocity;

    if(d == -1)
    {
        TBLTC = channel? UF1D_bltc : UF2D_bltc;
        TSSTC = channel? UF1D_sstc : UF2D_sstc;
        TBL = channel? UF1D_blm : UF2D_blm;
        TSS = channel? UF1D_ssm : UF2D_ssm;
        TSettlingVelocity = SettlingVelocity->Drc;
        TW = unity;
    }else
    {
        TBLTC = channel? UF1D_bltc_D.at(d) : UF2D_bltc_D.at(d);
        TSSTC = channel? UF1D_sstc_D.at(d) : UF2D_sstc_D.at(d);
        TBL = channel? UF1D_blm_D.at(d) : UF2D_blm_D.at(d);
        TSS = channel? UF1D_ssm_D.at(d) : UF2D_ssm_D.at(d);
        TW = W_D.at(d);
        TSettlingVelocity = settlingvelocities.at(d);
    }

    double surface = channel? DX->Drc * UF1D_LDDw->Drc : DX->Drc*ChannelAdj->Drc;
    double watervol = (channel? UF1D_f->Drc: UF2D_f->Drc);
    double blm = TBL->Drc;
    double ssm = TSS->Drc;
    double ssconc = watervol > UF_VERY_SMALL? ssm/watervol : 0.0;
    double blconc = watervol > UF_VERY_SMALL? blm/watervol : 0.0;
    double hf = (channel? UF1D_f->Drc: UF2D_f->Drc)/surface;
    double hs = (channel? UF1D_s->Drc: UF2D_s->Drc)/surface;
    double velocity = channel? (UF1D_fu->Drc):(sqrt(UF2D_fu->Drc*UF2D_fu->Drc + UF2D_fv->Drc*UF2D_fv->Drc));
    double discharge = hf * _dx * velocity;
    double blhf = std::min(0.05,hf);
    double blwatervol = blhf*surface;
    double bldischarge = blhf *  _dx * velocity;

    //calculate tranport capacity for bed load and suspended load
    TBLTC->Drc = UnifiedFlowTransportCapacity(r,c,d,channel,true);
    TSSTC->Drc = UnifiedFlowTransportCapacity(r,c,d,channel,false);


    if(!(hf > UF_VERY_SMALL) || surface < UF_VERY_SMALL)
    {
        //dump everything since there is hardly any water volume
        UF_SoilAdd(r,c,d,TBL->Drc,true);
        UF_SoilAdd(r,c,d,TSS->Drc,false);
        TBL->Drc=0;
        TSS->Drc=0;
        TBLTC->Drc=0;
        TSSTC->Drc=0;

    }else
    {
       ////Suspended sediment deposition
       //first check if sediment goes to suspended sediment layer or to bed layer
       double tobl = 0;
       double toss = 0;
       double TransportFactor;

       //deposition based on settling velocity
       TransportFactor = (1.0-exp(-dt*TSettlingVelocity/hf)) * watervol;

       double maxTC = std::max(TSSTC->Drc - ssconc,0.0) ;
       // positive difference: TC deficit becomes detachment (ppositive)
       double minTC = std::min(TSSTC->Drc - ssconc,0.0) ;

       tobl = TransportFactor * minTC;

       tobl = std::max(tobl,-TSS->Drc);
       TBL->Drc -= tobl;
       TSS->Drc += tobl;

       ////Suspended sediment detachment
       //erosion values based on settling velocity
       TransportFactor = dt*TSettlingVelocity * DX->Drc * SoilWidthDX->Drc;

       //correct detachment for grass strips, hard surfaces and houses
       double detachment = TW->Drc * maxTC * TransportFactor;
       if (GrassFraction->Drc > 0)
          detachment = (1-GrassFraction->Drc) * detachment;
       detachment = (1-StoneFraction->Drc) * detachment ;
       if (SwitchHardsurface)
          detachment = (1-HardSurface->Drc) * detachment;
       if (SwitchHouses)
          detachment = (1-HouseCover->Drc)* detachment;

       //check how much of the potential detachment can be detached from soil layer
       detachment = UF_SoilTake(r,c,d,detachment,channel,false);

       //### sediment balance
       TSS->Drc += detachment;
       double sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*hf;
       if(sssmax < TSS->Drc)
       {
           TBL->Drc +=(TSS->Drc - sssmax);
           TSS->Drc = sssmax;
       }

       ////recalculate bed load concentration
       blm = TBL->Drc;
       ssm = TSS->Drc;
       ssconc = watervol > UF_VERY_SMALL? ssm/watervol : 0.0;
       blconc = watervol > UF_VERY_SMALL? blm/watervol : 0.0;


       ////bed load deposition and detachment
       //### calc concentration and net transport capacity
       maxTC = std::max(TBLTC->Drc - blconc,0.0);
       // positive difference: TC deficit becomes detachment (ppositive)
       minTC = std::min(TBLTC->Drc - blconc,0.0);
       // negative difference: TC surplus becomes deposition (negative)
       // unit kg/m3

       ////Bed load sediment detachment
       //### detachment
       TransportFactor = dt*TSettlingVelocity * DX->Drc *SoilWidthDX->Drc;
       // detachment can only come from soil, not roads (so do not use flowwidth)
       // units s * m/s * m * m = m3

       detachment = TW->Drc * maxTC * TransportFactor;
       // unit = kg/m3 * m3 = kg
       detachment = std::min(detachment, maxTC * bldischarge*dt);
       // cannot have more detachment than remaining capacity in flow
       // use discharge because standing water has no erosion

       detachment = UF_SoilTake(r,c,d,detachment,channel,true);

       // IN KG/CELL

       //DETFlow->Drc = (1-Snowcover->Drc) * DETFlow->Drc ;
       // TODO: CHECK THIS no flow detachment on snow
       //is there erosion and sedimentation under the snowdeck?

       //### deposition
       if (BLDepthFlood->Drc > MIN_HEIGHT)
          TransportFactor = (1-exp(-dt*TSettlingVelocity/hf)) * watervol;
       else
          TransportFactor = 1*watervol;
       // if settl velo is very small, transportfactor is 0 and depo is 0
       // if settl velo is very large, transportfactor is 1 and depo is max

       //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
       // deposition can occur on roads and on soil (so use flowwidth)

       double deposition = minTC * TransportFactor;
       // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

       if (SwitchLimitDepTC)
          deposition = std::max(deposition, minTC *blwatervol);
       // cannot be more than sediment above capacity
       deposition = std::max(deposition, -TBL->Drc);
       // cannot have more depo than sediment presen

       UF_SoilAdd(r,c,d,deposition,channel);

       // IN KG/CELL
       TBL->Drc += detachment;
       TBL->Drc += deposition;

    }

}


//transport capacity
double TWorld::UnifiedFlowTransportCapacity(int r, int c, int _d, bool channel, bool bedload)
{
    double velocity = channel? (UF1D_fu->Drc):(sqrt(UF2D_fu->Drc*UF2D_fu->Drc + UF2D_fv->Drc*UF2D_fv->Drc));
    double surface = channel? DX->Drc * UF1D_LDDw->Drc : DX->Drc*ChannelAdj->Drc;
    double hf = (channel? UF1D_f->Drc: UF2D_f->Drc)/surface;
    double v = velocity;
    if(hf < UF_VERY_SMALL || v < UF_VERY_SMALL)
    {
        return 0.0;
    }

    //use overland flow equations
    if(hf < 0.25)
    {
        if(bedload)
        {
            return 0.0;
        }else
        {
            //use Govers transport capacity equation
            if(_d == -1)
            {
                double slope = channel? UF1D_LDDs->Drc : std::max(UF2D_SlopeX->Drc,UF2D_SlopeY->Drc);
                CG->Drc = pow((D50->Drc+5)/0.32, -0.6);
                DG->Drc = pow((D50->Drc+5)/300, 0.25);
                //### Calc transport capacity
                double omega = 0;

                    omega = 100.0* v*slope;

                // V in cm/s in this formula assuming grad is SINE
                double omegacrit = 0.4;
                // critical unit streampower in cm/s
                return std::min(MAXCONC, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
                // not more than 2650*0.32 = 848 kg/m3

            //use the Hairsine and Rose transport capacity equation
            }else
            {
                double slope = channel? UF1D_LDDs->Drc : std::max(UF2D_SlopeX->Drc,UF2D_SlopeY->Drc);
                double om =  100.0* v*slope;
                double omcr = 0.4;
                double tc =  W_D.at(_d)->Drc * (1.0/settlingvelocities.at(_d))*(1.0 * 0.013/9.81) * (2650.0/(2650.0 - 1000.0)) * ( std::max(0.0, (om - omcr))/WHrunoff->Drc) ;
                return std::min(MAXCONC,tc);
            }

        }
    }
    //use river-based equations
    if(_d == -1)
    {
        if(bedload)
        {
            //van rijn simple bed load
            double ps = 2400.0;
            double pw = 1000.0;
            double ucr;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            if( d50m < 0.005)
            {
               ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* hf/d90m);
            }else
            {
               ucr  = 0.19 * pow(d50m, 0.6) * log10(4.0* hf/d90m);
            }
            double me = std::max((v - ucr)/(sqrt(UF_Gravity * d50m * ((ps/pw)- 1.0))),0.0);
            double qs = 0.005 * ps*v *hf * pow(d50m/hf,1.2) * pow(me, 2.4);
            double tc =  qs/ (v * hf );
            return std::max(std::min( tc,MAXCONC ),0.0);

        }else
        {
            //van rijn simple suspended load
            double ucr;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400.0;
            double pw = 1000.0;
            double mu = 1.0;
            if( d50m < 0.0005)
            {
               ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* hf/d90m);
            }else
            {
               ucr  = 8.5 * pow(d50m, 0.6) * log10(4.0* hf/d90m);
            }
            double me = std::max((v - ucr)/sqrt(UF_Gravity * d50m * (ps/pw - 1)),0.0);
            double ds = d50m * UF_Gravity * ((ps/pw)-1)/(mu*mu);
            double qs = hf * 0.008 * ps*v * d50m * pow(me, 2.4) * pow(ds, -0.6);

            double tc =  qs/ (v * hf);
            return std::max(std::min(tc,MAXCONC),0.0);


        }

    }else
    {
        if(bedload)
        {
            //wu wang & Jia bed load
            double slope = channel? UF1D_LDDs->Drc : std::max(UF2D_SlopeX->Drc,UF2D_SlopeY->Drc);
            double ps = 2400.0;
            double pw = 1000.0;
            double h = hf;
            double n = std::max(N->Drc,0.001);
            double na = pow(graindiameters.at(_d)/100000.0,(1.0/6.0))/20.0;
            double phk = 0;
            double pek = 0;
            FOR_GRAIN_CLASSES
            {
                phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }
            double ppk = 1;
            ppk = pow(phk/pek,0.6);
            if(pek  == 0)
            {
                return 0;
            }

            double dh = (ChannelAdj->Drc *h)/(ChannelAdj->Drc + 2* h);
            double css = 0.03* (ps - pw) * (graindiameters.at(_d)/1000000.0) * ppk;

            double qs = 0.0053 *pow(std::max(pow(na/n,1.5)*((pw * dh * 9.81 * 0.1 * slope/css) -1.0 ), 0.0),2.2);
            qs = qs * sqrt((ps/pw - 1)*9.81*pow(graindiameters.at(_d)/1000000.0,3.0));

            double tc = ps * ChannelAdj->Drc * qs/ (v * hf*ChannelAdj->Drc);
            return std::max(std::min(tc,MAXCONCBL ),0.0);


        }else
        {
            //wu wang & Jia suspended load
            if(SSD_D.at(_d)->Drc < 0.004)
            {
                return 0;
            }
            double slope = channel? UF1D_LDDs->Drc : std::max(UF2D_SlopeX->Drc,UF2D_SlopeY->Drc);
            double ps = 2400.0;
            double pw = 1000.0;
            double h = hf;
            double phk = 0;
            double pek = 0;
            double sv = settlingvelocities.at(_d);
            double gd = graindiameters.at(_d)/1000000.0;
            FOR_GRAIN_CLASSES
            {
                phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }

            double ppk = 1;
            ppk = pow(phk/pek,0.6);
            if(pek == 0)
            {
                return 0;
            }

            double dh = (ChannelAdj->Drc *h)/(ChannelAdj->Drc + 2.0* h);

            double css = 0.03* (ps - pw) * (gd) * ppk;

            double qs = 0.0000262 *pow(std::max(( pw * 0.01 * h /css) - 1.0, 0.0)* v/(sqrt(sv)),2.2);
            qs =  qs * 1 * sqrt((ps/pw - 1)*9.81*pow(gd,3.0));

            double tc = ps * ChannelAdj->Drc * qs/ (v * hf * ChannelAdj->Drc);
            return std::max(std::min(tc,MAXCONC ),0.0);



        }
    }
}


//active entrainment
double TWorld::UnifiedFlowActiveEntrainment(double hf, double vf,double hs, double vs, double visc, int d)
{


}

void TWorld::UF_FlowEntrainment(double dt)
{





}
void TWorld::UF_FlowEntrainment(double dt, int r, int c,int d, bool channel)
{


}

double TWorld::UF_SoilTake(int r, int c, int d, double potential, bool channel,bool bedload)
{

    double detachment = potential;

    if (GrassFraction->Drc > 0)
       detachment = (1-GrassFraction->Drc) * detachment;
    // no flow detachment on grass strips

    // Detachment edxceptions:
    detachment = (1-StoneFraction->Drc) * detachment;
    // no flow detachment on stony surfaces

    if (SwitchHardsurface)
       detachment = (1-HardSurface->Drc) * detachment;
    // no flow detachment on hard surfaces

    if (SwitchHouses)
       detachment = (1-HouseCover->Drc)*detachment;
    // no flow det from house roofs

    //check how much of the potential detachment can be detached from soil layer
    detachment = DetachMaterial(r,c,d,channel,!channel, bedload, detachment);
    if(channel)
    {
        UF1D_Det->Drc += detachment;
    }else
    {
        UF2D_Det->Drc += detachment;
    }

}

void TWorld::UF_SoilAdd(int r, int c, int d, double deposition, bool channel)
{
    if(SwitchUseMaterialDepth)
    {
        StorageDep->Drc += deposition;
        if(SwitchUseGrainSizeDistribution)
        {
              StorageDep_D.Drcd += -deposition;
        }
    }

    if(channel)
    {
        UF1D_Dep->Drc += deposition;
    }else
    {
        UF2D_Dep->Drc += deposition;
    }

}