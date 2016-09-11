
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

#include <algorithm>
#include "operation.h"
#include "model.h"


void TWorld::SlopeStability()
{
    if(!SwitchErosion)
    {
        return;
    }
    if(!SwitchSlopeStability)
    {
        return;
    }if(this->InfilMethod == INFIL_NONE)
    {
        return;
    }

    SafetyFactor();

    return;
}

void TWorld::SlopeFailure()
{
    if(!SwitchSlopeStability)
    {
        return;
    }if(!SwitchSolidPhase)
    {
        return;
    }if(!SwitchSlopeFailure)
    {
        return;
    }if(!SwitchErosion)
    {
        return;
    }if(this->InfilMethod == INFIL_NONE)
    {
        return;
    }

    InitiateDebrisFlow();

    return;
}

void TWorld::SafetyFactor()
{

    FOR_ROW_COL_MV
    {
        tmb->Drc = 0;
        DFUnstable->Drc = 0;
        DFInitiationHeight->Drc = 0;
    }

    FOR_ROW_COL_MV
    {
        DFSFIterations->Drc = -1;
        DEMIterate->Drc = DEMOriginal->Drc;
        DFSoilDepth->Drc = SoilDepth1->Drc;
        if(SwitchTwoLayer)
        {
              DFSoilDepth->Drc += SoilDepth2->Drc;
        }
        DFSurfaceWaterHeight->Drc = UF2D_f->Drc/(_dx*_dx);
        DFSoilCohesion->Drc = Cohesion->Drc;
        if(this->InfilMethod != INFIL_NONE)
        {
            DFWaterHeight->Drc = L1->Drc;//ThetaI1->Drc * SoilDepth1->Drc;
            if(SwitchTwoLayer)
            {
                 DFWaterHeight->Drc += L2->Drc;//ThetaI2->Drc * SoilDepth2->Drc;
            }
        }
        DFWaterSuction->Drc = 0;
        DFPlantCohesion->Drc = RootCohesion->Drc;
        DFPlantPressure->Drc = 0.0;
        DFUnstable->Drc = 0;
        DFThreshold->Drc = SFMIN;
        DFThreshold1->Drc = SFMAX;
    }

    bool resolved = false;
    int iter = 0;
    while(!resolved)
    {
        FOR_ROW_COL_MV
        {
            DEMIterate->Drc -= tmb->Drc;
            DFWaterHeight->Drc = DFSoilDepth->Drc > 0? DFWaterHeight->Drc * (DFSoilDepth->Drc-tmb->Drc)/DFSoilDepth->Drc : 0.0;
            DFSoilDepth->Drc -= tmb->Drc;

            tma->Drc = DFUnstable->Drc;
            tmb->Drc = 0;
        }

        CalculateSafetyFactor(DEMIterate,DFSoilDepth,DFSurfaceWaterHeight,DFSoilCohesion,DFSoilInternalFrictionAngle,DFWaterHeight,DFWaterSuction,DFSoilDensity,DFPlantCohesion,DFPlantPressure,DFSafetyFactor,DFThreshold,DFThreshold1, tmb,DFInitiationHeight);


        resolved = true;
        int n_init = 0;
        double htotal = 0;
        FOR_ROW_COL_MV
        {
            if(tmb->Drc > 0)
            {
                DFUnstable->Drc = 1;
                resolved = false;
                n_init++;
                htotal += tmb->Drc;
            }

            if(tma->Drc == 1)
            {
                DFUnstable->Drc = 1;
            }

            if(DFUnstable->Drc == 1  && tma->Drc == 0)
            {
                DFSFIterations->Drc = iter + 1;
            }

            DFInitiationHeight->Drc += tmb->Drc;
        }

        if(iter == 0)
        {
            FOR_ROW_COL_MV
            {
                tmc->Drc = DFSafetyFactor->Drc;
            }
        }
        if(!SwitchSlopeFailure)
        {
            break;
        }
        iter ++;

    }

    FOR_ROW_COL_MV
    {
        DFSafetyFactor->Drc = tmc->Drc;
    }
}

bool TWorld::OUTORMV(int r, int c)
{
    if(r>=0 && r<_nrRows && c>=0 && c<_nrCols)
    {
        if(!pcr::isMV(LDD->data[r][c]))
        {
            return false;
        }
    }
    return true;

}

void TWorld::CalculateSafetyFactor(cTMap * _DEM,cTMap * _SoilDepth,
                                   cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                   cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                   cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                   cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                   cTMap * _SafetyFactor,cTMap * _Threshold,cTMap * _Threshold1,cTMap * _InititationHeight,cTMap * _Initiated)
{
    _SafetyFactor->setAllMV();
    _InititationHeight->setAllMV();
    FOR_ROW_COL_MV
    {

        double SlopeX = 0;
        double SlopeY = 0;
        double Slope = 0;

        //DEM
        double dem = _DEM->Drc;

        double nx = 0;
        double ny = 0;

        double demx = 0;
        double demy = 0;

        if(!OUTORMV(r+1,c))
        {
            nx += 1.0;
            demx += _DEM->data[r+1][c] - dem;
        }
        if(!OUTORMV(r-1,c))
        {
            nx += 1.0;
            demx += dem -_DEM->data[r-1][c];
        }

        if(!OUTORMV(r,c+1))
        {
            nx += 1.0;
            demx += _DEM->data[r][c+1] - dem;
        }
        if(!OUTORMV(r,c-1))
        {
            ny += 1.0;
            demx += dem -_DEM->data[r][c-1];
        }

        if(nx == 0)
        {
            SlopeX = 0;
        }else
        {
            SlopeX = (demx/nx)/_dx;
        }
        if(ny == 0)
        {
            SlopeY = 0;
        }else
        {
            SlopeY = (demy/ny)/_dx;
        }

        Slope = fabs(SlopeX) + fabs(SlopeY);

        DFSlope->Drc = Slope;
        double angle = atan(Slope);
        double cosa = cos(angle);
        double sina = sin(angle);
        double tanphi = tan(_InternalFrictionAngle->Drc);


        double t1 = (_SoilCohesion->Drc + _PlantCohesion->Drc + _SoilWaterSuction->Drc)
                    +
                    (
                        cosa * cosa*
                        (
                         _PlantPressure->Drc +
                               (
                                    _SoilDensity->Drc *_SoilDepth->Drc - 1000.0 *(_SoilWaterHeight->Drc )
                                )
                         )*tanphi
                     );

        double t2 = (
                    (_SoilDensity->Drc * _SoilDepth->Drc + _PlantPressure->Drc)
                *sina*cosa
                );

        double sf = (t2 > 0 ? t1/t2 : 1000.0);

        _SafetyFactor->Drc = sf;

    }


    FOR_ROW_COL_MV
    {
        bool nb_df = false;
        if(!OUTORMV(r-1,c)){if(_Initiated->data[r-1][c] > 0){nb_df = true;}}
        if(!OUTORMV(r+1,c)){if(_Initiated->data[r+1][c] > 0){nb_df = true;}}
        if(!OUTORMV(r,c-1)){if(_Initiated->data[r][c-1] > 0){nb_df = true;}}
        if(!OUTORMV(r,c+1)){if(_Initiated->data[r][c+1] > 0){nb_df = true;}}

        if(_SafetyFactor->Drc < _Threshold->Drc || (_SafetyFactor->Drc < _Threshold1->Drc && nb_df))
        {

            double cif = cos(_InternalFrictionAngle->Drc);
            double sif = sin(_InternalFrictionAngle->Drc);
            double d = _dx;
            double dx2 = _dx*_dx;
            double wf = _SoilWaterHeight->Drc/_SoilDepth->Drc;
            double wd = 1000.0;
            double pp = _PlantPressure->Drc;
            double hf = _SoilDepth->Drc - DFSlope->Drc * _dx;
            double pc = _PlantCohesion->Drc;
            double sd = _SoilDensity->Drc;
            double ws = _SoilWaterSuction->Drc;
            double sc = _SoilCohesion->Drc;
            double sf = _Threshold1->Drc;




            //solution for dynamic threshold value
            double t1 = (2.0 *hf *pc *cif + d* pp* sf*cif + 2.0 *hf *sc *cif -
                         d *hf *sd *sf * cif + 2.0 *hf *ws *cif -
                         dx2 *sd *sif + dx2 *wd *wf *sif);
            double t21 = (-2.0 *hf *pc *cif - d*pp*sf*cif - 2.0 *hf *sc *cif +
                    d *hf *sd *sf*cif - 2.0 *hf *ws *cif + dx2 *sd *sif -
                    dx2 *wd *wf *sif);
            double t22 = (pc*cif + sc*cif - d*sf*sd*cif +
                    ws*cif)*(dx2*pc*cif + hf*hf*pc*cif + d*hf*pp*sf*cif +
                    dx2*sc*cif + hf*hf*sc*cif + dx2*ws*cif +
                    hf*hf*ws*cif + dx2*pp*sif);

            double t2 = sqrt(t21*t21 - 4.0*t22);
            double t3 = (2.0*(pc*cif + sc*cif - d*sd*sf*cif + ws*cif));

            double h1 = (t1 + t2)/t3;
            double h2 = (t1 - t2)/t3;

            _InititationHeight->Drc = std::max(0.0,_SoilDepth->Drc - std::max(h1,h2));

            double DF_MinInitiationHeight = 0.1;

            if(_InititationHeight->Drc < DF_MinInitiationHeight)
            {
                _InititationHeight->Drc = 0;
            }
        }else
        {
            _InititationHeight->Drc = 0;
        }
    }


}

/*//Solution for threshold value of 1
double t1 = (2.0 *hf *pc *cif + d* pp *cif + 2.0 *hf *sc *cif -
             d *hf *sd *cif + 2.0 *hf *ws *cif -
             dx2 *sd *sif + dx2 *wd *wf *sif);
double t21 = (-2.0 *hf *pc *cif - d*pp*cif - 2.0 *hf *sc *cif +
        d *hf *sd *cif - 2.0 *hf *ws *cif + dx2 *sd *sif -
        dx2 *wd *wf *sif);
double t22 = (pc*cif + sc*cif - d*sd*cif +
        ws*cif)*(dx2*pc*cif + hf*hf*pc*cif + d*hf*pp*cif +
        dx2*sc*cif + hf*hf*sc*cif + dx2*ws*cif +
        hf*hf*ws*cif + dx2*pp*sif);*/

void TWorld::InitiateDebrisFlow()
{
    if((!(SwitchSolidPhase)) || (!(SwitchSlopeStability)))
    {
        return;
    }

    int n_init = 0;
    FOR_ROW_COL_MV
    {
        if(DFInitiationHeight->Drc > UF_VERY_SMALL)
        {
            n_init ++;
            DFTotalInitiationHeight->Drc += DFInitiationHeight->Drc;

            double h = DFInitiationHeight->Drc;
            double h_old = UF2D_s->Drc/(_dx*_dx);

            //change DEM (flow DEM is altered later through DEMChange)
            DEMOriginal->Drc -= h;
            //DEM->Drc -= h;
            DEMChange->Drc -= h;

            //fluid volume is added
            double fvol = std::min(SoilDepth1->Drc,h) * (ThetaS1->Drc) * _dx * _dx;
            if(SwitchTwoLayer)
            {
                fvol += std::min(SoilDepth2->Drc,std::max(h-SoilDepth1->Drc,0.0)) * (ThetaS2->Drc) * _dx * _dx;
            }

            UF2D_f->Drc += fvol;

            //solid phase volume is added
            double svol = (DFSoilRockFraction->Drc) * std::min(SoilDepth1->Drc,h) * (1.0-ThetaS1->Drc) * _dx * _dx;
            if(SwitchTwoLayer)
            {
                svol += (DFSoilRockFraction->Drc) * std::min(SoilDepth2->Drc,std::max(h-SoilDepth1->Drc,0.0)) * (1.0-ThetaS2->Drc) * _dx * _dx;
            }

            //and flow properties updated
            UF2D_su->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_su->Drc)/(svol + UF2D_s->Drc) : 0.0;
            UF2D_sv->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_sv->Drc)/(svol + UF2D_s->Drc) : 0.0;
            UF2D_d->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_d->Drc + DFSoilDensity->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_d->Drc;
            UF2D_ifa->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_ifa->Drc + DFSoilInternalFrictionAngle->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_ifa->Drc;
            UF2D_rocksize->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_rocksize->Drc + DFSoilRockSize->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_rocksize->Drc;

            UF2D_s->Drc += svol;

            //suspended sediment is added
            double SedMass = (1.0-DFSoilRockFraction->Drc) * std::min(SoilDepth1->Drc,h) * DFSoilDensity->Drc * (1.0-ThetaS1->Drc) * _dx * _dx;
            if(SwitchTwoLayer)
            {
                SedMass += (1.0-DFSoilRockFraction->Drc) * std::min(SoilDepth2->Drc,std::max(h-SoilDepth1->Drc,0.0)) * DFSoilDensity->Drc * (1.0-ThetaS2->Drc) * _dx * _dx;
            }

            double SedMassOld = UF2D_ssm->Drc;
            UF2D_ssm->Drc += SedMass;

            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    UF2D_ssm_D.Drcd += IW_D.Drcd * SedMass;
                }
            }

            //Depth of infiltration layers decreased
            double sd1 = SoilDepth1->Drc;
            SoilDepth1->Drc = std::max(0.0,SoilDepth1->Drc - h);
            if(SwitchTwoLayer)
            {
                 SoilDepth2->Drc = std::max(0.0,(SoilDepth2->Drc) - std::max(h - sd1,0.0));
            }

        }

    }
}
