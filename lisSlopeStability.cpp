
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
    sfset =true;

    FOR_ROW_COL_MV
    {
        tmb->Drc = 0;
        DFUnstable->Drc = 0;
        DFInitiationHeight->Drc = 0;

        if(SwitchBedrock)
        {
            tmd->Drc = 0;
            DFUnstable2->Drc = 0;
            DFInitiationHeight2->Drc = 0;
        }

    }



    //get all the nececcary input from the OpenLISEM model data, and put in in temporary maps
    //this way we can adapt stuff without breaking the rest of the model
    FOR_ROW_COL_MV
    {


        /*DFForcingDemand->Drc = 0.0;
        DFForcingCapacity->Drc = 0.0;
        DFForcing->Drc = 0.0;
        DFForcingAdded->Drc = 0.0;*/

        double area = _dx*_dx;
        DFSFIterations->Drc = -1;
        DEMIterate->Drc = DEMOriginal->Drc;
        DFSoilDepth->Drc = SoilDepth1->Drc;
        if(SwitchTwoLayer)
        {
              DFSoilDepth->Drc += SoilDepth2->Drc;
        }
        double soildepth0 = 0;
        double theta0 = 0;
        if(SwitchEntrainment && INCLUDE_DEPOSITS_STABILITY)
        {
              theta0 = (SoilRockMaterial->Drc)? SoilRockWater->Drc/SoilRockMaterial->Drc: 0.0;
              soildepth0 = SoilRockMaterial->Drc /area;
              DFSoilDepth->Drc += soildepth0;
        }

        //DFSoilDepth->Drc *= st_sdCalibration;

        DFSurfaceWaterHeight->Drc = 0;//UF2D_f->Drc/(_dx*_dx);
        DFSoilCohesion->Drc = DFSoilDepth->Drc > 0? (SwitchSeismic? std::min(1.0,std::max(0.0,(1.0 - StrengthLoss->Drc))) : 1.0) * 1000.0*((Cohesion->Drc * st_scCalibration * (DFSoilDepth->Drc-soildepth0)) + (10.0 * soildepth0))/(DFSoilDepth->Drc) : 0.0;
        if(SwitchBedrock)
        {
            DFSoilCohesion2->Drc = (SwitchSeismic? std::min(1.0,std::max(0.0,(1.0 - StrengthLoss2->Drc))) : 1.0) * DFSoilCohesion2Org->Drc;
        }

        DFSoilInternalFrictionAngle->Drc = (SwitchSeismic? std::min(1.0,std::max(0.0,(1.0 - StrengthLoss->Drc))) : 1.0) *DFSoilInternalFrictionAngleOrg->Drc;
        if(SwitchBedrock)
        {
            DFSoilInternalFrictionAngle2->Drc = (SwitchSeismic? std::min(1.0,std::max(0.0,(1.0 - StrengthLoss2->Drc))) : 1.0) *DFSoilInternalFrictionAngle2Org->Drc;
        }

        if(this->InfilMethod != INFIL_NONE)
        {
            DFWaterHeight->Drc = L1->Drc * (ThetaS1->Drc - ThetaI1->Drc) + ThetaI1->Drc * SoilDepth1->Drc;
            if(SwitchTwoLayer)
            {
                 DFWaterHeight->Drc += L2->Drc * (ThetaS2->Drc - ThetaI2->Drc) + ThetaI2->Drc * SoilDepth2->Drc;
            }
            if(SwitchEntrainment && INCLUDE_DEPOSITS_STABILITY)
            {
                DFWaterHeight->Drc += SoilRockWater->Drc /area;
            }
        }
        double saturation
                = std::max(0.0,std::min(1.0,(DFWaterHeight->Drc/DFSoilDepth->Drc)));
        double ts = 1;
        if(this->InfilMethod != INFIL_NONE)
        {
            ts =(ThetaS1->Drc);
            if(SwitchTwoLayer)
            {
                if((SoilDepth1->Drc + SoilDepth2->Drc) > 0)
                {
                    ts = (ThetaS1->Drc * SoilDepth1->Drc + ThetaS2->Drc * SoilDepth2->Drc)/(SoilDepth1->Drc + SoilDepth2->Drc);
                }
            }
            if(SwitchEntrainment && INCLUDE_DEPOSITS_STABILITY)
            {
                ts = DFSoilDepth->Drc > 0? ((ts * (DFSoilDepth->Drc-soildepth0)) + (UF_FLOWCOMPACTION_DEPOSITIONPOROSITY * soildepth0))/(DFSoilDepth->Drc) : 0.0;
            }
        }
        saturation = (ts > 0.1 ?saturation / ts : 0.0);

        if(this->InfilMethod != INFIL_NONE)
        {
            DFWaterSuction->Drc = Psi1->Drc;
            if(SwitchTwoLayer)
            {
                if((SoilDepth1->Drc + SoilDepth2->Drc) > 0)
                {
                    DFWaterSuction->Drc = (Psi1->Drc * SoilDepth1->Drc + Psi2->Drc * SoilDepth2->Drc)/(SoilDepth1->Drc + SoilDepth2->Drc);
                }
            }
            DFWaterSuction->Drc = DFWaterSuction->Drc * std::max(0.0,std::min(1.0,(3.0 - 3.0 *saturation)));
        }

        DFPlantCohesion->Drc = RootCohesion->Drc;
        DFPlantPressure->Drc = 0.0;
        DFUnstable->Drc = 0;
        DFThreshold->Drc = SFMIN;
        DFThreshold1->Drc = SFMAX;

        if(SF_Calibrate_First)
        {
            DFSFCalibration->Drc = 1.0;
            if(SwitchBedrock)
            {
                DFSFCalibration2->Drc = 1.0;
            }
        }
    }


    if(SwitchBedrock)
    {
        CalculateBedrockDepth(DEMIterate,DFSoilDepth, DFSoilDepth2);
    }


    if(SwitchUpslopeForcing)
    {
        CalculateSlopeForcing(DEMIterate,DFSoilDepth,DFSurfaceWaterHeight,DFSoilCohesion,DFSoilInternalFrictionAngle,DFWaterHeight,DFWaterSuction,DFSoilDensity,DFPlantCohesion,DFPlantPressure,DFSFCalibration, DFForcing, DFForcingUp,SwitchSeismic?PGACurrent : DFZERO);

        if(SwitchBedrock && SF_BedRock_First)
        {
            CalculateSlopeForcing(DEMIterate,DFSoilDepth2,DFSurfaceWaterHeight,DFSoilCohesion2,DFSoilInternalFrictionAngle2,DFZERO,DFZERO,DFSoilDensity2,DFZERO,DFAddedPressure,DFSFCalibration2, DFForcing2, DFForcingUp2,SwitchSeismic?PGACurrent : DFZERO );
        }

    }

    //store safety factor for display
    FOR_ROW_COL_MV
    {

        if(SwitchUpslopeForcing)
        {
            MaximumUpslopeForcing->Drc = std::max(MaximumUpslopeForcing->Drc,DFForcing->Drc + SwitchBedrock? DFForcing2->Drc : 0.0);
        }
        if(SwitchDownslopeForcing)
        {
            MinimumDownslopeForcing->Drc = std::min(MinimumDownslopeForcing->Drc,DFForcingUp->Drc + SwitchBedrock? DFForcingUp2->Drc : 0.0);
        }

    }


    //Since slope failure at one cell influences the next, we iterate untill we end up with a stable situation
    //the function CalculateSafetyFactor does the actual calculations for safety factor and possible slope failure depth
    //we then adapt the slope and do it again, untill nothing more fails

    bool resolved = false;
    int iter = 0;
    while(!resolved)
    {





        FOR_ROW_COL_MV
        {
            DEMIterate->Drc -= tmb->Drc;

            DFWaterHeight->Drc = DFSoilDepth->Drc > 0? DFWaterHeight->Drc * (DFSoilDepth->Drc-tmb->Drc)/DFSoilDepth->Drc : 0.0;
            DFSoilDepth->Drc -= tmb->Drc;
            //DFSoilDepth2->Drc -= tmd->Drc;

            tma->Drc = DFUnstable->Drc;
            tmb->Drc = 0;
            if(iter == 0)
            {
                tmc->Drc = 0;
            }

            if(SwitchBedrock)
            {
                DEMIterate->Drc -= tmd->Drc;
                tme->Drc = DFUnstable2->Drc;
            }
            tmd->Drc = 0;
        }

        CalculateSafetyFactors(DEMIterate,DFSoilDepth,DFSurfaceWaterHeight,DFSoilCohesion,DFSoilInternalFrictionAngle,DFWaterHeight,DFWaterSuction,DFSoilDensity,DFPlantCohesion,DFPlantPressure,DFSafetyFactor,DFThreshold,DFThreshold1, tmb,DFInitiationHeight,DFSFCalibration,DFForcing,DFForcingUp,SwitchSeismic?PGACurrent : DFZERO);

        if(SwitchBedrock  && SF_BedRock_First)
        {

            CalculateBedrockDepth(DEMIterate,DFSoilDepth, DFSoilDepth2);

            FOR_ROW_COL_MV
            {
                DFAddedPressure->Drc = DFPlantPressure->Drc +(DFSoilDensity->Drc *(DFSoilDepth->Drc - DFWaterHeight->Drc) + 1000.0 *(DFWaterHeight->Drc));
            }

            if(SwitchBedrock && SF_BedRock_First && (iter % 100) == 0)
            {


                CalculateSlopeForcing(DEMIterate,DFSoilDepth2,DFSurfaceWaterHeight,DFSoilCohesion2,DFSoilInternalFrictionAngle2,DFZERO,DFZERO,DFSoilDensity2,DFZERO,DFAddedPressure,DFSFCalibration2, DFForcing2, DFForcingUp2,SwitchSeismic?PGACurrent : DFZERO );
            }

            CalculateSafetyFactors(DEMIterate,DFSoilDepth2,DFSurfaceWaterHeight,DFSoilCohesion2,DFSoilInternalFrictionAngle2,DFZERO,DFZERO,DFSoilDensity2,DFZERO,DFAddedPressure,DFSafetyFactor2,DFThreshold,DFThreshold1, tmd,DFInitiationHeight2,DFSFCalibration2,DFForcing2,DFForcingUp2,SwitchSeismic?PGACurrent : DFZERO,false,0.8);
        }

        //CALIBRATE INITIAL STABILIsTY
        //if selected by user, all cells are forced to be at least stable
        //Safety factor of cells that were unstable increases up to the threshold plus a certain margin
        //stable cells get a calibration factor of 1, so nothing happens.

        if(SF_Calibrate_Initial)
        {
            if(SF_Calibrate_First)
            {
                SF_Calibrate_First = false;

                FOR_ROW_COL_MV
                {
                    if(DFInitiationHeight->Drc > 0 || DFSafetyFactor->Drc < DFThreshold->Drc * (1.0+SF_Calibrate_Margin))
                    {
                        DFSFCalibration->Drc = (DFThreshold->Drc * (1.0+SF_Calibrate_Margin *(DFSafetyFactor->Drc/(DFThreshold->Drc*(1.0+SF_Calibrate_Margin))) ))/DFSafetyFactor->Drc ;
                        tmb->Drc = 0;
                        DFSafetyFactor->Drc = DFSafetyFactor->Drc * DFSFCalibration->Drc;
                    }

                    if(SwitchBedrock)
                    {
                        if(DFInitiationHeight2->Drc > 0 || DFSafetyFactor2->Drc < DFThreshold->Drc * (1.0+SF_Calibrate_Margin))
                        {
                            DFSFCalibration2->Drc = (DFThreshold->Drc * (1.0+SF_Calibrate_Margin *(DFSafetyFactor2->Drc/(DFThreshold->Drc*(1.0+SF_Calibrate_Margin))) ))/DFSafetyFactor2->Drc ;
                            tmd->Drc = 0;
                            DFSafetyFactor2->Drc = DFSafetyFactor2->Drc * DFSFCalibration2->Drc;
                        }
                    }
                }

            }
        }


        //check if there were any slope failures in this iteration
        resolved = true;

        int n_init = 0;
        double htotal = 0;
        FOR_ROW_COL_MV
        {
            if(SwitchBedrock)
            {
                if(tmd->Drc > 0)
                {
                    //SF_BedRock_First = false;

                    tmb->Drc = DFSoilDepth->Drc;
                }
            }

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

        if(SwitchBedrock)
        {
            FOR_ROW_COL_MV
            {
                if(tmd->Drc > 0)
                {
                    DFUnstable2->Drc = 1;
                    resolved = false;
                    n_init++;
                    htotal += tmb->Drc;
                }

                if(tme->Drc == 1)
                {
                    DFUnstable2->Drc = 1;
                }

                if(DFUnstable2->Drc == 1  && tme->Drc == 0)
                {
                    DFSFIterations->Drc = iter + 1;
                }

                DFInitiationHeight2->Drc += tmd->Drc;
            }


        }
        if(iter == 0)
        {
            FOR_ROW_COL_MV
            {
                tmc->Drc = DFSafetyFactor->Drc;
                if(SwitchBedrock)
                {
                    tmf->Drc = DFSafetyFactor2->Drc;
                }
            }
        }


        //if slope failure is off anyway, we can get out of the loop anyway!
        if(!SwitchSlopeFailure)
        {
            SF_BedRock_First = true;
            break;
        }

        iter ++;

        if(iter > 200)
        {
            break;
        }

    }


    //store safety factor for display
    FOR_ROW_COL_MV
    {

        DFSafetyFactor->Drc = tmc->Drc;
        if(SwitchBedrock)
        {
            DFSafetyFactor2->Drc = tmf->Drc;
        }
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

void TWorld::CalculateBedrockDepth(cTMap * _DEM,cTMap * _SoilDepth,
                                   cTMap * _BedrockDepth)
{

    FOR_ROW_COL_MV
    {

        //get the actual slope for this cell

        double Slope = 0;
        double SlopeX = 0;
        double SlopeY = 0;
        //DEM
        double dem = _DEM->Drc;

        double nx = 0;
        double ny = 0;

        double demx = 0;
        double demy = 0;

        double min_surr = dem;

        if(!OUTORMV(r+1,c))
        {
            nx += 1.0;
            demx += _DEM->data[r+1][c] - dem;
            min_surr = std::min(min_surr,_DEM->data[r+1][c]);
        }
        if(!OUTORMV(r-1,c))
        {
            nx += 1.0;
            demx += dem -_DEM->data[r-1][c];
            min_surr = std::min(min_surr,_DEM->data[r-1][c]);
        }

        if(!OUTORMV(r,c+1))
        {
            ny += 1.0;
            demy += _DEM->data[r][c+1] - dem;
            min_surr = std::min(min_surr,_DEM->data[r][c+1]);
        }
        if(!OUTORMV(r,c-1))
        {
            ny += 1.0;
            demy += dem -_DEM->data[r][c-1];
            min_surr = std::min(min_surr,_DEM->data[r][c-1]);
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

        Slope = std::min(_dx,fabs(SlopeX) + fabs(SlopeY));

        _BedrockDepth->Drc = std::max(_dx* 0.5,0.5 * std::max(0.0,dem - min_surr));// + Slope * _dx;//std::max(_dx,Slope * 0.5 *_dx );

       }



}

void TWorld::CalculateSlopeForcing(cTMap * _DEM,cTMap * _SoilDepth,
                                   cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                   cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                   cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                   cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                   cTMap * _SFCalibration, cTMap * _DFForcing,
                                   cTMap * _DFForcingUp, cTMap * _GPA)
{
    FOR_ROW_COL_MV
    {

        DFForcingDemand->Drc = 0.0;
        DFForcingCapacity->Drc = 0.0;

        //get the actual slope for this cell
        double SlopeX = 0;
        double SlopeY = 0;
        double Slope = 0;

        //DEM
        double dem = _DEM->Drc;

        double nx = 0;
        double ny = 0;

        double demx = 0;
        double demy = 0;

        double min_surr = dem;

        if(!OUTORMV(r+1,c))
        {
            nx += 1.0;
            demx += _DEM->data[r+1][c] - dem;
            min_surr = std::min(min_surr,_DEM->data[r+1][c]);
        }
        if(!OUTORMV(r-1,c))
        {
            nx += 1.0;
            demx += dem -_DEM->data[r-1][c];
            min_surr = std::min(min_surr,_DEM->data[r-1][c]);
        }

        if(!OUTORMV(r,c+1))
        {
            ny += 1.0;
            demy += _DEM->data[r][c+1] - dem;
            min_surr = std::min(min_surr,_DEM->data[r][c+1]);
        }
        if(!OUTORMV(r,c-1))
        {
            ny += 1.0;
            demy += dem -_DEM->data[r][c-1];
            min_surr = std::min(min_surr,_DEM->data[r][c-1]);
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

        Slope = std::max(0.01,_DEM->Drc - min_surr)/_dx;

        DFSlope->Drc = Slope;
        DFSlopeX->Drc = SlopeX;
        DFSlopeY->Drc = SlopeY;

        if(!(std::fabs(SlopeX) + std::fabs(SlopeY) > 0))
        {
            DFSlopeX->Drc = 0.01;
            DFSlopeY->Drc = 0.01;

        }

        double angle = atan(Slope);
        double cosa = cos(angle);
        double sina = sin(angle);
        double tanphi = tan(_InternalFrictionAngle->Drc);

        DFForcingCapacity->Drc = _SFCalibration->Drc * ((_SoilCohesion->Drc + _PlantCohesion->Drc )
                    +
                    (-cosa * sina*(_GPA->Drc/9.81) * (_SoilDensity->Drc *(_SoilDepth->Drc - _SoilWaterHeight->Drc) + 1000.0 *(_SoilWaterHeight->Drc))
                        +
                     cosa * cosa*
                        (
                         _PlantPressure->Drc +
                               (
                                    _SoilDensity->Drc *(_SoilDepth->Drc - _SoilWaterHeight->Drc) + 1000.0 *(_SoilWaterHeight->Drc)

                                )
                         ))*tanphi
                     );

        DFForcingDemand->Drc = (_GPA->Drc/9.81)* (_SoilDensity->Drc * _SoilDepth->Drc + 1000.0 *(_SoilWaterHeight->Drc ))*cosa*cosa
               +(
                    (_SoilDensity->Drc * _SoilDepth->Drc + 1000.0 *(_SoilWaterHeight->Drc )+ _PlantPressure->Drc)
                *sina*cosa
                );

        _DFForcing->Drc = 0.0;
        _DFForcingUp->Drc = 0.0;

        DFForcingAdded->Drc = std::max(0.0,DFForcingDemand->Drc-DFForcingCapacity->Drc);
        DFForcingUpAdded->Drc = std::max(0.0,DFForcingCapacity->Drc-DFForcingDemand->Drc);
    }

    double scale = 9.81;

    bool stable = false;
    int iterf = 0;
    while(!stable)
    {

        stable = true;
        FOR_ROW_COL_MV
        {

            if(DFForcingAdded->Drc > 0 && DFSoilDepth->Drc > 0)
            {
                double dem = DEM->Drc;
                double demx1 = OUTORMV(r,c-1)? dem: DEM->data[r][c-1];
                double demx2 = OUTORMV(r,c+1)? dem: DEM->data[r][c+1];
                double demy1 = OUTORMV(r-1,c)? dem: DEM->data[r-1][c];
                double demy2 = OUTORMV(r+1,c)? dem: DEM->data[r+1][c];

                double forcingcapacity = std::max(0.0,(DFForcingCapacity->Drc-DFForcingDemand->Drc) -DFForcing->Drc);
                DFForcingAdded->Drc = 0.5*std::max(0.0,(forcingcapacity > 0?1.0:1.0) * (DFForcingAdded->Drc-forcingcapacity));

                double ls = std::sqrt(DFSlopeX->Drc *DFSlopeX->Drc + DFSlopeY->Drc *DFSlopeY->Drc);
                bool reverse = false;
                if((demx1 < demx2 && demx1 < dem) && !OUTORMV(r,c-1))
                {
                    stable = false;
                    double forceflux = DFForcingAdded->Drc *
                            std::max(0.0,std::min(1.0,(std::max(0.0,DFSlopeX->Drc * DFSlopeX->data[r][c-1] + DFSlopeY->Drc * DFSlopeY->data[r][c-1])/(ls*std::sqrt(DFSlopeX->data[r][c-1] *DFSlopeX->data[r][c-1]+ DFSlopeY->data[r][c-1] *DFSlopeY->data[r][c-1])))))*
                            std::max(0.0,std::min(1.0,  DFSlope->Drc)) * //std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                            (std::fabs(DFSlopeX->Drc)/DFSlope->Drc);//  * (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r][c-1])/DFSoilDepth->Drc);
                    DFForcingAdded->data[r][c-1] += forceflux;

                }
                if((demx1 > demx2 && demx2 < dem) && !OUTORMV(r,c+1))
                {
                    stable = false;
                    double forceflux = DFForcingAdded->Drc *
                            std::max(0.0,std::min(1.0,(std::fabs(DFSlopeX->Drc * DFSlopeX->data[r][c+1] + DFSlopeY->Drc * DFSlopeY->data[r][c+1])/(ls*std::sqrt(DFSlopeX->data[r][c+1] *DFSlopeX->data[r][c+1]+ DFSlopeY->data[r][c+1] *DFSlopeY->data[r][c+1])))))*
                            std::max(0.0,std::min(1.0,  DFSlope->Drc)) * //std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                            (std::fabs(DFSlopeX->Drc)/DFSlope->Drc);//  * (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r][c+1])/DFSoilDepth->Drc);
                    DFForcingAdded->data[r][c+1] +=  forceflux;

                }
                if(demy1 < demy2 && demy1 < dem && !OUTORMV(r-1,c))
                {
                    stable = false;
                    double forceflux = DFForcingAdded->Drc *
                            std::max(0.0,std::min(1.0,(std::fabs(DFSlopeX->Drc * DFSlopeX->data[r-1][c] + DFSlopeY->Drc * DFSlopeY->data[r-1][c])/(ls*std::sqrt(DFSlopeX->data[r-1][c] *DFSlopeX->data[r-1][c]+ DFSlopeY->data[r-1][c] *DFSlopeY->data[r-1][c])))))*
                            std::max(0.0,std::min(1.0, DFSlope->Drc)) * //std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                            (std::fabs(DFSlopeY->Drc)/DFSlope->Drc);//* (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r-1][c])/DFSoilDepth->Drc);
                    DFForcingAdded->data[r-1][c] +=  forceflux;

                }
                if(demy1 > demy2 && demy2 < dem && !OUTORMV(r+1,c))
                {
                    stable = false;
                    double forceflux = DFForcingAdded->Drc *
                            std::max(0.0,std::min(1.0,(std::fabs(DFSlopeX->Drc * DFSlopeX->data[r+1][c] + DFSlopeY->Drc * DFSlopeY->data[r+1][c])/(ls*std::sqrt(DFSlopeX->data[r+1][c] *DFSlopeX->data[r+1][c]+ DFSlopeY->data[r+1][c] *DFSlopeY->data[r+1][c])))))*
                            std::max(0.0,std::min(1.0, DFSlope->Drc));// * std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                            (std::fabs(DFSlopeY->Drc)/DFSlope->Drc);//  * (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r+1][c])/DFSoilDepth->Drc);
                    DFForcingAdded->data[r+1][c] +=  forceflux;

                }
                _DFForcing->Drc += DFForcingAdded->Drc ;
                DFForcingAdded->Drc = 0;
            }
        }

        iterf ++;
        if(iterf > 25)
        {
            break;
        }
    }

    if(SwitchDownslopeForcing)
    {
        stable = false;
        iterf = 0;
        while(!stable)
        {

            stable = true;
            FOR_ROW_COL_MV
            {

                if(DFForcingUpAdded->Drc > 0 && DFSoilDepth->Drc > 0)
                {
                    double dem = DEM->Drc;
                    double demx1 = OUTORMV(r,c-1)? dem: DEM->data[r][c-1];
                    double demx2 = OUTORMV(r,c+1)? dem: DEM->data[r][c+1];
                    double demy1 = OUTORMV(r-1,c)? dem: DEM->data[r-1][c];
                    double demy2 = OUTORMV(r+1,c)? dem: DEM->data[r+1][c];

                    double forcingcapacity = std::max(0.0,(DFForcingDemand->Drc-DFForcingCapacity->Drc) -DFForcingUp->Drc);
                    DFForcingUpAdded->Drc = 0.5* std::max(0.0,(forcingcapacity > 0? 1.0:1.0) * (DFForcingUpAdded->Drc-forcingcapacity));

                    double ls = std::sqrt(DFSlopeX->Drc *DFSlopeX->Drc + DFSlopeY->Drc *DFSlopeY->Drc);
                    bool reverse = false;
                    if((demx1 > demx2 && demx1 > dem) && !OUTORMV(r,c-1))
                    {
                        stable = false;
                        double forceflux = DFForcingUpAdded->Drc *
                                std::max(0.0,std::min(1.0,(std::max(0.0,DFSlopeX->Drc * DFSlopeX->data[r][c-1] + DFSlopeY->Drc * DFSlopeY->data[r][c-1])/(ls*std::sqrt(DFSlopeX->data[r][c-1] *DFSlopeX->data[r][c-1]+ DFSlopeY->data[r][c-1] *DFSlopeY->data[r][c-1])))))*
                                std::max(0.0,std::min(1.0, DFSlope->Drc)) * //std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                                (std::fabs(DFSlopeX->Drc)/DFSlope->Drc);//  * (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r][c-1])/DFSoilDepth->Drc);
                        DFForcingUpAdded->data[r][c-1] += forceflux;

                    }
                    if((demx1 < demx2 && demx2 > dem) && !OUTORMV(r,c+1))
                    {
                        stable = false;
                        double forceflux = DFForcingUpAdded->Drc *
                                std::max(0.0,std::min(1.0,(std::fabs(DFSlopeX->Drc * DFSlopeX->data[r][c+1] + DFSlopeY->Drc * DFSlopeY->data[r][c+1])/(ls*std::sqrt(DFSlopeX->data[r][c+1] *DFSlopeX->data[r][c+1]+ DFSlopeY->data[r][c+1] *DFSlopeY->data[r][c+1])))))*
                                std::max(0.0,std::min(1.0, DFSlope->Drc)) * //std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                                (std::fabs(DFSlopeX->Drc)/DFSlope->Drc);//  * (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r][c+1])/DFSoilDepth->Drc);
                        DFForcingUpAdded->data[r][c+1] +=  forceflux;

                    }
                    if(demy1 > demy2 && demy1 > dem && !OUTORMV(r-1,c))
                    {
                        stable = false;
                        double forceflux = DFForcingUpAdded->Drc *
                                std::max(0.0,std::min(1.0,(std::fabs(DFSlopeX->Drc * DFSlopeX->data[r-1][c] + DFSlopeY->Drc * DFSlopeY->data[r-1][c])/(ls*std::sqrt(DFSlopeX->data[r-1][c] *DFSlopeX->data[r-1][c]+ DFSlopeY->data[r-1][c] *DFSlopeY->data[r-1][c])))))*
                                std::max(0.0,std::min(1.0, DFSlope->Drc)) * //std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                                (std::fabs(DFSlopeY->Drc)/DFSlope->Drc);//* (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r-1][c])/DFSoilDepth->Drc);
                        DFForcingUpAdded->data[r-1][c] +=  forceflux;

                    }
                    if(demy1 < demy2 && demy2 > dem && !OUTORMV(r+1,c))
                    {
                        stable = false;
                        double forceflux = DFForcingUpAdded->Drc *
                                std::max(0.0,std::min(1.0,(std::fabs(DFSlopeX->Drc * DFSlopeX->data[r+1][c] + DFSlopeY->Drc * DFSlopeY->data[r+1][c])/(ls*std::sqrt(DFSlopeX->data[r+1][c] *DFSlopeX->data[r+1][c]+ DFSlopeY->data[r+1][c] *DFSlopeY->data[r+1][c])))))*
                                std::max(0.0,std::min(1.0, DFSlope->Drc));// * std::max(0.0, std::min(1.0, 1.0 - (1/(9.81 ) )*(_dx/DFSoilDepth->data[r][c]))) *
                                (std::fabs(DFSlopeY->Drc)/DFSlope->Drc);//  * (std::min(DFSoilDepth->Drc,DFSoilDepth->data[r+1][c])/DFSoilDepth->Drc);
                        DFForcingUpAdded->data[r+1][c] +=  forceflux;

                    }


                    _DFForcingUp->Drc += DFForcingUpAdded->Drc ;
                    DFForcingUpAdded->Drc = 0;
                }
            }

            iterf ++;
            if(iterf > 25)
            {
                return;
            }
        }
    }

}

void TWorld::CalculateSafetyFactors(cTMap * _DEM,cTMap * _SoilDepth,
                                   cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                   cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                   cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                   cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                   cTMap * _SafetyFactor,cTMap * _Threshold,
                                   cTMap * _Threshold1,cTMap * _InititationHeight,
                                   cTMap * _Initiated,cTMap * _SFCalibration,
                                   cTMap * _DFForcing, cTMap * _DFForcingUp,
                                   cTMap * _PGA, bool limit_slope, double factor_failmax)
{


    _SafetyFactor->setAllMV();
    _InititationHeight->setAllMV();
    FOR_ROW_COL_MV
    {

        //get the actual slope for this cell

        double SlopeX = 0;
        double SlopeY = 0;
        double Slope = 0;

        //DEM
        double dem = _DEM->Drc;

        double nx = 0;
        double ny = 0;

        double demx = 0;
        double demy = 0;

        double min_surr = dem;

        if(!OUTORMV(r+1,c))
        {
            nx += 1.0;
            demx += _DEM->data[r+1][c] - dem;
            min_surr = std::min(min_surr,_DEM->data[r+1][c]);
        }
        if(!OUTORMV(r-1,c))
        {
            nx += 1.0;
            demx += dem -_DEM->data[r-1][c];
            min_surr = std::min(min_surr,_DEM->data[r-1][c]);
        }

        if(!OUTORMV(r,c+1))
        {
            ny += 1.0;
            demy += _DEM->data[r][c+1] - dem;
            min_surr = std::min(min_surr,_DEM->data[r][c+1]);
        }
        if(!OUTORMV(r,c-1))
        {
            ny += 1.0;
            demy += dem -_DEM->data[r][c-1];
            min_surr = std::min(min_surr,_DEM->data[r][c-1]);
        }

        if(limit_slope)
        {
            FailLimit->Drc = std::max(0.0,dem - min_surr);
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

        Slope = std::max(0.0,_DEM->Drc - min_surr)/_dx;//  fabs(SlopeX) + fabs(SlopeY);

        DFSlope->Drc = Slope;
        double angle = atan(Slope);
        double cosa = cos(angle);
        double sina = sin(angle);
        double tanphi = tan(_InternalFrictionAngle->Drc);

        //std::max(1.0,-60.0 * (_SoilWaterHeight->Drc/_SoilDepth->Drc) + 19.0)
        //+_SoilWaterSuction->Drc
        //calculate safety factor

        //*(1 - std::min(0.0,std::max(1.0,_SoilWaterHeight->Drc/_SoilDepth->Drc)))

        double t1 =
                (_DFForcingUp->Drc + _SoilCohesion->Drc  + _PlantCohesion->Drc )
                    +
                    (
                        -cosa * sina*(_PGA->Drc/9.81)*
                        (
                            _SoilDensity->Drc *(_SoilDepth->Drc - _SoilWaterHeight->Drc) + 1000.0 *(_SoilWaterHeight->Drc)
                        )
                        +
                        cosa * cosa*
                        (

                               (
                                    _SoilDensity->Drc *(_SoilDepth->Drc - _SoilWaterHeight->Drc) + 1000.0 *(_SoilWaterHeight->Drc)
                                )
                         )
                     )*tanphi;

        double t2 = std::max(0.0,_DFForcing->Drc +
                (
                cosa*cosa*(_PGA->Drc/9.81)*(_SoilDensity->Drc * _SoilDepth->Drc + 1000.0 *(_SoilWaterHeight->Drc ))
                    +
                    ( _SoilDensity->Drc * _SoilDepth->Drc + 1000.0 *(_SoilWaterHeight->Drc ))
                *sina*cosa
                ));


        double sf;



        //above 1000 is irrelevant
        if(t2 < 1e-12)
        {

            sf = 100.999;
        }else
        {

            sf = std::max(0.0,std::min(t1/t2,1000.0));
        }


        //compensate for initial calibration
        _SafetyFactor->Drc = _SFCalibration->Drc * sf;

    }

    FOR_ROW_COL_MV
    {

        bool boundary = false;

        bool nb_df = false;
        if(!OUTORMV(r-1,c)){if(_Initiated->data[r-1][c] > 0){nb_df = true;}}else{boundary = true;}
        if(!OUTORMV(r+1,c)){if(_Initiated->data[r+1][c] > 0){nb_df = true;}}else{boundary = true;}
        if(!OUTORMV(r,c-1)){if(_Initiated->data[r][c-1] > 0){nb_df = true;}}else{boundary = true;}
        if(!OUTORMV(r,c+1)){if(_Initiated->data[r][c+1] > 0){nb_df = true;}}else{boundary = true;}

        if(boundary)
        {
            continue;
        }

        if((_SafetyFactor->Drc < _Threshold->Drc || (_SafetyFactor->Drc < _Threshold1->Drc && nb_df) ) && _SoilDepth->Drc > 0.1)
        {

            if(SF_Calibrate_LF)
            {
                if(DFFailureMask->Drc == 0)
                {
                    _InititationHeight->Drc = 0;
                    continue;
                }
            }

            double slope = DFSlope->Drc;
            double angle = atan(slope);
            double cosa = cos(angle);
            double sina = sin(angle);

            //get all relevant vvariables for calculation of stable depth
            double wf = _SoilWaterHeight->Drc/_SoilDepth->Drc;
            double wd = 1000.0;
            double hf = _SoilDepth->Drc - DFSlope->Drc * _dx;

            //threshold safety factor value
            //include sf calibration (compensation factor to make sure that the initial state is stable)
            double sf = _Threshold1->Drc;// / _SFCalibration->Drc;

            double h0 = hf;
            double a = _SoilCohesion->Drc + _PlantCohesion->Drc;
            double cc =  std::tan(_InternalFrictionAngle->Drc) *((1.0 - wf) *_SoilDensity->Drc-wd * wf);
            double b = (wd * wf + (1.0 - wf) *_SoilDensity->Drc);
            double d = 0;//(_DFForcing->Drc )/_SoilDepth->Drc + _PGA->Drc *(wd * wf + (1.0 - wf) *_SoilDensity->Drc)*cosa*cosa;
            double e = 0.0;
            double f = 0.0;
            double g = 0;//(_DFForcingUp->Drc)/_SoilDepth->Drc +_PGA->Drc * (wd * wf + (1.0 - wf) *_SoilDensity->Drc)*cosa*sina;
            double hdx = _dx;

            double soildepthn = _SoilDepth->Drc;
            double soildepths = _SoilDepth->Drc;



            //Use wolfram mathematica to find a analytical solution

            /*sol = Reduce[{sf == (a + (h - h1)*
            g + ((e + h*b)*Cos[ArcTan[(h - h0)/xd]]^2.0 ))/((h - h1)*
            d + ((f + h*c)*
             Sin[ArcTan[(h - h0)/xd]] Cos[ArcTan[(h - h0)/xd]])),
            sf != 0, xd != 0 }, h]*/

            //this is the root to the third power polynomial equation that is real

                        /*h == Root[a h0^2 - 1. g h0^2 h1 + d h0^2 h1 sf + f h0 sf xd + a xd^2 + e xd^2 - 1. g h1 xd^2 + d h1 sf xd^2
            + (-2. a h0 + g h0^2 + 2. g h0 h1 - 1. d h0^2 sf - 2. d h0 h1 sf - 1. f sf xd + c h0 sf xd + b xd^2 + g xd^2 - 1. d sf xd^2) #1
            + (a - 2. g h0 - 1. g h1 + 2. d h0 sf + d h1 sf - 1. c sf xd) #1^2
             + (g - 1. d sf) #1^3 &, 1] */

            if(SwitchSeismic || SwitchUpslopeForcing || SwitchDownslopeForcing)
            {



                //to get the third power root with real value
                //double c1 = 1.0 * a*h0*h0 - g*h0*h0*h0 + d * h0*h0*h0* sf+ f*h0*sf*hdx +a * hdx * hdx + e*hdx*hdx - g*h0*hdx*hdx + d*h0*sf*hdx*hdx;
                //double c2 = -2.0 * a * h0 + g * h0*h0 + 2.0 * g * h0*h0 - 1.0 * d * h0 * h0*sf - 2.0 * d*h0*h0*sf - 1.0 * f*sf*hdx + cc * h0*sf*hdx + b * hdx*hdx + g * hdx* hdx - d * sf*hdx*hdx;
                //double c3 = a - 2.0 * g * h0 - g * h0 + 2.0 * d * h0*sf + d * h0*sf - cc* sf*hdx;
                //double c4 = g - (d *sf);

                /*double c1 = a*h0*h0 + a * hdx*hdx;
                double c2 =-2.0*a*h0 + g * h0 * h0 - 1.0 * d * h0*h0*sf + c * h0*sf*hdx + b *hdx*hdx + g * hdx*hdx - 1.0 * d * sf * hdx*hdx ;
                double c3 = a - 2.0 * g * h0 + 2.0 * d * h0 * sf - 1.0 * cc * sf * hdx;
                double c4 = g - (d *sf);

                double *x = new double[3];
                x[0] = 0.0;
                x[1] = 0.0;
                x[2] = 0.0;

                int nsol = 0;
                if( std::fabs(c4) > 0.0)
                {
                    nsol = SolveP3(x,c3/c4,c2/c4,c1/c4);
                }else if( std::fabs(c3) > 0.0)
                {
                    nsol = SolveP2(x,c2/c3,c1/c3);
                }else if( std::fabs(c2) > 0.0)
                {
                    soildepths = -c1/c2;
                }else
                {
                    soildepths = _SoilDepth->Drc;
                }

                soildepths = (_SoilDepth->Drc);



                double ifzero = soildepths;
                if(nsol == 1)
                {
                    soildepths = x[0] < 0 ? ifzero : std::min(soildepths,x[0]);
                    soildepths = std::max(0.0,std::min(soildepths,_SoilDepth->Drc));

                }else if(nsol == 2)
                {
                    soildepths = x[0] < 0 ? ifzero : std::min(soildepths,x[0]);
                    soildepths = x[1] < 0 ? ifzero : std::min(soildepths,x[1]);
                    soildepths = std::max(0.0,std::min(soildepths,_SoilDepth->Drc));

                }else if(nsol == 3)
                {
                    soildepths = x[0] < 0 ? ifzero : std::min(soildepths,x[0]);
                    soildepths = x[1] < 0 ? ifzero : std::min(soildepths,x[1]);
                    soildepths = x[2] < 0 ? ifzero : std::min(soildepths,x[2]);
                    soildepths = std::max(0.0,std::min(soildepths,_SoilDepth->Drc));
                }

                if(soildepths > _SoilDepth->Drc - 1e-10)
                {
                    soildepths = 0;
                }*/

                /*double t1 = (-2.0*a*h0 - 2.0*g*h0 + b*hdx*hdx + 2.0*d*h0*sf +
                             c*h0*hdx*sf);

                double h1 = (0.5 *(2.0* a*h0 + 2.0*g*h0 - 1.0*b*hdx*hdx - 2.0*d*h0*sf - 1.0*c*h0*hdx*sf - 1.0 *
                                   std::sqrt(t1*t1 - 4.0 *(a + g - 1.0*d*sf - 1.0*c*hdx*sf)*(a*h0*h0 + g*h0*h0 + a*hdx*h0 + g*hdx*h0 - 1.0*d*h0*h0 *sf - 1.0*d*hdx*hdx*sf))))
                        /(a + g - 1.0*d*sf - 1.0* c*hdx*sf);
                double h2 = (0.5 *(2.0* a*h0 + 2.0*g*h0 - 1.0*b*hdx*hdx - 2.0*d*h0*sf - 1.0*c*h0*hdx*sf + 1.0 *
                                   std::sqrt(t1*t1 - 4.0 *(a + g - 1.0*d*sf - 1.0*c*hdx*sf)*(a*h0*h0 + g*h0*h0 + a*hdx*h0 + g*hdx*h0 - 1.0*d*h0*h0 *sf - 1.0*d*hdx*hdx*sf))))
                        /(a + g - 1.0*d*sf - 1.0* c*hdx*sf);;

                soildepths = std::max(0.0,std::min(_SoilDepth->Drc,std::min(_SoilDepth->Drc,std::max(0.0,std::min(h1,h2)))));

                qDebug() << _SoilDepth->Drc << h1 << h2;
                qDebug() << _SoilDepth->Drc << soildepths << _SafetyFactor->Drc;*/


            }

            double h1 = _SoilDepth->Drc;
            double h2 = _SoilDepth->Drc;
            double sc = 0;
            //if(!(SwitchSeismic || SwitchUpslopeForcing || SwitchDownslopeForcing))
            {
                //get all relevant vvariables for calculation of stable depth
                double cif = cos(_InternalFrictionAngle->Drc);
                double sif = sin(_InternalFrictionAngle->Drc);
                d = _dx;
                double dx2 = _dx*_dx;
                wf = _SoilWaterHeight->Drc/_SoilDepth->Drc;
                wd = 1000.0;
                double pp = _PlantPressure->Drc;
                hf = _SoilDepth->Drc - DFSlope->Drc * _dx;
                double pc = _PlantCohesion->Drc;
                double sd = _SoilDensity->Drc;
                double ws = _SoilWaterSuction->Drc;
                sc = _SoilCohesion->Drc - _DFForcing->Drc + _DFForcingUp->Drc -  (_PGA->Drc/9.81) *(wd * wf + (1.0 - wf) *_SoilDensity->Drc)*cosa*sina *std::tan(_InternalFrictionAngle->Drc) -  (_PGA->Drc/9.81) *(wd * wf + (1.0 - wf) *_SoilDensity->Drc)*cosa*cosa;

                //threshold safety factor value
                //include sf calibration (compensation factor to make sure that the initial state is stable)
                sf = _Threshold1->Drc;// / _SFCalibration->Drc;

                //solution for stable depth at which safety factor equals a threshold value
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

                //in the end, it is nothin more than a really complex implementation of the abc-formula
                // x1,2 = (-b +- Sqrt(b^2-4ac))/2a
                //Two solutions for this formula, in our case: 1 positive, 1 negative, so we pick the positive value

                double t2 = sqrt(t21*t21 - 4.0*t22);
                double t3 = (2.0*(pc*cif + sc*cif - d*sd*sf*cif + ws*cif));

                //positive and negative solution
                h1 = (t1 + t2)/t3;
                h2 = (t1 - t2)/t3;

                h1 = std::isnan(h1)? 0.0:h1;
                h2 = std::isnan(h2)? 0.0:h2;

                //qDebug() << h1 << h2 << _DFForcing->Drc << _DFForcingUp->Drc;
            }

            double hnew_select1 = std::max(h1,h2);
            double hnew_select2 = std::min(h1,h2);

            double hnew_select = hnew_select1;

            double hnew = std::max(_SoilDepth->Drc* factor_failmax,std::min(_SoilDepth->Drc,std::min(soildepths,std::min(soildepthn,std::min(_SoilDepth->Drc,std::max(0.0,hnew_select))))));


            double fheight = std::max(0.0,_SoilDepth->Drc - hnew);


            if(limit_slope)
            {
                fheight = std::min(fheight,FailLimit->Drc);
            }

            //final initiation height can not be more than Soildepth
            _InititationHeight->Drc =fheight;

            //minimum initiation height is 0.1 meters depth, otherwise not relevant
            double DF_MinInitiationHeight = 0.0;

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


double TWorld::CalculateSafetyFactorAt(int r, int c)
{
    return CalculateSafetyFactorAt(r,c,DFSlope->Drc,DFSoilDepth,DFSurfaceWaterHeight,DFSoilCohesion,DFSoilInternalFrictionAngle,DFWaterHeight,DFWaterSuction,DFSoilDensity,DFPlantCohesion,DFPlantPressure,DFSFCalibration);

}


double TWorld::CalculateSafetyFactorAt(int r, int c, double slope, cTMap * _SoilDepth,
                                   cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                   cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                   cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                   cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                   cTMap * _SFCalibration)
{


    return CalculateSafetyFactor(slope,_SoilDepth->Drc,_OverlandWater->Drc,
                                 _SoilCohesion->Drc,_InternalFrictionAngle->Drc,
                                 _SoilWaterHeight->Drc,_SoilWaterSuction->Drc,
                                 _SoilDensity->Drc,_PlantCohesion->Drc,
                                 _PlantPressure->Drc,_SFCalibration->Drc);
}

double TWorld::CalculateSafetyFactor(double slope, double _SoilDepth,
                                   double _OverlandWater, double _SoilCohesion,
                                   double _InternalFrictionAngle,double _SoilWaterHeight,
                                   double _SoilWaterSuction, double _SoilDensity,
                                   double _PlantCohesion, double _PlantPressure,
                                    double _SFCalibration)
{

    double angle = atan(slope);
    double cosa = cos(angle);
    double sina = sin(angle);
    double tanphi = tan(_InternalFrictionAngle);

    //calculate safety factor

    double t1 = (_SoilCohesion *(1 - std::min(0.0,std::max(1.0,_SoilWaterHeight/_SoilDepth))) + _PlantCohesion )
                +
                (
                    cosa * cosa*
                    (
                     _PlantPressure +
                           (
                                _SoilDensity *(_SoilDepth - _SoilWaterHeight) + 1000.0 *(_SoilWaterHeight)
                            )
                     )*tanphi
                 );

    double t2 = (
                (_SoilDensity * _SoilDepth + 1000.0 *(_SoilWaterHeight )+ _PlantPressure)
            *sina*cosa
            );

    //above 1000 is irrelevant
    double sf = (t2 > 0 ? std::max(0.0,std::min(_SFCalibration * t1/t2,1000.0)): 1000.0);
    return sf;
}

double TWorld::GetTotalSoilDepth(int r, int c)
{
    return DFSoilDepth->Drc;


}

double TWorld::SolveStableDepthAt(int r, int c)
{
    return SolveStableDepthAt(r,c,DFSlope->Drc,DFSoilDepth,DFSurfaceWaterHeight,DFSoilCohesion,DFSoilInternalFrictionAngle,DFWaterHeight,DFWaterSuction,DFSoilDensity,DFPlantCohesion,DFPlantPressure,DFThreshold,DFThreshold1,DFSFCalibration);


}

double TWorld::SolveStableDepthAt(int r, int c, double slope, cTMap * _SoilDepth,
                                   cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                   cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                   cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                   cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                   cTMap * _Threshold,cTMap * _Threshold1,cTMap * _SFCalibration)
{


    return SolveStableDepth(slope,_SoilDepth->Drc,_OverlandWater->Drc,
                            _SoilCohesion->Drc,_InternalFrictionAngle->Drc,
                            _SoilWaterHeight->Drc,_SoilWaterSuction->Drc,
                            _SoilDensity->Drc,_PlantCohesion->Drc,
                            _PlantPressure->Drc,_Threshold->Drc,
                            _Threshold1->Drc,_SFCalibration->Drc);
}

double TWorld::SolveStableDepth(double slope, double _SoilDepth,
                                   double _OverlandWater, double _SoilCohesion,
                                   double _InternalFrictionAngle,double _SoilWaterHeight,
                                   double _SoilWaterSuction, double _SoilDensity,
                                   double _PlantCohesion, double _PlantPressure,
                                   double _Threshold, double _Threshold1,double _SFCalibration)
{


    //get all relevant vvariables for calculation of stable depth
    double cif = cos(_InternalFrictionAngle);
    double sif = sin(_InternalFrictionAngle);
    double d = _dx;
    double dx2 = _dx*_dx;
    double wf = _SoilWaterHeight/_SoilDepth;
    double wd = 1000.0;
    double pp = _PlantPressure;
    double hf = _SoilDepth - slope * _dx;
    double pc = _PlantCohesion;
    double sd = _SoilDensity;
    double ws = _SoilWaterSuction;
    double sc = _SoilCohesion;

    //threshold safety factor value
    //include sf calibration (compensation factor to make sure that the initial state is stable)
    double sf = _Threshold1;// / _SFCalibration->Drc;


    //solution for stable depth at which safety factor equals a threshold value
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

    //in the end, it is nothin more than a really complex implementation of the abc-formula
    // x1,2 = (-b +- Sqrt(b^2-4ac))/2a
    //Two solutions for this formula, in our case: 1 positive, 1 negative, so we pick the positive value

    double t2 = sqrt(t21*t21 - 4.0*t22);
    double t3 = (2.0*(pc*cif + sc*cif - d*sd*sf*cif + ws*cif));

    //positive and negative solution
    double h1 = (t1 + t2)/t3;
    double h2 = (t1 - t2)/t3;

    //final stable depth
    return std::max(0.0,std::max(h1,h2));
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

            double h = !std::isnan(DFInitiationHeight->Drc)? DFInitiationHeight->Drc : 0.0;

            //change DEM (flow DEM is altered later through DEMChange)
            DEMOriginal->Drc -= h;
            //DEM->Drc -= h;
            DEMChange->Drc -= h;

            double area = _dx*_dx;
            if(SwitchEntrainment && INCLUDE_DEPOSITS_STABILITY)
            {
                  double depthdep = SoilRockMaterial->Drc /area;
                  double depthdepfail = std::min(h,depthdep);
                  h = std::max(0.0, h - depthdepfail);
                  double sdepvol = depthdepfail*area;
                  UF2D_f->Drc +=sdepvol;
                  UF2D_s->Drc +=sdepvol* 0.4;
                  SoilRockMaterial->Drc =std::max(0.0,SoilRockMaterial->Drc-sdepvol);
                  UF2D_su->Drc = (sdepvol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_su->Drc)/(sdepvol + UF2D_s->Drc) : 0.0;
                  UF2D_sv->Drc = (sdepvol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_sv->Drc)/(sdepvol + UF2D_s->Drc) : 0.0;
                  UF2D_d->Drc = (sdepvol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_d->Drc + DFSoilDensity->Drc * sdepvol)/(sdepvol + UF2D_s->Drc) : UF2D_d->Drc;
                  UF2D_ifa->Drc = (sdepvol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_ifa->Drc + DFSoilInternalFrictionAngle->Drc * sdepvol)/(sdepvol + UF2D_s->Drc) : UF2D_ifa->Drc;
                  UF2D_rocksize->Drc = (sdepvol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_rocksize->Drc + DFSoilRockSize->Drc * sdepvol)/(sdepvol + UF2D_s->Drc) : UF2D_rocksize->Drc;
            }


            //fluid volume is added

            double fh = 0;
            if(this->InfilMethod != INFIL_NONE)
            {
                fh = std::min(SoilDepth1->Drc,h) * (L1->Drc * (ThetaS1->Drc - ThetaI1->Drc) + ThetaI1->Drc * SoilDepth1->Drc)/SoilDepth1->Drc ;
                if(SwitchTwoLayer)
                {
                     fh += std::min(SoilDepth2->Drc,std::max(h-SoilDepth1->Drc,0.0)) * (L2->Drc * (ThetaS2->Drc - ThetaI2->Drc) + ThetaI2->Drc * SoilDepth2->Drc)/SoilDepth2->Drc;
                     L2->Drc = std::max(0.0,L2->Drc - std::max(0.0,h-SoilDepth1->Drc));
                }
                L1->Drc = std::max(0.0,L1->Drc - h);
            }

            double fvol = fh * _dx * _dx;
            fvol = !std::isnan(fvol)? fvol : 0.0;

            UF2D_f->Drc += fvol;

            UF_InitializedF += fvol;

            //solid phase volume is added
            double svol = std::min(SoilDepth1->Drc,h) * (1.0-ThetaS1->Drc) * _dx * _dx;
            if(SwitchTwoLayer)
            {
                svol += std::min(SoilDepth2->Drc,std::max(h-SoilDepth1->Drc,0.0)) * (1.0-ThetaS2->Drc) * _dx * _dx;
            }

            svol = !std::isnan(svol)? svol : 0.0;
            //and flow properties updated
            UF2D_su->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_su->Drc)/(svol + UF2D_s->Drc) : 0.0;
            UF2D_sv->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_sv->Drc)/(svol + UF2D_s->Drc) : 0.0;
            UF2D_d->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_d->Drc + DFSoilDensity->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_d->Drc;
            UF2D_ifa->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_ifa->Drc + DFSoilInternalFrictionAngle->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_ifa->Drc;
            UF2D_rocksize->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_rocksize->Drc + DFSoilRockSize->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_rocksize->Drc;

            UF2D_s->Drc += svol;

            //suspended sediment is added
            /*double SedMass = (1.0-DFSoilRockFraction->Drc) * std::min(SoilDepth1->Drc,h) * DFSoilDensity->Drc * (1.0-ThetaS1->Drc) * _dx * _dx;
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
            }*/

            //Depth of infiltration layers decreased
            double sd1 = SoilDepth1->Drc;
            SoilDepth1->Drc = std::max(0.0,SoilDepth1->Drc - h);
            if(SwitchTwoLayer)
            {
                 SoilDepth2->Drc = std::max(0.0,(SoilDepth2->Drc) - std::max(h - sd1,0.0));
            }

        }

        if(SwitchBedrock)
        {
            if(DFInitiationHeight2->Drc > UF_VERY_SMALL)
            {
                n_init ++;
                DFTotalInitiationHeight->Drc += DFInitiationHeight2->Drc;

                double h = DFInitiationHeight2->Drc;

                //change DEM (flow DEM is altered later through DEMChange)
                DEMOriginal->Drc -= h;
                //DEM->Drc -= h;
                DEMChange->Drc -= h;

                //solid phase volume is added
                double svol = DFInitiationHeight2->Drc * _dx * _dx;

                //and flow properties updated
                UF2D_su->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_su->Drc)/(svol + UF2D_s->Drc) : 0.0;
                UF2D_sv->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_sv->Drc)/(svol + UF2D_s->Drc) : 0.0;
                UF2D_d->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_d->Drc + DFSoilDensity->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_d->Drc;
                UF2D_ifa->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_ifa->Drc + DFSoilInternalFrictionAngle->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_ifa->Drc;
                UF2D_rocksize->Drc = (svol + UF2D_s->Drc)> 0? (UF2D_s->Drc *UF2D_rocksize->Drc + DFSoilRockSize->Drc * svol)/(svol + UF2D_s->Drc) : UF2D_rocksize->Drc;

                UF2D_s->Drc += svol;


            }

        }

    }

}
