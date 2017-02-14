
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
    }

    //calibration values for slope failure aspect
    st_scCalibration = getvaluedouble("Soil Cohesion Calibration");
    st_sifaCalibration = getvaluedouble("Soil Internal Friction Angle Calibration");
    st_sdCalibration = getvaluedouble("Soil Depth Calibration");
    st_csdCalibration = getvaluedouble("Create Stable Initial Safety Factor");
    st_csdsfCalibration = getvaluedouble("Minimum Safety Factor Calibration");

    //get all the nececcary input from the OpenLISEM model data, and put in in temporary maps
    //this way we can adapt stuff without breaking the rest of the model
    FOR_ROW_COL_MV
    {
        double area = DX->Drc*_dx;
        DFSFIterations->Drc = -1;
        DEMIterate->Drc = DEMOriginal->Drc;
        DFSoilDepth->Drc = SoilDepth1->Drc;
        if(SwitchTwoLayer)
        {
              DFSoilDepth->Drc += SoilDepth2->Drc;
        }
        double soildepth0 = 0;
        double theta0 = 0;
        if(SwitchEntrainment)
        {
              theta0 = (SoilRockMaterial->Drc)? SoilRockWater->Drc/SoilRockMaterial->Drc: 0.0;
              soildepth0 = SoilRockMaterial->Drc /area;
              DFSoilDepth->Drc += soildepth0;
        }
        DFSoilDepth->Drc *= st_sdCalibration;

        DFSurfaceWaterHeight->Drc = UF2D_f->Drc/(_dx*_dx);
        DFSoilCohesion->Drc = DFSoilDepth->Drc > 0? ((Cohesion->Drc * st_scCalibration * (DFSoilDepth->Drc-soildepth0)) + (10.0 * soildepth0))/(DFSoilDepth->Drc) : 0.0;
        if(this->InfilMethod != INFIL_NONE)
        {
            DFWaterHeight->Drc = L1->Drc * (ThetaS1->Drc - ThetaI1->Drc) + ThetaI1->Drc * SoilDepth1->Drc;
            if(SwitchTwoLayer)
            {
                 DFWaterHeight->Drc += L2->Drc * (ThetaS2->Drc - ThetaI2->Drc) + ThetaI2->Drc * SoilDepth2->Drc;
            }
            if(SwitchEntrainment)
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
            if(SwitchEntrainment)
            {
                ts = DFSoilDepth->Drc > 0? ((ts * (DFSoilDepth->Drc-soildepth0)) + (UF_FLOWCOMPACTION_DEPOSITIONPOROSITY * soildepth0))/(DFSoilDepth->Drc) : 0.0;
            }
        }
        saturation = saturation / ts;

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
            DFWaterSuction->Drc = DFWaterSuction->Drc * std::min(1.0,(3.0 - 3.0 *saturation));
        }

        DFPlantCohesion->Drc = RootCohesion->Drc;
        DFPlantPressure->Drc = 0.0;
        DFUnstable->Drc = 0;
        DFThreshold->Drc = SFMIN;
        DFThreshold1->Drc = SFMAX;

        if(SF_Calibrate_First)
        {
            DFSFCalibration->Drc = 1.0;
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

            tma->Drc = DFUnstable->Drc;
            tmb->Drc = 0;
        }

        CalculateSafetyFactors(DEMIterate,DFSoilDepth,DFSurfaceWaterHeight,DFSoilCohesion,DFSoilInternalFrictionAngle,DFWaterHeight,DFWaterSuction,DFSoilDensity,DFPlantCohesion,DFPlantPressure,DFSafetyFactor,DFThreshold,DFThreshold1, tmb,DFInitiationHeight,DFSFCalibration);

        //CALIBRATE INITIAL STABILITY
        //if selected by user, all cells are forced to be at leasat stable
        //Safety factor of cells that were unstable increases up to the threshold plus a certain margin
        //stable cells get a calibration factor of 1, so nothing happens.

        if(SF_Calibrate_Initial)
        {
            if(SF_Calibrate_First)
            {
                SF_Calibrate_First = false;

                FOR_ROW_COL_MV
                {

                    if(DFInitiationHeight->Drc > 0 || DFSafetyFactor->Drc < DFThreshold->Drc * (1+SF_Calibrate_Margin))
                    {

                        DFSFCalibration->Drc = (DFThreshold->Drc * (1.0+SF_Calibrate_Margin *(DFSafetyFactor->Drc/(DFThreshold->Drc*(1.0+SF_Calibrate_Margin))) ))/DFSafetyFactor->Drc ;
                        tmb->Drc = 0;
                        DFSafetyFactor->Drc = DFSafetyFactor->Drc * DFSFCalibration->Drc;

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

        //if slope failure is off anyway, we can get out of the loop anyway!
        if(!SwitchSlopeFailure)
        {
            break;
        }
        iter ++;

    }

    //store safety factor for display
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

void TWorld::CalculateSafetyFactors(cTMap * _DEM,cTMap * _SoilDepth,
                                   cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                   cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                   cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                   cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                   cTMap * _SafetyFactor,cTMap * _Threshold,cTMap * _Threshold1,cTMap * _InititationHeight,cTMap * _Initiated,cTMap * _SFCalibration)
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
            ny += 1.0;
            demy += _DEM->data[r][c+1] - dem;
        }
        if(!OUTORMV(r,c-1))
        {
            ny += 1.0;
            demy += dem -_DEM->data[r][c-1];
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

        //std::max(1.0,-60.0 * (_SoilWaterHeight->Drc/_SoilDepth->Drc) + 19.0)
        //+_SoilWaterSuction->Drc
        //calculate safety factor
        double t1 = (_SoilCohesion->Drc *(1 - std::min(0.0,std::max(1.0,_SoilWaterHeight->Drc/_SoilDepth->Drc))) + _PlantCohesion->Drc )
                    +
                    (
                        cosa * cosa*
                        (
                         _PlantPressure->Drc +
                               (
                                    _SoilDensity->Drc *(_SoilDepth->Drc - _SoilWaterHeight->Drc) + 1000.0 *(_SoilWaterHeight->Drc)
                                )
                         )*tanphi
                     );

        double t2 = (
                    (_SoilDensity->Drc * _SoilDepth->Drc + 1000.0 *(_SoilWaterHeight->Drc )+ _PlantPressure->Drc)
                *sina*cosa
                );
        /*double t1 = (_SoilCohesion->Drc + _PlantCohesion->Drc )
                            +
                            (
                                cosa * cosa*
                                (
                                 _PlantPressure->Drc +
                                       (
                                            _SoilDensity->Drc *(_SoilDepth->Drc) + 1000.0 *(_SoilWaterHeight->Drc)
                                        )
                                 )*tanphi
                             );

                double t2 = (
                            (_SoilDensity->Drc * _SoilDepth->Drc + 1000.0 *(_SoilWaterHeight->Drc )+ _PlantPressure->Drc)
                        *sina*cosa
                        );*/

        //above 1000 is irrelevant
        double sf = (t2 > 0 ? std::max(0.0,std::min(t1/t2,1000.0)): 1000.0);

        //compensate for initial calibration
        _SafetyFactor->Drc = _SFCalibration->Drc * sf;

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
            if(SF_Calibrate_LF)
            {
                if(DFFailureMask->Drc == 0)
                {
                    _InititationHeight->Drc = 0;
                    continue;
                }
            }

            //get all relevant vvariables for calculation of stable depth
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

            //threshold safety factor value
            //include sf calibration (compensation factor to make sure that the initial state is stable)
            double sf = _Threshold1->Drc;// / _SFCalibration->Drc;




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

            //final initiation height can not be more than Soildepth
            _InititationHeight->Drc = std::max(0.0,_SoilDepth->Drc - std::max(h1,h2));

            //minimum initiation height is 0.1 meters depth, otherwise not relevant
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

            double h = DFInitiationHeight->Drc;

            //change DEM (flow DEM is altered later through DEMChange)
            DEMOriginal->Drc -= h;
            //DEM->Drc -= h;
            DEMChange->Drc -= h;

            //fluid volume is added

            double fh = 0;
            if(this->InfilMethod != INFIL_NONE)
            {
                fh = std::min(SoilDepth1->Drc,h) * (L1->Drc * (ThetaS1->Drc - ThetaI1->Drc) + ThetaI1->Drc * SoilDepth1->Drc)/SoilDepth1->Drc ;
                if(SwitchTwoLayer)
                {
                     fh += std::min(SoilDepth2->Drc,std::max(h-SoilDepth1->Drc,0.0)) * (L2->Drc * (ThetaS2->Drc - ThetaI2->Drc) + ThetaI2->Drc * SoilDepth2->Drc)/SoilDepth2->Drc;
                     L2->Drc = std::max(0.0,L2->Drc - std::max(0.0,h-L1->Drc));
                }
                L1->Drc = std::max(0.0,L1->Drc - h);
            }

            double fvol = fh * _dx * _dx;

            UF2D_f->Drc += fvol;

            //solid phase volume is added
            double svol = std::min(SoilDepth1->Drc,h) * (1.0-ThetaS1->Drc) * _dx * _dx;
            if(SwitchTwoLayer)
            {
                svol += std::min(SoilDepth2->Drc,std::max(h-SoilDepth1->Drc,0.0)) * (1.0-ThetaS2->Drc) * _dx * _dx;
            }

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

    }
}
