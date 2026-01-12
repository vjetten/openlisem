/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file lisErosion.cpp
  \brief Flow and splash detachment functions for slopes and channels

functions: \n
- double TWorld::MaxConcentration(double watvol, double sedvol)\n
- void TWorld::SplashDetachment(void)\n
- double TWorld::GetTotalDW(int r, int c,QList<cTMap *> *M)\n
- double TWorld::GetDp(int r, int c,double p)\n
- double TWorld::GetDpMat(int r, int c,double p,QList<cTMap *> *M)\n
- double TWorld::GetMpMat(int r, int c,double p,QList<cTMap *> *M, QList<double> *V)\n
- double TWorld::GetSV(double d)\n
- void TWorld::SedimentSetMaterialDistribution(int r,int c)\n
- double TWorld::DetachMaterial(int r,int c, int d,bool channel, bool flood,bool bl,double detachment)\n
- void TWorld::FlowDetachment(void)\n
- void TWorld::ChannelFlowDetachment(int r, int c)\n
- void TWorld::RiverSedimentMaxC(int r, int c)\n
- void TWorld::RiverSedimentDiffusion(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)\n
- void TWorld::RiverSedimentLayerDepth(int r , int c)\n

*/

//#include <algorithm>
//#include "operation.h"
#include "model.h"

//---------------------------------------------------------------------------
//OBSOLETE!!!

/**
 * @fn void TWorld::FlowDetachment(void)
 * @brief Calculates flow detachment for overland flow in entire catchment
 *
 * This function uses the function for overland flow transport capacity to
 * calculate the potential detachment/deposition based on the settling velocity of the sediment.
 * When potential detachment is found, the fuction for taking soil
 * from the soil layer is used to find actual detachment.
 * When deposition is found, this sediment is added to the deposited soil layer.
 *
 * @see TWorld:OFTC
 * @see TWorld:GetSV
 * @see TWorld:DetachMaterial
 *
 */

// Overland flow erosion for 1D flow only
void TWorld::cell_FlowDetachment()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L  {
        double erosionwh = WHrunoff->Drc;
        double erosionwv = WHrunoff->Drc*CHAdjDX->Drc;

        //transport capacity
        DETFlow->Drc = 0;
        DEP->Drc = 0;
        TC->Drc = calcTCSuspended(r,c,-1, FS_SS_Method, WHrunoff->Drc, FlowWidth->Drc, V->Drc, 2);
        // trasnport capacity. 2 = kin wave. 1 = 2d flow and 0 is river

        if (erosionwh < HMIN) {
            if(DO_SEDDEP == 1) {
                DEP->Drc += -Sed->Drc;
                Sed->Drc = 0;
                Conc->Drc = 0;
                TC->Drc = 0;
            }
        } else {
            double maxTC = 0;
            double minTC = 0;

            double deposition = 0;
            double detachment = 0;
            double TransportFactor = 0;

            maxTC = qMax(TC->Drc - Conc->Drc,0.0);
            // positive difference: TC defi  cit becomes detachment (positive)
            minTC = qMin(TC->Drc - Conc->Drc,0.0);
            // negative difference: TC surplus becomes deposition (negative)
            // unit kg/m3

            //### deposition ###
            if (minTC < 0) {

                //TransportFactor = (1-exp(-_dt*SettlingVelocitySS->Drc/erosionwh)) * erosionwv;
                TransportFactor = _dt*SettlingVelocitySS->Drc * DX->Drc * ChannelAdj->Drc;
                // in m3
                // deposition can occur on roads and on soil (so use flowwidth)
                deposition = minTC * TransportFactor;
                // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

                deposition = qMax(deposition, -Sed->Drc);

                //if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                //    deposition = 0;
                // kin wave, there is no boundary

                if (SwitchSedtrap && SedMaxVolume->Drc == 0 && N->Drc == SedTrapN) {
                    N->Drc = Norg->Drc;
                }
                // mannings N becomes normal when sedtrap is full

                if (SwitchSedtrap && SedMaxVolume->Drc > 0) {
                    if (Sed->Drc > 0) {
                        double depvol = Sed->Drc/BulkDens; // m3
                        if (SedMaxVolume->Drc < depvol)
                            depvol = SedMaxVolume->Drc;
                        if (SedMaxVolume->Drc > 0){
                            deposition = -depvol*BulkDens;
                            maxTC = 0;
                        }
                        SedMaxVolume->Drc = SedMaxVolume->Drc - depvol;
                        SedimentFilter->Drc += depvol*BulkDens;
                    }
                }

                if(SwitchGridRetention) {
                    if (Sed->Drc > 0) {
                        double depvol = Sed->Drc/BulkDens; // sed in m3
                        if (GridRetention->Drc < depvol)
                            depvol = GridRetention->Drc;
                        if (GridRetention->Drc > 0){
                            deposition = -depvol*BulkDens;  // deposition is all that goes into trench
                            maxTC = 0;
                        }
                        GridRetention->Drc = GridRetention->Drc - depvol;
                    }
                }

            } else
              //### detachment ###
              if (maxTC > 0 && CohesionSoil->Drc > 0) {

                TransportFactor = _dt*SettlingVelocitySS->Drc * DX->Drc * SoilWidthDX->Drc;
                // soilwidth is erodible surface
                // TransportFactor = qMin(TransportFactor, Q->Drc*_dt);
                // detachment can only come from soil, not roads (so do not use flowwidth)
                // units s * m/s * m * m = m3

                detachment = Y->Drc * maxTC * TransportFactor;//qMin(TransportFactor, erosionwv);
                // unit = kg/m3 * m3 = kg (/cell)

                // exceptions
                if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                    detachment = 0;
                // VJ 190325 prevent any activity on the boundary!

                if (GrassFraction->Drc > 0)
                    detachment = (1-GrassFraction->Drc) * detachment;
                // no flow detachment on grass strips

                // Detachment edxceptions:
                detachment = (1-StoneFraction->Drc) * detachment;
                // no flow detachment on stony surfaces

                if (SwitchHouses)
                    detachment = (1-HouseCover->Drc)*detachment;
                // no flow det from house roofs

                if (SwitchSnowmelt)
                    detachment = (1-Snowcover->Drc) * detachment;

                detachment *= qBound(0.0,1.0 - (RoadWidthHSDX->Drc/_dx),1.0);
                // no flow detachment on hard surfaces, map is 0 is not selected

                if (SwitchSedtrap && SedMaxVolume->Drc >= 0)
                    detachment = 0;

                if (SwitchGridRetention && GridRetention->Drc >= 0)
                    detachment = 0;

                if(Sed->Drc+detachment > MAXCONC * erosionwv)
                    detachment = MAXCONC * erosionwv - Sed->Drc;
                // not more detachment then is possible to keep below diff(max concetrantion-sediment inf low)

              } // minv > 0

            //### sediment balance
            // add to sediment in flow (IN KG/CELL)
            Sed->Drc += detachment;
            Sed->Drc += deposition;
            Sed->Drc = qMax(0.0, Sed->Drc);
            DETFlow->Drc += detachment;
            DEP->Drc += deposition;
            Conc->Drc = MaxConcentration(erosionwv, Sed->Drc);

        }
    }}
}


// experimental nor used
void TWorld::cell_FlowDetachmentContinuous()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L  {
        double erosionwh = WHrunoff->Drc;
        double erosionwv = WHrunoff->Drc*CHAdjDX->Drc;

        //transport capacity
        DETFlow->Drc = 0;
        DEP->Drc = 0;
        TC->Drc = calcTCSuspended(r,c,-1, FS_SS_Method, WHrunoff->Drc, FlowWidth->Drc, V->Drc, 2);
        // trasnport capacity. 2 = kin wave. 1 = 2d flow and 0 is river

        if (erosionwh < HMIN) {
            if(DO_SEDDEP == 1) {
                DEP->Drc += -Sed->Drc;
                Sed->Drc = 0;
                Conc->Drc = 0;
                TC->Drc = 0;
            }
        } else {
            double deposition = 0;
            double detachment = 0;

            //### deposition ###
            deposition = _dt*qMin(1.0, SettlingVelocitySS->Drc/WH->Drc) * -Sed->Drc;
            // fraction of sediment always depostits
            deposition = qMax(deposition, -Sed->Drc);

            //if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
            //   deposition = 0;
            // prevent any activity on the boundary!

            if (SwitchSedtrap && SedMaxVolume->Drc == 0 && N->Drc == SedTrapN) {
                N->Drc = Norg->Drc;
            }
            // mannings N becomes normal when sedtrap is full

            if (SwitchSedtrap && SedMaxVolume->Drc > 0)
            {
                if (Sed->Drc > 0) {
                    double depvol = Sed->Drc * 1.0/BulkDens; // m3
                    if (SedMaxVolume->Drc < depvol)
                        depvol = SedMaxVolume->Drc;
                    if (SedMaxVolume->Drc > 0){
                        deposition = -depvol*BulkDens;
                    }
                }
            }

            if(SwitchGridRetention) {
                if (Sed->Drc > 0) {
                    double depvol = Sed->Drc/BulkDens; // sed in m3
                    if (GridRetention->Drc < depvol)
                        depvol = GridRetention->Drc;
                    if (GridRetention->Drc > 0){
                        deposition = -depvol*BulkDens;  // deposition is all that goes into trench

                    }
                    GridRetention->Drc = GridRetention->Drc - depvol;
                }
            }

            //### detachment ###
            if (CohesionSoil->Drc > 0) {
                double sed = qMax(0.0, Sed->Drc + deposition);
                double conc = MaxConcentration(WaterVolall->Drc, sed);
                detachment = Y->Drc * qMax(0.0, TC->Drc - conc) * _dt*SettlingVelocitySS->Drc * DX->Drc * SoilWidthDX->Drc;
                // unit = kg/m3 * m3 = kg (/cell)

                if (GrassFraction->Drc > 0)
                    detachment = (1-GrassFraction->Drc) * detachment;
                // no flow detachment on grass strips

                // Detachment edxceptions:
                detachment = (1-StoneFraction->Drc) * detachment;
                // no flow detachment on stony surfaces

                if (SwitchHouses)
                    detachment = (1-HouseCover->Drc)*detachment;
                // no flow det from house roofs

                if (SwitchSnowmelt)
                    detachment = (1-Snowcover->Drc) * detachment;

                detachment *= qBound(0.0,1.0 - (RoadWidthHSDX->Drc/_dx),1.0);
                // no flow detachment on hard surfaces, map is 0 is not selected

                if (SwitchSedtrap && SedMaxVolume->Drc > 0)
                    detachment = 0;

                if (SwitchGridRetention && GridRetention->Drc > 0)
                    detachment = 0;

                if(Sed->Drc+detachment > MAXCONC * erosionwv)
                    detachment = MAXCONC * erosionwv - Sed->Drc;
                // not more detachment then is possible to keep below diff(max concetrantion-sediment inf low)
            }

            //### sediment balance
            // add to sediment in flow (IN KG/CELL)
            Sed->Drc += detachment;
            Sed->Drc += deposition;
            Sed->Drc = qMax(0.0, Sed->Drc);
            DETFlow->Drc += detachment;
            DEP->Drc += deposition;
            Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);

        }
    }}
}


