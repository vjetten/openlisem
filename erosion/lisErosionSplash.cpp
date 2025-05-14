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


#include "model.h"


//---------------------------------------------------------------------------
void TWorld::cell_SplashDetachment()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L  {

        double _WH = FloodDomain->Drc == 0 ? WH->Drc : hmx->Drc;

        DETSplash->Drc = 0;

        if(_WH > HMIN && SplashStrength->Drc >= 0)
        {
            double DetDT1 = 0, DetDT2 = 0, DetLD1, DetLD2;
            double g_to_kg = 0.001;
            double Lc = SwitchLitter ? Litter->Drc : 0.0;
            double Cv = Cover->Drc;
            double strength = SplashStrength->Drc;
            double Int = Rain->Drc * 3600/_dt * 1000; // intensity in mm/h, Rain is in m
            double KE_DT = 0.0;
            double DETSplash_;

            switch (KEequationType)
            {
                case KE_EXPFUNCTION: KE_DT = KEParamater_a1*(1-(KEParamater_b1*exp(-KEParamater_c1*Int))); break;
                case KE_LOGFUNCTION: KE_DT = (Int > 1 ? KEParamater_a2 + KEParamater_b2*log10(Int) : 0); break;
                case KE_POWERFUNCTION: KE_DT = KEParamater_a3*pow(Int, KEParamater_b3); break;
                // kin energy in J/m2/mm
            }
            //VJ 110706  KE equations

            double directrain = (1-Lc) * (1-Cv)*Rainc->Drc * 1000;
            // Added litter also to directrain, assume it covers the entire cell, not only under the plant
            // rainfall between plants in mm

            double KE_LD = std::max(15.3*sqrt(PlantHeight->Drc)-5.87, 0.0);
            // kin energy in J/m2/mm
            double throughfall = (1-Lc) * Cv * LeafDrain->Drc * 1000;
            // leaf drip in mm, is calculated as plant leaf drip in interception function so mult cover

            double WH0 = exp(-1.48*_WH*1000);
            // water buffer effect on surface, WH in mm in this empirical equation from Torri ?

            // fraction ponded area
            double FPA = 1.0;
            if (RR->Drc > 0.1)
                FPA =  1-exp(-1.875*(_WH/(0.01*RR->Drc)));

            // Between plants, directrain is already with 1-cover
            DetDT1 = g_to_kg * FPA*strength*KE_DT*WH0 * directrain;
            //ponded areas, kg/m2/mm * mm = kg/m2
            DetDT2 = _WH > 0 ? g_to_kg * (1-FPA)*strength*KE_DT * directrain * SplashDelivery: 0.0;
            //dry areas, kg/m2/mm * mm = kg/m2


            if (SwitchKETimebased)
            {
                if (directrain > 0)
                {
                    DetDT1 = g_to_kg * FPA*strength*KE_DT*WH0 * _dt/3600;
                    //ponded areas, kg/m2/sec * sec = kg/m2
                    DetDT2 = g_to_kg * (1-FPA)*strength*KE_DT * _dt/3600 * SplashDelivery;
                    //dry areas, kg/m2/sec * sec = kg/m2
                }
            }
            //based on work by Juan Sanchez

            // Under plants, throughfall is already with cover
            DetLD1 = g_to_kg * FPA*(strength*KE_LD)*WH0 * throughfall;
            //ponded areas, kg/m2/mm * mm = kg/m2
            DetLD2 = g_to_kg * (1-FPA)*(strength*KE_LD) * throughfall * SplashDelivery;
            //dry areas, kg/m2/mm * mm = kg/m2

            DETSplash_ = DetLD1 + DetLD2 + DetDT1 + DetDT2;
            // Total splash kg/m2

            // Deal with all exceptions:

            DETSplash_ *= (SoilWidthDX->Drc*DX->Drc);
            // kg/cell, only splash over soilwidth, not roads/hardsurfaces and channels ! houses re not in soilwisth, need to be done here
            // FROM KG/M2 TO KG/CELL

            DETSplash_ = (1-StoneFraction->Drc) * DETSplash_;
            // no splash on stone surfaces

            if (SwitchGrassStrip)
                DETSplash_ = (1-GrassFraction->Drc) * DETSplash_;

            //      if(SwitchSedtrap)
            //          DETSplash->Drc = (1-SedimentFilter->Drc) * DETSplash->Drc;
            // assume sedtrap can have splash

            if (SwitchHouses)
                DETSplash_ = (1-HouseCover->Drc)*DETSplash_;
            //is already contained in soilwidth
            // no splash from house roofs

            if (SwitchSnowmelt)
                DETSplash_ = (1-Snowcover->Drc)*DETSplash_;
            // no splash on snow deck

            if(SwitchKinematic2D == K2D_METHOD_DYN) {
                SSFlood->Drc += DETSplash_;
                SSCFlood->Drc = MaxConcentration(WaterVolall->Drc, SSFlood->Drc);
            } else {
                if (FloodDomain->Drc > 0) {
                    SSFlood->Drc += DETSplash_;
                    SSCFlood->Drc = MaxConcentration(CHAdjDX->Drc * hmx->Drc, SSFlood->Drc);

                } else {
                    Sed->Drc += DETSplash_;
                    Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
                }
            }

            DETSplash->Drc = DETSplash_;
            // IN KG/CELL
        }
    }}
}
