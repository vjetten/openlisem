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

#include "lisemqt.h"
#include "global.h"
#include "model.h"

//---------------------------------------------------------------------------
void TWorld::InfilSwatre()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        // profile 0 is for impermeable surfaces
        if (ProfileID->Drc <= 0 || fractionImperm->Drc > 0.999) {
            InfilVol->Drc = 0;
            continue;
        }

        double tilevol = 0;

        double WHorig;
        if (FloodDomain->Drc == 0)
            WHorig = WH->Drc;
        else
            WHorig = hmx->Drc;

        SwatreSoilModel->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
        SwatreSoilModel->pixel[i_].tiledrain = 0;

        ComputeForPixel(i_, SwatreSoilModel);

        double WHN = SwatreSoilModel->pixel[i_].wh*0.01;

        Perc->Drc= SwatreSoilModel->pixel[i_].percolation*0.01;
        if (SwitchIncludeTile)
            tilevol = SwatreSoilModel->pixel[i_].tiledrain;  // is already in m3

        //TODO test infil swatre for crusts and compaction
        if (SwitchInfilCrust) {
            if (SwitchDynamicCrusting && ProfileIDCrust->Drc > 0) {
                CrustFraction->Drc = std::min(1.0, CrustFraction0->Drc + (1.0-exp(-0.2*std::max(0.0, RainCumCrust->Drc*1000))));
            }

            if (ProfileIDCrust->Drc > 0 && CrustFraction->Drc > 0) {
                SwatreSoilModelCrust->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
                SwatreSoilModelCrust->pixel[i_].tiledrain = 0;

                ComputeForPixel(i_, SwatreSoilModelCrust);

                double WHcrust = SwatreSoilModel->pixel[i_].wh*0.01;
                WHN = WHcrust*CrustFraction->Drc + WHN*(1-CrustFraction->Drc);
                // weighed average

                if (SwitchIncludeTile) {
                    tilevol = CrustFraction->Drc*SwatreSoilModelCrust->pixel[i_].tiledrain + tilevol*(1-CrustFraction->Drc);
                }
            }
        }

        if (SwitchInfilCompact) {
            if (ProfileIDCompact->Drc > 0 &&  CompactFraction->Drc > 0) {

                SwatreSoilModelCompact->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
                SwatreSoilModelCompact->pixel[i_].tiledrain = 0;

                ComputeForPixel(i_, SwatreSoilModelCompact);

                double WHcompact = SwatreSoilModelCompact->pixel[i_].wh*0.01;
                WHN = WHcompact*CompactFraction->Drc + WHN*(1-CompactFraction->Drc);
                // weighted average

                if (SwitchIncludeTile) {
                    tilevol = CompactFraction->Drc*SwatreSoilModelCompact->pixel[i_].tiledrain + tilevol*(1-CompactFraction->Drc);
                }
            }
        }

        if (SwitchGrassStrip) {
            if (ProfileIDGrass->Drc > 0 &&  GrassFraction->Drc > 0) {
                SwatreSoilModelGrass->pixel[i_].wh = WHorig*100;    // WH is in m, convert to cm
                SwatreSoilModelGrass->pixel[i_].tiledrain = 0;

                ComputeForPixel(i_, SwatreSoilModelGrass);

                double WHgrass = SwatreSoilModelCompact->pixel[i_].wh*0.01;
                WHN = WHgrass*GrassFraction->Drc + WHN*(1-GrassFraction->Drc);

                if (SwitchIncludeTile) {
                    tilevol = GrassFraction->Drc*SwatreSoilModelGrass->pixel[i_].tiledrain + tilevol*(1-GrassFraction->Drc);
                }
            }
        }

        if (SwitchIncludeTile)
            TileWaterVolSoil->Drc = tilevol;

        if (FloodDomain->Drc == 0)
            WH->Drc = WHN;
        else
            hmx->Drc = WHN;
        hmxWH->Drc = hmx->Drc + WH->Drc;
        WaterVolall->Drc = hmxWH->Drc*CHAdjDX->Drc;

        InfilVol->Drc = (WHorig - WHN) * FlowWidth->Drc * DX->Drc;
        // use flowwidth because impermeable is done separately

    }}

    //find depth wetting front, estimated at depth where h is initial value, very crude
    Fill(*Lwmm,0);
    for (int j = 0; j < SwatreSoilModel->pixel[0].profile->zone->nrNodes; j++) {
        cTMap *map = inith->at(j);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (j > 0 && SwatreSoilModel->pixel[i_].h[j] > map->Drc+1.0) {
                double l = SwatreSoilModel->pixel[i_].profile->zone->endComp[j-1]*10; // in mm
                double l1 = SwatreSoilModel->pixel[i_].profile->zone->endComp[j]*10; // in mm
                Lwmm->Drc = 0.5*(l+l1);
            }
        }}
    }

    // dump a map with h at every node
    if(SwitchDumphead) {
        for (int i = 0; i < SwatreSoilModel->pixel[0].profile->zone->nrNodes; i++) {

            QString dig = QString("%1").arg(i+1, 3, 10, QLatin1Char('0'));
            QString hname = QString("head0000.") + dig;
            QString tname = QString("theta000.") + dig;

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                hSwatre->Drc = SwatreSoilModel->pixel[i_].h[i];
                thetaSwatre->Drc = FindValue(hSwatre->Drc, SwatreSoilModel->pixel[i_].profile->horizon[i], H_COL, THETA_COL);
            }}
            report(*hSwatre, hname);
            report(*thetaSwatre, tname);
        }
    }
}
