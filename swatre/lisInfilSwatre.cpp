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
        double theta = 0;
        double perc = 0;
        double WHorig;
        if (FloodDomain->Drc == 0)
            WHorig = WH->Drc;
        else
            WHorig = hmx->Drc;

        PIXEL_INFO *pix = &SwatreSoilModel->pixel[i_];
        pix->wh = WHorig*100;    // WH is in m, convert to cm
        pix->tiledrain = 0;

        ComputeForPixel(pix);

        double WHN = pix->wh*0.01;
        theta = pix->thetaroot;
        perc = pix->percolation*0.01;
        if (SwitchIncludeTile)
            tilevol = pix->tiledrain;  // is already in m3

        //TODO test infil swatre for crusts and compaction
        if (SwitchInfilCrust) {
            if (SwitchDynamicCrusting && ProfileIDCrust->Drc > 0) {
                CrustFraction->Drc = std::min(1.0, CrustFraction0->Drc + (1.0-exp(-0.2*std::max(0.0, RainCumCrust->Drc*1000))));
            }

            if (ProfileIDCrust->Drc > 0 && CrustFraction->Drc > 0) {
                PIXEL_INFO *pixcr = &SwatreSoilModelCrust->pixel[i_];
                pixcr->wh = WHorig*100;    // WH is in m, convert to cm
                pixcr->tiledrain = 0;

                ComputeForPixel(pixcr);

                // weighed average
                WHN = pixcr->wh*0.01*CrustFraction->Drc + WHN*(1-CrustFraction->Drc);
                theta = pixcr->thetaroot*CrustFraction->Drc + theta*(1-CrustFraction->Drc);
                perc = pixcr->percolation*0.01*CrustFraction->Drc + perc*(1-CrustFraction->Drc);

                if (SwitchIncludeTile) {
                    tilevol = CrustFraction->Drc*pixcr->tiledrain + tilevol*(1-CrustFraction->Drc);
                }
            }
        }

        if (SwitchInfilCompact) {
            if (ProfileIDCompact->Drc > 0 &&  CompactFraction->Drc > 0) {
                PIXEL_INFO *pixcm = &SwatreSoilModelCompact->pixel[i_];
                pixcm->wh = WHorig*100;    // WH is in m, convert to cm
                pixcm->tiledrain = 0;

                ComputeForPixel(pixcm);

                WHN = pixcm->wh*0.01*CompactFraction->Drc + WHN*(1-CompactFraction->Drc);
                theta = pixcm->thetaroot*CompactFraction->Drc + theta*(1-CompactFraction->Drc);
                perc = pixcm->percolation*0.01*CompactFraction->Drc + perc*(1-CompactFraction->Drc);

                if (SwitchIncludeTile) {
                    tilevol = CompactFraction->Drc*pixcm->tiledrain + tilevol*(1-CompactFraction->Drc);
                }
            }
        }

        if (SwitchGrassStrip) {
            if (ProfileIDGrass->Drc > 0 &&  GrassFraction->Drc > 0) {
                PIXEL_INFO *pixgr = &SwatreSoilModelGrass->pixel[i_];
                pixgr->wh = WHorig*100;    // WH is in m, convert to cm
                pixgr->tiledrain = 0;

                ComputeForPixel(pixgr);

                WHN = pixgr->wh*0.01*GrassFraction->Drc + WHN*(1-GrassFraction->Drc);
                theta = pixgr->thetaroot*GrassFraction->Drc + theta*(1-GrassFraction->Drc);
                perc = pixgr->percolation*0.01*GrassFraction->Drc + perc*(1-GrassFraction->Drc);

                if (SwitchIncludeTile) {
                    tilevol = GrassFraction->Drc*pixgr->tiledrain + tilevol*(1-GrassFraction->Drc);
                }
            }
        }

        if (FloodDomain->Drc == 0)
            WH->Drc = WHN;
        else
            hmx->Drc = WHN;
        hmxWH->Drc = hmx->Drc + WH->Drc;
        WaterVolall->Drc = hmxWH->Drc*CHAdjDX->Drc;
        InfilVol->Drc = std::max(0.0, WHorig - WHN) * FlowWidth->Drc * DX->Drc;
        // use flowwidth because impermeable is done separately

        ThetaI1a->Drc = theta;
        Perc->Drc = perc;
        if (SwitchIncludeTile)
            TileWaterVolSoil->Drc = tilevol;

        //find depth wetting front, estimated at depth where h is initial value, very crude
        Lw->Drc = 0;
        for (int j = 0; j < pix->profile->zone->nrNodes; j++) {
//            if (j > 0 && (pix->h[j] > inith->at(j)->Drc+1.0 || pix->h[j] == 0)) {
              if (j > 0 && pix->h[j] > -10) {
                double l1 = pix->profile->zone->endComp[j-1]*0.01; // in m
                double l2 = pix->profile->zone->endComp[j]*0.01; // in m
                Lw->Drc = 0.5*(l1+l2);
            }
        }
    }}


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
