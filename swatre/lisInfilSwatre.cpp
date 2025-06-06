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
    // int numThreads = omp_get_max_threads();
    // QVector<NODES> threadBuffers;
    // threadBuffers.reserve(numThreads);

    // for (int i = 0; i < numThreads; ++i) {
    //     threadBuffers.append(NODES(MAX_NODES+3));
    // }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        // profile 0 is for impermeable surfaces
        if (ProfileID->Drc <= 0 || fractionImperm->Drc > 0.999) {
            InfilVol->Drc = 0;
            //continue;
        } else {
            int tid = omp_get_thread_num();
            NODES& local = threadBuffers[tid];

            // Reset all vectors to zero before use
            local.theta.fill(0.0);
            local.kavg.fill(0.0);
            local.k.fill(0.0);
            local.C.fill(0.0);
            local.thetaPrev.fill(0.0);
            local.h.fill(0.0);
            local.hPrev.fill(0.0);
            local.dz.fill(0.0);
            local.disZ.fill(0.0);
            local.S.fill(0.0);
            local.thoma.fill(0.0);
            local.thomb.fill(0.0);
            local.thomc.fill(0.0);
            local.thomf.fill(0.0);
            local.beta.fill(0.0);

            double tilevol = 0;
            double theta = 0;
            double perc = 0;
            double WHorig;
            if (FloodDomain->Drc == 0)
                WHorig = WH->Drc;
            else
                WHorig = hmx->Drc;

            SwatreSoilModel->pixel[i_].wh = WHorig;    // WH is in m, convert to cm
            SwatreSoilModel->pixel[i_].tiledrain = 0;

            ComputeForPixel(i_, SwatreSoilModel, local);

            double WHN = SwatreSoilModel->pixel[i_].wh;

            theta = SwatreSoilModel->pixel[i_].thetaroot;
            perc = SwatreSoilModel->pixel[i_].percolation;

            if (SwitchIncludeTile)
                tilevol = SwatreSoilModel->pixel[i_].tiledrain;  // is already in m3

            //TODO test infil swatre for crusts and compaction
            if (SwitchInfilCrust) {
                if (SwitchDynamicCrusting && ProfileIDCrust->Drc > 0) {
                    CrustFraction->Drc = std::min(1.0, CrustFraction0->Drc + (1.0-exp(-0.2*std::max(0.0, RainCumCrust->Drc*1000))));
                }

                if (ProfileIDCrust->Drc > 0 && CrustFraction->Drc > 0) {
                    SwatreSoilModelCrust->pixel[i_].wh = WHorig;    // WH is in m, convert to cm
                    SwatreSoilModelCrust->pixel[i_].tiledrain = 0;

                  //  ComputeForPixel(i_, SwatreSoilModelCrust);

                    double WHcrust = SwatreSoilModelCrust->pixel[i_].wh;
                    WHN = WHcrust*CrustFraction->Drc + WHN*(1-CrustFraction->Drc);
                    // weighed average

                    theta = SwatreSoilModelCrust->pixel[i_].thetaroot*CrustFraction->Drc + theta*(1-CrustFraction->Drc);
                    perc = SwatreSoilModelCrust->pixel[i_].percolation*CrustFraction->Drc + perc*(1-CrustFraction->Drc);

                    if (SwitchIncludeTile) {
                        tilevol = CrustFraction->Drc*SwatreSoilModelCrust->pixel[i_].tiledrain + tilevol*(1-CrustFraction->Drc);
                    }
                }
            }

            if (SwitchInfilCompact) {
                if (ProfileIDCompact->Drc > 0 &&  CompactFraction->Drc > 0) {

                    SwatreSoilModelCompact->pixel[i_].wh = WHorig;    // WH is in m, convert to cm
                    SwatreSoilModelCompact->pixel[i_].tiledrain = 0;

               //     ComputeForPixel(i_, SwatreSoilModelCompact);

                    double WHcompact = SwatreSoilModelCompact->pixel[i_].wh;
                    WHN = WHcompact*CompactFraction->Drc + WHN*(1-CompactFraction->Drc);
                    // weighted average
                    theta = SwatreSoilModelCompact->pixel[i_].thetaroot*CrustFraction->Drc + theta*(1-CrustFraction->Drc);
                    perc = SwatreSoilModelCompact->pixel[i_].percolation*CrustFraction->Drc + perc*(1-CrustFraction->Drc);
                    if (SwitchIncludeTile) {
                        tilevol = CompactFraction->Drc*SwatreSoilModelCompact->pixel[i_].tiledrain + tilevol*(1-CompactFraction->Drc);
                    }
                }
            }

            if (SwitchGrassStrip) {
                if (ProfileIDGrass->Drc > 0 &&  GrassFraction->Drc > 0) {
                    SwatreSoilModelGrass->pixel[i_].wh = WHorig;    // WH is in m, convert to cm
                    SwatreSoilModelGrass->pixel[i_].tiledrain = 0;

                 //   ComputeForPixel(i_, SwatreSoilModelGrass);

                    double WHgrass = SwatreSoilModelGrass->pixel[i_].wh;
                    WHN = WHgrass*GrassFraction->Drc + WHN*(1-GrassFraction->Drc);
                    theta = SwatreSoilModelGrass->pixel[i_].thetaroot*CrustFraction->Drc + theta*(1-CrustFraction->Drc);
                    perc = SwatreSoilModelGrass->pixel[i_].percolation*CrustFraction->Drc + perc*(1-CrustFraction->Drc);
                    if (SwitchIncludeTile) {
                        tilevol = GrassFraction->Drc*SwatreSoilModelGrass->pixel[i_].tiledrain + tilevol*(1-GrassFraction->Drc);
                    }
                }
            }

            if (FloodDomain->Drc == 0)
                WH->Drc = WHN;
            else
                hmx->Drc = WHN;
            hmxWH->Drc = /*hmx->Drc+*/  WH->Drc;

            WaterVolall->Drc = hmxWH->Drc*CHAdjDX->Drc;

            InfilVol->Drc = (WHorig - WHN) * FlowWidth->Drc * DX->Drc;
     //       if (WHorig - WHN < 0)
     //           qDebug() << r << c << WHorig << WHN << (WHorig - WHN) << fractionImperm->Drc << FlowWidth->Drc;
    //        InfilVol->Drc = std::max(0.0, WHorig - WHN) * FlowWidth->Drc * DX->Drc;
            // use flowwidth because impermeable is done separately

            ThetaI1a->Drc = theta;
            Perc->Drc = perc/_dt; //from m to m/sec
            if (SwitchIncludeTile)
                TileWaterVolSoil->Drc = tilevol;
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
