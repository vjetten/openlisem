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

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

#define QUNIT (QUnits == 1 ? 1.0 : 1000)

//---------------------------------------------------------------------------
/// Report maps for totals and mapseries (like report in PCRaster)
/// output filenames are fixed, cannot be changed by the user
/// outputnames that start with "out" are series
void TWorld::ReportMaps(void)
{
    if(SwitchInfiltration && InfilMethod != INFIL_SWATRE) {
        avgTheta();
    }
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        COMBO_V->Drc = V->Drc < 1e-5 ? 0 : V->Drc;
        VH->Drc = COMBO_V->Drc * hmxWH->Drc;
        Lwmm->Drc = Lw->Drc *1000;
    }}
report(*Lwmm,"LW");

    if(SwitchErosion)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            COMBO_SS->Drc = 0;
            COMBO_BL->Drc = 0;
            COMBO_TC->Drc = 0;

            COMBO_SS->Drc += SSFlood->Drc;
            COMBO_SS->Drc += Sed->Drc;

            COMBO_TC->Drc += SSTCFlood->Drc;
            COMBO_TC->Drc += TC->Drc;

            if (SwitchUse2Phase) {
                COMBO_BL->Drc += BLFlood->Drc;
                COMBO_TC->Drc += BLTCFlood->Drc;
            }

            if(SwitchIncludeChannel)
            {
                COMBO_SS->Drc += ChannelSSSed->Drc;
                if (SwitchUse2Phase)
                    COMBO_BL->Drc += ChannelBLSed->Drc;
                COMBO_TC->Drc += ChannelTC->Drc;
            }

            COMBO_SS->Drc = COMBO_SS->Drc  < 1e-6 ? 0 : COMBO_SS->Drc;
            COMBO_BL->Drc = COMBO_BL->Drc  < 1e-6 ? 0 : COMBO_BL->Drc;
        }}
    }

    // MAP DISPLAY VARIABLES
    if(SwitchInfiltration && InfilMethod != INFIL_SWATRE) {
        avgTheta();
    }


    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tm->Drc = (RainCumFlat->Drc)*1000.0;// + SnowmeltCum->Drc*DX->Drc/_dx) * 1000.0; // m to mm
    }}
    report(*tm, rainfallMapFileName);

    report(*InterceptionmmCum, interceptionMapFileName);

    report(*InfilmmCum, infiltrationMapFileName);

   // report(*runoffTotalCell, runoffMapFileName); // in mm, total runoff from cell (but there is also runon!)

    report(*Qm3total, runoffMapFileName); // in m3 total for this run

    report(*WHmax, floodWHmaxFileName);
    // report(*floodHmxMax, floodWHmaxFileName);  // BOTH overland flow and flood for all combinations

    report(*Qm3max,"qm3smax.map");

    if (SwitchGridRetention)
        report(*GridRetentionAct,"retentionm3.map");

    // max velocity on land in m/s
    report(*floodVMax, floodMaxVFileName);  // BOTH overland flow and flood for all combinations
    report(*floodVHMax, floodMaxVHFileName);  // momentum of all flow

    if (SwitchIncludeChannel)
    {
        report(*ChannelQntot, channelDischargeMapFileName);
        // total flow in river, cumulative during run, in m3 !!!

        report(*maxChannelflow, floodMaxQFileName);
        report(*maxChannelWH, floodMaxChanWHFileName);
    }

    if (SwitchIncludeStormDrains || SwitchIncludeTile)
    {
        report(*TileWaterVol, tileWaterVolfilename);
        // ADD SOIL
    }

    report(*floodTime, floodTimeFileName);
    report(*floodTimeStart, floodFEWFileName);

    if (SwitchGWflow)
        report(*GWWH,"groundwater.map");

    //===== SEDIMENT =====
    if(SwitchErosion)
    {
        double factor = 1.0;
        if(ErosionUnits == 2)
            factor = 1.0/(_dx*_dx);  //kg/m2
        else
            if (ErosionUnits == 0)
                factor = 10.0/(_dx*_dx); //ton/ha

        // all detachment combined
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            tm->Drc =qMax(0.0,TotalSoillossMap->Drc)*factor;
            tma->Drc =qMin(0.0,TotalSoillossMap->Drc)*factor;
        }}
        report(*tm, totalErosionFileName);
        // all deposition combined
        report(*tma, totalDepositionFileName);
        // all channel depostion combined

        if (SwitchIncludeChannel)
        {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (ChannelWidth->Drc > 0) {
                    tm->Drc =qMax(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc)*factor;
                    tma->Drc =qMin(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc)*factor;
                } else {
                    tm->Drc = 0;
                    tma->Drc = 0;
                }
            }}
            report(*tm, totalChanErosionFileName);
            report(*tma, totalChanDepositionFileName);
        }

        //copy(*tm, *TotalSoillossMap);
        //calcValue(*tm, factor, MUL);
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            tm->Drc = TotalSoillossMap->Drc  * factor;
        }}
        report(*tm, totalSoillossFileName);

        // total sediment

    }
}
//---------------------------------------------------------------------------
void TWorld::ReportMapSeries(void)
{
  //  qDebug << Outrunoff << Outwh << OutInt << Outvelo << Outinf << Outss;
    //discharge l/s or m3/s
    if (SwitchOutrunoff)
        report(*Qoutput, Outrunoff);
    // water height m
    if (SwitchOutwh)
        report(*hmxWH, Outwh);
    // interception mmtile
    if (SwitchOutInt)
        report(*InterceptionmmCum, OutInt);
    // velovity m/s
    if (SwitchOutvelo)
        report(*V, Outvelo);

    // infiltration mm
    if (SwitchOutinf)
        report(*InfilmmCum, Outinf);

    // surface storage (mm)
    if (SwitchOutss)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            tm->Drc = WHstore->Drc  * 1000;
        }}
        report(*tm, Outss);
    }

    if (SwitchIncludeTile|| SwitchIncludeStormDrains) {
        if (SwitchOutTiledrain) {
            if (QUnits == 1)
                report(*TileQn, OutTiledrain); //in m3/s
            else {
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    tm->Drc = TileQn->Drc  * 1000;
                }}
                report(*tm, OutTiledrain); //in l/s
            }
            if (SwitchOutTileVol) {
                report(*TileWaterVol, OutTileVol); //in m3
            }
        }
    }

    if (SwitchOutTheta) {
        if (SwitchInfiltration && InfilMethod != INFIL_SWATRE) { //InfilMethod != INFIL_NONE
            report(*ThetaI1a, OutTheta1);
            if (SwitchTwoLayer)
                report(*ThetaI2a, OutTheta2);
        }
    }

    if (SwitchOutGW && SwitchGWflow) {
            report(*GWWH, OutGW);
    }

    //===== SEDIMENT =====
    if(SwitchErosion)
    {
        double factor = 1.0;
        if(ErosionUnits == 2)
            factor = 1.0/(_dx*_dx);  //kg/m2
        else
            if (ErosionUnits == 0)
                factor = 10.0/(_dx*_dx); //ton/ha

        if (SwitchOutDet) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc =qMax(0.0,TotalSoillossMap->Drc)*factor;
            }}
            report(*tm, Outeros); // in units
        }

        // all deposition combined

        if (SwitchOutDep) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc =qMin(0.0,TotalSoillossMap->Drc)*factor;
            }}
            report(*tm, Outdepo); // in units
        }

        if (SwitchOutSL) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc =TotalSoillossMap->Drc*factor;
            }}
            report(*tm, OutSL);      // in user units
        }

        // total sediment
        if (SwitchOutSed) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc = (COMBO_SS->Drc + COMBO_BL->Drc)*factor;
            }}
            report(*tm, OutSed);      // in user units
        }
        if (SwitchOutConc) report(*TotalConc, Outconc);  // in g/l
        if (SwitchOutTC) report(*COMBO_TC, Outtc);      // in g/l

        if(SwitchUse2Phase) {
            if (SwitchOutSedSS) {
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    tm->Drc = COMBO_SS->Drc*factor;
                }}
            report(*tm, OutSedSS);      // in user units
            }
            if (SwitchOutSedBL) {
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    tm->Drc = COMBO_BL->Drc*factor;
                }}
                report(*tm, OutSedBL);      // in user units
            }
        }
    }
}
//---------------------------------------------------------------------------

