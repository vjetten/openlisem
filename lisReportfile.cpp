
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisReportfile.cpp
  \brief reporting maps, hydrographs, outlet/area totals and land unit stats

functions: \n
- void TWorld::OutputUI() fill output structure 'op' with results to talk to the interface.\n
- void TWorld::reportAll() \n
- void TWorld::ReportTimeseriesNew() report outlet data to textfile\n
- void TWorld::ReportTotalsNew() report totals to text file\n
- void TWorld::ReportMaps() report maps and mapseries\n
- void TWorld::CountLandunits() make a list of landunit numbers\n
- void TWorld::ReportLandunits() report text data per landunit\n
 */

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

//---------------------------------------------------------------------------
/// report to disk: timeseries at output points, totals, map series and land unit stats
void TWorld::reportAll(void)
{

    ReportTotalsNew();
    // report totals to a text file

    ReportTotalSeries();

    ReportTimeseriesNew();
    // report hydrographs ande sedigraphs at all points in outpoint.map

    if (!SwitchEndRun) {
        ReportMaps();
        ReportMapSeries();
    }
    // report all maps and mapseries

    ReportLandunits();
    // reportc stats per landunit class

    ChannelFloodStatistics();
    // report buildings submerged in flood level classes in 5cm intervals
}
//---------------------------------------------------------------------------
/** fill output structure 'op' with results to talk to the interface:
    report to screen, hydrographs and maps */
void TWorld::OutputUI(void)
{

    op.timestep = this->_dt/60.0;

    op.t = time_ms.elapsed()*0.001/60.0;
    op.t = omp_get_wtime()/60.0 - startTime;
    op.time = time/60;
    op.maxtime = op.t/runstep * op.maxstep;
    op.dx = _dx;
    op.runstep = runstep;
    op.maxstep = (int) ((EndTime-BeginTime)/_dt);
    op.EndTime = EndTime/60.0;
    op.CatchmentArea = CatchmentArea;

    op.Pmm = (RainAvgmm + SnowAvgmm)*3600/_dt;
    op.RainTotmm = RainTotmm + SnowTotmm;
    op.RainpeakTime = RainpeakTime/60;
    op.Rainpeak = Rainpeak;

    op.InfilTotmm = InfilTotmm;
    op.InfilKWTotmm = InfilKWTot; // infil part in kin wave not used

    op.SurfStormm = SurfStoremm;

    op.IntercTotmm = IntercTotmm;
    op.IntercLitterTotmm = IntercLitterTotmm;
    op.IntercHouseTotmm = IntercHouseTotmm;

    op.RunoffFraction = 0;
    if (op.RainTotmm > 0)
        op.RunoffFraction = std::max(0.0, (op.Qtotmm - op.BaseFlowtotmm)/op.RainTotmm);
    op.WaterVolTotmm = WaterVolRunoffmm;
    op.StormDrainTotmm = StormDrainTotmm;
    op.ChannelVolTotmm = ChannelVolTotmm;
    op.BaseFlowtotmm = BaseFlowTot * 1000.0/(_dx*_dx*nrCells);
    op.volFloodmm = floodVolTotmm;
    op.FloodTotMax = floodVolTotMax;
    op.FloodAreaMax = floodAreaMax;
    op.FloodArea = floodArea;

    op.Qtotmm = Qtotmm;
    op.Qtot = Qtot; // all outflow through channel and runoff for all open and outlets boundaries

    op.floodBoundaryTot = floodBoundaryTot;
    op.Qtile = QTiletot*1000.0/_dt;  //average tile output over all tile outlets as a flox in l/s
    op.Qtiletot = QTiletot;  //average tile output over all tile outlets as a flux in m3/s
    op.MB = MB;

    op.MBs = MBs;
    op.DetTotSplash = DetSplashTot*0.001; // convert from kg to ton per cell
    op.DetTotFlow = DetFlowTot*0.001;// + FloodDetTot*0.001; // convert from kg to ton
    op.DepTot = DepTot*0.001;// + FloodDepTot*0.001; // convert from kg to ton
    op.SedTot = SedTot*0.001;// + FloodSedTot*0.001; // convert from kg to ton

    op.ChannelDetTot = ChannelDetTot*0.001; // convert from kg to ton
    op.ChannelDepTot = ChannelDepTot*0.001; // convert from kg to ton
    op.ChannelSedTot = ChannelSedTot*0.001; // convert from kg to ton

    op.FloodDepTot = FloodDepTot*0.001;
    op.FloodDetTot = FloodDetTot*0.001;
    op.FloodSedTot = FloodSedTot*0.001;
    op.SoilLossTot = (SoilLossTot)*0.001; // convert from kg to ton
    op.floodBoundarySedTot = floodBoundarySedTot; // not used

//    if (noInterface)
//        return;
    //hydrographs

    op.OutletQ.at(0)->append(QtotT * 1000.0/_dt); //QtotT is in m3
    op.OutletQs.at(0)->append(SoilLossTotT);
    op.OutletC.at(0)->append(QtotT > MIN_FLUX? SoilLossTotT/QtotT : 0);
    op.OutletQtot.replace(0,Qtot);
    op.OutletQstot.replace(0,SoilLossTot/1000.0);


    double channelwh = 0;
    if(SwitchIncludeChannel) {
        FOR_ROW_COL_LDDCH5 {
                channelwh += ChannelWH->Drc;
        }}
    }
    op.OutletChannelWH.at(0)->append(channelwh);
    for(int j = 1; j < op.OutletIndices.length(); j++)
    {
        int r = op.OutletLocationX.at(j);
        int c = op.OutletLocationY.at(j);
        double discharge = Qoutput->Drc; //sum of current Qn, ChannelQn, Qflood in l/s, not Tile!
        double sedimentdischarge = SwitchErosion? Qsoutput->Drc  : 0.0; // in kg/s   * _dt
        double sedimentconcentration = SwitchErosion? TotalConc->Drc : 0.0;
        double channelwh = SwitchIncludeChannel? ChannelWH->Drc : 0.0;

        op.OutletQtot.replace(j,op.OutletQtot.at(j) + _dt * discharge/1000.0); //cumulative in m3/s
        op.OutletQstot.replace(j,op.OutletQstot.at(j) + sedimentdischarge/1000.0);
        op.OutletQ.at(j)->append(discharge);
        op.OutletQs.at(j)->append(sedimentdischarge);
        op.OutletC.at(j)->append(sedimentconcentration);
        op.OutletChannelWH.at(j)->append(std::isnan(channelwh)?0.0:channelwh);

    }

    for(int j = 0; j < op.OutletIndices.length(); j++)
    {
        bool peak = op.OutletQpeak.at(j) < op.OutletQ.at(j)->at(op.OutletQ.at(j)->length()-1);
        if(peak)
        {
            op.OutletQpeak.replace(j,op.OutletQ.at(j)->at(op.OutletQ.at(j)->length()-1));
            op.OutletQpeaktime.replace(j,time/60);

        }
    }

    //output maps
/*
    if(SwitchIncludeChannel)
    {
        if (SwitchChannelExtended)
        {
            fill(*extQCH,0.0);
            DistributeOverExtendedChannel(ChannelQn,extQCH);
            fill(*extWHCH,0.0);
            DistributeOverExtendedChannel(ChannelWH,extWHCH);
            fill(*extVCH,0.0);
            DistributeOverExtendedChannel(ChannelV,extVCH);
        } else {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                extQCH->Drc = ChannelQn->Drc;
                extWHCH->Drc = ChannelWH->Drc;
                extVCH->Drc = ChannelV->Drc;
            }}
        }
    }
*/
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
//        COMBO_VOFCH->Drc = V->Drc;
//        if (SwitchIncludeChannel) {
//            if (ChannelFlowWidth->Drc > 0)
//                COMBO_VOFCH->Drc = extVCH->Drc;
//        }
//        if(COMBO_VOFCH->Drc < 1e-6)
//            COMBO_VOFCH->Drc = 0;

        VH->Drc = V->Drc * hmxWH->Drc;
    }}

    if(SwitchErosion)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            COMBO_SS->Drc = 0;
            COMBO_BL->Drc = 0;
            COMBO_TC->Drc = TC->Drc;

            COMBO_SS->Drc += SSFlood->Drc;
            COMBO_TC->Drc += SSTCFlood->Drc;
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
            COMBO_SED->Drc = COMBO_SS->Drc + COMBO_BL->Drc;

        }}
    }

    //output maps for combo box
    for(int i = 0; i < op.ComboMapsSafe.length(); i++)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            op.ComboMapsSafe[i]->Drc = op.ComboMaps[i]->Drc * op.ComboScaling.at(i);
        }}
    }

    // ONLY ONCE
    if (runstep <= 1) {
        copy(*op.baseMap, *ShadeBW);
        copy(*op.baseMapDEM, *DEM);

        if (SwitchIncludeChannel) {
            copy(*op.channelMap, *LDDChannel);//*ChannelMaskExtended);
        }
        copy(*op.outletMap, *PointMap);

        if (SwitchRoadsystem) {
            copy(*op.roadMap, *RoadWidthDX);
        }
        if (SwitchHouses)
            copy(*op.houseMap, *HouseCover);

        if(SwitchFlowBarriers)
        {
            fill(*tma,0.0);
            FOR_ROW_COL_MV
            {
                tma->Drc = std::max(std::max(std::max(FlowBarrierN->Drc,FlowBarrierE->Drc),FlowBarrierW->Drc),FlowBarrierS->Drc);
            }
            copy(*op.flowbarriersMap,*tma);
        }
    }
    // MAP DISPLAY VARIABLES

}
//---------------------------------------------------------------------------
void TWorld::ReportTotalSeries(void)
{
    int DIG = ReportDigitsOut;
    QString newname1, pnr, sep = (SwitchWritePCRtimeplot ? " " : ",");
    int width = (!SwitchWritePCRtimeplot ? 0 : 3+DIG-3);

    newname1 = resultDir + totalSeriesFileName;
    // use simply resultdir + filename

    if (SwitchWriteHeaders) //  make file at first timestep
    {
        QFile fout(newname1);
        fout.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream out(&fout);

        out << "LISEM run - " << op.runfilename << "\n";
        out << "Time(min)" << sep << "P(mm)" << sep << "Ic(mm)";
        if (SwitchLitter)
            out << sep << "Ic(litter)(mm)";
        if (SwitchHouses)
            out << sep << "Ic(roof)(mm)";
        out << sep << "Inf(mm)";
        out << sep << "SS(mm)";
        if (SwitchIncludeStormDrains)
            out << sep << "StormDrain(mm)";
        out << sep << "Runoff(mm)";
        out << sep << "Flood(mm)";
        out << sep << "Flood area(m2)";
        out << sep << "Channels(mm)";
        out << sep << "Outflow(mm)";
        if (SwitchErosion) {
            out << sep << "Splash(ton)";
            out << sep << "FlowDet(ton)";
            out << sep << "Dep(ton)";
            out << sep << "Sed(ton)";
            out << sep << "ChanDet(ton)";
            out << sep << "ChanDep(ton)";
            out << sep << "ChanSed(ton)";
            out << sep << "FloodDet(ton)";
            out << sep << "FloodDep(ton)";
            out << sep << "FloodSed(ton)";
            out << sep << "SoilLoss(ton)";
        }
        out << "\n";
        fout.flush();
        fout.close();
    }


    QFile fout(newname1);
    fout.open(QIODevice::Append | QIODevice::Text);
    QTextStream out(&fout);
    out.setRealNumberPrecision(DIG);
    out.setFieldWidth(width);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    out << time/60;
    out << sep << op.RainTotmm;
    out << sep << op.IntercTotmm;
    if (SwitchLitter)
        out << sep << op.IntercLitterTotmm;
    if (SwitchHouses)
        out << sep << op.IntercHouseTotmm;
    out << sep << op.InfilTotmm;
    out << sep << op.SurfStormm;
    if (SwitchIncludeStormDrains)
    out << sep << op.StormDrainTotmm;
    out << sep << op.WaterVolTotmm;
    out << sep << op.volFloodmm;
    out << sep << op.FloodArea;
    out << sep << op.ChannelVolTotmm;
    out << sep << op.Qtotmm;
    if (SwitchErosion) {
        out << sep << op.DetTotSplash;
        out << sep << op.DetTotFlow;
        out << sep << op.DepTot;
        out << sep << op.SedTot;
        out << sep << op.ChannelDetTot;
        out << sep << op.ChannelDepTot;
        out << sep << op.ChannelSedTot;
        out << sep << op.FloodDetTot;
        out << sep << op.FloodDepTot;
        out << sep << op.FloodSedTot;
        out << sep << op.SoilLossTot;
    }
    out << "\n";

    /*
    if (SwitchErosion) {
        out << "\n";
        out << "\"Splash detachment (land) (ton):\"," << op.DetTotSplash<< "\n";
        out << "\"Flow detachment (land) (ton):\"," << op.DetTotFlow<< "\n";
        out << "\"Deposition (land) (ton):\"," << op.DepTot<< "\n";
        out << "\"Sediment (land) (ton):\"," << op.SedTot<< "\n";
        out << "\"Flow detachment (channels) (ton):\"," << op.ChannelDetTot<< "\n";
        out << "\"Deposition (channels) (ton):\"," << op.ChannelDepTot<< "\n";
        out << "\"Sediment (channels) (ton):\"," << op.ChannelSedTot<< "\n";
        out << "\"Flow detachment (flood) (ton):\"," << op.FloodDetTot<< "\n";
        out << "\"Deposition (flood) (ton):\"," << op.FloodDepTot<< "\n";
        out << "\"Susp. Sediment (flood) (ton):\"," << op.FloodSedTot<< "\n";
        out << "\"Total soil loss (ton):\"," << op.SoilLossTot<< "\n";
        out << "\"Average soil loss (kg/ha):\"," << (op.SoilLossTot*1000.0)/(op.CatchmentArea/10000.0)<< "\n";
        out << "\n";
    }
    */
    fout.flush();
    fout.close();


}
//---------------------------------------------------------------------------
/** reporting timeseries for every non zero point PointMap
 - 3 types of output: PCRaster timeplot format; SOBEK input format; flat comma delimited format
 - all points in one file or each point in a separate file
 - the types should be mututally exclusive in the interface and run file
*/
void TWorld::ReportTimeseriesNew(void)
{
    int nr = 0;

    int DIG = ReportDigitsOut;
    //int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
    double RainIntavg = RainAvgmm * 3600/_dt;
    double SnowIntavg = SnowAvgmm * 3600/_dt;
    QString newname1, pnr, sep = (SwitchWritePCRtimeplot ? " " : ",");
    int width = (!SwitchWritePCRtimeplot ? 0 : 3+DIG-3);
    // NOTE if SwitchWriteCommaDelimited = true then SwitchWritePCRtimeplot = false

    double QALL = QtotT * 1000.0/_dt; // total outflow on all sides in l/s, same as point 0 in interface
    double QSALL = SoilLossTotT;

    QFileInfo fi(resultDir + outflowFileName);

    //######  open files and write headers #####//

    //PCRaster and flat format are mutually exclusive
    if (SwitchWriteHeaders) //  make file at first timestep
    {
        SwitchWriteHeaders = false;
        if (SwitchSeparateOutput) // each point separate file
        {
            FOR_ROW_COL_MV
            {
                if ( PointMap->Drc > 0 )
                {
                    newname1 = fi.path() + "/" + fi.baseName() + "_" + QString::number((int)PointMap->Drc) + "." +  fi.suffix();

                    // make filename using point number

                    QFile fout(newname1);
                    fout.open(QIODevice::WriteOnly | QIODevice::Text);
                    QTextStream out(&fout);
                    out.setRealNumberPrecision(DIG);
                    out.setFieldWidth(width);
                    out.setRealNumberNotation(QTextStream::FixedNotation);

                    // HEADERS for the 3 types
                    if (SwitchWritePCRtimeplot)  //PCRaster timeplot format, cannot be SOBEK !
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << "#LISEM flow and sed output file for point #" << pnr << "\n";

                        // nr columns is time + rain (+ maybe snow) + Q + (maybe Qs + C)
                        int nrs = 4 + (SwitchErosion ? 3 : 0);
                        if (SwitchRainfall) nrs++;
                        if (SwitchSnowmelt) nrs++;
                        if (SwitchIncludeTile) nrs++;
                        pnr.setNum(nrs);
                        out << pnr << "\n";

                        out << "run step\n";
                        if (SwitchRainfall) out << "Pavg (mm/h)\n";
                        if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                        out << "Qall (l/s)\n";
                        if (SwitchIncludeChannel) out << "Qoutlet (l/s)\n" << "chanWH (m)\n";
                        if (SwitchIncludeTile) out << "Qdrain (l/s)\n";
                        if (SwitchErosion) out << "Qsall (kg/s)\n" << "Qs (kg/s)\n" << "C (g/l)\n";

                    }
                    else // flat format, comma delimited
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << "LISEM total flow and sed output file for point " << pnr << "\n";

                        out << "Time";
                        if (SwitchRainfall) out << ",Pavg";
                        if (SwitchSnowmelt) out << ",Snowavg";
                        out << ",Qall";
                        if (SwitchIncludeChannel) out << ",Q" << ",chanWH";
                        if (SwitchIncludeTile) out << ",Qtile";
                        if (SwitchErosion) out << ",Qsall,Qs,C";
                        out << "\n";

                        out << "min";
                        if (SwitchRainfall) out << ",mm/h";
                        if (SwitchSnowmelt) out << ",mm/h";
                        out << ",l/s";
                        if (SwitchIncludeChannel) out  << ",l/s" << ",m";
                        if (SwitchIncludeTile) out << ",l/s";
                        if (SwitchErosion) out << ",kg/s,kg/s,g/l";
                        out << "\n";
                    }
                    fout.close();
                }
            }
        }  // separate
        else   // HEADERS: all points in one file
        {

            newname1 = resultDir + outflowFileName;

            QFile fout(newname1);
            fout.open(QIODevice::WriteOnly | QIODevice::Text);
            QTextStream out(&fout);
            out.setRealNumberPrecision(DIG);
            out.setFieldWidth(width);
            out.setRealNumberNotation(QTextStream::FixedNotation);

            if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
            {
                // count nr of points
                FOR_ROW_COL_MV
                        if ( PointMap->Drc > 0 ) nr++;

                // nr columns is time + rain (+ maybe snow) + nr points*(Q + Qs + C)
                int nrs = 1; // runstep
                if (SwitchRainfall) nrs++;
                if (SwitchSnowmelt) nrs++;
                nrs += (1 + 2*nr); // all water points
                if (SwitchIncludeTile) nrs+=nr; //one tile per obs point
                if (SwitchErosion)
                    nrs += (1 + 2*nr); // all sed points

                pnr.setNum(nrs);

                out << "#LISEM total flow and sed output file for all reporting points\n";
                out <<  pnr << "\n";
                out << "Time (min)\n";
                if (SwitchRainfall) out << "Pavg (mm/h)\n";
                if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                out << "Qall (l/s)\n";
                FOR_ROW_COL_MV
                        if ( PointMap->Drc > 0 )
                {
                    pnr.setNum((int)PointMap->Drc);
                    out << "Q #" << pnr <<  "(l/s)\n";
                    if (SwitchIncludeChannel) out << "chanWH #" << pnr <<  "(m)\n";
                    if (SwitchIncludeTile) out << "Qtile #" << pnr <<  "(l/s)\n";
                }
                if (SwitchErosion)
                {
                    out << "Qsall (kg/s)\n";
                    FOR_ROW_COL_MV
                            if ( PointMap->Drc > 0 )
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << "Qs #"<< pnr << "(kg/s)\n" << "C #"<< pnr << "(g/l)\n";
                    }
                }
            }
            else // flat CSV format, comma delimited
            {
                out << "#LISEM total flow and sed output file for all reporting points in map\n";

                // line 1 variables
                out << "Time";
                // precipitation
                if (SwitchRainfall) out << ",P";
                if (SwitchSnowmelt) out << ",Snow";
                // total flow over the bundary
                out << ",Qall";
                // for all points in outpoint.map
                FOR_ROW_COL_MV
                        if ( PointMap->Drc > 0 )
                {
                    pnr.setNum((int)PointMap->Drc);
                    out << ",Q #" << pnr;
                    if (SwitchIncludeChannel) out << ",chanWH #" << pnr;
                    if (SwitchIncludeTile) out << ",Qtile #" << pnr;
                }
                // sediment
                if (SwitchErosion)
                {
                    out << ",Qsall";
                    FOR_ROW_COL_MV
                            if ( PointMap->Drc > 0 )
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << ",Qs #" << pnr << ",C #" << pnr;
                    }
                }
                out << "\n";

                //line 2 units
                out << "min";
                if (SwitchRainfall) out << ",mm/h";
                if (SwitchSnowmelt) out << ",mm/h";
                out << ",l/s";
                FOR_ROW_COL_MV
                        if ( PointMap->Drc > 0 )
                {
                    pnr.setNum((int)PointMap->Drc);
                    out << ",l/s #" << pnr;
                    if (SwitchIncludeChannel) out << ",m #" << pnr;
                    if (SwitchIncludeTile) out << ",l/s #" << pnr;
                }
                if (SwitchErosion)
                {
                    out << ",kg/s";
                    FOR_ROW_COL_MV
                            if ( PointMap->Drc > 0 )
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << ",kg/s #" << pnr << ",g/l #" << pnr;
                    }
                }
                out << "\n";
            }
            fout.close();
        } // all in one
    }  // opening files and writing header

    //######  open files and append values #####//

    if (SwitchSeparateOutput)
    {
        FOR_ROW_COL_MV
        {
            if ( PointMap->Drc > 0 ) // all points in separate files
            {
                newname1 = fi.path() + "/" + fi.baseName() + "_" + QString::number((int)PointMap->Drc) + "." +  fi.suffix();

                QFile fout(newname1);
                fout.open(QIODevice::Append | QIODevice::Text);

                QTextStream out(&fout);
                out.setRealNumberPrecision(DIG);
                out.setFieldWidth(width);
                out.setRealNumberNotation(QTextStream::FixedNotation);
                if (SwitchWritePCRtimeplot)
                    out << runstep;
                else
                    out << time/60;

                if (SwitchRainfall) out << sep << RainIntavg;
                if (SwitchSnowmelt) out << sep << SnowIntavg;

                out << sep << QALL << sep << Qoutput->Drc;  //Qoutput is sum channel, of, tile

                if (SwitchIncludeChannel) out << sep << ChannelWH->Drc;
                if (SwitchIncludeTile) out << sep << TileQn->Drc*1000;
                if (SwitchErosion) out << sep << QSALL << sep << Qsoutput->Drc << sep << TotalConc->Drc;
                out << "\n";
                fout.close();
            }  // if point
        }  //rows cols
    } //switch separate
    else
    {
        // all points in one file
        newname1 = resultDir + outflowFileName;
        // use simply resultdir + filename
        QFile fout(newname1);
        fout.open(QIODevice::Append | QIODevice::Text);
        QTextStream out(&fout);
        out.setRealNumberPrecision(DIG);
        out.setFieldWidth(width);
        out.setRealNumberNotation(QTextStream::FixedNotation);



        if (SwitchWritePCRtimeplot)
            out << runstep;
        else
            out << time/60;

        if (SwitchRainfall) out << sep << RainIntavg;
        if (SwitchSnowmelt) out << sep << SnowIntavg;

        out << sep << QALL;

        FOR_ROW_COL_MV
        {
            if ( PointMap->Drc > 0 )
            {
                out << sep << Qoutput->Drc;
                if (SwitchIncludeChannel) out << sep << ChannelWH->Drc;
                if (SwitchIncludeTile) out << sep << TileQn->Drc*1000;
            }
        }

        if (SwitchErosion)
        {
            out << sep << QSALL;
            FOR_ROW_COL_MV
            {
                if ( PointMap->Drc > 0 )
                    out << sep << Qsoutput->Drc << sep << TotalConc->Drc;
            }
        }
        out << "\n";
        fout.close();
    }
}
//---------------------------------------------------------------------------
/// Report totals of the main outlet nd general values for the catchment to a comma delimited text file
void TWorld::ReportTotalsNew(void)
{
    QFile fp(resultDir + resultFileName);
    if (!fp.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&fp);
    out.setRealNumberPrecision(9);
    out.setFieldWidth(16);
    out.setRealNumberNotation(QTextStream::FixedNotation);
    out << "\"LISEM run with:\"," << op.runfilename << "\n";
    out << "\"LISEM results at time (min):\"," << op.time <<"\n";
    out << "\"Catchment area (ha):\"," << op.CatchmentArea/10000.0<< "\n";
    out << "\"Total Precipitation (mm):\"," << op.RainTotmm<< "\n";
    out << "\"Total interception(mm):\"," << op.IntercTotmm<< "\n";
    out << "\"Total Litter interception (mm):\"," << op.IntercLitterTotmm<< "\n";
    out << "\"Total House interception (mm):\"," << op.IntercHouseTotmm<< "\n";
    out << "\"Total infiltration (mm):\"," << op.InfilTotmm<< "\n";
    out << "\"Surface storage (mm):\"," << op.SurfStormm<< "\n";
    out << "\"Storm Drain (mm):\"," << op.StormDrainTotmm<< "\n";
    if (SwitchKinematic2D == K2D_METHOD_KIN) {
        out << "\"Water in overland flow (mm):\"," << op.WaterVolTotmm<< "\n";
        out << "\"Water in flood (mm):\"," << 0.0 << "\n";
    } else {
       out << QString("\"Water in overland flow (h<%1)(mm)):\",%2\n").arg(minReportFloodHeight*1000).arg(op.WaterVolTotmm);
       out << QString("\"Water in flood (h>%1) (mm)):\",%2\n").arg(minReportFloodHeight*1000).arg(op.volFloodmm);
    }
    out << "\"Water in channels (mm):\"," << op.ChannelVolTotmm<< "\n";
    out << "\"Total outflow (all flows) (mm):\"," << op.Qtotmm<< "\n";
    out << "\n";
    out << "\"Total channel+OF discharg (m3):\"," << op.Qtot<< "\n";
    out << "\"Total flood discharge (m3):\"," << op.floodBoundaryTot<< "\n";
    out << "\"Total storm drain discharge (m3):\"," << op.Qtiletot<< "\n";
    out << "\"Peak time precipitation (min):\"," << op.RainpeakTime<< "\n";
    out << "\"Total discharge/Precipitation (%):\"," << op.RunoffFraction*100<< "\n";
    out << "\"Flood volume (max level) (m3):\"," << op.FloodTotMax<< "\n";
    out << "\"Flood area (max level) (m2):\"," << op.FloodAreaMax<< "\n";
    if (SwitchErosion) {
        out << "\n";
        out << "\"Splash detachment (land) (ton):\"," << op.DetTotSplash<< "\n";
        out << "\"Flow detachment (land) (ton):\"," << op.DetTotFlow+op.FloodDetTot<< "\n";
        out << "\"Deposition (land) (ton):\"," << op.DepTot+op.FloodDepTot<< "\n";
        out << "\"Sediment (land) (ton):\"," << op.SedTot+op.FloodSedTot<< "\n";
        out << "\"Flow detachment (channels) (ton):\"," << op.ChannelDetTot<< "\n";
        out << "\"Deposition (channels) (ton):\"," << op.ChannelDepTot<< "\n";
        out << "\"Sediment (channels) (ton):\"," << op.ChannelSedTot<< "\n";
    //    out << "\"Flow detachment (flood) (ton):\"," << op.FloodDetTot<< "\n";
    //    out << "\"Deposition (flood) (ton):\"," << op.FloodDepTot<< "\n";
    //    out << "\"Susp. Sediment (flood) (ton):\"," << op.FloodSedTot<< "\n";
        out << "\"Total soil loss (ton):\"," << op.SoilLossTot<< "\n";
        out << "\"Average soil loss (kg/ha):\"," << (op.SoilLossTot*1000.0)/(op.CatchmentArea/10000.0)<< "\n";
        out << "\n";
    }
    for(int i = 0; i< op.OutletQpeak.length();i++)
    {
        out << "\"Peak discharge for outlet " + QString::number(i) +" (l/s):\"," << op.OutletQpeak.at(i)<< "\n";
    }
    for(int i = 0; i< op.OutletQpeak.length();i++)
    {
        out << "\"Peak time discharge for outlet " + QString::number(i) +" (min):\"," << op.OutletQpeaktime.at(i)<< "\n";
    }
    fp.flush();
    fp.close();
}
//---------------------------------------------------------------------------
/// Report maps for totals and mapseries (like report in PCRaster)
/// output filenames are fixed, cannot be changed by the user
/// outputnames that start with "out" are series
void TWorld::ReportMaps(void)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tm->Drc = (RainCumFlat->Drc + SnowmeltCum->Drc*DX->Drc/_dx) * 1000.0; // m to mm
    }}
    report(*tm, rainfallMapFileName);

    report(*InterceptionmmCum, interceptionMapFileName);

    report(*InfilmmCum, infiltrationMapFileName);

    report(*runoffTotalCell, runoffMapFileName); // in mm, total runoff from cell (but there is also runon!)

    report(*WHmax, floodWHmaxFileName);
    // report(*floodHmxMax, floodWHmaxFileName);  // BOTH overland flow and flood for all combinations

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
        report(*TileQmax, tileQmaxfilename);
    }

    report(*floodTime, floodTimeFileName);
    report(*floodTimeStart, floodFEWFileName);

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
            tm->Drc =std::max(0.0,TotalSoillossMap->Drc)*factor;
            tma->Drc =std::min(0.0,TotalSoillossMap->Drc)*factor;
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
                    tm->Drc =std::max(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc)*factor;
                    tma->Drc =std::min(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc)*factor;
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
    //discharge l/s
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
        report(*COMBO_VOFCH, Outvelo);

    // infiltration mm
    if (SwitchOutinf)
        report(*InfilmmCum, Outinf);

    if (SwitchOutss)
    {
        //calcMapValue(*tm, *WHstore, 1000, MUL);// in mm
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            tm->Drc = WHstore->Drc  * 1000;
        }}
        report(*tm, Outss);

    }

    if (SwitchIncludeTile|| SwitchIncludeStormDrains)
    {
        if (SwitchOutTiledrain)
        {
           // calcMapValue(*tm, *TileQn, 1000, MUL);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc = TileQn->Drc  * 1000;
            }}
            report(*tm, OutTiledrain); //in l/s
        }
        if (SwitchOutTileVol)
        {
            // report(*TileV, "tilev"); //in m3/s
            report(*TileWaterVol, OutTileVol); //in m3
        }
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
            // all detachment combined
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc =std::max(0.0,TotalSoillossMap->Drc)*factor;
            }}
            report(*tm, Outeros); // in units
        }

        // all deposition combined

        if (SwitchOutDep) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc =std::min(0.0,TotalSoillossMap->Drc)*factor;
            }}
            report(*tm, Outdepo); // in units
        }

        if (SwitchOutSL) {
            //alcValue(*tm, factor, MUL);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc =TotalSoillossMap->Drc*factor;
            }}
            report(*tm, OutSL);      // in user units
        }

        // total sediment
        if (SwitchOutSed) {
//            copy(*tm, *COMBO_SED); //kg/cell
//            calcValue(*tm, factor, MUL);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                tm->Drc = COMBO_SED->Drc*factor;
            }}
            report(*tm, OutSed);      // in user units
        }
        if (SwitchOutConc) report(*TotalConc, Outconc);  // in g/l
        if (SwitchOutTC) report(*COMBO_TC, Outtc);      // in g/l

        if(SwitchUse2Phase) {
            if (SwitchOutSedSS) {
//                copy(*tm, *COMBO_SS); //kg/cell
//                calcValue(*tm, factor, MUL);
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    tm->Drc = COMBO_SS->Drc*factor;
                }}
            report(*tm, OutSedSS);      // in user units
            }
            if (SwitchOutSedBL) {
              //  copy(*tm, *COMBO_BL); //kg/cell
              //  calcValue(*tm, factor, MUL);
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
/// Land unit statistics: count nr land units in classifiedfile
// VJ 110110 count nr of land units in classified file
void TWorld::CountLandunits(void)
{
    if (!SwitchErosion)
        return;

    int i, j;
    for (i = 0; i < NRUNITS; i++)
    {
        unitList[i].nr = 0;
        unitList[i].var0 = 0;
        unitList[i].var1 = 0;
        unitList[i].var2 = 0;
        unitList[i].var3 = 0;
        unitList[i].var4 = 0;
        unitList[i].var5 = 0;
    }

    i = 0;
    FOR_ROW_COL_MV
    {
        bool found = false;

        for(j = 0; j <= i; j++)
            if ((long)LandUnit->Drc == unitList[j].nr)
                found = true;

        if(!found && i < NRUNITS)
        {
             unitList[i].nr = (long)LandUnit->Drc;
             i++;
        }
    }
    landUnitNr = i;
}
//---------------------------------------------------------------------------
//VJ 110110
/// Report the erosion totals per land unit
void TWorld::ReportLandunits(void)
{
    if (!SwitchErosion)
        return;

    #pragma omp parallel for num_threads(userCores)
    for (int i = 0; i < landUnitNr; i++)//landUnitNr; i++)
    {
        unitList[i].var0 = 0;
        unitList[i].var1 = 0;
        unitList[i].var2 = 0;
        unitList[i].var3 = 0;
    }

   #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        //variables are kg/cell convert to ton/cell
        for (int i = 0; i < landUnitNr; i++)
            if (unitList[i].nr == (int)LandUnit->Drc) {
                unitList[i].var0 += CellArea->Drc/10000;//ha
                unitList[i].var1 += TotalDetMap->Drc/1000; //ton/cell
                unitList[i].var2 += TotalDepMap->Drc/1000;
                unitList[i].var3 += TotalSoillossMap->Drc/1000;
            }
    }}


    QString name;
    name = resultDir + QFileInfo(totalLandunitFileName).baseName()+timestampRun+".csv";
    QFile fout(name);
    fout.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&fout);
    out.setRealNumberPrecision(3);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    out << "Landunit,Area,Detachment,Deposition,Soil Loss\n";
    out << "#,ha,ton,ton,ton\n";
    for (int i = 0; i < landUnitNr; i++)
        out << unitList[i].nr << ","
            << unitList[i].var0 << ","
            << unitList[i].var1 << ","
            << unitList[i].var2 << ","
            << unitList[i].var3 << "\n";
    fout.close();

}
//---------------------------------------------------------------------------
void TWorld::ChannelFloodStatistics(void)
{
    if(SwitchKinematic2D == K2D_METHOD_KIN)
        return;

    #pragma omp parallel for num_threads(userCores)
    for (int i = 0; i < NRUNITS; i++)
    {
        floodList[i].nr = i;
        floodList[i].var0 = 0.05*i; //depth
        floodList[i].var1 = 0;
        floodList[i].var2 = 0;
        floodList[i].var3 = 0;
        floodList[i].var4 = 0;
        floodList[i].var5 = 0;
        floodList[i].var6 = 0;
    }

    int nr = 0;
    FOR_ROW_COL_MV_L {
        double area = _dx*_dx;
        if(floodHmxMax->Drc > 0)  //floodHmxMax has zero under treshold
        {
            int i = 0;
            while (floodList[i].var0 < floodHmxMax->Drc && i < NRUNITS)
                i++;
            if (i > 0)
                i--;
            nr = std::max(nr, i);
            //qDebug() << nr << i << floodHmxMax->Drc;
            floodList[i].var1 += area; // area flooded in this class
            floodList[i].var2 += area*floodHmxMax->Drc; // vol flooded in this class
            floodList[i].var3 = std::max(floodTime->Drc/60.0,floodList[i].var3); // max time in this class
            floodList[i].var4 = std::max(floodTimeStart->Drc/60.0,floodList[i].var4); // max time in this class
            if (SwitchHouses)
                floodList[i].var5 += HouseCover->Drc*area;
            if (SwitchRoadsystem)
                floodList[i].var6 += RoadWidthDX->Drc*DX->Drc;
        }
    }}

    QString name;
    name = resultDir + QFileInfo(floodStatsFileName).baseName()+timestampRun+".csv";
    QFile fp(name);
    if (!fp.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    double totarea = 0;
    double totvol = 0;
    double totbuild = 0;
    double totroad = 0;
    for (int i = 1; i < nr+1; i++)
    {
        totarea += floodList[i].var1;
        totvol += floodList[i].var2;
        totbuild += floodList[i].var5;
        totroad += floodList[i].var6;
    }

    QTextStream out(&fp);
    out.setRealNumberPrecision(2);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    out << "\"LISEM run with:," << op.runfilename << "\"" << "\n";
    out << "\"results at time (min):\"" << op.time << "\n";
    out << "class,Depth,Area,Volume,Duration,Start,Structures,Roads\n";
    out << "#,m,m2,m3,h,h,m2,m2\n";
    out << "total" << ",>0.05," << totarea << "," << totvol << ",,," << totbuild << "," << totroad <<"\n";
    for (int i = 1; i < nr+1; i++)
        out << i << ","
            << floodList[i].var0 << ","
            << floodList[i].var1 << ","
            << floodList[i].var2 << ","
            << floodList[i].var3 << ","
            << floodList[i].var4 << ","
            << floodList[i].var5 << ","
            << floodList[i].var6
            << "\n";

    fp.flush();
    fp.close();

}
//---------------------------------------------------------------------------
// Make the maps to bedrawn in the interface as a copy in the op starcture
// reason is that all pointers are destroyed after the run so when lisem finishes
// the information on the output screen points to an empty pointer
// by copying the info remains available
// initialize maps for output to screen
// must be done after Initialize Data because then we know how large the map is
void TWorld::setupDisplayMaps()
{
    if (op.baseMap != 0)
    {
        delete op.baseMap;
        delete op.baseMapDEM;
        delete op.channelMap;
        delete op.outletMap;
        delete op.roadMap;
        delete op.houseMap;
        delete op.flowbarriersMap;
    }

    op.baseMap = new cTMap();
    op.baseMapDEM = new cTMap();
    op.channelMap = new cTMap();
    op.outletMap = new cTMap();
    op.roadMap = new cTMap();
    op.houseMap = new cTMap();
    op.flowbarriersMap = new cTMap();

    op.baseMap->MakeMap(LDD, 0);
    op.baseMapDEM->MakeMap(LDD, 0);
    op.channelMap->MakeMap(LDD, 0);
    op.outletMap->MakeMap(LDD, 0);
    op.roadMap->MakeMap(LDD, 0);
    op.houseMap->MakeMap(LDD, 0);
    op.flowbarriersMap->MakeMap(LDD, 0);
}
//---------------------------------------------------------------------------
void TWorld::setupHydrographData()
{
    ClearHydrographData();


    //get the sorted locations and index numbers of the outlet points
    QList<int> nr;
    //int maxnr = 0;

    //0 is reserved for total outflow (channel and overland flow)
    nr.append(0);
    op.OutletIndices.append(0);
    op.OutletLocationX.append(0);
    op.OutletLocationY.append(0);
    op.OutletQ.append(new QVector<double>);
    op.OutletQs.append(new QVector<double>);
    op.OutletC.append(new QVector<double>);
    op.OutletChannelWH.append(new QVector<double>);
    op.OutletQpeak.append(0);
    op.OutletQpeaktime.append(0);
    op.OutletQtot.append(0);
    op.OutletQstot.append(0);

    FOR_ROW_COL_MV
    {
        if(PointMap->Drc > 0)
        {
            nr.append((int)PointMap->Drc);
            op.OutletIndices.append((int)PointMap->Drc);
            op.OutletLocationX.append(r);
            op.OutletLocationY.append(c);
            op.OutletQ.append(new QVector<double>);
            op.OutletQs.append(new QVector<double>);
            op.OutletC.append(new QVector<double>);
            op.OutletChannelWH.append(new QVector<double>);
            op.OutletQpeak.append(0);
            op.OutletQpeaktime.append(0);
            op.OutletQtot.append(0);
            op.OutletQstot.append(0);
        }
    }

    QList<int> tx;
    QList<int> ty;
    tx.clear();
    tx.append(op.OutletLocationX);
    ty.clear();
    ty.append(op.OutletLocationY);
    op.OutletLocationX.clear();
    op.OutletLocationY.clear();

    std::sort(nr.begin(), nr.end());
    for(int i = 0; i < nr.length(); i++)
    {
        int j;
        for(j = 0; j < nr.length(); j++)
        {
            if(op.OutletIndices.at(j) == nr.at(i))
            {
                break;
            }
        }
        op.OutletLocationX.append(tx.at(j));
        op.OutletLocationY.append(ty.at(j));
    }
    op.OutletIndices.clear();
    op.OutletIndices.append(nr);

}

void TWorld::ClearHydrographData()
{
    for(int i =op.OutletIndices.length() - 1; i >-1 ; i--)
    {
        delete op.OutletQ.at(i);
        delete op.OutletQs.at(i);
        delete op.OutletC.at(i);
        delete op.OutletChannelWH.at(i);
    }

    op.OutletIndices.clear();
    op.OutletLocationX.clear();
    op.OutletLocationY.clear();
    op.OutletQ.clear();
    op.OutletQs.clear();
    op.OutletC.clear();
    op.OutletQpeak.clear();
    op.OutletQpeaktime.clear();
    op.OutletChannelWH.clear();
    op.OutletQtot.clear();
    op.OutletQstot.clear();
    if(SwitchImage)
    {
        op.has_image = true;
        op.Image = RGB_Image;
    }
}

//---------------------------------------------------------------------------
void TWorld::setColor(int i)
{
    if (i == 1){  //blue red
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.01);
        Colormap.append(0.1);
        Colormap.append(0.3);
        Colormap.append(1.0);

        Colors.clear();
        //        Colors.append("#dae3ff");
        //        Colors.append("#6daaff");
        //        Colors.append("#4c76e8");
        //        Colors.append("#4c57c3");
        Colors.append("#9eccee");//#bfdcf9");
        Colors.append("#427dc6");//#b1c0e9");
        Colors.append("#204ab5");//#7b94e7");
        Colors.append("#072a9c");//#2037b5");
        Colors.append("#df2a36");
    }
    if (i == 2){ // yellow red
        //        Colormap.clear();
        //        Colormap.append(0.0);
        //        Colormap.append(0.25);
        //        Colormap.append(0.75);
        //        Colormap.append(1.0);
        //        Colors.clear();
        //        Colors.append("#00FF00");
        //        Colors.append("#FFFF00");
        //        Colors.append("#FF0000");
        //        Colors.append("#A60000");
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);

        Colors.clear();
        Colors.append("#ffffb2");
        Colors.append("#fecc5c");
        Colors.append("#fd8d3c");
        Colors.append("#f03b20");
        Colors.append("#bd0026");


    }
    if (i == 3){ //blue
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);

        Colors.clear();
        Colors.append("#9eccee");//#bfdcf9");
        Colors.append("#427dc6");//#b1c0e9");
        Colors.append("#204ab5");//#7b94e7");
        Colors.append("#072a9c");//#2037b5");
        Colors.append("#07215e");//#183280");


    }
    if (i == 4) { // yellow blue
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);

        Colors.clear();
        Colors.append("#ffff61");
        Colors.append("#c7e55a");
        Colors.append("#32b1df");
        Colors.append("#3271ca");
        Colors.append("#2c3898");
    }
    if (i == 5) { // blue green

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#8888FF");
        Colors.append("#0000FF");
        Colors.append("#008800");
    }
    if (i == 6) { // red yellow green
        //        Colormap.clear();
        //        Colormap.append(0.0);
        //        Colormap.append(0.25);
        //        Colormap.append(0.5);
        //        Colormap.append(0.75);
        //        Colormap.append(1.0);
        //        Colors.clear();
        //        Colors.append("#A60000");
        //        Colors.append("#FF0000");
        //        Colors.append("#FFFF00");
        //        Colors.append("#00FF00");
        //        Colors.append("#007300");
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);

        Colors.clear();
        //        Colors.append("#1a9641");
        //        Colors.append("#a6d96a");
        //        Colors.append("#ffffa0");
        //        Colors.append("#fdae61");
        //        Colors.append("#e0181c");
        Colors.append("#d7191c");
        Colors.append("#fdae61");
        Colors.append("#fdfd7e");
        Colors.append("#abdda4");
        Colors.append("#2b83ba");


    }

    if (i == 7) { //blue green yellow red

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);
        //        Colors.clear();
        //        Colors.append("#007300");
        //        Colors.append("#00FF00");
        //        Colors.append("#FFFF00");
        //        Colors.append("#FF0000");
        //        Colors.append("#A60000");

        Colors.clear();
        //        Colors.append("#e0181c");
        //        Colors.append("#fdae61");
        //        Colors.append("#ffffa0");
        //        Colors.append("#a6d96a");
        //        Colors.append("#1a9641");
        Colors.append("#2b83ba");
        Colors.append("#a4ddd9");//abdda4");
        Colors.append("#ffffef");
        Colors.append("#d3b03e");//#fdae61");
        Colors.append("#d7191c");
    }
    if (i == 8) {  // green white yellow red
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.3);
        Colormap.append(0.5);
        Colormap.append(0.70);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#616ca2");
        Colors.append("#50B547");
        Colors.append("#FFFFFF");
        Colors.append("#ffff88");
        Colors.append("#FF0000");

    }
    if (i == 9) {  //YELLOW ORGANDE BROWN/RED

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#ffffd4");
        Colors.append("#fed98e");
        Colors.append("#fe9929");
        Colors.append("#d95f0e");
        Colors.append("#a94400");

    }
    if (i == 10) {  // greenish colors
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#f0f9e8");
        Colors.append("#bae4bc");
        Colors.append("#7bccc4");
        Colors.append("#43a2ca");
        Colors.append("#0868ac");
    }
}

void TWorld::GetComboMaps()
{
    //combo box maps
    ClearComboMaps();

    setColor(1);
    AddComboMap(0,"Total Discharge","l/s",Qoutput,Colormap,Colors,true,false,1.0, 1.0);
  //  if (FlowBoundaryType > 0)
  //  AddComboMap(0,"Boundary Discharge","l/s",K2DQ,Colormap,Colors,true,false,1000.0, 1.0);

    setColor(3);
    AddComboMap(0,"Water Height","m",hmxWH,Colormap,Colors,false,false,1.0,0.01);
    if (Switch2DDiagonalFlow)
       AddComboMap(0,"Diagonal Discharge","l/s",Qdiag,Colormap,Colors,false,false,1.0, 0.01);
    setColor(2);
    AddComboMap(0,"Flow Velocity","m/s",V /*COMBO_VOFCH*/,Colormap,Colors,false,false,1.0, 0.01);
    AddComboMap(0,"Flow Momentum","m2/s",VH,Colormap,Colors,false,false,1.0, 0.01); //VH

    if(SwitchIncludeChannel)
    {
//        setColor(1);
//        AddComboMap(0,"Channel Discharge","l/s",extQCH,Colormap,Colors,true,false,1000.0, 1.0);
//        setColor(3);
//        AddComboMap(0,"Channel Water Height","m",extWHCH,Colormap,Colors,false,false,1.0,0.01);
//        setColor(2);
//        AddComboMap(0,"Channel Velocity","m/s",extVCH,Colormap,Colors,false,false,1.0,0.01);
        setColor(1);
        AddComboMap(0,"Channel Discharge","l/s",ChannelQn,Colormap,Colors,true,false,1000.0, 1.0);
        setColor(3);
        AddComboMap(0,"Channel Water Height","m",ChannelWH,Colormap,Colors,false,false,1.0,0.01);
        setColor(2);
        AddComboMap(0,"Channel Velocity","m/s",ChannelV,Colormap,Colors,false,false,1.0,0.01);


    }

    if(SwitchIncludeTile || SwitchIncludeStormDrains) {
        setColor(1);
        AddComboMap(0,"Storm Drain Volume","m3",TileWaterVol,Colormap,Colors,false,false,1.0,1.0);
        AddComboMap(0,"Storm Drain Discharge","l/s",TileQn,Colormap,Colors,false,false,1000.0,1.0);
    }
    if(InfilMethod != INFIL_NONE)
    {
        setColor(4);
        AddComboMap(0,"Interception","mm",InterceptionmmCum,Colormap,Colors,false,false,1.0,1.0);
        AddComboMap(0,"Infiltration","mm",InfilmmCum,Colormap,Colors,false,false,1.0,1.0);

        if (InfilMethod != INFIL_SWATRE) {
            AddComboMap(0,"Moisture content 1","-",Thetaeff,Colormap,Colors,false,false,1.0,1.0);
            if (SwitchTwoLayer)
                AddComboMap(0,"Moisture content 2","-",ThetaI2,Colormap,Colors,false,false,1.0,1.0);
            if (SwitchPercolation)
                AddComboMap(0,"Percolation","mm",PercmmCum,Colormap,Colors,false,false,1.0,1.0);
        }
    }

    setColor(5);
    double factor = 3600000.0/_dt; //from m to mm/h

    AddComboMap(0,"Rainfall Cumulative","mm",RainCumFlat,Colormap,Colors,false,false,1000.0,0.1);
    AddComboMap(0,"Rainfall Intensity","mm/h",Rain,Colormap,Colors,false,false,factor,0.1);
  //  AddComboMap(0,"ETa cumulative","mm",ETa,Colormap,Colors,false,false,1000.0,0.1);

    if (SwitchKinematic2D == K2D_METHOD_DYN || SwitchKinematic2D == K2D_METHOD_KINDYN) {
        setColor(3);
      //  QString txt = QString("Flood Height");
      //  if (SwitchKinematic2D == K2D_METHOD_DYN)
        QString txt = QString("Max flood Height (h>%1m)").arg(minReportFloodHeight);

//        AddComboMap(0,txt,"m",hmxflood,Colormap,Colors,false,false,1.0,0.01);
//        setColor(3);
        AddComboMap(0,txt,"m",floodHmxMax,Colormap,Colors,false,false,1.0,0.01);
        setColor(6);
        AddComboMap(0,"Flood Start Time","min",floodTimeStart,Colormap,Colors,false,false,1.0,1.0);
        setColor(7);
        AddComboMap(0,"Flood duration","min",floodTime,Colormap,Colors,false,false,1.0,1.0);
        setColor(6);
        if (SwitchVariableTimestep) {
            AddComboMap(0,"Timestep","s",FloodDT,Colormap,Colors,false,false,1.0,0.01);
            AddComboMap(0,"Steps pr cell","-",FloodT,Colormap,Colors,false,false,1.0,1.0);
        }
    }

    if(SwitchErosion)
    {
        double step = 0.01;
        setColor(7);

        QString unit = "kg/cell";
        double factor = 1.0;
        if(ErosionUnits == 2)
        {
            factor = 1.0/(_dx*_dx);
            unit = "kg/m2";
        }else if (ErosionUnits == 0)
        {
            factor = 10.0/(_dx*_dx);
            unit = "t/ha";
        }
        AddComboMap(1,"Total Soil Loss",unit,TotalSoillossMap,Colormap,Colors,false,true,factor, step);

        setColor(9);
        AddComboMap(1,"Splash detachment",unit,DETSplashCum,Colormap,Colors,false,false,factor, step);
        AddComboMap(1,"Flow detachment",unit,DETFlowCum,Colormap,Colors,false,false,factor, step);
        AddComboMap(1,"Sed. Concentration","kg/m3",TotalConc,Colormap,Colors,false,false,1.0, step);
        if (SwitchSedtrap)
        AddComboMap(1,"Sed trap","kg/m3",SedMaxVolume,Colormap,Colors,false,false,1.0, step);

        if(SwitchUse2Phase) {
            AddComboMap(1,"Suspended sed.",unit,COMBO_SS/*SSFlood*/,Colormap,Colors,false,false,factor, step);
            AddComboMap(1,"Bedload sed.",unit,COMBO_BL /*BLFlood*/,Colormap,Colors,false,false,factor, step);
            AddComboMap(1,"TC suspended","kg/m3",SSTCFlood,Colormap,Colors,false,false,1.0, step);
            AddComboMap(1,"TC bedload","kg/m3",BLTCFlood,Colormap,Colors,false,false,1.0, step);
         //   AddComboMap(1,"SS depth","m",SSDepthFlood,Colormap,Colors,false,false,1.0, step);
         //   AddComboMap(1,"BL depth","m",BLDepthFlood,Colormap,Colors,false,false,1.0, step);
        } else {
            AddComboMap(1,"Sediment load",unit,COMBO_SED,Colormap,Colors,false,false,factor, step);
            AddComboMap(1,"Transport Capacity","kg/m3",COMBO_TC,Colormap,Colors,false,false,1.0, step);
        }

        setColor(10);
        AddComboMap(1,"Deposition",unit,DEPCum,Colormap,Colors,false,false,-factor, step);

        if(SwitchUseMaterialDepth) {
            AddComboMap(1,"Storage",unit,Storage,Colormap,Colors,false,false,-factor, step);
        AddComboMap(1,"Storage",unit,StorageDep,Colormap,Colors,false,false,-factor, step);
        }

/*
        if(SwitchUseGrainSizeDistribution)
        {
            setColor(9);

            FOR_GRAIN_CLASSES
            {
                AddComboMap(1,"Overland S.L. Grain Class " + QString::number(d),"kg/m2",Sed_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
            }
            if(SwitchIncludeChannel)
            {
                if(SwitchUse2Phase)
                {
                    FOR_GRAIN_CLASSES
                    {
                        AddComboMap(1,"Channel BL S.L. Grain Class " + QString::number(d),"kg/m2",RBL_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
                    }
                    FOR_GRAIN_CLASSES
                    {
                        AddComboMap(1,"Channel SS S.L. Grain Class " + QString::number(d),"kg/m2",RSS_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
                    }
                }else
                {
                    FOR_GRAIN_CLASSES
                    {
                        AddComboMap(1,"Channel S.L. Grain Class " + QString::number(d),"kg/m2",RBL_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
                    }
                }
            }
            if(SwitchChannelFlood)
            {
                FOR_GRAIN_CLASSES
                {
                    AddComboMap(1,"Flood BL S.L. Grain Class " + QString::number(d),"kg/m2",BL_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
                }
                FOR_GRAIN_CLASSES
                {
                    AddComboMap(1,"Flood SS S.L. Grain Class " + QString::number(d),"kg/m2",SS_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx),step);
                }
            }

        }
        */
    }

}
//---------------------------------------------------------------------------
void TWorld::ClearComboMaps()
{

    for(int i =op.ComboMapsSafe.length() - 1; i >-1 ; i--)
    {
        delete op.ComboMapsSafe.at(i);
    }
    op.ComboMapsSafe.clear();

    op.ComboLists.clear();
    op.ComboMaps.clear();
    op.ComboMapsSafe.clear();
    op.ComboColorMap.clear();
    op.ComboColors.clear();
    op.ComboLogaritmic.clear();
    op.ComboSymColor.clear();
    op.ComboMapNames.clear();
    op.ComboUnits.clear();
    op.ComboScaling.clear();
    op.userMinV.clear();
    op.userMaxV.clear();
    op.comboStep.clear();

    op.comboboxset = false;
}
//---------------------------------------------------------------------------
void TWorld::AddComboMap(int listn, QString name, QString unit,cTMap * map,QList<double> ColorMap, QList<QString> Colors,
                         bool log,bool symcol, double scale, double step)
{
    op.ComboLists.append(listn);
    op.ComboMaps.append(map);
    op.ComboMapsSafe.append(new cTMap());
    op.ComboMapsSafe.at(op.ComboMapsSafe.length()-1)->MakeMap(LDD,0.0);
    op.ComboColorMap.append(ColorMap);
    op.ComboColors.append(Colors);
    op.ComboLogaritmic.append(log);
    op.ComboSymColor.append(symcol);
    op.ComboMapNames.append(name);
    op.ComboUnits.append(unit);
    op.ComboScaling.append(scale);  // multiplier for display
    op.userMinV.append(0);  //initialize to 0, used to save users choice
    op.userMaxV.append(0);
    op.comboStep.append(step);

    op.comboboxset = false;
}
//---------------------------------------------------------------------------
void TWorld::run()
{
    QTimer::singleShot(0, this, SLOT(DoModel()));
    exec();
}
//---------------------------------------------------------------------------
void TWorld::stop()
{
    QMutexLocker locker(&mutex);
    stopRequested = true;
}
//---------------------------------------------------------------------------
