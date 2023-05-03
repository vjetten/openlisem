
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
**  website, information and code: https://github.com/vjetten/openlisem
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

#define QUNIT (QUnits == 1 ? 1.0 : 1000)

//---------------------------------------------------------------------------
/// report to disk: timeseries at output points, totals, map series and land unit stats
void TWorld::reportAll(void)
{
    ReportTotalsNew();
    // report totals to a text file

    ReportTimeseriesNew();
    // report hydrographs ande sedigraphs at all points in outpoint.map

    ReportTotalSeries();


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
    if (SwitchEventbased)
        op.Time.append(time/60);
    else
        op.Time.append(time/86400);
    op.maxtime = op.t/runstep * op.maxstep;
    op._dx = _dx;
    op._llx = _llx;
    op._lly = _lly;
    op._nrCols = _nrCols;
    op._nrRows = _nrRows;
    op.runstep = runstep;
    op.maxstep = (int) ((EndTime-BeginTime)/_dt);
    //op.EndTime = EndTime/60.0;
    op.CatchmentArea = CatchmentArea;

    op.Pmm.append((RainAvgmm + SnowAvgmm)*3600/_dt);
    op.RainTotmm = RainTotmm + SnowTotmm;
    op.ETaTotmm = ETaTotmm;
    op.GWlevel = GWlevel;
    op.RainpeakTime = RainpeakTime/60;
    op.Rainpeak = Rainpeak;

    op.InfilTotmm = InfilTotmm;
    op.InfilKWTotmm = InfilKWTot; // infil part in kin wave not used
    op.Theta1 = theta1tot;
    op.Theta2 = theta2tot;

    op.SurfStormm = SurfStoremm;

    op.IntercTotmm = IntercTotmm + IntercETaTotmm;
    op.IntercLitterTotmm = IntercLitterTotmm;
    op.IntercHouseTotmm = IntercHouseTotmm;

    op.RunoffFraction = 0;
    if (op.RainTotmm > 0)
        op.RunoffFraction = std::max(0.0, (op.Qtotmm - op.BaseFlowTotmm)/op.RainTotmm);
    op.WaterVolTotmm = WaterVolRunoffmm;
    op.StormDrainTotmm = StormDrainTotmm;
    op.ChannelVolTotmm = ChannelVolTotmm;
    op.BaseFlowTotmm = BaseFlowTotmm;

    op.volFloodmm = floodVolTotmm;
    op.FloodTotMax = floodVolTotMax;
    op.FloodAreaMax = floodAreaMax;
    op.FloodArea = floodArea;

    op.Qtotmm = Qtotmm;
    op.Qboundtotmm = FloodBoundarymm;
    op.Qtot = Qtot; // all outflow through channel and runoff for all open and outlets boundaries

    op.floodBoundaryTot = floodBoundaryTot;
    if (SwitchIncludeStormDrains || SwitchIncludeTile)
        op.Qtile.append(QTiletot*1000.0/_dt);  //average tile output over all tile outlets as a flox in l/s
    op.Qtiletot = QTiletot;  //average tile output over all tile outlets as a flux in m3/s
    op.MB = MB;

    if (SwitchErosion) {
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

        op.OutletQs.at(0)->append(SoilLossTot_dt); //timestep output in kg! SoilLossOutlet = sum of Qs*dt and channelQs*dt and boundaryQs
        op.OutletC.at(0)->append(Qtot_dt > MIN_FLUX? SoilLossTot_dt/Qtot_dt : 0);
        op.OutletQstot.replace(0,SoilLossTot*0.001);
    }


    //hydrographs

    // outlet 0 all flow
    op.OutletQ.at(0)->append(Qtot_dt * QUNIT/_dt); //Qtot_dt is in m3
    op.OutletQtot.replace(0,Qtot); // cumulative tot outflow
    op.OutletChannelWH.at(0)->append(0);

    for(int j = 1; j < op.OutletIndices.length(); j++)
    {
        int r = op.OutletLocationX.at(j);
        int c = op.OutletLocationY.at(j);

        double channelwh = SwitchIncludeChannel? ChannelWH->Drc : 0.0;
        op.OutletChannelWH.at(j)->append(std::isnan(channelwh)?0.0:channelwh); //? why nan

        if (SwitchIncludeChannel) {
            op.OutletQtot.replace(j,op.OutletQtot.at(j) + _dt * ChannelQn->Drc); //cumulative in m3/s
            op.OutletQ.at(j)->append(ChannelQn->Drc*QUNIT);
        } else {
            op.OutletQtot.replace(j,op.OutletQtot.at(j) + _dt * Qn->Drc); //cumulative in m3/s
            op.OutletQ.at(j)->append(Qn->Drc*QUNIT);
        }

        if (SwitchErosion) {
            if (SwitchIncludeChannel) {
                op.OutletQstot.replace(j,op.OutletQstot.at(j) + ChannelQsn->Drc*_dt/1000.0); // sum in kg of OF + channel
                op.OutletQs.at(j)->append(ChannelQsn->Drc);  // in kg/s
                op.OutletC.at(j)->append(ChannelConc->Drc);  // in kg/m3 or g/l
            } else {
                op.OutletQstot.replace(j,op.OutletQstot.at(j) + Qsn->Drc*_dt/1000.0); // sum in kg of OF + channel
                op.OutletQs.at(j)->append(Qsn->Drc);  // in kg/s
                op.OutletC.at(j)->append(Conc->Drc);  // in kg/m3 or g/l
            }
        }
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

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        COMBO_V->Drc = V->Drc < 1e-5 ? 0 : V->Drc;
        VH->Drc = COMBO_V->Drc * hmxWH->Drc;
    }}

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

    //output maps for combo box
    for(int i = 0; i < op.ComboMaps.length(); i++)
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            op.ComboMapsSafe[i]->Drc = op.ComboMaps[i]->Drc; // * op.ComboScaling.at(i); scaling is done filldrawmapdata
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
            FOR_ROW_COL_MV_L {
                if (RoadWidthDX->Drc > 0.2*_dx)
                    op.roadMap->Drc = RoadWidthDX->Drc;
                else
                    op.roadMap->Drc = 0;
                //copy(*op.roadMap, *RoadWidthDX);
            }}
        }
        if (SwitchHouses)
            copy(*op.houseMap, *HouseCover);

        if(SwitchHardsurface)
            copy(*op.hardsurfaceMap,*HardSurface);

        if(SwitchFlowBarriers)
        {
//            Fill(*tma,0.0);
//            FOR_ROW_COL_MV {
//                tma->Drc = std::max(std::max(std::max(FlowBarrierN->Drc,FlowBarrierE->Drc),FlowBarrierW->Drc),FlowBarrierS->Drc);
//            }
    //        copy(*op.flowbarriersMap,*tma);
        }
    }
    // MAP DISPLAY VARIABLES
    if(InfilMethod != INFIL_SWATRE && InfilMethod !=INFIL_NONE)
        avgTheta();
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
        SwitchWriteHeaders = false;
        QFile fout(newname1);
        fout.open(QIODevice::WriteOnly | QIODevice::Text);
        QTextStream out(&fout);

        out << "LISEM run - " << op.runfilename << "\n";
        out << "Time(min)" << sep << "P(mm)" << sep << "Ic(mm)";
        if (SwitchLitter)
            out << sep << "Ic(litter)(mm)";
        if (SwitchHouses)
            out << sep << "Ic(roof)(mm)";
        out << sep << "SS(mm)";
        if(SwitchIncludeET)
            out << sep << "ETa(mm)";
        out << sep << "Inf(mm)";
        out << sep << "Theta1 (-)";
        if (SwitchTwoLayer)
            out << sep << "Theta2 (-)";
        if (SwitchChannelBaseflow) {
            out << sep << "GWlevel (m)";
            out << sep << "Baseflow in (mm)";
        }
        if (SwitchIncludeStormDrains)
            out << sep << "StormDrain(mm)";
        out << sep << "Runoff(mm)";
        out << sep << "Flood(mm)";
        out << sep << "Flood area(m2)";
        out << sep << "Channels(mm)";
        out << sep << "Outflow(mm)";
        out << sep << "Boundary Outflow(mm)";
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
    out << sep << op.SurfStormm;
    if(SwitchIncludeET)
        out << sep << op.ETaTotmm;
    out << sep << op.InfilTotmm;
    out << sep << op.Theta1;
    if (SwitchTwoLayer)
        out << sep << op.Theta2;
    if (SwitchChannelBaseflow) {
        out << sep << op.GWlevel;
        out << sep << op.BaseFlowTotmm;
    }
    if (SwitchIncludeStormDrains)
        out << sep << op.StormDrainTotmm;
    out << sep << op.WaterVolTotmm;
    out << sep << op.volFloodmm;
    out << sep << op.FloodArea;
    out << sep << op.ChannelVolTotmm;
    out << sep << op.Qtotmm;
    out << sep << op.Qboundtotmm;
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

    double QALL = Qtot_dt * QUNIT/_dt; // total outflow for all outlets, same as point 0 in interface
    double QSALL = SoilLossTot_dt/_dt; //total sed loss in kg/s from all outlets, surface and boundary

    QFileInfo fi(resultDir + outflowFileName);

    //######  open files and write headers #####//

    QString unitS = "l/s";
    if (QUnits == 1)
        unitS = "m3/s";

    //PCRaster and flat format are mutually exclusive
    if (SwitchWriteHeaders) //  make file at first timestep
    {
        if (SwitchSeparateOutput) // each point separate file
        {
            FOR_ROW_COL_MV_OUTL {
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
                    //if (SwitchChannelBaseflow) nrs++;
                    if (SwitchIncludeTile) nrs++;
                    pnr.setNum(nrs);
                    out << pnr << "\n";

                    out << "run step\n";
                    if (SwitchRainfall) out << "Pavg (mm/h)\n";
                    if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                    out << "Qall " << unitS << "\n";
                    out << "QBound " << unitS << "\n";
                    if (SwitchIncludeChannel) {
                        out << "Q " << unitS << "\n" << "chanWH (m)\n";
                        // if (SwitchChannelBaseflow) out << ",Qb (l/s)";
                    } else {
                        out << "Q " << unitS << "\n";
                    }
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
                    out << ",Qbound";
                    if (SwitchIncludeChannel) {out << ",Q" << ",chanWH";}
                    else {out << ",Qrunoff";}
                   // if (SwitchChannelBaseflow) out << ",Qb";
                    if (SwitchIncludeTile) out << ",Qtile";
                    if (SwitchErosion) out << ",Qsall" << ",Qs" << ",C";
                    out << "\n";

                    out << "min";
                    if (SwitchRainfall) out << ",mm/h";
                    if (SwitchSnowmelt) out << ",mm/h";
                    out << "," << unitS;
                    out << "," << unitS;
                    if (SwitchIncludeChannel) out  << "," << unitS << ",m";
                  //  if (SwitchChannelBaseflowout << "," << unitS;
                    if (SwitchIncludeTile) out << "," << unitS;
                    if (SwitchErosion) out << ",kg/s"<< ",kg/s" << ",g/l";
                    out << "\n";
                }
                fout.close();
            }}
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
//                FOR_ROW_COL_MV
//                        if ( PointMap->Drc > 0 ) nr++;
                nr = crout_.size();

                // nr columns is time + rain (+ maybe snow) + nr points*(Q + Qs + C)
                int nrs = 1; // runstep
                if (SwitchRainfall) nrs++;
                if (SwitchSnowmelt) nrs++;
                nrs += (1 + 2*nr); // all water points
                if (SwitchIncludeTile) nrs+=nr;
                if (SwitchIncludeTile) nrs+=nr; //one tile per obs point
                if (SwitchErosion)
                    nrs += (1 + 2*nr); // all sed points

                pnr.setNum(nrs);  // nr of columnsin file

                out << "#LISEM total flow and sed output file for all reporting points\n";
                out <<  pnr << "\n";
                out << "Time (min)\n";
                if (SwitchRainfall) out << "Pavg (mm/h)\n";
                if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                out << "Qall (l/s)\n";
                out << "Qbound (l/s)\n";
                FOR_ROW_COL_MV_OUTL {
                    pnr.setNum(crout_[i_].nr);//(int)PointMap->Drc);
                    out << "Q #" << pnr <<  "(l/s)\n";
                    if (SwitchIncludeChannel) out << "chanWH #" << pnr <<  "(m)\n";
                 //   if (SwitchChannelBaseflow) out << "chanQb #" << pnr <<  "(l/s)\n";
                    if (SwitchIncludeTile) out << "Qtile #" << pnr <<  "(l/s)\n";
                }}
                if (SwitchErosion)
                {
                    out << "Qsall (kg/s)\n";
                    FOR_ROW_COL_MV_OUTL {
                        pnr.setNum(crout_[i_].nr);//(int)PointMap->Drc);
                        out << "Qs #"<< pnr << "(kg/s)\n" << pnr << "C #"<< pnr << "(g/l)\n";
                    }}
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
                out << ",Qbound";
                // for all points in outpoint.map
                FOR_ROW_COL_MV_OUTL {
                    pnr.setNum(crout_[i_].nr);//(int)PointMap->Drc);
                    out << ",Q #" << pnr;
                    if (SwitchIncludeChannel) out << ",chanWH #" << pnr;
                    if (SwitchIncludeTile) out << ",Qtile #" << pnr;
                }}
                // sediment
                if (SwitchErosion)
                {
                    out << ",Qsall";
                    FOR_ROW_COL_MV_OUTL {
                        pnr.setNum(crout_[i_].nr);//(int)PointMap->Drc);
                        out << ",Qs #" << pnr << ",C #" << pnr;
                    }}
                }
                out << "\n";

                //line 2 units
                out << "min";
                if (SwitchRainfall) out << ",mm/h";
                if (SwitchSnowmelt) out << ",mm/h";
                out << "," << unitS; // qall
                out << "," << unitS; // bound
                // for all outlets
                FOR_ROW_COL_MV_OUTL {
                    pnr.setNum(crout_[i_].nr);
                    out << "," << unitS << " #" << pnr;
                    if (SwitchIncludeChannel) out << ",m #" << pnr;
                    if (SwitchIncludeTile) out << "," << unitS << " #" << pnr;
                }}
                if (SwitchErosion)
                {
                    out << ",kg/s";
                    FOR_ROW_COL_MV_OUTL {
                        pnr.setNum(crout_[i_].nr);
                        out << ",kg/s #" << pnr << ",g/l #" << pnr;
                    }}
                }
                out << "\n";
            }
            fout.close();
        } // all in one
    }  // opening files and writing header

    //######  open files and append values #####//
    if (SwitchSeparateOutput)
    {
        // for all outlet points
        FOR_ROW_COL_MV_OUTL
        {
            newname1 = fi.path() + "/" + fi.baseName() + "_" + QString::number((int)PointMap->Drc) + "." +  fi.suffix();

            QFile fout(newname1);
            fout.open(QIODevice::Append | QIODevice::Text);

            QTextStream out(&fout);
            out.setFieldWidth(width);
            out.setRealNumberNotation(QTextStream::FixedNotation);
            out.setRealNumberPrecision(5);
            if (SwitchWritePCRtimeplot)
                out << runstep;
            else
                out << (time/60)/1440.0;

            out.setRealNumberPrecision(DIG);
            if (SwitchRainfall) out << sep << RainIntavg;
            if (SwitchSnowmelt) out << sep << SnowIntavg;

            out << sep << QALL << sep << BoundaryQ*QUNIT;

            if (SwitchIncludeChannel) {
                out << sep << ChannelQn->Drc*QUNIT;
                out << sep << ChannelWH->Drc;
                //if (SwitchChannelBaseflow) out << sep << Qbase->Drc*QUNIT;
            } else {
                out << sep << Qn->Drc*QUNIT;
            }

            if (SwitchIncludeTile) out << sep << TileQn->Drc*QUNIT;

            if (SwitchErosion) {
                out << sep << QSALL;
                if (SwitchIncludeChannel) {
                    out << sep << ChannelQsn->Drc;
                    out << sep << ChannelConc->Drc;
                } else {
                    out << sep << Qsn->Drc;
                    out << sep << TotalConc->Drc ;
                }
            }
            out << "\n";
            fout.close();
        }}
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
            out << (time/60)/1440.0;

        if (SwitchRainfall) out << sep << RainIntavg;
        if (SwitchSnowmelt) out << sep << SnowIntavg;

        out << sep << QALL << sep << BoundaryQ*QUNIT;

        FOR_ROW_COL_MV_OUTL
        {
            if (SwitchIncludeChannel) {
                out << sep << ChannelQn->Drc*QUNIT;
                out << sep << ChannelWH->Drc;
                //if (SwitchChannelBaseflow) out << sep << Qbase->Drc*QUNIT;
            } else {
                out << sep << Qn->Drc*QUNIT;
            }
            if (SwitchIncludeTile) out << sep << TileQn->Drc*QUNIT;
        }}

        if (SwitchErosion)
        {
            out << sep << QSALL;
            FOR_ROW_COL_MV_OUTL
            {
                if (SwitchIncludeChannel) {
                    out << sep << ChannelQsn->Drc;
                    out << sep << ChannelConc->Drc;
                } else {
                    out << sep << Qsn->Drc;
                    out << sep << TotalConc->Drc ;
                }
            }}
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
    out << "\"LISEM results at time (day:min):\"," << trunc(op.time/1440) << ":" << long(op.time) % 1440 <<"\n";
    if (op.CatchmentArea > 1e6)
        out << "\"Catchment area (km2):\"," << op.CatchmentArea/1e6<< "\n";
    else
        out << "\"Catchment area (m2):\"," << op.CatchmentArea<< "\n";
    out << "\"Total Precipitation (mm):\"," << op.RainTotmm<< "\n";
    out << "\"Total interception(mm):\"," << op.IntercTotmm<< "\n";
    out << "\"Total Litter interception (mm):\"," << op.IntercLitterTotmm<< "\n";
    out << "\"Total House interception (mm):\"," << op.IntercHouseTotmm<< "\n";
    out << "\"Surface storage (mm):\"," << op.SurfStormm<< "\n";
    out << "\"Total infiltration (mm):\"," << op.InfilTotmm<< "\n";
    out << "\"Total ETa (mm):\"," << op.ETaTotmm<< "\n";

    out << "\"Storm Drain (mm):\"," << op.StormDrainTotmm<< "\n";
    if (SwitchKinematic2D == K2D_METHOD_KIN) {
        out << "\"Water in overland flow (mm):\"," << op.WaterVolTotmm<< "\n";
        out << "\"Water in flood (mm):\"," << 0.0 << "\n";
    } else {
       out << QString("\"Water in overland flow (h<%1)(mm)):\",%2\n").arg(minReportFloodHeight*1000).arg(op.WaterVolTotmm);
       out << QString("\"Water in flood (h>%1) (mm)):\",%2\n").arg(minReportFloodHeight*1000).arg(op.volFloodmm);
    }
    out << "\"Water in channels (mm):\"," << op.ChannelVolTotmm<< "\n";
    out << "\"Water across boundary (mm):\"," << op.Qboundtotmm<< "\n";
    out << "\"Total baseflow inflow (mm):\"," << op.BaseFlowTotmm << "\n"; // not used!
    out << "\"Total outflow (all flows) (mm):\"," << op.Qtotmm+op.Qboundtotmm << "\n";
    out << "\n";
    out << "\"Total outflow (overland+channel+drains) (m3):\"," << op.Qtot<< "\n";
    out << "\"Total boundary outflow (m3):\"," << op.floodBoundaryTot<< "\n";
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
    for(int i = 1; i< op.OutletQpeak.length();i++)
    {
        out << "\"Peak discharge for outlet " + QString::number(i) +" (l/s):\"," << op.OutletQpeak.at(i)<< "\n";
    }
    for(int i = 1; i< op.OutletQpeak.length();i++)
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


    if (SwitchOutTheta) {
        if (InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE) {
            report(*ThetaI1a, OutTheta1);
            if (SwitchTwoLayer)
                report(*ThetaI2a, OutTheta2);
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
    name = resultDir + totalLandunitFileName;//QFileInfo(totalLandunitFileName).baseName()+"-"+op.timeStartRun+".csv";
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
        floodList[i].var0 = 0.05*i; //depth 5 cm intervals
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
                floodList[i].var6 += RoadWidthDX->Drc*DX->Drc; // WRONG: all road pixels is the surface, not the length
        }
    }}

    QFile fp(resultDir + floodStatsFileName);
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
    out << "\"results at time (day:min):\"" << op.time/1440 << ":" << long(op.time) % 1440 <<"\n";
    // "\"results at time (min):\"" << op.time << "\n";
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
  // op.OutletQb.append(new QVector<double>);
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
           // op.OutletQb.append(new QVector<double>);
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
//---------------------------------------------------------------------------
void TWorld::ClearHydrographData()
{
//    for(int i =op.OutletIndices.length() - 1; i >-1 ; i--)
//    {
//        delete op.OutletQ.at(i);
//        delete op.OutletQs.at(i);
//        delete op.OutletC.at(i);
//        delete op.OutletChannelWH.at(i);
//    }

    op.OutletIndices.clear();
    op.OutletLocationX.clear();
    op.OutletLocationY.clear();
    op.OutletQ.clear();
  //  op.OutletQb.clear();
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
// include all totals!!!
// trial, dump all values to restart run where it crashed
void TWorld::ReportDump(void)
{
    dumpDir = inputDir + "/Dump/";
    if (!QFileInfo(dumpDir).exists()){
        QDir().mkdir(dumpDir);
    }

    report(*RainCumFlat,dumpDir+"raincum.map"); //!
  //  report(*Rainpeak,dumpDir+"rainpeak.map");

    if (SwitchIncludeET) {
        report(*ETaCum,dumpDir+"ETaCum.map");
        report(*IntercETa,dumpDir+"IntercETa.map");
    }

    report(*Interc,dumpDir+"Interc.map");
    report(*CStor,dumpDir+"CStor.map");
    if (SwitchLitter)
        report(*LInterc,dumpDir+"LInterc.map");
    if (SwitchHouses)
        report(*IntercHouse,dumpDir+"IntercHouse.map");

    if(InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE) {
        report(*ThetaI1,dumpDir+"ThetaI1.map");
        if (SwitchTwoLayer)
            report(*ThetaI2,dumpDir+"ThetaI2.map");

       report(*InfilVol,dumpDir+"InfilVol.map");
    }

    report(*WHstore,dumpDir+"WHstore.map");
    report(*WH,dumpDir+"WH.map");
    report(*WHrunoff,dumpDir+"WHrunoff.map");
    report(*WHroad,dumpDir+"WHroad.map");

    report(*WaterVolall,dumpDir+"WaterVolall.map");
    report(*FloodWaterVol,dumpDir+"FloodWaterVol.map");


    if (SwitchIncludeChannel)
    {
        report(*ChannelWaterVol,dumpDir+"ChannelWaterVol.map");
        report(*ChannelWH,dumpDir+"ChannelWH.map");
        if (SwitchChannelBaseflow) {
            report(*Qbase,dumpDir+"Qbase.map");
            report(*GWWH,dumpDir+"GWWH.map");
        }
    }





}
