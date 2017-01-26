
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
    ReportTimeseriesNew();
    // report hydrographs ande sedigraphs at all points in outpoint.map

    ReportTotalsNew();
    // report totals to a text file

    ReportMaps();
    // report all maps and mapseries

    ReportLandunits();
    // reportc stats per landunit class

    ChannelFloodStatistics();
    // report buildings submerged in flood level classes
}
//---------------------------------------------------------------------------
/** fill output structure 'op' with results to talk to the interface:
    report to screen, hydrographs and maps */
void TWorld::OutputUI(void)
{

    //hydrographs
    op.timestep = this->_dt/60.0;
    op.OutletQ.at(0)->append((QtotT * 1000.0/_dt)); //QtotT is in m3
    op.OutletQs.at(0)->append(SoilLossTotT);
    op.OutletC.at(0)->append((QtotT) > MIN_FLUX? SoilLossTotT/(QtotT) : 0);
    op.OutletQtot.replace(0,Qtot);
    op.OutletQstot.replace(0,SoilLossTot/1000.0);

    double channelwh = 0;
    if(SwitchIncludeChannel)
    {

        FOR_ROW_COL_MV_CH
        {
            if(LDDChannel->Drc == 5)
            {
                channelwh += ChannelWH->Drc;
            }
        }
    }
    op.OutletChannelWH.at(0)->append(channelwh);
    for(int j = 1; j < op.OutletIndices.length(); j++)
    {
        int r = op.OutletLocationX.at(j);
        int c = op.OutletLocationY.at(j);

        double discharge = Qoutput->Drc; //sum of current Qn, ChannelQn, TileQn in l/s
        double sedimentdischarge = SwitchErosion? Qsoutput->Drc * _dt : 0.0;
        double sedimentconcentration = SwitchErosion? TotalConc->Drc : 0.0;
        double channelwh = SwitchIncludeChannel? ChannelWH->Drc : 0.0;

        op.OutletQtot.replace(j,op.OutletQtot.at(j) + _dt * discharge/1000.0); //cumulative in m3/s
        op.OutletQstot.replace(j,op.OutletQstot.at(j) + _dt * sedimentdischarge/1000.0);
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

    if(SwitchKinematic2D != K1D_METHOD)
    {
        FOR_ROW_COL_MV
        {
            if(K2DOutlets->Drc ==1)
            {
                V->Drc = 0;
            }
        }
    } //why set velocity at zero, and is also done

    //display maps
//    fill(*COMBO_QOFCH, 0.0);
//    calcMap(*COMBO_QOFCH, *Qn, ADD);
//    if (SwitchIncludeChannel)
//        calcMap(*COMBO_QOFCH, *ChannelQn, ADD);
//    if(SwitchChannelFlood)
//        calcMap(*COMBO_QOFCH, *Qflood, ADD);

    if (SwitchIncludeChannel)
    {
        fill(*tma,0.0);
        DistributeOverExtendedChannel(ChannelQn,tma);
        // VJ 161222 must be channelqn not channelq
    }

    FOR_ROW_COL_MV
    {
        COMBO_QOFCH->Drc = Qn->Drc;
        if(SwitchChannelFlood)
            COMBO_QOFCH->Drc += Qflood->Drc;
        if (SwitchIncludeChannel)
            if (ChannelWidthExtended->Drc > 0)
                COMBO_QOFCH->Drc = tma->Drc;
    }
    //VJ changed this to Qn and channel Qn

    FOR_ROW_COL_MV
    {
        COMBO_VOFCH->Drc = V->Drc;
        if(SwitchChannelFlood)
            COMBO_VOFCH->Drc += UVflood->Drc;
        if (SwitchIncludeChannel)
            if (ChannelWidthUpDX->Drc > 0)
                COMBO_VOFCH->Drc = ChannelV->Drc;
    }



    if(SwitchErosion)
    {
        fill(*COMBO_SS, 0.0);
        if(SwitchChannelFlood)
        {
            calcMap(*COMBO_SS, *BLFlood, ADD);
            calcMap(*COMBO_SS, *SSFlood, ADD);
        }
        if(SwitchIncludeChannel)
        {
            calcMap(*COMBO_SS, *ChannelBLSed, ADD);
            calcMap(*COMBO_SS, *ChannelSSSed, ADD);
        }
        calcMap(*COMBO_SS, *Sed, ADD);
    }

    //output maps for combo box
    for(int i = 0; i < op.ComboMapsSafe.length(); i++)
    {
        fill(*tma, 0.0);
        calcMapValue(*tma, *op.ComboMaps.at(i),op.ComboScaling.at(i), MUL);
        copy(*op.ComboMapsSafe.at(i), *tma);
    }

    //make sure sediment maps for all grain sizes are present
    if(SwitchErosion && SwitchUseGrainSizeDistribution)
    {
        FOR_GRAIN_CLASSES
        {
            if(op.graindiameters.length() < numgrainclasses + 1)
            {
                op.graindiameters.append(graindiameters.at(d));
            }else
            {
                break;
            }

        }
    }

    copy(*op.baseMap, *Shade);
    copy(*op.baseMapDEM, *DEM);

    if (SwitchIncludeChannel)
        copy(*op.channelMap, *ChannelWidthExtended);
    //BB 151118 might be better to draw LDD, since that is actually used to determine the presence of a channel
report(*ChannelWidthExtended,"cwe.map");
    if (SwitchRoadsystem)
    {
        copy(*op.roadMap, *RoadWidthDX);
       // calcMap(*op.roadMap, *HardSurface, ADD);
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

    // MAP DISPLAY VARIABLES

    op.dx = _dx;
    op.MB = MB;
    op.runstep = runstep;
    op.maxstep = (int) ((EndTime-BeginTime)/_dt);
    op.EndTime = EndTime/60.0;
    //	op.BeginTime = BeginTime/60.0;

    op.CatchmentArea = CatchmentArea;

    op.BaseFlowtotmm = BaseFlow * 1000.0/(_dx*_dx*nrCells);
    op.LitterStorageTotmm = mapTotal(*LInterc) *1000.0/(_dx*_dx*nrCells);
    op.ChannelVolTotmm = SwitchIncludeChannel? mapTotal(*ChannelWaterVol) * 1000.0/(_dx*_dx*nrCells) : 0.0;

    op.RainTotmm = RainTotmm + SnowTotmm;
    op.WaterVolTotmm = WaterVolRunoffmm;//WaterVolTotmm-SurfStoremm;
    op.Qtotmm = Qtotmm;
    op.Qtot = QtotOutlet;

    op.Qtile = 1000*TileQn->DrcOutlet;

    op.RainpeakTime = RainpeakTime/60;
    op.FloodTotMax = floodVolTotMax;
    op.FloodAreaMax = floodAreaMax;

    op.InfilTotmm = InfilTotmm;
    op.SurfStormm = SurfStoremm;
    op.IntercTotmm = IntercTotmm;// + IntercHouseTotmm;
    //houses
    op.IntercHouseTotmm = IntercHouseTotmm;
    op.InfilKWTotmm = InfilKWTot; // infil part in kin wave not used
    op.RunoffFraction = 0;
    if (op.RainTotmm > 0)
        op.RunoffFraction = std::max(0.0, (op.Qtotmm - op.BaseFlowtotmm)/op.RainTotmm);
    op.MBs = MBs;
    op.DetTotSplash=DetSplashTot*0.001; // convert from kg to ton per cell
    op.DetTotFlow=DetFlowTot*0.001; // convert from kg to ton
    op.DepTot=DepTot*0.001; // convert from kg to ton
    op.SedTot=SedTot*0.001; // convert from kg to ton

    op.ChannelDetTot=ChannelDetTot*0.001; // convert from kg to ton
    op.ChannelDepTot=ChannelDepTot*0.001; // convert from kg to ton
    op.ChannelSedTot=ChannelSedTot*0.001; // convert from kg to ton

    op.FloodDepTot = FloodDepTot*0.001;
    op.FloodDetTot = FloodDetTot*0.001;
    op.FloodSedTot = FloodSedTot*0.001;

    op.SoilLossTot=SoilLossTot/*Outlet*/ *0.001; // convert from kg to ton

    op.t = time_ms.elapsed()*0.001/60.0;
    op.time = time/60;
    op.maxtime = op.t/runstep * op.maxstep;

    op.Pmm = (RainAvgmm + SnowAvgmm)*3600/_dt;
    op.volFloodmm = floodTotmm;

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
    int hour = 0;
    int min = 0;
    int sec = 0;
    int DIG = ReportDigitsOut;
    //int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
    double RainIntavg = RainAvgmm * 3600/_dt;
    double SnowIntavg = SnowAvgmm * 3600/_dt;
    QString newname1, pnr, sep = (SwitchWritePCRtimeplot ? " " : ",");
    int width = (!SwitchWritePCRtimeplot ? 0 : 9+DIG-3);
    // NOTE if SwitchWriteCommaDelimited = true then SwitchWritePCRtimeplot = false

    double QALL = QtotT * 1000.0/_dt; // total outflow on all sides in l/s, same as point 0 in interface

    QFileInfo fi(resultDir + outflowFileName);

    if (SwitchSOBEKoutput)
        SwitchSeparateOutput = true;

    // - TODO: if runs are interrupted the number of lines win the SOBEK output will not be correct!
    if (SwitchSOBEKoutput && time >= 0)
    {
        sec = (int)time; //(time*60);
        hour = (int)(sec/3600);
        min = (int)(time/60 - hour*60);
        sec = (int)(sec - hour * 3600 - min * 60);
    }


    //######  open files and write headers #####//

    //SOBEK, PCRaster and flat format are mutually exclusive
    if (SwitchWriteHeaders) //  make file at first timestep
    {
        SwitchWriteHeaders = false;
        if (SwitchSeparateOutput) // each point separate file
        {
            FOR_ROW_COL_MV
            {
                if ( PointMap->Drc > 0 )
                {
                    newname1 = fi.path() + "/" + fi.baseName() + "_" +
                            QString::number((int)PointMap->Drc) + "." +  fi.suffix();
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
                        out << "#LISEM total flow and sed output file for point " << pnr << "\n";

                        // nr columns is time + rain (+ maybe snow) + Q + (maybe Qs + C)
                        int nrs = 3 + (SwitchErosion ? 2 : 0);
                        if (SwitchSnowmelt && SwitchRainfall) nrs++;
                        pnr.setNum(nrs);
                        out << pnr << "\n";

                        out << "run step\n";
                        if (SwitchRainfall) out << "Pavg (mm/h)\n";
                        if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                        out << "Qall (l/s)\n" << "Qoutlet (l/s)\n" << "chanWH (m)\n";
                        if (SwitchIncludeTile) out << "Qdrain (l/s)\n";
                        if (SwitchErosion) out << "Qs (kg/s)\n";
                        if (SwitchErosion) out << "C (g/l)\n";
                    }
                    else // SOBEK format
                        if (SwitchSOBEKoutput)
                        {
                            //NO HEADER
                            //                     pnr.setNum(SOBEKlines);
                            //                     out << "Q";
                            //                     if (SwitchErosion) out << " Qs C";
                            //                     out << "\n";
                            //                     out << "m3/s";
                            //                     if (SwitchErosion) out << " kg/s g/l";
                            //                     out << "\n";
                            //                     out << pnr << "\n";
                            out << "* " <<  (int)PointMap->Drc << "\n";
                        }
                        else // flat format, comma delimited
                        {
                            pnr.setNum((int)PointMap->Drc);
                            out << "LISEM total flow and sed output file for point " << pnr << "\n";
                            out << "Time";
                            if (SwitchRainfall) out << ",Pavg";
                            if (SwitchSnowmelt) out << ",Snowavg";
                            out << ",Q" << ",chanWH";
                            if (SwitchErosion) out << ",Qs,C";
                            out << "\n";
                            out << "min,mm/h";
                            if (SwitchSnowmelt) out << ",mm/h";
                            out << ",l/s" << ",l/s" << ",m";
                            if (SwitchIncludeTile) out << ",l/s";
                            if (SwitchErosion) out << ",kg/s,g/l";
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
                int nrs = 2+(1+(SwitchErosion ? 2 : 0))*nr;
                if (SwitchSnowmelt && SwitchRainfall) nrs++;
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
                    out << "chanWH #" << pnr <<  "(m)\n";
                    if (SwitchIncludeTile) out << "Qdr #" << pnr <<  "(l/s)\n";
                    if (SwitchErosion) out << "Qs #"<< pnr << "(kg/s)\n";
                    if (SwitchErosion) out << "C #"<< pnr << "(g/l)\n";
                }
                out << "\n";
            }
            // combination cannot occur: SOBEK is always separate output
            else // SOBEK format
                if (SwitchSOBEKoutput) //note: sobek input does not have rainfall
                {
                    FOR_ROW_COL_MV
                            if ( PointMap->Drc > 0 )
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << " Q #" << pnr;
                        if (SwitchErosion) out << " Qs #" << pnr;
                        if (SwitchErosion) out << " C #" << pnr;
                    }
                    out << "\n";
                    FOR_ROW_COL_MV
                            if ( PointMap->Drc > 0 )
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << " m3/s #" << pnr;
                        if (SwitchErosion) out << " kg/s #" << pnr;
                        if (SwitchErosion) out << " g/l #" << pnr;
                    }
                    out << "\n";
                    out << SOBEKdatestring << "\n";
                }
                else // flat CSV format, comma delimited
                {
                    out << "#LISEM total flow and sed output file for all reporting points in map\n";
                    out << "Time";
                    if (SwitchRainfall) out << ",P";
                    if (SwitchSnowmelt) out << ",Snow";
                    out << ",Qall";
                    FOR_ROW_COL_MV
                            if ( PointMap->Drc > 0 )
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << ",Q #" << pnr;
                        out << ",chanWH #" << pnr;
                        if (SwitchIncludeTile) out << ",Qdr #" << pnr;
                        if (SwitchErosion) out << ",Qs #" << pnr;
                        if (SwitchErosion) out << ",C #" << pnr;
                    }
                    out << "\n";
                    out << "min";
                    if (SwitchRainfall) out << ",mm/h";
                    if (SwitchSnowmelt) out << ",mm/h";
                    out << ",l/s";
                    FOR_ROW_COL_MV
                            if ( PointMap->Drc > 0 )
                    {
                        pnr.setNum((int)PointMap->Drc);
                        out << ",l/s #" << pnr;
                        out << ",m #" << pnr;
                        if (SwitchIncludeTile) out << ",l/s #" << pnr;
                        if (SwitchErosion) out << ",kg/s #" << pnr;
                        if (SwitchErosion) out << ",g/l #" << pnr;
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
           //     qDebug() << PointMap->Drc << r << c;
                newname1 = fi.path() + "/" + fi.baseName() + "_" +
                        QString::number((int)PointMap->Drc) + "." +  fi.suffix();
                QFile fout(newname1);
                fout.open(QIODevice::Append | QIODevice::Text);

                QTextStream out(&fout);
                out.setRealNumberPrecision(DIG);
                out.setFieldWidth(width);
                out.setRealNumberNotation(QTextStream::FixedNotation);

                if (!SwitchSOBEKoutput)   //PCRaster timeplot and flat CSV format
                {

                    if (SwitchWritePCRtimeplot)
                        out << runstep;
                    else
                        out << time/60;
                    if (SwitchRainfall) out << sep << RainIntavg;
                    if (SwitchSnowmelt) out << sep << SnowIntavg;
                    out << sep << QALL << Qoutput->Drc << sep << ChannelWH->Drc;
                    if (SwitchIncludeTile) out << sep << TileQn->Drc*1000;
                    if (SwitchErosion) out << sep << Qsoutput->Drc;
                    if (SwitchErosion) out << sep << TotalConc->Drc;
                    out << "\n";
                }
                else  //SOBEK format
                {
                    QString ss = QString("\"%1;%2:%3:%4\" %5 <\n").
                            arg(SOBEKdatestring).
                            arg((uint)hour,2,10,QLatin1Char('0')).
                            arg((uint)min,2,10,QLatin1Char('0')).
                            arg((uint)sec,2,10,QLatin1Char('0')).
                            arg(Qoutput->Drc/1000.0,0,'f',3);
                    out << ss;
                }
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

        if (!SwitchSOBEKoutput)
        {
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

                    out << sep << Qoutput->Drc << sep << ChannelWH->Drc;
                    if (SwitchIncludeTile) out << sep << TileQn->Drc*1000;
                    if (SwitchErosion) out << sep << Qsoutput->Drc;
                    if (SwitchErosion) out << sep << TotalConc->Drc;
                }
            }
            out << "\n";
        }
        else
            if (SwitchSOBEKoutput)
            {
                out.setFieldWidth(2);
                out << "\"" << SOBEKdatestring << ":" << hour << ":" <<  min << ":" <<  sec;
                out.setFieldWidth(8);
                FOR_ROW_COL_MV
                {
                    if ( PointMap->Drc > 0 )
                    {
                        out << " " << Qoutput->Drc/1000.0;
                        if (SwitchErosion) out << " " << Qsoutput->Drc;
                        if (SwitchErosion) out << " " << TotalConc->Drc;
                    }
                }
                out << " < \n";
            }
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
    out << "\"LISEM run with:," << op.runfilename << "\"\n";
    out << "\"LISEM results at time (min):," << op.time <<"\"\n";
    out << "\"Catchment area              (ha):\"," << op.CatchmentArea/10000.0<< "\n";
    out << "\"Total Precipitation         (mm):\"," << op.RainTotmm<< "\n";
    out << "\"Total discharge             (mm):\"," << op.Qtotmm<< "\n";
    out << "\"Total interception          (mm):\"," << op.IntercTotmm<< "\n";
    out << "\"Total House interception    (mm):\"," << op.IntercHouseTotmm<< "\n";
    out << "\"Total infiltration          (mm):\"," << op.InfilTotmm<< "\n";
    out << "\"Surface storage             (mm):\"," << op.SurfStormm<< "\n";
    out << "\"Water in runoff + channel   (mm):\"," << op.WaterVolTotmm<< "\n";
    out << "\"Total discharge             (m3):\"," << op.Qtot<< "\n";
    out << "\"Peak time precipitation    (min):\"," << op.RainpeakTime<< "\n";
    out << "\"Peak discharge/Precipitation (%):\"," << op.RunoffFraction*100<< "\n";
    out << "\"Flood volume (max level)    (m3):\"," << op.FloodTotMax<< "\n";
    out << "\"Flood area (max level)      (m2):\"," << op.FloodAreaMax<< "\n";
    out << "\"Splash detachment (land)   (ton):\"," << op.DetTotSplash<< "\n";
    out << "\"Flow detachment (land)     (ton):\"," << op.DetTotFlow<< "\n";
    out << "\"Deposition (land)          (ton):\"," << op.DepTot<< "\n";
    out << "\"Suspended Sediment (land)  (ton):\"," << op.SedTot<< "\n";
    out << "\"Flow detachment (channels) (ton):\"," << op.ChannelDetTot<< "\n";
    out << "\"Deposition (channels)      (ton):\"," << op.ChannelDepTot<< "\n";
    out << "\"Susp. Sediment (channels)  (ton):\"," << op.ChannelSedTot<< "\n";
    out << "\"Flow detachment (flood)    (ton):\"," << op.FloodDetTot<< "\n";
    out << "\"Deposition (flood)         (ton):\"," << op.FloodDepTot<< "\n";
    out << "\"Susp. Sediment (flood   )  (ton):\"," << op.FloodSedTot<< "\n";
    out << "\"Total soil loss            (ton):\"," << op.SoilLossTot<< "\n";
    out << "\"Average soil loss        (kg/ha):\"," << (op.SoilLossTot*1000.0)/(op.CatchmentArea/10000.0)<< "\n";
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
void TWorld::ReportMaps(void)
{
    FOR_ROW_COL_MV
    {
        tm->Drc = (RainCumFlat->Drc + SnowmeltCum->Drc*DX->Drc/_dx) * 1000.0; // m to mm
        tma->Drc = (Interc->Drc + IntercHouse->Drc)*1000.0/CellArea->Drc;
    }

    report(*tm, rainfallMapFileName);
    report(*tma, interceptionMapFileName);

    report(*InfilmmCum, infiltrationMapFileName);

    report(*runoffTotalCell, runoffMapFileName); // in m3, total runoff from cell (but there is also runon!)
    report(*runoffFractionCell, runoffFractionMapFileName);

    report(*WHmax, floodWHmaxFileName);

    if (SwitchIncludeChannel)
    {
        report(*ChannelQntot, channelDischargeMapFileName);
    }

    if(SwitchErosion)
    {
        // deal with erosion units, 0 = ton/ha, 1 = kg/m2, 2 = kg/cell
        if (ErosionUnits == 0)
            fill(*tma,10.0);
        else
            fill(*tma,1.0);
        if (ErosionUnits == 2 || ErosionUnits == 0)
            calcMap(*tma, *CellArea, DIV);


        // all detachment combined
        FOR_ROW_COL_MV
        {
            tm->Drc =std::max(0.0,TotalDetMap->Drc + TotalDepMap->Drc);
        }
        calcMap(*tm, *tma, MUL);
        report(*tm, totalErosionFileName);
        if (outputcheck[5].toInt() == 1)
            report(*tm, Outeros); // in units

        // all deposition combined
        FOR_ROW_COL_MV
        {
            tm->Drc =std::min(0.0,TotalDetMap->Drc + TotalDepMap->Drc);
        }
        calcMap(*tm, *tma, MUL);
        report(*tm, totalDepositionFileName);
        if (outputcheck[6].toInt() == 1)
            report(*tm, Outdepo); // in units

        // all channel depostion combined
        if (SwitchIncludeChannel)
        {
            FOR_ROW_COL_MV_CH
            {
                tm->Drc =std::min(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc);
            }
            calcMap(*tm, *tma, MUL);
            report(*tm, totalChanDepositionFileName);

            // all channel detachment combined
            FOR_ROW_COL_MV_CH
            {
                tm->Drc =std::max(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc);
            }
            calcMap(*tm, *tma, MUL);
            report(*tm, totalChanErosionFileName);
        }
        copy(*tm, *TotalSoillossMap); //kg/cell
        calcMap(*tm, *tma, MUL);
        report(*tm, totalSoillossFileName);
        if (outputcheck[16].toInt() == 1) report(*tm, OutSL);      // in user units

        // total sediment
        copy(*tm, *COMBO_SS); //kg/cell
        calcMap(*tm, *tma, MUL);
        if (outputcheck[17].toInt() == 1) report(*tm, OutSed);      // in user units

        if (outputcheck[1].toInt() == 1) report(*Conc, Outconc);  // in g/l
        if (outputcheck[4].toInt() == 1) report(*TC, Outtc);      // in g/l

    }

    if (SwitchChannelFlood)
    {
        report(*floodHmxMax, floodLevelFileName);
        report(*floodTime, floodTimeFileName);
        report(*floodTimeStart, floodFEWFileName);
        report(*maxChannelflow, floodMaxQFileName);
        report(*maxChannelWH, floodMaxChanWHFileName);
        report(*floodVMax, floodMaxVFileName);
    }

    if (outputcheck[0].toInt() == 1)
        report(*Qoutput, Outrunoff); // in l/s
    if (outputcheck[2].toInt() == 1)
    {
        calcMapValue(*tm, *WH, 1000, MUL);// WH in mm
        report(*tm, Outwh);
    }

    if (outputcheck[3].toInt() == 1)
        report(*runoffTotalCell, Outrwh); // in mm
    // changed to cum runoff in mm

    if (outputcheck[7].toInt() == 1)
        report(*V, Outvelo);

    if (outputcheck[8].toInt() == 1)
        report(*InfilmmCum, Outinf); // in mm

    if (outputcheck[9].toInt() == 1)
    {
        calcMapValue(*tm, *WHstore, 1000, MUL);// in mm
        report(*tm, Outss);
        /** TODO check this: surf store in volume m3 is multiplied by flowwidth? */
    }

    if (outputcheck[10].toInt() == 1) report(*ChannelWaterVol, Outchvol);


    if (SwitchIncludeTile && outputcheck.count() > 11)
    {
        if (outputcheck[11].toInt() == 1)
        {
            calcMapValue(*tm, *TileQn, 1000, MUL);// in mm
            report(*tm, OutTiledrain);
        }
    }
    if (SwitchChannelFlood)
    {
        if (outputcheck[12].toInt() == 1)
        {
            report(*hmx, OutHmx);
        }
        if (outputcheck[13].toInt() == 1)
        {
            report(*Qflood, OutQf);
        }
        if (outputcheck[14].toInt() == 1)
        {
            report(*UVflood, OutVf);
        }
        if (outputcheck[15].toInt() == 1)
        {
            report(*hmxWH, OutHmxWH);
        }

    }

}
//---------------------------------------------------------------------------
/// Land unit statistics: count nr land units in classifiedfile
// VJ 110110 count nr of land units in classified file
void TWorld::CountLandunits(void)
{
    int i, j;
    for (i = 0; i < 512; i++)
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

        if(!found && i < 512)
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
    QString newname1;

    if (!SwitchErosion)
        return;

    for (int i = 0; i < landUnitNr; i++)
    {
        unitList[i].var0 = 0;
        unitList[i].var1 = 0;
        unitList[i].var2 = 0;
        unitList[i].var3 = 0;
        unitList[i].var4 = 0;
        unitList[i].var5 = 0;
    }

    FOR_ROW_COL_MV
    {
        //variables are kg/cell convert to ton/cell
        for (int i = 0; i < landUnitNr; i++)
            if (unitList[i].nr == (long)LandUnit->Drc)
            {
                unitList[i].var0 += CellArea->Drc/10000;//ha
                unitList[i].var1 += TotalDetMap->Drc/1000; //ton/cell
                unitList[i].var2 += TotalDepMap->Drc/1000;
                unitList[i].var3 += TotalSoillossMap->Drc/1000;
            }
    }


    newname1 = resultDir + totalLandunitFileName;

    QFile fout(newname1);
    fout.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&fout);
    out.setRealNumberPrecision(3);
    out.setFieldWidth(12);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    //TODO: make this comma delimited?
    out << "    Landunit        Area  Detachment  Deposition   Soil Loss\n";
    out << "           #          ha         ton         ton         ton\n";
    for (int i = 0; i < landUnitNr; i++)
        out << unitList[i].nr
            << unitList[i].var0
            << unitList[i].var1
            << unitList[i].var2
            << unitList[i].var3
            << "\n";
    fout.close();

}
//---------------------------------------------------------------------------
void TWorld::ChannelFloodStatistics(void)
{
    if (!SwitchIncludeChannel)
        return;
    if (!SwitchChannelFlood)
        return;

    for (int i = 0; i < 256; i++)
    {
        floodList[i].nr = i;
        floodList[i].var0 = 0.1*i; //depth
        floodList[i].var1 = 0;
        floodList[i].var2 = 0;
        floodList[i].var3 = 0;
        floodList[i].var4 = 0;
        floodList[i].var5 = 0;
    }
    double area = _dx*_dx;
    int nr = 0;
    FOR_ROW_COL_MV
    {
        if(floodHmxMax->Drc > 0)//minReportFloodHeight)
        {
            int i = 0;
            while (floodList[i].var0 < floodHmxMax->Drc && i < 256)
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
        }
    }

    QString name;
    name = resultDir + QFileInfo(floodStatsFileName).baseName()+timestampRun+".csv";
    QFile fp(name);
    if (!fp.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&fp);
    out.setRealNumberPrecision(2);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    out << "\"LISEM run with:," << op.runfilename << "\"" << "\n";
    out << "\"results at time (min):\"" << op.time << "\n";
    out << "class,Depth,Area,Volume,Duration,Start,Structures\n";
    out << "#,m,m2,m3,h,h,m2\n";
    for (int i = 0; i < nr+1; i++)
        out << i << ","
            << floodList[i].var0 << ","
            << floodList[i].var1 << ","
            << floodList[i].var2 << ","
            << floodList[i].var3 << ","
            << floodList[i].var4 << ","
            << floodList[i].var5
            << "\n";

    double totarea = 0;
    double totvol = 0;
    double totbuild = 0;
    for (int i = 0; i < nr+1; i++)
    {
        totarea += floodList[i].var1;
        totvol += floodList[i].var2;
        totbuild += floodList[i].var5;
    }

    fp.flush();
    fp.close();

}//---------------------------------------------------------------------------
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
        delete op.roadMap;
        delete op.houseMap;
        delete op.flowbarriersMap;
    }

    op.baseMap = new cTMap();
    op.baseMapDEM = new cTMap();
    op.channelMap = new cTMap();
    op.roadMap = new cTMap();
    op.houseMap = new cTMap();
    op.flowbarriersMap = new cTMap();

    op.baseMap->MakeMap(LDD, 0);
    op.baseMapDEM->MakeMap(LDD, 0);
    op.channelMap->MakeMap(LDD, 0);
    op.roadMap->MakeMap(LDD, 0);
    op.houseMap->MakeMap(LDD, 0);
    op.flowbarriersMap->MakeMap(LDD, 0);
}
//---------------------------------------------------------------------------
void TWorld::setupHydrographData()
{
    // VJ 110630 show hydrograph for selected output point
    bool found = false;

    FOR_ROW_COL_MV
    {
        if(PointMap->Drc > 0)
        {
            found = true;
        }
    }
    if(!found)
    {
        ErrorString = QString("Outpoint.map has no values above 0");
        throw 1;
    }

    ClearHydrographData();


    //get the sorted locations and index numbers of the outlet points
    QList<int> nr;
    //int maxnr = 0;

    //0 is reserved for total outflow (channel and overland flow)
    nr.append(0);
    op.OutletIndices.append(0);
    op.OutletLocationX.append(0);
    op.OutletLocationY.append(0);
    op.OutletQ.append(new QList<double>);
    op.OutletQs.append(new QList<double>);
    op.OutletC.append(new QList<double>);
    op.OutletChannelWH.append(new QList<double>);
    op.OutletQpeak.append(0);
    op.OutletQpeaktime.append(0);
    op.OutletQtot.append(0);
    op.OutletQstot.append(0);

    FOR_ROW_COL_MV
    {
        if(PointMap->Drc > 0)
        {
           nr.append(PointMap->Drc);
           op.OutletIndices.append(PointMap->Drc);
           op.OutletLocationX.append(r);
           op.OutletLocationY.append(c);
           op.OutletQ.append(new QList<double>);
           op.OutletQs.append(new QList<double>);
           op.OutletC.append(new QList<double>);
           op.OutletChannelWH.append(new QList<double>);
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

    qSort(nr);
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
}

//---------------------------------------------------------------------------
void TWorld::GetComboMaps()
{
    //combo box maps
    ClearComboMaps();

    QList<double> Colormap;
    QList<QString> Colors;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.0005);
    Colormap.append(0.01);
    Colormap.append(0.05);
    Colormap.append(0.1);
    Colormap.append(0.9);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#8C8CFF");
    Colors.append("#8080FF");
    Colors.append("#4040ff");
    Colors.append("#0000FF");
    Colors.append("#00006F");
    Colors.append("#FF0000");
    Colors.append("#FF3300");
    AddComboMap(0,"Total Discharge","l/s",COMBO_QOFCH,Colormap,Colors,true,false,1000.0, 1.0);
  //  AddComboMap(0,"Overland Discharge","l/s",Qn,Colormap,Colors,true,false,1000.0, 1.0); //VJ changed to Qn instead of Q
    if(SwitchIncludeChannel)
        AddComboMap(0,"Channel Discharge","l/s",ChannelQn,Colormap,Colors,true,false,1000.0, 1.0); //Chnaged thhis to ChannelQn

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.25);
    Colormap.append(0.75);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#00FF00");
    Colors.append("#FFFF00");
    Colors.append("#FF0000");
    Colors.append("#A60000");
//    AddComboMap(0,"Overland Flow Velocity","m/s",V,Colormap,Colors,false,false,1.0, 0.01);
    AddComboMap(0,"Flow Velocity","m/s",COMBO_VOFCH,Colormap,Colors,false,false,1.0, 0.01);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"Channel Velocity","m/s",ChannelV,Colormap,Colors,false,false,1.0,0.01);
    }

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.5);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#7692FF");
    Colors.append("#0023b1");
    Colors.append("#001462");

    if(SwitchChannelFlood)
    {
        AddComboMap(0,"Water Height","m",hmxWH,Colormap,Colors,false,false,1.0,0.01);
        AddComboMap(0,"Flood Height","m",hmx,Colormap,Colors,false,false,1.0,0.01);
    }
    else
        AddComboMap(0,"Overland Flow Height","m",  /*K2DWHStore*/WHrunoff,Colormap,Colors,false,false,1.0, 0.01);

//    if(SwitchKinematic2D > 1)
//    {
//        AddComboMap(0,"Macro depression storage","m",K2DWHStore,Colormap,Colors,false,false,1.0, 0.01);
//    }

    if(InfilMethod != INFIL_NONE)
    {
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.1);
        Colormap.append(0.4);
        Colormap.append(1.0);

        Colors.clear();
        Colors.append("#f6f666");
        Colors.append("#FFFF55");
        Colors.append("#8080FF");
        Colors.append("#0000AA");
        AddComboMap(0,"Infiltration","mm",InfilmmCum,Colormap,Colors,false,false,1.0,1.0);
    }

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.5);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#8888FF");
    Colors.append("#0000FF");
    Colors.append("#008800");

    double factor = 3600000.0/_dt; //from m to mm/h

//    copy(*tm, *RainCumFlat );
//    calcMap(*tm,*SnowmeltCum, ADD);
//    copy(*tma, *Rain );
//    calcMap(*tma,*Snowmelt, ADD);
//    AddComboMap(0,"Precip. Cumulative","mm",tm,Colormap,Colors,false,false,1000.0,0.1);
//    AddComboMap(0,"Precip. Intensity","mm/h",tma,Colormap,Colors,false,false,factor,0.1);
    AddComboMap(0,"Rainfall Cumulative","mm",RainCumFlat,Colormap,Colors,false,false,1000.0,0.1);
    AddComboMap(0,"Rainfall Intensity","mm/h",Rain,Colormap,Colors,false,false,factor,0.1);

    if(SwitchChannelFlood)
    {
//        Colormap.clear();
//        Colormap.append(0.0);
//        Colormap.append(0.5);
//        Colormap.append(1.0);
//        Colors.clear();
//        Colors.append("#5477ff");
//        Colors.append("#0023b1");
//        Colors.append("#001462");
//   //     AddComboMap(0,"Flood Height","m",hmx,Colormap,Colors,false,false,1.0,0.01);
//        AddComboMap(0,"Water Height","m",hmxWH,Colormap,Colors,false,false,1.0,0.01);

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
//        AddComboMap(0,"Flood Velocity","m/s",UVflood,Colormap,Colors,false,false,1.0,0.01);

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.25);
        Colormap.append(0.5);
        Colormap.append(0.75);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#A60000");
        Colors.append("#FF0000");
        Colors.append("#FFFF00");
        Colors.append("#00FF00");
        Colors.append("#007300");

        AddComboMap(0,"Flood Start Time","min",floodTimeStart,Colormap,Colors,false,false,1.0,1.0);
        AddComboMap(0,"Flood duration","min",floodTime,Colormap,Colors,false,false,1.0,1.0);

    }

    if(SwitchErosion)
    {
        double step = 0.01;
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.3);
        Colormap.append(0.5);
        Colormap.append(0.70);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#616ca2");//#457A60");
        Colors.append("#50B547");//#96B547");
        Colors.append("#FFFFFF");
        Colors.append("#FFFF00");
        Colors.append("#FF0000");

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

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(0.9);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#FAFAD2");
        Colors.append("#FFFF66");
        Colors.append("#d47e17");//808000
        Colors.append("#804000");

        AddComboMap(1,"Sediment Load","kg/m2",COMBO_SS,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        AddComboMap(1,"Sed Concentration","kg/m3",TotalConc,Colormap,Colors,false,false,1.0, step);
        AddComboMap(1,"Splash detachment","kg/m2",DETSplashCum,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        AddComboMap(1,"Flow detachment","kg/m2",DETFlowCum,Colormap,Colors,false,false,1.0/(_dx*_dx), step);

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#EEEfcd");
        Colors.append("#50B547");//#96B547");
        Colors.append("#616ca2");//#457A60");
        AddComboMap(1,"Deposition","kg/m2",TotalDepMap,Colormap,Colors,false,false,-1.0/(_dx*_dx), step);

        if(SwitchUseGrainSizeDistribution)
        {
            Colormap.clear();
            Colormap.append(0.0);
            Colormap.append(0.5);
            Colormap.append(0.9);
            Colormap.append(1.0);
            Colors.clear();
            Colors.append("#FAFAD2");
            Colors.append("#FFFF66");
            Colors.append("#d47e17");
            Colors.append("#804000");

            FOR_GRAIN_CLASSES
            {
                AddComboMap(1,"Overland S.L. Grain Class " + QString::number(d),"kg/m2",Sed_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
            }
            if(SwitchIncludeChannel)
            {
                if(SwitchUse2Layer)
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
