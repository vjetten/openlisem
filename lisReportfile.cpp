
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
    FillTimeSeriesData();

    std::function<void(int)> freport = std::bind((&TWorld::reportWrapper),this,std::placeholders::_1);
    ThreadPool->RunReportFunction(freport);
}


void TWorld::reportWrapper(int not_used)
{

    ReportTimeseriesNew(not_used);
    // report hydrographs ande sedigraphs at all points in outpoint.map

    ReportTotalsNew(not_used);
    // report totals to a text file

    ReportMaps(not_used);
    // report all maps and mapseries

    ReportLandunits(not_used);
    // reportc stats per landunit class

    ChannelFloodStatistics(not_used);
    // report buildings submerged in flood level classes

}

void TWorld::FillTimeSeriesData()
{

    //save the timeseries data temporarily for multithreading!
    //only pointmap is threadsafe since it is not used!

    TSList_point.clear();
    TSList_rainav.clear();
    TSList_snowav.clear();
    TSList_q.clear();
    TSList_h.clear();
    TSList_qs.clear();
    TSList_c.clear();

    double RainIntavg = RainAvgmm * 3600/_dt;
    double SnowIntavg = SnowAvgmm * 3600/_dt;

    TSList_point.append(0);
    TSList_rainav.append(RainIntavg);
    TSList_snowav.append(SnowIntavg);


    TSList_q.append((QtotT * 1000.0/_dt));
    if(SwitchIncludeChannel)
    {
        double channelwh = 0;
        if(SwitchIncludeChannel)
        {

            FOR_ROW_COL_MV_CH
            {
                if(LDDChannel->Drc == 5)
                {
                    channelwh += UF1D_h->Drc;
                }
            }
        }
        TSList_h.append(channelwh);
    }else
    {
        TSList_h.append(0);
    }

    if(SwitchErosion)
    {
        TSList_qs.append(SoilLossTotT/_dt);
        TSList_c.append((QtotT) > 1e-6? SoilLossTotT/(QtotT) : 0);
    }else
    {
        TSList_qs.append(0.0);
        TSList_c.append(0.0);
    }

    FOR_ROW_COL_MV
    {
        if ( PointMap->Drc > 0 )
        {
            TSList_point.append(PointMap->Drc);
            TSList_rainav.append(RainIntavg);
            TSList_snowav.append(SnowIntavg);

            if(SwitchIncludeChannel)
            {
                TSList_q.append(Qoutput->Drc);
                TSList_h.append(UF1D_h->Drc);
                if(SwitchErosion)
                {
                    TSList_qs.append(Qsoutput->Drc);
                    TSList_c.append(TotalConc->Drc);
                }else
                {
                    TSList_qs.append(0.0);
                    TSList_c.append(0.0);
                }

            }else
            {
                TSList_q.append(0.0);
                TSList_h.append(0.0);
                if(SwitchErosion)
                {
                    TSList_qs.append(0.0);
                    TSList_c.append(0.0);
                }
            }


        }
    }


}

//---------------------------------------------------------------------------
/** fill output structure 'op' with results to talk to the interface:
    report to screen, hydrographs and maps */
void TWorld::OutputUI(void)
{

    //hydrographs
    op.timestep = this->_dt/60.0;
    op.OutletQ.at(0)->append((QtotT * 1000.0/_dt));
    op.OutletQs.at(0)->append(SoilLossTotT/_dt);
    op.OutletC.at(0)->append((QtotT) > 1e-6? SoilLossTotT/(QtotT) : 0);
    op.OutletQtot.replace(0,Qtot);
    op.OutletQstot.replace(0,SoilLossTot/1000.0);

    double channelwh = 0;
    if(SwitchIncludeChannel)
    {

        FOR_ROW_COL_MV_CH
        {
            if(LDDChannel->Drc == 5)
            {
                channelwh += UF1D_h->Drc;
            }
        }
    }
    op.OutletChannelWH.at(0)->append(channelwh);
    for(int j = 1; j < op.OutletIndices.length(); j++)
    {
        int r = op.OutletLocationX.at(j);
        int c = op.OutletLocationY.at(j);

        double discharge = Qoutput->Drc;
        double sedimentdischarge = SwitchErosion? Qsoutput->Drc * _dt : 0.0;
        double sedimentconcentration = SwitchErosion? UF1D_fsConc->Drc : 0.0;
        double channelwh = SwitchIncludeChannel? UF1D_h->Drc : 0.0;

        op.OutletQtot.replace(j,op.OutletQtot.at(j) + _dt * discharge/1000.0);
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

    //display maps
    fill(*COMBO_QOFCH, 0.0);
    calcMap(*COMBO_QOFCH, *UF2D_q, ADD);
    calcMap(*COMBO_QOFCH, *UF1D_q, ADD);

    //output maps for combo box
    for(int i = 0; i < op.ComboMaps.length(); i++)
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
        copy(*op.channelMap, *ChannelWidth);
    //BB 151118 might be better to draw LDD, since that is actually used to determine the presence of a channel
    if (SwitchRoadsystem)
        copy(*op.roadMap, *RoadWidthDX);
    if (SwitchHouses)
        copy(*op.houseMap, *HouseCover);

    copy(*op.vegcover,  *Cover);
    copy(*op.vegheight, *PlantHeight);
    copy(*op.randomroughness,*RR);

    copy(*op.gl_dem_change,*DEMChange);
    copy(*op.gl_flow_height, *UF2D_h);
    copy(*op.gl_flow_u, *UF2D_u);
    copy(*op.gl_flow_v, *UF2D_v);
    copy(*op.gl_flow_c, *UF2D_tConc);

    op.has_channel = SwitchIncludeChannel;
    if(SwitchIncludeChannel)
    {

        op.has_channel = SwitchChannelFlood;

        copy(*op.ch_w, *ChannelWidth);
        if(SwitchChannelFlood)
        {
            copy(*op.ch_d,*ChannelDepth);
        }else
        {
            fill(*op.ch_d,-1.0);
        }

        copy(*op.ch_ldd,*LDDChannel);
        copy(*op.gl_ch_flow_height, *UF1D_h);
        copy(*op.gl_ch_flow_v,*UF1D_velocity);
        copy(*op.gl_ch_flow_c,*UF1D_tConc);
    }

    // MAP DISPLAY VARIABLES

    op.dx = _dx;
    op.MB = MB;
    op.runstep = runstep;
    op.maxstep = (int) ((EndTime-BeginTime)/_dt);
    op.EndTime = EndTime/60.0;
    //	op.BeginTime = BeginTime/60.0;

    op.CatchmentArea = CatchmentArea;

    op.BaseFlowtot = BaseFlow * 1000.0/(_dx*_dx*nrCells);
    op.LitterStorageTot = mapTotal(*LInterc) *1000.0/(_dx*_dx*nrCells);
    op.ChannelVolTot = SwitchIncludeChannel? mapTotal(*UF1D_f) * 1000.0/(_dx*_dx*nrCells) : 0.0;

    op.RainTotmm = RainTotmm + SnowTotmm;
    op.WaterVolTotmm = WaterVolRunoffmm;//WaterVolTotmm-SurfStoremm;
    op.Qtotmm = Qtotmm;
    op.Qtot = QtotOutlet;

    op.RainpeakTime = RainpeakTime/60;
    op.FloodTotMax = floodVolTotMax;
    op.FloodAreaMax = floodAreaMax;

    op.InfilTotmm = InfilTotmm;
    op.SurfStormm = SurfStoremm;
    op.IntercTotmm = IntercTotmm;// + IntercHouseTotmm;
    //houses
    op.IntercHouseTotmm = IntercHouseTotmm;
    op.InfilKWTotmm = InfilKWTot; // infil part in kin wave not used
    op.RunoffFraction = (RainTotmm > 0 ? Qtotmm/RainTotmm : 0);

    op.MBs = MBs;
    op.DetTotSplash=DetSplashTot*0.001; // convert from kg to ton per cell
    op.DetTotFlow=DetFlowTot*0.001; // convert from kg to ton
    op.DepTot=DepTot*0.001; // convert from kg to ton
    op.SedTot=SedTot*0.001; // convert from kg to ton

    op.ChannelDetTot=ChannelDetTot*0.001; // convert from kg to ton
    op.ChannelDepTot=ChannelDepTot*0.001; // convert from kg to ton
    op.ChannelSedTot=ChannelSedTot*0.001; // convert from kg to ton

    op.FloodSed = FloodSedTot*0.001;
    op.FloodDepTot = FloodDepTot*0.001;
    op.FloodDetTot = FloodDetTot*0.001;

    op.SoilLossTot=SoilLossTot/*Outlet*/ *0.001; // convert from kg to ton

    op.t = time_ms.elapsed()*0.001/60.0;
    op.time = time/60;
    op.maxtime = op.t/runstep * op.maxstep;

    op.Pmm = (RainAvgmm + SnowAvgmm)*3600/_dt;
    op.volFloodmm = floodTotmm;

    op.rain_average = RainAvgmm*3600/_dt;
}
//---------------------------------------------------------------------------
/** reporting timeseries for every non zero point PointMap
 - 3 types of output: PCRaster timeplot format; SOBEK input format; flat comma delimited format
 - all points in one file or each point in a separate file
 - the types should be mututally exclusive in the interface and run file
*/
void TWorld::ReportTimeseriesNew(int not_used)
{
    int nr = 0;
    int hour = 0;
    int min = 0;
    int sec = 0;
    //int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
    double RainIntavg = RainAvgmm * 3600/_dt;
    double SnowIntavg = SnowAvgmm * 3600/_dt;
    QString newname1, pnr, sep = (SwitchWritePCRtimeplot ? " " : ",");
    int width = (!SwitchWritePCRtimeplot ? 0 : 9);
    // NOTE if SwitchWriteCommaDelimited = true then SwitchWritePCRtimeplot = false

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
            for(int p = 0; p < TSList_point.length(); p++)
                {
                    newname1 = fi.path() + "/" + fi.baseName() + "_" +
                            QString::number((int)TSList_point.at(p)) + "." +  fi.suffix();
                    // make filename using point number

                    QFile fout(newname1);
                    fout.open(QIODevice::WriteOnly | QIODevice::Text);
                    QTextStream out(&fout);
                    out.setRealNumberPrecision(3);
                    out.setFieldWidth(width);
                    out.setRealNumberNotation(QTextStream::FixedNotation);

                    // HEADERS for the 3 types
                    if (SwitchWritePCRtimeplot)  //PCRaster timeplot format, cannot be SOBEK !
                    {
                        pnr.setNum((int)TSList_point.at(p));
                        out << "#LISEM total flow and sed output file for point " << pnr << "\n";

                        // nr columns is time + rain (+ maybe snow) + Q + (maybe Qs + C)
                        int nrs = 3 + (SwitchErosion ? 2 : 0);
                        if (SwitchSnowmelt && SwitchRainfall) nrs++;
                        pnr.setNum(nrs);
                        out << pnr << "\n";

                        out << "run step\n";
                        if (SwitchRainfall) out << "Pavg (mm/h)\n";
                        if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                        out << "Qall (l/s)\n" << "chanWH (m)\n";
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
                            out << "* " <<  (int)TSList_point.at(p) << "\n";
                        }
                        else // flat format, comma delimited
                        {
                            pnr.setNum((int)TSList_point.at(p));
                            out << "LISEM total flow and sed output file for point " << pnr << "\n";
                            out << "Time";
                            if (SwitchRainfall) out << ",Pavg";
                            if (SwitchSnowmelt) out << ",Snowavg";
                            out << ",Q" << ",chanWH";
                            if (SwitchErosion) out << ",Qs,C";
                            out << "\n";
                            out << "min,mm/h";
                            if (SwitchSnowmelt) out << ",mm/h";
                            out << ",l/s" << ",m";
                            if (SwitchIncludeTile) out << ",l/s";
                            if (SwitchErosion) out << ",kg/s,g/l";
                            out << "\n";
                        }
                    fout.close();
                }
        }  // separate
        else   // HEADERS: all points in one file
        {

            newname1 = resultDir + outflowFileName;

            QFile fout(newname1);
            fout.open(QIODevice::WriteOnly | QIODevice::Text);
            QTextStream out(&fout);
            out.setRealNumberPrecision(3);
            out.setFieldWidth(width);
            out.setRealNumberNotation(QTextStream::FixedNotation);

            if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
            {
                // count nr of points
                nr = TSList_point.length();

                // nr columns is time + rain (+ maybe snow) + nr points*(Q + Qs + C)
                int nrs = 2+(1+(SwitchErosion ? 2 : 0))*nr;
                if (SwitchSnowmelt && SwitchRainfall) nrs++;
                pnr.setNum(nrs);

                out << "#LISEM total flow and sed output file for all reporting points\n";
                out <<  pnr << "\n";
                out << "Time (min)\n";
                if (SwitchRainfall) out << "Pavg (mm/h)\n";
                if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                for(int p = 0; p < TSList_point.length(); p++)
                    {
                    pnr.setNum((int)TSList_point.at(p));
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
                    for(int p = 0; p < TSList_point.length(); p++)
                        {
                        pnr.setNum((int)TSList_point.at(p));
                        out << " Q #" << pnr;
                        if (SwitchErosion) out << " Qs #" << pnr;
                        if (SwitchErosion) out << " C #" << pnr;
                    }
                    out << "\n";
                    for(int p = 0; p < TSList_point.length(); p++)
                        {
                        pnr.setNum((int)TSList_point.at(p));
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
                    for(int p = 0; p < TSList_point.length(); p++)
                        {
                        pnr.setNum((int)TSList_point.at(p));
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
                    for(int p = 0; p < TSList_point.length(); p++)
                        {
                        pnr.setNum((int)TSList_point.at(p));
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
        int point = 0;
        for(int p = 0; p < TSList_point.length(); p++)
            {
                point ++;

           //     qDebug() << PointMap->Drc << r << c;
                newname1 = fi.path() + "/" + fi.baseName() + "_" +
                        QString::number((int)TSList_point.at(p)) + "." +  fi.suffix();
                QFile fout(newname1);
                fout.open(QIODevice::Append | QIODevice::Text);

                QTextStream out(&fout);
                out.setRealNumberPrecision(3);
                out.setFieldWidth(width);
                out.setRealNumberNotation(QTextStream::FixedNotation);

                if (!SwitchSOBEKoutput)   //PCRaster timeplot and flat CSV format
                {

                    if (SwitchWritePCRtimeplot)
                        out << runstep;
                    else
                        out << time/60;
                    if (SwitchRainfall) out << sep << TSList_rainav.at(point);
                    if (SwitchSnowmelt) out << sep << TSList_snowav.at(point);
                    out << sep << TSList_q.at(point) << sep << TSList_h.at(point);
                    if (SwitchErosion) out << sep << TSList_qs.at(point);
                    if (SwitchErosion) out << sep << TSList_c.at(point);
                    out << "\n";
                }
                else  //SOBEK format
                {
                    QString ss = QString("\"%1;%2:%3:%4\" %5 <\n").
                            arg(SOBEKdatestring).
                            arg((uint)hour,2,10,QLatin1Char('0')).
                            arg((uint)min,2,10,QLatin1Char('0')).
                            arg((uint)sec,2,10,QLatin1Char('0')).
                            arg(TSList_q.at(point)/1000.0,0,'f',3);
                    out << ss;
                }
                fout.close();
            }  // if point
    } //switch separate
    else
    {
        // all points in one file
        newname1 = resultDir + outflowFileName;
        // use simply resultdir + filename
        QFile fout(newname1);
        fout.open(QIODevice::Append | QIODevice::Text);
        QTextStream out(&fout);
        out.setRealNumberPrecision(3);
        out.setFieldWidth(width);
        out.setRealNumberNotation(QTextStream::FixedNotation);

        if (!SwitchSOBEKoutput)
        {
            if (SwitchWritePCRtimeplot)
                out << runstep;
            else
                out << time/60;

            if(TSList_rainav.length() > 0)
            {
                if (SwitchRainfall) out << sep << TSList_rainav.at(0);

            }
            if(TSList_snowav.length() > 0)
            {
                if (SwitchSnowmelt) out << sep << TSList_snowav.at(0);
            }

            int point = 0;
            for(int p = 0; p < TSList_point.length(); p++)
                {
                    point ++;

                    out << sep << TSList_q.at(point) << sep <<TSList_h.at(point);
                    if (SwitchErosion) out << sep << TSList_qs.at(point);
                    if (SwitchErosion) out << sep << TSList_c.at(point);
                }

            out << "\n";
        }
        else
            if (SwitchSOBEKoutput)
            {
                out.setFieldWidth(2);
                out << "\"" << SOBEKdatestring << ":" << hour << ":" <<  min << ":" <<  sec;
                out.setFieldWidth(8);
                int point = 0;
                for(int p = 0; p < TSList_point.length(); p++)
                    {
                        point ++;
                        out << " " << TSList_q.at(point)/1000.0;
                        if (SwitchErosion) out << " " << TSList_qs.at(point);
                        if (SwitchErosion) out << " " << TSList_c.at(point);
                    }

                out << " < \n";
            }
        fout.close();
    }
}
//---------------------------------------------------------------------------
/// Report totals of the main outlet nd general values for the catchment to a comma delimited text file
void TWorld::ReportTotalsNew(int not_used)
{
    QFile fp(resultDir + resultFileName);
    if (!fp.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&fp);
    out.setRealNumberPrecision(5);
    out.setFieldWidth(12);
    out.setRealNumberNotation(QTextStream::FixedNotation);
    out << "\"LISEM run with:," << op.runfilename << "\"\n";
    out << "\"LISEM results at time (min):," << op.time <<"\"\n";
    out << "\"Catchment area              (ha):\"," << op.CatchmentArea/10000.0<< "\n";
    out << "\"Total rainfall              (mm):\"," << op.RainTotmm<< "\n";
    out << "\"Total discharge             (mm):\"," << op.Qtotmm<< "\n";
    out << "\"Total interception          (mm):\"," << op.IntercTotmm<< "\n";
    out << "\"Total House interception    (mm):\"," << op.IntercHouseTotmm<< "\n";
    out << "\"Total infiltration          (mm):\"," << op.InfilTotmm<< "\n";
    out << "\"Surface storage             (mm):\"," << op.SurfStormm<< "\n";
    out << "\"Water in runoff + channel   (mm):\"," << op.WaterVolTotmm<< "\n";
    out << "\"Total discharge             (m3):\"," << op.Qtot<< "\n";
    out << "\"Peak time rainfall         (min):\"," << op.RainpeakTime<< "\n";
    out << "\"Discharge/Rainfall           (%):\"," << op.RunoffFraction*100<< "\n";
    out << "\"Flood volume (max level)    (m3):\"," << op.FloodTotMax<< "\n";
    out << "\"Flood area (max level)      (m2):\"," << op.FloodAreaMax<< "\n";
    out << "\"Splash detachment (land)   (ton):\"," << op.DetTotSplash<< "\n";
    out << "\"Flow detachment (land)     (ton):\"," << op.DetTotFlow<< "\n";
    out << "\"Deposition (land)          (ton):\"," << op.DepTot<< "\n";
    out << "\"Suspended Sediment (land)  (ton):\"," << op.SedTot<< "\n";
    out << "\"Flow detachment (channels) (ton):\"," << op.ChannelDetTot<< "\n";
    out << "\"Deposition (channels)      (ton):\"," << op.ChannelDepTot<< "\n";
    out << "\"Susp. Sediment (channels)  (ton):\"," << op.ChannelSedTot<< "\n";
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
void TWorld::ReportMaps(int not_used)
{
    FOR_ROW_COL_MV
    {
        tm->Drc = RainCumFlat->Drc * 1000.0; // m to mm
        ThreadPool->tma->Drc = (Interc->Drc + IntercHouse->Drc)*1000.0/CellArea->Drc;
    }

    //creating a wig-shaped map
    /*for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            ThreadPool->tm->Drc = 1.0 * std::fabs(double(float(r) - 0.5 * _nrRows)) + 2.0 * double(_nrCols - float(c));

        }
    }*/

    //creating a test dem for a landslide dam
    /*for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            double x = std::fabs(double(float(r) - 0.5 * _nrRows)) < 3 ? std::pow((std::fabs(double(float(r) - 0.5 * _nrRows))/3.0),4.0)*3.0 : std::fabs(double(float(r) - 0.5 * _nrRows));

            ThreadPool->tm->Drc = 15.0 * x - 0.05 * double(_nrCols - float(c));

        }
    }*/

    //creating a flume
    /*double slope1 = std::tan(std::asin(10.0 *3.14/180));
    double slope2 = std::tan(std::asin(12.0 *3.14/180));
    double slope3 = std::tan(std::asin(14.0 *3.14/180));
    double slope4 = std::tan(std::asin(16.0 *3.14/180));
    double slope5 = std::tan(std::asin(18.0 *3.14/180));

    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            ThreadPool->tm->Drc = 1.0 + slope1 * std::fabs(double(float(_nrRows -r))) * _dx + 0.05 * _dx * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        }
    }
    for(int r = 1; r < _nrRows-1; r++)
    {
        for (int c = 1; c < _nrCols-1; c++)
        {
            ThreadPool->tm->Drc = (ThreadPool->tm->Drc * 2.0 + (!OUTORMV(r,c-1)? ThreadPool->tm->data[r][c-1] : ThreadPool->tm->Drc)+ (!OUTORMV(r,c+1)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r-1,c)? ThreadPool->tm->data[r-1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r+1,c)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc))/6.0;
        }
    }
    report(*ThreadPool->tm, QString("dem1"));
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            ThreadPool->tm->Drc = 1.0 + slope2 * std::fabs(double(float(_nrRows -r))) * _dx + 0.05 * _dx * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        }
    }
    for(int r = 1; r < _nrRows-1; r++)
    {
        for (int c = 1; c < _nrCols-1; c++)
        {
            ThreadPool->tm->Drc = (ThreadPool->tm->Drc * 2.0 + (!OUTORMV(r,c-1)? ThreadPool->tm->data[r][c-1] : ThreadPool->tm->Drc)+ (!OUTORMV(r,c+1)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r-1,c)? ThreadPool->tm->data[r-1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r+1,c)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc))/6.0;
        }
    }
    report(*ThreadPool->tm, QString("dem2"));
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            ThreadPool->tm->Drc = 1.0 + slope3 * std::fabs(double(float(_nrRows -r))) * _dx + 0.05 * _dx * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        }
    }
    for(int r = 1; r < _nrRows-1; r++)
    {
        for (int c = 1; c < _nrCols-1; c++)
        {
            ThreadPool->tm->Drc = (ThreadPool->tm->Drc * 2.0 + (!OUTORMV(r,c-1)? ThreadPool->tm->data[r][c-1] : ThreadPool->tm->Drc)+ (!OUTORMV(r,c+1)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r-1,c)? ThreadPool->tm->data[r-1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r+1,c)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc))/6.0;
        }
    }
    report(*ThreadPool->tm, QString("dem3"));
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            ThreadPool->tm->Drc = 1.0 + slope4 * std::fabs(double(float(_nrRows -r))) * _dx + 0.05 * _dx * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        }
    }
    for(int r = 1; r < _nrRows-1; r++)
    {
        for (int c = 1; c < _nrCols-1; c++)
        {
            ThreadPool->tm->Drc = (ThreadPool->tm->Drc * 2.0 + (!OUTORMV(r,c-1)? ThreadPool->tm->data[r][c-1] : ThreadPool->tm->Drc)+ (!OUTORMV(r,c+1)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r-1,c)? ThreadPool->tm->data[r-1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r+1,c)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc))/6.0;
        }
    }
    report(*ThreadPool->tm, QString("dem4"));
    for(int r = 0; r < _nrRows; r++)
    {
        for (int c = 0; c < _nrCols; c++)
        {
            ThreadPool->tm->Drc = 1.0 + slope5 * std::fabs(double(float(_nrRows -r))) * _dx + 0.05 * _dx * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

        }
    }
    for(int r = 1; r < _nrRows-1; r++)
    {
        for (int c = 1; c < _nrCols-1; c++)
        {
            ThreadPool->tm->Drc = (ThreadPool->tm->Drc * 2.0 + (!OUTORMV(r,c-1)? ThreadPool->tm->data[r][c-1] : ThreadPool->tm->Drc)+ (!OUTORMV(r,c+1)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r-1,c)? ThreadPool->tm->data[r-1][c] : ThreadPool->tm->Drc)+ (!OUTORMV(r+1,c)? ThreadPool->tm->data[r+1][c] : ThreadPool->tm->Drc))/6.0;
        }
    }
    report(*ThreadPool->tm, QString("dem5"));*/

    report(*ThreadPool->tm, rainfallMapFileName);
    report(*ThreadPool->tma, interceptionMapFileName);

    report(*InfilmmCum, infiltrationMapFileName);

    report(*runoffTotalCell, runoffMapFileName); // in m3, total runoff from cell (but there is also runon!)
    report(*runoffFractionCell, runoffFractionMapFileName);

    report(*WHmax, floodWHmaxFileName);


    if(SwitchErosion)
    {
        // deal with erosion units, 0 = ton/ha, 1 = kg/m2, 2 = kg/cell
        if (ErosionUnits == 0)
            fill(*ThreadPool->tma,10.0);
        else
            fill(*ThreadPool->tma,1.0);
        if (ErosionUnits == 2 || ErosionUnits == 0)
            calcMap(*ThreadPool->tma, *CellArea, DIV);


        // all detachment combined
        FOR_ROW_COL_MV
        {
            ThreadPool->tm->Drc =std::max(0.0,TotalDetMap->Drc + TotalDepMap->Drc);
        }
        calcMap(*(ThreadPool->tm), *(ThreadPool->tma), MUL);
        report(*(ThreadPool->tm), totalErosionFileName);
        if (outputcheck[5].toInt() == 1)
            report(*(ThreadPool->tm), Outeros); // in units

        // all deposition combined
        FOR_ROW_COL_MV
        {
            tm->Drc =std::min(0.0,TotalDetMap->Drc + TotalDepMap->Drc);
        }
        calcMap(*(ThreadPool->tm), *(ThreadPool->tma), MUL);
        report(*(ThreadPool->tm), totalDepositionFileName);
        if (outputcheck[6].toInt() == 1)
            report(*(ThreadPool->tm), Outdepo); // in units

        // all channel depostion combined
        if (SwitchIncludeChannel)
        {
            FOR_ROW_COL_MV_CH
            {
                ThreadPool->tm->Drc =std::min(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc);
            }
            calcMap(*(ThreadPool->tm), *(ThreadPool->tma), MUL);
            report(*(ThreadPool->tm), totalChanDepositionFileName);

            // all channel detachment combined
            FOR_ROW_COL_MV_CH
            {
                ThreadPool->tm->Drc =std::max(0.0,TotalChanDetMap->Drc + TotalChanDepMap->Drc);
            }
            calcMap(*(ThreadPool->tm), *(ThreadPool->tma), MUL);
            report(*(ThreadPool->tm), totalChanErosionFileName);
        }
        copy(*(ThreadPool->tm), *TotalSoillossMap); //kg/cell
        calcMap(*(ThreadPool->tm), *(ThreadPool->tma), MUL);
        report(*(ThreadPool->tm), totalSoillossFileName);
        if (outputcheck[16].toInt() == 1) report(*tm, OutSL);      // in user units

        // total sediment
        copy(*(ThreadPool->tm), *COMBO_SS); //kg/cell
        //calcMap(*(ThreadPool->tm), 1/(_dx*_dx), MUL);
        if (outputcheck[17].toInt() == 1) report(*(ThreadPool->tm), OutSed);      // in user units

        if (outputcheck[1].toInt() == 1) report(*UF2D_sConc, Outconc);  // in g/l
        if (outputcheck[4].toInt() == 1) report(*UF2D_sstc, Outtc);      // in g/l

    }

    report(*UF2D_f, FinalFluidPhaseFileName);

    if (SwitchChannelFlood)
    {
        report(*floodHmxMax, floodLevelFileName);
        report(*floodTime, floodTimeFileName);
        report(*floodTimeStart, floodFEWFileName);
        report(*maxChannelflow, floodMaxQFileName);
        report(*maxChannelWH, floodMaxChanWHFileName);
        report(*floodVMax, floodMaxVFileName);
    }

    if(SwitchSolidPhase)
    {
        report(*UF2D_s, FinalSolidPhaseFileName);
        report(*MaximumDebrisFlowHeight,FileName_MaxDebrisFlowHeight);
        report(*MaximumDebrisFlowVelocity,FileName_MaxDebrisFlowVelocity);
        report(*dfTimeStart,FileName_DebrisFlowStart);
    }

    if(SwitchSlopeFailure)
    {
        report(*TotalSlopeFailure,FileName_SlopeFailure);
    }
    if(SwitchSlopeStability)
    {
        report(*MinimumSafetyFactor,FileName_MinimumSafetyFactor);
    }
    if(SwitchEntrainment)
    {
        report(*TotalEntrainmentDet,FileName_EntrainmentDet);
        report(*TotalEntrainmentDep,FileName_EntrainmentDep);
    }


    if (outputcheck[0].toInt() == 1)
        report(*Qoutput, Outrunoff); // in l/s
    if (outputcheck[2].toInt() == 1)
    {
        calcMapValue(*(ThreadPool->tm), *WH, 1000, MUL);// WH in mm
        report(*(ThreadPool->tm), Outwh);
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
        calcMapValue(*(ThreadPool->tm), *WHstore, 1000, MUL);// in mm
        report(*(ThreadPool->tm), Outss);
        /** TODO check this: surf store in volume m3 is multiplied by flowwidth? */
    }

    if (outputcheck[10].toInt() == 1) report(*UF1D_f, Outchvol);


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
        if (outputcheck[16].toInt() == 1)
        {
            report(*hmxWH, OutHmxWH);
        }
    }

    if(SwitchSlopeStability)
    {

        if (outputcheck[18].toInt() == 1)
        {
            report(*DFSafetyFactor, OutSafetyFactor);
        }
    }

    if(SwitchSlopeFailure)
    {
        if (outputcheck[19].toInt() == 1)
        {
            report(*DFInitiationHeight, OutSlopeFailure);
        }
    }

    if(SwitchSolidPhase)
    {
        if (outputcheck[20].toInt() == 1)
        {
            report(*dfhmx, OutDFH);
        }
        if (outputcheck[21].toInt() == 1)
        {
            report(*dfUV, OutDFV);
        }
        if (outputcheck[22].toInt() == 1)
        {
            report(*UF2D_FPH, "fph");
        }
        if (outputcheck[23].toInt() == 1)
        {
            report(*UF2D_SPH, "sph");
        }

    }

    if(SwitchEntrainment)
    {
        if (outputcheck[24].toInt() == 1)
        {
            report(*TotalEntrainmentDet, "ent");
        }

    }

    if (outputcheck[25].toInt() == 1)
    {
        report(*UF2D_TimeStep, "ts");
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
void TWorld::ReportLandunits(int not_used)
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
void TWorld::ChannelFloodStatistics(int not_used)
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
    name = resultDir + QFileInfo(floodStatsFileName).baseName()+timestamp+".csv";
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
        delete op.vegcover;
        delete op.vegheight;
        delete op.randomroughness;
        delete op.gl_flow_height;
        delete op.gl_flow_u;
        delete op.gl_flow_v;
        delete op.gl_flow_c;
        delete op.gl_flow_height;
        delete op.ch_w;
        delete op.ch_d;
        delete op.ch_ldd;
        delete op.gl_ch_flow_height;
        delete op.gl_ch_flow_v;
        delete op.gl_ch_flow_c;
        delete op.gl_dem_change;
    }

    op.baseMap = new cTMap();
    op.baseMapDEM = new cTMap();
    op.channelMap = new cTMap();
    op.roadMap = new cTMap();
    op.houseMap = new cTMap();
    op.vegcover = new cTMap();
    op.vegheight = new cTMap();
    op.randomroughness = new cTMap();
    op.gl_flow_height = new cTMap();
    op.gl_flow_u = new cTMap();
    op.gl_flow_v = new cTMap();
    op.gl_flow_c = new cTMap();
    op.ch_w = new cTMap();
    op.ch_d = new cTMap();
    op.ch_ldd = new cTMap();
    op.gl_ch_flow_height = new cTMap();
    op.gl_ch_flow_v = new cTMap();
    op.gl_ch_flow_c = new cTMap();
    op.gl_dem_change = new cTMap();

    op.baseMap->MakeMap(LDD, 0);
    op.baseMapDEM->MakeMap(LDD, 0);
    op.channelMap->MakeMap(LDD, 0);
    op.roadMap->MakeMap(LDD, 0);
    op.houseMap->MakeMap(LDD, 0);
    op.vegcover->MakeMap(LDD, 0);
    op.vegheight->MakeMap(LDD, 0);
    op.randomroughness->MakeMap(LDD, 0);
    op.gl_flow_height->MakeMap(LDD, 0);
    op.gl_flow_u->MakeMap(LDD, 0);
    op.gl_flow_v->MakeMap(LDD, 0);
    op.gl_flow_c->MakeMap(LDD, 0);
    op.ch_w->MakeMap(LDD, 0);
    op.ch_d->MakeMap(LDD, 0);
    op.ch_ldd->MakeMap(LDD, 0);
    op.gl_ch_flow_height->MakeMap(LDD, 0);
    op.gl_ch_flow_v->MakeMap(LDD, 0);
    op.gl_ch_flow_c->MakeMap(LDD, 0);
    op.gl_dem_change->MakeMap(LDD, 0);
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
    int maxnr = 0;

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
    qDebug() << "dx read from map is " << _dx;

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
    AddComboMap(0,"Flow Height","m",UF2D_h,Colormap,Colors,false,false,1.0, 0.01);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"Channel Flow Height","m",UF1D_h,Colormap,Colors,false,false,1.0, 0.01);
    }

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
    AddComboMap(0,"Overland Flow Velocity","m/s",UF2D_velocity,Colormap,Colors,false,false,1.0, 0.01);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"Channel Velocity","m/s",UF1D_velocity,Colormap,Colors,false,false,1.0,0.01);
    }

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
    AddComboMap(0,"Overland Discharge","l/s",UF2D_q,Colormap,Colors,true,false,1000.0, 1.0);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"Channel Discharge","l/s",UF1D_q,Colormap,Colors,true,false,1000.0, 1.0);
        AddComboMap(0,"Total Discharge","l/s",COMBO_QOFCH,Colormap,Colors,true,false,1000.0, 1.0);
    }



    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.5);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#5477ff");
    Colors.append("#0023b1");
    Colors.append("#001462");

    AddComboMap(0,"Fluid Phase Volume","m3",UF2D_f,Colormap,Colors,false,false,1.0,1.0);

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.5);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#ffffff");
    Colors.append("#50B547");//#96B547");
    Colors.append("#616ca2");//#457A60");

    AddComboMap(0,"Solid Phase Volume","m3",UF2D_s,Colormap,Colors,false,false,1.0,1.0);
    //AddComboMap(0,"Solid-Fluid Drag Coefficient","-",UF2D_DC,Colormap,Colors,false,false,1.0,1.0);

    if(SwitchIncludeChannel)
    {

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#5477ff");
        Colors.append("#0023b1");
        Colors.append("#001462");
        AddComboMap(0,"Channel Fluid Phase Volume","m3",UF1D_f,Colormap,Colors,false,false,1.0,1.0);

        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#ffffff");
        Colors.append("#50B547");//#96B547");
        Colors.append("#616ca2");//#457A60");
        AddComboMap(0,"Channel Solid Phase Volume","m3",UF1D_s,Colormap,Colors,false,false,1.0,1.0);

    }

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
    Colors.append("#FF0000");

    AddComboMap(0,"Rainfall Cumulative","mm",RainCumFlat,Colormap,Colors,false,false,1.0,1.0);
    AddComboMap(0,"Rainfall Intensity","mm/h",Rain,Colormap,Colors,false,false,1.0,1.0);

    if(UF_DISPLAYFLOODMINIMUM > 0)
    {
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#5477ff");
        Colors.append("#0023b1");
        Colors.append("#001462");
        AddComboMap(0,"Flood Height","m",hmx,Colormap,Colors,false,false,1.0,0.01);

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
        AddComboMap(0,"Flood Velocity","m/s",UVflood,Colormap,Colors,false,false,1.0,0.01);

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
    }

    if(SwitchSolidPhase && UF_DISPLAYDEBRISFLOWMINIMUM > 0)
    {
        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#5477ff");
        Colors.append("#0023b1");
        Colors.append("#001462");
        AddComboMap(0,"Debris Flow Height","m",dfhmx,Colormap,Colors,false,false,1.0,0.01);

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
        AddComboMap(0,"Debris Flow Velocity","m/s",dfUV,Colormap,Colors,false,false,1.0,0.01);

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

        AddComboMap(0,"Debris Flow Start Time","min",dfTimeStart,Colormap,Colors,false,false,1.0,1.0);


    }

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

    AddComboMap(0,"Viscosity"," ",UF2D_visc,Colormap,Colors,false,false,1.0, 0.05);
    AddComboMap(0,"lax factor"," ",UF2D_Test,Colormap,Colors,false,false,1.0, 0.05);

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

    AddComboMap(0,"UF2D_DTStep" ,"kg/m2",UF2D_DTStep,Colormap,Colors,false,false,1.0,1.0);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"UF1D_DTStep" ,"kg/m2",UF1D_DTStep,Colormap,Colors,false,false,1.0,1.0);
    }
    AddComboMap(0,"UF2D_CoreMask" ,"kg/m2",ThreadPool->CoreMask2d,Colormap,Colors,false,false,1.0,1.0);
    if(SwitchIncludeChannel)
    {
        AddComboMap(0,"UF1D_CoreMask" ,"kg/m2",ThreadPool->CoreMask1d,Colormap,Colors,false,false,1.0,1.0);
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

        //AddComboMap(1,"Sediment Load","kg/m2",COMBO_SS,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        AddComboMap(1,"Sed Concentration","kg/m3",TotalConc,Colormap,Colors,false,false,1.0, step);
        if(SwitchIncludeChannel)
        {
            //AddComboMap(1,"Sediment Load","kg/m2",COMBO_SS,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
            AddComboMap(1,"Channel Sed Concentration","kg/m3",TotalConc,Colormap,Colors,false,false,1.0, step);

        }
        AddComboMap(1,"Splash detachment","kg/m2",DETSplashCum,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        AddComboMap(1,"Flow detachment","kg/m2",DETFlowCum,Colormap,Colors,false,false,1.0/(_dx*_dx), step);

        //AddComboMap(1,"UF1D_blm","kg/m2",UF1D_blm,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        //AddComboMap(1,"UF1D_ssm","kg/m2",UF1D_ssm,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        //AddComboMap(1,"UF2D_blm","kg/m2",UF2D_blm,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        //AddComboMap(1,"UF2D_ssm","kg/m2",UF2D_ssm,Colormap,Colors,false,false,1.0/(_dx*_dx), step);

        //AddComboMap(1,"UF1D_bltc","kg/m2",UF1D_bltc,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        //AddComboMap(1,"UF1D_sstc","kg/m2",UF1D_sstc,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        //AddComboMap(1,"UF2D_bltc","kg/m2",UF2D_bltc,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
        //AddComboMap(1,"UF2D_sstc","kg/m2",UF2D_sstc,Colormap,Colors,false,false,1.0/(_dx*_dx), step);


        Colormap.clear();
        Colormap.append(0.0);
        Colormap.append(0.5);
        Colormap.append(1.0);
        Colors.clear();
        Colors.append("#ffffff");
        Colors.append("#50B547");
        Colors.append("#616ca2");
        AddComboMap(1,"Deposition","kg/m2",TotalDepMap,Colormap,Colors,false,false,-1.0/(_dx*_dx), step);

        if(SwitchEntrainment)
        {
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

            AddComboMap(1,"Entrainment","m",TotalEntrainmentDet,Colormap,Colors,false,false,1.0/(_dx*_dx), step);

            if(SwitchDeposition || SwitchCompaction)
            {
                AddComboMap(1,"Deposition","m",TotalEntrainmentDep,Colormap,Colors,false,false,1.0/(_dx*_dx), step);
                //AddComboMap(1,"DepositionT","m",DepositionT,Colormap,Colors,false,false,1.0, step);
            }
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

            AddComboMap(1,"Entr. ShearStress","kg/m2",UF2D_ST,Colormap,Colors,false,false,1.0, step);
            //AddComboMap(1,"test2","kg/m2",Entrainmentshearstressc,Colormap,Colors,false,false,1.0, step);
            //AddComboMap(1,"test1","kg/m2",Entrainmentshearstress,Colormap,Colors,false,false,1.0, step);
            //AddComboMap(1,"sqx","kg/m2",UF2D_sqx,Colormap,Colors,false,false,1.0, step);
            //AddComboMap(1,"sqy","kg/m2",UF2D_sqy,Colormap,Colors,false,false,1.0, step);
        }

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

        AddComboMap(1,"DEMChange","m",DEMChange,Colormap,Colors,false,true,1.0, step);

        if(SwitchSlopeStability)
        {
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

            AddComboMap(1,"Safety Factor","-",DFSafetyFactor,Colormap,Colors,false,false,1.0,0.2);
            AddComboMap(1,"Forcing","-",DFForcing,Colormap,Colors,false,false,1.0,0.2);
            AddComboMap(1,"Forcing up","-",DFForcingUp,Colormap,Colors,false,false,1.0,0.2);

            if(SwitchBedrock)
            {
                AddComboMap(1,"BR Depth","m",DFSoilDepth2,Colormap,Colors,false,false,1.0,0.2);
                AddComboMap(1,"BR Safety Factor","-",DFSafetyFactor2,Colormap,Colors,false,false,1.0,0.2);
                AddComboMap(1,"BR Forcing","-",DFForcing2,Colormap,Colors,false,false,1.0,0.2);
                AddComboMap(1,"BR Forcing up","-",DFForcingUp2,Colormap,Colors,false,false,1.0,0.2);

            }
            if(SwitchSeismic)
            {


            }
        }
        if(SwitchSlopeFailure)
        {
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

            AddComboMap(1,"Failure Height","m",DFTotalInitiationHeight,Colormap,Colors,false,false,1.0,0.2);
        }
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
            Colors.append("#808000");
            Colors.append("#804000");

            FOR_GRAIN_CLASSES
            {
                AddComboMap(1,"Flood BL S.L. Grain Class " + QString::number(d),"kg/m2",UF2D_blm_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
            }
            FOR_GRAIN_CLASSES
            {
                AddComboMap(1,"Flood SS S.L. Grain Class " + QString::number(d),"kg/m2",UF2D_ssm_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx),step);
            }

            if(SwitchIncludeChannel)
            {

                FOR_GRAIN_CLASSES
                {
                    AddComboMap(1,"Channel BL S.L. Grain Class " + QString::number(d),"kg/m2",UF1D_blm_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
                }
                FOR_GRAIN_CLASSES
                {
                    AddComboMap(1,"Channel SS S.L. Grain Class " + QString::number(d),"kg/m2",UF1D_ssm_D.at(d),Colormap,Colors,false,false,1.0/(_dx*_dx), step);
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
    op.ComboScaling.append(scale);
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
