
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

#include "lisemqt.h"
#include "model.h"
#include "global.h"


//---------------------------------------------------------------------------
/** fill output structure 'op' with results to talk to the interface:
    report to screen, hydrographs and maps */
void TWorld::OutputUI(void)
{
    if (op.drawMapType == 1)
        op.DrawMap->copy(Qoutput);  //all output in m3/s
    if (op.drawMapType == 2)
        op.DrawMap->copy(InfilmmCum);  //infil in mm
    if (SwitchErosion && op.drawMapType == 3)
    {
        tmb->calc2Maps(TotalSoillossMap,CellArea, DIV);//from kg/cell to kg/m2
        tmb->calcValue(10, MUL); // from kg/m2 to ton/ha

        op.DrawMap->copy(tmb);  //soilloss in ton/ha
    }
    if (SwitchChannelFlood && op.drawMapType == 4)
    {
        FOR_ROW_COL_MV
                tmb->Drc = hmx->Drc < 0.01 ? 0 : hmx->Drc;
        op.DrawMap->copy(tmb);  //flood level in m
    }

    op.DrawMap1->copy(Qoutput);  //all output in m3/s
    op.DrawMap2->copy(InfilmmCum);  //infil in mm
    if (SwitchErosion)
    {
        tmb->calc2Maps(TotalSoillossMap,CellArea, DIV);
        tmb->calcValue(10, MUL); /* in kg/cell so div by area for kg/m2 and x10 for ton/ha */

        op.DrawMap3->copy(tmb);  //soilloss in ton/ha
    }
    if (SwitchChannelFlood)
    {
        FOR_ROW_COL_MV
                tmb->Drc = hmx->Drc < 0.01 ? 0 : hmx->Drc;
        op.DrawMap4->copy(tmb);  //flood level in m
    }

    op.baseMap->copy(Shade);
    if (SwitchIncludeChannel)
        op.channelMap->copy(ChannelWidth);
    if (SwitchRoadsystem)
        op.roadMap->copy(RoadWidthDX);

    op.dx = _dx;
    op.MB = MB;
    op.runstep = runstep;
    op.maxstep = (int) ((EndTime-BeginTime)/_dt);
    op.EndTime = EndTime/60.0;
    //	op.BeginTime = BeginTime/60.0;

    op.CatchmentArea = CatchmentArea;

    op.RainTotmm = RainTotmm + SnowTotmm;
    op.WaterVolTotmm = WaterVolTotmm-SurfStoremm;
    op.Qtotmm = Qtotmm;
    op.Qtot = QtotOutlet;
    op.QPlot = QPlot;  //VJ 110701
    op.QtotPlot = QtotPlot;  //VJ 110701
    op.QpeakPlot = QpeakPlot;  //VJ 110701
    op.SoilLossTotPlot = SoilLossTotPlot;
    op.Qpeak = Qpeak;
    op.QpeakTime = QpeakTime/60;
    op.RainpeakTime = RainpeakTime/60;

    op.InfilTotmm = InfilTotmm;
    op.SurfStormm = SurfStoremm;
    op.IntercTotmm = IntercTotmm;// + IntercHouseTotmm;
    //houses
    op.IntercHouseTotmm = IntercHouseTotmm;
    op.InfilKWTotmm = InfilKWTot; // infil part in kin wave not used
    op.RunoffFraction = (RainTotmm > 0 ? Qtotmm/RainTotmm : 0);

    op.MBs = MBs;
    op.DetTotSplash=DetSplashTot*0.001; // convert from kg to ton
    op.DetTotFlow=DetFlowTot*0.001; // convert from kg to ton
    op.DepTot=DepTot*0.001; // convert from kg to ton
    op.SedTot=SedTot*0.001; // convert from kg to ton

    op.ChannelDetTot=ChannelDetTot*0.001; // convert from kg to ton
    op.ChannelDepTot=ChannelDepTot*0.001; // convert from kg to ton
    op.ChannelSedTot=ChannelSedTot*0.001; // convert from kg to ton
    op.ChannelWH = ChannelWH->DrcPlot;

    op.SoilLossTot=SoilLossTotOutlet*0.001; // convert from kg to ton

    op.t = time_ms.elapsed()*0.001/60.0;
    op.time = time/60;
    op.maxtime = op.t/runstep * op.maxstep;

    op.Pmm = (RainAvgmm + SnowAvgmm)*3600/_dt;
    op.Q = Qoutput->DrcPlot; //Outlet;  //=> includes channel and tile
    op.Qs = Qsoutput->DrcPlot; //Outlet;
    op.C = TotalConc->DrcPlot; //Outlet;
    op.Qtile = 1000*TileQn->DrcPlot; //Outlet;
    // VJ 110630 show hydrograph for selected output point
    op.volFloodmm = floodTotmm;

    op.BufferVolTot = BufferVolin;//Tot;
    op.BufferSedTot = -BufferSedTot*0.001; // convert from kg to ton, negative beause is deposition
}
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
    // reportstats per landunit class
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
    int hour, min, sec;
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
                    out.setRealNumberPrecision(3);
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
                            out << ",l/s" << ",m";
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
            out.setRealNumberPrecision(3);
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
                newname1 = fi.path() + "/" + fi.baseName() + "_" +
                        QString::number((int)PointMap->Drc) + "." +  fi.suffix();
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
                    if (SwitchRainfall) out << sep << RainIntavg;
                    if (SwitchSnowmelt) out << sep << SnowIntavg;
                    out << sep << Qoutput->Drc << sep << ChannelWH->Drc;
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
        out.setRealNumberPrecision(3);
        out.setFieldWidth(width);
        out.setRealNumberNotation(QTextStream::FixedNotation);

        if (!SwitchSOBEKoutput)
        {
            if (SwitchWritePCRtimeplot)
                out << runstep;
            else
                out << time/60;

            FOR_ROW_COL_MV
            {
                if ( PointMap->Drc > 0 )
                {
                    if (SwitchRainfall) out << sep << RainIntavg;
                    if (SwitchSnowmelt) out << sep << SnowIntavg;
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
    out.setRealNumberPrecision(5);
    out.setFieldWidth(12);
    out.setRealNumberNotation(QTextStream::FixedNotation);
    out << "LISEM run with:," << op.runfilename << "\n";
    out << "LISEM results at time (min):," << op.time <<"\n";
    out << "---------------------------------------------\n";
    out << "Catchment area              (ha):," << op.CatchmentArea/10000.0<< "\n";
    out << "Total rainfall              (mm):," << op.RainTotmm<< "\n";
    out << "Total discharge             (mm):," << op.Qtotmm<< "\n";
    out << "Total interception          (mm):," << op.IntercTotmm<< "\n";
    out << "Total House interception    (mm):," << op.IntercHouseTotmm<< "\n";
    out << "Total infiltration          (mm):," << op.InfilTotmm<< "\n";
    out << "Surface storage             (mm):," << op.SurfStormm<< "\n";
    out << "Water in runoff + channel   (mm):," << op.WaterVolTotmm<< "\n";
    out << "Total discharge             (m3):," << op.Qtot<< "\n";
    out << "Peak discharge             (l/s):," << op.Qpeak<< "\n";
    out << "Peak time rainfall         (min):," << op.RainpeakTime<< "\n";
    out << "Peak time discharge        (min):," << op.QpeakTime<< "\n";
    out << "Discharge/Rainfall           (%):," << op.RunoffFraction*100<< "\n";
    out << "Splash detachment (land)   (ton):," << op.DetTotSplash<< "\n";
    out << "Flow detachment (land)     (ton):," << op.DetTotFlow<< "\n";
    out << "Deposition (land)          (ton):," << op.DepTot<< "\n";
    out << "Suspended Sediment (land)  (ton):," << op.SedTot<< "\n";
    out << "Flow detachment (channels) (ton):," << op.ChannelDetTot<< "\n";
    out << "Deposition (channels)      (ton):," << op.ChannelDepTot<< "\n";
    out << "Susp. Sediment (channels)  (ton):," << op.ChannelSedTot<< "\n";
    out << "Susp. Sediment (buffers)   (ton):," << op.BufferSedTot<< "\n";
    out << "Total soil loss            (ton):," << op.SoilLossTot<< "\n";
    out << "Average soil loss        (kg/ha):," << (op.SoilLossTot*1000.0)/(op.CatchmentArea/10000.0)<< "\n";

    fp.flush();
    fp.close();
}
//---------------------------------------------------------------------------
/// Report maps for totals and mapseries (like report in PCRaster)
/// output filenames are fixed, cannot be changed by the user
void TWorld::ReportMaps(void)
{
    //   if (units == 0)
    //      checkUnits_tonha->setChecked(true);
    //   if (units == 1)
    //      checkUnits_kgcell->setChecked(true);
    //   if (units == 2)
    //      checkUnits_kgm2->setChecked(true);
    if(SwitchErosion)
    {
        // VJ 110111 erosion units
        tm->copy(TotalDetMap); //kg/cell
        if (ErosionUnits == 2)  // in kg/m2
            tm->calcMap(CellArea, DIV);
        if (ErosionUnits == 0) // ton/ha
        {
            tm->calcMap(CellArea, DIV); //to kg/m2
            tm->calcValue(10, MUL); // * 0.001*10000 = ton/ha
        }

        tm->report(totalErosionFileName);
        if (outputcheck[5].toInt() == 1)
            tm->report(Outeros); // in units

        tm->copy(TotalDepMap); //kg/cell
        if (ErosionUnits == 2)  // in kg/m2
            tm->calcMap(CellArea, DIV);
        if (ErosionUnits == 0) // ton/ha
        {
            tm->calcMap(CellArea, DIV);
            tm->calcValue(10, MUL);
        }
        tm->report(totalDepositionFileName);
        if (outputcheck[6].toInt() == 1)
            tm->report(Outdepo); // in units

        tm->copy(TotalSoillossMap); //kg/cell
        if (ErosionUnits == 2)  // in kg/m2
            tm->calcMap(CellArea, DIV);
        if (ErosionUnits == 0) // ton/ha
        {
            tm->calcMap(CellArea, DIV);
            tm->calcValue(10, MUL);
        }
        tm->report(totalSoillossFileName);

        if (outputcheck[1].toInt() == 1) Conc->report(Outconc);  // in g/l
        if (outputcheck[4].toInt() == 1) TC->report(Outtc);      // in g/l
    }

    if (outputcheck[0].toInt() == 1)
        Qoutput->report(Outrunoff); // in l/s
    if (outputcheck[2].toInt() == 1)
    {
        tm->calcMapValue(WH, 1000, MUL);// WH in mm
        tm->report(Outwh);
    }
    if (outputcheck[3].toInt() == 1)
        WHrunoffCum->report(Outrwh); // in mm

    if (outputcheck[7].toInt() == 1) V->report(Outvelo);
    FOR_ROW_COL_MV
    {
        InfilVolCum->Drc += InfilVol->Drc + InfilVolKinWave->Drc;
        InfilmmCum->Drc = max(0, InfilVolCum->Drc*1000/CellArea->Drc);
    }
    if (outputcheck[8].toInt() == 1) InfilmmCum->report(Outinf); // in mm

    if (outputcheck[9].toInt() == 1)
    {
        tm->calcMapValue(WHstore, 1000, MUL);// in mm
        tm->report(Outss);
        /** TODO check this: surf store in volume m3 is multiplied by flowwidth? */
    }

    if (outputcheck[10].toInt() == 1) ChannelWaterVol->report(Outchvol);


    if (SwitchIncludeTile && outputcheck.count() > 11)
    {
        if (outputcheck[11].toInt() == 1)
        {
            tm->calcMapValue(TileQn, 1000, MUL);// in mm
            tm->report(OutTiledrain);
        }
    }
    if (SwitchChannelFlood)
    {
        if (outputcheck.count() > 12)
        {
            if (outputcheck[12].toInt() == 1)
                hmx->report(OutHmx);
            if (outputcheck[13].toInt() == 1)
            {
                Qflood->report(OutQf);
            }
            if (outputcheck[14].toInt() == 1)
            {
                UVflood->report(OutVf);
            }
        }
    }
    /* from old LISEM: order in run file
   char *q = strtok(p,",");SwitchMapoutRunoff= strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutConc  = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutWH    = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutWHC   = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutTC    = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutEros  = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutDepo  = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutV     = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutInf   = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutSs    = strcmp(q,"1") == 0;
   q = strtok(NULL,",");   SwitchMapoutChvol = strcmp(q,"1") == 0;
 */

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
        unitList[i].area = 0;
        unitList[i].totdet = 0;
        unitList[i].totdep = 0;
        unitList[i].totsl = 0;
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
        unitList[i].area = 0;
        unitList[i].totdet = 0;
        unitList[i].totdep = 0;
        unitList[i].totsl = 0;
    }

    FOR_ROW_COL_MV
    {
        //variables are kg/cell convert to ton/cell
        for (int i = 0; i < landUnitNr; i++)
            if (unitList[i].nr == (long)LandUnit->Drc)
            {
                unitList[i].area += CellArea->Drc/10000;//ha
                unitList[i].totdet += TotalDetMap->Drc/1000; //ton/cell
                unitList[i].totdep += TotalDepMap->Drc/1000;
                unitList[i].totsl += TotalSoillossMap->Drc/1000;
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
        out << unitList[i].nr << unitList[i].area << unitList[i].totdet
            << unitList[i].totdep << unitList[i].totsl << "\n";
    fout.close();

}
//---------------------------------------------------------------------------
