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
/// report to file: timeseries at output points, totals, map series and land unit stats
void TWorld::reportToFile(void)
{
    ReportTotalsNew();
    // report totals to a text file

    if (SwitchWritePCRtimeplot)
        ReportTimeseriesPCR();
    else
        ReportTimeseriesCSV();
    // report hydrographs ande sedigraphs at all points in outpoint.map

    ReportTotalSeries();
    // report catchment averages per timestep

    // spatial output, maps and mapseries
    // savemaptodisk reacts to printinterval
    if(savemaptodisk) {
        ReportMaps();
        ReportMapSeries();
    }
    // report all maps and mapseries

    ReportErosionLandunits();
    // report stats per landunit class

    FloodStatistics();
    // report buildings submerged in flood level classes in 5cm intervals
}
//---------------------------------------------------------------------------
void TWorld::setupHydrographData()
{
    // clear first
    op.OutletIndices.clear();
    op.OutletLocationX.clear();
    op.OutletLocationY.clear();
    op.OutletQ.clear();
    op.OutletQs.clear();
    op.OutletC.clear();
    op.Qbound.clear();
    op.OutletQpeak.clear();
    op.OutletQpeaktime.clear();
    op.OutletChannelWH.clear();
    op.OutletQtot.clear();
    op.OutletQstot.clear();

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
   // tx.clear();
    tx.append(op.OutletLocationX);
   // ty.clear();
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
/** fill output structure 'op' with results to talk to the interface:
    report to screen, hydrographs */
void TWorld::reportToUI(void)
{
    SwitchCorrectMB_WH = op.SwitchCorrectMB_WH;
    op.timestep = this->_dt/60.0;

    op.t = time_ms.elapsed()*0.001/60.0;
    op.t = omp_get_wtime()/60.0 - startTime;
    op.time = time/60.0; // current time in min
    if (SwitchEventbased)
        op.Time.append(time/60.0);  // vector of time in min
    else
        op.Time.append(time/86400.0); // vector of time in days
    op.maxtime = op.t/runstep * op.maxstep;
    op._dx = _dx;
    op._llx = _llx;
    op._lly = _lly;
    op._nrCols = _nrCols;
    op._nrRows = _nrRows;
    op.runstep = runstep;
    op.maxstep = (int) ((EndTime-BeginTime)/_dt_user);
    //op.EndTime = EndTime/60.0;
    op.CatchmentArea = CatchmentArea;

    op.RainTotmm = RainTotmm;// + SnowTotmm;
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
        op.RunoffFraction = qMax(0.0, (op.Qtotmm - op.BaseFlowTotmm)/op.RainTotmm);
    op.WaterVolTotmm = WaterVolRunoffmm;
    op.StormDrainTotmm = StormDrainTotmm;
    op.ChannelVolTotmm = ChannelVolTotmm;
    op.BaseFlowTotmm = BaseFlowTotmm;
    op.PeakFlowTotmm = PeakFlowTotmm;
    op.RetentionVolTot = RetentionVolTot;
    op.RetentionVolTotmm = RetentionVolTotmm;

    op.FloodVolmm = floodVolTotmm;
    op.FloodTotMax = floodVolTotMax;
    op.FloodAreaMax = floodAreaMax;
    op.FloodArea = floodArea;

    op.Qtotmm = Qtotmm;
    op.Qboundtotmm = Qboundtotmm;
    op.Qtot = Qtot; // all outflow through channel and runoff for all open and outlets boundaries

    op.QBoundaryTot = QBoundaryTot;
    op.Qtiletot = QTiletot;  //total volume of water in tiles in m3
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
        op.SoilLossTot = SoilLossTot*0.001; // convert from kg to ton
        op.floodBoundarySedTot = floodBoundarySedTot*0.001; // not used

        op.OutletQs.at(0)->append(SoilLossTot_dt); //timestep output in kg! SoilLossOutlet = sum of Qs*dt and channelQs*dt and QsBoundary
        op.OutletC.at(0)->append(Qtot_dt > MIN_FLUX? SoilLossTot_dt/Qtot_dt : 0);
        op.OutletQstot.replace(0,SoilLossTot*0.001);
    }
    if (SwitchPest) {
        op.PMOutW = PestOutW;
        op.PMerr = PMerr;
        op.PMinf = Pestinf;
        //op.PMperc = PestPerc;
        op.PestName = PestName;
        op.PMtotI = PMtotI;
        if (SwitchErosion) {
            op.PMOutS = PestOutS;
        }
    }


    //hydrographs
    op.Pmm.append((RainAvgmm)*3600/_dt); // + SnowAvgmm

    // outlet 0 all flow
    op.OutletQ.at(0)->append(Qtot_dt * QUNIT/_dt); //Qtot_dt is in m3

    op.Qbound.append(QBoundary*QUNIT);
    op.Qtile.append(QTile*QUNIT);  //average tile output over all tile outlets as a flux in l/s

    op.OutletQtot.replace(0,Qtot); // cumulative tot outflow
    op.OutletChannelWH.at(0)->append(0);

    for(int j = 1; j < op.OutletIndices.length(); j++)
    {
        int r = op.OutletLocationX.at(j);
        int c = op.OutletLocationY.at(j);

        double channelwh = SwitchIncludeChannel? ChannelWH->Drc : 0.0;
        op.OutletChannelWH.at(j)->append(std::isnan(channelwh)?0.0:channelwh); //? why nan

        if (SwitchIncludeChannel) {
            op.OutletQtot.replace(j,op.OutletQtot.at(j) + _dt * (ChannelQn->Drc + QBoundary)); //cumulative in m3/s
            op.OutletQ.at(j)->append(ChannelQn->Drc*QUNIT);
        } else {
            op.OutletQtot.replace(j,op.OutletQtot.at(j) + _dt * (Qn->Drc + QBoundary)); //cumulative in m3/s
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
        double p = op.OutletQpeak.at(j);
        double q = op.OutletQ.at(j)->last();  //at(op.OutletQ.at(j)->length()-1); // this point last in list

        if(p < q) {
            op.OutletQpeak.replace(j,q);
            if (SwitchEventbased)
                op.OutletQpeaktime.replace(j,time/60-op.BeginTime);
            else
                op.OutletQpeaktime.replace(j,time/60);
           // qDebug() << time << op.BeginTime;
        }
    }
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
        if (!fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
            ErrorString = "Cannot open the result file: "+totalSeriesFileName;
            throw 1;
        }
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
        if (SwitchChannelBaseflowStationary)
            out << sep << "Baseflow in (mm)";
         if (SwitchGWflow)
            out << sep << "GWlevel (m)";

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
        if (SwitchPest) {
            out << sep << "PMOutW";
            out << sep << "PMerr";
            out << sep << "PestPerc";
            out << sep << "Pestinf";
            if (SwitchErosion) {
                out << sep << "PMOutS";
            }
        }
        out << "\n";
        fout.flush();
        fout.close();
    }

    QFile fout(newname1);
    if (!fout.open(QIODevice::Append | QIODevice::Text)) {
        ErrorString = "Cannot open the result file: "+totalSeriesFileName;
        throw 1;
    }

    QTextStream out(&fout);
    out.setRealNumberPrecision(DIG);
    out.setFieldWidth(width);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    out << (time-BeginTime)/60;
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
    if (SwitchChannelBaseflowStationary)
        out << sep << op.BaseFlowTotmm;
    if (SwitchGWflow)
        out << sep << op.GWlevel;

    if (SwitchIncludeStormDrains)
        out << sep << op.StormDrainTotmm;
    out << sep << op.WaterVolTotmm;
    out << sep << op.FloodVolmm;
    out << sep << op.FloodArea;
    out << sep << op.ChannelVolTotmm;
    out << sep << op.Qtotmm;
    if (FlowBoundaryType > 0)
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
    if (SwitchPest) {
        out << sep << op.PMOutW;
        out << sep << op.PMerr;
        out << sep << op.PMperc;
        out << sep << op.PMinf;
        if (SwitchErosion) {
            out << sep << op.PMOutS;
        }

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
       out << QString("\"Water in overland flow (h<%1) (mm):\",%2\n").arg(minReportFloodHeight*1000).arg(op.WaterVolTotmm);
       out << QString("\"Water in flood (h>%1) (mm):\",%2\n").arg(minReportFloodHeight*1000).arg(op.FloodVolmm);
    }
    out << "\"Water in channels (mm):\"," << op.ChannelVolTotmm<< "\n";
    out << "\"Water across boundary (mm):\"," << op.Qboundtotmm<< "\n";
    out << "\"Water in rentention (m3):\"," << op.RetentionVolTot << "\n";
    out << "\"Total baseflow and GW inflow (mm):\"," << op.BaseFlowTotmm << "\n";
    out << "\"Total peakflow (mm):\"," << op.PeakFlowTotmm << "\n";
    out << "\"Total outflow (overland+channel) (mm):\"," << op.Qtotmm << "\n";
    out << "\"Total outflow (overland+channel) (m3):\"," << op.Qtot<< "\n";
    out << "\"Total boundary outflow (m3):\"," << op.QBoundaryTot<< "\n";
    out << "\"Total storm/tile drain discharge (m3):\"," << op.Qtiletot<< "\n";
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
    if (SwitchPest) {
        out << "\n";
        out << "\"Pesticide simulated:\"," << op.PestName << "\n";
        out << "\"Initial pesticide mass in system (mg) \", " << op.PMtotI << "\n";
        out << "\"Total dissolved pesticide transport (mg):\"," << op.PMOutW<< "\n";
        if (SwitchErosion) out << "\"Total particulate pesticide transport (mg):\"," << op.PMOutS<< "\n";
    }
    out << "\n";
    fp.flush();
    fp.close();
}

//---------------------------------------------------------------------------
void TWorld::ReportTimeseriesPCR(void)
{
    int nr = 0;

    int DIG = ReportDigitsOut;
    //int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
    double RainIntavg = RainAvgmm * 3600/_dt;
    //double SnowIntavg = SnowAvgmm * 3600/_dt;
    QString newname1, sep = " ";
    int width = 3+DIG-3;


    double QALL = Qtot_dt * QUNIT/_dt; // total outflow for all outlets, same as point 0 in interface, and all boundary, everuthing!
    double QSALL = SoilLossTot_dt/_dt; //total sed loss in kg/s from all outlets, surface and boundary

    QFileInfo fi(resultDir + outflowFileName);

    //######  open files and write headers #####//

    QString unitS = "l/s";
    if (QUnits == 1)
        unitS = "m3/s";

    // switchwriteheadersis done in report totals
    if (SwitchWriteHeaders) //  make file at first timestep
    {
        FOR_ROW_COL_MV_OUTL {
            newname1 = fi.path() + "/" + fi.baseName() + "_" + crout_[i_].code + "." +  fi.suffix();
            // make filename using point number

            QFile fout(newname1);
            if (!fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
                ErrorString = "Cannot open the result file: "+newname1;
                throw 1;
            }
            QTextStream out(&fout);
            out.setRealNumberPrecision(DIG);
            out.setFieldWidth(width);
            out.setRealNumberNotation(QTextStream::FixedNotation);

            // HEADERS for the 3 types
            if (SwitchWritePCRtimeplot)  //PCRaster timeplot format, cannot be SOBEK !
            {

                out << "#LISEM flow and sed output file for point #"+crout_[i_].code+"\n";

                // nr columns is time + rain + Q + (maybe Qs + C)
                int nrs = 5 + (SwitchErosion ? 3 : 0);
                if (SwitchRainfall) nrs++;
                if (SwitchSnowmelt) nrs++;
                //if (SwitchChannelBaseflowStationary || SwitchGWflow) nrs++;
                if (FlowBoundaryType > 0) nrs++;
                if (SwitchIncludeTile) nrs++;
                    out << nrs << "\n";

                out << "run step\n";
                out << "time (day)\n";
                if (SwitchRainfall) out << "Pavg (mm/h)\n";
                if (SwitchSnowmelt) out << "Snowavg (mm/h)\n";
                out << "Qall" << unitS << "\n";
                if (FlowBoundaryType > 0)
                    out << "QBound " << unitS << "\n";
                if (SwitchIncludeChannel) {
                    out << "Qchan"+crout_[i_].code << unitS;
                    out << "\n" << "WHchan" + crout_[i_].code + " (m)\n";
                } else {
                    out << "Qof " << unitS << "\n";
                }
                if (SwitchIncludeTile) out << "Qdrain (l/s)\n";
                if (SwitchErosion) {
                    out << "Qsall (kg/s)\n";
                    if (FlowBoundaryType > 0)
                        out << "QsBoundary (kg/s)";
                    if (SwitchIncludeChannel)
                        out << "Qschan%1"+crout_[i_].code +" (kg/s)\n";
                     else
                        out << "Qsof (kg/s)\n";
                    out << "C (g/l)\n";
                }

            }
             fout.close();
        }}
    }  // opening files and writing header

    //######  open files and append values #####//
    // for all outlet points
    FOR_ROW_COL_MV_OUTL
    {
        newname1 = fi.path() + "/" + fi.baseName() + "_" + crout_[i_].code + "." +  fi.suffix();

        QFile fout(newname1);
        if (!fout.open(QIODevice::Append | QIODevice::Text)) {
            ErrorString = "Cannot open the result file: "+newname1;
            throw 1;
        }

        QTextStream out(&fout);
        out.setFieldWidth(width);
        out.setRealNumberNotation(QTextStream::FixedNotation);
        out.setRealNumberPrecision(5);
        out << runstep;
        out << sep << (time/60)/1440.0;

        out.setRealNumberPrecision(DIG);
        if (SwitchRainfall) out << sep << RainIntavg;
        //if (SwitchSnowmelt) out << sep << SnowIntavg;

        out << sep << QALL;
        if (FlowBoundaryType > 0)
            out << sep << QBoundary*QUNIT;

        if (SwitchIncludeChannel) {
            out << sep << ChannelQn->Drc*QUNIT;
            out << sep << ChannelWH->Drc;
        } else {
            out << sep << Qn->Drc*QUNIT;
        }

        if (SwitchIncludeTile)
            out << sep << TileQn->Drc*QUNIT;

        if (SwitchErosion) {
            out << sep << QSALL;
            out << sep << QsBoundary;
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

}
//---------------------------------------------------------------------------
void TWorld::ReportTimeseriesCSV(void)
{
    int DIG = ReportDigitsOut;

    double RainIntavg = RainAvgmm * 3600/_dt;
    double SnowIntavg = SnowAvgmm * 3600/_dt;
    QString newname1, sep = ",";
    int width = 0;
    double QALL = Qtot_dt * QUNIT/_dt; // total outflow for all outlets, same as point 0 in interface
    double QSALL = SoilLossTot_dt/_dt; //total sed loss in kg/s from all outlets, surface and boundary

    QFileInfo fi(resultDir + outflowFileName);

    QString unitS = "l/s";
    if (QUnits == 1)
        unitS = "m3/s";

    //######  open files and write headers #####//

    //SwitchWriteHeaders is set to false in ReportTotalSeries(void)!

    if (SwitchWriteHeaders) {
        FOR_ROW_COL_MV_OUTL {
            newname1 = fi.path() + "/" + fi.baseName() + "_" + crout_[i_].code + "." +  fi.suffix();

            // make filename using point number

            QFile fout(newname1);
            if (!fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
                ErrorString = "Cannot open the result file: "+newname1;
                throw 1;
            }
            QTextStream out(&fout);
            out.setRealNumberPrecision(DIG);
            out.setFieldWidth(width);
            out.setRealNumberNotation(QTextStream::FixedNotation);

            out << "LISEM total flow and sed output file for point " + crout_[i_].code + "\n";
            // first row, variable names
            out << "Time";
            if (SwitchRainfall) out << ",Pavg";
            if (SwitchSnowmelt) out << ",Snowavg";
            out << ",Qall";
            if (FlowBoundaryType > 0)
                out << ",Qbound";
            if (SwitchIncludeChannel) {
                out << ",Qchan" + crout_[i_].code;
                out << ",WHchan" + crout_[i_].code;
            } else {
                out << ",Qrunoff";
            }
            if (SwitchIncludeTile) out << ",Qtile";
            if (SwitchErosion){
                    out << ",Qsall";
                if (FlowBoundaryType > 0)
                    out << ",Qsbound";
                if (SwitchIncludeChannel)
                    out << ",Qschan" + crout_[i_].code;
                else
                    out << QString(",Qsrunoff");
                out << ",Conc";
            }
            out << "\n";

            // second row, units
            out << "days"; //time
            if (SwitchRainfall) out << ",mm/h"; //rain
            if (SwitchSnowmelt) out << ",mm/h"; // snow
            out << "," << unitS; // qall
            if (FlowBoundaryType > 0)
                out << "," << unitS; //qbound
            if (SwitchIncludeChannel)
                out  << "," << unitS << ",m"; //qchannel
            else
                out  << "," << unitS; // Orunoff
            if (SwitchIncludeTile)
                out << "," << unitS;
            if (SwitchErosion) {
                out << ",kg/s";
                if (FlowBoundaryType > 0)
                    out << ",kg/s";
                out<< ",kg/s" << ",g/l";
            }
            out << "\n";
            fout.close();
        }}


    }  // opening files and writing header

    //######  open files and append values #####//
    // for all outlet points
    FOR_ROW_COL_MV_OUTL {
        newname1 = fi.path() + "/" + fi.baseName() + "_" + crout_[i_].code + "." +  fi.suffix();
        QFile fout(newname1);
        if (!fout.open(QIODevice::Append | QIODevice::Text)) {
            ErrorString = "Cannot open the result file: "+newname1;
            throw 1;
        }
        QTextStream out(&fout);
        out.setFieldWidth(width);
        out.setRealNumberNotation(QTextStream::FixedNotation);
        out.setRealNumberPrecision(7);

        out << (time/60)/1440.0;

        out.setRealNumberPrecision(DIG);
        if (SwitchRainfall) out << sep << RainIntavg;
        if (SwitchSnowmelt) out << sep << SnowIntavg;

        out << sep << QALL; // all water

        if (FlowBoundaryType > 0)
            out << sep << QBoundary*QUNIT;

        if (SwitchIncludeChannel) {
            out << sep << ChannelQn->Drc*QUNIT;
            out << sep << ChannelWH->Drc;
        } else {
            out << sep << Qn->Drc*QUNIT; // overlandflow
        }

        if (SwitchIncludeTile)
            out << sep << TileQn->Drc*QUNIT;

        if (SwitchErosion) {
            out << sep << QSALL;
            if (FlowBoundaryType > 0)
                out << sep << QsBoundary;
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
/// Report the erosion totals per land unit
void TWorld::ReportErosionLandunits(void)
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
             //   unitList[i].var1 += qMax(0.0,TotalSoillossMap->Drc/1000); //ton/cell
             //   unitList[i].var2 += qMin(0.0,TotalSoillossMap->Drc/1000);
                unitList[i].var1 += TotalSoillossMap->Drc/1000;
            }
    }}


    QString name;
    name = resultDir + totalLandunitFileName;//QFileInfo(totalLandunitFileName).baseName()+"-"+op.timeStartRun+".csv";
    QFile fout(name);
    if (!fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
        ErrorString = "Cannot open the result file: "+name;
        throw 1;
    }
    QTextStream out(&fout);
    out.setRealNumberPrecision(3);
    out.setRealNumberNotation(QTextStream::FixedNotation);

    // out << "Landunit,Area,Detachment,Deposition,Soil Loss\n";
    // out << "#,ha,ton,ton,ton\n";
    out << "Landunit,Area,Soil Loss\n";
    out << "#,ha,ton\n";
    for (int i = 0; i < landUnitNr; i++)
        out << unitList[i].nr << ","
            << unitList[i].var0 << ","
            << unitList[i].var1 << "\n";
          //  << unitList[i].var2 << ","
          //  << unitList[i].var3 << "\n";
    fout.close();

}
//---------------------------------------------------------------------------
void TWorld::FloodStatistics(void)
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
            nr = qMax(nr, i);
            //qDebug() << nr << i << floodHmxMax->Drc;
            floodList[i].var1 += area; // area flooded in this class
            floodList[i].var2 += area*floodHmxMax->Drc; // vol flooded in this class
            floodList[i].var3 = qMax(floodTime->Drc/60.0,floodList[i].var3); // max time in this class
            floodList[i].var4 = qMax(floodTimeStart->Drc/60.0,floodList[i].var4); // max time in this class
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
