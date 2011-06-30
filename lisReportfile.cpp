
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
- void TWorld::reportAll() \n
- void TWorld::ReportTimeseriesNew() \n
- void TWorld::ReportTotalsNew() \n
- void TWorld::ReportMaps() \n
- void TWorld::CountLandunits() \n
- void TWorld::ReportLandunits() \n
 */

#include "lisemqt.h"
#include "model.h"
#include "global.h"


//---------------------------------------------------------------------------
/// report to disk: timeseries at output points, totals, map series and land unit stats
void TWorld::reportAll()
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
void TWorld::ReportTimeseriesNew()
{
	int nr = 0;
	int hour, min, sec;
	int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
	double RainIntavg = RainAvgmm * 3600/_dt;
	double SnowIntavg = SnowAvgmm * 3600/_dt;
	QString newname1, pnr, sep = (SwitchWritePCRtimeplot ? " " : ",");
	int width = (!SwitchWritePCRtimeplot ? 0 : 9);

	QFileInfo fi(resultDir + outflowFileName);

   // - TODO: if runs are interrupted the number of lines win the SOBEK output will not be correct!
   if (SwitchSOBEKOutput && time > 0)
	{
		sec = (int)(time*60);
		hour = (int)(sec/3600);
		min = (int)(time - hour*60);
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
						out << "Q (l/s)\n";
						if (SwitchErosion) out << "Qs (kg/s)\n";
						if (SwitchErosion) out << "C (g/l)\n";
					}
					else // SOBEK format
						if (SwitchSOBEKoutput)
						{
                     pnr.setNum(SOBEKlines);
                     out << "Q";
                     if (SwitchErosion) out << " Qs C";
                     out << "\n";
                     out << "m3/s";
                     if (SwitchErosion) out << " kg/s g/l";
                     out << "\n";
                     out << pnr << "\n";
                  }
                  else // flat format, comma delimited
                  {
                     pnr.setNum((int)PointMap->Drc);
                     out << "LISEM total flow and sed output file for point " << pnr << "\n";
                     out << "Time";
                     if (SwitchRainfall) out << ",Pavg";
                     if (SwitchSnowmelt) out << ",Snowavg";
                     if (SwitchErosion) out << ",Qs,C";
                     out << "\n";
                     out << "min,mm/h";
                     if (SwitchSnowmelt) out << ",mm/h";
                     out << ",l/s";
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
					if (SwitchErosion) out << "Qs #"<< pnr << "(kg/s)\n";
					if (SwitchErosion) out << "C #"<< pnr << "(g/l)\n";
				}
				out << "\n";
			}
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
					out << sep << Qoutput->Drc;
					if (SwitchErosion) out << sep << Qsoutput->Drc;
					if (SwitchErosion) out << sep << TotalConc->Drc;
					out << "\n";
				}
				else  //SOBEK format
				{
					out.setFieldWidth(2);
					out << "\"" << SOBEKdatestring << ":" << hour << ":" <<  min << ":" <<  sec;
					out.setFieldWidth(8);
					out << Qoutput->Drc;
					if (SwitchErosion) out << " " << Qsoutput->Drc;
					if (SwitchErosion) out << " " << TotalConc->Drc;
					out << " <\n";
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
					out << sep << Qoutput->Drc;
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
                  out << " " << Qoutput->Drc;
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
void TWorld::ReportTotalsNew()
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
	out << "Total soil loss            (ton):," << op.SoilLossTot<< "\n";
	out << "Average soil loss        (kg/ha):," << op.SoilLossTot/op.CatchmentArea*10000.0*1000<< "\n";

	fp.flush();
	fp.close();
}
//---------------------------------------------------------------------------
/// Report maps for totals and mapseries (like report in PCRaster)
/// output filenames are fixed, cannot be changed by the user
void TWorld::ReportMaps()
{
	if(SwitchErosion)
	{
      // VJ 110111 erosion units
      tm->copy(TotalDetMap); //kg/cell
      if (ErosionUnits < 2)  // in kg/m2
         tm->calc(CellArea, DIV);
      if (ErosionUnits == 0) // ton/ha
         tm->calcV(10, MUL);
      tm->mwrite(totalErosionFileName);
      if (outputcheck[5].toInt() == 1)
         tm->report(Outeros); // in units

      tm->copy(TotalDepMap); //kg/cell
      if (ErosionUnits < 2)  // in kg/m2
         tm->calc(CellArea, DIV);
      if (ErosionUnits == 0) // ton/ha
         tm->calcV(10, MUL);
      tm->mwrite(totalDepositionFileName);
      if (outputcheck[6].toInt() == 1)
         tm->report(Outdepo); // in units

      tm->copy(TotalSoillossMap); //kg/cell
      if (ErosionUnits < 2)  // in kg/m2
         tm->calc(CellArea, DIV);
      if (ErosionUnits == 0) // ton/ha
         tm->calcV(10, MUL);
      tm->mwrite(totalSoillossFileName);

      if (outputcheck[1].toInt() == 1) Conc->report(Outconc);  // in g/l
		if (outputcheck[4].toInt() == 1) TC->report(Outtc);      // in g/l
	}

	if (outputcheck[0].toInt() == 1) Qoutput->report(Outrunoff); // in l/s
	if (outputcheck[2].toInt() == 1)
	{
		tm->calc2V(WH, 1000, MUL);// WH in mm
		tm->report(Outwh);
	}
	if (outputcheck[3].toInt() == 1)	WHrunoffCum->report(Outrwh); // in mm

	if (outputcheck[7].toInt() == 1) V->report(Outvelo);
	FOR_ROW_COL_MV
	{
		InfilVolCum->Drc += InfilVol->Drc + InfilVolKinWave->Drc;
		tm->Drc = max(0, InfilVolCum->Drc*1000/CellArea->Drc);
	}
	if (outputcheck[8].toInt() == 1) tm->report(Outinf);

	if (outputcheck[9].toInt() == 1)
	{
		tm->calc2V(WHstore, 1000, MUL);// in mm
		tm->report(Outss);
      /** TODO check this: surf store in volume m3 is multiplied by flowwidth? */
	}

	if (outputcheck[10].toInt() == 1) ChannelWaterVol->report(Outchvol);

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
void TWorld::CountLandunits()
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
void TWorld::ReportLandunits()
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
      for (int i = 0; i < landUnitNr; i++)
         if (unitList[i].nr == (long)LandUnit->Drc)
         {
            unitList[i].area += CellArea->Drc/10000;
            unitList[i].totdet += TotalDetMap->Drc * CellArea->Drc/1000;
            unitList[i].totdep += TotalDepMap->Drc * CellArea->Drc/1000;
            unitList[i].totsl += TotalSoillossMap->Drc * CellArea->Drc/1000;
         }
   }
   // in ton

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
