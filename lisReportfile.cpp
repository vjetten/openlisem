/*
 * reportfile.cpp
 *
 *  Created on: Mar 2, 2010
 *      Author: jetten
 */

#include "ifacebasic.h"
#include "model.h"
#include "global.h"
//#include <qstring.h>


//---------------------------------------------------------------------------
// reporting timeseries for every non zero point PointMap
// - 3 types of output: PCRaster timeplot format; SOBEK input format; flat comma delimited format
// - all points in one file or each point in a separate file
// - the types should be mututally exclusive in the interface and run file
// - TO BE DONE: if runs are interrupted the number of lines win the SOBEK output will not be correct!
void TWorld::ReportTimeseriesNew()
{
   int nr = 0;
   int hour, min, sec;
   int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
   QString newname1, pnr;
   QFileInfo fi(resultDir + outflowFileName);

   if (SwitchSOBEKOutput && time > 0)
   {
      sec = (int)(time*60);
      hour = (int)(sec/3600);
      min = (int)(time - hour*60);
      sec = (int)(sec - hour * 3600 - min * 60);
   }

   //SOBEK, PCRaster and flat format are mutually exclusive
   if (time <= BeginTime) //  make file at first timestep
   {
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
             out.setFieldWidth(8);
             out.setRealNumberNotation(QTextStream::FixedNotation);

             // HEADERS for the 3 types
             if (SwitchWritePCRtimeplot)  //PCRaster timeplot format, cannot be SOBEK !
             {
            	 pnr.setNum((int)PointMap->Drc);
                 out << "#LISEM total flow and sed output file for point " << pnr << "\n";
                 SwitchErosion ? out << "5\n" : out << "3\n";
                 out << "run step\n";
                 out << "Pavg (mm)\n";
                 out << "Q (l/s)\n";
                 if (SwitchErosion) out << "Qs (kg/s)\n";
                 if (SwitchErosion) out << "C conc (g/l)\n";
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
                 out << "#LISEM total flow and sed output file for point " << pnr << "\n";
                 out << "Time,Pavg";
                 if (SwitchErosion) out << ",Qs,C";
                 out << "\n";
                 out << "min,mm,l/s";
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
        out.setFieldWidth(8);
        out.setRealNumberNotation(QTextStream::FixedNotation);

        if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
        {
            // count nr of points
            FOR_ROW_COL_MV
              if ( PointMap->Drc > 0 )
                nr++;
            int nrs = 2+(1+(SwitchErosion ? 2 : 0))*nr;
			pnr.setNum(nrs);
            out << "#LISEM total flow and sed output file for all reporting points\n";
            out <<  pnr << "\n";
            out << "Time (min)\n";
            out << "Pavg (mm)\n";
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
           out << ",P";
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
           out << ",mm";
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

   //######  write all the values #####//
   //######  write all the values #####//
   //######  write all the values #####//

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
        out.setFieldWidth(8);
        out.setRealNumberNotation(QTextStream::FixedNotation);

        if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
        {
        	out << runstep;
        	out << Rain->Drc;
        	out << Qn->Drc*1000;
        	if (SwitchErosion) out << Qs->Drc;
        	if (SwitchErosion) out << SedVol->Drc;
        	out << "\n";
        }
        else  //SOBEK format
        if (SwitchSOBEKoutput)   //PCRaster timeplot and flat format
        {
            out.setFieldWidth(2);
        	out << "\"" << SOBEKdatestring << ":" << hour << ":" <<  min << ":" <<  sec;
            out.setFieldWidth(8);
        	out << Qn->Drc;
        	if (SwitchErosion) out << " " << Qs->Drc;
        	if (SwitchErosion) out << " " << SedVol->Drc;
        	out << " <\n";
        }
        else
        {
        	out << time/60;
        	out << "," << Rain->Drc;
        	out << "," << Qn->Drc*1000;
        	if (SwitchErosion) out << "," << Qs->Drc;
        	if (SwitchErosion) out << "," << SedVol->Drc;
        	out << "\n";
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
      out.setFieldWidth(8);
      out.setRealNumberNotation(QTextStream::FixedNotation);

      if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
      {
          out << runstep;
          FOR_ROW_COL_MV
          {
            if ( PointMap->Drc > 0 )
            {
               out << Rain->Drc;
               out << Qn->Drc*1000;
               if (SwitchErosion) out << Qs->Drc;
               if (SwitchErosion) out << SedVol->Drc;
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
               out << " " << Qn->Drc;
               if (SwitchErosion) out << " " << Qs->Drc;
               if (SwitchErosion) out << " " << SedVol->Drc;
            }
          }
          out << " < \n";
      }
      else // flat comma delimited
      {
          out << time/60;
          FOR_ROW_COL_MV
          {
            if ( PointMap->Drc > 0 )
            {
               out << "," << Rain->Drc;
               out << "," << Qn->Drc*1000;
               if (SwitchErosion) out << "," << Qs->Drc;
               if (SwitchErosion) out << "," << SedVol->Drc;
            }
          }
          out << "\n";
      }
      fout.close();
   }
}
//---------------------------------------------------------------------------
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
    out << "Average surface storage     (mm):," << op.SurfStorTotmm<< "\n";
    out << "Water in runoff             (mm):," << op.WaterVolTotmm<< "\n";
    out << "Total discharge             (m3):," << op.Qtot<< "\n";
    out << "Peak discharge             (l/s):," << op.Qpeak<< "\n";
    out << "Peak time rainfall         (min):," << op.RainpeakTime<< "\n";
    out << "Peak time discharge        (min):," << op.QpeakTime<< "\n";
    out << "Discharge/Rainfall           (%):," << op.RunoffFraction<< "\n";
    out << "Splash detachment (land)   (ton):," << op.DetTotSplash<< "\n";
    out << "Flow detachment (land)     (ton):," << op.DetTotFlow<< "\n";
    out << "Deposition (land)          (ton):," << op.DepTot<< "\n";
    out << "Suspended Sediment (land)  (ton):," << op.SedVolTot<< "\n";
    out << "Flow detachment (channels) (ton):," << op.ChannelDetTot<< "\n";
    out << "Deposition (channels)      (ton):," << op.ChannelDepTot<< "\n";
    out << "Susp. Sediment (channels)  (ton):," << op.ChannelSedTot<< "\n";
    out << "Total soil loss            (ton):," << op.SoilLossTot<< "\n";
    out << "Average soil loss        (kg/ha):," << op.SoilLossTot/op.CatchmentArea*10000.0*1000<< "\n";

    fp.flush();
    fp.close();
}
//---------------------------------------------------------------------------
void TWorld::ReportMaps()
{
    if (SwitchErosion)
    {
    	FOR_ROW_COL_MV
    	{
    		TotalDetMap->Drc += DETSplash->Drc + DETFlow->Drc;
    		TotalDepMap->Drc += DEP->Drc;
    		if (SwitchIncludeChannel)
    		{
    			TotalDetMap->Drc += ChannelDetFlow->Drc;
    			TotalDepMap->Drc += ChannelDep->Drc;
    		}
    		TotalSoillossMap->Drc = TotalDetMap->Drc + TotalDepMap->Drc;
    	}

    	TotalDetMap->mwrite(totalErosionFileName);
    	TotalDepMap->mwrite(totalDepositionFileName);
    	TotalSoillossMap->mwrite(totalSoillossFileName);
    }

    //Qoutput->copy(Qn);
    //Qoutput->calc(ChannelQn, ADD);
    //Qoutput->calcV(1000, MUL); // sum Qn and channel in l/s
    //Qoutput->report(Outrunoff);
    FOR_ROW_COL_MV
    {

        Qoutput->Drc = 1000*(Qn->Drc + ChannelQn->Drc);
        if (Outlet->Drc == 1)
        {
        	double oldpeak = Qpeak;
        	Qpeak = max(Qpeak, Qoutput->Drc);
        	if (oldpeak < Qpeak)
        		QpeakTime = time;
        }

    }
 //   Qoutput->report(Outrunoff);

//SedVol->report("sedvol");


    if (outputcheck[0].toInt() == 1) Qn->report(Outrunoff);
    if (outputcheck[1].toInt() == 1) Conc->report(Outconc); // in overland flow, not channels
    if (outputcheck[2].toInt() == 1) WH->report(Outwh);
    //if (outputcheck[3].toInt() == 1) WHrunoff->report(Outrwh); // CUMULATIVE RUNOFF?
    if (outputcheck[4].toInt() == 1) TC->report(Outtc);
    if (outputcheck[5].toInt() == 1) TotalDetMap->report(Outeros);
    if (outputcheck[6].toInt() == 1) TotalDepMap->report(Outdepo);
    if (outputcheck[7].toInt() == 1) V->report(Outvelo);
    if (outputcheck[8].toInt() == 1) fact->report(Outinf);
    if (outputcheck[9].toInt() == 1) WHstore->report(Outss);
    if (outputcheck[10].toInt() == 1) ChannelWaterVol->report(Outchvol);
    /*
    OUTRUNOFF=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\ro
    OUTCONC=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\conc
    OUTWH=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\wh
    OUTRWH=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\whc
    OUTTC=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\tc
    OUTEROS=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\det
    OUTDEPO=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\depo
    OUTVELO=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\velo
    OUTINF=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\inf
    OUTSS=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\sstor
    OUTCHVOL=D:\data\jantiene\Trier_17Oct03\Results_Trier_17Oct03\chanvol
*/


	//fpot->report("fpot",runstep);
		//fact->report("fact",runstep);
		//RainCum->report("rainc", runstep);
	    //Fcum->report("fcum", runstep);
	    //WHstore->report("sstor", runstep);
	    //Qn->report("Qn", runstep);
		//LeafDrain->report("ld", runstep);
		//SettlingVelocity->report("sv", runstep);

}
//---------------------------------------------------------------------------
/*
 * void TWorld::ReportTimeseries()
{
   FILE *fileout;
   int nr = 0;

   char SOBEKstr[11];
   int hour, min, sec;
   int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
   QString newname1;
   QFileInfo fi(resultDir + outflowFileName);

   if (SwitchSOBEKOutput && time > 0)
   {
      sec = (int)(time*60);
      hour = (int)(sec/3600);
      min = (int)(time - hour*60);
      sec = (int)(sec - hour * 3600 - min * 60);
   }
   if (SwitchSOBEKOutput)
       strcpy(SOBEKstr, SOBEKdatestring.toLatin1());

   //SOBEK, PCRaster and flat format are mutually exclusive
   if (time <= BeginTime) //  make file at first timestep
   {
      if (SwitchSeparateOutput) // each point separate file
      {
        FOR_ROW_COL_MV
        {
         if ( PointMap->Drc > 0 )
         {
//             sprintf(newname1,"%s_%d.%s",(const char*)fi.baseName().toLatin1(), (int)PointMap->Drc, (const char*)fi.suffix().toLatin1());
             newname1 = fi.path() + "/" + fi.baseName() + "_" +
            		    QString::number((int)PointMap->Drc) + "." +  fi.suffix();
             // make filename using point number

             //fileout = fopen(name_s, "w");
             QFile fout(newname1);
             fout.open(QIODevice::WriteOnly | QIODevice::Text);
             QTextStream out(&fout);
             out.setRealNumberPrecision(5);
             out.setFieldWidth(12);
             out.setRealNumberNotation(QTextStream::FixedNotation);
             //fileout = fopen(newname1.toAscii().constData(), "w");

             // HEADERS for the 3 types
             if (SwitchWritePCRtimeplot)  //PCRaster timeplot format, cannot be SOBEK !
             {
                 fprintf(fileout,"#LISEM total flow and sed output file for point %d\n", (int)PointMap->Drc);

                 if (SwitchErosion)
                   fprintf(fileout,"5\n");
                 else
                   fprintf(fileout,"3\n");
                 fprintf(fileout,"run step\n");
                 fprintf(fileout,"Pavg (mm)\n");
                 fprintf(fileout,"Q (l/s)\n");
                 if (SwitchErosion) fprintf(fileout,"Qs (kg/s)\n");
                 if (SwitchErosion) fprintf(fileout,"C conc (g/l)\n");
             }
             else // SOBEK format
             if (SwitchSOBEKoutput)
             {
                 fprintf(fileout,"Q");
                 if (SwitchErosion) fprintf(fileout," Qs");
                 if (SwitchErosion) fprintf(fileout," C");
                 fprintf(fileout,"\n");
                 fprintf(fileout,"m3/s");
                 if (SwitchErosion) fprintf(fileout," kg/s");
                 if (SwitchErosion) fprintf(fileout," g/l");
                 fprintf(fileout,"\n");
                 fprintf(fileout,"*%d\n",SOBEKlines);
             }
             else // flat format, comma delimited
             {
                 fprintf(fileout,"#LISEM total flow and sed output file for point %d\n", (int)PointMap->Drc);
                 fprintf(fileout,"Time");
                 fprintf(fileout,",Pavg");
                 if (SwitchErosion) fprintf(fileout,",Qs");
                 if (SwitchErosion) fprintf(fileout,",C");
                 fprintf(fileout,"\n");
                 fprintf(fileout,"min");
                 fprintf(fileout,",mm");
                 fprintf(fileout,",l/s");
                 if (SwitchErosion) fprintf(fileout,",kg/s");
                 if (SwitchErosion) fprintf(fileout,",g/l");
                 fprintf(fileout,"\n");
             }
             fclose(fileout);
           }
         }
      }  // separate
      else   // HEADERS: all points in one file
      {

        newname1 = resultDir + outflowFileName;
        // use simply resultdir + filename
//    	char *name_s = new char[newname1.length()+1];
//        strcpy(name_s, newname1.toAscii().constData());

//        fileout = fopen(name_s, "w");
        fileout = fopen(newname1.toAscii().constData(), "w");

//        delete[] name_s;

        if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
        {
            // count nr of points
            FOR_ROW_COL_MV
              if ( PointMap->Drc > 0 )
                nr++;

            fprintf(fileout,"#LISEM total flow and sed output file for all reporting points\n");
            fprintf(fileout,"%d\n",2+(1+(SwitchErosion ? 2 : 0))*nr);
            // number of variables = 2 + (1 or 2) * number of points
            fprintf(fileout,"Time (min)\n");
            fprintf(fileout,"Pavg (mm)\n");
            FOR_ROW_COL_MV
              if ( PointMap->Drc > 0 )
              {
                 fprintf(fileout,"Q #%d (l/s)\n",(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"Qs #%d (kg/s)\n",(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"C #%d (g/l)\n",(int)PointMap->Drc);
              }
            fprintf(fileout,"\n");
        }
        else // SOBEK format
        if (SwitchSOBEKoutput) //note: sobek input does not have rainfall
        {
          FOR_ROW_COL_MV
            if ( PointMap->Drc > 0 )
            {
               fprintf(fileout," Q #%d",(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout," Qs #%d",(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout," C #%d",(int)PointMap->Drc);
            }
          fprintf(fileout,"\n");
          FOR_ROW_COL_MV
            if ( PointMap->Drc > 0 )
            {
               fprintf(fileout," m3/s #%d",(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout," kg/s #%d",(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout," g/l #%d",(int)PointMap->Drc);
            }
           fprintf(fileout,"\n");
           fprintf(fileout,"*%d\n",SOBEKlines);
        }
        else // flat CSV format, comma delimited
        {
           fprintf(fileout,"#LISEM total flow and sed output file for all reporting points in map\n");
           fprintf(fileout,"Time");
           fprintf(fileout,",P");
           FOR_ROW_COL_MV
             if ( PointMap->Drc > 0 )
             {
                 fprintf(fileout,",Q #%d",(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,",Qs #%d",(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,",C #%d",(int)PointMap->Drc);
             }
           fprintf(fileout,"\n");
           fprintf(fileout,"min");
           fprintf(fileout,",mm");
           FOR_ROW_COL_MV
             if ( PointMap->Drc > 0 )
             {
                 fprintf(fileout,",l/s #%d",(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,",kg/s #%d",(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,",g/l #%d",(int)PointMap->Drc);
             }
           fprintf(fileout,"\n");
         }
         fclose(fileout);
      } // all in one
   }  // opening files and writing header

   //######  write all the values #####//
   //######  write all the values #####//
   //######  write all the values #####//
   if (SwitchSeparateOutput)
   {
     FOR_ROW_COL_MV
     {
      if ( PointMap->Drc > 0 ) // all points in separate files
      {
          newname1 = fi.path() + "/" + fi.baseName() + "_" +
         		    QString::number((int)PointMap->Drc) + "." +  fi.suffix();
          // make filename using point number

//          char *name_s = new char[newname1.length()+1];
//          strcpy(name_s, newname1.toAscii().constData());

//          fileout = fopen(name_s, "w");
          fileout = fopen(newname1.toAscii().constData(), "w");

          //delete[] name_s;

        if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
        {
//        	fprintf(fileout,"%8.3f",time/60);
        	fprintf(fileout,"%8ld",runstep);
        	fprintf(fileout," %8.3f", Rain->Drc);
        	fprintf(fileout," %8.3g", Qn->Drc*1000);
        	if (SwitchErosion) fprintf(fileout," %8.3g", Qs->Drc);
        	if (SwitchErosion) fprintf(fileout," %8.3g", SedVol->Drc);
        	fprintf(fileout,"\n");
        }
        else  //SOBEK format
        if (SwitchSOBEKoutput)   //PCRaster timeplot and flat format
        {
        	fprintf(fileout,"\"%s; %02d:%02d:%02d\"", SOBEKstr, hour, min, sec);
        	fprintf(fileout," %8.3g", Qn->Drc);
        	if (SwitchErosion) fprintf(fileout," %8.3g", Qs->Drc);
        	if (SwitchErosion) fprintf(fileout," %8.3g", SedVol->Drc);
        	fprintf(fileout," <\n");
        }
        else
        {
        	fprintf(fileout,"%8.3f",time/60);
        	fprintf(fileout,",%8.3f", Rain->Drc);
        	fprintf(fileout,",%8.3g", Qn->Drc*1000);
        	if (SwitchErosion) fprintf(fileout,",%8.3g", Qs->Drc);
        	if (SwitchErosion) fprintf(fileout,",%8.3g", SedVol->Drc);
        	fprintf(fileout,"\n");
        }
      }  // if point
     }  //rows cols
   } //switch separate
   else
   {
	   // all points in one file
       newname1 = resultDir + outflowFileName;
       // use simply resultdir + filename
//   	   char *name_s = new char[newname1.length()+1];
//       strcpy(name_s, newname1.toAscii().constData());

//       fileout = fopen(name_s, "w");
       fileout = fopen(newname1.toAscii().constData(), "w");
//       //fileout = fopen((const char*)newname1.toLatin1(), "w");

//       delete[] name_s;


      if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
      {
          fprintf(fileout,"%8ld",runstep);
//          fprintf(fileout,"%8.3f",time);
          FOR_ROW_COL_MV
          {
            if ( PointMap->Drc > 0 )
            {
               fprintf(fileout," %8.3f",Rain->Drc);
               fprintf(fileout," %8.3g",Qn->Drc*1000);
               if (SwitchErosion) fprintf(fileout," %8.3g", Qs->Drc);
               if (SwitchErosion) fprintf(fileout," %8.3g", SedVol->Drc);
            }
          }
          fprintf(fileout,"\n");
      }
      else
      if (SwitchSOBEKoutput)
      {
          fprintf(fileout,"\"%s; %02d:%02d:%02d\"", SOBEKstr, hour, min, sec);
          FOR_ROW_COL_MV
          {
            if ( PointMap->Drc > 0 )
            {
               fprintf(fileout," %8.3g",Qn->Drc);
               if (SwitchErosion) fprintf(fileout," %8.3g",Qs->Drc);
               if (SwitchErosion) fprintf(fileout," %8.3g",SedVol->Drc);
            }
          }
          fprintf(fileout," < \n");
      }
      else // flat comma delimited
      {
          fprintf(fileout,"%8.3f",time/60);
          FOR_ROW_COL_MV
          {
            if ( PointMap->Drc > 0 )
            {
               fprintf(fileout,",%8.3f",Rain->Drc);
               fprintf(fileout,",%8.3g",Qn->Drc*1000);
               if (SwitchErosion) fprintf(fileout,",%8.3g", Qs->Drc);
               if (SwitchErosion) fprintf(fileout,",%8.3g", SedVol->Drc);
            }
          }
          fprintf(fileout,"\n");
      }
   }

   fclose(fileout);
}
//---------------------------------------------------------------------------
void TWorld::ReportTotals()
{
	QString newname1 = resultDir + resultFileName;
    // use simply resultdir + filename
	char *name_s = new char[newname1.length()+1];
    strcpy(name_s, newname1.toAscii().constData());

    FILE *fp = fopen(name_s,"w");

    delete[] name_s;

  //  if (fp == NULL)
    //   qFatal(QString("error while opening file: " + resultFileName).toLatin1());

    fprintf(fp,"LISEM run with:,%s\n",(const char *)op.runfilename.toAscii().constData());
    fprintf(fp,"LISEM results at time (min):,%.3f\n",op.time);
    fprintf(fp,"---------------------------------------------\n");
    fprintf(fp,"Catchment area              (ha):,  %12.5f\n",op.CatchmentArea/10000.0);
    fprintf(fp,"Total rainfall              (mm):,  %12.5f\n",op.RainTotmm);
    fprintf(fp,"Total discharge             (mm):,  %12.5f\n",op.Qtotmm);
    fprintf(fp,"Total interception          (mm):,  %12.5f\n",op.IntercTotmm);
    fprintf(fp,"Total infiltration          (mm):,  %12.5f\n",op.InfilTotmm);
    fprintf(fp,"Average surface storage     (mm):,  %12.5f\n",op.SurfStorTotmm);
    fprintf(fp,"Water in runoff             (mm):,  %12.5f\n",op.WaterVolTotmm);
    fprintf(fp,"Total discharge             (m3):,  %12.5f\n",op.Qtot);
    fprintf(fp,"Peak discharge             (l/s):,  %12.5f\n",op.Qpeak);
    fprintf(fp,"Peak time rainfall         (min):,  %12.5f\n",op.RainpeakTime);
    fprintf(fp,"Peak time discharge        (min):,  %12.5f\n",op.QpeakTime);
    fprintf(fp,"Discharge/Rainfall           (%%):,  %12.5f\n",op.RunoffFraction);
    fprintf(fp,"Splash detachment (land)   (ton):,  %12.5f\n",op.DetTotSplash);
    fprintf(fp,"Flow detachment (land)     (ton):,  %12.5f\n",op.DetTotFlow);
    fprintf(fp,"Deposition (land)          (ton):,  %12.5f\n",op.DepTot);
    fprintf(fp,"Suspended Sediment (land)  (ton):,  %12.5f\n",op.SedVolTot);
    fprintf(fp,"Flow detachment (channels) (ton):,  %12.5f\n",op.ChannelDetTot);
    fprintf(fp,"Deposition (channels)      (ton):,  %12.5f\n",op.ChannelDepTot);
    fprintf(fp,"Susp. Sediment (channels)  (ton):,  %12.5f\n",op.ChannelSedTot);
    fprintf(fp,"Total soil loss            (ton):,  %12.5f\n",op.SoilLossTot);
    fprintf(fp,"Average soil loss        (kg/ha):,  %12.5f\n",op.SoilLossTot/op.CatchmentArea*10000.0*1000);
    fclose(fp);

  if (SwitchMulticlass || SwitchNutrients)
  {
    fprintf(fp,"Soil loss mu class 0      (kg):,  %12.5f\n",TotalSedDischargeMu0);
    fprintf(fp,"Soil loss mu class 1      (kg):,  %12.5f\n",TotalSedDischargeMu1);
    fprintf(fp,"Soil loss mu class 2      (kg):,  %12.5f\n",TotalSedDischargeMu2);
    fprintf(fp,"Soil loss mu class 3      (kg):,  %12.5f\n",TotalSedDischargeMu3);
    fprintf(fp,"Soil loss mu class 4      (kg):,  %12.5f\n",TotalSedDischargeMu4);
    fprintf(fp,"Soil loss mu class 5      (kg):,  %12.5f\n",TotalSedDischargeMu5);
  }
  if (SwitchNutrients)
  {
    fprintf(fp,"P loss solution           (kg):,  %12.5f\n",TotalPSolution);
    fprintf(fp,"NO3 loss solution         (kg):,  %12.5f\n",TotalNO3Solution);
    fprintf(fp,"NH4 loss solution         (kg):,  %12.5f\n",TotalNH4Solution);
    fprintf(fp,"P loss suspension         (kg):,  %12.5f\n",TotalPSuspension);
    fprintf(fp,"NO3 loss suspension       (kg):,  %12.5f\n",TotalNO3Suspension);
    fprintf(fp,"NH4 loss suspension       (kg):,  %12.5f\n",TotalNH4Suspension);
    fprintf(fp,"P loss infiltration       (kg):,  %12.5f\n",TotalPInfiltration);
    fprintf(fp,"NO3 loss infiltration     (kg):,  %12.5f\n",TotalNO3Infiltration);
    fprintf(fp,"NH4 loss infiltration     (kg):,  %12.5f\n",TotalNH4Infiltration);
    fprintf(fp,"P loss Deposition         (kg):,  %12.5f\n",TotalPDeposition);
    fprintf(fp,"NO3 loss Deposition       (kg):,  %12.5f\n",TotalNO3Deposition);
    fprintf(fp,"NH4 loss Deposition       (kg):,  %12.5f\n",TotalNH4Deposition);
    fprintf(fp,"P loss Detachment         (kg):,  %12.5f\n",TotalPDetachment);
    fprintf(fp,"NO3 loss Detachment       (kg):,  %12.5f\n",TotalNO3Detachment);
    fprintf(fp,"NH4 loss Detachment       (kg):,  %12.5f\n",TotalNH4Detachment);
  }
    //fprintf(fp,"---------------------------------------------\n");

}

 */

