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
// - TO BE DONE: if runs are interrupted the number of lines win the SOBEK output wioll not be correct!
void TWorld::ReportTimeseries()
{
   FILE *fileout;
   int nr = 0;
   //char *ext, base[128], newname[128];
   char SOBEKstr[11];
   char sep[4], after[4];
   int hour, min, sec;
   int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;
   QString newname1;

   if (SwitchSOBEKOutput && time > 0)
   {
      sec = (int)(time*60);
      hour = (int)(sec/3600);
      min = (int)(time - hour*60);
      sec = (int)(sec - hour * 3600 - min * 60);
   }
   if (SwitchSOBEKOutput)
       strcpy(SOBEKstr, SOBEKdatestring.toLatin1());

   if (SwitchWritePCRtimeplot) strcpy(sep, ",\0");
     else sep[0] = '\0';
   if (SwitchSOBEKoutput) strcpy(sep, " \0");

   if (SwitchWritePCRtimeplot) strcpy(after,"\n\0");
     else after[0] = '\0';

   //base[0]='\0';
   //newname[0]='\0';

   QFileInfo fi(outflowFileName);

   //strcpy(base, outflowFileName.toLatin1());
   //ext = strrchr(outflowFileName.toLatin1(),'.');
   //base[outflowFileName.length()- strlen(ext)]= '\0';

   //SOBEK, PCRaster and flat format are mutually exclusive
   if (time <= BeginTime) //  make file at first timestep
   {
      if (SwitchSeparateOutput) // each point separate file
      {
        FOR_ROW_COL_MV
         if ( PointMap->Drc > 0 )
         {

           //  sprintf(newname,"%s_%d.%s",fi.baseName().toLatin1(), (int)PointMap->Drc, fi.suffix().toLatin1());
             newname1 = fi.path() + "/" + fi.baseName() + "_" + QString::number((int)PointMap->Drc) + "." +  fi.suffix();
             // make filename using point number
             fileout = fopen(newname1.toLatin1(), "w");

             // HEADERS for the 3 types
             if (SwitchWritePCRtimeplot)  //PCRaster timeplot format, cannot be SOBEK !
             {
                 fprintf(fileout,"#LISEM total flow and sed output file for point %d\n", (int)PointMap->Drc);
                 if (SwitchErosion)
                   fprintf(fileout,"5\n");
                 else
                   fprintf(fileout,"3\n");
                 fprintf(fileout,"Time (min)");
                 fprintf(fileout,"%s%sPavg (mm)",after,sep);
                 fprintf(fileout,"%s%sQ (l/s)",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sQs (kg/s)",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sC conc (g/l)",after,sep);
                 fprintf(fileout,"\n");
             }
             else // SOBEK format
             if (SwitchSOBEKoutput)
             {
                 fprintf(fileout,"%s%sQ",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sQs",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sC",after,sep);
                 fprintf(fileout,"\n");
                 fprintf(fileout,"%s%sm3/s",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%skg/s",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sg/l",after,sep);
                 fprintf(fileout,"\n");
                 fprintf(fileout,"*%d\n",SOBEKlines);
             }
             else // flat format, comma delimited
             {
                 fprintf(fileout,"#LISEM total flow and sed output file for point %d\n", (int)PointMap->Drc);
                 fprintf(fileout,"Time");
                 fprintf(fileout,"%s%sPavg",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sQs",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sC",after,sep);
                 fprintf(fileout,"\n");
                 fprintf(fileout,"min");
                 fprintf(fileout,"%s%smm",after,sep);
                 fprintf(fileout,"%s%sl/s",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%skg/s",after,sep);
                 if (SwitchErosion) fprintf(fileout,"%s%sg/l",after,sep);
                 fprintf(fileout,"\n");
             }
             fclose(fileout);
         }
      }  // separate
      else   // HEADERS: all points in one file
      {
        // count nr of points
        FOR_ROW_COL_MV
          if ( PointMap->Drc > 0 )
            nr++;

        fileout = fopen(outflowFileName.toLatin1(), "w");
        // open file using filename, all in one file
        //fprintf(fileout,"#LISEM total flow and sed output file for all reporting points in map\n");

        if (SwitchWritePCRtimeplot)   //PCRaster timeplot format
        {
            fprintf(fileout,"#LISEM total flow and sed output file for all reporting points in map\n");
            fprintf(fileout,"%d\n",2+(1+(SwitchErosion ? 2 : 0))*nr);
            // number of variables = 2 + (1 or 2) * number of points
            fprintf(fileout,"Time (min)%s",after);
            fprintf(fileout,"%sPavg (mm)",sep);
            FOR_ROW_COL_MV
              if ( PointMap->Drc > 0 )
              {
                 fprintf(fileout,"%s%sQ #%d (l/s)",after,sep,(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"%s%sQs #%d (kg/s)",after,sep,(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"%s%sC #%d (g/l)",after,sep,(int)PointMap->Drc);
              }
            fprintf(fileout,"\n");
        }
        else // SOBEK format
        if (SwitchSOBEKoutput) //note: sobek input does not have rainfall
        {
          FOR_ROW_COL_MV
            if ( PointMap->Drc > 0 )
            {
               fprintf(fileout,"%s%sQ #%d",after,sep,(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout,"%s%sQs #%d",after,sep,(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout,"%s%sC #%d",after,sep,(int)PointMap->Drc);
            }
          fprintf(fileout,"\n");
          FOR_ROW_COL_MV
            if ( PointMap->Drc > 0 )
            {
               fprintf(fileout,"%s%sm3/s #%d",after,sep,(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout,"%s%skg/s #%d",after,sep,(int)PointMap->Drc);
               if (SwitchErosion) fprintf(fileout,"%s%sg/l #%d",after,sep,(int)PointMap->Drc);
            }
           fprintf(fileout,"\n");
           fprintf(fileout,"*%d\n",SOBEKlines);
        }
        else // flat CSV format
        {
           fprintf(fileout,"#LISEM total flow and sed output file for all reporting points in map\n");
           fprintf(fileout,"Time%s",after);
           fprintf(fileout,"%sPavg",sep);
           FOR_ROW_COL_MV
             if ( PointMap->Drc > 0 )
             {
                 fprintf(fileout,"%s%sQ #%d",after,sep,(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"%s%sQs #%d",after,sep,(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"%s%sC #%d",after,sep,(int)PointMap->Drc);
             }
           fprintf(fileout,"\n");
           fprintf(fileout,"min%s",after);
           fprintf(fileout,"%smm",sep);
           FOR_ROW_COL_MV
             if ( PointMap->Drc > 0 )
             {
                 fprintf(fileout,"%s%sl/s #%d",after,sep,(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"%s%skg/s #%d",after,sep,(int)PointMap->Drc);
                 if (SwitchErosion) fprintf(fileout,"%s%sg/l #%d",after,sep,(int)PointMap->Drc);
             }
           fprintf(fileout,"\n");
         }
         fclose(fileout);
      } // all in one
   }  // opening files and writing header

   // write all the values
   if (SwitchSeparateOutput)
   {
     FOR_ROW_COL_MV
      if ( PointMap->Drc > 0 )
      {
        newname1 = fi.path() + "/" + fi.baseName() + "_" + QString::number((int)PointMap->Drc) + "." +  fi.suffix();
       // sprintf(newname,"%s_%d.%s",base, (int)PointMap->Drc, ext);
        fileout = fopen(newname1.toLatin1(), "a");

        if (!SwitchSOBEKoutput)   //PCRaster timeplot and flat format
        {
          fprintf(fileout,"%8.3f",time);
          fprintf(fileout,"%s%8.3f", sep, Rain->Drc);
          fprintf(fileout,"%s%8.3g", sep, Q->Drc);
          if (SwitchErosion) fprintf(fileout,"%s%8.3g", sep, Qs->Drc);
          if (SwitchErosion) fprintf(fileout,"%s%8.3g", sep, SedVol->Drc);
          fprintf(fileout,"\n");
        }
        else  //SOBEK format
        {
          fprintf(fileout,"\"%s; %02d:%02d:%02d\"", SOBEKstr, hour, min, sec);
          fprintf(fileout,"%s%8.3g", sep, Q->Drc);
          if (SwitchErosion) fprintf(fileout,"%s%8.3g", sep, Qs->Drc);
          if (SwitchErosion) fprintf(fileout,"%s%8.3g", sep, SedVol->Drc);
          fprintf(fileout," <\n");
        }
      }
   }
   else
   {
      after[0] = '\0';
      fileout = fopen(outflowFileName.toLatin1(), "a");

      if (!SwitchSOBEKoutput)
      {
        fprintf(fileout,"%8.3f",time);
        FOR_ROW_COL_MV
          if ( PointMap->Drc > 0 )
          {
             fprintf(fileout,"%s%s%8.3f",after, sep, Rain->Drc);
             fprintf(fileout,"%s%s%8.3g", after, sep, Q->Drc);
             if (SwitchErosion) fprintf(fileout,"%s%s%8.3g", after, sep, Qs->Drc);
             if (SwitchErosion) fprintf(fileout,"%s%s%8.3g", after, sep, SedVol->Drc);
          }
        fprintf(fileout,"\n");
      }
      else
      {
        fprintf(fileout,"\"%s; %02d:%02d:%02d\"", SOBEKstr, hour, min, sec);
        FOR_ROW_COL_MV
          if ( PointMap->Drc > 0 )
          {
             fprintf(fileout,"%s%s%8.3g", after, sep, Q->Drc);
             if (SwitchErosion) fprintf(fileout,"%s%s%8.3g", after, sep, Qs->Drc);
             if (SwitchErosion) fprintf(fileout,"%s%s%8.3g", after, sep, SedVol->Drc);
          }
        fprintf(fileout," < \n");
      }
   }

   fclose(fileout);
}
//---------------------------------------------------------------------------
void TWorld::ReportTotals()
{
    FILE *fp = fopen(resultFileName.toLatin1(),"w");
 //   if (fp == NULL)
   //    DEBUG(QString("error while opening file: ") + resultFileName);

    fprintf(fp,"LISEM run with: %s\n",(const char *)op.runfilename.toLatin1());
    fprintf(fp,"LISEM results at time (min): %.3f\n",op.time);
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
    fprintf(fp,"Average soil loss        (kg/ha):,  %12.5f\n",op.SoilLossTot/CatchmentArea*10000.0);
/*
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
 */
    //fprintf(fp,"---------------------------------------------\n");
    fclose(fp);

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

    	TotalDetMap->report(totalErosionFileName);
    	TotalDepMap->report(totalDepositionFileName);
    	TotalSoillossMap->report(totalSoillossFileName);
    }

    if (outputcheck[0].toInt() == 1) Qn->report(Outrunoff);


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


