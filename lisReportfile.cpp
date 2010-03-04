/*
 * reportfile.cpp
 *
 *  Created on: Mar 2, 2010
 *      Author: jetten
 */

#include "model.h"
#include <qstring.h>


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
   char *ext, base[128], newname[128], SOBEKstr[11];
   char sep[4], after[4];
   int hour, min, sec;
   int SOBEKlines = (int) (EndTime-BeginTime)/_dt+1;

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

   base[0]='\0';
   newname[0]='\0';

   strcpy(base, outflowFileName.toLatin1());
   ext = strrchr(outflowFileName.toLatin1(),'.');
   base[outflowFileName.length()- strlen(ext)]= '\0';

   //SOBEK, PCRaster and flat format are mutually exclusive
   if (time <= BeginTime) //  make file at first timestep
   {
      if (SwitchSeparateOutput) // each point separate file
      {
        FOR_ROW_COL_MV
         if ( PointMap->Drc > 0 )
         {
             sprintf(newname,"%s_%d%s",base, (int)PointMap->Drc, ext);
             // make filename using point number
             fileout = fopen(newname, "w");

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
        sprintf(newname,"%s_%d%s",base, (int)PointMap->Drc, ext);
        fileout = fopen(newname, "a");

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
