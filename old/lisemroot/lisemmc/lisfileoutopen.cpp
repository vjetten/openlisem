
// *****************************************************************************
// ******** OPEN and print first line output files *****************************
// *****************************************************************************

//       reportPointMap(outPointFileName, OutPoint, OutPoint, -1);
//       reportMap(outflowFileName, Outlet, -1, NULL_MAP, NULL_MAP, NULL_MAP, NULL_MAP, SwitchSeparateOutput);
//       reportPointMap(outflowFileName, OutPoint, -1, OutPoint, OutPoint, OutPoint, SwitchSeparateOutput);

/*
       if (SwitchOutlet1)
       {
         fpoutflow1 = fopen(outflowFileName,"w");
         if (fpoutflow1 == NULL)
           LisemError("error while opening file output main outlet: ",outflowFileName);
          if (!SwitchMulticlass)
          {
            if (SwitchWritePCRtimeplot)
             fprintf(fpoutflow1,timeplotoutputstr);
            else
             fprintf(fpoutflow1,commaoutputstr);
            if (SwitchWritePCRtimeplot)
             fprintf(fpoutflow1,ncommaformat,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0);
            else
             fprintf(fpoutflow1,commaformat,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0);
          }
          else
          {
            if (SwitchNutrients)
            {
              if (SwitchWritePCRtimeplot)
               fprintf(fpoutflow1,timeplotoutputstrNUT);
              else
               fprintf(fpoutflow1,commaoutputstrNUT);
              if (SwitchWritePCRtimeplot)
                fprintf(fpoutflow1,ncommaformatNUT,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
              else
                fprintf(fpoutflow1,commaformatNUT,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            }
            else
            {
              if (SwitchWritePCRtimeplot)
                fprintf(fpoutflow1,timeplotoutputstrMC);
              else
                fprintf(fpoutflow1,commaoutputstrMC);
              if (SwitchWritePCRtimeplot)
                fprintf(fpoutflow1,ncommaformatMC,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
              else
                fprintf(fpoutflow1,commaformatMC,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            }
          }
         fflush(fpoutflow1);
         fclose(fpoutflow1);
       }
       if (SwitchOutlet2)
       {
         fpoutflow2 = fopen(outflowFileName2,"w");
         if (fpoutflow2 == NULL)
            LisemError("error while opening file for ouput subcatchment 1: ",outflowFileName2);
         if (SwitchWritePCRtimeplot)
           fprintf(fpoutflow2,timeplotoutputstr);
         else
           fprintf(fpoutflow2,commaoutputstr);
         if (SwitchWritePCRtimeplot)
           fprintf(fpoutflow2,ncommaformat,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0);
         else
           fprintf(fpoutflow2,commaformat,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0);
         fclose(fpoutflow2);
       }
       if (SwitchOutlet3)
       {
         fpoutflow3 = fopen(outflowFileName3,"w");
         if (fpoutflow3 == NULL)
            LisemError("error while opening file for ouput subcatchment 2: ",outflowFileName3);
         if (SwitchWritePCRtimeplot)
           fprintf(fpoutflow3,timeplotoutputstr);
         else
           fprintf(fpoutflow3,commaoutputstr);
         if (SwitchWritePCRtimeplot)
           fprintf(fpoutflow3,ncommaformat,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0);
         else
           fprintf(fpoutflow3,commaformat,STARTINTERVAL, 0.0, 0.0, 0.0, 0.0);
         fclose(fpoutflow3);
       }
 */
// *****************************************************************************
// ******** OPEN and print first line output files *****************************
// *****************************************************************************

//VJ 040823 add buffer output
       if (BufferFileName[0] != '\0')
       {
           BufferFout = fopen(BufferFileName,"w");
           if (BufferFout == NULL)
             LisemError("error while opening file output buffer(s)",BufferFileName);

           if (SwitchWritePCRtimeplot)
               fprintf(BufferFout,timeplotoutputstrBuffer);
           else
               fprintf(BufferFout,commaoutputstrBuffer);

           if (SwitchWritePCRtimeplot)
               fprintf(BufferFout,ncommaformat4,STARTINTERVAL, 0.0, 0.0, 0.0);
           else
               fprintf(BufferFout,commaformat4,STARTINTERVAL, 0.0, 0.0, 0.0);
       }



//VJ 050913 add pest output option: file pestout.txt in result dir
      if (SwitchPestout)
      {
          strcpy(PestoutFileName, CatPath("pestout.txt", RESPATH));
          fpestout = fopen(PestoutFileName,"w");
          fprintf(fpestout,"# pest result for lisem, discharge main outlet\n");
          fprintf(fpestout,"2\n");
          fprintf(fpestout,"time (min)\n");
          fprintf(fpestout,"Q (l/s)\n");
          fprintf(fpestout," %.3f %.3f\n",0,0);
          fclose(fpestout);
       }



