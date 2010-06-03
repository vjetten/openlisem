// *****************************************************************************
// ******** OUTPUT maps ********************************************************
// *****************************************************************************
       _spatial(REAL4, totwater);
       _spatial(REAL4, totsediment);
       calc(" totwater = WaterHVolout ");
       calc(" totsediment = Sedin ");
       //when MC or Nutrients Sedin is summed in LisPMCkinematic
       //in all functions concentration is calculated as sed/volume, NOT Qsed/Q
       //this assumes instant mixing of all sediment in total volume in the timestep

       if (SwitchIncludeChannel)
       {
         calc(" ChannelSedin = cover(ChannelSedin, 0)");
         calc(" ChannelVolout = cover(ChannelVolout, 0)");
         calc(" totwater += ChannelVolout ");
         calc(" totsediment += ChannelSedin ");
       }
       if (SwitchWheelPresent)
       {
         calc(" WheelSedin = cover(WheelSedin, 0)");
         calc(" WheelVolout = cover(WheelVolout, 0)");
         calc(" totwater += WheelVolout ");
         calc(" totsediment += WheelSedin ");
       }
       _spatial(REAL4, sedconcfile);
       calc("sedconcfile = mif(totwater gt 1e-6, min(850,totsediment/totwater), 0)");
//VJ 090214 limit to 1e-6 (=0.001 liter)
//VJ 090214 limit to max 850
       //in all functions concentration is calculated as sed/volume, NOT Qsed/Q

       if ((PERIODMAPS && PERIODCOUNT == PERIOD) ||
           ((WRITETIME[timeIndex]<= timestepindex) && !(WRITETIME[timeIndex] == 0))
          )
        {
           LastPCRTimestep++;
             // increase sequential timestepnr for pcraster output
          //VJ
          // toggle between output in discharge l/s or unit discharge l/s/m
           _spatial(REAL4, QOVERLANDPLUSCHANNEL);
           if(SwitchRunoffPerM)
                calc(" QOVERLANDPLUSCHANNEL = Qout/DX*1000 ");
           else
                calc(" QOVERLANDPLUSCHANNEL = Qout*1000 ");

           if (SwitchMapoutRunoff) writeTimeseries(QOVERLANDPLUSCHANNEL, mapname("outrunoff"));
           if (SwitchMapoutWH) writeTimeseries(WH, mapname("outwh"));
           if (SwitchMapoutRWH) writeTimeseries(RunoffMeanHin, mapname("outrwh"));

           if (SwitchMapoutTC) writeTimeseries(TransportCapOUT, mapname("outtc"));

           if (SwitchMapoutV) writeTimeseries(V, mapname("outvelo"));
           _spatial(REAL4, InfilOutput);
           calc(" InfilOutput = InfilVol/(DX*DXc)*1000 ");
           if (SwitchMapoutInf) writeTimeseries(InfilOutput, mapname("outinf"));
           if (SwitchMapoutSs) writeTimeseries(SurfStorH, mapname("outss"));

           if (SwitchChannelBaseflow)
             if (SwitchMapoutChvol) writeTimeseries(ChannelVolout, mapname("outchvol"));
             //??? not sure why this?

           if (!SwitchNoErosion)
           {
               if (SwitchMapoutConc)
                  writeTimeseries(sedconcfile, mapname("outconc"));
                  //in all functions concentration is calculated as sed/volume, NOT Qsed/Q

           	//VJ 050822 three units to choose from
               _spatial(REAL4, unitfactor);
               if (ErosionUnits == 0)
	               calc(" unitfactor = 1.0 ");
               if (ErosionUnits == 1)
                 	calc(" unitfactor = 1.0/(DXc*DX) ");
               if (ErosionUnits == 2)
                 	calc(" unitfactor = 10.0/(DXc*DX) ");

               _spatial(REAL4, Erosion);
               calc(" Erosion = SumTotErosion * unitfactor");

               write(Erosion, totalErosionFileName);

               if (SwitchMapoutEros)
                    writeTimeseries(Erosion, mapname("outeros"));

               _spatial(REAL4, Deposition);
               calc(" Deposition = SumTotDeposition * unitfactor ");

               write(Deposition, totalDepositionFileName);

               if (SwitchMapoutDepo)
                    writeTimeseries(Deposition, mapname("outdepo"));

               if (SwitchMulticlass || SwitchNutrients)
                  mcalc(" SedConc# = mif(WaterHVolin gt 0 ,SedinMu#/WaterHVolin,0) ",6);

               if (SwitchMulticlass)
               {
                 if (SwitchMapoutMC0) writeTimeseries(SedConc0, mapname("outmu0"));
                 if (SwitchMapoutMC1) writeTimeseries(SedConc1, mapname("outmu1"));
                 if (SwitchMapoutMC2) writeTimeseries(SedConc2, mapname("outmu2"));
                 if (SwitchMapoutMC3) writeTimeseries(SedConc3, mapname("outmu3"));
                 if (SwitchMapoutMC4) writeTimeseries(SedConc4, mapname("outmu4"));
                 if (SwitchMapoutMC5) writeTimeseries(SedConc5, mapname("outmu5"));
                 if (SwitchMapoutMC6) writeTimeseries(D50susp, mapname("outD50susp"));
               }
               if (SwitchNutrients)
               {
                 if (SwitchMapoutPsol) writeTimeseries(NutPSolution, mapname("outpsolut"));
                 if (SwitchMapoutPsus) writeTimeseries(NutPSuspension, mapname("outpsus"));
                 if (SwitchMapoutPinf) writeTimeseries(NutPInfiltration, mapname("outpinf"));
                 if (SwitchMapoutNH4sol) writeTimeseries(NutNH4Solution, mapname("outnh4solut"));
                 if (SwitchMapoutNH4sus) writeTimeseries(NutNH4Suspension, mapname("outnh4sus"));
                 if (SwitchMapoutNH4inf) writeTimeseries(NutNH4Infiltration, mapname("outnh4inf"));
                 if (SwitchMapoutNO3sol) writeTimeseries(NutNO3Solution, mapname("outno3solut"));
                 if (SwitchMapoutNO3sus) writeTimeseries(NutNO3Suspension, mapname("outno3sus"));
                 if (SwitchMapoutNO3inf) writeTimeseries(NutNO3Infiltration, mapname("outno3inf"));
						//timeseries
						
                 if (SwitchMapoutPdep)   write(SumTotNutPdep, mapname("outpdep"));
                 if (SwitchMapoutNH4dep) write(SumTotNutNH4dep, mapname("outnh4dep"));
                 if (SwitchMapoutNO3dep) write(SumTotNutNO3dep, mapname("outno3dep"));
                 if (SwitchMapoutPdet)   write(SumTotNutPdet, mapname("outpdet"));
                 if (SwitchMapoutNH4det) write(SumTotNutNH4det, mapname("outnh4det"));
                 if (SwitchMapoutNO3det) write(SumTotNutNO3det, mapname("outno3det"));
                 // total for the area in kg per m2 (conversion is done in NUTmassbalsed
               }
               if (SwitchGullies)
               {
                 if (SwitchMapoutGul0) writeTimeseries(GullyDepth, mapname("outguld"));
                 if (SwitchMapoutGul1) writeTimeseries(GullyWidthDX, mapname("outgulw"));
                 if (SwitchMapoutGul2) writeTimeseries(GullyArea, mapname("outgula"));
                 if (SwitchMapoutGul3) writeTimeseries(GullyCSFraction, mapname("outgulf"));
                 if (SwitchMapoutGul4) writeTimeseries(DEM, mapname("outguldem"));
               }
           } // if !noerosion

           timeIndex++;

       } //if periodmaps

//VJ 031218 set timeseries also when not PCRaster output
       if(timestepindex >= ENDINTERVAL)
       {
           LisIFace->Messages->Lines->Append("Setting timeseries min and max, please wait ...");

           if (SwitchMapoutRunoff) 	SetTimeseriesMinmax(mapname("outrunoff"));
           if (SwitchMapoutWH) 		SetTimeseriesMinmax(mapname("outwh"));
           if (SwitchMapoutRWH) 	SetTimeseriesMinmax(mapname("outrwh"));
           if (SwitchMapoutTC) 		SetTimeseriesMinmax(mapname("outtc"));
           if (SwitchMapoutV) 		SetTimeseriesMinmax(mapname("outvelo"));
           if (SwitchMapoutInf) 	SetTimeseriesMinmax(mapname("outinf"));
           if (SwitchMapoutSs) 		SetTimeseriesMinmax(mapname("outss"));
//?????????
           if (SwitchChannelBaseflow)
            if (SwitchMapoutChvol) SetTimeseriesMinmax(mapname("outchvol"));

           if (!SwitchNoErosion)
           {
               if (SwitchMapoutConc)    SetTimeseriesMinmax(mapname("outconc"));
               //in all functions concentration is calculated as sed/volume, NOT Qsed/Q

               if (SwitchMapoutEros)    SetTimeseriesMinmax(mapname("outeros"));
               if (SwitchMapoutDepo)    SetTimeseriesMinmax(mapname("outdepo"));

               if (SwitchMulticlass)
               {
                 if (SwitchMapoutMC0) SetTimeseriesMinmax(mapname("outmu0"));
                 if (SwitchMapoutMC1) SetTimeseriesMinmax(mapname("outmu1"));
                 if (SwitchMapoutMC2) SetTimeseriesMinmax(mapname("outmu2"));
                 if (SwitchMapoutMC3) SetTimeseriesMinmax(mapname("outmu3"));
                 if (SwitchMapoutMC4) SetTimeseriesMinmax(mapname("outmu4"));
                 if (SwitchMapoutMC5) SetTimeseriesMinmax(mapname("outmu5"));
                 if (SwitchMapoutMC6) SetTimeseriesMinmax(mapname("outD50susp"));
               }
               if (SwitchNutrients)
               {
                 if (SwitchMapoutPsol)   SetTimeseriesMinmax(mapname("outpsolut"));
                 if (SwitchMapoutPsus)   SetTimeseriesMinmax(mapname("outpsus"));
                 if (SwitchMapoutPinf)   SetTimeseriesMinmax(mapname("outpinf"));
                 if (SwitchMapoutNH4sol) SetTimeseriesMinmax(mapname("outnh4solut"));
                 if (SwitchMapoutNH4sus) SetTimeseriesMinmax(mapname("outnh4sus"));
                 if (SwitchMapoutNH4inf) SetTimeseriesMinmax(mapname("outnh4inf"));
                 if (SwitchMapoutNO3sol) SetTimeseriesMinmax(mapname("outno3solut"));
                 if (SwitchMapoutNO3sus) SetTimeseriesMinmax(mapname("outno3sus"));
                 if (SwitchMapoutNO3inf) SetTimeseriesMinmax(mapname("outno3inf"));
               }
               if (SwitchGullies)
               {
                 if (SwitchMapoutGul0) SetTimeseriesMinmax(mapname("outguld"));
                 if (SwitchMapoutGul1) SetTimeseriesMinmax(mapname("outgulw"));
                 if (SwitchMapoutGul2) SetTimeseriesMinmax(mapname("outgula"));
                 if (SwitchMapoutGul3) SetTimeseriesMinmax(mapname("outgulf"));
                 if (SwitchMapoutGul4) SetTimeseriesMinmax(mapname("outguldem"));
               }
           } // if !noerosion
           LisIFace->Messages->Lines->Append("Done!");
       }//laststep
       // ******************************************************************

       PERIODCOUNT++;
       if (PERIODCOUNT > PERIOD)
          PERIODCOUNT = 1;



// *****************************************************************************
// ******** OUTPUT to Outlet files *********************************************
// *****************************************************************************



//NOTE: ??????????? QOUT AND QSEDOUT ARE ACCUMULATED WITH CHANNEL AND WHEEL AND GULLY IN MASSBAL AND MASSBALSED
//NOTE: SED CONC TO FILE IS VOL/VOL WHILE TO COLUMN FILE IS FLUX/FLUX BECAUSE THAT IS WHAT IS LEAVING

//VJ 040216 if more outlet 1 then recalc conc as sum sediment /sum water
/*
       _nonspatial(REAL4, SCOutlet1);
       calc(" SCOutlet1 = sum(mif(Outlet1 eq 1, Qout, 0))");   //qout in m3/s
       if (SCOutlet1 > 1e-6)
         calc(" SCOutlet1 = sum(mif(Outlet1 eq 1, Qsedout, 0))/SCOutlet1 "); //qsedout in kg/s
       else
          SCOutlet1 = 0;   //SCOutlet in kg/s / m3/s = kg/m3 = g/l
       SCOutlet1 = MIN(850, SCOutlet1);
*/

//VJ 080619  SED CONCENTRATION AS weight SED / VOLUME water at outlet
//VJ 090214 simplified to sedconcfile
//OULET 1

       char *pstr, *pstrNUT, *pstrMC;
       int oR;

       _nonspatial(REAL4, SCOutlet1);
       calc(" SCOutlet1 = sum(mif(Outlet1 eq 1, sedconcfile, 0)) ");
       _nonspatial(REAL4, QOutlet1);
       calc(" QOutlet1 = sum(mif(Outlet1 eq 1, Qout, 0)) ");
       _nonspatial(REAL4, QSEDOutlet1);
       QSEDOutlet1 = QOutlet1*SCOutlet1 ;
/*
//OULET 2
       _nonspatial(REAL4, SCOutlet2);
       calc(" SCOutlet2 = sum(mif(Outlet2 eq 1, sedconcfile, 0)) ");
       _nonspatial(REAL4, QOutlet2);
       calc(" QOutlet2 = sum(mif(Outlet2 eq 1, Qout, 0)) ");
       _nonspatial(REAL4, QSEDOutlet2);
       QSEDOutlet2 = QOutlet2*SCOutlet2 ;

//OULET 3
       _nonspatial(REAL4, SCOutlet3);
       calc(" SCOutlet3 = sum(mif(Outlet2 eq 1, sedconcfile, 0)) ");
       _nonspatial(REAL4, QOutlet3);
       calc(" QOutlet3 = sum(mif(Outlet3 eq 1, Qout, 0)) ");
       _nonspatial(REAL4, QSEDOutlet3);
       QSEDOutlet3 = QOutlet3*SCOutlet3 ;
*/

//=========------- output to files-------===========//

       if (SwitchWritePCRtimeplot)
       {
         pstr = ncommaformat;
         pstrNUT = ncommaformatNUT;
         pstrMC = ncommaformatMC;
       }
       else
       {
         pstr = commaformat;
         pstrNUT = commaformatNUT;
         pstrMC = commaformatMC;
       }
       //defined in lisheadin.h

       _spatial(REAL4, Qoutr);
       calc(" Qoutr = Qout*sedconcfile");

       _spatial(REAL4, QSout);
       calc(" QSout = Qout*sedconcfile");

       _spatial(REAL4, RainAvgOut);
       calc(" RainAvgOut = RainfallAverageH*3600.0/DTSEC ");

// VJ 091211 NEW REPORTING
 //      reportPointMap(outflowFileName, Outlet, timestepindex, RainAvgOut, Qoutr, QSout, sedconcfile, SwitchSeparateOutput);

/*
       if (SwitchOutlet1)
       {
           fpoutflow1 = fopen(outflowFileName,"a");
           if (fpoutflow1 == NULL)
              LisemError("error while opening file: ",outflowFileName);

           if (!SwitchMulticlass)
           {
             oR = fprintf(fpoutflow1,pstr,
                     timestepindex,
                     RainfallAverageH*3600.0/DTSEC,
//VJ 080620 changed
                     QOutlet1*1000,
                     QSEDOutlet1,
                     SCOutlet1 );
           }
           else  //MULTICLASS
           {
             _nonspatial(REAL4, SCOutlet1Mu0);
             _nonspatial(REAL4, SCOutlet1Mu1);
             _nonspatial(REAL4, SCOutlet1Mu2);
             _nonspatial(REAL4, SCOutlet1Mu3);
             _nonspatial(REAL4, SCOutlet1Mu4);
             _nonspatial(REAL4, SCOutlet1Mu5);

             if (QOutflow > 1e-6)
             {
                 SCOutlet1Mu0 = SedOutflow0/QOutflow;
                 SCOutlet1Mu1 = SedOutflow1/QOutflow;
                 SCOutlet1Mu2 = SedOutflow2/QOutflow;
                 SCOutlet1Mu3 = SedOutflow3/QOutflow;
                 SCOutlet1Mu4 = SedOutflow4/QOutflow;
              	  SCOutlet1Mu5 = SedOutflow5/QOutflow;
             }
             else
             {
                 SCOutlet1Mu0 = 0;
                 SCOutlet1Mu1 = 0;
                 SCOutlet1Mu2 = 0;
                 SCOutlet1Mu3 = 0;
                 SCOutlet1Mu4 = 0;
             	  SCOutlet1Mu5 = 0;
             }
             if (!SwitchNutrients)
             {
                 oR = fprintf(fpoutflow1,pstrMC,
                         timestepindex,
                         RainfallAverageH*3600.0/DTSEC,
//VJ 080620 changed
                         QOutlet1*1000,
                         QSEDOutlet1,
                         SCOutlet1,
                         SCOutlet1Mu0,SCOutlet1Mu1,SCOutlet1Mu2,
                         SCOutlet1Mu3,SCOutlet1Mu4,SCOutlet1Mu5);
             }
             else
             {
//VJ 030415 changed output of nutrients, conc as flux nut/flux water
                 _nonspatial(REAL4, NutPConc1);
                 _nonspatial(REAL4, NutNH4Conc1);
                 _nonspatial(REAL4, NutNO3Conc1);
                 // solution
                 _nonspatial(REAL4, NutPConc2);
                 _nonspatial(REAL4, NutNH4Conc2);
                 _nonspatial(REAL4, NutNO3Conc2);
                 // suspension
                 
                 if (QOutflow > 1e-8)
                 {
                     NutPConc1 = NutPOutflow/QOutflow;
                     NutNH4Conc1 = NutNH4Outflow/QOutflow;
                     NutNO3Conc1 = NutNO3Outflow/QOutflow;
                     //solution concentration kg/s / m3/s = g/l
                     NutPConc2 = NutPOutflows/QOutflow;
                     NutNH4Conc2 = NutNH4Outflows/QOutflow;
                     NutNO3Conc2 = NutNO3Outflows/QOutflow;
                     //suspension concentration kg/s / m3/s = g/l
                 }
                 else
                 {
                     NutPConc1   = 0;
                     NutNH4Conc1 = 0;
                     NutNO3Conc1 = 0;                     
                     NutPConc2   = 0;
                     NutNH4Conc2 = 0;
                     NutNO3Conc2 = 0;
                 }
                 oR = fprintf(fpoutflow1,pstrNUT,
                         timestepindex,
                         RainfallAverageH*3600.0/DTSEC,
//VJ 080620 changed
                         QOutlet1*1000,
                         QSEDOutlet1,
                         SCOutlet1,
                         SCOutlet1Mu0,SCOutlet1Mu1,SCOutlet1Mu2,
                         SCOutlet1Mu3,SCOutlet1Mu4,SCOutlet1Mu5,
                         NutPConc1, NutNH4Conc1, NutNO3Conc1,
                         NutPConc2, NutNH4Conc2, NutNO3Conc2);
             }
           }
           fclose(fpoutflow1);
       }


       // ***** outlet subcatchment 2 ****/

 /*
       _nonspatial(REAL4, SCOutlet2);
       calc(" SCOutlet2 = sum(mif(Outlet2 eq 1, Qout, 0))");
       if (SCOutlet2 > 1e-6)
         calc(" SCOutlet2 = sum(mif(Outlet2 eq 1, Qsedout, 0))/SCOutlet2 ");
       else
          SCOutlet2 = 0;
       SCOutlet2 = MIN(850, SCOutlet2);

//VJ 040216 if more outlet 2 then recalc conc as sum sediment /sum water
       _nonspatial(REAL4, SCOutlet2);
//       calc(" SCOutlet2 = sum(mif(Outlet2 eq 1, sedconcfile, 0)) ");
       calc(" SCOutlet2 = sum(mif(Outlet2 eq 1, totwater, 0))");

//VJ 050808 fixed a typing error here
       if (SCOutlet2 > 0)
         calc(" SCOutlet2 = sum(mif(Outlet2 eq 1, totsediment, 0))/SCOutlet2 ");
       else
          SCOutlet2 = 0;
*/
/*
       if (SwitchOutlet2)
       {
           fpoutflow2 = fopen(outflowFileName2,"a");
           if (fpoutflow2 == NULL)
              LisemError("error while opening file output subcatchment 1: ",outflowFileName2);

           oR = fprintf(fpoutflow2,pstr,
                   timestepindex,
                   RainfallAverageH*3600.0/DTSEC,
                   QOutlet2*1000,
                   QSEDOutlet2,
                   SCOutlet2 );
           fclose(fpoutflow2);
           POSTCOND(oR > 5);
       }
*/

       // ***** outlet subcatchment 3 ****/

//VJ 040216 if more outlet 3 then recalc conc as sum sediment /sum water
/*
       _nonspatial(REAL4, QOutlet3);
       calc(" QOutlet3 = sum(mif(Outlet3 eq 1, Qout, 0)) ");
       calc(" SCOutlet3 = sum(mif(Outlet3 eq 1, totwater, 0))");
//VJ 050808 fixed a typing error here
       if (SCOutlet3 > 0)
         calc(" SCOutlet3 = sum(mif(Outlet3 eq 1, totsediment, 0))/SCOutlet3 ");
       else
          SCOutlet3 = 0;
//       calc(" SCOutlet3 = sum(mif(Outlet3 eq 1, sedconcfile, 0)) ");
       _nonspatial(REAL4, SCOutlet3);
       calc(" SCOutlet3 = sum(mif(Outlet3 eq 1, Qout, 0))");
       if (SCOutlet3 > 1e-6)
         calc(" SCOutlet3 = sum(mif(Outlet3 eq 1, Qsedout, 0))/SCOutlet3 ");
       else
          SCOutlet3 = 0;

*/
 /*
       if (SwitchOutlet3)
       {
           fpoutflow3 = fopen(outflowFileName3,"a");
           if (fpoutflow3 == NULL)
              LisemError("error while opening file output subcatchment 2: ",outflowFileName3);

           oR = fprintf(fpoutflow3,pstr,
                   timestepindex,
                   RainfallAverageH*3600.0/DTSEC,
                   QOutlet3*1000,
                   QSEDOutlet3,
                   SCOutlet3 );
           fclose(fpoutflow3);
           POSTCOND(oR > 5);
       }
*/

//VJ 050913 add pestout possibility
      if (SwitchPestout)
      {
		    pestcounter++;
   	    if (pestcounter+0.1 > pestreport)
      	 {
             fpestout = fopen(PestoutFileName,"a");
	          fprintf(fpestout," %.3f %.3f\n",timestepindex,QOutflow*1000);
   	       fclose(fpestout);
      	    pestcounter = 0;
	       }
      }

//VJ 080217 moved bufferoutput to here, experimental, not well tested
      if (SwitchBuffers){
           timeoutput("buffervol.csv",BufferID,BufferVolumeCurrent,timestepindex);
      }

        //********TOTALS FILE UPDATE *************


        FILE *fp = fopen(resultFileName,"w");
        if (fp == NULL)
           LisemError("error while opening file: ",resultFileName);

        fprintf(fp,"LISEM run with: %s\n",LisIFace->RunFilename.c_str());
        fprintf(fp,"LISEM results at time (min): %.3f\n",(timestepindex));
        fprintf(fp,"---------------------------------------------\n");
        fprintf(fp,"Catchment area            (ha):,  %f\n",CatchmentArea/10000.0);
        fprintf(fp,"Total rainfall            (mm):,  %12.5f\n",TotalRainVol/CatchmentArea*1000.0);
        fprintf(fp,"Total discharge           (mm):,  %12.5f\n",TotalDischarge/CatchmentArea*1000.0);
        fprintf(fp," Total interception       (mm):,  %12.5f\n",TotalInterceptionVol/CatchmentArea*1000.0);
        fprintf(fp," Total infiltration       (mm):,  %12.5f\n",TotalInfilVol/CatchmentArea*1000.0);
        fprintf(fp," Average surface storage  (mm):,  %12.5f\n",SumDeprStoreVol/CatchmentArea*1000.0);
        fprintf(fp," Water in runoff          (mm):,  %12.5f\n",(SumWaterHVolout-SumDeprStoreVol)/CatchmentArea*1000.0);
        fprintf(fp," Mass balance error (water)(%%):,  %12.5f\n",MassBalanceError);
        fprintf(fp,"Total discharge           (m3):, %12.5f\n",TotalDischarge);
        fprintf(fp,"Peak discharge           (l/s):,  %12.5f\n",PeakDischarge);
//VJ 080217 add baseflow
//        fprintf(fp,"Baseflow discharge           (l/s):,  %12.5f\n",BaseflowDischarge);
        fprintf(fp,"Peak time                (min):,  %12.5f\n",PeakTime);
      if (TotalRainVol > 0)
        fprintf(fp,"Discharge/Rainfall         (%%):,  %10.3f  \n",TotalDischarge/TotalRainVol*100);
      else
        fprintf(fp,"Discharge/rainfall         (%%):,  %10.3f  \n",0.0);
        fprintf(fp,"Splash detachment        (ton):,  %12.5f\n",TotalSplashDetachment/1000.0);
        fprintf(fp,"Flow detachment (land)   (ton):,  %12.5f\n",TotalFlowDetachment/1000.0);
        fprintf(fp,"Deposition (land)        (ton):,  %12.5f\n",TotalDeposition/1000.0);
        //VJ 050822 what about gullies?????
        fprintf(fp,"Erosion channel/wheeltr. (ton):,  %12.5f\n",(TotalChannelFlowDetachment+TotalWheelFlowDetachment)/1000.0);
        fprintf(fp,"Deposition channel/whlt  (ton):,  %12.5f\n",(TotalChannelDeposition+TotalWheelDeposition)/1000.0);
        fprintf(fp,"Suspended Sediment       (ton):,  %12.5f\n",OutputSedSusp/1000.0);
        fprintf(fp,"Susp. Sediment chan/whlt (ton):,  %12.5f\n",OutputChanSedSusp/1000.0);
        fprintf(fp,"Total soil loss          (ton):,  %12.5f\n",TotalSedDischarge/1000.0);
        fprintf(fp,"Average soil loss      (kg/ha):,  %12.5f\n",TotalSedDischarge/(CatchmentArea/10000.0));

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
//VJ 030415 added nutrient det and dep
      }
        //fprintf(fp,"---------------------------------------------\n");
        fclose(fp);


