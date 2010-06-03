//these variables link the output of lisem to the screen, in units that are showed on the screen
//the purpose is to make the lisem model and the interfce as seperate as possible so that a
//change of model does not affect the interface too much and vice versa 

       StartTime = STARTINTERVAL;
       EndTime = ENDINTERVAL;
       CatchmentAreaHa = CatchmentArea/10000.0;

       CurrentTime = timestepindex;
       CellSize = DX;
       AVGRainMM = RainfallAverageH*3600.0/DTSEC;
       // used for line graph
       TOTRainMM = TotalRainVol/CatchmentArea*1000.0;
       TOTInterceptionMM = TotalInterceptionVol/CatchmentArea*1000.0;
       TOTSurfStorMM = SumDeprStoreVol/CatchmentArea*1000.0;
       TOTInfilMM = TotalInfilVol/CatchmentArea*1000.0;
       TOTRunoffMM = (SumWaterHVolout+SumChannelVolout+SumWheelVolout)/CatchmentArea*1000.0;
       // SumWaterHVolout includes depression storage
       TOTDischargeMM = TotalDischarge/CatchmentArea*1000.0;
       TOTDischargeM3 = TotalDischarge;
       if (TotalRainVol > 0)
          P_QPercentage = TotalDischarge/TotalRainVol*100;
       else
          P_QPercentage = 0;
       PeakDischargeQ = PeakDischarge;
       //baseflow estimate cannot be done because all is mixed, inflection point analysis!
       BaseflowDischargeQ = BaseflowDischarge;
       PeakDischargeTime = PeakTime;
       PeakRainfallTime = PeakRainTime;
       TOTBufferVolume = SumBufferVolume + SumBufferVolumeChannel;
       TOTBufferSedVolume = SumBufferSedVolume/1000.0 + SumBufferSedVolumeChannel/1000.0;
       MassBalance = MassBalanceError;
       SedMassBalance = SedMassBalanceError;
       OutputSedSusp = SumSedout/1000.0;
       //VJ 050822 what about gullies
       OutputChanSedSusp = (SumWheelSedout+SumChannelSedout)/1000.0;
       TOTSplashErosion = TotalSplashDetachment/1000.0;
       TOTFlowErosion = TotalFlowDetachment/1000.0;
       TOTDeposition = TotalDeposition/1000.0;
       TOTChanFlowErosion = (TotalWheelFlowDetachment+TotalChannelFlowDetachment)/1000.0;
       TOTChanDeposition = (TotalWheelDeposition+TotalChannelDeposition)/1000.0;
       TOTSoilLoss = TotalSedDischarge/1000.0;
       AVGSoilLoss = TotalSedDischarge/(CatchmentArea/10000.0);
       OutputDischarge = QOutflow*1000.0;
//       OutputDischarge1 = QOutlet2*1000.0;
//       OutputDischarge2 = QOutlet3*1000.0;

       if (QOutflow > 0)
          OutputSedConc = SCOutlet1;
       else
          OutputSedConc = 0;
       if (QOutflow > 0)
          OutputSediment = SedOutflow;
          //in multised this is the total calc in sedmassball
       else
          OutputSediment = 0;

       StepCounter = stepNr;

       //functions called by synchronize, safe method to call a function from a thread
       //functions InitializeTotalsSync(void) and  UpdateTotalsSync(void) are in lisscreenoutput.cpp

       if (!InitDone)
       {
          Mout = &QoutOut;
          Synchronize(DrawMapInitSync);
       }
       else
       {
          switch (LisIFace->GroupDisplayMap->ItemIndex)
          {
           case 0: Mout = &QoutOut;  break;
           case 1: Mout = &WHOut;  break;
           case 2: Mout = &VOut;  break;
           case 3: Mout = &InfilOut;  break;
           case 4: Mout = &SoilLossOut; break;
           case 5: Mout = &DetOut; break;
           case 6: Mout = &DepOut; break;
           case 7: Mout = &TCOut; break;
           case 8: Mout = &SCOut; break;
         }
         Synchronize(DrawMapSync);
       }

       if (!InitDone)
          Synchronize(InitializeTotalsSync);
       Synchronize(UpdateTotalsSync);




