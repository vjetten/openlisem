     /**** vars for TOTALS, syntax:
     Totalxxx        is used for cumulative non spatial total for whole run, totaldischarge
     Sumxxx          is used for spatial total in one timestep, e,g, sumwatervolout
     SumTotxxx       is used for spatial total of run, e.g. erosion.map 
     */

     _nonspatial(REAL4, MinimumHeight);
     MinimumHeight = 0.001;
        //0.01 mm is minimum water height acceptable to avoid divisaion by near zero

     _nonspatial(REAL4, CatchmentArea);   // total catchment area in m2
     calc(" CatchmentArea = sum(DXc * DX) ");
     // using real slope

     _spatial(REAL4, RainHCum);
     calc(" RainHCum = 0.0 ");
     // cumulative rainfall per cell in mm
//VJ030415 changed names to cum

     _nonspatial(REAL4, TotalRainVol);
     TotalRainVol = 0.0;
        // total rainfall in total catchment area, in m3

     _spatial(REAL4, InterceptionHCum);
     calc("InterceptionHCum = 0.0");
     //cumulative interception per cell in mm
//VJ030415 changed names to cum

     _spatial(REAL4, WHCum);
     calc("WHCum = 0.0");
//VJ 100115 cumulative runoff


     _nonspatial(REAL4, TotalInterceptionVol);
     TotalInterceptionVol = 0.0;
        // total interception in total catchment area, in m3

     _nonspatial(REAL4, SumDeprStoreVol);
        // non spatial sum of water in surface storage in one DT

     _nonspatial(REAL4, TotalInfilVol);
     TotalInfilVol = 0.0;
        // total m3 of water infiltrated, nonspatial

     _spatial(REAL4, InfilVol);
     calc(" InfilVol = 0.0 ");
      // total m3 of water infiltrated, spatial

     _spatial(REAL4, InfilKinWave);
     //VJ 050831 aux variable needed for calc of infil in kin wave

     _nonspatial(REAL4, SumWaterHVolout);
     SumWaterHVolout = 0;

     _nonspatial(REAL4, SumWaterHVolin);
     SumWaterHVolin = 0;

     _nonspatial(REAL4, SumSedout);
     SumSedout = 0;

     _nonspatial(REAL4, SumChannelVolout);
     SumChannelVolout = 0;

     _nonspatial(REAL4, SumChannelSedout);
     SumChannelSedout = 0;

     _nonspatial(REAL4, SumWheelVolout);
     SumWheelVolout = 0.0;

     _nonspatial(REAL4, SumWheelSedout);
     SumWheelSedout = 0.0;

     _nonspatial(REAL4, SumGullyVolout);
     SumGullyVolout = 0;

     _nonspatial(REAL4, SumGullySedout);
     SumGullySedout = 0;


       // sum of map water and sed in channel before and after kin wave
       // used in mass balance error and corrections

     _nonspatial(REAL4, TotalDischarge);
     TotalDischarge = 0.0;
        // total discharge from catchment, in m3

     _nonspatial(REAL4, PeakDischarge);
     PeakDischarge = 0.0;
        // peak discharge from catchment, in l/s

//VJ 080217 add baseflow        
     _nonspatial(REAL4, BaseflowDischarge);
     BaseflowDischarge = 0.0;
        // baseflow discharge from catchment, in l/s

     _nonspatial(REAL4, TotalBaseflowVol);
     if (SwitchChannelBaseflow)
      calc(" TotalBaseflowVol = sum(ChannelVolin) ");
     else
      calc(" TotalBaseflowVol = 0 ");
     // total initial volume in channel, ChannelVolin is calculated in lisparamchannel.cpp

     _nonspatial(REAL4, PeakTime);
     PeakTime = 0.0;
        // time of peak discharge, in min

     _nonspatial(REAL4, RainfallAverageH);
     RainfallAverageH = 0.0;
     _nonspatial(REAL4, PeakRainfall);
     PeakRainfall = 0.0;
     _nonspatial(REAL4, PeakRainTime);
     PeakRainTime = 0.0;
        // time of peak rainfall, in min

     _nonspatial(REAL4, TotalSedDischarge);
     TotalSedDischarge = 0.0;
        // total soil loss from catchment, in kg

     _nonspatial(REAL4, TotalSplashDetachment);
     TotalSplashDetachment = 0.0;
        // total rainfall detachment on land

     _nonspatial(REAL4, TotalFlowDetachment);
     TotalFlowDetachment = 0.0;
        // total flow detachment on land

     _nonspatial(REAL4, TotalDeposition);
     TotalDeposition = 0.0;
        // total deposition on land

     _nonspatial(REAL4, TotalChannelFlowDetachment);
     TotalChannelFlowDetachment = 0.0;
        // total flow detachment in channel

     _nonspatial(REAL4, TotalChannelDeposition);
     TotalChannelDeposition = 0.0;
        // total deposition in channel

     _nonspatial(REAL4, TotalWheelFlowDetachment);
     TotalWheelFlowDetachment = 0.0;
        // total flow detachment in wheel track

     _nonspatial(REAL4, TotalWheelDeposition);
     TotalWheelDeposition = 0.0;
        // total deposition in wheel track

     _nonspatial(REAL4, TotalGullyFlowDetachment);
     TotalGullyFlowDetachment = 0.0;
        // total flow detachment in gullies

     _nonspatial(REAL4, TotalGullyDeposition);
     TotalGullyDeposition = 0.0;
        // total deposition in gullies

     _spatial(REAL4, SumTotErosion);
     calc(" SumTotErosion = 0.0 ");
     // spatial sum of all detachment

     _spatial(REAL4, SumTotDeposition);
     calc(" SumTotDeposition = 0.0 ");
     // spatial sum of all deposition in one timstep, used for outputmap

     _nonspatial(REAL4, MassBalanceError);
     MassBalanceError = 0.0;
        // water mass balance

     _nonspatial(REAL4, SedMassBalanceError);
     SedMassBalanceError = 0.0;
        // sediment mass balance

//############  TOTALS for multiclass ###########     
 
     _nonspatial(REAL4, TotalSedinMu0);
     TotalSedinMu0 = 0.0;
     _nonspatial(REAL4, TotalSedinMu1);
     TotalSedinMu1 = 0.0;
     _nonspatial(REAL4, TotalSedinMu2);
     TotalSedinMu2 = 0.0;
     _nonspatial(REAL4, TotalSedinMu3);
     TotalSedinMu3 = 0.0;
     _nonspatial(REAL4, TotalSedinMu4);
     TotalSedinMu4 = 0.0;
     _nonspatial(REAL4, TotalSedinMu5);
     TotalSedinMu5 = 0.0;
     //map totals 
     
     _nonspatial(REAL4, TotalSedDischargeMu0);
     TotalSedDischargeMu0 = 0.0;
     _nonspatial(REAL4, TotalSedDischargeMu1);
     TotalSedDischargeMu1 = 0.0;
     _nonspatial(REAL4, TotalSedDischargeMu2);
     TotalSedDischargeMu2 = 0.0;
     _nonspatial(REAL4, TotalSedDischargeMu3);
     TotalSedDischargeMu3 = 0.0;
     _nonspatial(REAL4, TotalSedDischargeMu4);
     TotalSedDischargeMu4 = 0.0;
     _nonspatial(REAL4, TotalSedDischargeMu5);
     TotalSedDischargeMu5 = 0.0;
     // discharge totals

//############  TOTALS for nutrients ###########
// solution and suspendion are leaving the catchment at the outlet (sum of discharge)
// infiltration, deposition and detachment are totals of the cumulative maps, kg/m2
     _nonspatial(REAL4, TotalPSolution);
     TotalPSolution = 0.0;
     _nonspatial(REAL4, TotalPSuspension);
     TotalPSuspension = 0.0;
     _nonspatial(REAL4, TotalPInfiltration);
     TotalPInfiltration = 0.0;
     _nonspatial(REAL4, TotalPDeposition);
     TotalPDeposition = 0.0;
     _nonspatial(REAL4, TotalPDetachment);
     TotalPDetachment = 0.0;
 
     _nonspatial(REAL4, TotalNO3Solution);
     TotalNO3Solution = 0.0;
     _nonspatial(REAL4, TotalNO3Suspension);
     TotalNO3Suspension = 0.0;
     _nonspatial(REAL4, TotalNO3Infiltration);
     TotalNO3Infiltration = 0.0;
     _nonspatial(REAL4, TotalNO3Deposition);
     TotalNO3Deposition = 0.0;
     _nonspatial(REAL4, TotalNO3Detachment);
     TotalNO3Detachment = 0.0;
 
     _nonspatial(REAL4, TotalNH4Solution);
     TotalNH4Solution = 0.0;
     _nonspatial(REAL4, TotalNH4Suspension);
     TotalNH4Suspension = 0.0;
     _nonspatial(REAL4, TotalNH4Infiltration);
     TotalNH4Infiltration = 0.0;
     _nonspatial(REAL4, TotalNH4Deposition);
     TotalNH4Deposition = 0.0;
     _nonspatial(REAL4, TotalNH4Detachment);
     TotalNH4Detachment = 0.0;
     // map totals
     //VJ 030415 added erosion and deposition
     
		_spatial(REAL4, SumTotNutPdep  );     
		_spatial(REAL4, SumTotNutNH4dep);     
		_spatial(REAL4, SumTotNutNO3dep);     
		_spatial(REAL4, SumTotNutPdet  );     
		_spatial(REAL4, SumTotNutNH4det);     
		_spatial(REAL4, SumTotNutNO3det);     
		calc("SumTotNutPdep   = 0 ");
		calc("SumTotNutNH4dep = 0 ");
		calc("SumTotNutNO3dep = 0 ");
		calc("SumTotNutPdet   = 0 ");
		calc("SumTotNutNH4det = 0 ");
		calc("SumTotNutNO3det = 0 ");
      //sum of the run for output

//VJ 040823 include buffers      
     _nonspatial(REAL4, SumBufferVolume);
     SumBufferVolume = 0;
     _nonspatial(REAL4, SumBufferVolumeChannel);
     SumBufferVolumeChannel = 0;
     _nonspatial(REAL4, SumBufferSedVolume);
     SumBufferSedVolume = 0;
     _nonspatial(REAL4, SumBufferSedVolumeChannel);
     SumBufferSedVolumeChannel = 0;

