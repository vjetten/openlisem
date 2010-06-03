// *********************************************************
// ******** MULTICLASS SEDIMENT ****************************
// *********************************************************

    if (!SwitchNoErosion)
    {

         // *****************************************************************************
         // ******** Transport capacity and settling velocity ****************************
         // *****************************************************************************

           // **** transport capacity of overland flow ****
           _spatial(REAL4, TCMu0);
           _spatial(REAL4, TCMu1);
           _spatial(REAL4, TCMu2);
           _spatial(REAL4, TCMu3);
           _spatial(REAL4, TCMu4);
           _spatial(REAL4, TCMu5);



           _spatial(REAL4, StreamPower);
           calc(" StreamPower = Gradient*V*100 ");

           mcalc(" TCMu# = mif(StreamPower gt CriticalStreamPower, "
                 " CG#*(StreamPower - CriticalStreamPower)**DG#, 0) ",6);
           mcalc(" TCMu# = TCCAL * 2650*min(TCMu#,0.32) ",6);
           calc(" TransportCapOUT = 0 ");
           mcalc(" TransportCapOUT += TCMu# * FractionMu# ",6);
           //sum TC based on soil fraction
          // copied here for timeseries output

           _spatial(REAL4, TCSoilMu0);
           _spatial(REAL4, TCSoilMu1);
           _spatial(REAL4, TCSoilMu2);
           _spatial(REAL4, TCSoilMu3);
           _spatial(REAL4, TCSoilMu4);
           _spatial(REAL4, TCSoilMu5);
          mcalc(" TCSoilMu# = TCMu# * FractionMu# ",6);

           _spatial(REAL4, TCSuspMu0);
           _spatial(REAL4, TCSuspMu1);
           _spatial(REAL4, TCSuspMu2);
           _spatial(REAL4, TCSuspMu3);
           _spatial(REAL4, TCSuspMu4);
           _spatial(REAL4, TCSuspMu5);
           _spatial(REAL4, SedinMuSum);
           calc(" SedinMuSum = 0 ");
           mcalc(" SedinMuSum += SedinMu#",6);
           mcalc(" TCSuspMu# = mif(SedinMuSum gt 0, TCMu# * SedinMu#/SedinMuSum, 0) ",6);

           mcalc(" SedConc# = mif(WaterHVolin gt 0 ,SedinMu#/WaterHVolin,0) ",6);
//           calc(" SedinMuSum = 0 ");
//           mcalc(" SedinMuSum += SedConc#",6);
  //BULLSHIT

// *****************************************************************************
// ******** FLOW detachment and deposition *************************************
// *****************************************************************************

           _spatial(REAL4, RoadWidthFactor);
           calc(" RoadWidthFactor = mif(FlowWidth gt 0,1.0-RoadWidthDX/FlowWidth,0) ");
              // roadwidthfactor: no flow detachment on roads
           _spatial(REAL4, DF0);
           _spatial(REAL4, DF1);
           _spatial(REAL4, DF2);
           _spatial(REAL4, DF3);
           _spatial(REAL4, DF4);
           _spatial(REAL4, DF5);
           mcalc(" DF# = max(Y*(TCSoilMu#-SedConc#),0)",6);
//           mcalc(" DF# *= mif(WH gt 0, (1-exp(-SV#/(0.001*WH)*DTSEC))*WaterHVolin*RoadWidthFactor, 0)",6);
	       if (SwitchAltErosion)
    	   {
	    	   mcalc(" DF# *= Qin*DTSEC ");
	       }
    	   else
           {
	           mcalc(" DF# *= SV#*DTSEC*DX*FlowWidth*RoadWidthFactor",6);
           }

           mcalc(" DF# = min(DF#,max(0,(TCSoilMu#-SedConc#))*WaterHVolin)",6);

           _spatial(REAL4, Dep0);
           _spatial(REAL4, Dep1);
           _spatial(REAL4, Dep2);
           _spatial(REAL4, Dep3);
           _spatial(REAL4, Dep4);
           _spatial(REAL4, Dep5);
//           mcalc(" Dep# = mif(TCSoilMu#+TCSuspMu# le SedConc#, (TCSuspMu#-SedConc#), 0)",6);
//           mcalc(" Dep# = mif(SumTCMu le SedinMuSum, (TCSuspMu#-SedConc#), 0)",6);
           mcalc(" Dep# = min(TCSuspMu#-SedConc#, 0)",6);
           mcalc(" Dep# *= mif(WH gt 0, (1-exp(-SV#/(0.001*WH)*DTSEC))*WaterHVolin, 0)",6);
//           mcalc(" Dep# *= SV#*DTSEC*DX*FlowWidth",6);
           mcalc(" Dep# = mif(DF# gt 0, 0, Dep#) ",6);
             // no deposition is there is detachment

           mcalc(" Dep# = max(Dep#, min(0,TCSuspMu#-SedConc#)*WaterHVolin) ",6);

          if (SwitchGrassPresent)
          {
             mcalc(" DF# = mif(GrassWidth gt 0, (1-GrassFraction)*DF#, DF#) ",6);
                //no erosion on grasstrip
             mcalc(" Dep# = mif(GrassWidth gt 0, (1-GrassFraction)*Dep#-GrassFraction*SedinMu#,Dep#) ",6);
               // total deposition on grasstrip

          }

   //       mcalc(" DF# = mif(RunoffMeanHin gt MinimumHeight, DF#, 0) ",6);
   //       mcalc(" Dep# = mif(RunoffMeanHin gt MinimumHeight, Dep#, -SedinMu#) ",6);


//          calc(" DepositionSplash =  ");
          //mcalc(" DepositionSplash += mif(TCMu# == 0, -DetachmentSplash*FractionMu#, 0) ",6);
          calc(" DepositionSplash = mif(TransportCapOUT == 0, -DetachmentSplash, 0) ");

          mcalc(" SedinMu# += (DetachmentSplash+DepositionSplash)*FractionMu# ",6);

          mcalc(" SedinMu# += DF# ",6);
          mcalc(" SedinMu# += Dep# ",6);

          mcalc(" SedinMu# = max(SedinMu#, 0)",6);

          }//not noerosion
