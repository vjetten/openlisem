// *********************************************************
// ******** MULTICLASS SEDIMENT ****************************
// *********************************************************

    if (!SwitchNoErosion)
    {

       // ********************************************************************
       // ******** Transport capacity  ***************************************
       // ********************************************************************

        // **** transport capacity of the soil material ****
        _spatial(REAL4, TCsoil);
        _spatial(REAL4, StreamPower);
        calc(" StreamPower = Gradient*V*100 ");

        calc(" TCsoil = mif(StreamPower gt CriticalStreamPower, "
               " CGsoil*(StreamPower - CriticalStreamPower)**DGsoil, 0) ");
        calc(" TCsoil = 2650*min(TCsoil,0.32) ");
        calc(" TransportCapOUT = TCsoil ");
         //sum TC based on soil fraction
         // copied here for timeseries output

        _spatial(REAL4, TCSoilMu0);
        _spatial(REAL4, TCSoilMu1);
        _spatial(REAL4, TCSoilMu2);
        _spatial(REAL4, TCSoilMu3);
        _spatial(REAL4, TCSoilMu4);
        _spatial(REAL4, TCSoilMu5);
        mcalc(" TCSoilMu# = TCsoil * FractionMu# ",6);


        // ********************************************************************
        // ******** Add splash ************************************************
        // ********************************************************************

        calc(" DepositionSplash = mif(TCsoil == 0, -DetachmentSplash, 0) ");
        mcalc(" SedinMu# += (DetachmentSplash+DepositionSplash)*FractionMu# ",6);


        // *******************************************************************
        // ******** concentration                  ***************************
        // *******************************************************************

         mcalc(" SedConc# = mif(WaterHVolin gt 0 ,SedinMu#/WaterHVolin,0) ",6);

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
         mcalc(" DF# = Y*max(TCSoilMu#-SedConc#,0)",6);
   //      mcalc(" DF# *= mif(WH gt 0, (1-exp(-SV#/(0.001*WH)*DTSEC))*WaterHVolin*RoadWidthFactor, 0)",6);
         mcalc(" DF# *= SV#*DTSEC*DXc*FlowWidth*RoadWidthFactor",6);
         mcalc(" DF# = min(DF#,max(0,TCSoilMu#-SedConc#)*WaterHVolin)",6);

         if (SwitchGrassPresent)
         {
            mcalc(" DF# = mif(GrassWidth gt 0, (1-GrassFraction)*DF#, DF#) ",6);
                //no erosion on grasstrip
         }

        // **** new concentration ****

//         mcalc(" SedinMu# += DF# ",6);
//         mcalc(" SedConc# = mif(WaterHVolin gt 0 ,SedinMu#/WaterHVolin,0) ",6);

        //================
        //***DEPOSITION***
        //================

        // **** transport capacity of the suspended material ****
        // *** based on D50 of the suspended mat.

//        _spatial(REAL4, D50susp);   moved to paraminput for output

        findD50(D50susp,SedinMu0,SedinMu1,SedinMu2,SedinMu3,SedinMu4,SedinMu5,*mu_cl);

        _spatial(REAL4, CGsusp);
        _spatial(REAL4, DGsusp);
        calc(" CGsusp = ((D50susp+5)/0.32)**-0.6 ");
        calc(" DGsusp = ((D50susp+5)/300)**0.25 ");

        _spatial(REAL4, TCsusp);
        calc(" TCsusp = mif(StreamPower gt CriticalStreamPower, "
              " CGsusp*(StreamPower - CriticalStreamPower)**DGsusp, 0) ");
        calc(" TCsusp = 2650*min(TCsusp,0.32) ");
        calc(" TCsusp = mif(D50susp == 0, 0, TCsusp)");

        _spatial(REAL4, TCSuspMu0);
        _spatial(REAL4, TCSuspMu1);
        _spatial(REAL4, TCSuspMu2);
        _spatial(REAL4, TCSuspMu3);
        _spatial(REAL4, TCSuspMu4);
        _spatial(REAL4, TCSuspMu5);
        _spatial(REAL4, SedinMuSum);
        calc(" SedinMuSum = 0 ");
        calc(" SedinMuSum += SedinMu0");
        calc(" SedinMuSum += SedinMu1");
        calc(" SedinMuSum += SedinMu2");
        calc(" SedinMuSum += SedinMu3");
        calc(" SedinMuSum += SedinMu4");
        calc(" SedinMuSum += SedinMu5");

        mcalc(" TCSuspMu# = mif(SedinMuSum gt 0,TCsusp*SedinMu#/SedinMuSum,0) ",6);

        _spatial(REAL4, Dep0);
        _spatial(REAL4, Dep1);
        _spatial(REAL4, Dep2);
        _spatial(REAL4, Dep3);
        _spatial(REAL4, Dep4);
        _spatial(REAL4, Dep5);
        mcalc(" Dep# = min(TCSuspMu#-SedConc#, 0)",6);
        mcalc(" Dep# *= mif(WH gt 0, (1-exp(-DTSEC*SV#/(0.001*WH)))*WaterHVolin, 1)",6);
      //        mcalc(" Dep# *= SV#*DTSEC*DX*FlowWidth",6);
        mcalc(" Dep# = max(Dep#, min(0,TCSuspMu#-SedConc#)*WaterHVolin) ",6);
           // no more deposition then TC surplus

        if (SwitchGrassPresent)
        {
           mcalc(" Dep# = mif(GrassWidth gt 0, (1-GrassFraction)*Dep#-GrassFraction*SedinMu#,Dep#) ",6);
             // total deposition on grasstrip
        }


        mcalc(" DF# = mif(RunoffMeanHin gt MinimumHeight, DF#, 0) ",6);
        mcalc(" Dep# = mif(RunoffMeanHin gt MinimumHeight, Dep#, -SedinMu#) ",6);

        mcalc(" SedinMu# += DF# ",6);
        mcalc(" SedinMu# += Dep# ",6);

        mcalc(" SedinMu# = max(SedinMu#, 0)",6);

        if (SwitchNutrients)
        {
            // deposition nut as fraction of total sediment
            // erosion as uptake * local content 
            _spatial(REAL4, NutDepo);
            calc(" NutDepo = mif(SedinMu0 gt 0.001, NutPSuspension*Dep0/SedinMu0, 0) ");
            calc(" NutPSuspension = max(0, NutPSuspension + DF0*NutPContent*NutPConversion + NutDepo) ");

            calc(" NutDepo = mif(SedinMu0 gt 0.001, NutNO3Suspension*Dep0/SedinMu0, 0) ");
            calc(" NutNO3Suspension = max(0,NutNO3Suspension+ DF0*NutNO3Content*NutNO3Conversion + NutDepo) ");

            calc(" NutDepo = mif(SedinMu0 gt 0.001, NutNH4Suspension*Dep0/SedinMu0, 0) ");
            calc(" NutNH4Suspension = max(0, NutNH4Suspension+DF0*NutNH4Content*NutNH4Conversion + NutDepo) ");
        }

    }//not noerosion
