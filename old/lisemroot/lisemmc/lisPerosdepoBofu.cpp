// *****************************************************************************
// ******** FLOW detachment **************************************************
// *****************************************************************************


    if (!SwitchNoErosion)
    {
        _spatial(REAL4, FluidDensity);
        calc(" FluidDensity = 1000 ");
        _spatial(REAL4, Omega);
        _spatial(REAL4, OmegaCrit);
        _spatial(REAL4, R);
        calc(" R = mif(FlowWidth gt 0, FlowWidth*(RunoffMeanHin/1000)/(FlowWidth+2*RunoffMeanHin/1000),0) ");
        calc(" Omega = 9.8*Gradient*FluidDensity*V*R ");
        calc(" OmegaCrit = 9.8*FluidDensity*R*0.004 ");
//writeTimeseries(OmegaCrit,"omega");

       // **** sediment concentration ****
       calc(" SedConcentration = mif( WaterHVolin gt 0, Sedin/WaterHVolin, 0) ");
          // sediment concentration in kg/m3
         _spatial(REAL4, F);
         calc(" F = 0.05");

         _spatial(REAL4, c_t);
        calc(" c_t = mif(RunoffMeanHin gt 0, "
             "  2650 * F * max(Omega-OmegaCrit,0)/((2650-FluidDensity)*"
             "  9.87*RunoffMeanHin/1000*SettlingVelocity*10),0) ");
 //       writeTimeseries(c_t,"ct");
 //       writeTimeseries(SedConcentration,"c");

        _spatial(REAL4, H);
         calc(" H = 0.1 ");
         _spatial(REAL4, term1);
         calc(" term1 = (1-H)*Y*F*max(Omega-OmegaCrit,0)");
         _spatial(REAL4, term2);
         calc(" term2 = mif(RunoffMeanHin gt 0,F*H*2650/((2650-FluidDensity)*9.8*RunoffMeanHin/1000)*max(Omega-OmegaCrit,0),0)");


        _spatial(REAL4, dSedin);
//        calc(" dSedin = mif(c_t gt 0, DXc * Y* F* max(Omega-OmegaCrit,0)*(1-SedConcentration/c_t), 0)");
        calc(" dSedin = DXc*(term1 + term2 - SedConcentration*SettlingVelocity)");



 //       writeTimeseries(dSedin,"ds");

       _spatial(REAL4, DetachmentFlow);
       _spatial(REAL4, DepositionFlow);
       calc(" DetachmentFlow = max(dSedin, 0)");
       calc(" DepositionFlow = min(dSedin, 0)");

       calc(" DepositionSplash = mif(Omega le OmegaCrit, -DetachmentSplash, 0)");
          // no streampower no splash removal
          
       calc(" Sedin += DetachmentSplash ");
       calc(" Sedin += DepositionSplash ");
       calc(" Sedin += DetachmentFlow ");
       calc(" Sedin += DepositionFlow ");


/*
       _spatial(REAL4, streampower);
       calc(" streampower = 100*Gradient*V ");
       _spatial(REAL4, TransportCapacity);
       calc(" TransportCapacity = mif(streampower gt CriticalStreamPower, "
            " CGovers*(streampower-CriticalStreamPower)**(DGovers), 0) ");
          // TransportCapacity is the volumetric transport capacity of sediment
          // in (cm3 soil)/(cm3 water)
          // V is flow rate in m/s; formula requires cm/s, thus factor 100

       calc(" TransportCapacity = 2650*min(TransportCapacity,0.32) ");
          // 2650 is the particle density in kg/m3
          // the unit of TransportCapacity is now in kg/m3

       calc(" TransportCapOUT = TransportCapacity ");
          // copied here for timeseries output

       // **** sediment concentration ****
       calc(" SedConcentration = mif( WaterHVolin gt 0, Sedin/WaterHVolin, 0) ");
          // sediment concentration in kg/m3

       // **** settling velocity ****

       _spatial(REAL4, RoadWidthFactor);
       calc(" RoadWidthFactor = mif(FlowWidth gt 0,1.0-RoadWidthDX/FlowWidth,0) ");
          // roadwidthfactor: no flow detachment on roads

       _spatial(REAL4, TransportFactor);
       calc("TransportFactor = DTSEC*SettlingVelocity*DX*FlowWidth ");
       calc("TransportFactor = mif(Qin gt 0,(1-exp(-SettlingVelocity/Qin*DXc*FlowWidth)),0)");

       _spatial(REAL4, DetachmentFlow);
       calc(" DetachmentFlow = Y*max(TransportCapacity-SedConcentration,0)*TransportFactor*RoadWidthFactor");


       if (SwitchGrassPresent)
         calc(" DetachmentFlow = mif(GrassWidth gt 0, (1-GrassFraction)*DetachmentFlow,DetachmentFlow) ");
          // no flow detachment on field strips and waterways

       calc(" DetachmentFlow = min(DetachmentFlow, max(0,TransportCapacity-SedConcentration)*WaterHVolin)");
          // no more detachment than there is capacity to transport

          // the unit of TransportCapacity is now in kg/m3

//       calc("TransportFactor = mif(WH gt 0,(1-exp(-DTSEC*SettlingVelocity/(0.001*WH)))*WaterHVolin,0)");
       calc("TransportFactor = mif(Qin gt 0,(1-exp(-SettlingVelocity/Qin*DXc*FlowWidth)),0)");
writeTimeseries(TransportFactor,"tf"); 
       _spatial(REAL4, DepositionFlow);
       calc(" DepositionFlow = min(0, TransportCapacity-SedConcentration)*TransportFactor");

       if (SwitchGrassPresent)
         calc(" DepositionFlow = mif(GrassWidth gt 0, "
              "-Sedin*GrassFraction+(1-GrassFraction)*DepositionFlow, DepositionFlow) ");
         // total deposition on grassstrips

//          calc(" DepositionFlow = max(DepositionFlow, -Sedin) ");
       calc(" DepositionFlow = max(DepositionFlow, min(0,TransportCapacity-SedConcentration)*WaterHVolin) ");
          // Deposition amount (kg)
          // (positive = supply to the flow; negative = deposition)
          // DepositionFlow is the amount of deposition (negative!)
          // deposition cannot be larger than the amount of available sediment
          // =max because of negative values
          // if the flow is 0, than all sediment is deposited


// *****************************************************************************
// ******** Combine all fluxes for Sedin *********************************
// *****************************************************************************

       calc(" DetachmentFlow = mif(RunoffMeanHin gt MinimumHeight, DetachmentFlow, 0) ");
       calc(" DepositionFlow = mif(RunoffMeanHin gt MinimumHeight, DepositionFlow, -Sedin) ");

       calc(" DepositionSplash = mif(TransportCapacity == 0, -DetachmentSplash, 0) ");
          // none delivered splash is counted with deposition

       calc(" DetachmentSplash = 0 ");
       calc(" DepositionSplash = 0 ");

       calc(" DetachmentFlow = 0 ");
       calc(" DepositionFlow = 0");
*/

          // Sedin is the amount of sediment in kg available for transport
          // the available sediment cannot be smaller than 0
          // if sedin < 0 add to deposition
    } //not no erosion


