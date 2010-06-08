// *****************************************************************************
// ******** FLOW detachment **************************************************
// *****************************************************************************


    if (!SwitchNoErosion)
    {
       _spatial(REAL4, streampower);
       calc(" streampower = 100*Gradient*V ");
//VJ 050303 V was not recalculated in lisPkinematic so questionable if streampower
//is good in channel cells

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

//       calc(" DepositionSplash = mif(TransportCapacity == 0, -DetachmentSplash, 0) ");
//==> VJ 100508 waterheight minimal not TC ???
// VJ 100109 replaced by sed concentration check
       calc(" Sedin += DetachmentSplash ");
//       calc(" Sedin += DepositionSplash ");
       _spatial(REAL4, sedvoltemp);
       calc(" SedConcentration = mif( WaterHVolin gt 0, Sedin/WaterHVolin, 1000) ");
       calc(" SedConcentration = min(SedConcentration, 848) ");
       calc(" sedvoltemp = SedConcentration*WaterHVolin ");
       calc(" DepositionSplash += min(0, sedvoltemp-Sedin) ");
       calc(" Sedin = WaterHVolin * SedConcentration ");

       // **** sediment concentration ****
//       calc(" SedConcentration = mif( WaterHVolin gt 0, Sedin/WaterHVolin, 0) ");
          // sediment concentration in kg/m3
//==> VJ 100508 check max conc < 848 ???

       // **** settling velocity ****

 //      _spatial(REAL4, RoadWidthFactor);
 //      calc(" RoadWidthFactor = mif(FlowWidth gt 0,1.0-RoadWidthDX/FlowWidth,0) ");
          // roadwidthfactor: no flow detachment on roads

       _spatial(REAL4, TransportFactor);
      // calc("TransportFactor = DTSEC*SettlingVelocity*DXc*FlowWidth ");
       calc("TransportFactor = DTSEC*SettlingVelocity*DXc*SoilWidthDX*PondAreaFract ");
       //VJ 100603 changed: this way roads are ignored, like in openlisem 

       _spatial(REAL4, DetachmentFlow);
       if (SwitchAltErosion)
       {
	       calc(" DetachmentFlow = Y*max(0,TransportCapacity-SedConcentration)*Qin*DTSEC ");
       }
       else
       {
	       calc(" DetachmentFlow = Y*max(TransportCapacity-SedConcentration,0)*TransportFactor");//*RoadWidthFactor");
       }
//==> VJ 100508 do roadwidth directly ???


   	   calc(" DetachmentFlow = min(DetachmentFlow, max(0,TransportCapacity-SedConcentration)*WaterHVolin)");
          // no more detachment than there is capacity to transport

       if (SwitchGrassPresent)
         calc(" DetachmentFlow = mif(GrassWidth gt 0, (1-GrassFraction)*DetachmentFlow,DetachmentFlow) ");
          // no flow detachment on field strips and waterways

       if (SwitchNoErosionOutlet)
       {
          calc(" DetachmentFlow = mif(Outlet1 eq 1, 0,DetachmentFlow) ");
       }
//VJ 040224 added if outlet no detachment

       calc("TransportFactor = mif(WH gt 0,(1-exp(-DTSEC*SettlingVelocity/(0.001*WH)))*WaterHVolin,WaterHVolin)");
       //MOET 1 ZIJN ALS WH = 0 -> voledige depositie!
//==> VJ 100508 not 1 but WaterHVolin ???
       

       _spatial(REAL4, DepositionFlow);
       calc(" DepositionFlow = min(0, TransportCapacity-SedConcentration)*TransportFactor");

       //calc(" DepositionFlow = max(DepositionFlow, min(0,TransportCapacity-SedConcentration)*WaterHVolin) ");
       calc(" DepositionFlow = max(DepositionFlow, -Sedin) ");
          // Deposition amount (kg)
          // (positive = supply to the flow; negative = deposition)
          // DepositionFlow is the amount of deposition (negative!)
          // deposition cannot be larger than the amount of available sediment
          // =max because of negative values
          // if the flow is 0, than all sediment is deposited

 /* discussie met Mike Kirkby
       _spatial(REAL4, DepositionFlow);
       _spatial(REAL4, newconc);
       calc(" newconc = mif(Qin gt 0, SedConcentration * exp(-1.0*SettlingVelocity/(Qin/DX)*DXc),0) ");
       writeTimeseries(newconc,"nc");
       calc(" DepositionFlow = -max(0,SedConcentration-newconc)* WaterHVolin ");
 */
       if (SwitchNoErosionOutlet)
       {
          calc(" DepositionFlow = mif(Outlet1 eq 1, 0,DepositionFlow) ");
       }
//VJ 040224 added if outlet no deposition

       if (SwitchGrassPresent)
         calc(" DepositionFlow = mif(GrassWidth gt 0, "
              "-Sedin*GrassFraction+(1-GrassFraction)*DepositionFlow, DepositionFlow) ");
         // total deposition on grassstrips

        calc(" DetachmentFlow *= (1-StoneFraction) ");
        //VJ 031112 no erosion where stones

//VJ 080423 Snowmelt
       if (SwitchSnowmelt) //if is actually not necessary, snow cover is 0 when no snow or switched off
       {
           calc(" DetachmentFlow *= (1-SnowCover) ");
           calc(" DepositionFlow *= (1-SnowCover) ");
       }


// *****************************************************************************
// ******** Combine all fluxes for Sedin *********************************
// *****************************************************************************

//       calc(" DetachmentFlow = mif(RunoffMeanHin gt MinimumHeight, DetachmentFlow, 0) ");
//       calc(" DepositionFlow = mif(RunoffMeanHin gt MinimumHeight, DepositionFlow, -Sedin) ");
 
       calc(" Sedin += DetachmentFlow ");
       calc(" Sedin += DepositionFlow ");

       calc(" SedConcentration = mif( WaterHVolin gt 0, Sedin/WaterHVolin, 1000) ");
       calc(" SedConcentration = min(SedConcentration, 848) ");
       calc(" sedvoltemp = SedConcentration*WaterHVolin ");
       calc(" DepositionFlow += min(0, sedvoltemp-Sedin) ");
       calc(" Sedin = WaterHVolin * SedConcentration ");


    } //not no erosion


