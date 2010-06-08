// *****************************************************************************
// ******** kinematic wave calculation for overland flow ***********************
// *****************************************************************************
//VJ 030701 changed definition, wheelwidth and giullywidth are 0 anyway when not chosen
// avoid confusion
//      calc(" StoneWidthDX = max(DX-ChannelWidthUpDX-GullyWidthDX-RoadWidthDX-WheelWidthDX,0)*StoneFraction ");
//      calc(" SoilWidthDX = max(DX-ChannelWidthUpDX-GullyWidthDX-RoadWidthDX-WheelWidthDX-StoneWidthDX,0) ");
//      calc(" SoilWidthDX = max(DX-ChannelWidthUpDX-GullyWidthDX-RoadWidthDX-WheelWidthDX,0) ");
//      calc(" FlowWidth = RoadWidthDX+(SoilWidthDX)*PondAreaFract ");
       // adjust because of possible changes in gully width and channel width

      // **** recalc Alpha and Qin after tochannel, towheel or togully
      calc(" Perimeter = FlowWidth+2*RunoffMeanHin/1000 ");
      calc(" HydrRadius = mif(Perimeter gt 0,(FlowWidth*RunoffMeanHin/1000)/Perimeter, 0) ");

      calc(" V = HydrRadius**(2.0/3.0) * sqrt(Gradient) / N ");
//VJ 050303 V was not recalculated to questionable if streampower is good in channel cells

      calc(" Alpha = ( N/sqrt(Gradient) *(Perimeter**(2.0/3.0)) )**Beta ");
      calc(" Qin = mif(Alpha gt 0 , ((FlowWidth * RunoffMeanHin/1000)/Alpha)**(1/Beta) , 0) ");

 //    _spatial(REAL4, DeprSedin);
 //    calc(" DeprSedin = mif(WH gt 0, Sedin*SurfStorH/WH,0) ");
 //    calc(" Sedin -= DeprSedin ");
      // sediment in depr storage

     calc(" WH = SurfStorH ");
       // the waterheight is reduced to the retention
       // the new detention RunoffMeanHout will be calculated, based on the new
       // Qout, and will be added after that to WH

//     calc(" SumDeprStoreVol = sum(SurfStorH*(SoilWidthDX+StoneWidthDX)*DXc/1000.0) ");
     calc(" SumDeprStoreVol = sum(SurfStorH*SoilWidthDX*DXc/1000.0) ");
        // spatial total depression storage for screen output, for one DT

     calc(" WaterHVolin = (FlowWidth*RunoffMeanHin + SurfStorH*SoilWidthDX)*DXc/1000.0 ");
//     calc(" SumWaterHVolin = sum(WaterHVolin) ");
        // needed for infil in kin wave

     calc(" WHRoad = 0 ");
       // all water on roads contributes to flow

     _spatial(REAL4, Qsedin);
     calc(" Qsedin = mif(WaterHVolin gt 0, Qin * Sedin/WaterHVolin, 0) ");
     // Qsedin is the sediment discharge in kg/s (kg*m3/s/m3)
     // calculated as concentraion * waterdischarge )

     _nonspatial(REAL4, SumSedin);
     calc(" SumSedin = sum(Sedin) ");
     //calc total for sed mass balance error correction, kg

     _spatial(REAL4, Qout);
     _spatial(REAL4, Qsedout);
     _spatial(REAL4, q);
        //infiltration surplus in m2/s
     if (SwitchKinwaveInfil)
     //VJ 050704 changed DXc to DX
     //VJ 051005 changed back to DXc!
       calc(" q = min((InfilSurplus)/(1000*DTSEC)*DX,0) ");
       // unit now m2/s multiply with DXc to get m3/s
     else
       calc(" q = 0 ");

     _spatial(REAL4, qs);
     calc(" qs = 0");
     
//VJ 040823 include buffers, needed to correct infil in kin wave
     if (SwitchBuffers){
        calc(" buffersize = sum(BufferVolumeCurrent) ");
        calc(" sedbuffersize = sum(BufferSedVolumeCurrent) ");
     }

     kineDX(Qout, Qin, q,
            Qsedout, Qsedin, qs,
            BufferVolumeCurrent, BufferSedVolumeCurrent,
            LDD, Alpha, Beta, DTSEC, DXc,
            int(SwitchSedtrap), BufferSedBulkDensity);
            //NOTE: q and qs now contain all incoming water and sediment of a cell (Qin and Sin) !!!

//VJ 050831 CHANGED: infil changes in the kin wave:
//          it is the sum of all fluxes in the cell in m3/s. Needed for infiltration calc

     calc(" Qsedout = min(Qsedout, qs + Sedin/DTSEC) ");
     // cannot have more Qsedout than sediment present in cell

     calc(" QOutflow = sum(mif(Outlet1 eq 1, Qout, 0))");
     calc(" SedOutflow = sum(mif(Outlet1 eq 1, Qsedout, 0))");
     // kin wave can handle more pits now so make sure QOutflow
     // is Qout at outlet 1

     _spatial(REAL4, RunoffAout);
     calc(" RunoffAout = Alpha*(Qout**Beta)");
     // RunoffAout is the cross section of the flow after the KW

     _spatial(REAL4, DXcorr);
     calc(" DXcorr = DX ");

     if (SwitchWheelPresent)
        calc(" DXcorr -= cover(WheelWidthDX,0) ");

     if (SwitchIncludeChannel)
        calc(" DXcorr -= cover(ChannelWidthUpDX,0) ");

     if (SwitchGullies)
        calc(" DXcorr -= cover(GullyWidthDX,0) ");

   //  calc(" DXcorr = max(DXcorr, DX*0.1) ");
        // to avoid a dxcorr of zero

    // **** calc new waterheights by adding new detention to retention

     _spatial(REAL4, RunoffHout);
     calc(" RunoffHout = 1000 * RunoffAout/DXcorr ");
     // RunoffHout is the average water depth in the pixel, in mm

     calc(" WH += RunoffHout");  // depr storage is added here
     calc(" WHRoad = RunoffHout");
       // add new WH to surfaces
     _spatial(REAL4, WaterHVolout);
//     calc(" WaterHVolout = (RoadWidthDX*WHRoad+(SoilWidthDX+StoneWidthDX)*WH)*DXc/1000 ");
     calc(" WaterHVolout = (RoadWidthDX*WHRoad + SoilWidthDX*WH)*DXc/1000 ");
       // new volume
     calc(" SumWaterHVolout = sum(WaterHVolout) ");
     // spatial total amount of water on surface in m3 after kin wave
     // used in water mass balance error, surface depression storage included

     if (SwitchKinwaveInfil){     //always true
        calc(" InfilKinWave =(q*DTSEC + WaterHVolin - WaterHVolout - Qout*DTSEC) ");
        // q is all incoming Q in a cell, changed in kin wave
        calc(" InfilVol += InfilKinWave ");
        // spatial cumulative infil for output
        calc(" TotalInfilVol += sum(InfilKinWave)");
        // add volume infiltrated in kin wave to suminfiltration (NOTE, this includes possible errors)
     }
     //VJ 100509 changed to decrease mbe



//VJ 040823 include buffers, correct infil in kin wave which has the buffer effect
// (all losses end up in infil)
     if (SwitchBuffers){
        calc(" buffersize =  buffersize - sum(BufferVolumeCurrent) ");
        calc(" SumBufferVolume = sum(BufferVolumeInit) - sum(BufferVolumeCurrent) ");
        calc(" TotalInfilVol -= buffersize ");
     }

     if (!SwitchNoErosion)
     {
           // **** Calculation New Sediment amount Sedout

           if (SwitchBuffers){
              calc(" SumBufferSedVolume = sum(BufferSedVolumeCurrent) ");
              calc(" sedbuffersize = SumBufferSedVolume - sedbuffersize");
           }

           _spatial(REAL4, Sedout);
           calc(" Sedout = max(0, qs*DTSEC + Sedin - Qsedout*DTSEC) ");
           //sediment = all incoming + what was there - all outgoing
//           calc(" Sedout += DeprSedin ");
/*
           if (SwitchCorrectMassSED) // default true !
           {
               calc(" SumSedout = sum(Sedout) ");
               if (SumSedout < -1e-6 || SumSedout > 1e-6)
                 calc(" Sedout += Sedout*(SumSedin-SumSedout-SedOutflow*DTSEC)/(SumSedout) ");
           }
 */
           calc(" Sedin = Sedout");
           _spatial(REAL4, sedvoltemp);
           calc(" SedConcentration = mif( WaterHVolout gt 0, Sedin/WaterHVolout, 1000) ");
           calc(" SedConcentration = min(SedConcentration, 848) ");
           calc(" sedvoltemp = SedConcentration*WaterHVolout ");
           calc(" DepositionFlow += min(0, sedvoltemp-Sedin) ");
           calc(" Sedin = WaterHVolout * SedConcentration ");
           // VJ 100501 added conc check

           calc(" SumSedout = sum(Sedin) ");
               // needed in overall mass balance error

     }//no erosion


