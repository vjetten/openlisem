// *****************************************************************************
// ******** kinematic wave calculation for overland flow ***********************
// *****************************************************************************
//VJ 030701 changed definition, wheelwidth and giullywidth are 0 anyway when not chosen
// avoid confusion
//??WAAROM HERHALING HIER?
      calc(" StoneWidthDX = max(DX-ChannelWidthUpDX-GullyWidthDX-RoadWidthDX-WheelWidthDX,0)*StoneFraction ");
      calc(" SoilWidthDX = max(DX-ChannelWidthUpDX-GullyWidthDX-RoadWidthDX-WheelWidthDX-StoneWidthDX,0) ");
      calc(" FlowWidth = RoadWidthDX+(SoilWidthDX+StoneWidthDX)*PondAreaFract ");

         // adjust because of possible changes in gully width and channel width

      // **** recalc Alpha and Qin after tochannel, towheel or togully
      calc(" Perimeter = FlowWidth+2*RunoffMeanHin/1000 ");
      calc(" HydrRadius = mif(Perimeter gt 0,(FlowWidth*RunoffMeanHin/1000)/Perimeter, 0) ");

// just some tryouts      
//      calc(" Perimeter = FlowWidth ");
//      calc(" HydrRadius = RunoffMeanHin/1000 ");
//      calc(" HydrRadius = mif(PondAreaFract gt 0, 0.001*(-0.0084*RR*10 + 0.7966)/PondAreaFract "
//                "* (RunoffMeanHin**(0.0027*RR*10 + 1.0124)),0) ");


      calc(" V = HydrRadius**(2.0/3.0) * sqrt(Gradient) / N ");
//VJ 050303 V was not recalculated to questionable if streampower is good in channel cells

      calc(" Alpha = ( N/sqrt(Gradient) *(Perimeter**(2.0/3.0)) )**Beta ");
      calc(" Qin = mif(Alpha gt 0 , ((FlowWidth * RunoffMeanHin/1000)/Alpha)**(1/Beta) , 0) ");
  write(Qin,"qin.map");
  write(Alpha,"alpha.map");

     _spatial(REAL4, DeprSedin);
     calc(" DeprSedin = mif(WH gt 0, Sedin*SurfStorH/WH,0) ");
      // sediment in depr storage

//     _spatial(REAL4, SwitchTill);
//     calc("SwitchTill = mif(WH gt 1000, 1, 0) ");

     calc(" WH = SurfStorH ");
       // the waterheight is reduced to the retention
       // the new detention RunoffMeanHout will be calculated, based on the new
       // Qout, and will be added after that to WH

     calc(" SumDeprStoreVol = sum(SurfStorH*(SoilWidthDX+StoneWidthDX)*DXc/1000.0) ");
        // spatial total depression storage for mass balance, for one DT

     calc(" SumWaterHVolin = sum(WaterHVolin) ");
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

 //    _spatial(REAL4, QoutTill);
 //   _spatial(REAL4, QsedoutTill);

     _spatial(REAL4, Qout);
     _spatial(REAL4, Qsedout);
     _spatial(REAL4, infil);
        //infiltration surplus in m2/s
     if (SwitchKinwaveInfil)
     //VJ 050704 changed DXc to DX
     //VJ 051005 changed back to DXc!
       calc(" infil = min((InfilSurplus)/(1000*DTSEC)*DX,0) ");
       // unit now m2/s multiply with DXc to get m3/s
     else
       calc(" infil = 0 ");

     calc(" InfilVolinKinWave = max(-InfilSurplus/(1000)*DX*DXc,0) ");
     //VJ 050831 volume available for infil in m3

     _spatial(REAL4, stemp);
     calc(" stemp = 0");
     
//VJ 040823 include buffers, needed to correct infil in kin wave
     if (SwitchBuffers){
        calc(" buffersize = sum(BufferVolumeCurrent) ");
        calc(" sedbuffersize = sum(BufferSedVolumeCurrent) ");
     }

     kineDX(Qout, Qin, infil,
            Qsedout, Qsedin, stemp,
            BufferVolumeCurrent, BufferSedVolumeCurrent,
            LDD, Alpha, Beta, DTSEC, DXc);
//VJ 050831 CHANGED: infil changes in the kin wave:
//          it is the sum of all fluxes in the cell in m3/s. Needed for infiltration calc

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

     calc(" DXcorr = max(DXcorr, DX*0.1) ");
        // to avoid a dxcorr of zero

    // **** calc new waterheights by adding new detention to retention

     _spatial(REAL4, RunoffHout);
     calc(" RunoffHout = 1000 * RunoffAout/DXcorr ");
     // RunoffHout is the average water depth in the pixel, in mm

     calc(" WH += RunoffHout");
     calc(" WHRoad = RunoffHout");
       // add new WH to surfaces
     _spatial(REAL4, WaterHVolout);
     calc(" WaterHVolout = (RoadWidthDX*WHRoad+(SoilWidthDX+StoneWidthDX)*WH)*DXc/1000 ");
       // new volume
     calc(" SumWaterHVolout = sum(WaterHVolout) ");
     // spatial total amount of water on surface in m3 after kin wave
     // used in water mass balance error, surface depression storage included

     //VJ 050831 made changes here
     if (SwitchKinwaveInfil){     //always true
//        calc(" TotalInfilVol += sum(min(infil*DTSEC+WaterHVolin, InfilVolinKinWave)) ");
        //nonspatial infil adjustment in m3
        calc(" InfilVol += min(infil*DTSEC+WaterHVolout, InfilVolinKinWave) ");
        //spatial infil adjustment in m3

       // old stuff
       calc(" TotalInfilVol += max(0,SumWaterHVolin - SumWaterHVolout - QOutflow*DTSEC) ");
       // add volume infiltrated in kin wave to suminfiltration (NOTE, this includes possible errors)
     }

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

           _spatial(REAL4, Sedout);

//           calc(" Sedout = mif(Qout gt 0, Qsedout/Qout * WaterHVolout, 0) ");
              // sediment load is volume * concentration (=Qs/Q)
//of
           if (SwitchBuffers){
              calc(" SumBufferSedVolume = sum(BufferSedVolumeCurrent) ");
              calc(" sedbuffersize = SumBufferSedVolume - sedbuffersize");
           }


           calc(" Sedout = mif(Qout gt 0, Qsedout/Qout * RunoffAout*DXc, 0) ");
           calc(" Sedout += DeprSedin ");

           if (SwitchCorrectMassSED) // default true !
           {
               calc(" SumSedout = sum(Sedout) ");
               if (SumSedout < -1e-6 || SumSedout > 1e-6)
                 calc(" Sedout += Sedout*(SumSedin-SumSedout-SedOutflow*DTSEC)/(SumSedout) ");
           }

           calc(" Sedin = Sedout");
       //    calc(" DepositionFlow += mif( Sedin lt -0.01, Sedin, 0)");

           calc(" Sedin = max(Sedin, 0) ");
           calc(" SumSedout = sum(Sedin) ");
           // needed in mass balance error

     }//no erosion


