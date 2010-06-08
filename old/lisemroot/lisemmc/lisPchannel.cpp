// *****************************************************************************
// ******** channel flow calculations ******************************************
// *****************************************************************************

   /*
   <-  chan widthUp  ->
   \  |            |  /
    \ |      chan h| /
     \|____________|/ ->side angle, 0 = rectangular
        chan width
   */


     calc(" ChannelVolin += RunoffVolinToChannel ");
        // add overland flow in channel cells

//VJ 040823: bug fix : rainfall volume added over projected cell DX and not DXc
     calc(" ChannelVolin += mif(ChannelWidthDX gt 0, (RainH*ChannelWidthUpDX*DX)/1000, 0) ");
     // add rainfall on channel, no interception

//VJ 080217 add baseflow
     if (SwitchChannelBaseflow)
     {
         _nonspatial(REAL4, checksum);
         calc(" checksum = sum(RunoffVolinToChannel)");
         if (checksum > 0)
            startbaseflowincrease = true;
         // check if is somewhere overland flow going into the channel
         //when this is true the channel baseflow is rising with a fraction until the end of the run
         startbaseflowincrease = checksum > 0;
         _spatial(REAL4, ChannelBasevolume);
         calc(" ChannelBasevolume = ChannelBaseflux * DTSEC ");

         if (startbaseflowincrease)
         {
            calc(" AddVolume += (ChannelBasevolume*ChannelBaseincrease) ");
            //addvolume declared in lisparamchannel.cpp
         }
         else
         {
            calc(" AddVolume = ChannelBasevolume ");
         }
         calc(" ChannelVolin += AddVolume ");

         calc(" TotalBaseflowVol += sum(AddVolume) ");
     }

//VJ 040823: no rainfall on buffer, already done in Prainfall
     if (SwitchBuffers)
        calc(" ChannelVolin -= mif(BufferID gt 0, (RainH*ChannelWidthUpDX*DX)/1000, 0) ");
        // ADD rainfall, no interception

     // channel height from volume can be solved by a squareroot equation:
     begin{
       _spatial(REAL4, a);
       calc(" a = mif(ChannelVolin > 0, ChannelSide*DXc/ChannelVolin, 0 ) ");

       _spatial(REAL4, b);
       calc(" b = mif(ChannelVolin > 0, ChannelWidthDX*DXc/ChannelVolin, 0 ) ");

       _spatial(REAL4, c);
       calc(" c = -1.0 ");

        calc(" ChannelWH = 1000*mif(ChannelVolin gt 0 and a gt 0,              "
             "          ( -b + sqrt(b*b-4*a*c) ) / (2*a) , 0)     ");
     }end
//VJ 030412 factor 1000 in wrong place
     calc(" ChannelWH = mif(ChannelSide eq 0, ChannelVolin/(DXc*ChannelWidthDX)*1000, ChannelWH) ");
       // ChannelWH is the water depth in the channel in mm

     calc(" ChannelWidthUpDX = mif(ChannelWidthDX gt 0, ChannelWidthDX+2*ChannelSide*ChannelWH/1000,0) ");
     // up width is width + 2*tan angle * wh
     calc(" ChannelWidthUpDX = cover(min(ChannelWidthUpDX, 0.9*DX-RoadWidthDX), 0)");
      // channelwidth cannot be larger than 90% of the grid cell
      //VJ 100206 added roadwidth here

     // **** flow rate and flux in channel
     _spatial(REAL4, cosside);
     calc(" cosside = cos(atan(ChannelSide))");
     _spatial(REAL4, ChannelPerimeter);
     calc(" ChannelPerimeter = ChannelWidthDX+2*0.001*ChannelWH/cosside ");

     _spatial(REAL4, ChannelCSArea);
//     calc(" ChannelCSArea = (ChannelWidthDX+0.5*(ChannelWidthUpDX-ChannelWidthDX)) * ChannelWH/1000");
     calc(" ChannelCSArea = (ChannelWidthDX+(ChannelWidthUpDX-ChannelWidthDX)) * ChannelWH/1000");
     //BUG NOT 0.5! but 2* 0.5
//VJ 100206 changed calculations of channel A and P
     _spatial(REAL4, ChannelV);
     calc(" ChannelV = mif(ChannelPerimeter gt 0, sqrt(ChannelGradient)/ChannelN * (ChannelCSArea/ChannelPerimeter)**(2.0/3.0), 0) ");
        // ChannelV is the channel flow rate, in meters per second
        // ChannelGradient is the slope in the channel, in percentage, 10% = 0.10
        // ChannelN is Mannings n in channel
        // ChannelWH is the water height in channel, in mm
        // ChannelWidthDX is the channel width, in meters
        // ChannelSide is the tangent of the channel edge angle,
        // ChannelSide = 0 for rectangular channel, 1 for 45 degrees

     _spatial(REAL4, ChannelBeta);
     calc(" ChannelBeta = 0.6 ");
     _spatial(REAL4, ChannelAlpha);
     calc(" ChannelAlpha = (ChannelN/sqrt(ChannelGradient)*ChannelPerimeter**(2.0/3.0))**ChannelBeta ");
     _spatial(REAL4, ChannelQin);
     calc(" ChannelQin = mif(ChannelAlpha gt 0 , (ChannelCSArea/ChannelAlpha)**(1/ChannelBeta), 0 ) ");

     if (SwitchBuffers)
        calc(" SumBufferVolume = sum(BufferVolumeInitChannel) - sum(BufferVolumeCurrentChannel) ");


// *****************************************************************************
// **** channel detachment and transport calculations **************************
// *****************************************************************************

     if (!SwitchNoErosion)
     {
         calc(" ChannelSedin += SedinToChannel ");
            // the overland flow sediment is added to the channel system

         // **** sediment concentration
         _spatial(REAL4, ChannelSedConcentration);
         calc(" ChannelSedConcentration = mif( ChannelVolin gt 0, ChannelSedin/ChannelVolin, 1000) ");
         calc(" ChannelSedConcentration = min(ChannelSedConcentration, 848) ");
          _spatial(REAL4, sedvoltemp);
         calc(" sedvoltemp = ChannelSedConcentration*ChannelVolin ");
         _spatial(REAL4, ChannelDeposition);
         calc(" ChannelDeposition = min(0, sedvoltemp-ChannelSedin) ");
         calc(" ChannelSedin = ChannelVolin * ChannelSedConcentration ");

         // **** transport capacity
         _spatial(REAL4, ChannelTransportCapacity);
         calc(" ChannelTransportCapacity = mif(ChannelGradient*ChannelV*100 gt CriticalStreamPower,   "
              "   CGovers*(ChannelGradient*ChannelV*100-CriticalStreamPower)**(DGovers),0)    ");
         calc(" ChannelTransportCapacity = 2650*min(ChannelTransportCapacity,0.32) ");
            // TransportCapacity is the volumetric transport capacity of sediment
            // in (cm3 soil)/(cm3 water)
            // ChannelV is the channel flow rate in m/s formula requires cm/s, thus factor 100
            // CriticalStreamPower is the critical unit stream power = 0.4 cm/s
            // 2650 is the particle density in kg/m3, the unit of TransportCapacity is now in kg/m3

         // **** detachment in channels

         _spatial(REAL4, TransportFactor);
         calc(" TransportFactor = SettlingVelocity*DTSEC*DXc*ChannelWidthUpDX ");

         _spatial(REAL4, ChannelDetachmentFlow);
	       if (SwitchAltErosion)
    	   {
	    	   calc(" ChannelDetachmentFlow = ChannelY*max(0, ChannelTransportCapacity-ChannelSedConcentration)*ChannelQin*DTSEC ");
	      }
    	   else
         {
	         calc(" ChannelDetachmentFlow = ChannelY * "
                 " max(0, ChannelTransportCapacity-ChannelSedConcentration)*TransportFactor ");
         }
            // ChannelDetachmentFlow is detachment by overland flow (kg)
            // ChannelY is the flow detachment efficiency coefficient (constant map)
            // ChannelSedin is the amount of sediment in kg
            // ChannelVolin is the amount of channel water, in m3
            // which can transport an amount of sediment in m3

         calc(" ChannelDetachmentFlow = min(ChannelDetachmentFlow, "
              " max(0, ChannelTransportCapacity-ChannelSedConcentration)*ChannelVolin) ");
            // detachment cannot be larger than the remaining transp.capacit

//         calc(" TransportFactor = mif(ChannelWH gt 0, (1-exp(-DTSEC*SettlingVelocity/(0.001*ChannelWH)))*ChannelVolin, ChannelVolin) ");

         calc(" ChannelDeposition = min(0, ChannelTransportCapacity-ChannelSedConcentration)*TransportFactor");
//         calc(" ChannelDeposition = max(ChannelDeposition, min(0,ChannelTransportCapacity-ChannelSedConcentration)*ChannelVolin) ");
         calc(" ChannelDeposition = max(ChannelDeposition, -ChannelSedin) ");
            // Deposition amount (kg)
            // (positive = supply to the flow; negative = deposition)
            // deposition cannot be larger than the amount of available sediment

          calc(" ChannelSedin += ChannelDetachmentFlow ");
          calc(" ChannelSedin += ChannelDeposition ");
            // ChannelSedin is the amount of sediment in kg available for transport
            // ChannelDeposition  = negative

          calc(" ChannelSedConcentration = mif(ChannelVolin gt 0, ChannelSedin/ChannelVolin, 1000) ");
          calc(" ChannelSedConcentration = min(ChannelSedConcentration, 848) ");
          calc(" sedvoltemp = ChannelSedConcentration*ChannelVolin ");
          calc(" ChannelDeposition += min(0, sedvoltemp-ChannelSedin) ");
          calc(" ChannelSedin = ChannelVolin * ChannelSedConcentration ");
          //VJ 100509 channel sed concentration check
     }// switch no erosion

// *****************************************************************************
// ******** kinematic wave calculation for CHANNEL flow ************************
// *****************************************************************************

//     _nonspatial(REAL4, SumChannelVolin);
//     calc(" SumChannelVolin = sum(ChannelVolin) ");
     // used in mass balance error

     _spatial(REAL4, ChannelQsedin);
//     calc(" ChannelQsedin = mif(ChannelVolin gt 0, ChannelSedin*ChannelQin/ChannelVolin, 0) ");
     calc(" ChannelQsedin = ChannelQin * ChannelSedConcentration ");
        // ChannelQsedin is the sediment discharge in kg/s (m3/s*kg*m3)

     _nonspatial(REAL4, SumChannelSedin);
     calc(" SumChannelSedin = sum(ChannelSedin) ");
     // used in mass balance error

//VJ 050704 added channel infil in m2/s
     _spatial(REAL4, chq);
     if (SwitchChannelInfil)
       calc(" chq = -(ChannelKsat *  ChannelPerimeter /3600000.0) ");
       //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
     else
       calc(" chq = 0 ");

     _spatial(REAL4, chqs);
     calc(" chqs = 0");//(DepositionFlow+DetachmentFlow+DetachmentSplash+DepositionSplash)/DXc ");

     _spatial(REAL4, ChannelQout);
     _spatial(REAL4, ChannelQsedout);

     kineDX(ChannelQout,ChannelQin, chq,
               ChannelQsedout,ChannelQsedin, chqs,
               BufferVolumeCurrentChannel, BufferSedVolumeCurrentChannel,
               ChannelLDD, ChannelAlpha, ChannelBeta, DTSEC, DXc,
               int(SwitchSedtrap), BufferSedBulkDensity);
//VJ 050831 CHANGED: infil map changes during run: after the kin wave
// it contains the sum of all incoming fluxes in the cell! Needed for infiltration calc
//VJ 100206 CHANGED again, q contains Qin in KineDX, only incoming fluxes

     calc(" ChannelQsedout = min(ChannelQsedout, ChannelSedin/DTSEC + chqs)");
     // no more sedoutflow than sum inflow in a cell plus what is in the cell/dt
     calc(" ChannelQOutflow = sum(mif(Outlet1 eq 1,ChannelQout, 0))");
     calc(" ChannelSedOutflow = sum(mif(Outlet1 eq 1, ChannelQsedout, 0))");
     // recalc outlet discharge, kin wave can handle more pits

     calc(" ChannelCSArea = ChannelAlpha*(ChannelQout**ChannelBeta)");
     // ChannelCSArea is A (cross section) in m2

     //moved to lisparaminit.cpp
     //_spatial(REAL4, ChannelVolout);
     calc(" ChannelVolout = ChannelCSArea*DXc ");
     // m3

     calc(" SumChannelVolout = sum(ChannelVolout) ");
     // spatial total of water in channel, m3 , after kin wave
     // needed in screen output

//VJ 050704 added channel infil in m2/s
//VJ 050831 corrected channel infil calc
     if (SwitchKinwaveInfil)
     {
//VJ 100206 changed to infil per cell = sum inflow (q) + (volume before - after) - outflow)
        calc(" InfilVolinKinWave = chq*DTSEC + ChannelVolin - ChannelVolout - ChannelQout*DTSEC ");
        calc(" InfilVol += cover(InfilVolinKinWave,0) ");
	     calc(" TotalInfilVol += sum(cover(InfilVolinKinWave,0)) ");
        // add to totals
	   }

     calc(" ChannelVolout = max(ChannelVolout, 0) ");
     calc(" ChannelVolin = ChannelVolout ");

     if (SwitchBuffers)
        calc(" SumBufferVolumeChannel = sum(BufferVolumeInitChannel) - sum(BufferVolumeCurrentChannel) ");

     if (!SwitchNoErosion)
     {
         // **** correct mass balance error SED
         _spatial(REAL4, ChannelSedout);
         calc(" ChannelSedout = max(0, chqs*DTSEC + ChannelSedin - ChannelQsedout*DTSEC) ");

/*
         calc(" ChannelSedout = mif(ChannelQout gt 0, ChannelQsedout/ChannelQout*ChannelVolout, 0) ");
         // sediment discharge
         if (SwitchCorrectMassSED)
         {
             calc(" SumChannelSedout = sum(ChannelSedout) ");
             if (SumChannelSedout > 0)
                calc(" ChannelSedout += (SumChannelSedin-SumChannelSedout-ChannelSedOutflow*DTSEC)*ChannelSedout/SumChannelSedout ");
         }
*/
         calc(" ChannelSedin = mif(ChannelWidthUpDX gt 0,ChannelSedout, 0) ");
         calc(" ChannelSedConcentration = mif( ChannelVolout gt 0, ChannelSedin/ChannelVolout, 1000) ");
         calc(" ChannelSedConcentration = min(ChannelSedConcentration, 848) ");
         calc(" sedvoltemp = ChannelSedConcentration*ChannelVolout ");
         calc(" ChannelDeposition += min(0, sedvoltemp-ChannelSedin) ");
         calc(" ChannelSedin = ChannelVolout * ChannelSedConcentration ");

         calc(" SumChannelSedout = sum(ChannelSedin) ");

         if (SwitchBuffers)
            calc(" SumBufferSedVolumeChannel = sum(BufferSedVolumeCurrentChannel) ");
     }// switch no erosion



