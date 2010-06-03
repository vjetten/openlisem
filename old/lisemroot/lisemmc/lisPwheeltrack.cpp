// *****************************************************************************
// ******** wheeltrack flow calculations ***************************************
// ******** last changes: VJ 26/02/2001 ****************************************
// *****************************************************************************

//VJ 030701 corrected all DX to DXc
//VJ 030702:
//##################################################################
//!!!!!!!!!!!!! ALL THIS IS FOR ONE TRACK WIH A COMPLEX PERIMETER!
//##################################################################

       _spatial(REAL4, WheelVolin);
       calc(" WheelVolin = WHWheelTrack*WheelWidthDX*DXc/1000 ");

         // volume in wheeltracks (m3)
       calc(" WheelVolin += cover(RunoffVolinToWheel,0) ");
       // divide by wheelnumber because calc for one track

         // add to flow from gridcells
//       calc(" WHWheelTrack = mif(WheelWidthDX gt 0, WheelVolin/(DXc*WheelWidthDX/WheelNumber) * 1000, 0) ");
       calc(" WHWheelTrack = mif(WheelWidthDX gt 0, WheelVolin/(DXc*WheelWidthDX) * 1000, 0) ");
         // recalc wheeltrack water height (mm)

       // **** flow rate and flux in channel
/*
       _spatial(REAL4, WheelV);
       _spatial(REAL4, WheelP);
       //calc(" WheelP = mif(WheelNumber gt 0, WheelNumber*(WheelWidthDX/WheelNumber + 2*WHWheelTrack/1000),0) ");
       calc(" WheelP = WheelWidthDX + 2*WheelNumber*WHWheelTrack/1000 ");
       // perimeter of wheelnumber tracks
       calc(" WheelV = mif(WheelNumber gt 0 and WheelP gt 0, ((WheelWidthDX/WheelNumber*WHWheelTrack/1000)/WheelP)**(2/3)"
            "           * sqrt(WheelGradient)/WheelN, 0) ");

       // speed in ONE wheel track
       _spatial(REAL4, WheelBeta);
       calc(" WheelBeta = 0.6 ");
       _spatial(REAL4, WheelAlpha);
       calc(" WheelAlpha = (WheelN/sqrt(WheelGradient)*WheelP**(2/3))**WheelBeta ");
         // Alpha for one of the wheeltracks
       _spatial(REAL4, WheelQin);
       calc(" WheelQin = mif(WheelAlpha le 0 , 0, "
            " (WheelWidthDX/WheelNumber*WHWheelTrack/1000)/WheelAlpha)**(1/WheelBeta) ");
       // Qin is for ONE wheeltrack

       // Ingrid: Voor kinematic wave totale Qin gebruiken en ook alpha aanpassen voor totaal!!
       //  calc(" WheelQin = cover(WheelNumber*WheelQin,0) ");
       //  calc(" WheelAlpha = cover(WheelNumber**(1-WheelBeta) * WheelAlpha,0) ");
       //  discharge as sum of wheelnumber discharges
*/

       _spatial(REAL4, WheelV);
       _spatial(REAL4, WheelP);
       calc(" WheelP = WheelWidthDX + 2*WheelNumber*WHWheelTrack/1000 ");
       calc(" WheelV = mif(WheelP gt 0, ((WheelWidthDX*WHWheelTrack/1000)/WheelP)**(2/3)"
            "           * sqrt(WheelGradient)/WheelN, 0) ");
       _spatial(REAL4, WheelBeta);
       calc(" WheelBeta = 0.6 ");
       _spatial(REAL4, WheelAlpha);
       calc(" WheelAlpha = (WheelN/sqrt(WheelGradient)*WheelP**(2/3))**WheelBeta ");
         // Alpha for one of the wheeltracks
       _spatial(REAL4, WheelQin);
       calc(" WheelQin = mif(WheelAlpha le 0 , 0, "
            " (WheelWidthDX*WHWheelTrack/1000)/WheelAlpha)**(1/WheelBeta) ");
       // Qin is for ONE wheeltrack


// *****************************************************************************
// **** Wheel detachment and transport calculations **************************
// *****************************************************************************

       if (!SwitchNoErosion)
       {
//           calc(" WheelSedin += SedinToWheel/WheelNumber	 ");
           calc(" WheelSedin += SedinToWheel	 ");
              // the overland flow sediment is added to the Wheel system
				  //ALL THIS IS FOR ONE WHEEL TRACK!

           _spatial(REAL4, WheelTransportCapacity);
           calc(" WheelTransportCapacity = mif(WheelGradient*WheelV*100 gt CriticalStreamPower,   "
                "   CGovers*(WheelGradient*WheelV*100-CriticalStreamPower)**(DGovers),0)    ");
           calc(" WheelTransportCapacity = 2650*min(WheelTransportCapacity,0.32) ");

           // **** sediment concentration

           _spatial(REAL4, WheelSedConcentration);
           calc(" WheelSedConcentration = mif(WheelVolin gt 0, WheelSedin/WheelVolin, 0) ");

           // **** detachment in Wheels
           _spatial(REAL4, TCAux);
           calc(" TCAux = SettlingVelocity*DTSEC*DXc*WheelWidthDX ");

           _spatial(REAL4, WheelDetachmentFlow);
	       if (SwitchAltErosion)
    	   {
	    	   calc(" WheelDetachmentFlow = WheelY*max(0, WheelTransportCapacity-WheelSedConcentration)*WheelQin*DTSEC ");
	       }
    	   else
           {
        	   calc(" WheelDetachmentFlow = WheelY*max(0, WheelTransportCapacity-WheelSedConcentration)*TCAux ");
           }
                // detachment in kg
           calc(" WheelDetachmentFlow = min(WheelDetachmentFlow, max(0, WheelTransportCapacity-WheelSedConcentration)*WheelVolin) ");


          calc(" TCAux = mif(WHWheelTrack gt 0, (1-exp(-DTSEC*SettlingVelocity/(0.001*WHWheelTrack)))*WheelVolin, 1) ");
          _spatial(REAL4, WheelDeposition);
          calc(" WheelDeposition = min(0, WheelTransportCapacity-WheelSedConcentration)*TCAux");
          calc(" WheelDeposition = max(WheelDeposition, min(0,WheelTransportCapacity-WheelSedConcentration)*WheelVolin) ");
                // deposition in kg

          calc(" WheelDetachmentFlow = mif(WHWheelTrack gt MinimumHeight, WheelDetachmentFlow, 0)");
          calc(" WheelDeposition = mif(WHWheelTrack gt MinimumHeight, WheelDeposition, -WheelSedin)");

//ALL THIS IS FOR ONE WHEEL TRACK!

          calc(" WheelSedin += WheelDetachmentFlow ");
          calc(" WheelSedin += WheelDeposition ");
               // wheelsedin is amount of sediment in suspension (kg)

       }// switch no erosion

// *****************************************************************************
// ******** kinematic wave calculation for Wheel flow ************************
// *****************************************************************************

       _spatial(REAL4, WheelQsedin);
       calc(" WheelQsedin = mif(WheelVolin gt 0, WheelSedin*WheelQin/WheelVolin, 0) ");
          // WheelQsedin is the sediment discharge in kg/s (kg*m3/s/m3)

       _nonspatial(REAL4, SumWheelVolin);
       calc(" SumWheelVolin = sum(cover(WheelVolin,0)) ");
       // used in mass balance error

       _nonspatial(REAL4, SumWheelSedin);
       calc(" SumWheelSedin = sum(cover(WheelSedin,0)) ");
       // used in mass balance error

       _spatial(REAL4, WheelQout);
       _spatial(REAL4, WheelQsedout);

      _spatial(REAL4, infilWT);
      if (SwitchKinwaveInfil)
     //VJ 050704 changed DXc to DX
     //VJ 051005 changed back to DXc !!
        calc(" infilWT = InfilSurplusWT/(1000*DTSEC)*DX ");
      else
        calc(" infilWT = 0 ");
     _spatial(REAL4, stempWT);
     calc(" stempWT = 0");

     calc(" InfilVolinKinWave = max(-InfilSurplusWT/(1000)*DX*DXc,0) ");
     //VJ 050831 volume available for infil in m3


     //VJ 040823 include buffers, needed to correct infil in kin wave
     if (SwitchBuffers)
        calc(" buffersize = sum(BufferVolumeCurrent) ");

     kineDX(WheelQout,WheelQin, infilWT,
            WheelQsedout,WheelQsedin, stempWT,
            BufferVolumeCurrent, BufferSedVolumeCurrent,
            WheelLDD, WheelAlpha, WheelBeta, DTSEC, DXc, int(SwitchSedtrap));
//VJ 050831 CHANGED: infil changes: it is the sum of all fluxes in the cell. Needed for infiltration calc

       calc("WheelQOutflow = sum(mif(Outlet1 eq 1,WheelQout,0))");
       calc("WheelSedOutflow = sum(mif(Outlet1 eq 1,WheelQsedout,0))");

       // kin wave can handle more pits now, so ensure output catchment is outlet1

       _spatial(REAL4, WheelCSArea);
       calc(" WheelCSArea = WheelAlpha*(WheelQout**WheelBeta) ");
       // WheelCSArea is A (total cross section of two tracks) in m2

       _spatial(REAL4, WheelVolout);
       calc(" WheelVolout = WheelCSArea*DXc ");
       // m3

       calc(" SumWheelVolout = sum(cover(WheelVolout,0)) ");
       // spatial total of water in Wheel, m3 , after kin wave

       calc(" pitsout = cover(mif(WheelLDD eq 5, WheelQout),0)");

       if (SwitchCorrectMass)
       {
           if (SumWheelVolout > 0)
           calc(" WheelVolout += WheelVolout*"
                "(SumWheelVolin-SumWheelVolout-sum(pitsout)*DTSEC)/SumWheelVolout ");
           // divide error caused by pits over network
           calc(" SumWheelVolout = sum(WheelVolout) ");
       }

       calc(" WheelVolout = max(WheelVolout, 0) ");
       calc(" WheelVolin = WheelVolout ");
       _spatial(REAL4, WHWheelin);
       calc("WHWheelin = WHWheelTrack ");
//       calc(" WHWheelTrack = mif(WheelWidthDX gt 0, WheelVolout*1000/(DXc*WheelWidthDX/WheelNumber), 0)");
       calc(" WHWheelTrack = mif(WheelWidthDX gt 0, WheelVolout*1000/(DXc*WheelWidthDX), 0)");
            // recalc WH from volume

     //VJ 050831 made changes here
     if (SwitchKinwaveInfil){
        //calc(" TotalInfilVol += sum(min(infil*DTSEC+WaterHVolin, InfilVolinKinWave)) ");
        //nonspatial infil adjustment in m3
        calc(" InfilVol += cover(min(infil*DTSEC+WaterHVolin, InfilVolinKinWave),0) ");
        //spatial infil adjustment in m3

        calc(" TotalInfilVol += max(0,SumWheelVolin - SumWheelVolout - sum(pitsout*DTSEC)) ");
       // old stuff
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
         // **** correct mass balance error SED
         _spatial(REAL4, WheelSedout);
         calc(" WheelSedout = mif(WheelQout gt 0, WheelQsedout/WheelQout*WheelVolout, 0) ");

         // sediment discharge

         calc(" pitsoutsed = cover(mif(WheelLDD eq 5, WheelQsedout),0)");
         if (SwitchCorrectMassSED)
         {
             calc(" SumWheelSedout = sum(WheelSedout) ");
             if (SumWheelSedout > 0)
                calc(" WheelSedout += (SumWheelSedin-SumWheelSedout-pitsoutsed*DTSEC)*WheelSedout/SumWheelSedout ");
         }

         calc(" WheelSedout = max(WheelSedout, 0) ");
         calc(" WheelSedin = mif(WheelWidthDX gt 0,WheelSedout, 0) ");

         calc(" SumWheelSedout = sum(WheelSedin) ");

             // recalc sediment outlet
     }// switch no erosion


// *****************************************************************************
// ******** Overflow wheeltracks ***********************************************
// *****************************************************************************


/* trial doesn't work
/*  THROWING Q ON THE LDD NETWORK DOESN'T WORK
  //     _spatial(REAL4, Diffvol);
  //     calc(" Diffvol = pitsout*DTSEC ");
  //     calc(" WheelVolout -= mif(WheelLDD eq 5, Diffvol, 0) ");
//       calc(" WheelVolin = WheelVolout ");
    //   calc(" SumWheelVolout = sum(cover(WheelVolout,0)) ");
  //     calc(" WHWheelTrack = mif(WheelWidthDX gt 0,WheelVolout/(WheelWidthDX*DXc)*1000,0)");

    //   calc(" WheelSedin -= mif(WheelLDD eq 5, pitsoutsed*DTSEC, 0) ");
//       calc(" SumWheelSedout = sum(WheelSedin) ");
*/

//NOTE : FractionToWheeltrack determines mass balance error important!!!

       _spatial(REAL4, DiffWH);
//       calc(" DiffWH = mif(WheelLDD eq 5, max(0, WHWheelTrack-WheelDepth),0)");
       calc(" DiffWH = cover(max(WHWheelTrack-WheelDepth,0),0)");
         // calc overflow water difference

       if (!SwitchNoErosion)
       {
          calc(" WheelSedConcentration = mif(WheelVolin gt 0, WheelSedin/WheelVolin, 0) ");
         _spatial(REAL4, DiffSed);
          calc(" DiffSed = cover(DiffWH/1000*WheelWidthDX*DXc * WheelSedConcentration, 0) ");
          calc(" Sedin += DiffSed ");
          calc(" WheelSedin -= DiffSed ");
          calc(" SumWheelSedout = sum(WheelSedin) ");
       }

       _spatial(REAL4, overflowfraction);
       calc(" overflowfraction = WheelWidthDX/(DX-RoadWidthDX-ChannelWidthUpDX) ");
//       calc(" overflowfraction = WheelWidthDX/(DX-RoadWidthDX-ChannelWidthUpDX-WheelWidthDX) ");
       //fraction of channel in cell

       calc(" DiffWH = cover(DiffWH*overflowfraction,0) ");

       calc(" WHWheelTrack = min(WheelDepth+DiffWH, WHWheelTrack)");
//VJ this works

       calc(" WH += DiffWH ");

       calc(" WheelVolout = WHWheelTrack/1000*WheelWidthDX*DXc ");
       calc(" WheelVolin = WheelVolout ");
       calc(" SumWheelVolout = sum(cover(WheelVolout,0)) ");





