// *****************************************************************************
// ******** channel flow calculations ******************************************
// *****************************************************************************

       calc(" ChannelVolin += RunoffVolinToChannel ");
         // add overland flow in channel cells
         
//VJ 040823: bug fix : rainfall volume added over projected cell DX and nor DXc
       calc(" ChannelVolin += mif(ChannelWidthDX gt 0, (RainH*ChannelWidthUpDX*DX)/1000, 0) ");
         // add rainfall, no interception

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

       // can be solved by a squareroot equation:
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
       calc(" ChannelWidthUpDX = cover(min(ChannelWidthUpDX, 0.9*DX), 0)");

       _nonspatial(REAL4,MAXWIDTH);
       calc(" MAXWIDTH = maxcell(ChannelWidthUpDX+RoadWidthDX) ");
       if (MAXWIDTH > DX)
       {
          LisemError("TOTAL ChannelWidthDX + RoadWidthDX exceeds cell width\
                        error in CHANWIDT.MAP or ROADWIDT.MAP");
       }

       // **** flow rate and flux in channel
       _spatial(REAL4, ChannelPerimeter);
       calc(" ChannelPerimeter = ChannelWidthDX+2*ChannelWH/1000*sqrt(ChannelSide*ChannelSide+1) ");
       _spatial(REAL4, ChannelCSArea);
       calc(" ChannelCSArea = ChannelWidthDX*ChannelWH/1000+ChannelSide*sqr(ChannelWH/1000)");
       _spatial(REAL4, ChannelV);
       calc(" ChannelV = mif(ChannelPerimeter gt 0, sqrt(ChannelGradient)/ChannelN * (ChannelCSArea/ChannelPerimeter)**(2.0/3.0), 0) ");

       _spatial(REAL4, ChannelBeta);
       calc(" ChannelBeta = 0.6 ");
       _spatial(REAL4, ChannelAlpha);
       calc(" ChannelAlpha = (ChannelN/sqrt(ChannelGradient)*ChannelPerimeter**(2.0/3.0))**ChannelBeta ");
       _spatial(REAL4, ChannelQin);
       calc(" ChannelQin = mif(ChannelAlpha gt 0 , (ChannelCSArea/ChannelAlpha)**(1/ChannelBeta), 0 ) ");

// *****************************************************************************
// **** channel detachment and transport calculations **************************
// *****************************************************************************

       if (!SwitchNoErosion)
       {
           mcalc(" ChannelSedinMu# += SedinToChannel# ", 6);

           // **** sediment concentration *****

           _spatial(REAL4, ChannelSedConc0);
           _spatial(REAL4, ChannelSedConc1);
           _spatial(REAL4, ChannelSedConc2);
           _spatial(REAL4, ChannelSedConc3);
           _spatial(REAL4, ChannelSedConc4);
           _spatial(REAL4, ChannelSedConc5);
           mcalc(" ChannelSedConc# = mif(ChannelVolin gt 0,ChannelSedinMu#/ChannelVolin,0) ",6);

           // **** transport capacity detachment *****

           _spatial(REAL4, ChannelTC);
           calc(" ChannelTC = CGovers*max(0,(ChannelGradient*ChannelV*100-CriticalStreamPower)**DGovers,0) ");
           calc(" ChannelTC = TCCAL * 2650*min(ChannelTC,0.32) ");

           _spatial(REAL4, ChannelTCSoil0);
           _spatial(REAL4, ChannelTCSoil1);
           _spatial(REAL4, ChannelTCSoil2);
           _spatial(REAL4, ChannelTCSoil3);
           _spatial(REAL4, ChannelTCSoil4);
           _spatial(REAL4, ChannelTCSoil5);
           mcalc(" ChannelTCSoil# = ChannelTC * FractionMu# ",6);

           // **** detachment in channels *******

           _spatial(REAL4, ChannelDF0);
           _spatial(REAL4, ChannelDF1);
           _spatial(REAL4, ChannelDF2);
           _spatial(REAL4, ChannelDF3);
           _spatial(REAL4, ChannelDF4);
           _spatial(REAL4, ChannelDF5);
//VJ 031211 fixed error, changed Y to ChannelY
           mcalc(" ChannelDF# = ChannelY*max(0,ChannelTCSoil#-ChannelSedConc#)",6);
           mcalc(" ChannelDF# *= mif(ChannelWH gt 0, (1-exp(-SV#/(0.001*ChannelWH)*DTSEC))*ChannelVolin, 0) ",6);
           mcalc(" ChannelDF# = min(ChannelDF#,max(0,ChannelTCSoil#-ChannelSedConc#)*Channelvolin",6);

           // **** NEW sediment concentration *****

           mcalc(" ChannelSedinMu# += ChannelDF# ",6);

           mcalc(" ChannelSedConc# = mif(ChannelVolin gt 0,ChannelSedinMu#/ChannelVolin,0) ",6);

           // **** transport capacity deposition *****

           findD50(D50susp,ChannelSedinMu0,ChannelSedinMu1,ChannelSedinMu2,
                   ChannelSedinMu3,ChannelSedinMu4,ChannelSedinMu5,*mu_cl);
      //     _spatial(REAL4, CGsusp);
      //     _spatial(REAL4, DGsusp);
           calc(" CGsusp = ((D50susp+5)/0.32)**-0.6 ");
           calc(" DGsusp = ((D50susp+5)/300)**0.25 ");

           _spatial(REAL4, ChannelTCsusp);
           calc(" ChannelTCsusp = mif(StreamPower gt CriticalStreamPower, "
                " CGsusp*(StreamPower - CriticalStreamPower)**DGsusp, 0) ");
           calc(" ChannelTCsusp = 2650*min(ChannelTCsusp,0.32) ");
           calc(" ChannelTCsusp = mif(D50susp == 0, 0, ChannelTCsusp)");

           _spatial(REAL4, ChannelTCSusp0);
           _spatial(REAL4, ChannelTCSusp1);
           _spatial(REAL4, ChannelTCSusp2);
           _spatial(REAL4, ChannelTCSusp3);
           _spatial(REAL4, ChannelTCSusp4);
           _spatial(REAL4, ChannelTCSusp5);
           
           _spatial(REAL4, ChannelSedinMuSum);
           calc(" ChannelSedinMuSum = ChannelSedinMu0");
           calc(" ChannelSedinMuSum += ChannelSedinMu1");
           calc(" ChannelSedinMuSum += ChannelSedinMu2");
           calc(" ChannelSedinMuSum += ChannelSedinMu3");
           calc(" ChannelSedinMuSum += ChannelSedinMu4");
           calc(" ChannelSedinMuSum += ChannelSedinMu5");

           mcalc(" ChannelTCSusp# = mif(ChannelSedinMuSum gt 0, "
                 "ChannelTCsusp * ChannelSedinMu#/ChannelSedinMuSum, 0) ",6);

           // **** deposition in channels *****

           _spatial(REAL4, ChannelDep0);
           _spatial(REAL4, ChannelDep1);
           _spatial(REAL4, ChannelDep2);
           _spatial(REAL4, ChannelDep3);
           _spatial(REAL4, ChannelDep4);
           _spatial(REAL4, ChannelDep5);
           mcalc(" ChannelDep# = min(0,ChannelTCSusp#-ChannelSedConc#)",6);
           mcalc(" ChannelDep# *= mif(ChannelWH gt 0, (1-exp(-SV#/(0.001*ChannelWH)*DTSEC))*ChannelVolin, 0) ",6);
           mcalc(" ChannelDep# = max(ChannelDep#,min(0,ChannelTCSusp#-ChannelSedConc#)*ChannelVolin),0) ",6);

//           mcalc(" ChannelSedin# += ChannelDF# ",6);
           mcalc(" ChannelSedinMu# += ChannelDep# ",6);
           mcalc(" ChannelSedinMu# = max(ChannelSedinMu#, 0) ",6);

       }// switch no erosion

// *****************************************************************************
// ******** kinematic wave calculation for CHANNEL flow ************************
// *****************************************************************************

       _nonspatial(REAL4, SumChannelVolin);
       calc(" SumChannelVolin = sum(ChannelVolin) ");
       // used in mass balance error

       _nonspatial(REAL4, ChannelSumSedin0);
       _nonspatial(REAL4, ChannelSumSedin1);
       _nonspatial(REAL4, ChannelSumSedin2);
       _nonspatial(REAL4, ChannelSumSedin3);
       _nonspatial(REAL4, ChannelSumSedin4);
       _nonspatial(REAL4, ChannelSumSedin5);
       mcalc(" ChannelSumSedin# = sum(ChannelSedinMu#) ",6);

      _spatial(REAL4, ChannelQsedin0);
      _spatial(REAL4, ChannelQsedin1);
      _spatial(REAL4, ChannelQsedin2);
      _spatial(REAL4, ChannelQsedin3);
      _spatial(REAL4, ChannelQsedin4);
      _spatial(REAL4, ChannelQsedin5);
      mcalc(" ChannelQsedin# = mif(ChannelVolin gt 0, ChannelSedin#*ChannelQin/ChannelVolin) ",6);

      _spatial(REAL4, ChannelQsedout0);
      _spatial(REAL4, ChannelQsedout1);
      _spatial(REAL4, ChannelQsedout2);
      _spatial(REAL4, ChannelQsedout3);
      _spatial(REAL4, ChannelQsedout4);
      _spatial(REAL4, ChannelQsedout5);

      _spatial(REAL4, ChannelQout);
      _spatial(REAL4, chinfil);
     if (SwitchChannelInfil)
       calc(" chinfil = -(ChannelKsat *  ChannelPerimeter /3600000.0) ");
       //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
     else
       calc(" chinfil = 0 ");

     //calc(" InfilVolinKinWave = -chinfil * DTSEC * DXc ");
     //VJ 050831 calc basic input infil vol kin wave
     //VJ 060817 added DXc because it is a volume m2/s * m * s = m3     

      minematic(ChannelQout, ChannelQin, chinfil,
            ChannelQsedout0, ChannelQsedout1, ChannelQsedout2, ChannelQsedout3, ChannelQsedout4, ChannelQsedout5,
            ChannelQsedin0, ChannelQsedin1, ChannelQsedin2, ChannelQsedin3, ChannelQsedin4, ChannelQsedin5,
            ChannelLDD, ChannelAlpha, ChannelBeta, DTSEC, DXc);
//VJ 050831 CHANGED: infil map changes during run: after the kin wave
// it contains the sum of all fluxes in the cell! Needed for infiltration calc

       calc(" ChannelQOutflow = sum(mif(Outlet1 eq 1, ChannelQout, 0))");
       mcalc(" ChannelSedOutflow# = sum(mif(Outlet1 eq 1, ChannelQsedout#, 0))",6);
       // kin wave can handle more pits now, so ensure output catchment is outlet1

       calc(" ChannelCSArea = ChannelAlpha*(ChannelQout**ChannelBeta)");
      // ChannelCSArea is A (cross section) in m2

      //moved to lisparaminit.cpp
      //       _spatial(REAL4, ChannelVolout);
       calc(" ChannelVolout = ChannelCSArea*DXc ");
       // m3

       calc(" SumChannelVolout = sum(ChannelVolout) ");
       // spatial total of water in channel, m3 , after kin wave

/*  VJ 051005 blocked this, works for channel so why not here
       if (SwitchCorrectMass)
       begin{
           _nonspatial(REAL4, ChannelError);
            calc(" ChannelError = SumChannelVolin-SumChannelVolout-ChannelQOutflow*DTSEC");
           if (SumChannelVolout > 0)
              calc(" ChannelVolout += ChannelError/SumChannelVolout*ChannelVolout ");
           calc(" ChannelVolout = max(ChannelVolout, 0) ");
           // divide error caused by pits over network
           calc(" SumChannelVolout = sum(ChannelVolout) ");
       }end
*/

//VJ 051005 added corrected channel infil calc
   if (SwitchKinwaveInfil)
   {
      //nonspatial infil adjustment in m3
      calc(" InfilVol += cover(min(chinfil*DTSEC+ChannelVolout, InfilVolinKinWave),0) ");
      //spatial infil adjustment in m3
	   calc(" TotalInfilVol += max(0,SumChannelVolin - SumChannelVolout - ChannelQOutflow*DTSEC) ");
	 }

   calc(" ChannelVolout = max(ChannelVolout, 0) ");
   calc(" ChannelVolin = ChannelVolout ");

   if (SwitchBuffers)
      calc(" SumBufferVolumeChannel = sum(BufferVolumeInitChannel) - sum(BufferVolumeCurrentChannel) ");

   if (!SwitchNoErosion)
   {
     // **** Calculation New Sediment amount Sedout
     _spatial(REAL4, Channelqfactor);
     calc(" Channelqfactor = mif(ChannelQout gt 0.00001, ChannelVolout/ChannelQout, 0) ");
     mcalc(" ChannelSedinMu# = Channelqfactor*ChannelQsedout# ",6);

     // reuse channelsedin
     if (SwitchCorrectMassSED) // default true !
     {
         _nonspatial(REAL4, ChannelSumSedout0);
         _nonspatial(REAL4, ChannelSumSedout1);
         _nonspatial(REAL4, ChannelSumSedout2);
         _nonspatial(REAL4, ChannelSumSedout3);
         _nonspatial(REAL4, ChannelSumSedout4);
         _nonspatial(REAL4, ChannelSumSedout5);
         mcalc(" ChannelSumSedout# = sum(ChannelSedinMu#) ",6);

         if (ChannelSumSedout0 > 0)
          calc(" ChannelSedinMu0 += ChannelSedinMu0*(ChannelSumSedin0-ChannelSumSedout0-ChannelSedOutflow0*DTSEC)/ChannelSumSedout0 ");
         if (ChannelSumSedout1 > 0)
          calc(" ChannelSedinMu1 += ChannelSedinMu1*(ChannelSumSedin1-ChannelSumSedout1-ChannelSedOutflow1*DTSEC)/ChannelSumSedout1 ");
         if (ChannelSumSedout2 > 0)
          calc(" ChannelSedinMu2 += ChannelSedinMu2*(ChannelSumSedin2-ChannelSumSedout2-ChannelSedOutflow2*DTSEC)/ChannelSumSedout2 ");
         if (ChannelSumSedout3 > 0)
          calc(" ChannelSedinMu3 += ChannelSedinMu3*(ChannelSumSedin3-ChannelSumSedout3-ChannelSedOutflow3*DTSEC)/ChannelSumSedout3 ");
         if (ChannelSumSedout4 > 0)
          calc(" ChannelSedinMu4 += ChannelSedinMu4*(ChannelSumSedin4-ChannelSumSedout4-ChannelSedOutflow4*DTSEC)/ChannelSumSedout4 ");
         if (ChannelSumSedout5 > 0)
          calc(" ChannelSedinMu5 += ChannelSedinMu5*(ChannelSumSedin5-ChannelSumSedout5-ChannelSedOutflow5*DTSEC)/ChannelSumSedout5 ");
      }


      mcalc(" ChannelSedinMu# = max(ChannelSedinMu#, 0) ",6);

      calc(" ChannelSedin = max(ChannelSedinMu0, 0) ");
      calc(" ChannelSedin += max(ChannelSedinMu1, 0) ");
      calc(" ChannelSedin += max(ChannelSedinMu2, 0) ");
      calc(" ChannelSedin += max(ChannelSedinMu3, 0) ");
      calc(" ChannelSedin += max(ChannelSedinMu4, 0) ");
      calc(" ChannelSedin += max(ChannelSedinMu5, 0) ");

      calc(" SumChannelSedout = sum(ChannelSedin) ");
         // needed in mass balance error

      if (SwitchBuffers)
         calc(" SumBufferSedVolumeChannel = sum(BufferSedVolumeCurrentChannel) ");
}//no erosion


