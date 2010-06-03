// *****************************************************************************
// ******** kinematic wave calculation for overland flow ***********************
// *****************************************************************************

     // **** recalc Alpha and Qin after tochannel

      calc(" Alpha = (N/sqrt(Gradient)*(FlowWidth+2*RunoffMeanHin/1000)**(2.0/3.0))**Beta ");
      calc(" Qin = mif(Alpha gt 0 , ((FlowWidth * RunoffMeanHin/1000)/Alpha)**(1/Beta), 0) ");
        // recalc because of fraction to wheel

       //calc(" WH -= RunoffHin ");
       calc(" WH = SurfStorH ");
          // the waterheight is reduced to the retention
          // the new detention RunoffMeanHout will be calculated, based on the new
          // Qout, and will be added after that to WH

//       calc(" SumDeprStoreVol = sum(WH*(SoilWidthDX+StoneWidthDX)*DXc/1000.0) ");
       calc(" SumDeprStoreVol = sum(SurfStorH*(SoilWidthDX+StoneWidthDX)*DXc/1000.0) ");
          // spatial total depression storage for mass balance, for one DT

     _spatial(REAL4, DeprSedin0);
     _spatial(REAL4, DeprSedin1);
     _spatial(REAL4, DeprSedin2);
     _spatial(REAL4, DeprSedin3);
     _spatial(REAL4, DeprSedin4);
     _spatial(REAL4, DeprSedin5);
     mcalc(" DeprSedin# = mif(WH gt 0, SedinMu#*SurfStorH/WH,0) ",6);


       calc(" SumWaterHVolin = sum(WaterHVolin) ");
        // needed for infil in kin wave

       calc(" WHRoad = 0 ");
          // all water on roads contributes to flow

       _nonspatial(REAL4, SumSedin0);
       _nonspatial(REAL4, SumSedin1);
       _nonspatial(REAL4, SumSedin2);
       _nonspatial(REAL4, SumSedin3);
       _nonspatial(REAL4, SumSedin4);
       _nonspatial(REAL4, SumSedin5);
       mcalc(" SumSedin# = sum(SedinMu#) ",6);

        calc(" SumSedout = SumSedin0 ");
        calc(" SumSedout += SumSedin1 ");
        calc(" SumSedout += SumSedin2 ");
        calc(" SumSedout += SumSedin3 ");
        calc(" SumSedout += SumSedin4 ");
        calc(" SumSedout += SumSedin5 ");

      _spatial(REAL4, Qsedin0);
      _spatial(REAL4, Qsedin1);
      _spatial(REAL4, Qsedin2);
      _spatial(REAL4, Qsedin3);
      _spatial(REAL4, Qsedin4);
      _spatial(REAL4, Qsedin5);
      mcalc(" Qsedin# = mif(WaterHVolin gt 0, Qin * SedinMu#/WaterHVolin, 0) ",6);

      _spatial(REAL4, Qsedout0);
      _spatial(REAL4, Qsedout1);
      _spatial(REAL4, Qsedout2);
      _spatial(REAL4, Qsedout3);
      _spatial(REAL4, Qsedout4);
      _spatial(REAL4, Qsedout5);
      mcalc(" Qsedout# = 0",6);

      _spatial(REAL4, Qout);

      _spatial(REAL4, infil);
      if (SwitchKinwaveInfil)
     //VJ 050704 changed DXc to DX
        calc(" infil = InfilSurplus/(1000*DTSEC)*DX ");
      else
        calc(" infil = 0 ");

      minematic(Qout, Qin, infil,
                Qsedout0, Qsedout1, Qsedout2, Qsedout3, Qsedout4, Qsedout5,
                Qsedin0, Qsedin1, Qsedin2, Qsedin3, Qsedin4, Qsedin5,
                LDD, Alpha, Beta, DTSEC, DXc);

      calc(" QOutflow = sum(mif(Outlet1 eq 1, Qout, 0))");
      mcalc(" SedOutflow# = sum(mif(Outlet1 eq 1, Qsedout#, 0))",6);
        //outlet

      _spatial(REAL4, RunoffAout);
      calc(" RunoffAout = Alpha*(Qout**Beta)");


      _spatial(REAL4, DXcorr);
      calc(" DXcorr = DX ");

//      if (SwitchWheelPresent)
//        calc(" DXcorr -= WheelWidthDX ");

      if (SwitchIncludeChannel)
        calc(" DXcorr -= cover(ChannelWidthUpDX,0) ");

      calc(" DXcorr = max(DXcorr, DX*0.1) ");

      _spatial(REAL4, RunoffHout);
      calc(" RunoffHout = 1000 * RunoffAout/DXcorr ");
      // RunoffHout is the average water depth in the pixel, in mm

      calc(" WH += RunoffHout");
      calc(" WHRoad = RunoffHout");
      calc(" WHRoad = RunoffHout");
       // add new WH to surfaces
      //hier


      _spatial(REAL4, WaterHVolout);
      calc(" WaterHVolout = (RoadWidthDX*WHRoad+(SoilWidthDX+StoneWidthDX)*WH)*DXc/1000 ");

      calc(" SumWaterHVolout = sum(WaterHVolout) ");
      // spatial total amount of water on surface in m3 after kin wave
      // used in water mass balance error

     //VJ 051005 made changes here
     if (SwitchKinwaveInfil){
//        calc(" TotalInfilVol += sum(min(infil*DTSEC+WaterHVolin, InfilVolinKinWave)) ");
        //nonspatial infil adjustment in m3
        calc(" InfilVol += min(infil*DTSEC+WaterHVolout, InfilVolinKinWave) ");
        //spatial infil adjustment in m3
       // old stuff
       calc(" TotalInfilVol += max(0,SumWaterHVolin - SumWaterHVolout - QOutflow*DTSEC) ");
       // add volume infiltrated in kin wave to suminfiltration (NOTE, this includes possible errors)
     }

       if (!SwitchNoErosion)
       {
         // **** Calculation New Sediment amount Sedout
         _spatial(REAL4, qfactor);
       calc(" qfactor = mif(Qout gt 0, WaterHVolout/Qout, 0) ");
        // klopt beter
      //     calc(" qfactor = mif(Qout gt 0,RunoffAout*DXc/Qout, 0) ");

         mcalc(" SedinMu# = qfactor*Qsedout# ",6);
         mcalc(" SedinMu# += DeprSedin# ",6);

         _nonspatial(REAL4, SumSedout0);
         _nonspatial(REAL4, SumSedout1);
         _nonspatial(REAL4, SumSedout2);
         _nonspatial(REAL4, SumSedout3);
         _nonspatial(REAL4, SumSedout4);
         _nonspatial(REAL4, SumSedout5);
         mcalc(" SumSedout# = sum(SedinMu#) ",6);

         if (SwitchCorrectMassSED) // default true !
         {
             if (SumSedout0 > 0)
              calc(" SedinMu0 += SedinMu0*(SumSedin0-SumSedout0-SedOutflow0*DTSEC)/SumSedout0 ");
             if (SumSedout1 > 0)
              calc(" SedinMu1 += SedinMu1*(SumSedin1-SumSedout1-SedOutflow1*DTSEC)/SumSedout1 ");
             if (SumSedout2 > 0)
              calc(" SedinMu2 += SedinMu2*(SumSedin2-SumSedout2-SedOutflow2*DTSEC)/SumSedout2 ");
             if (SumSedout3 > 0)
              calc(" SedinMu3 += SedinMu3*(SumSedin3-SumSedout3-SedOutflow3*DTSEC)/SumSedout3 ");
             if (SumSedout4 > 0)
              calc(" SedinMu4 += SedinMu4*(SumSedin4-SumSedout4-SedOutflow4*DTSEC)/SumSedout4 ");
             if (SumSedout5 > 0)
              calc(" SedinMu5 += SedinMu5*(SumSedin5-SumSedout5-SedOutflow5*DTSEC)/SumSedout5 ");
         }
         mcalc(" SedinMu# = max(SedinMu#, 0)",6);

          calc(" Sedin = SedinMu0 ");
          calc(" Sedin += SedinMu1 ");
          calc(" Sedin += SedinMu2 ");
          calc(" Sedin += SedinMu3 ");
          calc(" Sedin += SedinMu4 ");
          calc(" Sedin += SedinMu5 ");
          // needed for screen output

          calc(" SumSedout = sum(Sedin) ");
          // needed for mass balance
  }//no erosion


