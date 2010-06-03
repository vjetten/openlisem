// *****************************************************************************
// ****** INFILTRATION *********************************************************
// *****************************************************************************

      _spatial(REAL4, WHinf);
      calc(" WHinf = WH ");
        // WH before infiltration

      if (SwitchWheelPresent)
      {
          _spatial(REAL4, WHWTinf);
          calc(" WHWTinf = WHWheelTrack ");
      }

     /* total volume before infiltration */
//VJ 030701 no difference WT vars are 0 anyway
     _spatial(REAL4, InfilVolin);
//     if (SwitchWheelPresent)
      calc(" InfilVolin = (RoadWidthDX*WHRoad+WheelWidthDX*WHWheelTrack+(SoilWidthDX+StoneWidthDX)*WH)*DXc/1000 ");
//     else
//      calc(" InfilVolin = (RoadWidthDX*WHRoad+(SoilWidthDX+StoneWidthDX)*WH)*DXc/1000 ");
        // rainfall/water on channels is treated in the channel section of LISEM

     calc(" WHCompact = WH ");
     calc(" WHCrust = WH ");
     calc(" WHGrass = WH ");
     //N.B. WHWheeltrack is done separately only when version is wheeltracks

     // impermeable catchment, for testing
     if (INFIL_METHOD == INFIL_NONE)
        ReportV("no infiltration calculated, 100 percent runoff possible");

     // **** SWATRE all soil types
     if (INFIL_METHOD == INFIL_SWATRE)
     {
         ReportV("using SWATRE infiltration model");
         #include "lisPinfswa.cpp"

         if (SwitchNutrients)
            swatretheta(soilModelcrack, thetatop, 1, 1);
           // average humidity from first 3 layers
     }

     // **** Holtan/Overton infiltration
     if (INFIL_METHOD == INFIL_HOLTAN)
     {
         ReportV("using Holtan infiltration model");
         #include "lisPinfhol.cpp"
         if (SwitchNutrients)
            calc(" thetatop = ThetaS1 ");
     }


     // **** Green/Ampt infiltration

     if (INFIL_METHOD == INFIL_GREENAMPT || INFIL_METHOD == INFIL_GREENAMPT2)
     {
         ReportV("using Green & Ampt infiltration model");

         #include "lisPinfGA.cpp"
         if (SwitchGrassPresent)
         {
             #include "lisPinfGAgr.cpp"
         }
         if (SwitchWheelPresent)
         {
             #include "lisPinfGAwt.cpp"
         }

         if (SwitchNutrients)
            calc(" thetatop = ThetaS1 ");
           // assumption that moisture content for nutrients equals saturation
//         DTHOUR = temp;
     }
     // **** simple subtraction of Ksat!, for testing
     if (INFIL_METHOD == INFIL_KSAT)
     {
          ReportV("using Ksat subtraction");

          _spatial(REAL4, Ksateff);
          calc(" Ksateff = Ksat1*(1-CrustFraction-GrassFraction-CompactFraction)");
            //spatial average of ksat, all fractions are zero when not existant !!!

          if (SwitchCrustPresent)
             calc(" Ksateff += KsatCR*CrustFraction ");

          if (SwitchCompactPresent)
             calc(" Ksateff += KsatComp*CompactFraction ");

          if (SwitchGrassPresent)
             calc(" Ksateff += KsatGR*GrassFraction ");

          calc(" InfilPot =  - Ksateff * DTHOUR");
          calc(" WH = max(WH - Ksateff * DTHOUR, 0) ");

     }


// WH is now the water height after infiltration !!!

//VJ 040823 include buffers, no infil in buffer cells
     if (SwitchBuffers){
        calc(" WH = mif(BufferID gt 0, WHinf, WH) ");
        calc(" WHCompact = mif(BufferID gt 0, WHinf, WHCompact) ");
        calc(" WHCrust = mif(BufferID gt 0, WHinf, WHCrust) ");
        calc(" WHGrass = mif(BufferID gt 0, WHinf, WHGrass) ");
     }

//we now have potentially five waterheights: WH, WHGrass, WHCrust, WHCompact, WHWheel
//in case of WHWheel, it is a separate channel and already accounted for in the fractions
// if wheeltrack version is not used it is simply a fraction compacted, so wheels are kept separate

//calculate weigthed average WH


     calc(" WH = mif(hardsurface le 0, \
            WH*(1-CrustFraction-GrassFraction-CompactFraction) + WHGrass*GrassFraction \
            + WHCrust*CrustFraction + WHCompact*CompactFraction, 0) ");

        // new gridcell WH as weighted average of surface type fractions and WHs

//MAXIMIZE infiluence of grassstrip, in the sense that when a gstrip is present, its potential infil
// determines that of the cell
     if (SwitchGrassPresent)
     begin{
        _spatial(REAL4, oldWH);
        calc(" oldWH = WH ");
        calc(" WH = mif(GrassWidth gt 0, max(0, WH+InfilPotGrass), WH) ");
        calc(" InfilPotGrass = mif(GrassWidth gt 0, InfilPotGrass+oldWH-WH,0) ");
     }end

//calc average potential infiltration needed for surplus kin wave
     calc(" InfilPot = InfilPot*(1-CrustFraction-GrassFraction-CompactFraction) + InfilPotGrass*GrassFraction "
            "         + InfilPotCrust*CrustFraction + InfilPotCompact*CompactFraction");
        // NB infilpot is negative, needed in kin wave infiltration surplus


     calc(" WHinf -= WH ");
          // WHinf is now average infiltrated water, needed also for nutrients

     calc(" InfilSurplus = WHinf + min(InfilPot, 0) ");
          // remaining AVERAGE infil needed for kin wave: actual (WHinf) - potential infiltration (InfilPot)
          // infilpot and infilsurplus are NEGATIVE values

      if (SwitchWheelPresent)
      {
          calc(" WHWTinf -= WHWheelTrack ");
          calc(" InfilSurplusWT = WHWTinf + min(InfilPotWheel, 0) ");
      }
      //wheeltracks separate

     _spatial(REAL4, InfilVolout);
     calc(" InfilVolout = (RoadWidthDX*WHRoad+WheelWidthDX*WHWheelTrack+(SoilWidthDX+StoneWidthDX)*WH)*DXc/1000 ");
     //volume in cell after infiltration, in m3

     calc(" TotalInfilVol += sum(InfilVolin - InfilVolout) "); /* + ChannelInfilVol */
     // non spatial infil total. This has to adjusted for infil in kin wave
     // needed for mass balance and screen output

     calc(" InfilVol += InfilVolin - InfilVolout");// + ChannelInfilVol) ");
     //spatial infiltration volume needed for map output

       // rainfall/water on channels is treated in the channel section of LISEM
       // volume infiltrated water in m3
       // used in mass balance check


//CHECK where is channel infil accounted for?
