// METHOD 2
// **** Swatre model Soil Water using tables
// **** (C) Winand Staring Centre, modified by Utrecht University

     static SOIL_MODEL *soilModelcrust;
     static SOIL_MODEL *soilModelcrack;
     static SOIL_MODEL *soilModelwheel;
     static SOIL_MODEL *soilModelgrass;
     static SOIL_MODEL *soilModelcompact;
/*
     if (SwitchCompactPresent)
        calc(" WHCompact = WH ");
     if (SwitchCrustPresent)
        calc(" WHCrust = WH ");
     if (SwitchGrassPresent)
        calc(" WHGrass = WH ");
     //N.B. WHWheeltrack is done separately only when version is weeltracks

     _spatial(REAL4, InfilPotCrust);
     calc(" InfilPotCrust = 0 ");
     _spatial(REAL4, InfilPotCompact);
     calc(" InfilPotCompact = 0 ");
     _spatial(REAL4, InfilPotGrass);
     calc(" InfilPotGrass = 0 ");
*/

     /* INITIALIZE FIRST STEP */
     if (stepNr == 0)
     {
	  // check input
          LisIFace->Messages->Lines->Append("Preparing SWATRE soil profiles ...");
        _spatial_input(REAL4, PROFILE, mapname("profmap"));
        celltest(LDD, PROFILE);

        _spatial_input_if(REAL4, PROFWHEEL, mapname("profwltr"),SwitchWheelPresent);
        _spatial_input_if(REAL4, PROFCRUST, mapname("profcrst"),SwitchCrustPresent);
        _spatial_input_if(REAL4, PROFGRASS, mapname("profgras"),SwitchGrassPresent);
        _spatial_input_if(REAL4, PROFCOMP, mapname("profcomp"),SwitchCompactPresent);

        initsoil(mapname("profinp"), tableDir);
               // reading in soil physics tables

         SwatreInitialized = true;

	 timestepindex -= DTMIN; // adjust index for a moment

	 soilModelcrack = initswatre(PROFILE, mapname("inithead"),  HeadOut, RESPATH,minDtSwatre, precisionSwatre, timestepindex);
	     // PROFILE is the soil profile map;
	     // inithead.001 etc are the initial pressure heads, in path CatPath
	     // inithead are negative numbers, pF 1 = -10.000
	     // RESPATH and index are added in initswatre call

     if (SwitchCrustPresent)
     {
         celltest(LDD, PROFCRUST);
         calc(" PROFCRUST = mif(CrustFraction gt 0, PROFCRUST, -1) ");
         soilModelcrust = initswatre(PROFCRUST, mapname("inithead"),  NULL_MAP, RESPATH,
                                     minDtSwatre, precisionSwatre, timestepindex);
     }

     if (SwitchCompactPresent)
     {
         celltest(LDD, PROFCOMP);
         calc(" PROFCOMP = mif(CompactFraction gt 0, PROFCOMP, -1) ");
         soilModelgrass = initswatre(PROFCOMP, mapname("inithead"),  NULL_MAP, RESPATH,
                                     minDtSwatre, precisionSwatre, timestepindex);
     }

     if (SwitchGrassPresent)
     {
         celltest(LDD, PROFGRASS);
         calc(" PROFGRASS = mif(GrassWidth gt 0, PROFGRASS, -1) ");
         soilModelgrass = initswatre(PROFGRASS, mapname("inithead"),  NULL_MAP, RESPATH,
                                     minDtSwatre, precisionSwatre, timestepindex);
     }

     if (SwitchWheelPresent)
     {
         celltest(LDD, PROFWHEEL);
         calc(" PROFWHEEL = mif(WheelWidthDX gt 0, PROFWHEEL, -1) ");
             soilModelwheel = initswatre(PROFWHEEL, mapname("inithead"),  NULL_MAP, RESPATH,
                                     minDtSwatre, precisionSwatre, timestepindex);
     }

	 timestepindex += DTMIN; // undo index adjustment

     }

     swatrestep(soilModelcrack, WH, InfilPot, HeadOut, RESPATH, DTDAY, timestepindex,
               ksatCalibration, SwitchGeometricMean);
     if (SwitchCrustPresent)
        swatrestep(soilModelcrust, WHCrust, InfilPotCrust, NULL_MAP, RESPATH, DTDAY, timestepindex,
                   ksatCalibration, SwitchGeometricMean);
     if (SwitchCompactPresent)
        swatrestep(soilModelcompact, WHCompact, InfilPotCompact, NULL_MAP, RESPATH, DTDAY, timestepindex,
                   ksatCalibration, SwitchGeometricMean);
     if (SwitchGrassPresent)
        swatrestep(soilModelgrass, WHGrass, InfilPotGrass, NULL_MAP, RESPATH, DTDAY, timestepindex,
                   ksatCalibration, SwitchGeometricMean);

     if (SwitchWheelPresent)
        swatrestep(soilModelwheel, WHWheelTrack, InfilPotWheel, NULL_MAP, RESPATH, DTDAY, timestepindex,
                   ksatCalibration, SwitchGeometricMean);


/*
MOVED to LisPInfiltration
     calc(" WH = WH*(1-CrustFraction-GrassFraction-CompactFraction) + WHGrass*GrassFraction "
           "      + WHCrust*CrustFraction + WHCompact*CompactFraction ");

     if (SwitchGrassPresent)
     begin{
        _spatial(REAL4, oldWH);
        calc(" oldWH = WH ");
        calc(" WH = mif(GrassWidth gt 0, max(0, WH+InfilPotGrass), WH) ");
        calc(" InfilPotGrass += mif(GrassWidth gt 0, oldWH-WH,0) ");
     }end

        // adjestment of WH
     calc(" InfilPot = InfilPot*(1-CrustFraction-GrassFraction-CompactFraction) + InfilPotGrass*GrassFraction "
            "         + InfilPotCrust*CrustFraction + InfilPotCompact*CompactFraction");
        // NB infilpot is negative
*/
