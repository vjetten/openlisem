//******************************************************************************
//****** INPUT of runfile VARIABLES AND MAPS ***********************************
//******************************************************************************


     //-------------------------------------------------------------------------
     //******************variable declarations *********************************
     //-------------------------------------------------------------------------

     _nonspatial(REAL4, STARTINTERVAL);   // in minutes
     _nonspatial(REAL4, ENDINTERVAL);     // in minutes
     _nonspatial(REAL4, DTSEC);           // timeinterval in seconds

     _nonspatial(REAL4, DX);              // elementsize in meters
     _nonspatial(REAL4, DTMIN);           // timeinterval in minutes
     _nonspatial(REAL4, DTHOUR);           // timeinterval in days
     _nonspatial(REAL4, DTDAY);           // timeinterval in days

     _nonspatial(REAL4, NR_VALS_START);  // number of cells that have not MV
     _nonspatial(REAL4, NR_VALS_NOW);

     //-------------------------------------------------------------------------
     //******************get directories and paths******************************
     //-------------------------------------------------------------------------

     if (SwitchOutputTimeStep)
     {
        PERIOD = GetInt("Output interval");
        PERIODMAPS = TRUE;
     }

     if (SwitchOutputTimeUser)
     {
        for (int i = 0; i < LisIFace->E_OutputTimes->Lines->Count; i++)
             WRITETIME[i] = LisIFace->E_OutputTimes->Lines->Strings[i].ToDouble();
     }




     //-------------------------------------------------------------------------
     //******************get run time and timestep *****************************
     //-------------------------------------------------------------------------

     STARTINTERVAL = (REAL4)GetFloat("Begin time");
     ENDINTERVAL =  (REAL4)GetFloat("End time");
     DTSEC = (REAL4)GetFloat("Timestep");
     timestep = DTSEC/60.0;
     // needed for output window

     DTMIN = DTSEC/60.0;
     DTHOUR = DTSEC/3600.0;
     DTDAY = DTSEC/(24.0*60.0*60.0);

     //VJ 080522 multi timestep GA solution
      _nonspatial(REAL4, DTHOURGA);           // timeinterval in seconds for GA infiltration
      _nonspatial(INT4, nrGAsteps);
      nrGAsteps = DTSEC > 10 ? (INT4) DTSEC/10.0 + 1 : 1;
      DTHOURGA = DTHOUR/(REAL4)nrGAsteps;

     //-------------------------------------------------------------------------
     //******************get other parameters***********************************
     //-------------------------------------------------------------------------

     _nonspatial(REAL4, CriticalStreamPower);
     CriticalStreamPower = 0.4;
        // critical unit stream power, in cm/s!
        // 0.4, according to Govers, 1990

     _nonspatial(REAL4, SplashDelivery);
     SplashDelivery = (REAL4)GetFloat("Splash Delivery Ratio");
//     SplashDelivery = 0.1;
     // Splash delivery ratio, determining the fraction of splash
     // transported into the flow system

     _nonspatial(REAL4, StripN);
     StripN = (REAL4)GetFloat("Grassstrip Mannings n");
        // Manning's n for grass-strips on slopes as a soil conservation tool


     //-------------------------------------------------------------------------
     // ******** initialize mask and DX    *************
     //-------------------------------------------------------------------------

//     InitMask(mapname("area"));
     InitMask(mapname("LDD"));
     //VJ 090521 changed this from ID map to LDD

     DX = RgiveCellSizeX();
     // initialize immediately after InitMask

     _nonspatial(REAL4, TMP_RANGE);
     // parameter used for special range checks
     _spatial(REAL4, TMP_VALUE);
     // parameter used for special range checks


     //-------------------------------------------------------------------------
     // ******** drainage basin morphology *************
     //-------------------------------------------------------------------------

       // VJ-> NOT ANY MORE: type LDD has a check on nr. of pits
       _spatial_input(LDD, LDD, mapname("LDD"));
       calc(" NR_VALS_START = count(mif(LDD gt 0, 1, 0))");
       //VJ 090208 most important map, all is checked agains LDD or chanLDD

//VJ 090208
// AREA no longer needed, function replaced by LDD
//     _spatial_input(UINT1, AREA, mapname("ID"));//Area"));
//     calc(" NR_VALS_START = count(AREA)");

      _spatial_input(UINT1, RainGaugeID, mapname("ID"));
      celltest(LDD, RainGaugeID);

      //     calc(" NR_VALS_START = count(ID)");//AREA)");
 //VJ 080423 snowmelt
      _spatial_input_if(UINT1, SnowID, mapname("SnowID"),SwitchSnowmelt);
      _spatial(REAL4, SnowCover);
      calc(" SnowCover = 0 ");
      if (SwitchSnowmelt)
      {
          celltest(LDD, SnowID);
          calc(" SnowCover = mif(SnowID gt 0, 1.0, 0.0) ");
      }

      _spatial_input(REAL4, Gradient,mapname("Grad"));
      celltest(LDD, Gradient);
      rangetest(Gradient,R_GE_LE,0.0001,1,"Slope gradient values must be sin(angle)");

      _spatial(REAL4, DXc);
      calc(" DXc = DX/cos(asin(Gradient)) ");

      _spatial(REAL4, DXx);
      calc(" DXx = DX ");
      //spatial variant of DX for testing

     _spatial_input(REAL4, Outlet, mapname("Outlet"));
     celltest(LDD, Outlet);

     _spatial_input(REAL4, OutPoint, mapname("outpoint"));
     celltest(LDD, OutPoint);

      /*
     {
       char p[256];
       strcpy(outPointFileName, CatPath(p, RESPATH));
       strcat(outPointFileName, mapnameonly("outpoint"));
       strcat(outPointFileName, ".csv");
     }
     */

     _spatial_input(REAL4, RoadWidthDX, mapname("Road"));
     celltest(LDD, RoadWidthDX);
     rangetest(RoadWidthDX,R_GE_LT, 0, DX,"@@Width of roads");
     calc(" NR_VALS_NOW = count(RoadWidthDX)");
     if (NR_VALS_NOW < NR_VALS_START)
        LisemError("wrong number of pixels in Roadwidt.map, non roads must have value 0!");

     _spatial(REAL4, hardsurface);
     calc (" hardsurface = 0.0 ");

     if (SwitchHardsurface)
     {
        _spatial_input(REAL4, hardsurf,mapname("hardsurf"));
        celltest(LDD, hardsurf);
        rangetest(hardsurf,R_GE_LE,0,1, "hard surface");
        calc (" hardsurface = 1.0*hardsurf ");
     }
       //VJ 080613 include hard surfaces
       //VJ 091216 moved to general, influences infiltration too!!!


     //-------------------------------------------------------------------------
     // ******** land use variables ********************
     //-------------------------------------------------------------------------

//VJ 100116 Interception
     // Van Hoyningen-Huene (1981), p.46; maximum interception (mm)
/*
Orig. LISEM (crops): S = 0.935+0.498*LAI-0.00575*LAI^2
Pinus:         S = 0.2331*LAI          (n=12,R2=0.88)
Douglas Fir:   S = 0.3165*LAI          (n=4, R2=0.83)
Olive:         S = 1.46 * LAI^0.56     (n=5, R2=0.87)
Eucalypt:      S = 0.0918*LAI^1.04     (n=8, R2=0.51)
Rainforest:    S = 0.2856*LAI          (n=5, R2=0.60)
Bracken:       S = 0.1713*LAI          (n=8, R2=0.98)
Clumped grass: S = 0.59 * LAI^0.88     (n=6, R2=0.82)
*/

     _spatial_input(REAL4, LAI,mapname("LAI"));
     celltest(LDD, LAI);
     rangetest(LAI,R_GE_LE,0,12, "Leaf area index");
         // LAI should be entered for the crop/vegetation,
         // not as a pixel average, VegetatFraction deals with that

     _spatial(REAL4, CanopyStorage);
     calc(" CanopyStorage = 0");
     InterceptionLAIType = 0;

     if (!SwitchInterceptionLAI)
     {
        _spatial_input(REAL4, Smax,mapname("smax"));
        celltest(LDD, Smax);
        rangetest(Smax,R_GE_LE,0,100, "Max Canopy Storage");
        calc(" CanopyStorage = Smax ");
     }
     else
     {
         InterceptionLAIType = GetInt("Canopy storage equation");
     }

     if (SwitchInterceptionLAI)
     {

         switch (InterceptionLAIType)
         {
          case 0: calc(" CanopyStorage = 0.935 + (0.498*LAI) - (0.00575 * sqr(LAI) ) ");break;
          case 1: calc(" CanopyStorage = 0.2331 * LAI "); break;
          case 2: calc(" CanopyStorage = 0.3165 * LAI "); break;
          case 3: calc(" CanopyStorage = 1.46 * (LAI^0.56)"); break;
          case 4: calc(" CanopyStorage = 0.0918 * (LAI^1.04) "); break;
          case 5: calc(" CanopyStorage = 0.2856 * LAI "); break;
          case 6: calc(" CanopyStorage = 0.1713 * LAI "); break;
          case 7: calc(" CanopyStorage = 0.59 * (LAI^0.88) "); break;
         }
     }

     _nonspatial(REAL4, StemflowFraction);
     StemflowFraction = GetFloat("Stemflow fraction");

     _spatial_input(REAL4, VegetatFraction, mapname("cover"));
     celltest(LDD, VegetatFraction);
     rangetest(VegetatFraction,R_GE_LE,0,1, "Soil coverage");

     _spatial_input(REAL4, CropHeight,mapname("CH"));
     celltest(LDD, CropHeight);
     rangetest(CropHeight,R_GE_LE,0,30, "Crop height (m)");

     //-------------------------------------------------------------------------
     // ******** soil surface variables ****************
     //-------------------------------------------------------------------------

     _spatial_input(REAL4, RR,mapname("RR"));
     celltest(LDD, RR);
     rangetest(RR,R_GE_LE,0,10, "Random roughness");

     _spatial_input(REAL4, N,mapname("manning"));
     celltest(LDD, N);
     rangetest(N,R_GE_LE, 0.00001,10, "Manning's n");

     nCalibration = -1.0;
     nCalibration = GetFloat("N calibration");
     _nonspatial(REAL4, n_cal);
     n_cal = nCalibration;
     if (nCalibration >= 0)
     {
       calc(" N = N * n_cal ");
     }

     _spatial_input(REAL4, StoneFraction,mapname("Stonefrc"));
     celltest(LDD, StoneFraction);
     rangetest(StoneFraction,R_GE_LE, 0, 1, "Stone cover fraction");
         // enter the fraction of the soil surface covered with stones
         // value between 0 and 1; 1 is fully covered with stones

     _spatial(REAL4, CrustFraction);
     calc(" CrustFraction = 0 ");
     if (SwitchInfilCrust)
     {
       _spatial_input(REAL4, CrustFraction1, mapname("Crustfrc"));
       celltest(LDD, CrustFraction1);
       rangetest(CrustFraction1,R_GE_LE, 0, 1, "Fraction of crusted soils");
       calc(" CrustFraction1 = min(max(CrustFraction1, 0),1) ");
       calc(" TMP_RANGE = sum(CrustFraction1) ");
       SwitchCrustPresent = (TMP_RANGE > 0) && SwitchInfilCrust;
       if (!SwitchCrustPresent)
          LisemError("Cannot simulate crusts, crust fraction map empty");
       else
       {
          calc("CrustFraction = CrustFraction1 ");
       }
     }

     _spatial(REAL4, CompactFraction);
     calc(" CompactFraction = 0 ");
     if (SwitchInfilCompact)
     {
       _spatial_input(REAL4, CompactFraction1, mapname("Compfrc"));
       celltest(LDD, CompactFraction1);
       rangetest(CompactFraction1,R_GE_LE, 0, 1, "Fraction of compacted soils");
       calc(" CompactFraction1 = min(max(CompactFraction1, 0),1) ");
       calc(" TMP_RANGE = sum(CompactFraction1) ");
       SwitchCompactPresent = (TMP_RANGE > 0) && SwitchInfilCompact;
       if (!SwitchCompactPresent)
          LisemError("Cannot simulate compacted areas, compact fraction map empty");

       if (!SwitchCompactPresent)
          LisemError("Cannot simulate compacted soil, Compacted fraction map empty");
       else
       {
          calc("CompactFraction = CompactFraction1 ");
       }
     }

     _spatial(REAL4, GrassWidth);
     calc("GrassWidth = 0");
     if (SwitchInfilGrass)
     {
       _spatial_input(REAL4, GrassWidth1,mapname("grasswidth"));
       celltest(LDD, GrassWidth1);
       rangetest(GrassWidth1,R_GE_LE, 0, DX, "@@Grass-strips");
       // enter the width of field strips or grassed waterways in m
       calc(" TMP_RANGE = sum(GrassWidth1)");
       SwitchGrassPresent = (TMP_RANGE > 0) && SwitchInfilGrass;
       if (!SwitchGrassPresent)
          LisemError("Cannot simulate grass strips, grass width map empty");
       else
       {
          calc("GrassWidth = GrassWidth1 ");
       }

     }

       // needed for swatre
     _spatial(REAL4, GrassFraction);
     calc(" GrassFraction = min(1, GrassWidth/DX) ");
     // VJ 040218: added max 1

     // ************************************************************************
     // ******** GRASS-STRIPS **************************************************
     // ************************************************************************

     _spatial(REAL4, N_org);
     //VJ 090529 needed for sed trap, Nstrip stops when trap is full
     calc("N_org = N");
     if (SwitchGrassPresent)
        calc("N = mif(GrassWidth gt 0, N*(1-GrassFraction)+StripN*GrassFraction, N)");

     //-------------------------------------------------------------------------
     // ******** Buffers ********************
     //-------------------------------------------------------------------------
//VJ 040823 include buffers
     _spatial_input_if(UINT1, BufferID, mapname("bufferID"), SwitchBuffers);
     _spatial_input_if(REAL4, BufferVolume, mapname("bufferVolume"), SwitchBuffers);
//     _spatial_input_if(REAL4, BufferArea, mapname("bufferArea"), SwitchBuffers);
//     _spatial_input_if(REAL4, BufferDis, mapname("bufferdischarge"), SwitchBuffers);
     _nonspatial(REAL4, BufferSedBulkDensity);
     BufferSedBulkDensity = GetFloat("Sediment bulk density");

     if (SwitchBuffers)
     {
        celltest(LDD, BufferID);
        celltest(LDD, BufferVolume);
  //      celltest(LDD, BufferArea);
  //      celltest(LDD, BufferDis);
        calc(" BufferID = cover(BufferID,0) ");
        calc(" BufferVolume = cover(BufferVolume,0) ");
    //    calc(" BufferArea = cover(BufferArea,0) ");
   //     calc(" BufferDis = cover(BufferDis,0) ");
     }

     //-------------------------------------------------------------------------
     // ******** Erosion deposition variables **********
     //-------------------------------------------------------------------------

    if (!SwitchNoErosion)
    {
       _spatial_input(REAL4, AggregateStab,mapname("AggrStab"));
       celltest(LDD, AggregateStab);
       rangetest(AggregateStab,R_LE, R_DUMMY, 200, "Aggregate stability index");
       calc("TMP_RANGE=mincell(AggregateStab)");
       if (TMP_RANGE <= 0)
          LisemWarning("aggregate stability values <= 0.0 found, for those cells COHESION is used for splash detachment");

       _spatial_input(REAL4, D50,mapname("D50"));
       celltest(LDD, D50);
       calc("TMP_RANGE=maxcell(D50)");
       if (TMP_RANGE > 300)
           LisemWarning("Transport equations not suited for D50 > 300 mu, check D50.MAP");
       calc(" TMP_RANGE = mincell(D50) ");
       if (TMP_RANGE < 5)
           LisemWarning("Transport equations not suited for D50 < 5 mu, check D50.MAP");

       _spatial_input(REAL4, CohesionSoil,mapname("Coh"));
       celltest(LDD, CohesionSoil);
       _spatial_input(REAL4, CohesionRoot,mapname("CohAdd"));
       celltest(LDD, CohesionRoot);

       _spatial(REAL4, CohesionTotal);
       calc(" CohesionTotal = CohesionSoil+CohesionRoot ");
     //  rangetest(CohesionTotal,R_GE, 0.196, R_DUMMY,
     //      "@Cohesion values must be > 0.196 (detachment efficiency coef. Y < 1) "
     //      " error in CHANCOH.MAP (enter large values, e.g. 9999 for non-erodible surfaces)");
    }



