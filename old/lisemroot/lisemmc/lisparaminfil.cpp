     //******************************************************
     //***** INFILTRATION INPUT AND INIT ********************
     //******************************************************

     //-------------------------------------------------------------------------
     //******************get SWATRE infil parameters****************************
     //-------------------------------------------------------------------------

     INFIL_METHOD = GetInt("Infil Method");

     ksatCalibration = GetFloat("Ksat calibration"); ///100.0;
  //   if (!runv3)
    //    ksatCalibration /= 100;

     _nonspatial(REAL4, precisionSwatre);
     precisionSwatre = 5;

     _nonspatial(REAL4, minDtSwatre);
        // minimum timestep for Swatre submodel, in days
     if (INFIL_METHOD == INFIL_SWATRE)
     {
        minDtSwatre = GetFloat("SWATRE internal minimum timestep"); 
        if (minDtSwatre > DTSEC)
           minDtSwatre = DTSEC;
        minDtSwatre/=86400;
     }

     //-------------------------------------------------------------------------

     _spatial(REAL4, InfilSurplus);
     calc(" InfilSurplus = 0 ");
       // surplus for kin wave

     _spatial(REAL4, InfilSurplusWT);
     calc(" InfilSurplusWT = 0 ");

     _spatial(REAL4, InfilPot);
     calc(" InfilPot = 0 ");
     _spatial(REAL4, InfilPotCrust);
     calc(" InfilPotCrust = 0 ");
     _spatial(REAL4, InfilPotCompact);
     calc(" InfilPotCompact = 0 ");
     _spatial(REAL4, InfilPotGrass);
     calc(" InfilPotGrass = 0 ");
       // potential infiltration needed for infilsurplus

     if (SwitchDumphead && (INFIL_METHOD != INFIL_SWATRE))
        SwitchDumphead = false;
     _spatial_input_if(REAL4, HeadOut,mapname("headout"),SwitchDumphead);
     if (SwitchDumphead)
     {
         celltest(LDD, HeadOut);
     }

     _nonspatial(REAL4, ksatcal);
     ksatcal = ksatCalibration;

     if (INFIL_METHOD != INFIL_SWATRE && INFIL_METHOD != INFIL_NONE)
     {
        _spatial_input(REAL4, Ksat1, mapname("ksat1"));
        celltest(LDD, Ksat1);
        rangetest(Ksat1,R_GT, 0, R_DUMMY,"@Ksat layer 1 must be > 0 mm/h");
        calc("Ksat1 *= ksatcal ");

        _spatial_input_if(REAL4, KsatWT, mapname("ksatwt"), SwitchWheelPresent);
        if (SwitchWheelPresent)
        {
           celltest(LDD, KsatWT);
           calc("KsatWT *= ksatcal ");
        }

        _spatial_input_if(REAL4, KsatComp, mapname("ksatcomp"), SwitchCompactPresent);
        if (SwitchCompactPresent)
        {
           celltest(LDD, KsatComp);
           calc("KsatComp *= ksatcal ");
        }

       _spatial_input_if(REAL4, KsatGR, mapname("ksatgras"), SwitchGrassPresent);
       if (SwitchGrassPresent)
       {
           celltest(LDD, KsatGR);
           calc("KsatGR *= ksatcal ");
       }

        _spatial_input_if(REAL4, KsatCR, mapname("ksatcrst"), SwitchCrustPresent);
        if (SwitchCrustPresent)
        {
           celltest(LDD, KsatCR);
           calc("KsatCR *= ksatcal ");
        }
     }

     if (INFIL_METHOD == INFIL_MOREL || INFIL_METHOD == INFIL_GREENAMPT ||
         INFIL_METHOD == INFIL_GREENAMPT2 || INFIL_METHOD == INFIL_SMITH)
     {
        _spatial_input(REAL4, ThetaS1, mapname("thetas1"));
        celltest(LDD, ThetaS1);
        rangetest(ThetaS1,R_GE_LE,0,1,"Soil moisture at saturation");

        _spatial_input(REAL4, ThetaI1, mapname("thetai1"));
        celltest(LDD, ThetaI1);
        rangetest(ThetaI1,R_GE_LE,0,1,"Initial Soil Moisture");

        _spatial (REAL4, THETADIF);
        calc(" THETADIF = ThetaI1-ThetaS1 ");
        rangetest(THETADIF,R_LE, R_DUMMY, 0,"@Initial Soil Moisture cannot be larger than Saturated value,layer 1");
         // THETAS and THETAI as a fraction, e.g. 0.45

        _spatial_input_if(REAL4, PSI1, mapname("psi1"), (bool)(INFIL_METHOD == INFIL_GREENAMPT || INFIL_METHOD == INFIL_GREENAMPT2));
        if((bool)(INFIL_METHOD == INFIL_GREENAMPT || INFIL_METHOD == INFIL_GREENAMPT2))
        {
            celltest(LDD, PSI1);
        }
         // PSI in cm and positive!

        /* not needd for now
        _spatial_input_if(REAL4, GCapDrive, mapname("gcapdrive"), (bool)(INFIL_METHOD == INFIL_MOREL || INFIL_METHOD == INFIL_SMITH));
        */
        _spatial_input_if(REAL4, SoilDepth1, mapname("soildep1"), (bool)(INFIL_METHOD == INFIL_MOREL || INFIL_METHOD == INFIL_SMITH || INFIL_METHOD == INFIL_GREENAMPT || INFIL_METHOD == INFIL_GREENAMPT2 ));
        if((bool)(INFIL_METHOD == INFIL_MOREL || INFIL_METHOD == INFIL_SMITH || INFIL_METHOD == INFIL_GREENAMPT || INFIL_METHOD == INFIL_GREENAMPT2 ))
        {
            celltest(LDD, SoilDepth1);
        }
         // SoilDepth1 in mm

        // depth of wetting front in 1st layer
        _spatial(REAL4, L1);
        calc(" L1 = 1e-6 ");

        _spatial(REAL4, L1Grass);
        calc(" L1Grass = 1e-6 ");

        _spatial(REAL4, L1Wheel);
        calc(" L1Wheel = 1e-6 ");

          //VJ 050812 drainage exponent as in K = Ksat * (thetai/thetas)** DrFactor
        _spatial_input_if(REAL4, DrFactor, mapname("drfactor"), SwitchDrainage);
        if (SwitchDrainage)
        {
           celltest(LDD, DrFactor);
	        rangetest(DrFactor,R_GE_LE,1,10,"Drainage exponent");
        }

     }
     if (INFIL_METHOD == INFIL_GREENAMPT2)
     {
//VJ 060109 bug fix: include this switch here
        SwitchInfilGA2 = true;

        _spatial_input(REAL4, Ksat2, mapname("ksat2"));
        rangetest(Ksat2,R_GT, 0, R_DUMMY,"@Ksat layer 2 must be > 0 mm/h");
        calc("Ksat2 *= ksatcal ");
         // Ksat1 in mm/h

        _spatial_input(REAL4, ThetaS2, mapname("thetas2"));
        rangetest(ThetaS2,R_GE_LE,0,1,"Soil moisture at saturation");

        _spatial_input(REAL4, ThetaI2, mapname("thetai2"));
        rangetest(ThetaI2,R_GE_LE,0,1,"Initial Soil Moisture");

        _spatial (REAL4, THETADIF);
        calc(" THETADIF = ThetaI2-ThetaS2 ");
        rangetest(THETADIF,R_LE, R_DUMMY, 0,"@Initial Soil Moisture cannot be larger than Saturated value,layer 2");
         // THETAS and THETAI as a fraction, e.g. 0.45

        _spatial_input_if(REAL4, PSI2, mapname("psi2"), (bool)(INFIL_METHOD == INFIL_GREENAMPT2));
         // PSI in cm and positive!

        _spatial_input_if(REAL4, SoilDepth2, mapname("soildep2"), (bool)(INFIL_METHOD == INFIL_GREENAMPT2 ));
         // SoilDepth1 in mm

        if ((bool)(INFIL_METHOD == INFIL_GREENAMPT2))
        {
           celltest(LDD, Ksat2);
           celltest(LDD, ThetaS2);
           celltest(LDD, ThetaI2);
           celltest(LDD, PSI2);
           celltest(LDD, SoilDepth2);
        }


        // depth of wetting front in 2nd layer
        _spatial(REAL4, L2);
        calc(" L2 = 0 ");

        _spatial(REAL4, L2Grass);
        calc(" L2Grass = 0 ");

        _spatial(REAL4, L2Wheel);
        calc(" L2Wheel = 0 ");

     }



