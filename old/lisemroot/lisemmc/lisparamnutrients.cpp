   //---------------------------------------------------------------------------
   //********************* Nutrients *******************************************
   //---------------------------------------------------------------------------

   // INPUT
   _spatial_input_if(REAL4, NutPContent, mapname("pcont"), SwitchNutrients); //XSOIL
   _spatial_input_if(REAL4, NutP_Qinit, mapname("psolute"), SwitchNutrients); //Q
   _spatial_input_if(REAL4, NutPEfficiency, mapname("pefficiency"), SwitchNutrients); //epsilon
   _spatial_input_if(REAL4, NutPSorption, mapname("psorp"), SwitchNutrients);  // kd
   _spatial_input_if(REAL4, NutPConversion, mapname("pconv"), SwitchNutrients); //ER

   _spatial_input_if(REAL4, NutNO3Content, mapname("no3cont"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNO3_Qinit, mapname("no3solute"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNO3Efficiency, mapname("no3efficiency"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNO3Sorption, mapname("no3sorp"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNO3Conversion, mapname("no3conv"), SwitchNutrients);

   _spatial_input_if(REAL4, NutNH4Content, mapname("nh4cont"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNH4_Qinit, mapname("nh4solute"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNH4Efficiency, mapname("nh4efficiency"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNH4Sorption, mapname("nh4sorp"), SwitchNutrients);
   _spatial_input_if(REAL4, NutNH4Conversion, mapname("nh4conv"), SwitchNutrients);

   _spatial_input_if(REAL4, BulkDensity, mapname("bulk"), SwitchNutrients);
   if (SwitchNutrients)
   {
      celltest(LDD, NutPContent      );
      celltest(LDD, NutP_Qinit       );
      celltest(LDD, NutPEfficiency   );
      celltest(LDD, NutPSorption     );
      celltest(LDD, NutPConversion   );
      celltest(LDD, NutNO3Content    );
      celltest(LDD, NutNO3_Qinit     );
      celltest(LDD, NutNO3Efficiency );
      celltest(LDD, NutNO3Sorption   );
      celltest(LDD, NutNO3Conversion );
      celltest(LDD, NutNH4Content    );
      celltest(LDD, NutNH4_Qinit     );
      celltest(LDD, NutNH4Efficiency );
      celltest(LDD, NutNH4Sorption   );
      celltest(LDD, NutNH4Conversion );
   }

   _spatial(REAL4, thetatop);
   // moisture content surface layers (calculated in SWATRE)

   _spatial_if(REAL4, NutPSolution, SwitchNutrients);
   _spatial_if(REAL4, NutNH4Solution, SwitchNutrients);
   _spatial_if(REAL4, NutNO3Solution, SwitchNutrients);

   _spatial_if(REAL4, NutPSuspension, SwitchNutrients);
   _spatial_if(REAL4, NutNH4Suspension, SwitchNutrients);
   _spatial_if(REAL4, NutNO3Suspension, SwitchNutrients);

   _spatial_if(REAL4, NutPInfiltration, SwitchNutrients);
   _spatial_if(REAL4, NutNO3Infiltration, SwitchNutrients);
   _spatial_if(REAL4, NutNH4Infiltration, SwitchNutrients);
   // NUT solution and suspension and infiltration maps (kg)

   _spatial_if(REAL4, ChanNutPSolution, SwitchNutrients);
   _spatial_if(REAL4, ChanNutNH4Solution, SwitchNutrients);
   _spatial_if(REAL4, ChanNutNO3Solution, SwitchNutrients);
   _spatial_if(REAL4, ChanNutPSuspension, SwitchNutrients);
   _spatial_if(REAL4, ChanNutNH4Suspension, SwitchNutrients);
   _spatial_if(REAL4, ChanNutNO3Suspension, SwitchNutrients);
//VJ 030414 added channel solution and suspension (kg)

   _spatial_if(REAL4, NutWaterVolin, SwitchNutrients);
   // water volume excluding roads and channels

   //_spatial(REAL4, NutQout);
   //calc("NutQout = 0");
   // Qout of kin wave unsed in ninematic.... not used anymore

   _nonspatial(REAL4, NutPOutflow);
   _nonspatial(REAL4, NutNH4Outflow);
   _nonspatial(REAL4, NutNO3Outflow);
   _nonspatial(REAL4, NutPOutflows);
   _nonspatial(REAL4, NutNH4Outflows);
   _nonspatial(REAL4, NutNO3Outflows);
   NutPOutflow = 0;
   NutNH4Outflow = 0;
   NutNO3Outflow = 0;
   NutPOutflows = 0;
   NutNH4Outflows = 0;
   NutNO3Outflows = 0;
     // discharge of nutrients from the catchments in kg/s

   _nonspatial(REAL4, ChanNutPOutflow);
   _nonspatial(REAL4, ChanNutNH4Outflow);
   _nonspatial(REAL4, ChanNutNO3Outflow);
   _nonspatial(REAL4, ChanNutPOutflows);
   _nonspatial(REAL4, ChanNutNH4Outflows);
   _nonspatial(REAL4, ChanNutNO3Outflows);
   ChanNutPOutflow = 0;
   ChanNutNH4Outflow = 0;
   ChanNutNO3Outflow = 0;
   ChanNutPOutflows = 0;
   ChanNutNH4Outflows = 0;
   ChanNutNO3Outflows = 0;
     // channel discharge of nutrients from the catchments in kg/s

//VJ 030627 fixed typo xxxPdep = xxxPDep in all the code
 	_spatial_if(REAL4, NutPDep  , SwitchNutrients);
	_spatial_if(REAL4, NutNH4Dep, SwitchNutrients);
	_spatial_if(REAL4, NutNO3Dep, SwitchNutrients);
	_spatial_if(REAL4, NutPDet  , SwitchNutrients);
	_spatial_if(REAL4, NutNH4Det, SwitchNutrients);
	_spatial_if(REAL4, NutNO3Det, SwitchNutrients);
   // total nutrients in suspension eroded or deposited kg
   //VJ 030415 added here
   // initialisations
   if (SwitchNutrients)
   {
      _spatial(REAL4, NutPExtRate);
      calc(" NutPExtRate = 0 ");
      _spatial(REAL4, NutNH4ExtRate);
      calc(" NutNH4ExtRate = 0 ");
      _spatial(REAL4, NutNO3ExtRate);
      calc(" NutNO3ExtRate = 0 ");
        // nutrients extracted from soil in kg

      calc(" NutWaterVolin = 0");
        // volume used to calculate concentrations

      calc(" NutPSuspension = 0 ");
      calc(" NutNH4Suspension = 0 ");
      calc(" NutNO3Suspension = 0 ");
        // nutrients in suspension in surface water in kg
        // related to Mu0 fractionof

      calc(" NutPSolution = 0 ");
      calc(" NutNH4Solution = 0 ");
      calc(" NutNO3Solution = 0 ");
        // nutrients in solution in surface water in kg
        // "in" means before kinematic wave

      calc(" NutPInfiltration = 0 ");
      calc(" NutNO3Infiltration = 0 ");
      calc(" NutNH4Infiltration = 0 ");
        // MinfiltX: nutrients infiltrated in kg

      calc(" ChanNutPSuspension = 0 ");
      calc(" ChanNutNH4Suspension = 0 ");
      calc(" ChanNutNO3Suspension = 0 ");
      calc(" ChanNutPSolution = 0 ");
      calc(" ChanNutNH4Solution = 0 ");
      calc(" ChanNutNO3Solution = 0 ");
        // idem channels (kg)
//VJ 030414 added channel solution and suspension

		calc(" NutPDep   = 0 ");
		calc(" NutNH4Dep = 0 ");
		calc(" NutNO3Dep = 0 ");
		calc(" NutPDet   = 0 ");
		calc(" NutNH4Det = 0 ");
		calc(" NutNO3Det = 0 ");
		//total eros/depo kg
//VJ 030415 added total eros/depo suspension kg
		
      _spatial(REAL4, NutPConcentration);
      calc(" NutPConcentration = 0 ");
      _spatial(REAL4, NutNH4Concentration);
      calc(" NutNH4Concentration = 0 ");
      _spatial(REAL4, NutNO3Concentration);
      calc(" NutNO3Concentration = 0 ");
        // CLX: nutrients in solution in surface water in kg/m3

      _spatial(REAL4, activesurface);
      if (SwitchIncludeChannel)
       calc(" activesurface = DXc*(DX-ChannelWidthDX-RoadWidthDX)");
      else
       calc(" activesurface = DXc*(DX-RoadWidthDX)");

       // Q: solute store in soil water per cell
       // in kg/m2 of surface, no source from wheeltracks, channels and roads
      _spatial(REAL4, NutP_Q);
      calc(" NutP_Q = NutP_Qinit*activesurface ");
      _spatial(REAL4, NutNH4_Q);
      calc(" NutNH4_Q = NutNH4_Qinit*activesurface ");
      _spatial(REAL4, NutNO3_Q);
      calc(" NutNO3_Q = NutNO3_Qinit*activesurface ");
      //Q in kg = kg/m2 * m2

        // soil moisture average of layers 1-3
      calc(" thetatop = 0 ");
 }

