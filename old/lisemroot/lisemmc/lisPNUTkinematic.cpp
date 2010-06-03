// *****************************************************************************
// **** kinematic wave nutrients in solution and suspension of overland flow ***
// *****************************************************************************

         //NOTE:
         //CALL AFTER KIN WAVE MULTICLASS BECAUSE ALPHA, QIN ETC ARE CALCULATED THERE


		 // algorithm: MLFXin and   MLFXout in kg/m3/s
		 // note the volume of water on  the gridcell is taken so MLFX is in kg/s
		 // for  the amount of water m3 on the surface of the gridcell

		 // make a  flux of the solved matter in kg/s
		_spatial(REAL4, NutPQin);
		_spatial(REAL4, NutNH4Qin);
		_spatial(REAL4, NutNO3Qin);
		_spatial(REAL4, NutPQins);
		_spatial(REAL4, NutNH4Qins);
		_spatial(REAL4, NutNO3Qins);

		_spatial(REAL4, qfactor);
		calc("qfactor = mif(NutWaterVolin gt 0.00001, Qin/NutWaterVolin, 0) ");
		calc("NutPQin    = NutPSolution     * qfactor ");
		calc("NutNH4Qin  = NutNH4Solution   * qfactor ");
		calc("NutNO3Qin  = NutNO3Solution   * qfactor ");
		calc("NutPQins   = NutPSuspension   * qfactor ");
		calc("NutNH4Qins = NutNH4Suspension * qfactor ");
		calc("NutNO3Qins = NutNO3Suspension * qfactor ");
		// nutrient fluxes between cells in m3/s  * kg/m3 = kg/s
		// NutWaterVolin BASED ON WH BEFORE KIN WAVE

		_nonspatial(REAL4, SumNutP);
		_nonspatial(REAL4, SumNutNH4);
		_nonspatial(REAL4, SumNutNO3);
		_nonspatial(REAL4, SumNutPs);
		_nonspatial(REAL4, SumNutNH4s);
		_nonspatial(REAL4, SumNutNO3s);
		calc(" SumNutP = sum(NutPSolution) ");
		calc(" SumNutNH4 = sum(NutNH4Solution) ");
		calc(" SumNutNO3 = sum(NutNO3Solution) ");
		calc(" SumNutPs = sum(NutPSuspension) ");
		calc(" SumNutNH4s = sum(NutNH4Suspension) ");
		calc(" SumNutNO3s = sum(NutNO3Suspension) ");
		// sums  needed for mass balance

		_spatial(REAL4, NutPQout);
		_spatial(REAL4, NutNH4Qout);
		_spatial(REAL4, NutNO3Qout);
		_spatial(REAL4, NutPQouts);
		_spatial(REAL4, NutNH4Qouts);
		_spatial(REAL4, NutNO3Qouts);
		_spatial(REAL4, Qout);

		_spatial(REAL4, infil);
		if (SwitchKinwaveInfil)
     //VJ 050704 changed DXc to DX
		calc(" infil = InfilSurplus/(1000*DTSEC)*DX ");
		else
		calc(" infil = 0 ");

		// fluxes after kin  wave
		minematic(Qout, Qin, infil,
		  NutPQout, NutNH4Qout, NutNO3Qout, NutPQouts, NutNH4Qouts, NutNO3Qouts,
		  NutPQin, NutNH4Qin, NutNO3Qin, NutPQins, NutNH4Qins, NutNO3Qins,
		  LDD, Alpha, Beta, DTSEC, DXc);
//NOTE infil is just for the correct kin wave here, it is already accounted for in MCkinwave
          // fluxes after kin wave
          
		calc(" NutPOutflow = sum(mif(Outlet1 eq 1, NutPQout, 0)) ");
		calc(" NutNH4Outflow = sum(mif(Outlet1 eq 1, NutNH4Qout, 0)) ");
		calc(" NutNO3Outflow = sum(mif(Outlet1 eq 1, NutNO3Qout, 0)) ");
		calc(" NutPOutflows = sum(mif(Outlet1 eq 1, NutPQouts, 0)) ");
		calc(" NutNH4Outflows = sum(mif(Outlet1 eq 1, NutNH4Qouts, 0)) ");
		calc(" NutNO3Outflows = sum(mif(Outlet1 eq 1, NutNO3Qouts, 0)) ");


		calc(" NutWaterVolin = WH*(SoilWidthDX+StoneWidthDX)*DXc/1000 ");
		// recalc volume  in cell because after kinematic wave for water
		//THIS IS WH AFTER THE MC KIN WAVE WATER

		calc("qfactor = mif(Qout gt 0.00001, NutWaterVolin/Qout,0)"); 
		calc("NutPSolution = NutPQout * qfactor "); 
		calc("NutNH4Solution = NutNH4Qout * qfactor ");
		calc("NutNO3Solution = NutNO3Qout *	qfactor ");
		calc("NutPSuspension = NutPQouts * qfactor ");
		calc("NutNH4Suspension = NutNH4Qouts * qfactor ");
		calc("NutNO3Suspension = NutNO3Qouts * qfactor ");
		// new nutrients in solution  after kin wave in kg

		if (SwitchCorrectMassSED) // default true !
		{
			_nonspatial(REAL4, SumNutPout);
			_nonspatial(REAL4, SumNutNH4out);
			_nonspatial(REAL4, SumNutNO3out);
			_nonspatial(REAL4, SumNutPouts);
			_nonspatial(REAL4, SumNutNH4outs);
			_nonspatial(REAL4, SumNutNO3outs);
			calc(" SumNutPout = sum(NutPSolution) ");
			calc(" SumNutNH4out = sum(NutNH4Solution) ");
			calc(" SumNutNO3out = sum(NutNO3Solution) ");
			calc(" SumNutPouts = sum(NutPSuspension) ");
			calc(" SumNutNH4outs = sum(NutNH4Suspension) ");
			calc(" SumNutNO3outs = sum(NutNO3Suspension) ");

			if (fabs(SumNutPout) > 1e-6)
				calc(" NutPSolution += NutPSolution*(SumNutP-SumNutPout-NutPOutflow*DTSEC)/SumNutPout ");
			if (fabs(SumNutNH4out) > 1e-6)
				calc(" NutNH4Solution += NutNH4Solution*(SumNutNH4-SumNutNH4out-NutNH4Outflow*DTSEC)/SumNutNH4out ");
			if (fabs(SumNutNO3out) > 1e-6)
				calc(" NutNO3Solution += NutNO3Solution*(SumNutNO3-SumNutNO3out-NutNO3Outflow*DTSEC)/SumNutNO3out ");
			if (fabs(SumNutPouts) > 1e-6)
				calc(" NutPSuspension += NutPSuspension*(SumNutPs-SumNutPouts-NutPOutflows*DTSEC)/SumNutPouts ");
			if (fabs(SumNutNH4outs) > 1e-6)
				calc(" NutNH4Suspension += NutNH4Suspension*(SumNutNH4s-SumNutNH4outs-NutNH4Outflows*DTSEC)/SumNutNH4outs ");
			if (fabs(SumNutNO3outs) > 1e-6)
				calc(" NutNO3Suspension += NutNO3Suspension*(SumNutNO3s-SumNutNO3outs-NutNO3Outflows*DTSEC)/SumNutNO3outs ");
				// correct mass balance
		}


