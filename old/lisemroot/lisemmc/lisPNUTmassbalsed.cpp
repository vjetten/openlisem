//======================================================================
// ********************** NUTRIENT TOTALS ******************************
//======================================================================
//VJ 030415 redone completely

		 // ***** NON SPATIAL *****
		 
		if (SwitchIncludeChannel)
		{
			NutPOutflow    += ChanNutPOutflow   ;
			NutNH4Outflow  += ChanNutNH4Outflow ;
			NutNO3Outflow  += ChanNutNO3Outflow ;
			NutPOutflows   += ChanNutPOutflows  ;
			NutNH4Outflows += ChanNutNH4Outflows;
			NutNO3Outflows += ChanNutNO3Outflows;                         
		}
			 
		TotalPSolution     += NutPOutflow   *DTSEC;
		TotalNH4Solution   += NutNH4Outflow *DTSEC;
		TotalNO3Solution   += NutNO3Outflow *DTSEC;
		TotalPSuspension   += NutPOutflows  *DTSEC;
		TotalNH4Suspension += NutNH4Outflows*DTSEC;
		TotalNO3Suspension += NutNO3Outflows*DTSEC;
		// total nutrients at outlet in solution and suspension (kg)
				
		calc(" TotalPInfiltration   += sum(NutPInfiltration  ) ");
		calc(" TotalNH4Infiltration += sum(NutNO3Infiltration) ");
		calc(" TotalNO3Infiltration += sum(NutNH4Infiltration) ");
		// total nutrients lost through infiltratrion (kg)
		 

		calc(" TotalPDeposition   += sum(NutPDep)  /CatchmentArea");
		calc(" TotalNH4Deposition += sum(NutNH4Dep)/CatchmentArea");
		calc(" TotalNO3Deposition += sum(NutNO3Dep)/CatchmentArea");
		
		calc(" TotalPDetachment   += sum(NutPDet)  /CatchmentArea");
		calc(" TotalNH4Detachment += sum(NutNH4Det)/CatchmentArea");
		calc(" TotalNO3Detachment += sum(NutNO3Det)/CatchmentArea");
		// total erosion or deposition of suspended material, average in kg/m2
		

			// ***** SPATIAL *****
		 
		calc("SumTotNutPdep   += NutPDep  /(DX*DX) ");
		calc("SumTotNutNH4dep += NutNH4Dep/(DX*DX) ");
		calc("SumTotNutNO3dep += NutNO3Dep/(DX*DX) ");
		calc("SumTotNutPdet   += NutPDet  /(DX*DX) ");
		calc("SumTotNutNH4det += NutNH4Det/(DX*DX) ");
		calc("SumTotNutNO3det += NutNO3Det/(DX*DX) ");
		// calculate totals for the run for output kg per m2
