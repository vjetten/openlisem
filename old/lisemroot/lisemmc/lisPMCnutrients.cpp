// *****************************************************************************
// ****** nutrient balance in SOLUTION *****************************************
// *****************************************************************************

      // IMPORTANT:
      // IT IS ASSUMED THAT NUTRIENT SOURCES ARE ONLY ON FIELD SURFACES
      // SO WHEEL TRACKS AND CHANNELS DO NOT HAVE NUTRIENTS
      // "FARMER ONLY SPRAYS ON FIELDS !!!"
      // ALL WATER VOLUME PERTAIN TO FIELD SURFACE ONLY (SOILWIDTHDX+STONEWIDTHDX)

      // NutXSolution is the variable that loops (MLXt)
      // new fluxes are added and subtracted here
      // than it is put into the kin wave to be routed:
      // 1) MLXt = MLXt + dt(EXT - MinfiltX)
      //    new estimate of solved nuts with extraction (+) and infiltration (-)
      // 2) MLXt = MLXt + dt(SMLFXin - SMLFXout)
      // 3) MLXt+dt = MLXt
      //    lateral done in the kin wave


// NutX_Q = the amount of nutrient in the soil in kg
// NutXSolution = the amount of nutrient in solution in the runoff water (above ground) in kg
// NutXInfiltration = the nutrient infiltrating from the runoff into the soil in kg
// NutX

//------------------------------------------------------------------------------
// ****** solved nutrient concentrations ***************************************
//------------------------------------------------------------------------------

         // algorithm: CLX = MLX / zl
         // assumed that zl is thickness of water layer * water surface = m3
      calc(" NutWaterVolin = WH*(SoilWidthDX+StoneWidthDX)*DXc/1000 ");
         //NOTE also used before kin wave, this is WH before kin wave

      calc(" NutPConcentration = mif (WH gt 0, NutPSolution/NutWaterVolin, 0) ");
      calc(" NutNH4Concentration = mif (WH gt 0, NutNH4Solution/NutWaterVolin, 0) ");
      calc(" NutNO3Concentration = mif (WH gt 0, NutNO3Solution/NutWaterVolin, 0) ");
         // recalc concentration in kg/m3
         // "in" means before kin wave

//------------------------------------------------------------------------------
// ****** nutrient loss by infiltration ****************************************
//------------------------------------------------------------------------------

         // algorithm: MinfiltX = CLX*q
      calc(" NutPInfiltration = NutPConcentration*WHinf*(SoilWidthDX+StoneWidthDX)*DXc/1000 ");
      calc(" NutNH4Infiltration = NutNH4Concentration*WHinf*(SoilWidthDX+StoneWidthDX)*DXc/1000 ");
      calc(" NutNO3Infiltration = NutNO3Concentration*WHinf*(SoilWidthDX+StoneWidthDX)*DXc/1000 ");
        // infiltrated nutrients in solution: kg/m3 concentration * m3 infiltrated water
        // in this timestep: WHInf/1000*DX*DX = m3/s*dtsec

      calc(" NutPInfiltration = min(NutPInfiltration, NutPSolution)");
      calc(" NutNH4Infiltration = min(NutNH4Infiltration, NutNH4Solution)");
      calc(" NutNO3Infiltration = min(NutNO3Infiltration, NutNO3Solution)");

      calc(" TotalPInfiltration += sum(NutPInfiltration/CatchmentArea)");
      calc(" TotalNH4Infiltration += sum(NutNH4Infiltration/CatchmentArea)");
      calc(" TotalNO3Infiltration += sum(NutNO3Infiltration/CatchmentArea)");
      //total nuts infiltrated in the area

        // not more infiltration then there is in solution
//------------------------------------------------------------------------------
// ****** nutrient exchange in active layer solution ***************************
//------------------------------------------------------------------------------
// these are nutrients that enter the solution from the soil storage

         // algorithm: EXT = Q*epsilon/(theta+gamma*kd)
     calc(" NutPExtRate = mif (NutPSorption gt 0, NutP_Q*NutPEfficiency/(thetatop+BulkDensity*NutPSorption), 0)");
     calc(" NutNH4ExtRate = mif (NutNH4Sorption gt 0, NutNH4_Q*NutNH4Efficiency/(thetatop+BulkDensity*NutNH4Sorption), 0)");
     calc(" NutNO3ExtRate = mif (NutNO3Sorption gt 0, NutNO3_Q*NutNO3Efficiency/(thetatop+BulkDensity*NutNO3Sorption), 0)");
      //in paramnutrients: calc(" NutP_Q = NutP_Qinit*activesurface ");
      //Q in kg/m2 * m2 = kg
      // Extraction rate in
      // kg*1/s / (-+kg/m3*m3/kg) = kg/s

        //EXT: extraction rate in kg/s
     calc(" NutPExtRate = min(NutPExtRate, NutP_Q/DTSEC) ");
     calc(" NutNH4ExtRate = min(NutNH4ExtRate, NutNH4_Q/DTSEC) ");
     calc(" NutNO3ExtRate = min(NutNO3ExtRate, NutNO3_Q/DTSEC) ");

        // algorithm: Q(t+dt) = Q(t) - dt*EXT
        //no more extraction than there is present
     calc(" NutP_Q = max(0, NutP_Q - DTSEC*NutPExtRate)");
     calc(" NutNH4_Q = max(0, NutNH4_Q - DTSEC*NutNH4ExtRate)");
     calc(" NutNO3_Q = max(0, NutNO3_Q - DTSEC*NutNO3ExtRate)");
       // adjust solute store in active layer
       // units: kg = kg - s*kg/s


//------------------------------------------------------------------------------
// ****** adjust nutrients in surface water in solution ************************
//------------------------------------------------------------------------------

       // algorithm: MLX(t+dt) = MLX(t) + dt(MLFxin-MLFXout + EXT - MinfiltX)
     calc(" NutPSolution = max(0, NutPSolution - NutPInfiltration + NutPExtRate*DTSEC) ");
     calc(" NutNH4Solution = max(0, NutNH4Solution - NutNH4Infiltration + NutNH4ExtRate*DTSEC) ");
     calc(" NutNO3Solution = max(0, NutNO3Solution - NutNO3Infiltration + NutNO3ExtRate*DTSEC) ");
       // nutrients in solution before kinematic wave in kg
       // is recalculated after the "ninematic"



//------------------------------------------------------------------------------
// ****** nutrients in surface water in suspension *****************************
//------------------------------------------------------------------------------

     //SUSPENSION CALCUALTED AFTER THE MULTICLASS KINEMATIC WAVE
