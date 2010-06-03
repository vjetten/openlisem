/**************************************************/
/* GREEN AND AMPT MODULE FOR GRASSSTRIPS          */
/* all vars in mm                                 */
/**************************************************/

//NOTE regular infil GA is done before this

      _spatial(REAL4, DeltaF);
      calc(" DeltaF = DTHOUR*KsatGR*(1+(PSI1*10+WHGrass)/(L1Grass+0.001)) ");
         // infiltration according to G&A solution of Darcy, kutilek en Nielsen, p 138

      calc(" InfilPotGrass = -DeltaF ");
           //potential infil needed for surplus kin wave

      if (!SwitchInfilGA2)
      {
         // not more than water on surface
         // if no second layer than do check with WH now else at the end of second layer
         calc(" DeltaF = min(WHGrass, DeltaF) ");
      }

      _spatial(REAL4, dL1);
      calc(" dL1 = mif(ThetaS1 gt ThetaI1, DeltaF/(ThetaS1-ThetaI1), SoilDepth1-L1Grass) ");
          // increase of wetting front 1st layer
      calc(" dL1 = max(0, dL1) ");
      //VJ 120805 if not impermeable and soildepth1 limited, dL1 could become negative
      //now set to 0 if negative which means that the infiltration becomes constant (except for WH)
      //soildepth1 is a control variable


      if (!SwitchInfilGA2)
      {
         if (SwitchImpermeable)
         {
            //cut off dL1
            calc(" dL1 = min(dL1, SoilDepth1-L1Grass) ");
            //VJ 050812 correct DeltaF in case of impermeable
            calc(" DeltaF = min(SoilDepth1 - L1Grass, DeltaF ");
         }

         calc(" L1Grass += dL1 ");
         // depth of wetting front

         //VJ 050812 Adjust moisture content if not impermeable
         if (SwitchDrainage)
         {
              _spatial(REAL4, Kd);
              calc(" Kd = KsatGR*(ThetaI1/ThetaS1)**DrFactor ");
              calc(" ThetaI1 -= mif(ThetaI1 > 0.02*ThetaS1, (Kd*DTHOUR)/SoilDepth1, 0) ");
         }
      }


      //********* SECOND LAYER GA ******************

      if (SwitchInfilGA2)
      {
         _spatial(UINT2, check);
         calc(" check = mif(L2 eq 0 and dL1+L1Grass gt SoilDepth1, 1, 0) ");
         calc(" check = mif(L2 gt 0, 2, check) ");
            //check: 0 first layer, 1: water entering 2nd; 2: water in second

         _spatial(REAL4, dL2);
         calc(" dL2 = mif(check eq 1, SoilDepth1-(L1Grass+dL1)*(ThetaS1-ThetaI1)/(ThetaS2-ThetaI2), 0) ");
           // if water just enters 2nd layer make a guess for L2

//VJ 050812 recalc dL2 in case water just enters layer 2
         calc(" dL2 = 0 ");
         // set dL2 to 0 (e.g. check is 0)
         calc(" L1Grass += dL1 ");
         // increase water layer in first
         calc(" dL2 = mif(check eq 1, L1Grass - SoilDepth1, dL2 ");
         // dL2 is overshoot of L1
         calc(" dL2 = mif(check eq 1 and ThetaS2 gt ThetaI2, dL2 * (ThetaS1-ThetaI1)/(ThetaS2-ThetaI2), 1) ");
         // correct dL2 for porosity in second layer
         calc(" L1Grass = min(L1Grass, SoilDepth1) ");
           // limit L1 to Soildepth1

         calc(" L2 = mif(check eq 1, dL2, L2) ");
           // if water just enters 2nd layer make a guess for L2

         _spatial(REAL4, Reff);
         calc(" Reff = mif(KsatGR eq 0 and Ksat2 eq 0, 1e20, L1Grass/KsatGR + L2/Ksat2 ");
           // effective hydraulic conductivity as resistance

         calc(" DeltaF = mif(check gt 0, DTHOUR*((PSI2*10+WHGrass)+L1Grass+L2)/Reff), DeltaF) ");
           // if L2 gt 0, change deltaf else leave it at 1st layer deltaf

         calc(" L2 = mif(check eq 1, 0, L2) ");
           // put L2 back to add the real dL2

         calc(" DeltaF = min(WHGrass, DeltaF) ");
            // not more than water on surface

         calc(" dL2 = mif(ThetaS2 gt ThetaI2, DeltaF/(ThetaS2-ThetaI2), SoilDepth2-L2) ");
            // calc real dL2
         calc(" dL2 = max(dL2, 0) ");
         //VJ 050812 cannot be negative, can happen when soildepth2 is limited
         calc(" dL2 = mif(check gt 0, dL2, 0) ");
			//if not in layer 2 then set back to zero, just to be sure
			
         if (SwitchImpermeable)
         {
            calc(" dL2 = min(dL2, SoilDepth2-L2) ");
            //VJ 050812 correct DeltaF in case of impermeable
            calc(" DeltaF = min(SoilDepth2 - L2, DeltaF ");
            // adjust DeltaF if this is the case
         }

         calc(" L2 += dL2 ");
         // depth of wetting front
         //VJ 050812 Adjust moisture content if not impermeable

         if (SwitchDrainage)
         {
              _spatial(REAL4, Kd);
              calc(" Kd = Ksat2*(ThetaI2/ThetaS2)**DrFactor ");
              calc(" ThetaI2 -= mif(ThetaI2 > 0.02*ThetaS2, (Kd*DTHOUR)/SoilDepth2, 0) ");
         }
      }

      calc(" WHGrass -= DeltaF ");
        // adjust WH with new infil

