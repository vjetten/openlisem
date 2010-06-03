/**************************************************/
/* GREEN AND AMPT MODULE FOR GRASSSTRIP          */
/* all vars in mm                                 */
/* VJ 060109: number of bug fixes (brackets) version 2.392
/**************************************************/

//NOTE regular infil GA is done before this

      calc( " L1Grass = mif( ThetaI1 ge ThetaS1, SoilDepth1, L1Grass )");
        //VJ 051006 check for saturated soil

      _spatial(REAL4, DeltaF);
      calc(" DeltaF = DTHOUR*KsatGR*(1+(PSI1*10+WHGrass)/(L1Grass+0.001)) ");
         // infiltration according to G&A solution of Darcy, kutilek en Nielsen, p 138

      calc(" InfilPotGrass = -DeltaF ");
           //potential infil needed for surplus kin wave

		calc(" DeltaF = min(WHGrass, DeltaF) ");
			// no more infil than water on surface

		_spatial(REAL4, dL1);
		calc(" dL1 = mif(ThetaI1 lt ThetaS1, DeltaF/(ThetaS1-ThetaI1), 0) ");
			 // increase of wetting front 1st layer

      if (!SwitchInfilGA2)
      {
        if (SwitchImpermeable)
         {
            //cut off dL1
            calc(" dL1 = min(dL1, SoilDepth1-L1Grass) ");
            //VJ 050812 correct DeltaF in case of impermeable
            calc(" DeltaF = min(SoilDepth1 - L1Grass, DeltaF ");

//VJ 050829 fixed, infilpot should become 0 when impermeable
				calc(" InfilPotGrass = -DeltaF ");
         }

         calc(" L1Grass += dL1 ");
         // depth of wetting front

         //VJ 050812 Adjust moisture content if not impermeable
         if (SwitchDrainage)
         {
				//VJ 050829 changed drainage: decreases thetai or decreases L1
				  _spatial(REAL4, Kd);
				  calc(" Kd = mif(L1Grass ge SoilDepth1, KsatGR, KsatGR*(ThetaI1/ThetaS1)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = KsatGR

				  calc(" ThetaI1 -= mif(L1Grass lt SoilDepth1, (Kd*DTHOUR)/(SoilDepth1-L1Grass),  0) ");
				  // adjust thetai if soil is not full
				  calc(" ThetaI1 = min(max(ThetaI1, 0.02*ThetaS1), ThetaS1-0.01) ");
				  // limit thetai between 0.02*pore and pore-0.01

				  calc(" L1Grass -= mif(L1Grass ge SoilDepth1, Kd*DTHOUR, 0) ");
				  // if profile full decrease L1 with drainage
         }
      }

      //********* SECOND LAYER GA ******************

      if (SwitchInfilGA2)
      {
         calc( " L2Grass = mif( ThetaI2 lt ThetaS2, SoilDepth2, L2Grass) ");

         _spatial(UINT1, check);
         calc(" check = mif(L2Grass eq 0 and dL1+L1Grass gt SoilDepth1, 1, 0) ");
         calc(" check = mif(L2Grass gt 0, 2, check) ");
            //check: 0 first layer, 1: water entering 2nd; 2: water in second

//VJ 050812 recalc dL2 in case water just enters layer 2
         _spatial(REAL4, dL2);
			calc(" dL2 = 0 ");
			// set dL2 to 0 (e.g. check is 0)

         calc(" L1Grass += dL1 ");
         // increase water layer in first

         calc(" dL2 = mif(check eq 1 and ThetaI2 lt ThetaS2, (L1Grass - SoilDepth1)*(ThetaS1-ThetaI1)/(ThetaS2-ThetaI2), dL2) ");
			// dL2 is overshoot of L1, corrected for porosity in second layer, S always > I

         calc(" L1Grass = min(L1Grass, SoilDepth1) ");
           // limit L1 to Soildepth1

         calc(" L2Grass = mif(check eq 1, dL2, L2Grass) ");
           // if water just enters 2nd layer make a guess for L2

         _spatial(REAL4, Reff);
         calc(" Reff = mif(KsatGR eq 0 or Ksat2 eq 0, 1e20, L1Grass/KsatGR + L2Grass/Ksat2) ");
           // effective hydraulic conductivity as resistance
			  // if both zero assume a very high resistance

         calc(" DeltaF = mif(check gt 0, DTHOUR*((PSI2*10+WHGrass)+L1Grass+L2Grass)/Reff), DeltaF) ");
           // if L2 gt 0, change deltaf else leave it at 1st layer deltaf

			calc(" InfilPotGrass = -DeltaF ");
			  //adjust potential infil for kin wave when 2nd layer

         calc(" L2Grass = mif(check eq 1, 0, L2Grass) ");
           // put L2 back to add the real dL2

         calc(" DeltaF = min(WHGrass, DeltaF) ");
            // not more than water on surface

         calc(" dL2 = mif(ThetaI2 lt ThetaS2, DeltaF/(ThetaS2-ThetaI2), 0) ");
            // calc real dL2
         calc(" dL2 = mif(check gt 0, dL2, 0) ");
			//if not in layer 2 then set back to zero, just to be sure
			//in that case DeltaF only applies to L1

         if (SwitchImpermeable)
         {
            calc(" dL2 = min(dL2, SoilDepth2-L2Grass) ");
            //VJ 050812 correct DeltaF in case of impermeable
            calc(" DeltaF = min(SoilDepth2 - L2Grass, DeltaF) ");
            // adjust DeltaF if this is the case
//VJ 050829 fixed, infilpot should become 0 when impermeable
				calc(" InfilPotGrass = -DeltaF ");
         }

         calc(" L2Grass += dL2 ");
         // depth of wetting front
         //VJ 050812 Adjust moisture content if not impermeable

         if (SwitchDrainage)
         {
              _spatial(REAL4, Kd);
				  calc(" Kd = mif(L2Grass ge SoilDepth2, Ksat2, Ksat2*(ThetaI2/ThetaS2)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = ksat

				  calc(" ThetaI2 -= mif(L2Grass lt SoilDepth2, (Kd*DTHOUR)/(SoilDepth2-L2Grass),  0) ");
				  // adjust thetai if soil is not full
				  calc(" ThetaI2 = min(max(ThetaI2, 0.02*ThetaS2), ThetaS2-0.01) ");
				  // limit thetai between 0.02*pore and pore-0.01

				  calc(" L2Grass -= mif(L2Grass ge SoilDepth2, Kd*DTHOUR, 0) ");
				  // if profile full decrease L with drainage
         }
      }

      calc(" WHGrass -= DeltaF ");
        // adjust WH with new infil

