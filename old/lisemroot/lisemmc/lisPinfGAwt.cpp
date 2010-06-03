/**************************************************/
/* GREEN AND AMPT MODULE FOR WheelSTRIP          */
/* all vars in mm                                 */
/* VJ 060109: number of bug fixes (brackets) version 2.392
/**************************************************/

//NOTE regular infil GA is done before this

      calc( " L1Wheel = mif( ThetaI1 ge ThetaS1, SoilDepth1, L1Wheel )");
        //VJ 051006 check for saturated soil

      _spatial(REAL4, DeltaF);
      calc(" DeltaF = DTHOUR*KsatGR*(1+(PSI1*10+WHWheel)/(L1Wheel+0.001)) ");
         // infiltration according to G&A solution of Darcy, kutilek en Nielsen, p 138

      calc(" InfilPotWheel = -DeltaF ");
           //potential infil needed for surplus kin wave

		calc(" DeltaF = min(WHWheel, DeltaF) ");
			// no more infil than water on surface

		_spatial(REAL4, dL1);
		calc(" dL1 = mif(ThetaI1 lt ThetaS1, DeltaF/(ThetaS1-ThetaI1), 0) ");
			 // increase of wetting front 1st layer

      if (!SwitchInfilGA2)
      {
        if (SwitchImpermeable)
         {
            //cut off dL1
            calc(" dL1 = min(dL1, SoilDepth1-L1Wheel) ");
            //VJ 050812 correct DeltaF in case of impermeable
            calc(" DeltaF = min(SoilDepth1 - L1Wheel, DeltaF ");

//VJ 050829 fixed, infilpot should become 0 when impermeable
				calc(" InfilPotWheel = -DeltaF ");
         }

         calc(" L1Wheel += dL1 ");
         // depth of wetting front

         //VJ 050812 Adjust moisture content if not impermeable
         if (SwitchDrainage)
         {
				//VJ 050829 changed drainage: decreases thetai or decreases L1
				  _spatial(REAL4, Kd);
				  calc(" Kd = mif(L1Wheel ge SoilDepth1, KsatGR, KsatGR*(ThetaI1/ThetaS1)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = KsatGR

				  calc(" ThetaI1 -= mif(L1Wheel lt SoilDepth1, (Kd*DTHOUR)/(SoilDepth1-L1Wheel),  0) ");
				  // adjust thetai if soil is not full
				  calc(" ThetaI1 = min(max(ThetaI1, 0.02*ThetaS1), ThetaS1-0.01) ");
				  // limit thetai between 0.02*pore and pore-0.01

				  calc(" L1Wheel -= mif(L1Wheel ge SoilDepth1, Kd*DTHOUR, 0) ");
				  // if profile full decrease L1 with drainage
         }
      }

      //********* SECOND LAYER GA ******************

      if (SwitchInfilGA2)
      {
         calc( " L2Wheel = mif( ThetaI2 lt ThetaS2, SoilDepth2, L2Wheel )");

         _spatial(UINT1, check);
         calc(" check = mif(L2Wheel eq 0 and dL1+L1Wheel gt SoilDepth1, 1, 0) ");
         calc(" check = mif(L2Wheel gt 0, 2, check) ");
            //check: 0 first layer, 1: water entering 2nd; 2: water in second

//VJ 050812 recalc dL2 in case water just enters layer 2
         _spatial(REAL4, dL2);
			calc(" dL2 = 0 ");
			// set dL2 to 0 (e.g. check is 0)

         calc(" L1Wheel += dL1 ");
         // increase water layer in first

         calc(" dL2 = mif(check eq 1 and ThetaI2 lt ThetaS2, (L1Wheel - SoilDepth1)*(ThetaS1-ThetaI1)/(ThetaS2-ThetaI2), dL2) ");
			// dL2 is overshoot of L1, corrected for porosity in second layer, S always > I

         calc(" L1Wheel = min(L1Wheel, SoilDepth1) ");
           // limit L1 to Soildepth1

         calc(" L2Wheel = mif(check eq 1, dL2, L2Wheel) ");
           // if water just enters 2nd layer make a guess for L2

         _spatial(REAL4, Reff);
         calc(" Reff = mif(KsatGR eq 0 or Ksat2 eq 0, 1e20, L1Wheel/KsatGR + L2Wheel/Ksat2) ");
           // effective hydraulic conductivity as resistance
			  // if both zero assume a very high resistance

         calc(" DeltaF = mif(check gt 0, DTHOUR*((PSI2*10+WHWheel)+L1Wheel+L2Wheel)/Reff), DeltaF) ");
           // if L2 gt 0, change deltaf else leave it at 1st layer deltaf

			calc(" InfilPotWheel = -DeltaF ");
			  //adjust potential infil for kin wave when 2nd layer

         calc(" L2Wheel = mif(check eq 1, 0, L2Wheel) ");
           // put L2 back to add the real dL2

         calc(" DeltaF = min(WHWheel, DeltaF) ");
            // not more than water on surface

         calc(" dL2 = mif(ThetaI2 lt ThetaS2, DeltaF/(ThetaS2-ThetaI2), 0) ");
            // calc real dL2
         calc(" dL2 = mif(check gt 0, dL2, 0) ");
			//if not in layer 2 then set back to zero, just to be sure
			//in that case DeltaF only applies to L1

         if (SwitchImpermeable)
         {
            calc(" dL2 = min(dL2, SoilDepth2-L2Wheel) ");
            //VJ 050812 correct DeltaF in case of impermeable
            calc(" DeltaF = min(SoilDepth2 - L2Wheel, DeltaF ");
            // adjust DeltaF if this is the case
//VJ 050829 fixed, infilpot should become 0 when impermeable
				calc(" InfilPotWheel = -DeltaF ");
         }

         calc(" L2Wheel += dL2 ");
         // depth of wetting front
         //VJ 050812 Adjust moisture content if not impermeable

         if (SwitchDrainage)
         {
              _spatial(REAL4, Kd);
				  calc(" Kd = mif(L2Wheel ge SoilDepth2, Ksat2, Ksat2*(ThetaI2/ThetaS2)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = ksat

				  calc(" ThetaI2 -= mif(L2Wheel lt SoilDepth2, (Kd*DTHOUR)/(SoilDepth2-L2Wheel),  0) ");
				  // adjust thetai if soil is not full
				  calc(" ThetaI2 = min(max(ThetaI2, 0.02*ThetaS2), ThetaS2-0.01) ");
				  // limit thetai between 0.02*pore and pore-0.01

				  calc(" L2Wheel -= mif(L2Wheel ge SoilDepth2, Kd*DTHOUR, 0) ");
				  // if profile full decrease L with drainage
         }
      }

      calc(" WHWheel -= DeltaF ");
        // adjust WH with new infil

