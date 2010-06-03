/**************************************************/
/* GREEN AND AMPT MODULE FOR WHEELTRACKS          */
/* all vars in mm                                 */
/**************************************************/

//NOTE regular infil GA is done before this

      _spatial(REAL4, DeltaF);
      calc(" DeltaF = DTHOUR*KsatWT*(1+(PSI1*10+WHWheelTrack)/(L1Wheel+0.001)) ");
         // infiltration according to G&A solution of Darcy, kutilek en Nielsen, p 138

      calc(" InfilPotWheel = -DeltaF ");
           //potential infil needed for surplus kin wave

		calc(" DeltaF = min(WHWheel, DeltaF) ");
			// no more infil than water on surface

		_spatial(REAL4, dL1);
		calc(" dL1 = DeltaF/(ThetaS1-ThetaI1) ");
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
				  calc(" Kd = mif(L1Wheel ge SoilDepth1, KsatWT, KsatWT*(ThetaI1/ThetaS1)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = ksatWT

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
         _spatial(UINT2, check);
         calc(" check = mif(L2 eq 0 and dL1+L1Wheel gt SoilDepth1, 1, 0) ");
         calc(" check = mif(L2 gt 0, 2, check) ");
            //check: 0 first layer, 1: water entering 2nd; 2: water in second
         
//VJ 050812 recalc dL2 in case water just enters layer 2
         _spatial(REAL4, dL2);
			calc(" dL2 = 0 ");
			// set dL2 to 0 (e.g. check is 0)

         calc(" L1Wheel += dL1 ");
         // increase water layer in first

         calc(" dL2 = mif(check eq 1, L1Wheel - SoilDepth1, dL2 ");
         // dL2 is overshoot of L1

         calc(" L1Wheel = min(L1Wheel, SoilDepth1) ");
           // limit L1 to Soildepth1

         calc(" L2 = mif(check eq 1, dL2, L2) ");
           // if water just enters 2nd layer make a guess for L2

         _spatial(REAL4, Reff);
         calc(" Reff = mif(KsatWT eq 0 or Ksat2 eq 0, 1e20, L1Wheel/KsatWT + L2/Ksat2 ");
           // effective hydraulic conductivity as resistance
			  // if both zero assume a very high resistance
			  
         calc(" DeltaF = mif(check gt 0, DTHOUR*((PSI2*10+WHWheelTrack)+L1Wheel+L2)/Reff), DeltaF) ");
           // if L2 gt 0, change deltaf else leave it at 1st layer deltaf

			calc(" InfilPotWheel = -DeltaF ");
			  //adjust potential infil for kin wave when 2nd layer

         calc(" L2 = mif(check eq 1, 0, L2) ");
           // put L2 back to add the real dL2

         calc(" DeltaF = min(WHWheelTrack, DeltaF) ");
            // not more than water on surface

         calc(" dL2 = DeltaF/(ThetaS2-ThetaI2) ");
            // calc real dL2
         calc(" dL2 = mif(check gt 0, dL2, 0) ");
			//if not in layer 2 then set back to zero, just to be sure
			//in that case DeltaF only applies to L1

         if (SwitchImpermeable)
         {
            calc(" dL2 = min(dL2, SoilDepth2-L2) ");
            //VJ 050812 correct DeltaF in case of impermeable
            calc(" DeltaF = min(SoilDepth2 - L2, DeltaF ");
            // adjust DeltaF if this is the case
//VJ 050829 fixed, infilpot should become 0 when impermeable
				calc(" InfilPotWheel = -DeltaF ");
         }

         calc(" L2 += dL2 ");
         // depth of wetting front
         //VJ 050812 Adjust moisture content if not impermeable
         
         if (SwitchDrainage)
         {
              _spatial(REAL4, Kd);
				  calc(" Kd = mif(L2 ge SoilDepth2, Ksat2, Ksat2*(ThetaI2/ThetaS2)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = ksat

				  calc(" ThetaI2 -= mif(L2 lt SoilDepth2, (Kd*DTHOUR)/(SoilDepth2-L2),  0) ");
				  // adjust thetai if soil is not full
				  calc(" ThetaI2 = min(max(ThetaI2, 0.02*ThetaS2), ThetaS2-0.01) ");
				  // limit thetai between 0.02*pore and pore-0.01

				  calc(" L2 -= mif(L2 ge SoilDepth2, Kd*DTHOUR, 0) ");
				  // if profile full decrease L with drainage
         }
      }

      calc(" WHWheelTrack -= DeltaF ");
        // adjust WH with new infil

