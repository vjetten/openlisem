/**************************************************/
/* GREEN AND AMPT MODULE */
/* all vars in mm */
/* VJ 050812: fixed a number of errors, version 2.36
/* VJ 050829: second revision, version 2.36
/* VJ 060109: number of bug fixes (brackets) version 2.392
/**************************************************/

// WHEELTRACK AFTER GRASSSTRIPS SEPARATE, CRUST AND COMPACTION AVERAGED FOR PIXEL
//
//NOTE: thetai < thetas is checked when reading in, it is always smaller
//NOTE: psi is checked for positive when reading in

      calc( " L1 = mif( ThetaI1 ge ThetaS1, SoilDepth1, L1) ");
        //VJ 051006 check for saturated soil

		//Calculate effective Ksat
		_spatial(REAL4, Ksateff);
		calc(" Ksateff = Ksat1*(1-CrustFraction-CompactFraction)");
		  //spatial average of ksat, all fractions are zero when not existant !!!

		if (SwitchCrustPresent)
			calc(" Ksateff += KsatCR*CrustFraction ");

		if (SwitchCompactPresent)
			calc(" Ksateff += KsatComp*CompactFraction ");

  		calc(" InfilPot = 0 ");
      //init infilpot for this timestep

	_spatial(REAL4, DeltaF);
	_spatial(REAL4, dL1);
   _spatial(REAL4, Kd);
   _spatial(REAL4, dL2);


//   float temp = DTHOUR;
//   DTHOUR = DTHOURGA;
    // nrGAsteps = 1;
//   for (int step = 0; step < nrGAsteps; step ++)
//   {

		calc(" DeltaF = DTHOUR*Ksateff*(1+(PSI1*10+WH)/(L1)) ");
			// potential infiltration according to G&A solution of Darcy, kutilek en Nielsen, p 138

		calc(" InfilPot -= DeltaF ");
			//potential infil needed for surplus kin wave

		calc(" DeltaF = min(WH, DeltaF) ");
			// no more infil than water on surface

		calc(" dL1 = mif(ThetaI1 lt ThetaS1, DeltaF/(ThetaS1-ThetaI1), 0) ");
			 // increase of wetting front 1st layer

		if (!SwitchInfilGA2)
		{
			if (SwitchImpermeable)
			{
				//cut off dL1
				calc(" dL1 = min(dL1, SoilDepth1-L1) ");
				//VJ 050812 correct DeltaF in case of impermeable
				calc(" DeltaF = min(SoilDepth1 - L1, DeltaF) ");

//VJ 050829 fixed, infilpot should become 0 when impermeable
				calc(" InfilPot = -DeltaF ");
			}

			calc(" L1 += dL1 ");
			// increase depth of wetting front

			if (SwitchDrainage)
			{
				//VJ 050829 changed drainage: decreases thetai or decreases L1

				  calc(" Kd = mif(L1 ge SoilDepth1, Ksateff, Ksateff*(ThetaI1/ThetaS1)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = ksat
              //VJ 051005 L2 cannot become greater than soildepth!, changed eq to ge

				  calc(" ThetaI1 -= mif(L1 lt SoilDepth1, (Kd*DTHOUR)/(SoilDepth1-L1),  0) ");
				  // adjust thetai if soil is not full
				  calc(" ThetaI1 = min(max(ThetaI1, 0.02*ThetaS1), ThetaS1-0.01) ");
				  // limit thetai between 0.02*pore and pore-0.01

				  calc(" L1 -= mif(L1 ge SoilDepth1, Kd*DTHOUR, 0) ");
				  // if profile full decrease L1 with drainage

			}
		}

		//********* SECOND LAYER GA ******************

		if (SwitchInfilGA2)
		{

			_spatial(UINT1, check);
			calc(" check = mif(L2 eq 0 and dL1+L1 gt SoilDepth1, 1, 0) ");
			calc(" check = mif(L2 gt 0, 2, check) ");
				//check: 0 water in first layer, 1: water entering 2nd; 2: water in second layer

//VJ 050812 recalc dL2 in case water just enters layer 2

			calc(" dL2 = 0 ");
			// set dL2 to 0 (e.g. check is 0)

			calc(" L1 += dL1 ");
			// increase water layer in first
			calc(" dL2 = mif(check eq 1 and ThetaI2 lt ThetaS2, (L1 - SoilDepth1)*(ThetaS1-ThetaI1)/(ThetaS2-ThetaI2), dL2) ");
			// dL2 is overshoot of L1, corrected for porosity in second layer, S always > I

			calc(" L1 = min(L1, SoilDepth1) ");
			// limit L1 to Soildepth1

         calc(" L2 = mif(check eq 0, 0, L2) ");
			  // if water is in the 1st layer L2 = 0
			calc(" L2 = mif(check eq 1, dL2, L2) ");
			  // if water just enters 2nd layer make a guess for L2
         calc(" L2 = mif(ThetaI2 ge ThetaS2, SoilDepth2, L2) ");
           // if layer allready full then L2 is layer depth

			_spatial(REAL4, Reff);
			calc(" Reff = mif(Ksateff eq 0 or Ksat2 eq 0, 1e20, L1/Ksateff + L2/Ksat2) ");
			  // effective hydraulic conductivity as resistance
			  // if either one is zero assume a very high resistance

          //average K = K1*L1/(L1+L2) + K2*L2/(L1+L2) = (K1L1+K2L2)/L1+L2)

			calc(" DeltaF = mif(check ne 0, DTHOUR*((PSI2*10+WH)+L1+L2)/Reff, DeltaF) ");
			  // if L2 gt 0, change deltaf else leave it at 1st layer deltaf

			calc(" InfilPot -= DeltaF ");
			  //adjust potential infil for kin wave when 2nd layer

			calc(" L2 = mif(check eq 1, 0, L2) ");
			  // put L2 back to 0 when just entering layer 2, to add the real dL2

			calc(" DeltaF = min(WH, DeltaF) ");
				// not more than water on surface

			calc(" dL2 = mif(ThetaI2 lt ThetaS2, DeltaF/(ThetaS2-ThetaI2), 0) ");
			// calc real dL2
			calc(" dL2 = mif(check gt 0, dL2, 0) ");
			//if not in layer 2 then set back to zero, just to be sure
			//in that case DeltaF only applies to L1

			if (SwitchImpermeable)
			{
				calc(" dL2 = min(dL2, SoilDepth2-L2) ");
				//VJ 050812 correct DeltaF in case of impermeable
				calc(" DeltaF = min(SoilDepth2 - L2, DeltaF) ");
				// adjust DeltaF if this is the case
//VJ 050829 fixed, infilpot should become 0 when impermeable
				calc(" InfilPot = -DeltaF ");
			}

			calc(" L2 += dL2 ");
			// depth of wetting front

			//VJ 050812 Adjust moisture content if not impermeable
			if (SwitchDrainage)
			{
				//VJ 050829 changed drainage: decreases thetai or decreases L1
				  calc(" Kd = mif(L2 ge SoilDepth2, Ksat2, Ksat2*(ThetaI2/ThetaS2)**DrFactor) ");
				  //if profile full then thetai = thetas, kd = ksat
              //VJ 051005 L2 can become greater than soildepth!, changed eq to ge

				  calc(" ThetaI2 -= mif(L2 lt SoilDepth2, (Kd*DTHOUR)/(SoilDepth2-L2),  0) ");
				  // adjust thetai if soil is not full
				  calc(" ThetaI2 = min(max(ThetaI2, 0.02*ThetaS2), ThetaS2-0.01) ");
				  // limit thetai between 0.02*pore and pore-0.01

				  calc(" L2 -= mif(L2 ge SoilDepth2, Kd*DTHOUR, 0) ");
				  // if profile full decrease L with drainage

			}
		}

		calc(" WH -= DeltaF ");
		  // adjust WH with new infil
  //}//10 sec steps

  //   DTHOUR = temp;

     calc(" WHCompact = WH ");
	  calc(" WHCrust = WH ");
		 // omslachtig maar nodig ivm SWATRE compatibiliteit

	  calc(" InfilPotCrust = InfilPot ");
	  calc(" InfilPotCompact = InfilPot ");


