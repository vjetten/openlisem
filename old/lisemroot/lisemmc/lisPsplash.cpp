// *****************************************************************************
// ******** SPLASH detachment **************************************************
// *****************************************************************************
  if (!SwitchNoErosion)
  {

     // **** calculation of kinetic energy

     _spatial(REAL4, KELeafDrainage);
     calc(" KELeafDrainage = max(15.8*sqrt(CropHeight)-5.87,0) ");
        // KELeafDrainage is the Kinetic Energy from Leaf Drainage (J/m2/mm)
        // CropHeight is the effective height of the plant canopy (m)
        // cannot be less than 0

     _spatial(REAL4, KEDirectThroughfall);
     calc(" KEDirectThroughfall = mif(RainIntensity > 1, 8.95+8.44*log10(RainIntensity), 0)");
        // KEDirectThroughfall is the Kinetic Energy from Direct Throughfall (J/m2/mm)
        // RainIntensity is the rainfall intensity (mm/h)
        // cannot be less than 0
        // prevent the case when RainIntensity equals 0: log10(0)!!

     // **** splash detachment


     // KEDirectThroughfall is the Kinetic Energy from Direct Throughfall (J/m2/mm)
     // KELeafDrainage is the Kinetic Energy from Leaf Drainage (J/m2/mm)
     // VegetatFraction is the soil cover by vegetation (fraction)
     // AggregateStab is the aggregate stability index
     // CohesionSoil is the soil cohesion in kPa
     // if aggregate stability could not be measured, the model
     // uses COH. The Aggregate Stability map should contain 0 or negative values!
     // /1000 is to go from grams per timestep to kg per timestep
     // from KE = 0 to KE = 10 the relationship is linear!,
     //    based on field experiments
     // these formulas are from EUROSEM, but modified and calibrated
     //    for the Limburg soils
     // NO SPLASH ON STONE COVERED SOILS
     // splash is calculated seperately for ponded and non ponded areas;
     // 1 - for ponded areas: WHall is all WH concentrated on the ponded surface

     _spatial(REAL4, WHall);
/*
     if (SwitchWheelPresent)
       calc(" WHall = mif(PondAreaFract gt 0, (WH/PondAreaFract*(SoilWidthDX+StoneWidthDX)"
            "+WHWheelTrack*WheelWidthDX)/(StoneWidthDX+SoilWidthDX+WheelWidthDX), 0) ");
     else
       calc(" WHall = mif(PondAreaFract gt 0, WH/PondAreaFract, 0) ");
*/
//VJ 030701 wat een bullshit, WH is gewwon gemiddelde water hoogte in de cell waar er water is
//dus:
     calc(" WHall = WH ");


     _spatial(REAL4, DetachmentSplash);
     _spatial(REAL4, DepositionSplash);
     calc(" DetachmentSplash = 0 ");
     calc(" DepositionSplash = 0 ");

     _spatial(REAL4, DetachmentDirectThroughfall);
        // splash detachment direct rainfall
     _spatial(REAL4, DetachmentLeafDrainage);
        // splash detachment leaf drip

     _spatial(REAL4, WH0);
     calc(" WH0 = exp(-1.48*WHall) ");
        // water height factor

//hier
     _spatial(REAL4, PondedAreaSplash);
     if (SwitchWheelPresent)
       calc(" PondedAreaSplash = (WheelWidthDX+PondAreaFract*SoilWidthDX)*DXc/1000 ");
     else
       calc(" PondedAreaSplash = PondAreaFract*SoilWidthDX*DXc/1000 ");
        // wet splash area

     _spatial(REAL4, RainfallDirectH);
        // direct rainfall
//VJ 040823 use corrected rainfall height here
     calc(" RainfallDirectH = RainHc ");

     _spatial(REAL4, ThroughfallH);
     calc(" ThroughfallH = (RainHc-InterceptionH)*(1-StemflowFraction) ");
//VJ 100116 CHANGED TO STEMFLOW FRACTION     
        // factor 0.6 = stemflow 40% !
        // this corresponds to an leaf to ground surface angle of 36.87 degrees
        // the effect of leaf drainage is neglegible when CH<0.15 m.
        // InterceptionH is already taking VegetatFraction into account


// WAT EEN ROTZOOI IS DIT !!!!!!
// GEEN DOCUMENTATIE

     // **** 1) direct rainfall on ponded areas

     calc(" DetachmentDirectThroughfall = mif(KEDirectThroughfall gt 10 and AggregateStab gt 0,  "
          "(2.82/AggregateStab*KEDirectThroughfall*WH0 + 2.96)*RainfallDirectH*PondedAreaSplash,"
          "KEDirectThroughfall/10*(2.82/AggregateStab*10*WH0 + 2.96)*RainfallDirectH*PondedAreaSplash) ");

     calc(" DetachmentDirectThroughfall = mif(KEDirectThroughfall gt 10 and AggregateStab le 0,"
          "(0.1033/CohesionSoil*KEDirectThroughfall * WH0 + 3.58)  "
          " * RainfallDirectH * PondedAreaSplash, DetachmentDirectThroughfall) ");

     calc(" DetachmentDirectThroughfall = mif(KEDirectThroughfall le 10 and AggregateStab le 0, "
          "    KEDirectThroughfall/10* (0.1033/CohesionSoil*10 * WH0 + 3.58) "
          "    * RainfallDirectH * PondedAreaSplash, DetachmentDirectThroughfall) ");

     // **** 2) throughfall on ponded areas,

     calc(" DetachmentLeafDrainage = mif(KELeafDrainage gt 10 and AggregateStab gt 0,  "
          "(2.82/AggregateStab*KELeafDrainage*WH0 + 2.96)*ThroughfallH*PondedAreaSplash,"
          "KELeafDrainage/10*(2.82/AggregateStab*10 * WH0 + 2.96)*ThroughfallH*PondedAreaSplash) ");

     calc(" DetachmentLeafDrainage = mif(KELeafDrainage gt 10 and AggregateStab le 0,  "
          "(0.1033/CohesionSoil*KELeafDrainage*WH0 + 3.58)*ThroughfallH*PondedAreaSplash, "
          "DetachmentLeafDrainage) ");

     calc(" DetachmentLeafDrainage = mif(KELeafDrainage le 10 and AggregateStab le 0,  "
          "KELeafDrainage/10*(0.1033/CohesionSoil*10*WH0 + 3.58)*ThroughfallH*PondedAreaSplash, "
          "DetachmentLeafDrainage) ");

     calc(" DetachmentSplash = (1-VegetatFraction)*DetachmentDirectThroughfall "
                       " + VegetatFraction*DetachmentLeafDrainage ");
// the types of splash detachment are added

     // **** 3 - direct rainfall for non-ponded areas: WH0 = 1

     calc(" PondedAreaSplash = (1-PondAreaFract)*SoilWidthDX*DX/1000 ");
     // dry area

     calc(" DetachmentDirectThroughfall = mif(KEDirectThroughfall gt 10 and AggregateStab gt 0,     "
          "    (2.82/AggregateStab*KEDirectThroughfall+2.96) * RainfallDirectH * PondedAreaSplash, "
          "    KEDirectThroughfall/10*(2.82/AggregateStab*10+2.96) * RainfallDirectH * PondedAreaSplash) ");
     calc(" DetachmentDirectThroughfall = mif(KEDirectThroughfall gt 10 and AggregateStab le 0,     "
          "    (0.1033/CohesionSoil*KEDirectThroughfall+3.58) * RainfallDirectH * PondedAreaSplash, DetachmentDirectThroughfall) ");
     calc(" DetachmentDirectThroughfall = mif(KEDirectThroughfall le 10 and AggregateStab le 0,     "
          "    KEDirectThroughfall/10*(0.1033/CohesionSoil*10+3.58) * RainfallDirectH * PondedAreaSplash, DetachmentDirectThroughfall) ");


     // **** 4 - throughfall for non-ponded areas: WH0 = 1

     calc(" DetachmentLeafDrainage = mif(KELeafDrainage gt 10 and AggregateStab gt 0,     "
          "    (2.82/AggregateStab*KELeafDrainage+2.96) * ThroughfallH * PondedAreaSplash,                            "
          "    KELeafDrainage/10*(2.82/AggregateStab*10+2.96) * ThroughfallH * PondedAreaSplash ) ");
     calc(" DetachmentLeafDrainage = mif(KELeafDrainage gt 10 and AggregateStab le 0,            "
          "    (0.1033/CohesionSoil*KELeafDrainage+3.58) * ThroughfallH * PondedAreaSplash, DetachmentLeafDrainage) ");
     calc(" DetachmentLeafDrainage = mif(KELeafDrainage le 10 and AggregateStab le 0,            "
          "    KELeafDrainage/10*(0.1033/CohesionSoil*10+3.58) * ThroughfallH * PondedAreaSplash, DetachmentLeafDrainage) ");


     calc(" DetachmentSplash += SplashDelivery*((1-VegetatFraction)*DetachmentDirectThroughfall "
                             "+ VegetatFraction*DetachmentLeafDrainage) ");
    // the types of splash detachment are added

     if (SwitchGrassPresent)
        calc(" DetachmentSplash = mif(GrassWidth gt 0, (1-GrassFraction)*DetachmentSplash,DetachmentSplash) ");
        // no splash detachment on field strips and waterways

     calc("DetachmentSplash *= (1 - StoneFraction) ");
     // NO SPLASH ON STONES

//VJ 080423 Snowmelt, no splash activity on snowcover     
     calc(" DetachmentSplash *= (1 - SnowCover) ");

//VJ 080614 no splash on hard surfaces
     calc(" DetachmentSplash = mif(hardsurface gt 0, 0, DetachmentSplash)");

 } //!noerosion
