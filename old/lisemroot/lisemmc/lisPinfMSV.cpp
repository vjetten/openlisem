/**************************************************/
/* Morel-Seytoux and Verdin MODULE */
/* all vars in mm */
/**************************************************/
/*
    float time1 = (index-STARTINTERVAL)/60; // time in hours

    _spatial(REAL4, Bparam);
    calc(" Bparam = (WH+GCapDrive)*(ThetaS1-ThetaI1) ");
       // auxillary parameter, room in soil

    _spatial(REAL4, ksattemp);
    calc(" ksattemp = Ksat1 ");
    _spatial(REAL4, icumtemp);
    calc(" icumtemp = Icum ");

    _spatial(REAL4, tpond);
    calc(" tpond = mif(RainIntensity ne ksattemp, max(0, (ksattemp*Bparam)/(RainIntensity**2 - ksattemp*RainIntensity)),-999) ");
      // time to ponding

    _spatial(INT1, pondflag);
    calc(" pondflag = mif((time1 ge tpond and tpond gt 0) or (WH gt 0), 1, 0) ");
      // ponding if t > tp AND tp > 0 OR as long as there is water at the surface

    _spatial(REAL4, Aparam);
    calc(" Aparam = mif(pondflag eq 1,((B+icumtemp)**2)/(2*ksattemp*Bparam*((RainIntensity/ksattemp-1+0.01)**2)),0) ");

    _spatial(REAL4, infilrate);
    calc(" infilrate = mif(RainIntensity gt ksattemp,"
         " ksattemp + 0.5*sqrt(2/B*ksattemp*(Bparam+icumtemp)**2)/sqrt(time1-tpond+Aparam), RainIntensity)");
     // infiltration rate if p > ks (mm/h) equals rainfall or M-S & V

     _spatial(REAL4, infil);
    calc(" infil = min(WH, infilrate*DTHOUR) ");
     // cannot infiltrate more than there is
    calc(" infil = mif(Bparam eq 0, 0, infil) ");
     // no infiltration if soil at saturation (full)

    if (SwitchImpermeable)
    {
       _spatial(REAL4, L1);
       calc(" L1 = mif(ThetaS1-ThetaI1 gt 0, (InfilCum+infil)/(ThetaS1-ThetaI1), 0) ");
         // param for soil depth check
       calc(" infil = mif(L1 ge SoilDepth, 0, infil) ");
       calc(" InfilCum += infil ");
       // overall cumulative infil
    }

    calc(" Icum += mif(time1 lt tpond, infil, 0) ");
    // cumulative infiltration until ponding

    calc(" WH -= infil ");
     // update water height

    if (SwitchCompactPresent)
    {
        _spatial(REAL4, ksattemp);
        calc(" ksattemp = Ksat1 ");
        _spatial(REAL4, icumtemp);
        calc(" icumtemp = Icum ");

        _spatial(REAL4, tpond);
        calc(" tpond = mif(RainIntensity ne ksattemp, max(0, (ksattemp*Bparam)/(RainIntensity**2 - ksattemp*RainIntensity)),-999) ");

        _spatial(REAL4, Aparam);
        calc(" Aparam = ((B+icumtemp)**2)/(2*ksattemp*Bparam*((RainIntensity/ksattemp-1)**2)) ");

        _spatial(REAL4, infilrate);
        calc(" infilrate = mif(RainIntensity gt ksattemp,"
             " ksattemp + 0.5*sqrt(2/B*ksattemp*(Bparam+icumtemp)**2)/sqrt(time1-tpond+Aparam), RainIntensity)");
         // infiltration rate (mm/h) equals rainfall or M-S & V
*/
