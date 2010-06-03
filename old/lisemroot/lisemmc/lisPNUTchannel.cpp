// *****************************************************************************
// ******** kinematic wave calculation for CHANNEL flow ************************
// *****************************************************************************
//VJ added channel nutrients 030414

      calc(" ChanNutPSuspension   += NutPSusToChannel  ");
      calc(" ChanNutNH4Suspension += NutNO3SusToChannel");
      calc(" ChanNutNO3Suspension += NutNH4SusToChannel");
      calc(" ChanNutPSolution     += NutPSolToChannel  ");
      calc(" ChanNutNH4Solution   += NutNO3SolToChannel");
      calc(" ChanNutNO3Solution   += NutNH4SolToChannel");
      // add overland flow nutrients (kg)
      
      _spatial(REAL4, ChanNutPSusQin  );
      _spatial(REAL4, ChanNutNH4SusQin);
      _spatial(REAL4, ChanNutNO3SusQin);
      _spatial(REAL4, ChanNutPSolQin  );
      _spatial(REAL4, ChanNutNH4SolQin);
      _spatial(REAL4, ChanNutNO3SolQin);
      // nutrient fluxes (kg/s)

      _nonspatial(REAL4, SumChannelVolin);
      calc(" SumChannelVolin = sum(ChannelVolin) ");
      // used in mass balance error

      _nonspatial(REAL4, ChanSuminNutPSus  );
      _nonspatial(REAL4, ChanSuminNutNH4Sus);
      _nonspatial(REAL4, ChanSuminNutNO3Sus);
      _nonspatial(REAL4, ChanSuminNutPSol  );
      _nonspatial(REAL4, ChanSuminNutNH4Sol);
      _nonspatial(REAL4, ChanSuminNutNO3Sol);
      calc("ChanSuminNutPSus   = sum(ChanNutPSuspension  )");
      calc("ChanSuminNutNH4Sus = sum(ChanNutNH4Suspension)");
      calc("ChanSuminNutNO3Sus = sum(ChanNutNO3Suspension)");
      calc("ChanSuminNutPSol   = sum(ChanNutPSolution    )");
      calc("ChanSuminNutNH4Sol = sum(ChanNutNH4Solution  )");
      calc("ChanSuminNutNO3Sol = sum(ChanNutNO3Solution  )");
      // totals for mass balance error

     _spatial(REAL4, Channelqfactor);
      calc(" Channelqfactor = mif(ChannelVolin gt 0.00001, ChannelQin/ChannelVolin, 0) ");
      calc(" ChanNutPSusQin   = ChanNutPSuspension  *Channelqfactor ");
      calc(" ChanNutNH4SusQin = ChanNutNH4Suspension*Channelqfactor ");
      calc(" ChanNutNO3SusQin = ChanNutNO3Suspension*Channelqfactor ");
      calc(" ChanNutPSolQin   = ChanNutPSolution    *Channelqfactor ");
      calc(" ChanNutNH4SolQin = ChanNutNH4Solution  *Channelqfactor ");
      calc(" ChanNutNO3SolQin = ChanNutNO3Solution  *Channelqfactor ");
      // calc fluxes in kg * m3/s/m3 = kg/s

      _spatial(REAL4, ChanNutPSusQout  );
      _spatial(REAL4, ChanNutNH4SusQout);
      _spatial(REAL4, ChanNutNO3SusQout);
      _spatial(REAL4, ChanNutPSolQout  );
      _spatial(REAL4, ChanNutNH4SolQout);
      _spatial(REAL4, ChanNutNO3SolQout);
      // fluxes out kg/s
      
      _spatial(REAL4, ChannelQout);
      // discharge out m3, alread calculated in lisPMCchannel but done again here

      _spatial(REAL4, chinfil);
      if (SwitchChannelInfil)
         calc(" chinfil = -(ChannelKsat *  ChannelPerimeter /3600000.0) ");
       //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
      else
         calc(" chinfil = 0 ");

      minematic(ChannelQout, ChannelQin, chinfil,
      		      ChanNutPSusQout,ChanNutNH4SusQout,ChanNutNO3SusQout,ChanNutPSolQout,ChanNutNH4SolQout,ChanNutNO3SolQout,
      		      ChanNutPSusQin,ChanNutNH4SusQin,ChanNutNO3SusQin,ChanNutPSolQin,ChanNutNH4SolQin,ChanNutNO3SolQin,
           	    ChannelLDD, ChannelAlpha, ChannelBeta, DTSEC, DXc);

      calc(" ChanNutPOutflow    = sum(mif(Outlet1 eq 1, ChanNutPSusQout  , 0))");  
      calc(" ChanNutNH4Outflow  = sum(mif(Outlet1 eq 1, ChanNutNH4SusQout, 0))"); 
      calc(" ChanNutNO3Outflow  = sum(mif(Outlet1 eq 1, ChanNutNO3SusQout, 0))"); 
      calc(" ChanNutPOutflows   = sum(mif(Outlet1 eq 1, ChanNutPSolQout  , 0))"); 
      calc(" ChanNutNH4Outflows = sum(mif(Outlet1 eq 1, ChanNutNH4SolQout, 0))");
      calc(" ChanNutNO3Outflows = sum(mif(Outlet1 eq 1, ChanNutNO3SolQout, 0))");
       // channel nutrient flux at outlet       
      

      // **** Calculation New Sediment amount Sedout
      calc(" Channelqfactor = mif(ChannelQout gt 0.00001, ChannelVolout/ChannelQout, 0) ");
      calc(" ChanNutPSuspension   = Channelqfactor*ChanNutPSusQout  ");
      calc(" ChanNutNH4Suspension = Channelqfactor*ChanNutNH4SusQout");
      calc(" ChanNutNO3Suspension = Channelqfactor*ChanNutNO3SusQout");
      calc(" ChanNutPSolution     = Channelqfactor*ChanNutPSolQout  ");
      calc(" ChanNutNH4Solution   = Channelqfactor*ChanNutNH4SolQout");
      calc(" ChanNutNO3Solution   = Channelqfactor*ChanNutNO3SolQout");
      
      
      if (SwitchCorrectMassSED) // default true !
      {
          _nonspatial(REAL4, ChanSumoutNutPSus  );               
          _nonspatial(REAL4, ChanSumoutNutNH4Sus);               
          _nonspatial(REAL4, ChanSumoutNutNO3Sus);               
          _nonspatial(REAL4, ChanSumoutNutPSol  );               
          _nonspatial(REAL4, ChanSumoutNutNH4Sol);               
          _nonspatial(REAL4, ChanSumoutNutNO3Sol);               
          calc("ChanSumoutNutPSus   = sum(ChanNutPSuspension  )");
          calc("ChanSumoutNutNH4Sus = sum(ChanNutNH4Suspension)");
          calc("ChanSumoutNutNO3Sus = sum(ChanNutNO3Suspension)");
          calc("ChanSumoutNutPSol   = sum(ChanNutPSolution    )");
          calc("ChanSumoutNutNH4Sol = sum(ChanNutNH4Solution  )");
          calc("ChanSumoutNutNO3Sol = sum(ChanNutNO3Solution  )");
      
          if (ChanSumoutNutPSus > 0)
                calc(" ChanNutPSuspension += ChanNutPSuspension*(ChanSuminNutPSus-ChanSumoutNutPSus-ChanNutPOutflow*DTSEC)/ChanSumoutNutPSus ");
          if (ChanSumoutNutNH4Sus > 0)
             calc(" ChanNutNH4Suspension += ChanNutNH4Suspension*(ChanSuminNutNH4Sus-ChanSumoutNutNH4Sus-ChanNutNH4Outflow*DTSEC)/ChanSumoutNutNH4Sus ");
          if (ChanSumoutNutNO3Sus > 0)
                calc(" ChanNutNO3Suspension += ChanNutNO3Suspension*(ChanSuminNutNO3Sus-ChanSumoutNutNO3Sus-ChanNutNO3Outflow*DTSEC)/ChanSumoutNutNO3Sus ");
          if (ChanSumoutNutPSol > 0)
             calc(" ChanNutPSolution += ChanNutPSolution*(ChanSuminNutPSol-ChanSumoutNutPSol-ChanNutPOutflow*DTSEC)/ChanSumoutNutPSol ");
          if (ChanSumoutNutNH4Sol > 0)
             calc(" ChanNutNH4Solution += ChanNutNH4Solution*(ChanSuminNutNH4Sol-ChanSumoutNutNH4Sol-ChanNutNH4Outflow*DTSEC)/ChanSumoutNutNH4Sol ");
          if (ChanSumoutNutNO3Sol > 0)
             calc(" ChanNutNO3Solution += ChanNutNO3Solution*(ChanSuminNutNO3Sol-ChanSumoutNutNO3Sol-ChanNutNO3Outflow*DTSEC)/ChanSumoutNutNO3Sol ");
       }
      
      
