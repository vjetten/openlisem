// *****************************************************************************
// ******** channel flow calculations ******************************************
// *****************************************************************************


     calc(" ChannelVolin += RunoffVolinToChannel ");
     
//VJ 040823: bug fix : rainfall volume added over projected cell DX and nor DXc
     calc(" ChannelVolin += mif(ChannelWidthDX gt 0, (RainH*ChannelWidthUpDX*DX)/1000, 0) ");
        // ADD rainfall, no interception

//     calc(" ChannelVolin -= ChannelInfilVol ");
//     calc(" ChannelVolin = max(ChannelVolin, 0)");
        // SUBTRACT the infiltration in channels, max 0, never happens, but never say never

     // ChannelVolin = (ChannelWidthDX*ChannelWH + ChannelSide*sqr(ChannelWH) )*DX

     // can be solved by a squareroot equation:

     begin{
       _spatial(REAL4, a);
       calc(" a = mif(ChannelVolin > 0, ChannelSide*DXc/ChannelVolin, 0 ) ");

       _spatial(REAL4, b);
       calc(" b = mif(ChannelVolin > 0, ChannelWidthDX*DXc/ChannelVolin, 0 ) ");

       _spatial(REAL4, c);
       calc(" c = -1.0 ");

        calc(" ChannelWH = mif(ChannelVolin gt 0 and a gt 0,              "
             "          ( -b + sqrt(b*b-4*a*c) ) / (2*a) , 0)     ");
     }end

     calc(" ChannelWH = mif(ChannelSide eq 0, ChannelVolin/DXc/ChannelWidthDX*1000, ChannelWH*1000) ");
       // ChannelWH is the water depth in the channel in mm

     calc(" ChannelWidthUpDX = mif(ChannelWidthDX gt 0, ChannelWidthDX+2*ChannelSide*ChannelWH/1000,0) ");
     calc(" ChannelWidthUpDX = cover(min(ChannelWidthUpDX, 0.9*DX), 0)");
      // channelwidth cannot be larger than 90% of the grid cell

     _nonspatial(REAL4,MAXWIDTH);
     calc(" MAXWIDTH = maxcell(ChannelWidthUpDX+RoadWidthDX) ");
     if (MAXWIDTH > DX)
     {
        LisemError("TOTAL ChannelWidthDX + RoadWidthDX exceeds cell width\
                      error in CHANWIDT.MAP or ROADWIDT.MAP");
     }

     // **** flow rate and flux in channel
     _spatial(REAL4, ChannelPerimeter);
     calc(" ChannelPerimeter = ChannelWidthDX+2*ChannelWH/1000*sqrt(ChannelSide*ChannelSide+1) ");
     _spatial(REAL4, ChannelCSArea);
     calc(" ChannelCSArea = ChannelWidthDX*ChannelWH/1000+ChannelSide*sqr(ChannelWH/1000)");
     _spatial(REAL4, ChannelV);
     calc(" ChannelV = mif(ChannelPerimeter gt 0, sqrt(ChannelGradient)/ChannelN * (ChannelCSArea/ChannelPerimeter)**(2.0/3.0), 0) ");
        // ChannelV is the channel flow rate, in meters per second
        // ChannelGradient is the slope in the channel, in percentage, 10% = 0.10
        // ChannelN is Mannings n in channel
        // ChannelWH is the water height in channel, in mm
        // ChannelWidthDX is the channel width, in meters
        // ChannelSide is the tangent of the channel edge angle,
        // ChannelSide = 0 for rectangular channel, 1 for 45 degrees

// ADR: changed ChannelQin calculation, based on eq. 9.3.3 from Chow:
     _spatial(REAL4, ChannelBeta);
     calc(" ChannelBeta = 0.6 ");
     _spatial(REAL4, ChannelAlpha);
     calc(" ChannelAlpha = (ChannelN/sqrt(ChannelGradient)*ChannelPerimeter**(2.0/3.0))**ChannelBeta ");
     _spatial(REAL4, ChannelQin);
     calc(" ChannelQin = mif(ChannelAlpha gt 0 , (ChannelCSArea/ChannelAlpha)**(1/ChannelBeta), 0 ) ");

// *****************************************************************************
// **** channel detachment and transport calculations **************************
// *****************************************************************************

     if (!SwitchNoErosion)
     {
         calc(" ChannelSedin += SedinToChannel ");
            // the overland flow sediment is added to the channel system

         // **** sediment concentration
         _spatial(REAL4, ChannelSedConcentration);
         calc(" ChannelSedConcentration = mif(ChannelVolin gt 0 ,ChannelSedin/ChannelVolin, 0) ");


        _spatial(REAL4, FluidDensity);
        calc(" FluidDensity = 1000 ");
        _spatial(REAL4, ChannelOmega);
        _spatial(REAL4, ChannelOmegaCrit);
        _spatial(REAL4, ChannelR);
        calc(" ChannelR = mif(ChannelPerimeter gt 0, ChannelCSArea/ChannelPerimeter,0) ");
        calc(" ChannelOmega = 9.8*ChannelGradient*FluidDensity*ChannelV*ChannelR ");
        calc(" ChannelOmegaCrit = 9.8*FluidDensity*ChannelR*0.004 ");
//writeTimeseries(OmegaCrit,"omega");

         _spatial(REAL4, ChannelF);
         calc(" ChannelF = 0.05");

        _spatial(REAL4, ChannelH);
         calc(" ChannelH = 0.1 ");
         _spatial(REAL4, term1);
         calc(" term1 = (1-ChannelH)*ChannelY*ChannelF*max(ChannelOmega-ChannelOmegaCrit,0)");
         _spatial(REAL4, term2);
         calc(" term2 = mif(ChannelWH gt 0,ChannelF*ChannelH*2650/((2650-FluidDensity)*9.8*ChannelWH/1000)*max(ChannelOmega-ChannelOmegaCrit,0),0)");


        _spatial(REAL4, dChannelSedin);
//        calc(" dSedin = mif(c_t gt 0, DXc * Y* F* max(Omega-OmegaCrit,0)*(1-SedConcentration/c_t), 0)");
        calc(" dChannelSedin = DXc*(term1 + term2 - ChannelSedConcentration*SettlingVelocity)");

       _spatial(REAL4, ChannelDetachmentFlow);
       _spatial(REAL4, ChannelDeposition);
       calc(" ChannelDetachmentFlow = max(dChannelSedin, 0)");
       calc(" ChannelDeposition = min(dChannelSedin, 0)");


          calc(" ChannelSedin += ChannelDetachmentFlow ");
          calc(" ChannelSedin += ChannelDeposition ");
          //calc(" ChannelSedin = max(ChannelSedin, 0) ");
            // ChannelSedin is the amount of sediment in kg available for transport
            // ChannelDeposition (negative value) is added, so ChannelSedin decreases!
     }// switch no erosion

// *****************************************************************************
// ******** kinematic wave calculation for CHANNEL flow ************************
// *****************************************************************************

     _spatial(REAL4, ChannelQsedin);
     calc(" ChannelQsedin = mif(ChannelVolin gt 0, ChannelSedin*ChannelQin/ChannelVolin, 0) ");
        // ChannelQsedin is the sediment discharge in kg/s (kg*m3/s/m3)

     _nonspatial(REAL4, SumChannelVolin);
     calc(" SumChannelVolin = sum(ChannelVolin) ");
     // used in mass balance error

     _nonspatial(REAL4, SumChannelSedin);
     calc(" SumChannelSedin = sum(ChannelSedin) ");
     // used in mass balance error
//VJ 050704 added channel infil in m2/s
     _spatial(REAL4, chinfil);
     if (SwitchChannelInfil)
       calc(" chinfil = -(ChannelKsat *  ChannelPerimeter /3600000.0) ");
       //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
     else
       calc(" chinfil = 0 ");

//     calc(" InfilVolinKinWave = chinfil*DTSEC*DXc ");
     //VJ 050831 calc basic input infil vol kin wave
     //VJ 060817 added DXc because it is a volume m2/s * m * s = m3     

     _spatial(REAL4, chtemp);
     calc(" chtemp = 0");//(DepositionFlow+DetachmentFlow+DetachmentSplash+DepositionSplash)/DXc ");

     _spatial(REAL4, ChannelQout);
     _spatial(REAL4, ChannelQsedout);

     kineDX(ChannelQout,ChannelQin,chinfil,
               ChannelQsedout,ChannelQsedin, chtemp,
               ChannelLDD, ChannelAlpha, ChannelBeta, DTSEC,DXc);
//VJ 050831 CHANGED: infil changes: it is the sum of all fluxes in the cell. Needed for infiltration calc

     calc(" ChannelQOutflow = sum(mif(Outlet1 eq 1,ChannelQout, 0))");
     calc(" ChannelSedOutflow = sum(mif(Outlet1 eq 1, ChannelQsedout, 0))");
     // recalc outlet discharge, kin wave can handle more pits

     calc(" ChannelCSArea = ChannelAlpha*(ChannelQout**ChannelBeta)");
     // ChannelCSArea is A (cross section) in m2

//hier
   //  _spatial(REAL4, ChannelVolout);
     calc(" ChannelVolout = ChannelCSArea*DXc ");
     // m3

     calc(" SumChannelVolout = sum(ChannelVolout) ");
     // spatial total of water in channel, m3 , after kin wave

 /* VJ 080217 not needed for water
     if (SwitchCorrectMass)
     {
         if (SumChannelVolout > 0)
            calc(" ChannelVolout += ChannelVolout*"
                 "(SumChannelVolin-SumChannelVolout-ChannelQOutflow*DTSEC)/SumChannelVolout ");
            // divide error caused by pits over network
         calc(" SumChannelVolout = sum(ChannelVolout) ");
     }
 */
 
     calc(" ChannelVolout = max(ChannelVolout, 0) ");
     calc(" ChannelVolin = ChannelVolout ");

//VJ 050831 corrected channel infil calc
     if (SwitchKinwaveInfil){
        calc(" TotalInfilVol += sum(min(chinfil*DTSEC+ChannelVolin, InfilVolinKinWave)) ");
        //nonspatial infil adjustment in m3
        calc(" InfilVol += cover(min(chinfil*DTSEC+ChannelVolin, InfilVolinKinWave),0) ");
        //spatial infil adjustment in m3
	     calc(" TotalInfilVol += max(0,SumChannelVolin - SumChannelVolout - ChannelQOutflow*DTSEC) ");
	  }

     if (!SwitchNoErosion)
     {
         // **** correct mass balance error SED
         _spatial(REAL4, ChannelSedout);
         calc(" ChannelSedout = mif(ChannelQout gt 0, ChannelQsedout/ChannelQout*ChannelVolout, 0) ");

         // sediment discharge

         if (SwitchCorrectMassSED)
         {
             calc(" SumChannelSedout = sum(ChannelSedout) ");
             if (SumChannelSedout > 0)
                calc(" ChannelSedout += (SumChannelSedin-SumChannelSedout-ChannelSedOutflow*DTSEC)*ChannelSedout/SumChannelSedout ");
         }

         calc(" ChannelSedout = max(ChannelSedout, 0) ");
         calc(" ChannelSedin = mif(ChannelWidthUpDX gt 0,ChannelSedout, 0) ");

         calc(" SumChannelSedout = sum(ChannelSedin) ");

             // recalc sediment outlet
     }// switch no erosion



