// *****************************************************************************
// ******** gully flow calculations ********************************************
// *****************************************************************************


           calc(" GullyVolin += RunoffVolinToGully ");
             // add overland flow volume to gully

//VJ 040823: bug fix : rainfall volume added over projected cell DX and nor DXc
           calc(" GullyVolin += mif(GullyWidthDX gt 0, (RainH*GullyWidthDX*DX)/1000, 0) ");
             // add rainfall (no interception) to gully volume (m3)

           calc(" GullyWH = mif(GullyWidthDX gt 0, GullyVolin/(DXc*GullyWidthDX)*1000,0) ");
           // recalc WH in gully in mm

           // **** flow rate and flux in Gully
           _spatial(REAL4, GullyPerimeter);
           calc(" GullyPerimeter = GullyWidthDX+2*GullyWH/1000 ");

           _spatial(REAL4, GullyCSArea);
           calc(" GullyCSArea = GullyWidthDX*GullyWH/1000 ");
           _spatial(REAL4, GullyV);
           calc(" GullyV = mif(GullyPerimeter gt 0, sqrt(Gradient)/GullyN * (GullyCSArea/GullyPerimeter)**(2.0/3.0), 0) ");

           _spatial(REAL4, GullyBeta);
           calc(" GullyBeta = 0.6 ");
           _spatial(REAL4, GullyAlpha);
           calc(" GullyAlpha = (GullyN/sqrt(Gradient)*GullyPerimeter**(2.0/3.0))**GullyBeta ");
           calc(" GullyAlpha = mif(FCrit gt 0, GullyAlpha, Alpha) ");
           // gullyout can go further than gullies alone, enable kin wave to use all of ldd
           _spatial(REAL4, GullyQin);
           calc(" GullyQin = mif(GullyAlpha gt 0 , (GullyCSArea/GullyAlpha)**(1/GullyBeta), 0 ) ");

           calc(" GullyQin *= FCrit ");
           calc(" GullyVolin *= FCrit ");
             // only gully flow where FCrit eq 1

// *****************************************************************************
// **** Gully detachment and transport calculations **************************
// *****************************************************************************

           calc(" GullySedin += SedinToGully ");
              // the overland flow sediment is added to the Gully system

           _spatial(REAL4, GullyTransportCapacity);
           calc(" GullyTransportCapacity = mif(Gradient*GullyV*100 gt CriticalStreamPower,   "
                "   CGovers*(Gradient*GullyV*100-CriticalStreamPower)**(DGovers),0)    ");
           calc(" GullyTransportCapacity = 2650*min(GullyTransportCapacity,0.32) ");

         // **** sediment concentration

           _spatial(REAL4, GullySedConcentration);
           calc(" GullySedConcentration = mif(GullyVolin gt 0 ,GullySedin/GullyVolin,0) ");

         // **** detachment in Gullys

           _spatial(REAL4, CTCAux);
           calc(" CTCAux = SettlingVelocity*DTSEC*DXc*GullyWidthDX ");

//VJ 040329 changed cohesion calculation
           _spatial(REAL4, Ytemp);
           calc(" Ytemp = mif(GullyDepth gt GullyDepth2, "
                  "(2*GullyDepth2*Y+2*(GullyDepth-GullyDepth2)*GullyY)/"
                  "(2*GullyDepth2*Y+2*(GullyDepth-GullyDepth2)*GullyY+GullyWidthDX*GullyY), Y) ");
            // apparent Y is Y distributed along the perimeter

           _spatial(REAL4, GullyDetachmentFlow);
	       if (SwitchAltErosion)
    	   {
	    	   calc(" GullyDetachmentFlow = Ytemp* max(0, GullyTransportCapacity-GullySedConcentration)*GullyQin*DTSEC ");
	       }
    	   else
           {
	           calc(" GullyDetachmentFlow = Ytemp*CTCAux* max(0, GullyTransportCapacity-GullySedConcentration) ");
           }
           calc(" GullyDetachmentFlow = min(max(0,GullyTransportCapacity*GullyVolin-GullySedin),GullyDetachmentFlow) ");
              // detachment cannot be larger than the remaining transp.capacit

//           calc(" CTCAux = mif(GullyWH gt 0, (1-exp(-SettlingVelocity*DTSEC/(0.001*GullyWH)))*GullyVolin,0) ");
           calc(" CTCAux = mif(GullyWH gt 0, (1-exp(-SettlingVelocity*DTSEC/(0.001*GullyWH)))*GullyVolin,1) ");
           _spatial(REAL4, GullyDeposition);
           calc(" GullyDeposition = CTCAux * min(0,GullyTransportCapacity-GullySedConcentration)");
           calc(" GullyDeposition = max(GullyDeposition, min(0,GullyTransportCapacity*GullyVolin-GullySedin)) ");
              // deposition cannot be larger than the transp.capacit surplus


           calc(" GullyDetachmentFlow *= FCrit ");
           calc(" GullyDeposition *= FCrit ");
           calc(" GullySedConcentration *= FCrit ");


         calc(" GullyDetachmentFlow = mif(GullyWH gt MinimumHeight, GullyDetachmentFlow, 0)");
         calc(" GullyDeposition = mif(GullyWH gt MinimumHeight, GullyDeposition, -GullySedin)");
             //!!!!!!!!!!!!!!!!!
         calc(" GullySedin += GullyDetachmentFlow ");
         calc(" GullySedin += GullyDeposition ");
         calc(" GullySedin = max(GullySedin, 0) ");
              // GullySedin is the amount of sediment in kg available for transport
              // GullyDeposition (negative value) is added, so GullySedin decreases!

// *****************************************************************************
// ******** kinematic wave calculation for Gully flow ************************
// *****************************************************************************

           _spatial(REAL4, GullyQsedin);
           calc(" GullyQsedin = mif(GullyVolin gt 0, GullySedin*GullyQin/GullyVolin, 0) ");
              // GullyQsedin is the sediment discharge in kg/s (kg*m3/s/m3)

           _nonspatial(REAL4, SumGullyVolin);
           calc(" SumGullyVolin = sum(GullyVolin) ");
           // used in mass balance error

           _nonspatial(REAL4, SumGullySedin);
           calc(" SumGullySedin = sum(GullySedin) ");
           // used in mass balance error

//VJ 060404 added Gully infil in m2/s
           _spatial(REAL4, ginfil);
           if (SwitchGullyInfil)
           {
              _spatial(REAL4, GullyKsat );
              calc(" GullyKsat = GullyKsat1 ");
              calc(" GullyKsat = mif(GullyDepth gt SoilDepth2, GullyKsat2, GullyKsat) ");
              calc(" ginfil = -(GullyKsat *  GullyPerimeter /3600000.0) ");
           }
          //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
           else
             calc(" ginfil = 0 ");

           //ginfil is < 0 if calculated, corrected 060508  
           calc(" InfilVolinKinWave = -ginfil*DTSEC*DXc ");
           //VJ 060817 added DXc because it is a volume m2/s * m * s = m3

           _spatial(REAL4, gtemp);
           calc(" gtemp = 0");

           _spatial(REAL4, GullyQout);
           _spatial(REAL4, GullyQsedout);
           kineDX(GullyQout,GullyQin,ginfil,
                  GullyQsedout,GullyQsedin, gtemp,
                  Bufferdummy, Bufferdummy,
                  LDD, GullyAlpha, GullyBeta, DTSEC,DXc,int(SwitchSedtrap));
//VJ 050831 CHANGED: infil changes: it is the sum of all fluxes in the cell. Needed for infiltration calc
              // N.B. gullyout not restricted to gullies, outflow over normal LDD
           calc(" GullyQOutflow = sum(mif(Outlet1 eq 1,GullyQout, 0))");
           calc(" GullySedOutflow = sum(mif(Outlet1 eq 1, GullyQsedout, 0))");

           // kin wave can handle more pits now, so ensure output catchment is outlet1

     // ADR M3 calculation changed, now using eq. 9.3.3 of Chow:
           calc(" GullyCSArea = mif(GullyWidthDX gt 0, GullyAlpha*(GullyQout**GullyBeta),0)");
           // GullyCSArea is A (cross section) in m2

           _spatial(REAL4, GullyVolout);
           calc(" GullyVolout = GullyCSArea*DXc ");
           // m3

           calc(" SumGullyVolout = sum(GullyVolout) ");
           // spatial total of water in Gully, m3 , after kin wave
           if (SwitchCorrectMass)
           {
               _nonspatial(REAL4, GullyError);
               GullyError = SumGullyVolin-SumGullyVolout-GullyQOutflow*DTSEC;
               if (SumGullyVolout > 0)
                  calc(" GullyVolout += GullyError/SumGullyVolout*GullyVolout ");
               calc(" GullyVolout = max(GullyVolout, 0) ");
               // divide error caused by pits over network
           }

           calc(" SumGullyVolout = sum(GullyVolout) ");
           calc(" GullyVolin = GullyVolout ");

//VJ 060404 added gully infil in m2/s
     if (SwitchKinwaveInfil){
        calc(" InfilVol += cover(min(ginfil*DTSEC+GullyVolout, InfilVolinKinWave),0) ");
        //spatial infil adjustment in m3
	     calc(" TotalInfilVol += max(0,SumGullyVolin - SumGullyVolout - GullyQOutflow*DTSEC) ");
	  }

           // **** correct mass balance error SED
           _spatial(REAL4, GullySedout);
           calc(" GullySedout = mif(GullyQout gt 0, GullyQsedout*GullyVolout/GullyQout, 0) ");
           calc(" SumGullySedout = sum(GullySedout) ");
           // sediment discharge


           if (SwitchCorrectMassSED)
           {
               _nonspatial(REAL4, GullySedError);
               GullySedError = SumGullySedin-SumGullySedout-GullySedOutflow*DTSEC;
               if (SumGullySedout > 0)
                  calc(" GullySedout += GullySedError*GullySedout/SumGullySedout ");
               calc(" GullySedout = max(GullySedout, 0) ");

           }
           calc(" GullySedin = mif(GullyWidthDX gt 0,GullySedout, 0) ");
           calc(" SumGullySedout = sum(GullySedin) ");



// *****************************************************************************
// ******** ADJUST GULLY DIMENSIONS ********************************************
// *****************************************************************************
//VJ 040329 moved factor decl up, used now above
            // recalculate thse variables with new values

           _spatial(REAL4, oldwidth);
           calc(" oldwidth = GullyWidthDX ");
           
           calc(" GullyPerimeter = GullyWidthDX+2*GullyDepth ");

           _spatial(REAL4, factor);
           calc(" factor = mif(GullyPerimeter gt 0, 2*GullyDepth2/GullyPerimeter, 1) ");

           _spatial(REAL4, BD);
           calc(" BD = mif(GullyDepth gt GullyDepth2, factor*BulkDensity1+(1-factor)*BulkDensity2, BulkDensity1) ");
           //if bd1 is deeper than gullydepth2 then bd2
           _spatial(REAL4, volCS);
//           calc(" volCS = (GullyDetachmentFlow + GullyDeposition)/(BD*DXc) ");
           calc(" volCS = (GullyDetachmentFlow)/(BD*DXc) ");
              // cross section gully, kan negatief zijn
//no depo for the time being

//VJ 040329 added switch equal erosion over width and depth in gully, set in paramgully
          calc(" factor = mif(GullyDepth gt GullyDepth2, "
                  "(2*GullyDepth2*Y+2*(GullyDepth-GullyDepth2)*GullyY)/"
                  "(2*GullyDepth2*Y+2*(GullyDepth-GullyDepth2)*GullyY+GullyWidthDX*GullyY), 0) ");
//           calc(" factor = mif(GullyDepth gt GullyDepth2, "
//                  "(2*GullyDepth2/CohesionTotal+2*(GullyDepth-GullyDepth2)/GullyCohesion)/"
//                  "(2*GullyDepth2/CohesionTotal+2*(GullyDepth-GullyDepth2)/GullyCohesion+GullyWidthDX/GullyCohesion), 0) ");

           //cohesion total is coh layer 1 + plant coh

           _spatial(REAL4, deltaw);
           calc(" deltaw = 0");
           calc(" deltaw = mif(GullyDepth gt GullyDepth2 and GullyDepth gt 0, volCS*min(1.0,factor/(2*GullyDepth)), 0) ");
           //cannot be more width increase then volcs (min(1.0, ....)
           calc(" deltaw = mif(GullyDepth lt GullyDepth2 and GullyPerimeter gt 0, volCS/GullyPerimeter, deltaw) ");
//homogeneous situation:
              // vol divided over width and depth after ratio of 2*d/P and w/P
              // delta vol depth = deltad*w = w/P*vol -> deltad = vol/P
              // delta vol width = 2*deltaw*d = 2*d/P*vol -> deltaw = vol/P

           _spatial(REAL4, calcwidth);
           calc(" calcwidth = GullyWidthDX ");
           if (!SwitchGullyEqualWD)
           {
               calc(" calcwidth = qwa*GullyQin**qwb ");
//               calc(" GullyWidthDX = mif(GullyDepth lt GullyDepth2, qwa*GullyQin**qwb, GullyWidthDX) ");
                // if in first layer use Q-w relation
           }
//         calc(" GullyWidthDX += deltaw ");
//         calc(" GullyWidthDX = max(oldwidth, GullyWidthDX)*FCrit ");

           calc(" GullyWidthDX = (max(GullyWidthDX, calcwidth)+deltaw)*FCrit ");
           calc(" GullyWidthDX = min(DX*0.9,GullyWidthDX) ");
           // GULLY CANNOT HAVE THE WIDTH OF THE PIXEL!

           calc(" GullyDepth = mif(GullyWidthDX gt 0, (GullyDepth*oldwidth+ volCS)/GullyWidthDX, 0) ");
           calc(" GullyDepth = cover(GullyDepth, 0) ");

           calc(" DEM = DEMinit-GullyDepth ");
           // adjust slope angle
           calc(" GullyArea=GullyDepth*GullyWidthDX ");
           calc(" GullyCSFraction=GullyArea/0.0929 ");


