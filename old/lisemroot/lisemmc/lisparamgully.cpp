
//---------------------------------------------------------------------------
//********************* GULLIES *********************************************
//---------------------------------------------------------------------------

   _spatial(REAL4, GullyWidthDX);
   calc("GullyWidthDX = 0");
   _nonspatial(REAL4, GullyQOutflow);
   _nonspatial(REAL4, GullySedOutflow);
   GullyQOutflow = 0;
   GullySedOutflow = 0;
     _spatial(REAL4, GullyPits);

   // INPUT
   _spatial_input_if(REAL4, DEM, mapname("dem"), SwitchGullies);
//   _spatial_input_if(REAL4, GullyWidthInit, mapname("gullyinit"), SwitchGullies);
//SET AT 10 CM!!!
   _spatial_input_if(REAL4, GullyN, mapname("gullyn"), SwitchGullies);
   _spatial_input_if(REAL4, GullyCohesion, mapname("gullycoh"), SwitchGullies);
   _spatial_input_if(REAL4, GullyDepth2, mapname("gullydep"), SwitchGullies);
   _spatial_input_if(REAL4, BulkDensity1, mapname("bulkdens1"), SwitchGullies);
   _spatial_input_if(REAL4, BulkDensity2, mapname("bulkdens2"), SwitchGullies);
//VJ 040331 Included init gully dimensions
   _spatial_input_if(REAL4, GullyWidthInit, mapname("gulwinit"), (SwitchGullies && SwitchGullyInit));
   _spatial_input_if(REAL4, GullyDepthInit, mapname("guldinit"), (SwitchGullies && SwitchGullyInit));
//VJ 060404 added gully infil, 060508 added switchgiullyinfil, corrected gulksat1/2 name
   _spatial_input_if(REAL4, GullyKsat1, mapname("gulksat1"),SwitchGullies & SwitchGullyInfil);
   _spatial_input_if(REAL4, GullyKsat2, mapname("gulksat2"),SwitchGullies & SwitchGullyInfil);
     // in mm/h

    if (SwitchGullies)
    {
       celltest(LDD, DEM);
       celltest(LDD, GullyN);
       celltest(LDD, GullyCohesion);
       celltest(LDD, GullyDepth2);
       celltest(LDD, BulkDensity1);
       celltest(LDD, BulkDensity2);
       celltest(LDD, GullyWidthInit);
       celltest(LDD, GullyDepthInit);
       celltest(LDD, GullyKsat1);
       celltest(LDD, GullyKsat2);
    }

   _spatial_if(REAL4, DEMdown, SwitchGullies);
   _spatial_if(REAL4, DEMinit, SwitchGullies);
   _spatial_if(REAL4, GullyDepth, SwitchGullies);
   // gully depth changes
   _spatial_if(REAL4, GullyCSFraction, SwitchGullies);
   // gully fraction of square foot CSA
   _spatial_if(REAL4, GullyArea, SwitchGullies);
   // gully fraction of square foot CSA
   _spatial_if(REAL4, GullyAlpha, SwitchGullies);
   _spatial_if(REAL4, GullyQin, SwitchGullies);

   if (SwitchGullies)
   {
     calc(" DEMinit = DEM ");
     // to calc gully depth

//VJ 040331 Included init gully dimensions
     if (SwitchGullyInit)
     {
        calc(" GullyDepth = GullyDepthInit ");
     }
     else
     {
        calc(" GullyDepth = 0 ");
     }
     calc(" GullyCSFraction = 0 ");
     calc(" GullyArea = 0 ");

     _spatial(REAL4, UPSarea);
     _spatial(REAL4, UPSareain);
     calc("UPSareain = DX*DXc");
     accuflux(UPSarea,UPSareain, LDD);

//?????????     calc(" GullyWidthDX = GullyWidthInit ");

     _spatial(REAL4, FCrit);
     if (LisIFace->E_Fcritical->ItemIndex == 0)
        calc(" FCrit = mif((Gradient*(UPSarea/DX)**0.4 gt 0.5),1,0)");
//        calc(" FCrit = mif((Gradient*(UPSarea/DX)**0.4 gt 0.5) and (Gradient gt 0.04),1,0)");
     if (LisIFace->E_Fcritical->ItemIndex == 1)
        calc(" FCrit = mif((Gradient*(UPSarea/DX) gt 18) and (ln(UPSarea/Gradient) gt 6.8),1,0)");
     if (LisIFace->E_Fcritical->ItemIndex == 2)
        calc(" FCrit = mif(Gradient gt 0.025*(UPSarea*0.0001)**-0.4,1,0)");
     if (LisIFace->E_Fcritical->ItemIndex == 3)
        calc(" FCrit = mif(Gradient*(UPSarea/DX)**0.4 gt 0.72,1,0)");
     if (LisIFace->E_Fcritical->ItemIndex == 4)
        calc(" FCrit = mif((Gradient*(UPSarea/DX) gt 40) and (ln(UPSarea/Gradient) gt 9.8),1,0)");
     if (LisIFace->E_Fcritical->ItemIndex == 5)
        calc(" FCrit = 1 ");

     _nonspatial(REAL4, tgrad);
     tgrad = atof(LisIFace->E_ThresholdGrad->Text.c_str());
     if (LisIFace->E_Fcritical->ItemIndex == 5)
        tgrad = 0;
     _spatial(REAL4, FCrittemp1);
     _spatial(REAL4, FCrittemp2);
     calc(" FCrittemp1 = FCrit*mif(Gradient ge tgrad,1,0) ");

     REAL4 v = 0;
     edge(FCrittemp2,FCrittemp1,v);
     accuflux(FCrit, FCrittemp2, LDD);
     calc(" FCrit = mif(FCrit gt 0,1,FCrit) ");

     upstream(FCrittemp1, FCrit, LDD);
     downstream(FCrittemp2, FCrit, LDD);
     calc(" FCrit = mif(FCrittemp1 eq 1 and FCrittemp2 eq 1, 1, FCrit) ");
     calc(" FCrit = mif(FCrittemp1 eq 0 and FCrittemp2 eq 0, 0, FCrit) ");
     upstream(FCrittemp1, FCrit, LDD);
     downstream(FCrittemp2, FCrit, LDD);
     calc(" FCrit = mif(FCrittemp1 eq 1 and FCrittemp2 eq 1, 1, FCrit) ");
     calc(" FCrit = mif(FCrittemp1 eq 0 and FCrittemp2 eq 0, 0, FCrit) ");
     calc(" FCrittemp1 = FCrit ");
     edge(FCrit,FCrittemp1,v);

     _spatial_input(REAL4, nonFCrit, mapname("nonfcrit"));
     calc(" FCrit = mif(nonFCrit == 1,0,FCrit) ");
     writepath(FCrit,"fcrit.map");

     _nonspatial(REAL4, qwa);
     _nonspatial(REAL4, qwb);
     if (LisIFace->E_QWrelation->ItemIndex == 0)
     {
         SwitchGullyEqualWD = true;
//VJ 040329 added equal erosion over width and depth in gully, set in paramgully
     }
     if (LisIFace->E_QWrelation->ItemIndex == 1)
     {
          qwa = 12.85; qwb=0.59;
     }
     if (LisIFace->E_QWrelation->ItemIndex == 2)
     {
          qwa = 2.65; qwb=0.37;
     }
     if (LisIFace->E_QWrelation->ItemIndex == 3)
     {
            qwa = atof(LisIFace->E_QWparama->Text.c_str());
            qwb = atof(LisIFace->E_QWparamb->Text.c_str());
     }

       _spatial(REAL4,GullyVolin);
       calc("GullyVolin = 0");
       _spatial(REAL4,GullySedin);
       calc("GullySedin = 0");
       _spatial(REAL4, GullyWH);
       calc(" GullyWH = 0");

       _spatial(REAL4, GullyY);
       calc(" GullyY = mif(GullyCohesion gt 0.2, 1/(0.89+0.56*GullyCohesion), 1)");
       //     calc(" GullyY = 0.79*exp(-0.85*GullyCohesion) ");

       calc(" GullyDepth2 /= 100 ");
       // convert to meters


       _spatial(REAL4, DEM2);
       calc(" DEM2 = DEM-GullyDepth2 ");
/*
     _spatial(REAL4, MaxDepressionStorage);
           calc(" MaxDepressionStorage = max(0, 0.243*RR + 0.010*RR*RR - 0.012*RR*Gradient*100) ");
           _spatial(REAL4, PAFcoeff);
           calc(" PAFcoeff = 1.4057*(RR/10)**-0.9422 ");
           calc(" GullyWidthInit = mif(RR gt 0.013, 1-exp(-PAFcoeff*MaxDepressionStorage/10), 1)");
*/
//VJ 040331 Included init gully dimensions
     if (SwitchGullyInit)
     {
       calc("GullyWidthDX = GullyWidthInit*FCrit ");
     }
     else
     {
       calc("GullyWidthDX = 0.2*FCrit");
       // init at 0.2 m
     }
   }

