
     //-------------------------------------------------------------------------
     //********************* Multiclass *******************
     //-------------------------------------------------------------------------
     _nonspatial(REAL4, TCCAL);
     TCCAL = 1.0;//atof(LisIFace->E_TCCalibration->Text.c_str());
     // calibration for TC usd in multiclas sediment

     _nonspatial(REAL4, SedOutflow0 );
     _nonspatial(REAL4, SedOutflow1 );
     _nonspatial(REAL4, SedOutflow2 );
     _nonspatial(REAL4, SedOutflow3 );
     _nonspatial(REAL4, SedOutflow4 );
     _nonspatial(REAL4, SedOutflow5 );

     SedOutflow0 = 0.0;
     SedOutflow1 = 0.0;
     SedOutflow2 = 0.0;
     SedOutflow3 = 0.0;
     SedOutflow4 = 0.0;
     SedOutflow5 = 0.0;

     _nonspatial(REAL4, ChannelSedOutflow0);
     _nonspatial(REAL4, ChannelSedOutflow1);
     _nonspatial(REAL4, ChannelSedOutflow2);
     _nonspatial(REAL4, ChannelSedOutflow3);
     _nonspatial(REAL4, ChannelSedOutflow4);
     _nonspatial(REAL4, ChannelSedOutflow5);

     ChannelSedOutflow0 = 0.0;
     ChannelSedOutflow1 = 0.0;
     ChannelSedOutflow2 = 0.0;
     ChannelSedOutflow3 = 0.0;
     ChannelSedOutflow4 = 0.0;
     ChannelSedOutflow5 = 0.0;

     _spatial_input_if(REAL4, FractionMu0, mapname("fractionmu0"), SwitchMulticlass);
     _spatial_input_if(REAL4, FractionMu1, mapname("fractionmu1"), SwitchMulticlass);
     _spatial_input_if(REAL4, FractionMu2, mapname("fractionmu2"), SwitchMulticlass);
     _spatial_input_if(REAL4, FractionMu3, mapname("fractionmu3"), SwitchMulticlass);
     _spatial_input_if(REAL4, FractionMu4, mapname("fractionmu4"), SwitchMulticlass);
     _spatial_input_if(REAL4, FractionMu5, mapname("fractionmu5"), SwitchMulticlass);
     if (SwitchMulticlass)
     {
        celltest(LDD,FractionMu0);
        celltest(LDD,FractionMu1);
        celltest(LDD,FractionMu2);
        celltest(LDD,FractionMu3);
        celltest(LDD,FractionMu4);
        celltest(LDD,FractionMu5);
     }

     _spatial_if(REAL4, SedinMu0, SwitchMulticlass);
     _spatial_if(REAL4, SedinMu1, SwitchMulticlass);
     _spatial_if(REAL4, SedinMu2, SwitchMulticlass);
     _spatial_if(REAL4, SedinMu3, SwitchMulticlass);
     _spatial_if(REAL4, SedinMu4, SwitchMulticlass);
     _spatial_if(REAL4, SedinMu5, SwitchMulticlass);

     _spatial_if(REAL4, SedConc0, SwitchMulticlass);
     _spatial_if(REAL4, SedConc1, SwitchMulticlass);
     _spatial_if(REAL4, SedConc2, SwitchMulticlass);
     _spatial_if(REAL4, SedConc3, SwitchMulticlass);
     _spatial_if(REAL4, SedConc4, SwitchMulticlass);
     _spatial_if(REAL4, SedConc5, SwitchMulticlass);

     _spatial_if(REAL4, D50susp, SwitchMulticlass);

     _spatial_if(REAL4, ChannelSedinMu0, SwitchMulticlass);
     _spatial_if(REAL4, ChannelSedinMu1, SwitchMulticlass);
     _spatial_if(REAL4, ChannelSedinMu2, SwitchMulticlass);
     _spatial_if(REAL4, ChannelSedinMu3, SwitchMulticlass);
     _spatial_if(REAL4, ChannelSedinMu4, SwitchMulticlass);
     _spatial_if(REAL4, ChannelSedinMu5, SwitchMulticlass);

     if (SwitchMulticlass)
     {

       _spatial(REAL4, mutot);
       calc(" mutot = 0 ");
       mcalc(" mutot += FractionMu# ",6);
       mcalc(" FractionMu# /= mutot ",6);

       _nonspatial(REAL4, mu0);
       _nonspatial(REAL4, mu1);
       _nonspatial(REAL4, mu2);
       _nonspatial(REAL4, mu3);
       _nonspatial(REAL4, mu4);
       _nonspatial(REAL4, mu5);
/*
        mu0=2;
        mu1=16.0;
        mu2=31.9;
        mu3=52.9;
        mu4=74.8;
        mu5=105.1;
*/
        mu0=LisIFace->ClassMu[0];
        mu1=LisIFace->ClassMu[1];
        mu2=LisIFace->ClassMu[2];
        mu3=LisIFace->ClassMu[3];
        mu4=LisIFace->ClassMu[4];
        mu5=LisIFace->ClassMu[5];

        mu_cl[0] = mu0;
        mu_cl[1] = mu1;
        mu_cl[2] = mu2;
        mu_cl[3] = mu3;
        mu_cl[4] = mu4;
        mu_cl[5] = mu5;
   /*
       _nonspatial(REAL4, CG0);
       _nonspatial(REAL4, CG1);
       _nonspatial(REAL4, CG2);
       _nonspatial(REAL4, CG3);
       _nonspatial(REAL4, CG4);
       _nonspatial(REAL4, CG5);
       _nonspatial(REAL4, DG0);
       _nonspatial(REAL4, DG1);
       _nonspatial(REAL4, DG2);
       _nonspatial(REAL4, DG3);
       _nonspatial(REAL4, DG4);
       _nonspatial(REAL4, DG5);
       calc(" CG0 = ((mu0+5)/0.32)**-0.6 ");
       calc(" CG1 = ((mu1+5)/0.32)**-0.6 ");
       calc(" CG2 = ((mu2+5)/0.32)**-0.6 ");
       calc(" CG3 = ((mu3+5)/0.32)**-0.6 ");
       calc(" CG4 = ((mu4+5)/0.32)**-0.6 ");
       calc(" CG5 = ((mu5+5)/0.32)**-0.6 ");
       calc(" DG0 = ((mu0+5)/300)**0.25 ");
       calc(" DG1 = ((mu1+5)/300)**0.25 ");
       calc(" DG2 = ((mu2+5)/300)**0.25 ");
       calc(" DG3 = ((mu3+5)/300)**0.25 ");
       calc(" DG4 = ((mu4+5)/300)**0.25 ");
       calc(" DG5 = ((mu5+5)/300)**0.25 ");
*/
       mcalc(" SedinMu# = 0 ",6);
       mcalc(" ChannelSedinMu# = 0.0", 6);

       _nonspatial(REAL4, SV0);
       _nonspatial(REAL4, SV1);
       _nonspatial(REAL4, SV2);
       _nonspatial(REAL4, SV3);
       _nonspatial(REAL4, SV4);
       _nonspatial(REAL4, SV5);
       calc(" SV0 = 2*(2650-1000)*9.80*sqr(mu0/2000000)/(9*0.001) ");
       calc(" SV1 = 2*(2650-1000)*9.80*sqr(mu1/2000000)/(9*0.001) ");
       calc(" SV2 = 2*(2650-1000)*9.80*sqr(mu2/2000000)/(9*0.001) ");
       calc(" SV3 = 2*(2650-1000)*9.80*sqr(mu3/2000000)/(9*0.001) ");
       calc(" SV4 = 2*(2650-1000)*9.80*sqr(mu4/2000000)/(9*0.001) ");
       calc(" SV5 = 2*(2650-1000)*9.80*sqr(mu5/2000000)/(9*0.001) ");

       _spatial(REAL4, D50soil);
       findD50(D50soil,FractionMu0,FractionMu1,FractionMu2,FractionMu3,FractionMu4,FractionMu5,*mu_cl);
       write(D50soil,"d50soil.map");

       _spatial(REAL4, CGsoil);
       _spatial(REAL4, DGsoil);
       calc(" CGsoil = ((D50soil+5)/0.32)**-0.6 ");
       calc(" DGsoil = ((D50soil+5)/300)**0.25 ");
   }

