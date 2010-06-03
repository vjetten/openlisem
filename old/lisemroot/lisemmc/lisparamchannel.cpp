     //-------------------------------------------------------------------------
     //********************* CHANNELS **********************
     //-------------------------------------------------------------------------

     // declare channel maps
     _nonspatial(REAL4, ChannelQOutflow );
     _nonspatial(REAL4, ChannelSedOutflow);
//VJ 080217 add baseflow at outlet
     _nonspatial(REAL4, ChannelQBaseflow );


     _spatial_input_if(LDD, ChannelLDD, mapname("lddchan"),SwitchIncludeChannel);
     _spatial_input_if(REAL4, ChannelGradient, mapname("ChanGrad"),SwitchIncludeChannel);
     _spatial_input_if(REAL4, ChannelN, mapname("Chanman"),SwitchIncludeChannel);
     _spatial_input_if(REAL4, ChannelCohesion,mapname("ChanCoh"),SwitchIncludeChannel);
     _spatial_input_if(REAL4, ChannelWidthDX,mapname("ChanWidth"),SwitchIncludeChannel);
     _spatial_input_if(REAL4, ChannelSide  ,mapname("ChanSide"),SwitchIncludeChannel);
         // ChannelSide is the tangent of the channel edge angle,
         // ChannelSide = 0 for rectangular channel, 1 for 45 degrees

     if (SwitchIncludeChannel)
     {
        celltest(ChannelLDD, ChannelGradient);
        celltest(ChannelLDD, ChannelN);
        celltest(ChannelLDD, ChannelCohesion);
        celltest(ChannelLDD, ChannelWidthDX);
        celltest(ChannelLDD, ChannelSide);
     }

//VJ 050704 added channel infil here
     _spatial_input_if(REAL4, ChannelKsat, mapname("chanksat"),SwitchChannelInfil);
     // in mm/h
     if (SwitchChannelInfil){
        celltest(ChannelLDD, ChannelKsat);
        }

//VJ 080214 baseflow variables
     _spatial_input_if(REAL4, ChannelBaseflux, mapname("chanbaseflux"),SwitchChannelBaseflow);
     _spatial_input_if(REAL4, ChannelBaseincrease, mapname("chanincrease"),SwitchChannelBaseflow);
     if (SwitchChannelBaseflow){
        celltest(ChannelLDD, ChannelBaseflux);
        celltest(ChannelLDD, ChannelBaseincrease);
        }

     _spatial(REAL4, ChannelWidthUpDX);
     calc(" ChannelWidthUpDX = 0 ");
       // help var
     _spatial_if(REAL4, ChannelVolin, SwitchIncludeChannel);
     _spatial_if(REAL4, ChannelWH, SwitchIncludeChannel);

      // initalise channel maps
     if (SwitchIncludeChannel)
     {
         ChannelQOutflow = 0.0;
         ChannelSedOutflow = 0.0;

         rangetest(ChannelGradient,R_GE_LE,0.0001,10,"Slope gradient values in channels");
         rangetest(ChannelN,R_GE_LE, 0.001,10, "Manning's n channel");

         ChnCalibration = -1.0;
         ChnCalibration = GetFloat("Channel N calibration");
         _nonspatial(REAL4, chn_cal);
         chn_cal = ChnCalibration;
         if (ChnCalibration >= 0)
         {
           calc(" ChannelN = ChannelN * chn_cal ");
         }

         ChKsatCalibration = -1.0;
         ChKsatCalibration = GetFloat("Channel Ksat calibration");
         _nonspatial(REAL4, chksat_cal);
         chksat_cal = ChnCalibration;
         if (ChKsatCalibration >= 0 && SwitchChannelInfil)
         {
           calc(" ChannelKsat = ChannelKsat * chn_cal ");
         }


         rangetest(ChannelCohesion, R_GE, 0.196, R_DUMMY,
             "@Cohesion values must be > 0.196 (detachment efficiency coef. Y < 1)"
             " error in CHANCOH.MAP (enter large values, e.g. 9999 for non-erodible surfaces)");

        // first check if pixels are in [0,DX]
//         rangetest(ChannelWidthDX,R_LE, R_DUMMY, 0.9*DX, "@@Width of channel must be larger than 90% of DX");

         calc(" TMP_RANGE = sum(mif(ChannelWidthDX gt 0.9*DX,1,0)) ");
         if (TMP_RANGE > 0)
            LisemError(" Channel Width >= 0.9*CellWidth");
//         calc(" TMP_RANGE = sum(mif(ChannelWidthDX le 0.01,1,0)) ");
//         if (TMP_RANGE > 0)
//            LisemError(" Channel Width < 1 cm");
               // now check if they are in the ChannelLDD > 0
//         calc(" TMP_VALUE = mif(ChannelLDD,ChannelWidthDX) ");
//         rangetest(TMP_VALUE,R_GE_LE, 0.001, DX,"@@Channel bottom widths (map: CHANWIDT.MAP)");
                //    tmp variable like CHANWIDTH1 is in error message not specified by map name
                //    put that in the message like here: (map: CHANWIDTH.MAP)
                // @@ creates msg about cellwidth

         //VJ cut down channel maps to ldd mask
         calc(" ChannelCohesion = mif(ChannelLDD, ChannelCohesion) ");
         calc(" ChannelWidthDX = mif(ChannelLDD, ChannelWidthDX) ");
         calc(" ChannelSide = mif(ChannelLDD, ChannelSide) ");
         calc(" ChannelGradient = mif(ChannelLDD, ChannelGradient) ");
         calc(" ChannelN = mif(ChannelLDD, ChannelN) ");

         //VJ check nr channel cells, many users make errors!
         calc(" NR_VALS_START = count(ChannelLDD)");
         calc(" NR_VALS_NOW = sum(mif(ChannelGradient gt 0,1,0))");
         if (NR_VALS_NOW != NR_VALS_START)
            LisemError("Channel LDD and Channel gradient do not have the same nr of cells!");
         calc(" NR_VALS_NOW = sum(mif(ChannelCohesion gt 0,1,0))");
         if (NR_VALS_NOW != NR_VALS_START)
            LisemError("Channel LDD and Channel cohesion do not have the same nr of cells!");
         calc(" NR_VALS_NOW = sum(mif(ChannelN gt 0,1,0))");
         if (NR_VALS_NOW != NR_VALS_START)
            LisemError("Channel LDD and Channel Manning's N do not have the same nr of cells!");
         calc(" NR_VALS_NOW = sum(mif(ChannelWidthDX gt 0,1,0))");
         if (NR_VALS_NOW != NR_VALS_START)
            LisemError("Channel LDD and Channel width do not have the same nr of cells!");


         calc(" ChannelWidthUpDX = cover(ChannelWidthDX,0) ");
         // Width of channel with respect to channel precipitation,
         // since ChannelWidthDX is only the bottom width

         calc(" TMP_VALUE = ChannelWidthUpDX+RoadWidthDX ");
         rangetest(TMP_VALUE, R_GE_LE,0, DX,
             "@@Total of channel and road width in some cells!");

//VJ 031211 removed < 100 check for cohesion
         _spatial(REAL4, ChannelY);
         calc(" ChannelY = 1.0/(0.89+0.56*ChannelCohesion) ");
//         calc(" ChannelY = mif(ChannelCohesion <100, 1.0/(0.89+0.56*ChannelCohesion), 0) ");
           // enter high ChannelCohesion (9999) for non-erodible surfaces
         calc(" ChannelY = min(ChannelY,1.0) ");

         calc(" ChannelWH = mif(ChannelLDD, 0.0)");
         calc(" ChannelVolin = mif(ChannelLDD, 0.0)");

//VJ 080217 initialize channel volume to initial baseflow volume
         _spatial_if(REAL4, AddVolume, SwitchChannelBaseflow);
         if (SwitchChannelBaseflow)
         {
            ChannelQBaseflow = 0.0;
            
            calc(" AddVolume = 0 ");
            _spatial_input(REAL4, ChannelVolInit, mapname("chanvolini"));
            celltest(ChannelLDD, ChannelVolInit);

            calc(" ChannelVolin = mif(ChannelLDD, ChannelVolInit)");
         }

         _spatial(REAL4, ChannelSedin);
         calc(" ChannelSedin = mif(ChannelLDD, 0.0)");

}



