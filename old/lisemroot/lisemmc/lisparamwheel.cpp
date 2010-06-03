
     //-------------------------------------------------------------------------
     //********************* Wheeltracks *******************
     //-------------------------------------------------------------------------

     _nonspatial(REAL4, WheelQOutflow );
     WheelQOutflow = 0.0;
     _nonspatial(REAL4, WheelSedOutflow);
     WheelSedOutflow = 0.0;
     _spatial(REAL4, WheelSedin);
     calc(" WheelSedin = 0.0");
     _spatial(REAL4, InfilPotWheel);
     calc(" InfilPotWheel = 0 ");
      // global vars needed anyway

     _spatial_input_if(REAL4, WheelWidth,mapname("Wheelwidth"),SwitchWheelAsChannel);
     _spatial_input_if(LDD, WheelLDD, mapname("lddwheel"), SwitchWheelAsChannel);
     _spatial_input_if(REAL4, WheelGradient, mapname("wheelgradient"), SwitchWheelAsChannel);
     _spatial_input_if(REAL4, WheelN, mapname("wheelman"), SwitchWheelAsChannel);
     _spatial_input_if(REAL4, WheelCohesion,mapname("wheelcohesion"), SwitchWheelAsChannel);
     _spatial_input_if(REAL4, WheelDepth,mapname("wheeldepth"), SwitchWheelAsChannel);
     _spatial_input_if(REAL4, WheelNumber,mapname("wheelnbr"), SwitchWheelAsChannel);

     _spatial_if(REAL4, RunoffVolinToWheel, SwitchWheelAsChannel);
     _spatial_if(REAL4, SedinToWheel, SwitchWheelAsChannel);

     if (SwitchWheelAsChannel)
     {
//VJ 080622 moved to here
         _spatial(REAL4, pitsout);     //water out of pits
         calc(" pitsout = 0.0");
         _spatial(REAL4, pitsoutsed);    //sediment out of pts
         calc(" pitsoutsed = 0.0");

         calc(" RunoffVolinToWheel = 0 ");
         calc(" SedinToWheel = 0");


         rangetest(WheelNumber,R_GE_LE, 1, 10,"@@Wheelnumber: nr of tracks in a cell");

         rangetest(WheelWidth,R_GE_LE, 0, DX,"@@Wheeltrack width");
         // wheel track width in m
         calc(" TMP_RANGE = sum(WheelWidth) ");
         SwitchWheelPresent = (TMP_RANGE > 0);
         if(!SwitchWheelPresent)
            LisemError("Cannot simulate wheel tracks, wheel width map is empty");

         calc(" WheelWidth = cover(WheelWidth, 0) ");

         rangetest(WheelGradient,R_GE_LE,0.0001,1,"@@Slope gradient values in wheeltracks (must be sine)");

         rangetest(WheelN,R_GE_LE, 0.0001,10, "@@Manning's n wheeltracks");
         calc(" WheelN = mif(WheelLDD eq 5, 10, WheelN) ");
//VJ 030702: give pits in wheeltracks a high value to simulate still standing water
// obstruction against the end of the track          

         rangetest(WheelCohesion, R_GE, 0.196, R_DUMMY,
             "@Cohesion values must be > 0.196 (detachment efficiency coef. Y < 1)"
             " error in WheelCOH.MAP (enter large values, e.g. 9999 for non-erodible surfaces)");

         // depth of wheeltracks in cm
         calc(" WheelDepth = cover(WheelDepth*10, 0) ");
         // convert to mm
         calc(" WheelNumber = cover(WheelNumber, 0)");
         calc(" WheelGradient = cover(WheelGradient, 0)");
         calc(" WheelN = cover(WheelN, 0)");

//         calc(" NR_VALS_START = count(AREA)");
//         calc(" NR_VALS_NOW = count(WheelWidth)");
//         if (NR_VALS_NOW < NR_VALS_START)
//           LisemError("wrong number of pixels in Wheelwid.map, non wheeltracks must have value 0!");

         _spatial(REAL4, WheelY);
         calc(" WheelY = mif(WheelCohesion <100, 1.0/(0.89+0.56*WheelCohesion), 0) ");
           // enter high ChannelCohesion (9999) for non-erodible surfaces
         calc(" WheelY = min(WheelY,1.0) ");
         calc(" WheelY = cover(WheelY, 0)");

//         _spatial(REAL4, L1Wheel);
//         calc(" L1Wheel = 0 ");
//????????????????????????????
//is alred created in lisparaminfil
     }

