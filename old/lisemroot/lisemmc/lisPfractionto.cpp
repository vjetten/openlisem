// *****************************************************************************
// ******** transport to the channel system ************************************
// *****************************************************************************


           _spatial_if(REAL4, RunoffVolinToChannel, SwitchIncludeChannel);
           _spatial_if(REAL4, SedinToChannel, SwitchIncludeChannel);

           _spatial_if(REAL4, SedinToChannel0, (SwitchIncludeChannel && SwitchMulticlass));
           _spatial_if(REAL4, SedinToChannel1, (SwitchIncludeChannel && SwitchMulticlass));
           _spatial_if(REAL4, SedinToChannel2, (SwitchIncludeChannel && SwitchMulticlass));
           _spatial_if(REAL4, SedinToChannel3, (SwitchIncludeChannel && SwitchMulticlass));
           _spatial_if(REAL4, SedinToChannel4, (SwitchIncludeChannel && SwitchMulticlass));
           _spatial_if(REAL4, SedinToChannel5, (SwitchIncludeChannel && SwitchMulticlass));

           _spatial_if(REAL4, NutPSolToChannel, (SwitchIncludeChannel && SwitchNutrients));
           _spatial_if(REAL4, NutNO3SolToChannel, (SwitchIncludeChannel && SwitchNutrients));
           _spatial_if(REAL4, NutNH4SolToChannel, (SwitchIncludeChannel && SwitchNutrients));
           _spatial_if(REAL4, NutPSusToChannel, (SwitchIncludeChannel && SwitchNutrients));
           _spatial_if(REAL4, NutNO3SusToChannel, (SwitchIncludeChannel && SwitchNutrients));
           _spatial_if(REAL4, NutNH4SusToChannel, (SwitchIncludeChannel && SwitchNutrients));
//VJ 030414 added solution and suspension

           _spatial(REAL4, FractionToChannel);
           calc(" FractionToChannel = 0 ");
             // needed for gully
           if (SwitchIncludeChannel) // if channel included
           begin{
               calc(" FractionToChannel = mif(ChannelWidthUpDX gt 0 and DX gt ChannelWidthUpDX, "
                    " min(DTSEC*V/(0.5*(DX-ChannelWidthUpDX)),1), 0) ");
               calc(" FractionToChannel = mif(DX eq ChannelWidthUpDX, 1, FractionToChannel)");
               calc(" FractionToChannel = mif(ChannelLDD, FractionToChannel)");
//VJ 040823 if buffer in channel than all overlandflow in channel
               if (SwitchBuffers)
                  calc(" FractionToChannel = mif(BufferID gt 0 and ChannelWidthUpDX gt 0, 1.0, FractionToChannel) ");

//VJ 090225 throw all water and sed in channel cell at main outlet, assume there is no overland flow adjacent to the channel
               if (SwitchAllinChannel)
                  calc(" FractionToChannel = mif(Outlet1 eq 1,1.0, FractionToChannel)");


//VJ 060505 bug fix: cover with 0 else problems with wheeltrack as channel
               calc(" FractionToChannel = cover(FractionToChannel, 0) ");

                 // depends on overall velocity in cell

               calc(" RunoffMeanHin -= FractionToChannel*RunoffMeanHin ");
               // fraction of runoff water flowing into channel
               //calc(" WH -= cover(FractionToChannel*RunoffMeanHin,0)");
               calc(" RunoffVolinToChannel = FractionToChannel*RunoffVolin ");
               calc(" RunoffVolin -= cover(RunoffVolinToChannel,0)");
               //calc(" WaterHVolin -= cover(RunoffVolinToChannel,0)");

               if (!SwitchMulticlass)
               {
                  calc(" SedinToChannel = FractionToChannel*Sedin ");
                  calc(" Sedin -= cover(SedinToChannel,0)");
               }

               if (SwitchMulticlass)
               {
                  mcalc(" SedinToChannel# = FractionToChannel*SedinMu# ",6);
                  mcalc(" SedinMu# -= cover(SedinToChannel#,0)", 6);
               }
               if (SwitchNutrients)
               {
                  calc(" NutPSolToChannel = FractionToChannel*NutPSolution ");
                  calc(" NutPSolution -= cover(NutPSolToChannel,0)");
                  calc(" NutNO3SolToChannel = FractionToChannel*NutNO3Solution ");
                  calc(" NutNO3Solution -= cover(NutNO3SolToChannel,0)");
                  calc(" NutNH4SolToChannel = FractionToChannel*NutNH4Solution ");
                  calc(" NutNH4Solution -= cover(NutNH4SolToChannel,0)");

                  calc(" NutPSusToChannel = FractionToChannel*NutPSuspension ");
                  calc(" NutPSuspension -= cover(NutPSusToChannel,0)");
                  calc(" NutNO3SusToChannel = FractionToChannel*NutNO3Suspension ");
                  calc(" NutNO3Suspension -= cover(NutNO3SusToChannel,0)");
                  calc(" NutNH4SusToChannel = FractionToChannel*NutNH4Suspension ");
                  calc(" NutNH4Suspension -= cover(NutNH4SusToChannel,0)");
//VJ 030414 added solution and suspension, 030415: corrected mcalc to calc
               }
    	   }end// if channel included

// *****************************************************************************
// ******** transport to the wheeltrack system *********************************
// *****************************************************************************
//VJ 030701 moved to paramwheels
//           _spatial_if(REAL4, RunoffVolinToWheel, SwitchWheelPresent);
//           _spatial_if(REAL4, SedinToWheel, SwitchWheelPresent);
           if (SwitchWheelPresent)
           begin{
               _spatial(REAL4, FractionToWheel);
               //calc(" FractionToWheel = cover(WheelWidthDX/DX, 0) ");
               // fraction is geometric only

//               calc(" FractionToWheel = mif(WheelWidthDX gt 0, "
//                    " min(DTSEC*V/(0.5*(DX-WheelWidthDX)/WheelNumber),1), 0) ");
//VJ 030701 no more account for wheelnumber, assumed in the middel with width wheelwidth

               calc(" FractionToWheel = mif(WheelWidthDX gt 0, "
                    " min(DTSEC*V/(0.5*(DX-WheelWidthDX)),1), 0) ");
               calc(" FractionToWheel = mif(DX eq WheelWidthDX, 1, FractionToWheel)");
//               calc(" FractionToWheel = mif(WheelLDD eq 5, 0, FractionToWheel)");
               calc(" FractionToWheel = mif(WHWheelTrack ge WheelDepth, 0, FractionToWheel)");  \
//VJ 030702: THIS FIXES A LOT OF PROBLEMS:no back flow into whelltrack if water higher or equal depth

                 // no inflow on pits otherwise circular flow
               calc(" FractionToWheel = cover(FractionToWheel, 0) ");
               // fraction depends on overall velocity

               calc(" FractionToWheel = mif(FractionToChannel+FractionToWheel gt 1, "
                    "  1-FractionToChannel, FractionToWheel)");
               // sum fraction le 1.0

               calc(" RunoffMeanHin -= FractionToWheel*RunoffMeanHin");
//               calc(" WH -= FractionToWheel*RunoffMeanHin");

               calc(" RunoffVolinToWheel = FractionToWheel*RunoffVolin ");
               calc(" RunoffVolin -= RunoffVolinToWheel ");
//               calc(" WaterHVolin -= RunoffVolinToWheel");

               calc(" SedinToWheel = FractionToWheel*Sedin ");
               calc(" Sedin -= SedinToWheel ");
           }end

    // -------------------------------------------------------------------------
    //***** transport to GULLY system ******************************************
    // -------------------------------------------------------------------------

            _spatial(REAL4, FractionToGully);
            calc(" FractionToGully = 0 ");
            _spatial_if(REAL4, RunoffVolinToGully, SwitchGullies);
            _spatial_if(REAL4, SedinToGully, SwitchGullies);

            if (SwitchGullies)
            begin{
                calc(" FractionToGully = min(1, mif(GullyWidthDX gt 0 and DX gt GullyWidthDX, "
                     " DTSEC*V/(0.5*(DX-GullyWidthDX)), 0)) ");
               calc(" FractionToGully = mif(DX eq GullyWidthDX, 1, FractionToGully)");
               calc(" FractionToGully = mif(FractionToChannel+FractionToGully gt 1, "
                    "  1-FractionToChannel, FractionToGully)");
               calc(" FractionToGully = cover(FractionToGully, 0)");

               calc(" RunoffVolinToGully = FractionToGully*RunoffVolin ");
               calc(" RunoffMeanHin -= cover(FractionToGully*RunoffMeanHin,0)");
//               calc(" WH -= cover(FractionToGully*RunoffMeanHin,0)");
               calc(" RunoffVolin -= cover(RunoffVolinToGully,0)");
//               calc(" WaterHVolin -= cover(RunoffVolinToGully,0)");

               calc(" SedinToGully = FractionToGully*Sedin ");
               calc(" Sedin -= SedinToGully ");
            }end


