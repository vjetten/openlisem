//==============================================================================
// ********************** SEDIMENT BALANCE TOTALS ******************************
//==============================================================================
     if (!SwitchNoErosion)
    {
           // ***** NON SPATIAL *****

           if (SwitchIncludeChannel)
              SedOutflow += ChannelSedOutflow;
           if (SwitchWheelPresent)
              SedOutflow += WheelSedOutflow;
           if (SwitchGullies)
              SedOutflow += GullySedOutflow;

           TotalSedDischarge += SedOutflow*DTSEC;
             // vj : sediment combined for total outflow through outlet channel,
             // before: was added to ChannelVolin and ChannelSedin before kinematic wave !
            // total soil loss outlet

           // *****  NON SPATIAL totals:
           // overland:
           calc(" TotalSplashDetachment += sum(DetachmentSplash) ");
           calc(" TotalFlowDetachment += sum(DetachmentFlow) ");
           calc(" TotalDeposition += sum(DepositionSplash)+sum(DepositionFlow)");
           // channels:
           if (SwitchIncludeChannel)
           {
              calc(" TotalChannelFlowDetachment += sum(ChannelDetachmentFlow) ");
              calc(" TotalChannelDeposition += sum(ChannelDeposition) ");
           }
           // wheeltracks:
           if (SwitchWheelPresent)
           {
              calc(" TotalWheelFlowDetachment += sum(WheelDetachmentFlow) ");
              calc(" TotalWheelDeposition += sum(WheelDeposition) ");
           }
           // Gullies:
           if (SwitchGullies)
           {
              calc(" TotalGullyFlowDetachment += sum(GullyDetachmentFlow) ");
              calc(" TotalGullyDeposition += sum(GullyDeposition) ");
           }
             // NB do not count this for mass balance! is eroded again

           // overland + channel suspended:
           _nonspatial(REAL4, SumSedInFlow);
           SumSedInFlow = SumSedout + SumChannelSedout + SumWheelSedout + SumGullySedout;
//           -SumBufferSedVolume-SumBufferSedVolumeChannel;

                //*******  mass balance error sediment

           // NB deposition = negative !!!
           _nonspatial(REAL4, allerosion);
           allerosion = TotalFlowDetachment + TotalSplashDetachment
                        + TotalChannelFlowDetachment
                        + TotalWheelFlowDetachment
                        + TotalGullyFlowDetachment;

           _nonspatial(REAL4, alldeposition);

           alldeposition = TotalDeposition
                           +TotalChannelDeposition
                           +TotalWheelDeposition
                           +TotalGullyDeposition;

           if (allerosion > 0)
              SedMassBalanceError =
               100.0*(allerosion + alldeposition - SumSedInFlow - TotalSedDischarge)/allerosion;
           else
              SedMassBalanceError = 0;


           // ***** SPATIAL *****

           calc(" SumTotErosion += DetachmentSplash+DetachmentFlow");

           if (SwitchIncludeChannel)
              calc(" SumTotErosion += mif(ChannelWidthUpDX gt 0, ChannelDetachmentFlow, 0) ");
           if (SwitchWheelPresent)
              calc(" SumTotErosion += mif(WheelWidthDX gt 0, WheelDetachmentFlow, 0) ");
           if (SwitchGullies)
              calc(" SumTotErosion += mif(GullyWidthDX gt 0, GullyDetachmentFlow, 0) ");

           calc(" SumTotDeposition += (DepositionSplash+DepositionFlow)");//+DepositionGrassStrip)");

           if (SwitchIncludeChannel)
              calc(" SumTotDeposition += mif(ChannelWidthUpDX gt 0, ChannelDeposition, 0) ");
           if (SwitchWheelPresent)
              calc(" SumTotDeposition += mif(WheelWidthDX gt 0, WheelDeposition, 0) ");
           if (SwitchGullies)
              calc(" SumTotDeposition += mif(GullyWidthDX gt 0, GullyDeposition, 0) ");

           // ***** SPATIAL *****
           //VJ 050831  added spatial totals
           if (SwitchIncludeChannel)
              calc(" Qsedout += cover(ChannelQsedout,0) ");
           if (SwitchWheelPresent)
              calc(" Qsedout += cover(WheelQsedout,0) ");
           if (SwitchGullies)
              calc(" Qsedout += cover(GullyQsedout,0) ");
           // add for map output
 }
