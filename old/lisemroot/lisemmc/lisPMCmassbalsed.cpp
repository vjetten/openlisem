//==============================================================================
// ********************** SEDIMENT BALANCE TOTALS ******************************
//==============================================================================

           // ***** NON SPATIAL *****

           if (SwitchIncludeChannel)
              mcalc("SedOutflow# += ChannelSedOutflow#",6);

           calc("SedOutflow = SedOutflow0 ");
           calc("SedOutflow += SedOutflow1 ");
           calc("SedOutflow += SedOutflow2 ");
           calc("SedOutflow += SedOutflow3 ");
           calc("SedOutflow += SedOutflow4 ");
           calc("SedOutflow += SedOutflow5 ");

           calc(" TotalSedDischarge += SedOutflow0*DTSEC ");
           calc(" TotalSedDischarge += SedOutflow1*DTSEC ");
           calc(" TotalSedDischarge += SedOutflow2*DTSEC ");
           calc(" TotalSedDischarge += SedOutflow3*DTSEC ");
           calc(" TotalSedDischarge += SedOutflow4*DTSEC ");
           calc(" TotalSedDischarge += SedOutflow5*DTSEC ");

          mcalc(" TotalSedDischargeMu# += SedOutflow#*DTSEC ",6);

//LisemWarning("totdis ",TotalSedDischarge);
            // total soil loss outlet

           // *****  NON SPATIAL totals:
           // overland:
           calc(" TotalSplashDetachment += sum(DetachmentSplash) ");
           calc(" TotalFlowDetachment += sum(DF0) ");
           calc(" TotalFlowDetachment += sum(DF1) ");
           calc(" TotalFlowDetachment += sum(DF2) ");
           calc(" TotalFlowDetachment += sum(DF3) ");
           calc(" TotalFlowDetachment += sum(DF4) ");
           calc(" TotalFlowDetachment += sum(DF5) ");
//LisemWarning("totdf2 ",TotalFlowDetachment);
           calc(" TotalDeposition += sum(DepositionSplash)");//+sum(DepositionGrassStrip) ");
           calc(" TotalDeposition += sum(Dep0)");
           calc(" TotalDeposition += sum(Dep1)");
           calc(" TotalDeposition += sum(Dep2)");
           calc(" TotalDeposition += sum(Dep3)");
           calc(" TotalDeposition += sum(Dep4)");
           calc(" TotalDeposition += sum(Dep5)");

           // overland + channel suspended:
           _nonspatial(REAL4, SumSedInFlow);
           SumSedInFlow = SumSedout;
           
           // channels:
           if (SwitchIncludeChannel)
           {
              mcalc(" TotalChannelFlowDetachment += sum(cover(ChannelDF#,0)) ",6);
              mcalc(" TotalChannelDeposition += sum(cover(ChannelDep#,0)) ",6);
              SumSedInFlow += SumChannelSedout;
           }

                //*******  mass balance error sediment

           // NB deposition = negative !!!
           _nonspatial(REAL4, allerosion);
           allerosion = TotalFlowDetachment + TotalSplashDetachment
                        + TotalChannelFlowDetachment;

           _nonspatial(REAL4, alldeposition);
            alldeposition = TotalDeposition+TotalChannelDeposition;

           if (allerosion > 0)
              SedMassBalanceError =
               100.0*(allerosion + alldeposition - SumSedInFlow- TotalSedDischarge)/allerosion;
           else
              SedMassBalanceError = 0;

           // ***** SPATIAL *****

             calc(" SumTotErosion += DetachmentSplash");
             calc(" SumTotErosion += DF0");
             calc(" SumTotErosion += DF1");
             calc(" SumTotErosion += DF2");
             calc(" SumTotErosion += DF3");
             calc(" SumTotErosion += DF4");
             calc(" SumTotErosion += DF5");

             if (SwitchIncludeChannel)
                mcalc(" SumTotErosion += cover(ChannelDF#, 0) ",6);
//             if (SwitchWheelAsChannel)
//                calc(" SumTotErosion += mif(WheelWidthDX gt 0, WheelDetachmentFlow, 0) ");
                // the erosion amount is all the detachment.

             calc(" SumTotDeposition += DepositionSplash");
             calc(" SumTotDeposition += Dep0");
             calc(" SumTotDeposition += Dep1");
             calc(" SumTotDeposition += Dep2");
             calc(" SumTotDeposition += Dep3");
             calc(" SumTotDeposition += Dep4");
             calc(" SumTotDeposition += Dep5");

             if (SwitchIncludeChannel)
                mcalc(" SumTotDeposition += cover(ChannelDep#, 0) ",6);



