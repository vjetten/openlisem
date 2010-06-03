//==============================================================================
// ********************** WATER BALANCE TOTALS *********************************
//==============================================================================


            // ***** NON SPATIAL *****

           if (SwitchIncludeChannel)
               QOutflow += ChannelQOutflow;
           if (SwitchWheelPresent)
               QOutflow += WheelQOutflow;
           if (SwitchGullies)
               QOutflow += GullyQOutflow;

           TotalDischarge += QOutflow*DTSEC;
              // VJ: water combined for total outflow through outlet point,
              // before: was added to ChannelVolin and ChannelSedin before kinematic wave !
              // total discharge outlet

           // **** peak flow detection
           if ((QOutflow*1000) > PeakDischarge)      //qoutflow contains all outflows, channel surface gully etc
           {
                PeakTime = timestepindex ;
                PeakDischarge = QOutflow*1000;
              // calculate peak discharge, in l/s
              // calculate time of peak discharge, in min
           }

//VJ 080217 ADD BASEFLOW HERE
         //BaseflowDischarge = ChannelQBaseflow*1000;

           //*******  mass balance error water
//VJ 040823 include buffers, add sum volume map
//VJ 080217 add baseflow initial total as input
           if (TotalRainVol > 0)
              MassBalanceError = 100/TotalRainVol * (TotalRainVol + TotalBaseflowVol -
                  TotalInfilVol-TotalInterceptionVol-TotalDischarge
                  -SumWaterHVolout-SumChannelVolout-SumWheelVolout-SumGullyVolout
                  -SumBufferVolume-SumBufferVolumeChannel);
           else
              MassBalanceError = 0;
           // note: SumWaterHVolout has all the water, incl depression storage


           // ***** SPATIAL *****

           if (SwitchIncludeChannel)
              calc(" Qout += cover(ChannelQout,0) ");
           if (SwitchWheelPresent)
              calc(" Qout += cover(WheelQout,0) ");
           if (SwitchGullies)
              calc(" Qout += cover(GullyQout,0) ");
           // sum for map output



