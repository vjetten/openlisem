// *****************************************************************************
// ****** rainfall & interception **********************************************
// *****************************************************************************

         _spatial(REAL4, RainIntensity);
         calc("RainIntensity = 0");
//VJ 080628 try until now the main loop. If interval exceeds the input then simply stop, not go on with 0         
     //    try{
             if(timestepindex >= ENDINTERVAL+DTMIN)
                timestepindex--;

             timeinput(RainIntensity,RainGaugeID,rainFileName,timestepindex,DTMIN,firsttime, 1);
             //VJ 080628 firsttime causes to open de file for reading

             if(timestepindex >= ENDINTERVAL+DTMIN)
                timestepindex++;
       //  }
       /*
         catch(Exception &E)
         {
               LisemWarning("no information on rainfall interval, Rain Intensity set to zero");
         }
         */
           _spatial(REAL4, RainH);
           calc(" RainH = RainIntensity * DTHOUR ");
               // RainIntensity = rainfall intensity (in mm/h).
               // multiplied by the time interval (in sec)
               // divided by 3600 (to go to hours) gives
               // the amount of rainfall in each time interval, in mm
//VJ 040823 added corrected rainfall height
           _spatial(REAL4, RainHc);
           calc(" RainHc = RainH * DX/DXc ");
               // height corrected for surface because water spreads out

           //_nonspatial(REAL4, RainfallAverageH); moved to lisparamtotals.cpp
           calc(" RainfallAverageH = sum(RainH)/count(RainH) ");
           // used in file and screen output, correct back to dx*dx for good intensities

           // **** peak rainfall detection
           if (RainfallAverageH >= PeakRainfall)
           {
                PeakRainTime = timestepindex ;
                PeakRainfall = RainfallAverageH;
           }

//VJ 040823 add real rainfall height here
           calc(" TotalRainVol += sum((RainH*DX*DX)/1000) ");
           // rainfall cum spatial in m3

//- interception seen as rigid storage SMax filling up and overflowing
//- overflow flux is identical to rainfall flux in intensity
//- SMax is the storage of the plants inside the gridcell, not the average storage of the gridcell
// so if a single tree inside a cell has an SMax of 2mm even if it covers 10%, the Smax of that cell is 2
// - therefore the same goes for LAI: the LAI of the plants inside the gridcell
// this is also easier to observe. The LAI from a satellite image is the average LAI of a cell, must be divided by Cover

           if (SwitchBuffers)
              calc(" CanopyStorage = mif(BufferID gt 0, 0, CanopyStorage) ");
//VJ 040823 include buffers, no interception in buffers

           calc(" CanopyStorage = mif(hardsurface gt 0, 0, CanopyStorage) ");
// VJ 100131 include no interception on hard surfaces

           _spatial(REAL4, InterceptionH);
           calc(" InterceptionH = InterceptionHCum ");
           // save momentarily because equation is for cumulative interception

//VJ 040823: use corrected rainfall height
           calc(" RainHCum += RainHc ");
           // total rainfall depth
           calc(" InterceptionHCum = mif(CanopyStorage gt 0, VegetatFraction*CanopyStorage* "
                 "(1-exp(-0.0653*LAI*RainHCum/CanopyStorage) ),0) ");
             //cumulative interception for grid cell,
             // CanopyStorage is calculated from LAI
             // canopy openess according to Aston (1979), based on Merriam (1960/1973)
             // 0.046*LAI = k = (1-p); p = 1-0.046*LAI ; Aston (1979)
             // REEVALUATED BY DE JONG AND JETTEN
             // note: LAI is not a pixel average, but the average for PER!

           calc(" InterceptionH = InterceptionHCum - InterceptionH ");
             // interception in this timestep
//           calc(" TotalInterceptionVol = sum((InterceptionHCum*(SoilWidthDX+StoneWidthDX)*DXc)/1000) ");
           calc(" TotalInterceptionVol = sum(InterceptionHCum*SoilWidthDX*DXc/1000) ");
               // total momentary interception in m3, no interception on roads
//VJ030415 changed names to cum


           _spatial(REAL4, ThroughfallH);
           calc(" ThroughfallH = (VegetatFraction * RainHc - InterceptionH)*(1-StemflowFraction) ");
//VJ 100510 if LAI already is the avertage LAI for a gridcell with cover included, the part falling on the LAI is the rainfall*cover
//VJ 100510 moved veg fraction to this formula else drainage is all rain - fraction of interception, overestimate!
//VJ 100116 CHANGED TO STEMFLOW FRACTION
        // factor 0.6 = stemflow 40% !
        // this corresponds to an leaf to ground surface angle of 36.87 degrees
        // the effect of leaf drainage is neglegible when CH<0.15 m.
        // InterceptionH is already taking VegetatFraction into account


// *****************************************************************************
// ****** snowmelt *************************************************************
// *****************************************************************************
//VJ 080423 Snowmelt
//after interception, add snowmelt to net rainfall

  if (SwitchSnowmelt)
  {
         _spatial(REAL4, SnowmeltIntensity);
         calc("SnowmeltIntensity = 0");
//VJ 080628 try until now the main loop. If interval exceeds the input then simply stop, not go on with 0
//         try{
            if(timestepindex >= ENDINTERVAL+DTMIN)
                timestepindex--;

             timeinput(SnowmeltIntensity,SnowID,snowmeltFileName,timestepindex,DTMIN, firsttime, 2);

             if(timestepindex >= ENDINTERVAL+DTMIN)
                timestepindex++;
/*
         }
         catch(Exception &E)
         {
               LisemWarning("no information on snowmelt interval, snowmelt Intensity set to zero");
         }
*/
         _spatial(REAL4, SnowmeltH);
         calc(" SnowmeltH = SnowmeltIntensity * DTHOUR ");

         calc(" RainHc += SnowmeltH * DX/DXc ");
             // add coorected snowmelt height to rainfall

         _nonspatial(REAL4, SnowmeltAverageH);
         calc(" SnowmeltAverageH = sum(SnowmeltH)/count(SnowmeltH) ");
         // used in file and screen output, correct back to dx*dx for good intensities

         calc(" RainfallAverageH += SnowmeltAverageH ");

         calc(" TotalRainVol += sum((SnowmeltH*DX*DX)/1000) ");
         // rainfall cum spatial in m3
  }


// *****************************************************************************
// ****** Add to waterheight ***************************************************
// *****************************************************************************


           // **** increase all area water heights with net rainfall
           calc(" WH += RainHc - InterceptionH ");

//VJ 040823 include buffers, add direct rainfall to buffers
           if (SwitchBuffers) {
              // subtract rainfall from still available volume
              calc(" BufferVolumeCurrent -= mif(BufferVolumeCurrent gt 0, (RainH*DX*DX)/1000, 0) ");
              // if still room in buffer water height is zero
              calc(" WH = mif(BufferVolumeCurrent gt 0, 0, WH) ");
              // if buffer is negative (overflow), then WH is overflow height in mm
              calc(" WH = mif(BufferVolumeCurrent lt 0, -(BufferVolumeCurrent*1000)/(DX*DXc), WH) ");
           }

           if (SwitchWheelPresent)
              calc(" WHWheelTrack += RainHc ");
              // no interception on wheeltracks

           calc(" WHRoad += RainHc ");
               // rainfall directly in channels is treated in the channel section below

