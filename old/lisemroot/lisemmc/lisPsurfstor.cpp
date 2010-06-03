// *****************************************************************************
// ******* SURFACE storage in micro-depressions ********************************
// *****************************************************************************

           _spatial(REAL4, RunoffHin);
           calc(" RunoffHin = 0 ");
           calc(" RunoffHin = mif(WH gt DeprStoreStartH and WH le DeprStoreRainH and (DeprStoreRainH-DeprStoreStartH) gt 0, "
                " (DeprStoreRainH-DeprStoreMaxH)*(WH-DeprStoreStartH)/(DeprStoreRainH-DeprStoreStartH),RunoffHin)");
           calc(" RunoffHin = mif((DeprStoreRainH-DeprStoreStartH) le 0, WH, RunoffHin)");
           calc(" RunoffHin = mif(WH gt DeprStoreRainH, WH-DeprStoreMaxH, RunoffHin)");
              /* DeprStoreMaxH = MDS die bereikt wordt bij een WH van DeprStoreRainH
                 en de runoff tot dat moment is dus (DeprStoreRainH-DeprStoreMaxH) */
              // RunoffHin is the detention (potential overland flow), in mm
              // which is equal to the sum of all water minus the storage in depressions
              // WH is the available amount of water in mm
              // DeprStoreStartH is the net amount of water needed for starting overland flow
              // DeprStoreRainH is the total amount of net rainfall needed to fill all depressions
              // DeprStoreMaxH is the maximum amount of depression storage
              // between DeprStoreStartH and DeprStoreRainH a linear relationship is assumed

           _spatial(REAL4, PondAreaFract);
//          calc(" PondAreaFract = mif(DeprStoreMaxH gt 0, "
//                "   1-exp(-5*PondAreaFractMax*(WH/DeprStoreMaxH)), 1) ");
//          calc(" PondAreaFract = mif(RR gt 0.013, 1-exp(-PAFcoeff*WH/10), 1)");

          calc(" PondAreaFract = 1-exp(-1.88*WH/(RR*10)) ");
          //VJ 050302 vervangen pond area formule
          //NB WH en RR moeten dezelfde waardes hebben

/* OLD linear way of doing ponded area
           calc(" PondAreaFract = 0 ");
           calc(" PondAreaFract = mif(WH gt DeprStoreStartH and WH le DeprStoreMaxH, "
                " PondAreaFractMax*(WH-DeprStoreStartH)/(DeprStoreMaxH-DeprStoreStartH), PondAreaFract) ");
           calc(" PondAreaFract = mif(WH gt DeprStoreMaxH and WH le DeprStoreRainH, "
                "   (1-PondAreaFractMax)*(WH-DeprStoreMaxH)/(DeprStoreRainH-DeprStoreMaxH), PondAreaFract) ");
           calc(" PondAreaFract = mif(WH gt DeprStoreRainH, 1, PondAreaFract)");
*/
              // NOTE DeprStoreMaxH etc can become 0 !!!!

           _spatial(REAL4, PondedIsolatedH);
/*           calc(" PondedIsolatedH = mif(DeprStoreMaxH gt 0, "
                " min(0.2*(1- ((WH-RunoffHin)/DeprStoreMaxH-0.75)/0.25),0.2),0) ");
           calc(" PondedIsolatedH = max(PondedIsolatedH, 0) ");
*/
           calc(" PondedIsolatedH = 0 ");
           //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111

           // summary DEPRESSION STORAGE: defined are:
           // - area with roads: RoadWidthDX (m)
           // infiltration is 0, no splash and flow detachment
           // - area with channels: ChannelWidthUpDX (m)
           // infiltration is 0
           // - area with wheeltracks: WheelWidthDX (m)
           // different infiltration characteristics, larger flow rate
           // - area with ponding: PondAreaFract*(SoilWidthDX+StoneWidthDX) (m)
           // PondAreaFract*SoilWidthDX (m) : area with waterlayer influencing e.g. splash detachment
           // - area with no ponding: (1-PondAreaFract)*(SoilWidthDX+StoneWidthDX) (m)
           // (1-PondAreaFract)*SoilWidthDX (m) : area with splash directly on soil surface
           // - area with ponding in isolated depressions: (PondedIsolatedH*PondAreaFract)*(SoilWidthDX+StoneWidthDX) (m)
           // no flow detachment
           // - area with ponding in non-isolated depressions: (PondAreaFract*(1-PondedIsolatedH))*(SoilWidthDX+StoneWidthDX) (m)
           // PondAreaFract*(1-PondedIsolatedH)*SoilWidthDX (m) : flow detachment
           // - in general: no detachment on stone covered surfaces

           // **** hydr radius (m), total overland flow (m3), avg detention (m)
           //VJ added for FAIR
           _spatial(REAL4, FlowWidth);
           _spatial(REAL4, RunoffVolin);
           _spatial(REAL4, WaterHVolin);

           if (SwitchWheelAsChannel)
           {
             calc(" FlowWidth = RoadWidthDX+(SoilWidthDX+StoneWidthDX)*PondAreaFract*(1-PondedIsolatedH) ");

             calc(" RunoffVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + "
                  " RunoffHin*(1-PondedIsolatedH)*(SoilWidthDX+StoneWidthDX))*DX/1000 , 0) ");
                 // RunoffVolin is the amount of runoff water in a cell, m3
                 // on roads, in wheeltracks, and non-isolated depressions
                 // the water in isolated depressions is kept (see kinematic wave)
                 // PondAreaFract should not be used here!!!!
             calc(" WaterHVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + "
                  "       WH*(SoilWidthDX+StoneWidthDX) )*DX/1000 , 0) ");
           }
           else
           {
             calc(" FlowWidth = RoadWidthDX+WheelWidthDX+(SoilWidthDX+StoneWidthDX)*PondAreaFract*(1-PondedIsolatedH) ");

             calc(" RunoffVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + WHWheelTrack*WheelWidthDX + "
                  "       RunoffHin*(1-PondedIsolatedH)*(SoilWidthDX+StoneWidthDX) )*DX/1000 , 0) ");
                 // RunoffVolin is the amount of runoff water in a cell, m3
                 // on roads, in wheeltracks, and non-isolated depressions
                 // the water in isolated depressions is kept (see kinematic wave)
                 // PondAreaFract should not be used here!!!!
             calc(" WaterHVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + WHWheelTrack*WheelWidthDX + "
                  "       WH*(SoilWidthDX+StoneWidthDX) )*DX/1000 , 0) ");
           }

           if (SwitchGrassPresent)
              calc(" FlowWidth = mif(GrassWidth gt 0, GrassFraction*DX+(1-GrassFraction)*FlowWidth, FlowWidth) ");

//           calc(" FlowWidth = mif(GrassWidth gt 0, DX, FlowWidth) ");
//           calc(" FlowWidth = mif(GrassWidth gt 9.0, GrassWidth, FlowWidth) ");
               // flowwidth is modified on grassstrips and grassed waterways
               // flow is considered as a thin layer of water here
               // a strip width > 9.0 m is considered a grassed waterway

           _spatial(REAL4, RunoffMeanHin);
           calc(" RunoffMeanHin = mif(FlowWidth gt 0, RunoffVolin*1000/(FlowWidth*DX), 0) ");
//VJ 030701 not sure if this is still needed
//           _spatial(REAL4, TotalMeanHin);
//           calc(" TotalMeanHin = mif(FlowWidth gt 0, WaterHVolin*1000/(FlowWidth*DX), 0) ");
               // mean runoff height (in mm) in wheeltracks, roads and
               // non-isolated depressions, water IN FLOW




