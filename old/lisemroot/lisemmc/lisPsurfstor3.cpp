// *****************************************************************************
// ******* SURFACE storage in micro-depressions ********************************
// *****************************************************************************

           _spatial(REAL4, RunoffHin);
             // water depth of moving water
           calc(" RunoffHin = mif(MDS-SDS gt 0, (WH-SDS) * (1-exp(-WH*(WH-SDS)/(MDS-SDS))), WH)");
           calc(" RunoffHin = mif(WH le SDS, 0, RunoffHin)");

           if(SwitchSimpleDepression)
           {
              calc(" RunoffHin = mif(WH le MDS, 0, WH-MDS)");
           }

           _spatial(REAL4, PondAreaFract);
//           calc(" PondAreaFract = mif(RR gt 0.013, 1-exp(-PAFcoeff*WH), 1)");
//VJ 040224 BUG FIX WH was given as WH/10 which is 0.1mm!!!!!!!!
            //-PAFcoef to initial section
            //WH in mm here

//VJ 050506 new results from analysis of hydr radius and ponded area fraction
           calc(" PondAreaFract = 1-exp(-1.875*WH/(RR*10)) ");
          //NB WH en RR moeten dezelfde eenheden hebben, bijv allebei in mm

           _spatial(REAL4, SurfStorH);
           calc(" SurfStorH = WH-RunoffHin");

           //calc(" PondAreaFract = 1.0 ");
//VJ 050302 vervangen pond area formule

//VJ 040823 include buffers, no depr stor in buffers
           if (SwitchBuffers)
              calc(" PondAreaFract = mif(BufferID gt 0, 1, PondAreaFract)");

           // **** hydr radius (m), total overland flow (m3), avg detention (m)
           //VJ added for FAIR
           _spatial(REAL4, FlowWidth);
           _spatial(REAL4, RunoffVolin);
           _spatial(REAL4, WaterHVolin);

//VJ 030701 if wheeltracks as channels then not in flowwidth (same as Channels),
//if not wheeltracks then only compacted areas so only in infiltration and affecting
           calc(" FlowWidth = RoadWidthDX+(SoilWidthDX+StoneWidthDX)*PondAreaFract ");
           calc(" RunoffVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + "
                " RunoffHin*(SoilWidthDX+StoneWidthDX))*DXc/1000 , 0) ");
               // RunoffVolin is the amount of runoff water in a cell, m3
               // on roads, in wheeltracks
           calc(" WaterHVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + "
                "       WH*(SoilWidthDX+StoneWidthDX) )*DXc/1000 , 0) ");
/* OBSOLETE
           if (!SwitchWheelAsChannel)
           {
             calc(" FlowWidth = RoadWidthDX+(SoilWidthDX+StoneWidthDX)*PondAreaFract ");
             calc(" RunoffVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + "
                  " RunoffHin*(SoilWidthDX+StoneWidthDX))*DXc/1000 , 0) ");
                 // RunoffVolin is the amount of runoff water in a cell, m3
                 // on roads, in wheeltracks
             calc(" WaterHVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + "
                  "       WH*(SoilWidthDX+StoneWidthDX) )*DXc/1000 , 0) ");
           }
           else
           {
             calc(" FlowWidth = RoadWidthDX+(SoilWidthDX+StoneWidthDX)*PondAreaFract ");
             calc(" RunoffVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + WHWheelTrack*WheelWidthDX + "
                  "       RunoffHin*(SoilWidthDX+StoneWidthDX) )*DXc/1000 , 0) ");
                 // RunoffVolin is the amount of runoff water in a cell, m3
                 // on roads, in wheeltracks
             calc(" WaterHVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + WHWheelTrack*WheelWidthDX + "
                  "       WH*(SoilWidthDX+StoneWidthDX) )*DXc/1000 , 0) ");
           }
*/

           if (SwitchGrassPresent)
              calc(" FlowWidth = mif(GrassWidth gt 0, GrassFraction*DX+(1-GrassFraction)*FlowWidth, FlowWidth) ");

           _spatial(REAL4, RunoffMeanHin);
           calc(" RunoffMeanHin = mif(FlowWidth gt 0, RunoffVolin*1000/(FlowWidth*DXc), 0) ");

           _spatial(REAL4, Perimeter);
           calc(" Perimeter = FlowWidth+2*RunoffMeanHin/1000 ");

           _spatial(REAL4, HydrRadius);
           calc(" HydrRadius = mif(Perimeter gt 0,(FlowWidth*RunoffMeanHin/1000)/Perimeter, 0) ");
// just some tryouts
//           calc(" Perimeter = FlowWidth ");
//           calc(" HydrRadius = RunoffMeanHin/1000 ");
//           calc(" HydrRadius = mif(PondAreaFract gt 0, 0.001*(-0.0084*RR*10 + 0.7966)/PondAreaFract "
//                "* (RunoffMeanHin**(0.0027*RR*10 + 1.0124)),0) ");
//           writeTimeseries(PondAreaFract , "fpa");
// *****************************************************************************
// ******** overland flow rate calculation *************************************
// *****************************************************************************

           _spatial(REAL4, V);
           _spatial(REAL4, Beta);
           _spatial(REAL4, Alpha);
           _spatial(REAL4, Qin);
          //VJ NB:
          //RunoffVolin = amount of water available for runoff (potential) in m3
          //Qin = amount of water flowing between cells in m3/sec
          //so Qin*DTSEC = NOT EQUAL to RunoffVolin

           calc(" V = HydrRadius**(2.0/3.0) * sqrt(Gradient) / N ");
              // overland flow rate V, in meters per second
              // RunoffMeanHin is potential overland flow, in mm
              // Gradient is the slope in percentage, 10% = 0.10
              // N is Manning's n

           calc(" Beta = 0.6 ");
           calc(" Alpha = ( N/sqrt(Gradient) *(Perimeter**(2.0/3.0)) )**Beta ");
              // Alpha is the alpha as used in the kinematic wave formula
              // In: Applied Hydrology, p.283

           calc(" Qin = mif(Alpha gt 0 , ((FlowWidth * RunoffMeanHin/1000)/Alpha)**(1/Beta) , 0) ");
              // Qin, Overland Flow Discharge between cells, in m/s*m*m = m3/s
              // formula 9.3.3 from Chow Applied Hydrology


