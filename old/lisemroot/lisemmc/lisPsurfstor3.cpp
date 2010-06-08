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

// VJ 100131 include no surface storage on hard surfaces
           calc(" RunoffHin = mif(hardsurface gt 0, WH, RunoffHin) ");

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
//           calc(" FlowWidth = max(0.01*DX, RoadWidthDX+(SoilWidthDX+StoneWidthDX)*PondAreaFract) ");
           calc(" FlowWidth = max(0.01*DX, RoadWidthDX + SoilWidthDX*PondAreaFract) ");
           //VJ 100502 added min 0.01DX here
//           calc(" RunoffVolin = (WHRoad*RoadWidthDX + RunoffHin*(SoilWidthDX+StoneWidthDX))*DXc/1000 ");
           calc(" RunoffVolin = (WHRoad*RoadWidthDX + RunoffHin*SoilWidthDX)*DXc/1000 ");
               // RunoffVolin is the amount of runoff water in a cell, m3
               // on roads, in wheeltracks
//           calc(" WaterHVolin = mif(FlowWidth gt 0, (WHRoad*RoadWidthDX + "
//                "       WH*(SoilWidthDX+StoneWidthDX) )*DXc/1000 , 0) ");
           calc(" WaterHVolin = (WHRoad*RoadWidthDX + WH*SoilWidthDX)*DXc/1000 ");
               // watervolin is total water (including storage)

           if (SwitchGrassPresent)
              calc(" FlowWidth = mif(GrassWidth gt 0, GrassFraction*DX+(1-GrassFraction)*FlowWidth, FlowWidth) ");

// *****************************************************************************
// ******** overland flow rate calculation *************************************
// *****************************************************************************

           _spatial(REAL4, RunoffMeanHin);
           calc(" RunoffMeanHin = mif(FlowWidth gt 0, RunoffVolin*1000/(FlowWidth*DXc), 0) ");
           // average out WH with all surface types in a cell, roads etc
           // so here the WH is increased because in fact it is devided by flowwith causing all the water to flow on a fpa surface

           _spatial(REAL4, Perimeter);
           calc(" Perimeter = FlowWidth+2*RunoffMeanHin/1000 ");

           _spatial(REAL4, HydrRadius);
           calc(" HydrRadius = mif(Perimeter gt 0,(FlowWidth*RunoffMeanHin/1000)/Perimeter, 0) ");

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


