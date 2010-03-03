#include "model.h"


//---------------------------------------------------------------------------
void TWorld::GridCell(void)
{
   FOR_ROW_COL_MV
   {
      if (SwitchIncludeChannel)
      if (RoadWidthDX->Drc > 0)
          ChannelWidthUpDX->Drc = min(0.9*_dx-RoadWidthDX->Drc, ChannelWidthUpDX->Drc);
      // channel cannot be wider than _dx-road

      WheelWidthDX->Drc = 0;
      if (SwitchWheelPresent)
         WheelWidthDX->Drc = max(0,_dx-RoadWidthDX->Drc-ChannelWidthUpDX->Drc)/_dx * WheelWidth->Drc;
      // adjust wheelwidth in cells with other surfaces

      SoilWidthDX->Drc = max(0, _dx - ChannelWidthUpDX->Drc
                                    - GullyWidthDX->Drc
                                    - RoadWidthDX->Drc
                                    - WheelWidthDX->Drc);
   }
}
//---------------------------------------------------------------------------
void TWorld::SurfaceStorage(void)
{
   FOR_ROW_COL_MV
   {
       double RRm = 0.01*RR->Drc; // assume RR in cm convert to m
       double wh = WH->Drc, whflow = 0;
       double SDS;
       double mds = MDS->Drc;

       SDS = 0.1*MDS->Drc;

       fpa->Drc = 1-exp(-1.875*(wh/RRm));
       // fraction ponded area of a gridcell
       FlowWidth->Drc = max(0.01*_dx, fpa->Drc*SoilWidthDX->Drc + RoadWidthDX->Drc);
       // calculate flowwidth by fpa*surface + road, excludes channel already

       if (mds > 0)
          whflow = (wh-SDS) * (1-exp(-1000*wh*(wh-SDS)/(mds-SDS)));
          // non-linear release fo water from depression storage
       else
          whflow = wh;
       if (wh < SDS)
          whflow = 0;
       // subtract surface storage and calc water available for runoff, in m
       // assumed on soilsurface because there is the roughness

       WHstore->Drc = wh - whflow;
       // average water stored on flowwidth and not available for flow, in m

       //WHrunoff->Drc *= ((SoilWidthDX->Drc + RoadWidthDX->Drc)/FlowWidth->Drc);
       // WHrunoff = average water height incl roads and all water is now in FlowWidth,
       //so excluding non-ponded (dry) areas, thereby increasing the average water level

       WaterVolRunoff->Drc = DX->Drc*( whflow*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
       // total volume available for flow, surface + road

       if (GrassPresent->Drc == 1)
         FlowWidth->Drc = GrassWidthDX->Drc + (1-GrassFraction->Drc)*FlowWidth->Drc;
       // assume grassstrip spreads water over entire width

       if (FlowWidth->Drc > 0)
         WHrunoff->Drc = WaterVolRunoff->Drc/(DX->Drc*FlowWidth->Drc);
       else
         WHrunoff->Drc = 0;
       // average WHrunoff from soil surface + roads, because kin wave can only do one discharge
   }

   CalcVelDisch();
   // calculate Alpha, Q and V from WHrunoff and flowwidth
}
//---------------------------------------------------------------------------

