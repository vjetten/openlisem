/*---------------------------------------------------------------------------
project: openLISEM
name: lisOverland.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisOverland.cpp:
- calculate velocity and discharge
- calculate overland flow, calls kinematic wave
---------------------------------------------------------------------------*/

#include "model.h"

//---------------------------------------------------------------------------
void TWorld::CalcVelDisch(void)
{
   FOR_ROW_COL_MV
   {
       double Perim, beta = 0.6;
       double _23 = 2.0/3.0;
       double beta1 = 1/beta;

       // avg WH from soil surface and roads, over width FlowWidth

       Perim = 2*WHrunoff->Drc+FlowWidth->Drc;
       if (Perim > 0)
          R->Drc = WHrunoff->Drc*FlowWidth->Drc/Perim;
       else
          R->Drc = 0;

       Alpha->Drc = pow(N->Drc/sqrt(Grad->Drc) * pow(Perim, _23),beta);

       if (Alpha->Drc > 0)
          Q->Drc = pow((FlowWidth->Drc*WHrunoff->Drc)/Alpha->Drc, beta1);
       else
          Q->Drc = 0;

       V->Drc = pow(R->Drc, _23)*sqrt(Grad->Drc)/N->Drc;
    }
}
//---------------------------------------------------------------------------
void TWorld::OverlandFlow(void)
{
   ToChannel();
   // calculate what goes in the channel is used

   /*---- Water ----*/
   FOR_ROW_COL_MV
   {
       // recalculate after subtractions in "to channel"
       WaterVolin->Drc = DX->Drc * (WHrunoff->Drc*FlowWidth->Drc + WHstore->Drc*SoilWidthDX->Drc);
       // WaterVolin total water volume in m3 before kin wave, runoff may be adjusted in tochannel
       q->Drc = FSurplus->Drc*_dx/_dt;
       // infil flux in kin wave <= 0, in m2/s, use _dx bexcause in kiv wave DX is used
       Qoutflow->Drc = 0;
      // init current outflow in all pits, rst is 0,  in m3
    }

    /*---- Sediment ----*/
    if (SwitchErosion)
    {
       FOR_ROW_COL_MV
       {

           Qs->Drc =  Q->Drc * Conc->Drc;
           // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
           Qsoutflow->Drc = 0;
           // init outflow of sed in pits
       }
    }


   Qn->setMV();
   // flag all Qn gridcell with MV for in kin wave

   FOR_ROW_COL_MV
   {
     if (LDD->Drc == 5) // if outflow point, pit
     {
    	 //TODO: WHEN MORE PITS QPEAK IS FIRST INSTEAD OF MAIN PIT
        Kinematic(r,c, LDD, Q, Qn, Qs, Qsn, q, Alpha, DX, WaterVolin, SedVol, BufferVol, BufferSedVol);

     	  Qoutflow->Drc = Qn->Drc * _dt;
        	if (SwitchErosion)
        		Qsoutflow->Drc = Qsn->Drc * _dt;
			// these maps now contain m3 and kg per timestep in pit cells
     }
   }

   FOR_ROW_COL_MV
   {
      double WHoutavg = (Alpha->Drc*pow(Qn->Drc, 0.6))/(_dx-ChannelWidthUpDX->Drc);
      // WH based on A/dx = alpha Q^beta / dx
      //TODO _dx also needs to be corrected for wheeltracks and gullies

      WHroad->Drc = WHoutavg;
      // set road to average outflowing wh, no surface storage.

      WH->Drc = WHstore->Drc + WHoutavg;
      // add new average waterlevel (A/dx) to stored water

      //V->Drc = Qn->Drc/(WHoutavg*(_dx-ChannelWidthUpDX->Drc));
      // recalc velocity for output to map ????

      WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc );
       // new water volume after kin wave, all water incl depr storage

      InfilVolKinWave->Drc = q->Drc*_dt + WaterVolin->Drc - WaterVolall->Drc - Qn->Drc*_dt;
      //infiltrated volume is sum of incoming fluxes+volume before - outgoing flux - volume after

      if (SwitchErosion)
      {
         Conc->Drc = MaxConcentration(WaterVolall->Drc, SedVol->Drc, DEP->Drc);
         // correct for very high concentrations, 850 after Govers et al
         // recalc sediment volume
      }
   }
}
//---------------------------------------------------------------------------
