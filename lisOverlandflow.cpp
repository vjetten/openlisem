/*---------------------------------------------------------------------------
project: openLISEM
name: lisOverland.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisOverland.cpp:
- calculate velocity and discharge
- calculate overland flow, calls kinematic wave
---------------------------------------------------------------------------*/

#include "model.h"

//---------------------------------------------------------------------------
//fraction of water and sediment flowing into the channel
void TWorld::ToChannel(void)
{
	if (SwitchIncludeChannel)
	{
		FOR_ROW_COL_MV_CH
		{
			double fractiontochannel = min(_dt*V->Drc/(0.5*(_dx-ChannelWidthUpDX->Drc)), 1.0);
			double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

			if (SwitchAllinChannel)
				if (Outlet->Drc == 1)
					fractiontochannel = 1.0;
			// in catchment outlet cell, throw everything in channel

			if (SwitchBuffers)
				if (BufferID->Drc > 0)
					fractiontochannel = 1.0;
			// where there is a buffer in the channel, all goes in the channel

			RunoffVolinToChannel->Drc = fractiontochannel*Volume;
			// water diverted to the channel
			WHrunoff->Drc *= (1-fractiontochannel);
			// adjust water height
			if (SwitchErosion)
			{
				SedToChannel->Drc = fractiontochannel*Sed->Drc;
				//sediment diverted to the channel
				Sed->Drc -= SedToChannel->Drc;
				// adjust sediment in suspension
			}
		}
		CalcVelDisch();
		// recalc velocity and discharge
	}
}
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

	/*---- Water ----*/
	// recalculate water vars after subtractions in "to channel"
	FOR_ROW_COL_MV
	{
		WaterVolin->Drc = DX->Drc * (WHrunoff->Drc*FlowWidth->Drc + WHstore->Drc*SoilWidthDX->Drc);
		// WaterVolin total water volume in m3 before kin wave, WHrunoff may be adjusted in tochannel
		q->Drc = FSurplus->Drc*_dx/_dt;
		// infil flux in kin wave <= 0, in m2/s, use _dx bexcause in kiv wave DX is used
		Qoutflow->Drc = 0;
		// init current outflow in all pits, rst is 0,  in m3
	}

	/*---- Sediment ----*/
	if (SwitchErosion)
	{
		// calc seediment flux going in kin wave as Qs = Q*C
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

	// do kin wave for all pits
	FOR_ROW_COL_MV
	{
		if (LDD->Drc == 5) // if outflow point, pit
		{
			//TODO: WHEN MORE PITS QPEAK IS FIRST INSTEAD OF MAIN PIT
			Kinematic(r,c, LDD, Q, Qn, Qs, Qsn, q, Alpha, DX, WaterVolin, Sed, BufferVol, BufferSed);

			Qoutflow->Drc = Qn->Drc * _dt;
			if (SwitchErosion)
				Qsoutflow->Drc = Qsn->Drc * _dt;
			// these maps now contain m3 and kg per timestep in pit cells
		}
	}

	// calculate resulting flux Qn back to water height on surface
	FOR_ROW_COL_MV
	{
		double WHoutavg = (Alpha->Drc*pow(Qn->Drc, 0.6))/(_dx-ChannelWidthUpDX->Drc);
		// WH based on A/dx = alpha Q^beta / dx
		//TODO _dx also needs to be corrected for wheeltracks and gullies

		WHroad->Drc = WHoutavg;
		// set road to average outflowing wh, no surface storage.

		WH->Drc = WHoutavg + WHstore->Drc;
		// add new average waterlevel (A/dx) to stored water

		//V->Drc = Qn->Drc/(WHoutavg*(_dx-ChannelWidthUpDX->Drc));
		// recalc velocity for output to map ????

		WaterVolall->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc );
		// new water volume after kin wave, all water incl depr storage

		double diff = q->Drc*_dt + WaterVolin->Drc - WaterVolall->Drc - Qn->Drc*_dt;
		if (InfilMethod == INFIL_NONE)
		{
			WaterVolall->Drc = q->Drc*_dt + WaterVolin->Drc - Qn->Drc*_dt;
			InfilVolKinWave->Drc = diff;
		}
		else
			InfilVolKinWave->Drc = diff;
		//infiltrated volume is sum of incoming fluxes+volume before - outgoing flux - volume after

		if (SwitchBuffers && BufferVol->Drc > 0)
		{
			InfilVolKinWave->Drc = 0;
			WaterVolall->Drc = 0;
		}
		// ass long as the cell has a buffer and it is not full there
		// is not infil and no normal watervol


		if (SwitchErosion)
		{
			Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc, DEP->Drc);
			// correct for very high concentrations, 850 after Govers et al
			// recalc sediment volume
		}
	}
}
//---------------------------------------------------------------------------
