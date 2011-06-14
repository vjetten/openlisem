
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisErosion.cpp
  \brief Flow and splash detachment functions for slopes and channels

functions: \n
- double TWorld::MaxConcentration(double watvol, double sedvol, double dep) \n
- void TWorld::SplashDetachment(void) \n
- void TWorld::FlowDetachment(void) \n
- void TWorld::ChannelFlowDetachment(void) \n
 */

#include "model.h"


//---------------------------------------------------------------------------
double TWorld::MaxConcentration(double watvol, double sedvol, double dep)
{
	double conc = (watvol > 0 ? sedvol/watvol : 1000);
	if (conc > MAXCONC)
	{
		dep += min(0, MAXCONC*watvol - sedvol);
		conc = MAXCONC;
	}
	sedvol = conc*watvol;

	return(conc);
}
//---------------------------------------------------------------------------
void TWorld::SplashDetachment(void)
{
	if (!SwitchErosion)
		return;

	FOR_ROW_COL_MV
	{
		double b, strength, DetDT1, DetDT2, DetLD1, DetLD2;
		double g_to_kg = 0.001;

		double Int = Rain->Drc * 3600/_dt * 1000;
		// intensity in mm/h, Rain is in m

		double KE_DT = 28.3*(1-(0.52*exp(-0.042*Int)));
		// kin energy in J/m2/mm, Van DIjk general equation 2002
      KE_DT = 7.60*pow(Int, 0.22);
      // kin energy in J/m2/mm, cape verde Sanchez et al in prep 2011



      /* TODO: allow many different equations here, interface choice */
		/* equation in LISEM, based on Eurosem, Morgan 1998
       if (Int > 1)
           KE_DT = 8.95+8.44*log10(Int);
       else
    	   KE_DT = 0;
		 */

		double directrain = (1-Cover->Drc)*Rainc->Drc * 1000; //
		// rainfall between plants in mm

		double KE_LD = max(15.3*sqrt(PlantHeight->Drc)-5.87, 0);
		// kin energy in J/m2/mm
//		double throughfall = LeafDrain->Drc * (1-StemflowFraction) * 1000;
      double throughfall = LeafDrain->Drc * 1000;
      // leaf drip in mm is already calculated with plant covern in interception function
      // VJ 110206 stemflow is also accounted for

		double WH0 = exp(-1.48*WH->Drc*1000);
		// water buffer effect on surface, WH in mm in this empirical equation from Torri ?

		if (AggrStab > 0)
		{
			strength = 2.82/AggrStab->Drc; b = 2.96;
		}
		else
		{
			strength = 0.1033/CohesionSoil->Drc; b = 3.58;
		}

		// Between plants, directrain is already with 1-cover
		DetDT1 = g_to_kg * fpa->Drc*(strength*KE_DT*WH0+b) * directrain;
		//ponded areas, kg/m2/mm * mm = kg/m2
		DetDT2 = g_to_kg * (1-fpa->Drc)*(strength*KE_DT+b) * directrain * SplashDelivery;
		//dry areas, kg/m2/mm * mm = kg/m2

		// Under plants, throughfall is already with cover
		DetLD1 = g_to_kg * fpa->Drc*(strength*KE_LD*WH0+b) * throughfall;
		//ponded areas, kg/m2/mm * mm = kg/m2
		DetLD2 = g_to_kg * (1-fpa->Drc)*(strength*KE_LD+b) * throughfall * SplashDelivery;
		//dry areas, kg/m2/mm * mm = kg/m2

		DETSplash->Drc = DetLD1 + DetLD2 + DetDT1 + DetDT2;
		// Total splash kg/m2

		// Deal with all exceptions:

		DETSplash->Drc *= (SoilWidthDX->Drc*DX->Drc);
		// kg/cell, only splash over soilwidth, not roads and channels

		DETSplash->Drc = (1-StoneFraction->Drc) * DETSplash->Drc;
		// no splash on stone surfaces

      if (GrassFraction->Drc > 0)
			DETSplash->Drc = (1-GrassFraction->Drc) * DETSplash->Drc;
		// no splash on grass strips

		if (SwitchBuffers && !SwitchSedtrap)
			if(BufferID->Drc > 0)
				DETSplash->Drc = 0;
		// no splash in buffers, but sedtrap can have splash

      if (SwitchHardsurface)
         DETSplash->Drc = (1-HardSurface->Drc)*DETSplash->Drc;
		// no splash on hard surfaces
		DETSplash->Drc = (1-Snowcover->Drc)*DETSplash->Drc;
		// no splash on snow deck
	}
}
//---------------------------------------------------------------------------
void TWorld::FlowDetachment(void)
{
	if (!SwitchErosion)
		return;

	FOR_ROW_COL_MV
	{
		//### Calc transport capacity
		double omega = 100*V->Drc*Grad->Drc;
		// V in cm/s in this formula assuming grad is SINE
		double omegacrit = 0.4;
		// critical unit streampower in cm/s
		TC->Drc = min(MAXCONC, 2650 * CG->Drc * pow(max(0, omega - omegacrit), DG->Drc));
		// not more than 2650*0.32 = 848 kg/m3

		//### Add splash to sediment
		Sed->Drc += DETSplash->Drc;
		// add splash to sed volume

		//### calc concentration and net transport capacity
		DEP->Drc = 0;
		// init deposition for this timestep
		Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc, DEP->Drc);
		// limit sed concentration to max

		double maxTC = max(TC->Drc - Conc->Drc,0);
		double minTC = min(TC->Drc - Conc->Drc,0);
		// unit kg/m3

		//### detachment
		double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
		// detachment can only come from soil, not roads (so do not use flowwidth)
		// units s * m/s * m * m = m3

		DETFlow->Drc = Y->Drc * maxTC * TransportFactor;
		// unit = kg/m3 * m3 = kg
		DETFlow->Drc = min(DETFlow->Drc, maxTC * WaterVolall->Drc);
		// cannot have more detachment than remaining capacity in flow

      if (GrassFraction->Drc > 0)
			DETFlow->Drc = (1-GrassFraction->Drc) * DETFlow->Drc;
		// no flow detachment on grass strips

		// Detachment edxceptions:
		DETFlow->Drc = (1-StoneFraction->Drc) * DETFlow->Drc ;
		// no flow detachment on stony surfaces

      if (SwitchHardsurface)
         DETFlow->Drc = (1-HardSurface->Drc) * DETFlow->Drc ;
		// no flow detachment on hard surfaces

		//DETFlow->Drc = (1-Snowcover->Drc) * DETFlow->Drc ;
      /* TODO: CHECK THIS no flow detachment on snow */
        //is there erosion and sedimentation under the snowdeck?

		//### deposition
		if (WH->Drc > MIN_HEIGHT)
			TransportFactor = (1-exp(-_dt*SettlingVelocity->Drc/WH->Drc)) * WaterVolall->Drc;
		else
			TransportFactor = WaterVolall->Drc;
		//   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
		// deposition can occur on roads and on soil (so use flowwidth)

		double deposition = minTC * TransportFactor;
		// max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

		//deposition = max(deposition, minTC * WaterVol->Drc);
		// cannot be more than sediment above capacity
		deposition = max(deposition, -Sed->Drc);
		// cannot have more depo than sediment present
      /* TODO what about this: which one to choose */

		//deposition = (1-Snowcover->Drc) * deposition;
      /* TODO: TRUE??? no deposition on snow */

      if (GrassFraction->Drc > 0)
			deposition = -Sed->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
		// generate 100% deposition on grassstrips
		//? bit tricky, maximizes effect on grassstrips ?

		//### sediment balance
		DEP->Drc += deposition;
		Sed->Drc += DETFlow->Drc;
		Sed->Drc += deposition;

		Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc, DEP->Drc);
		// limit concentration to 850 and throw rest in deposition
	}
}
//---------------------------------------------------------------------------
void TWorld::ChannelFlowDetachment(void)
{
	if (!SwitchErosion)
		return;

	FOR_ROW_COL_MV_CH
	{
		double omega = 100*ChannelV->Drc*ChannelGrad->Drc;
		double omegacrit = 0.4;

		ChannelTC->Drc = min(MAXCONC, 2650 * CG->Drc * pow(max(0, omega - omegacrit), DG->Drc));
		// Channel transport capacity, ot more than 2650*0.32 = 848 kg/m3

		ChannelDep->Drc = 0;
		ChannelDetFlow->Drc = 0;

		ChannelSed->Drc += SedToChannel->Drc;
		// add sed flow into channel from slope

		ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc, ChannelDep->Drc);
		// set conc to max and add surplus sed to ChannelDep

		double maxTC = max(ChannelTC->Drc - ChannelConc->Drc,0);
		double minTC = min(ChannelTC->Drc - ChannelConc->Drc,0);
		// basic TC deficit, posiutive if detachment, negative if depositioin, units kg/m3

		double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * ChannelWidthUpDX->Drc;
		// units s * m/s * m * m = m3
		ChannelDetFlow->Drc = ChannelY->Drc * maxTC * TransportFactor;
		// unit kg/m3 * m3 = kg
		ChannelDetFlow->Drc = min(ChannelDetFlow->Drc, maxTC * ChannelWaterVol->Drc);
		// cannot have more detachment than remaining capacity in flow

		//or:? ChannelDetFlow->Drc = ChannelY->Drc * maxTC * Q->Drc*_dt;
		// simple erosion formula

		double deposition = minTC * TransportFactor;
		// max deposition in kg/s  < 0
		deposition = max(deposition, -ChannelSed->Drc);
		// cannot be more than sediment above capacity
		//or:?     deposition = max(deposition, minTC * ChannelWaterVol->Drc);

		ChannelDep->Drc += deposition;
		ChannelSed->Drc += deposition;
		ChannelSed->Drc += ChannelDetFlow->Drc;

		ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc, ChannelDep->Drc);
	}
}
//---------------------------------------------------------------------------

