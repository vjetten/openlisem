/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * Infiltration: all 1 and 2 layer infiltration functions: G&A, S&P, ksat
 */

#include "model.h"
#include "swatre_g.h"

//NOTE fact and fpot have a unit of m (not m/s)

#define tiny 1e-8

//---------------------------------------------------------------------------
// DOESN'T WORK YET
void TWorld::InfilSwatre(void)
{
	fact->copy(WH);

	SwatreStep(SwatreSoilModel, WH, fpot);
	// WH and fpot done in swatrestep
	FOR_ROW_COL_MV
			fact->Drc = (fact->Drc - WH->Drc);

	if (SwitchInfilCrust)
	{
		tm->copy(WH);
		tma->fill(0);
		SwatreStep(SwatreSoilModelCrust, tm, tma);
		FOR_ROW_COL_MV
		{
			//tm = WHcrust and tma = fpot crust
			fact->Drc = (tm->Drc - WH->Drc)*CrustFraction->Drc + fact->Drc*(1-CrustFraction->Drc);
			WH->Drc = tm->Drc*CrustFraction->Drc + WH->Drc*(1-CrustFraction->Drc);
			fpot->Drc = tma->Drc*CrustFraction->Drc + fpot->Drc*(1-CrustFraction->Drc);
		}
	}

	if (SwitchInfilCompact)
	{
		tm->copy(WH);
		tma->fill(0);
		SwatreStep(SwatreSoilModelCompact, tm, tma);
		FOR_ROW_COL_MV
		{
			fact->Drc = (tm->Drc - WH->Drc)*CompactFraction->Drc + fact->Drc*(1-CompactFraction->Drc);
			WH->Drc = tm->Drc*CompactFraction->Drc + WH->Drc*(1-CompactFraction->Drc);
			fpot->Drc = tma->Drc*CompactFraction->Drc + fpot->Drc*(1-CompactFraction->Drc);
		}
	}

	if (SwitchInfilGrass)
	{
		factgr->copy(WHGrass);
		SwatreStep(SwatreSoilModelGrass, WHGrass, fpotgr);
		FOR_ROW_COL_MV
				factgr->Drc = (factgr->Drc - WHGrass->Drc);
	}
}
//---------------------------------------------------------------------------
// DOESN'T WORK YET
void TWorld::InfilMorelSeytoux1(void)
{
	FOR_ROW_COL_MV
	{
		double fact1;
		double Ks = Ksateff->Drc/3600000.0;  //in m/s
		double fwh = WH->Drc; // in m, in WH is old WH + net rainfall
		double rt = fwh/_dt; // m/s pseudo rainfall, all water/dt = rate
		double A = 0, B, tp;
		double tt=time-BeginTime;

		if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
			Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
		// if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2

		B = (fwh+Psi1->Drc*0.01)*(ThetaS1->Drc-ThetaI1->Drc); //m
		tp = max(0, Ks*B/(rt*rt - Ks*rt)); //sec

		if (Ks < rt )
		{
			if (tp < tt+_dt)
				A = pow(B+Fcum->Drc,2)/(2*Ks*B*pow(rt/Ks-1,2));
			fpot->Drc = 0.5*sqrt(2*Ks*pow(B+Fcum->Drc,2)/B)*(1/sqrt(tt-tp+A))+Ks;
		}
		else
			fpot->Drc = Ks;
		// potential infiltration in m
		// psi input map is in cm so multiply 0.01 for m

		fact1 = min(fpot->Drc, fwh);
		// actual infil in m, cannot have more infil than water on the surface

		fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc);

	}
}
//---------------------------------------------------------------------------
// Solution Eurosem v2 manual page 10, Morgan et al 1998
void TWorld::InfilSmithParlange1(void)
{
	FOR_ROW_COL_MV
	{
		double fact1;
		double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
		double fwh = WH->Drc; // in m, in WH is old WH + net rainfall
		double Cdexp, B = (fwh + Psi1->Drc*0.01)*max(ThetaS1->Drc-ThetaI1->Drc, tiny);
		// psi input map is in cm so multiply 0.01 for m

		if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
			Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
		// smallest of the two ksats in two laters blocks the flow

		Cdexp = exp(Fcum->Drc/B);
		fpot->Drc = Ks*Cdexp/(Cdexp-1);
		// potential infiltration in m

		fact1 = min(fpot->Drc, fwh);
		// actual infil in m, cannot have more infil than water on the surface

		fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc);
		// adjust fact->Drc for twolayer, impermeable etc

		if(SwitchInfilGrass)
		{
			fwh = WHGrass->Drc;
			Ks = KsatGrass->Drc*_dt/3600000.0;  //in m

			if (SwitchTwoLayer && L1gr->Drc > SoilDepth1->Drc - tiny)
				Ks = min(KsatGrass->Drc, Ksat2->Drc)*_dt/3600000.0;

			B = (fwh + Psi1->Drc*0.01)*max(ThetaS1->Drc-ThetaI1->Drc, tiny);

			Cdexp = exp(Fcum->Drc/B);
			fpotgr->Drc = Ks*Cdexp/(Cdexp-1);
			// potential infiltration in m

			fact1 = min(fpotgr->Drc, fwh);

			factgr->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc);
		}
	}
}
//---------------------------------------------------------------------------
// Solution Kutilek and Nielsen 2004 pag 138
void TWorld::InfilGreenAmpt1(void)
{
	FOR_ROW_COL_MV
	{
		double fact1;
		double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
		double fwh = WH->Drc; // in m, WH is old WH + net rainfall

		if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
			Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;
		// if wetting front > layer 1 than ksat is determined by smallest ksat1 and ksat2

		fpot->Drc = Ks*(1+(Psi1->Drc*0.01+fwh)/(L1->Drc+L2->Drc));
		// potential infiltration in m, Darcy : Q = K * (dh/dz + 1)
		// L1 initialized at 1e-10, psi in cm so multiply 0.01 for m

		fact1 = min(fpot->Drc, fwh);
		// actual infil in m, cannot have more infil than water on the surface

		fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc);
		// adjust fact for twolayer, impermeable etc

		if(SwitchInfilGrass)
		{
			fwh = WHGrass->Drc; // in m, WH is old WH + net rainfall
			Ks = KsatGrass->Drc*_dt/3600000.0;  //in m

			if (SwitchTwoLayer && L1gr->Drc > SoilDepth1->Drc - tiny)
				Ks = min(KsatGrass->Drc, Ksat2->Drc)*_dt/3600000.0;

			fpotgr->Drc = Ks*(1+(Psi1->Drc*0.01+fwh)/(L1gr->Drc+L2->Drc));

			fact1 = min(fpotgr->Drc, fwh);

			factgr->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc);
		}
	}
}
//---------------------------------------------------------------------------
// Direct subtraction of Ksat, added for testing purposes!
void TWorld::InfilKsat(void)
{
	FOR_ROW_COL_MV
	{
		double fact1;
		double Ks = Ksateff->Drc*_dt/3600000.0;  //in m
		double fwh = WH->Drc; // in m, WH is old WH + net rainfall

		if (SwitchTwoLayer && L1->Drc > SoilDepth1->Drc - tiny)
			Ks = min(Ksateff->Drc, Ksat2->Drc)*_dt/3600000.0;

		fpot->Drc = Ks;
		// potential infil equals Ksat, unit is m

		fact1 = min(fpot->Drc, fwh);
		// actual infil in m, cannot have more infil than water on the surface

		fact->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1->Drc, &L2->Drc);
		// adjust fact for twolayer, impermeable etc

		if(SwitchInfilGrass)
		{
			fwh = WHGrass->Drc; // in m, WH is old WH + net rainfall
			Ks = KsatGrass->Drc*_dt/3600000.0;  //in m

			if (SwitchTwoLayer && L1gr->Drc > SoilDepth1->Drc - tiny)
				Ks = min(KsatGrass->Drc, Ksat2->Drc)*_dt/3600000.0;

			fpotgr->Drc = Ks;

			fact1 = min(fpotgr->Drc, fwh);

			factgr->Drc = IncreaseInfiltrationDepth(r, c, fact1, &L1gr->Drc, &L2gr->Drc);
		}
	}
}
//---------------------------------------------------------------------------
// function to increase wetting front and deal with 2nd layer
// returns actual infiltration
// this function is called form all infiltration functions
double TWorld::IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p)
{
	double dL1, dL2;
	double L1, L2;
	L1 = *L1p;
	L2 = *L2p;

	dL1 = fact/max(tiny, ThetaS1->Drc-ThetaI1->Drc);
	// increase in depth is infiltration/available porespace

	if (SwitchImpermeable && !SwitchTwoLayer)
	{
		// reaches bottom in this timestep, fact is remaining water
		if (L1+dL1 > SoilDepth1->Drc && L1 < SoilDepth1->Drc)
		{
			dL1 = max(0, SoilDepth1->Drc - L1);
			// if L1 = soildepth1, increase will be 0
			fact = dL1 * (ThetaS1->Drc-ThetaI1->Drc);
		}
		// soil is full, fact = 0, no more increase
		if (L1 > SoilDepth1->Drc - tiny)
			fact = 0;
	}

	if (!SwitchTwoLayer && (L1+dL1 > SoilDepth1->Drc - tiny))
		dL1 = SoilDepth1->Drc - L1;
	L1 += dL1;
	// increase infiltration depth L1  = fact/avail pore space
	// if soil is full and not impermeable, infiltration continues as steady state

	//infiltration can go on in the second layer & first layer is full
	if (SwitchTwoLayer && (L1 > SoilDepth1->Drc - tiny))
	{
		dL2 = fact/max(tiny, (ThetaS2->Drc-ThetaI2->Drc));
		// water enters into second layer
		if (SwitchImpermeable)
		{
			// reaches bottom in this timestep, fact is remaining water
			if (L2+dL2 > SoilDepth2->Drc && L2 < SoilDepth2->Drc)
			{
				dL2 = max(0, SoilDepth2->Drc - L2);
				// if L2 = soildepth2, increase will be 0
				fact = dL2 * (ThetaS2->Drc-ThetaI2->Drc);
			}
			// soil is full, fact = 0, no more increase
			if (L2 > SoilDepth2->Drc - tiny)
				fact = 0;
		}
		// soil is full, but not impermeable, infiltration continues as steady state
		if (L2+dL2 > SoilDepth2->Drc - tiny)
			dL2 = SoilDepth2->Drc - L2;
		L2 += dL2;
		// increase infiltration depth L2  = fact/avail pore space
	}

	*L1p = (REAL8)L1;
	*L2p = (REAL8)L2;

	return fact; //m
}
//---------------------------------------------------------------------------
/*
 * main function:
 * add rain to WH, calc effective Ksat
 */
void TWorld::Infiltration(void)
{
	FOR_ROW_COL_MV
	{
		WH->Drc += RainNet->Drc;
		// add net to water rainfall on soil surface (in m)

		if (GrassPresent->Drc > 0)
			WHGrass->Drc += RainNet->Drc;
		// net rainfall on grass strips, infil is calculated separately for grassstrips

		if (RoadWidthDX->Drc > 0)
			WHroad->Drc += Rainc->Drc;
		// assume no interception and infiltration on roads, gross rainfall

		InfilVol->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
		// vol before infil

		if (InfilMethod != INFIL_SWATRE)
		{
			Ksateff->Drc = Ksat1->Drc*(1-CrustFraction->Drc-CompactFraction->Drc);
			// avg ksat of "normal" surface with crusting and compaction fraction, fractions are 0 when there is none

			if (SwitchInfilCrust)
				Ksateff->Drc += KsatCrust->Drc*CrustFraction->Drc;
			if (SwitchInfilCompact)
				Ksateff->Drc += KsatCompact->Drc*CompactFraction->Drc;
			// adjust effective infil for crusting and compaction

			Ksateff->Drc *= ksatCalibration;
			// apply runfile/iface calibration factor
			if (SwitchBuffers && !SwitchSedtrap)
				if(BufferID->Drc > 0)
					Ksateff->Drc = 0;
			//VJ 1000608 no infil in buffers, , but sedtrap can have infil
		}
	} //row col

	//select an infiltration type
	switch (InfilMethod)
	{
	case INFIL_NONE : fact->fill(0); fpot->fill(0);break;
	case INFIL_SWATRE : InfilSwatre(); break;
	case INFIL_HOLTAN : break;
	case INFIL_GREENAMPT : InfilGreenAmpt1(); break;
	case INFIL_GREENAMPT2 : InfilGreenAmpt1(); break;
	case INFIL_KSAT : InfilKsat(); break;
	case INFIL_MOREL : InfilMorelSeytoux1(); break; //TODO: DOESN'T WORK YET
	case INFIL_SMITH : InfilSmithParlange1(); break;
	}
	// these function result in an actual infiltration "fact" (in m)
	// and potential infiltration "fpot" (in m)
	// each function deals with grass strips as a separate infiltration process

	FOR_ROW_COL_MV
	{
		if (SwitchBuffers && !SwitchSedtrap)
			if(BufferID->Drc > 0 && BufferVol->Drc > 0)
				WH->Drc = 0;
		//VJ 100608 no infil in buffers until it is full
		//TODO NOTE CORRECT FOR RAINFALL IF WH IS 0

		if (InfilMethod != INFIL_SWATRE)
		{
			WH->Drc -= fact->Drc;
			if (WH->Drc < 0) // in case of rounding of errors
			{
				fact->Drc += WH->Drc;
				WH->Drc = 0;
			}
			// subtract fact->Drc from WH, cannot be more than WH

			Fcum->Drc += fact->Drc;
			// cumulative infil in m

			if (GrassPresent->Drc > 0)
			{
				WHGrass->Drc -= factgr->Drc;
				if (WHGrass->Drc < 0) // in case of rounding of errors
				{
					factgr->Drc += WHGrass->Drc;
					WHGrass->Drc = 0;
				}
				Fcumgr->Drc += factgr->Drc;
			}
			// calculate and correct water height on grass strips
		}

		if (GrassPresent->Drc > 0)
			WH->Drc = WH->Drc*(1-GrassFraction->Drc) + GrassFraction->Drc * WHGrass->Drc;
		// average water height if grasstrip present

		FSurplus->Drc = min(0, fact->Drc - fpot->Drc);
		// negative surplus of infiltration in m for kinematic wave in m

		if (GrassPresent->Drc > 0)
			FSurplus->Drc = min(0, factgr->Drc - fpotgr->Drc);
		// if grasstrip present use grasstrip surplus as entire surplus

		InfilVol->Drc -= DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
		//	InfilVol->Drc = max(0, InfilVol->Drc);
		// infil volume is WH before - water after
	}
}
//---------------------------------------------------------------------------
void TWorld::SoilWater()
{
	if (!SwitchSoilwater)
		return;

	FOR_ROW_COL_MV
	{
		double Percolation = 0;
		if (!SwitchImpermeable || SwitchTwoLayer)
		{
			if (SwitchInfilCompact)
				Percolation = Ksateff->Drc * pow(ThetaI1->Drc/ThetaS1->Drc, 9);
			else
				Percolation = Ksat1->Drc * pow(ThetaI1->Drc/ThetaS1->Drc, 9);
		}
		//NOTE it is assumed that crusting does NOT but compaction DOES affect the percolation of the first layer

		Soilwater->Drc = L1->Drc*ThetaS1->Drc + (SoilDepth1->Drc - L1->Drc)*ThetaI1->Drc;
		Soilwater->Drc = max(0, Soilwater->Drc - Percolation);
		Soilwater->Drc = min(Soilwater->Drc, SoilDepth1->Drc*ThetaS1->Drc);

		ThetaI1->Drc = Soilwater->Drc/SoilDepth1->Drc;

		// if evaporation decrease L1

		if (SwitchTwoLayer)
		{
			double Percolation2 = 0;
			if (!SwitchImpermeable)
				Percolation2 = Ksat2->Drc * pow(ThetaI2->Drc/ThetaS2->Drc, 9);

			Soilwater2->Drc = L2->Drc*ThetaS2->Drc + (SoilDepth2->Drc - L2->Drc)*ThetaI2->Drc;
			Soilwater2->Drc = max(0, Soilwater2->Drc - Percolation2);
			Soilwater2->Drc = min(Soilwater2->Drc, SoilDepth2->Drc*ThetaS2->Drc);
			ThetaI2->Drc = Soilwater2->Drc/SoilDepth2->Drc;
		}

	}
}

