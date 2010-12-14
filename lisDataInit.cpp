/*---------------------------------------------------------------------------
project: openLISEM
name: lisDatainit.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisDatainit.cpp:
- newmap, readmap, destoy data
 - get inputdata, init data, init boolean options (SwitchSomething)
---------------------------------------------------------------------------*/


#include "model.h"
#include "swatre_g.h"
#include <qstring.h>

//char SwatreErrorString[256];
//---------------------------------------------------------------------------
void TWorld::InitMapList(void)
{

	maplistnr = 0;
	for (int i = 0; i < NUMNAMES; i++)
	{
		maplist[i].m = NULL;
	}

}
//---------------------------------------------------------------------------
TMMap *TWorld::NewMap(double value)
{
	TMMap *_M = new TMMap();

	_M->_MakeMap(Mask, value);

	if (_M)
	{
		maplist[maplistnr].m = _M;
		maplistnr++;
	}

	return(_M);
}
//---------------------------------------------------------------------------
TMMap *TWorld::ReadMapMask(QString name)
{

	TMMap *_M = new TMMap();

	_M->PathName = name;
	bool res = _M->LoadFromFile();
	if (!res)
	{
		ErrorString = "Cannot find map " +_M->PathName;
		throw 1;

	}

	for (int r = 0; r < nrRows; r++)
		for (int c = 0; c < nrCols; c++)
			if (!IS_MV_REAL8(&Mask->Drc) && IS_MV_REAL8(&_M->Drc))
			{
            ErrorString =
                  QString("Missing value in map %1 where LDD has a value,\n at row=%2 and col=%3").arg(name).arg(r).arg(c);

            throw 1;
         }

	if (_M)
	{
		maplist[maplistnr].m = _M;
		maplistnr++;
	}

	//msleep(100);
	//emit debug(_M->PathName);

	return(_M);

}
//---------------------------------------------------------------------------
TMMap *TWorld::ReadMap(cTMap *Mask, QString name)
{

	TMMap *_M = new TMMap();

	_M->PathName = /*inputdir + */name;
	//	DEBUG(_M->PathName);

	bool res = _M->LoadFromFile();
	if (!res)
	{
		ErrorString = "Cannot find map " +_M->PathName;
		throw 1;
	}

	for (int r = 0; r < nrRows; r++)
		for (int c = 0; c < nrCols; c++)
			if (!IS_MV_REAL8(&Mask->Drc) && IS_MV_REAL8(&_M->Drc))
			{
            QString sr, sc;
            sr.setNum(r); sc.setNum(c);
            ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+name;
            throw 1;
         }

	if (_M)
	{
		maplist[maplistnr].m = _M;
		maplistnr++;
	}

	//msleep(100);
	//emit debug(_M->PathName);

	return(_M);

}
//---------------------------------------------------------------------------
void TWorld::DestroyData(void)
{
	for (int i = 0; i < maplistnr; i++)
	{
		if (maplist[i].m != NULL)
		{
			maplist[i].m->KillMap();
			maplist[i].m = NULL;
		}
	}
	if (nrrainfallseries > 1)
	{
		for (int r=0; r < nrrainfallseries; r++)
			delete[] RainfallSeries[r];
		delete[] RainfallSeries;
	}

	if (InfilMethod == INFIL_SWATRE && initSwatreStructure)
	{
		FreeSwatreInfo();
		if (SwatreSoilModel)
			CloseSwatre(SwatreSoilModel);
		if (SwatreSoilModelCrust)
			CloseSwatre(SwatreSoilModelCrust);
		if (SwatreSoilModelCompact)
			CloseSwatre(SwatreSoilModelCompact);
		if (SwatreSoilModelGrass)
			CloseSwatre(SwatreSoilModelGrass);
	}
}
//---------------------------------------------------------------------------
TMMap *TWorld::InitMask(QString name)
{
	// read map and make a mask map

	TMMap *_M = new TMMap();

	_M->PathName = /*inputdir + */name;
	bool res = _M->LoadFromFile();
	if (!res)
	{
		ErrorString = "Cannot find map " +_M->PathName;
		throw 1;
	}

	if (_M)
	{
		maplist[maplistnr].m = _M;
		maplistnr++;
	}

	Mask = new TMMap();
	Mask->_MakeMap(_M, 1.0);
	maplist[maplistnr].m = Mask;
	maplistnr++;
	_dx = Mask->MH.cellSizeX*1.0000000;
	nrRows = Mask->nrRows;
	nrCols = Mask->nrCols;

	//msleep(100);
	//emit debug(_M->PathName);
	return(_M);

}
//---------------------------------------------------------------------------
TMMap *TWorld::InitMaskChannel(QString name)
{

	TMMap *_M = new TMMap();

	_M->PathName = /*inputdir + */name;
	bool res = _M->LoadFromFile();
	if (!res)
	{
		ErrorString = "Cannot find map " +_M->PathName;
		throw 1;
	}

	if (_M)
	{
		maplist[maplistnr].m = _M;
		maplistnr++;
	}

	MaskChannel = new TMMap();
	MaskChannel->_MakeMap(_M, 1.0);
	maplist[maplistnr].m = MaskChannel;
	maplistnr++;

	return(_M);

}
//---------------------------------------------------------------------------
void TWorld::GetInputData(void)
{
   LDD = InitMask(getvaluename("ldd"));
	// LDD is also mask and reference file, everthiung has to fit LDD
	// chanels use channel LDD as mask

	Grad = ReadMap(LDD,getvaluename("grad"));  // must be SINE of the slope angle !!!
	Outlet = ReadMap(LDD,getvaluename("outlet"));
   Outlet->cover(0);
   // fill outlet with zero, some users have MV where no outlet

	FOR_ROW_COL_MV
	{
		if (Outlet->Drc == 1)
		{
			if (LDD->Drc != 5)
			{
				ErrorString = "Main outlet gridcell does not coincide with pit in LDD";
				throw 1;
			}
			else
			{
				c_outlet = c;
				r_outlet = r;
			}
		}
	}

	if (SwitchRainfall)
   {
      RainZone = ReadMap(LDD,getvaluename("id"));
   }
	Snowcover = NewMap(0);
	if (SwitchSnowmelt)
	{
		SnowmeltZone = ReadMap(LDD,getvaluename("SnowID"));
		FOR_ROW_COL_MV
		{
			Snowcover->Drc = (SnowmeltZone->Drc == 0 ? 0 : 1.0);
		}
	}
	N = ReadMap(LDD,getvaluename("manning"));
	RR = ReadMap(LDD,getvaluename("RR"));
	PlantHeight = ReadMap(LDD,getvaluename("CH"));
	LAI = ReadMap(LDD,getvaluename("lai"));
	Cover = ReadMap(LDD,getvaluename("cover"));

	GrassPresent = NewMap(0);
	if (SwitchInfilGrass)
	{
		GrassWidthDX = ReadMap(LDD,getvaluename("grasswidth"));
		GrassFraction->copy(GrassWidthDX);
		GrassFraction->calcV(_dx, DIV);
		StripN = getvaluedouble("Grassstrip Mannings n");
		FOR_ROW_COL_MV
		{
			if (GrassWidthDX > 0)
			{
				GrassPresent->Drc = 1;
				N->Drc = N->Drc*(1-GrassFraction->Drc)+StripN*GrassFraction->Drc;
			}
		}
	}
	else
		GrassFraction = NewMap(0);

	if(InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE)
	{
      Ksat1 = ReadMap(LDD,getvaluename("ksat1"));
      SoilDepth1 = ReadMap(LDD,getvaluename("soildep1"));
      SoilDepth1->calcV(1000, DIV);
      //VJ 101213 convert from mm to m
      // can be zero for outcrops
		FOR_ROW_COL_MV
		{
			bool check = false;
			if (SoilDepth1->Drc < 0)
				check = true;
			if (check)
			{
				ErrorString = QString("SoilDepth1 values < 0 at row %1, col %2").arg(r).arg(c);
				throw 1;
			}
		}

		if(InfilMethod != INFIL_KSAT)
		{
         ThetaS1 = ReadMap(LDD,getvaluename("thetas1"));
         ThetaI1 = ReadMap(LDD,getvaluename("thetai1"));
         Psi1 = ReadMap(LDD,getvaluename("psi1"));
			if (SwitchTwoLayer)
			{
            ThetaS2 = ReadMap(LDD,getvaluename("thetaS2"));
            ThetaI2 = ReadMap(LDD,getvaluename("thetaI2"));
            Psi2 = ReadMap(LDD,getvaluename("psi2"));
            Ksat2 = ReadMap(LDD,getvaluename("ksat2"));
            SoilDepth2 = ReadMap(LDD,getvaluename("soilDep2"));
            SoilDepth2->calcV(1000, DIV);
            //VJ 101213 convert from mm to m
            FOR_ROW_COL_MV
				{
					bool check = false;
               if (SoilDepth1->Drc < 0)
						check = true;
					if (check)
					{
						ErrorString = QString("SoilDepth1 values <= 0 at row %1, col %2").arg(r).arg(c);
						throw 1;
					}
				}
			}
		}
		if (SwitchInfilCrust)
		{
			CrustFraction = ReadMap(LDD,getvaluename("crustfrc"));
			KsatCrust = ReadMap(LDD,getvaluename("ksatcrst"));
		}
		else
			CrustFraction = NewMap(0);
		if (SwitchInfilCompact)
		{
			CompactFraction = ReadMap(LDD,getvaluename("compfrc"));
			KsatCompact = ReadMap(LDD,getvaluename("ksatcomp"));
		}
		else
			CompactFraction = NewMap(0);
	}
	if (InfilMethod == INFIL_SWATRE)
	{
		ProfileID = ReadMap(LDD,getvaluename("profmap"));

		if (SwitchInfilGrass)
			ProfileIDGrass = ReadMap(LDD,getvaluename("profgrass"));

		if (SwitchInfilCrust)
		{
			CrustFraction = ReadMap(LDD,getvaluename("crustfrc"));
			ProfileIDCrust = ReadMap(LDD,getvaluename("profcrust"));
		}
		else
			CrustFraction = NewMap(0);

		if (SwitchInfilCompact)
		{
			CompactFraction = ReadMap(LDD,getvaluename("compfrc"));
			ProfileIDCompact = ReadMap(LDD,getvaluename("profcomp"));
		}
		else
			CompactFraction = NewMap(0);

		int res = ReadSwatreInput(SwatreTableName, SwatreTableDir);
		if (res)
			throw res;
	}

	StoneFraction  = ReadMap(LDD,getvaluename("stonefrc"));
	// WheelWidth  = ReadMap(LDD,getvaluename("wheelwidth"));
	RoadWidthDX  = ReadMap(LDD,getvaluename("road"));
	HardSurface = ReadMap(LDD,getvaluename("hardsurf"));

	if (SwitchErosion)
	{
		Cohesion = ReadMap(LDD,getvaluename("coh"));
		RootCohesion = ReadMap(LDD,getvaluename("cohadd"));
		AggrStab = ReadMap(LDD,getvaluename("AggrStab"));
		D50 = ReadMap(LDD,getvaluename("D50"));
	}

	if (SwitchIncludeChannel)
	{
		LDDChannel = InitMaskChannel(getvaluename("lddchan"));

		ChannelWidth = ReadMap(LDDChannel, getvaluename("chanwidth"));
		ChannelSide = ReadMap(LDDChannel, getvaluename("chanside"));
		ChannelGrad = ReadMap(LDDChannel, getvaluename("changrad"));
		ChannelGrad->calcV(0.001, ADD);
		ChannelN = ReadMap(LDDChannel, getvaluename("chanman"));
		ChannelCohesion = ReadMap(LDDChannel, getvaluename("chancoh"));
		if (SwitchChannelInfil)
			ChannelKsat = ReadMap(LDDChannel, getvaluename("chanksat"));

		ChannelGrad->cover(0);
		ChannelSide->cover(0);
		ChannelWidth->cover(0);
		ChannelN->cover(0);
	}

	PointMap = ReadMap(LDD,getvaluename("outpoint"));

	if (SwitchBuffers || SwitchSedtrap)
	{
		BufferID = ReadMap(LDD,getvaluename("bufferID"));
		BufferVol = ReadMap(LDD,getvaluename("bufferVolume"));
		BulkDens = getvaluedouble("Sediment bulk density");
		// also sed trap use bufffervol to calculate the max sed store

		FOR_ROW_COL_MV
		{
			if (SwitchBuffers && BufferID->Drc > 0)
			{
				Grad->Drc = 0.001;
				RR->Drc = 0.01;
				N->Drc = 0.25;
				// very arbitrary!!!
				Cover->Drc = 0;
				if (SwitchIncludeChannel && ChannelGrad->Drc > 0)
				{
					ChannelGrad->Drc = 0.001;
					ChannelN->Drc = 0.25;
				}
			}
			// adjust soil and surface [roperties for buffercells, not sed traps
		}
	}
	else
	{
		BufferID = NewMap(0);
		BufferVol = NewMap(0);
	}
}
//---------------------------------------------------------------------------
void TWorld::IntializeData(void)
{

	//TO DO add units and descriptions

	//totals for mass balance and file/screen output
	Qtot = 0;
	QtotOutlet = 0;
	Qtotmm = 0;
	Qpeak = 0;
	QpeakTime = 0;
	MB = 0;
	InfilTot = 0;
	InfilTotmm = 0;
	InfilKWTot = 0;
	IntercTot = 0;
	IntercTotmm = 0;
	WaterVolTot = 0;
	WaterVolTotmm = 0;
	RainTot = 0;
	RainTotmm = 0;
	Rainpeak = 0;
	RainpeakTime = 0;
   RainAvgmm = 0;
   SnowAvgmm = 0;
	SnowTot = 0;
	SnowTotmm = 0;
	Snowpeak = 0;
	SnowpeakTime = 0;
	DetSplashTot = 0;
	DetFlowTot = 0;
	DepTot = 0;
	DetTot = 0;
	DepTot = 0;
	SoilLossTot = 0;
	SoilLossTotOutlet = 0;
	SedTot = 0;
	MBs = 0;
	ChannelVolTot=0;
	ChannelSedTot = 0;
	ChannelDepTot = 0;
	ChannelDetTot = 0;
	BufferVolTot = 0;
	BufferVolTotInit = 0;
	BufferSedTot = 0;

	tm = NewMap(0); // temp map for aux calculations
	tma = NewMap(0); // temp map for aux calculations
	nrCells = Mask->MapTotal();

	//terrain maps
	DX = NewMap(0);
	CellArea = NewMap(0);
	FOR_ROW_COL_MV
	{
		Grad->Drc = max(0.0001, Grad->Drc);
		DX->Drc = _dx/cos(asin(Grad->Drc));
		CellArea->Drc = DX->Drc * _dx;
	}
	CatchmentArea = CellArea->MapTotal();

	WheelWidth = NewMap(0);
	WheelWidthDX = NewMap(0);
	SoilWidthDX = NewMap(0);
	GullyWidthDX = NewMap(0);
	MDS = NewMap(0);
	FOR_ROW_COL_MV
	{
		double RRmm = 10* RR->Drc;
		MDS->Drc = max(0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*tan(asin(Grad->Drc))*100);
		MDS->Drc /= 1000; // convert to m
	}

	// rainfall and interception maps
	Rain = NewMap(0);
	Rainc = NewMap(0);
	RainCum = NewMap(0);
	RainNet = NewMap(0);
	LeafDrain = NewMap(0);
	RainIntensity = NewMap(0);
	RainM3 = NewMap(0);
	CStor = NewMap(0);
	Interc = NewMap(0);
	
   Snowmelt = NewMap(0);
	Snowmeltc = NewMap(0);
	SnowmeltCum = NewMap(0);
	//SnowmeltNet = NewMap(0);

	// infiltration maps
	InfilVolKinWave = NewMap(0);
	InfilVol = NewMap(0);
	InfilVolCum = NewMap(0);
	fact = NewMap(0);
	fpot = NewMap(0);
	factgr = NewMap(0);
	fpotgr = NewMap(0);
	Ksateff = NewMap(0);
	FSurplus = NewMap(0);

	if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
	{
		Fcum = NewMap(1e-10);
		L1 = NewMap(1e-10);
		L2 = NewMap(1e-10);
		Fcumgr = NewMap(1e-10);
		L1gr = NewMap(1e-10);
		L2gr = NewMap(1e-10);
      if (InfilMethod != INFIL_KSAT)
      {
         Soilwater = NewMap(0);
         Soilwater2 = NewMap(0);
         Soilwater->calc2(ThetaI1, SoilDepth1, MUL);
         if (SwitchTwoLayer)
         {
            Soilwater2->calc2(ThetaI2, SoilDepth2, MUL);
         }
      }
	}

	// runoff maps
	WH = NewMap(0);
	WHrunoff = NewMap(0);
	WHrunoffCum = NewMap(0);
	WHstore = NewMap(0);
	WHroad = NewMap(0);
	WHGrass = NewMap(0);
	FlowWidth = NewMap(0);
	fpa = NewMap(0);
	V = NewMap(0);
	Alpha = NewMap(0);
	Q = NewMap(0);
	Qn = NewMap(0);
	Qoutput = NewMap(0);
	Qsoutput = NewMap(0);
	Qoutflow = NewMap(0); // value of Qn*dt in pits only
	q = NewMap(0);
	R = NewMap(0);
	Perim = NewMap(0);
	WaterVolrunoff = NewMap(0);
	WaterVolin = NewMap(0);
	WaterVolall = NewMap(0);

	// calibration
	ksatCalibration = getvaluedouble("Ksat calibration");
	nCalibration = getvaluedouble("N calibration");
	ChnCalibration = getvaluedouble("Channel Ksat calibration");
	ChKsatCalibration = getvaluedouble("Channel N calibration");
	SplashDelivery = getvaluedouble("Splash Delivery Ratio");
	StemflowFraction = getvaluedouble("Stemflow fraction");
	N->calcV(nCalibration, MUL);

	if (SwitchIncludeChannel)
	{
		ChannelN->calcV(ChnCalibration, MUL);
		if (SwitchChannelInfil)
			ChannelKsat->calcV(ChKsatCalibration, MUL);

	}

	
   SwatreSoilModel = NULL;
	SwatreSoilModelCrust = NULL;
	SwatreSoilModelCompact = NULL;
	SwatreSoilModelGrass = NULL;
	if (InfilMethod == INFIL_SWATRE)
	{
		double precision = 5.0;
		// note "5" is a precision factor dewtermining next timestep, set to 5 in old lisem
		SwatreSoilModel = InitSwatre(ProfileID, initheadName, swatreDT, precision,
                                   ksatCalibration, SwitchGeometric, SwitchImpermeable);
		if (SwatreSoilModel == NULL)
			throw 3;

		if (SwitchInfilCrust)
		{
			SwatreSoilModelCrust = InitSwatre(ProfileIDCrust, initheadName, swatreDT, precision,
                                           ksatCalibration, SwitchGeometric, SwitchImpermeable);
			if (SwatreSoilModelCrust == NULL)
				throw 3;
		}
		if (SwitchInfilCompact)
		{
			SwatreSoilModelCompact = InitSwatre(ProfileIDCompact, initheadName, swatreDT, precision,
                                             ksatCalibration, SwitchGeometric, SwitchImpermeable);
			if (SwatreSoilModelCompact == NULL)
				throw 3;
		}
		if (SwitchInfilGrass)
		{
			SwatreSoilModelGrass = InitSwatre(ProfileIDGrass, initheadName, swatreDT, precision,
                                           ksatCalibration, SwitchGeometric, SwitchImpermeable);
			if (SwatreSoilModelGrass == NULL)
				throw 3;
		}
      initSwatreStructure = true;
      // flag: structure is created and can be destroyed in function destroydata
	}

	InterceptionLAIType = getvalueint("Canopy storage equation");
	if (InterceptionLAIType == 8)
		SwitchInterceptionLAI = false;
	if (SwitchInterceptionLAI)
	{
		CanopyStorage = NewMap(0); //m
		InterceptionLAIType = getvalueint("Canopy storage equation");
		FOR_ROW_COL_MV
		{
			switch (InterceptionLAIType)
			{
			case 0: CanopyStorage->Drc = 0.935+0.498*LAI->Drc-0.00575*(LAI->Drc * LAI->Drc);break;
			case 1: CanopyStorage->Drc = 0.2331 * LAI->Drc; break;
			case 2: CanopyStorage->Drc = 0.3165 * LAI->Drc; break;
			case 3: CanopyStorage->Drc = 1.46 * pow(LAI->Drc,0.56); break;
			case 4: CanopyStorage->Drc = 0.0918 * pow(LAI->Drc,1.04); break;
			case 5: CanopyStorage->Drc = 0.2856 * LAI->Drc; break;
			case 6: CanopyStorage->Drc = 0.1713 * LAI->Drc; break;
			case 7: CanopyStorage->Drc = 0.59 * pow(LAI->Drc,0.88); break;
			}
		}
	}
	else
		CanopyStorage = ReadMap(LDD,getvaluename("smax"));
	CanopyStorage->calcV(0.001, MUL); // to m
	//NOTE: LAI is still needed for canopy openness, can be circumvented with cover

	// erosion maps
	Qs = NewMap(0);
	Qsn = NewMap(0);
	Qsoutflow = NewMap(0);
	DETSplash = NewMap(0);
	DETFlow = NewMap(0);
	DEP = NewMap(0);
	Sed = NewMap(0);
	TC = NewMap(0);
	Conc = NewMap(0);
	CG = NewMap(0);
	DG = NewMap(0);
	SettlingVelocity = NewMap(0);
	CohesionSoil = NewMap(0);
	Y = NewMap(0);

	TotalDetMap = NewMap(0);
	TotalDepMap = NewMap(0);
	TotalSoillossMap = NewMap(0);
	TotalSed = NewMap(0);
	TotalWatervol = NewMap(0);
	TotalConc = NewMap(0);


	if (SwitchErosion)
	{
		FOR_ROW_COL_MV
		{
			CG->Drc = pow((D50->Drc+5)/0.32, -0.6);
			DG->Drc = pow((D50->Drc+5)/300, 0.25);
			SettlingVelocity->Drc = 2*(2650-1000)*9.80*pow(D50->Drc/2000000, 2)/(9*0.001);
			CohesionSoil->Drc = Cohesion->Drc + Cover->Drc*RootCohesion->Drc;
			// soil cohesion everywhere, plantcohesion only where plants
			Y->Drc = min(1.0, 1.0/(0.89+0.56*CohesionSoil->Drc));
		}
	}

	// channel maps that must be there even if cxhannel is switched off
	SedToChannel = NewMap(0);
	ChannelWidthUpDX = NewMap(0);
	ChannelWaterVol = NewMap(0);
	ChannelQoutflow = NewMap(0);
	RunoffVolinToChannel = NewMap(0);
	ChannelQsoutflow = NewMap(0);
	ChannelQ = NewMap(0);
	ChannelQn = NewMap(0);
	ChannelQs = NewMap(0);
	ChannelQsn = NewMap(0);
	ChannelV = NewMap(0);
	ChannelWH = NewMap(0);
	Channelq = NewMap(0);
	ChannelAlpha = NewMap(0);
	ChannelMask = NewMap(0);
	ChannelDX = NewMap(0);
	ChannelDetFlow = NewMap(0);
	ChannelDep = NewMap(0);
	ChannelSed = NewMap(0);
	ChannelConc = NewMap(0);
	ChannelTC = NewMap(0);
	ChannelY = NewMap(0);

	if (SwitchIncludeChannel)
	{

		ChannelWidthUpDX->copy(ChannelWidth);
		ChannelWidthUpDX->cover(0);
		FOR_ROW_COL_MV
		{
			if (!IS_MV_REAL8(&LDDChannel->Data[r][c]))
				ChannelMask->Drc = 1;
			else
				SET_MV_REAL8(&ChannelMask->Data[r][c]);
		}
		FOR_ROW_COL_MV_CH
		{
			ChannelDX->Drc = _dx/ChannelGrad->Drc;
			ChannelY->Drc = min(1.0, 1.0/(0.89+0.56*ChannelCohesion->Drc));
		}
	}

	//VJ 100514 buffer and sedtrap maps
   /** TODO how calculate max sed store with only sed traps? */
	// use slope of cell:        | /
	//                           |/
	// then max store is _dx/cos = DX*height fence * bulk dens?
	if (SwitchBuffers)
	{
		BufferVolInit = NewMap(0);
		ChannelBufferVolInit = NewMap(0);

		if (SwitchIncludeChannel)
		{
			ChannelBufferVol = NewMap(0);
			FOR_ROW_COL_MV_CH
               if (BufferID->Drc > 0)
			{
				ChannelBufferVol->Drc = BufferVol->Drc;
				BufferVol->Drc = 0;
			}
			//split buffers in channel buffers and slope buffers
			// in "ToCHannel" all flow in a buffer is dumped in the channel
			ChannelBufferVolInit->copy(ChannelBufferVol);
			// copy initial max volume of buffers in channels
		}
		BufferVolInit->copy(BufferVol);
		// copy initial max volume of remaining buffers on slopes
		BufferVolTotInit = BufferVolInit->MapTotal() + ChannelBufferVolInit->MapTotal();
		// sum up total initial volume available in buffers
		//BufferVol->fill(0);
		//ChannelBufferVol->fill(0);
		// rset to zero to fill up

	}


	BufferSed = NewMap(0);
	if (SwitchBuffers || SwitchSedtrap)
	{
		BufferSedInit = NewMap(0);
		ChannelBufferSedInit = NewMap(0);

		BufferSed->calc2V(BufferVol, BulkDens, MUL);
		//NOTE: buffer sed vol is maximum store in kg and will decrease while it
		// fills up. It is assumed that the sedimented part contains a pore volume
		// that can contain water, like  loose soil. Thsi is determined by the bulkdensity
		if (SwitchIncludeChannel)
		{
			ChannelBufferSed = NewMap(0);
			FOR_ROW_COL_MV_CH
               if (BufferID->Drc > 0)
			{
				ChannelBufferSed->Drc = BufferSed->Drc;
				BufferSed->Drc = 0;
			}
			//split buffers in channel buffers and slope buffers
			// in "ToCHannel" all flow in a buffer is dumped in the channel
			ChannelBufferSedInit->copy(ChannelBufferSed);
			// copy initial max volume of buffers in channels
		}
		BufferSedInit->copy(BufferSed);
		// copy initial max volume of remaining buffers on slopes
		BufferSedTotInit = BufferSedInit->MapTotal() + ChannelBufferSedInit->MapTotal();
		// sum up total initial volume available in buffers
		//BufferSed->fill(0);
		//ChannelBufferSed->fill(0);
		// rset to zero to fill up
	}
}
//---------------------------------------------------------------------------
void TWorld::IntializeOptions(void)
{
	nrrainfallseries = 0;

	//dirs and names
	resultDir.clear();
	inputDir.clear();
	outflowFileName = QString("totals.txt");//.clear();
	totalErosionFileName = QString("erososion.map");//.clear();
	totalDepositionFileName = QString("deposition.map");//.clear();
	totalSoillossFileName = QString("soilloss.map");//.clear();
	outflowFileName = QString("hydrohtaph.csv");//.clear();
	rainFileName.clear();
	rainFileDir.clear();
	snowmeltFileName.clear();
	snowmeltFileDir.clear();
	SwatreTableDir.clear();
	SwatreTableName = QString("profile.inp");//.clear();
	resultFileName.clear();

	SwitchHardsurface = false;
	SwatreInitialized = false;
	SwitchInfilGA2 = false;
	SwitchCrustPresent = false;
	SwitchWheelPresent = false;
	SwitchCompactPresent = false;
	SwitchIncludeChannel = false;
	SwitchChannelBaseflow = false;
	startbaseflowincrease = false;
	SwitchChannelInfil = false;
	SwitchAllinChannel = false;
	SwitchErosion = false;
	SwitchAltErosion = false;
	SwitchSimpleDepression = false;
	SwitchBuffers = false;
	SwitchSedtrap = false;
	SwitchSnowmelt = false;
	SwitchRunoffPerM = false;
	SwitchInfilCompact = false;
	SwitchInfilCrust = false;
	SwitchInfilGrass = false;
	SwitchImpermeable = false;
	SwitchDumphead = false;
	SwitchWheelAsChannel = false;
	SwitchMulticlass = false;
	SwitchNutrients = false;
	SwitchGullies = false;
	SwitchGullyEqualWD = false;
	SwitchGullyInfil = false;
	SwitchGullyInit = false;
	SwitchOutputTimeStep = false;
	SwitchOutputTimeUser = false;
	SwitchMapoutRunoff = false;
	SwitchMapoutConc = false;
	SwitchMapoutWH = false;
	SwitchMapoutWHC = false;
	SwitchMapoutTC = false;
	SwitchMapoutEros = false;
	SwitchMapoutDepo = false;
	SwitchMapoutV = false;
	SwitchMapoutInf = false;
	SwitchMapoutSs = false;
	SwitchMapoutChvol = false;
	SwitchWritePCRnames = false;
	SwitchWritePCRtimeplot = false;
	SwitchNoErosionOutlet = false;
	SwitchDrainage = false;
	SwitchPestout = false;
	SwitchSeparateOutput = false;
	SwitchSOBEKOutput = false;
	SwitchInterceptionLAI = false;
	SwitchTwoLayer = false;
	SwitchSimpleSedKinWave = false;
	SwitchSOBEKoutput = false;
	SwitchPCRoutput = false;
	SwitchSoilwater = false;
	SwitchGeometric = true;

	SwitchWriteHeaders = true; // write headers in output files in first timestep

   initSwatreStructure = false;
   // check to flag when swatre 3D structure is created, needed to clean up data
}
//---------------------------------------------------------------------------

