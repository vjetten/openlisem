/*---------------------------------------------------------------------------
project: openLISEM
name: lisDatainit.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisDatainit.cpp:
- newmap, readmap, destoy data
- get inputdata, init data, init boolean options (SwitchSomething)
---------------------------------------------------------------------------*/

#include "model.h"
#include <qstring.h>

//---------------------------------------------------------------------------
void TWorld::InitMapList(void)
{

  maplistnr = 0;
  for (int i = 0; i < 200; i++)
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

    _M->PathName = /*inputdir + */name;
    bool res = _M->LoadFromFile();
    if (!res)
    {
      ErrorString = "Cannot find map " +_M->PathName;
      throw 1;

    }

    for (int r = 0; r < nrRows; r++)
     for (int c = 0; c < nrCols; c++)
     if (!IS_MV_REAL4(&Mask->Drc) && IS_MV_REAL4(&_M->Drc))
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
TMMap *TWorld::ReadMap(cTMap *Mask, QString name)
{

    TMMap *_M = new TMMap();

    _M->PathName = /*inputdir + */name;
    bool res = _M->LoadFromFile();
    if (!res)
    {
      ErrorString = "Cannot find map " +_M->PathName;
      throw 1;
    }

    for (int r = 0; r < nrRows; r++)
     for (int c = 0; c < nrCols; c++)
     if (!IS_MV_REAL4(&Mask->Drc) && IS_MV_REAL4(&_M->Drc))
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
  baseNameMap = getvaluename("grad");
  Grad = ReadMap(LDD,getvaluename("grad"));  // must be SINE of the slope angle !!!
  Outlet = ReadMap(LDD,getvaluename("outlet"));

  FOR_ROW_COL_MV
  {
	  if (Outlet->Drc == 1 && LDD->Drc != 5)
	  {
	      ErrorString = "Main outlet gridcell does not coincide with pit in LDD";
	      throw 1;
	  }
  }
  RainZone = ReadMap(LDD,getvaluename("id"));
  N = ReadMap(LDD,getvaluename("manning"));
  RR = ReadMap(LDD,getvaluename("RR"));
  PlantHeight = ReadMap(LDD,getvaluename("CH"));
  LAI = ReadMap(LDD,getvaluename("lai"));
  Cover = ReadMap(LDD,getvaluename("cover"));

  ThetaS1 = ReadMap(LDD,getvaluename("ThetaS1"));
  ThetaI1 = ReadMap(LDD,getvaluename("ThetaI1"));
  Psi1 = ReadMap(LDD,getvaluename("Psi1"));
  Ksat1 = ReadMap(LDD,getvaluename("Ksat1"));
  SoilDepth1 = ReadMap(LDD,getvaluename("SoilDep1"));
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

  GrassPresent = NewMap(0);
  if (SwitchInfilGrass)
  {
      KsatGrass = ReadMap(LDD,getvaluename("ksatgras"));
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

  if (SwitchTwoLayer)
  {
      ThetaS2 = ReadMap(LDD,getvaluename("ThetaS2"));
      ThetaI2 = ReadMap(LDD,getvaluename("ThetaI2"));
      Psi2 = ReadMap(LDD,getvaluename("Psi2"));
      Ksat2 = ReadMap(LDD,getvaluename("Ksat2"));
      SoilDepth2 = ReadMap(LDD,getvaluename("SoilDep2"));
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

}
//---------------------------------------------------------------------------
void TWorld::IntializeData(void)
{

  //TO DO add units and descriptions

      //totals for mass balance
      Qtot = 0;
      Qtotmm = 0;
      Qpeak = 0;
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
      DetTotSplash = 0;
      DetTotFlow = 0;
      DepTot = 0;
      DetTot = 0;
      DepTot = 0;
      SoilLossTot = 0;
      SedVolTot = 0;
      MBs = 0;
      ChannelVolTot=0;
      ChannelSedTot = 0;
      ChannelDepTot = 0;
      ChannelDetTot = 0;

      tm = NewMap(0); // temp map for aux calculations
      nrCells = Mask->MapTotal();

      //terrain maps
//      Grad->calcV(0.0001, ADD);
      DX = NewMap(0);
      FOR_ROW_COL_MV
      {
    	  Grad->Drc = max(0.0001, Grad->Drc);
          DX->Drc = _dx/cos(asin(Grad->Drc));
      }
      CatchmentArea = DX->MapTotal() * _dx;

      WheelWidthDX = NewMap(0);
      SoilWidthDX = NewMap(0);
      GullyWidthDX = NewMap(0);
      StoneWidthDX = NewMap(0);
      MDS = NewMap(0);
      FOR_ROW_COL_MV
      {
         double RRmm = 10* RR->Drc;
         MDS->Drc = max(0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*tan(asin(Grad->Drc))*100);
         MDS->Drc /= 1000; // convert to m
      }

      // rainfall and interception maps
      Rain = NewMap(0);
      RainCum = NewMap(0);
      RainNet = NewMap(0);
      LeafDrain = NewMap(0);
      RainIntensity = NewMap(0);
      RainM3 = NewMap(0);
      CStor = NewMap(0);
      Interc = NewMap(0);

      // infiltration maps
      InfilMethod = getvalueint("Infil Method");
      if (InfilMethod == INFIL_GREENAMPT2)
    	  SwitchTwoLayer = true;
      InfilVolKinWave = NewMap(0);
      InfilVol = NewMap(0);
      Fcum = NewMap(1e-10);
      L1 = NewMap(1e-10);
      L2 = NewMap(1e-10);
      FSurplus = NewMap(0);
      fact = NewMap(0);
      fpot = NewMap(0);
      Fcumgr = NewMap(1e-10);
      L1gr = NewMap(1e-10);
      L2gr = NewMap(1e-10);
      factgr = NewMap(0);
      fpotgr = NewMap(0);
      Ksateff = NewMap(0);
      Soilwater = NewMap(0);
      Soilwater2 = NewMap(0);
      Soilwater->calc2(ThetaI1, SoilDepth1, MUL);
      if (SwitchTwoLayer)
      {
          Soilwater2->calc2(ThetaI2, SoilDepth2, MUL);
      }
      // runoff maps
      WH = NewMap(0);
      WHrunoff = NewMap(0);
      WHstore = NewMap(0);
      WHroad = NewMap(0);
      WHinf = NewMap(0);
      WHGrass = NewMap(0);
      FlowWidth = NewMap(0);
      fpa = NewMap(0);
      V = NewMap(0);
      Alpha = NewMap(0);
      Q = NewMap(0);
      Qn = NewMap(0);
      Qoutput = NewMap(0);
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

      QString output = getvaluename("CheckOutputMaps");

      N->calcV(nCalibration, MUL);
      if (SwitchIncludeChannel)
      {
          ChannelN->calcV(ChnCalibration, MUL);
          if (SwitchChannelInfil)
             ChannelKsat->calcV(ChKsatCalibration, MUL);

      }

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
      SedVol = NewMap(0);

      TC = NewMap(0);
      Conc = NewMap(0);
      CG = NewMap(0);
      DG = NewMap(0);
      SettlingVelocity = NewMap(0);
      CohesionSoil = NewMap(0);
      Y = NewMap(0);

      if (SwitchErosion)
      {
        FOR_ROW_COL_MV
        {
           CG->Drc = pow((D50->Drc+5)/0.32, -0.6);
           DG->Drc = pow((D50->Drc+5)/300, 0.25);
           SettlingVelocity->Drc = 2*(2650-1000)*9.80*pow(D50->Drc/2000000, 2)/(9*0.001);
           //           CohesionSoil->Drc = Cohesion->Drc + Cover->Drc*RootCohesion->Drc;
           // soil cohesion everywhere, plantcohesion only where plants
           CohesionSoil->Drc = Cohesion->Drc + RootCohesion->Drc;
           Y->Drc = min(1.0, 1.0/(0.89+0.56*CohesionSoil->Drc));
        }
      }

      // channel maps
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
      ChannelSedVol = NewMap(0);
      ChannelConc = NewMap(0);
      ChannelTC = NewMap(0);
      ChannelY = NewMap(0);

      if (SwitchIncludeChannel)
      {

           ChannelWidthUpDX->copy(ChannelWidth);
           ChannelWidthUpDX->cover(0);
           FOR_ROW_COL_MV
           {
                if (!IS_MV_REAL4(&LDDChannel->Data[r][c]))
                   ChannelMask->Drc = 1;
                else
                   SET_MV_REAL4(&ChannelMask->Data[r][c]);
           }
           FOR_ROW_COL_MV_CH
           {
               ChannelDX->Drc = _dx/ChannelGrad->Drc;
               ChannelY->Drc = min(1.0, 1.0/(0.89+0.56*ChannelCohesion->Drc));
           }
      }

      TotalDetMap = NewMap(0);
      TotalDepMap = NewMap(0);
      TotalSoillossMap = NewMap(0);
}
//---------------------------------------------------------------------------
void TWorld::IntializeOptions(void)
{

  //dirs and names
  resultDir.clear();
  inputDir.clear();
  outflowFileName.clear();
  totalErosionFileName.clear();
  totalDepositionFileName.clear();
  totalSoillossFileName.clear();
  outflowFileName.clear();
  rainFileName.clear();
  rainFileDir.clear();
  snowmeltFileName.clear();
  snowmeltFileDir.clear();
  tableDir.clear();
  resultFileName.clear();

    SwitchHardsurface =
    SwatreInitialized =
    SwitchInfilGA2 =
    SwitchCrustPresent =
    SwitchWheelPresent =
    SwitchCompactPresent =
    SwitchIncludeChannel =
    SwitchChannelBaseflow =
    startbaseflowincrease =
    SwitchChannelInfil =
    SwitchAllinChannel =
    SwitchErosion =
    SwitchAltErosion =
    SwitchSimpleDepression =
    SwitchBuffers =
    SwitchSedtrap =
    SwitchSnowmelt =
    SwitchRunoffPerM =
    SwitchInfilCompact =
    SwitchInfilCrust =
    SwitchInfilGrass =
    SwitchImpermeable =
    SwitchDumphead =
    SwitchGeometricMean =
    SwitchWheelAsChannel =
    SwitchMulticlass =
    SwitchNutrients =
    SwitchGullies =
    SwitchGullyEqualWD =
    SwitchGullyInfil =
    SwitchGullyInit =
    SwitchOutputTimeStep =
    SwitchOutputTimeUser =
    SwitchMapoutRunoff =
    SwitchMapoutConc =
    SwitchMapoutWH =
    SwitchMapoutWHC =
    SwitchMapoutTC =
    SwitchMapoutEros =
    SwitchMapoutDepo =
    SwitchMapoutV =
    SwitchMapoutInf =
    SwitchMapoutSs =
    SwitchMapoutChvol =
    SwitchWritePCRnames =
    SwitchWritePCRtimeplot =
    SwitchNoErosionOutlet =
    SwitchDrainage =
    SwitchPestout =
    SwitchSeparateOutput =
    SwitchSOBEKOutput =
    SwitchInterceptionLAI =
    SwitchTwoLayer =
    SwitchSimpleSedKinWave =
    SwitchSOBEKoutput =
    SwitchPCRoutput =
    SwitchSoilwater = false;
}
//---------------------------------------------------------------------------
/*
void TWorld::InitMask(cTMap *M)
{
    Mask = new TMMap();
    Mask->_MakeMap(M, 1.0);
    maplist[maplistnr].m = Mask;
    maplistnr++;
    _dx = Mask->MH.cellSizeX*1.0000000;
    nrRows = Mask->nrRows;
    nrCols = Mask->nrCols;
}
*/