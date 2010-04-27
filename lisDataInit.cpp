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
TMMap *TWorld::ReadMap(QString name)
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
TMMap *TWorld::ReadMapMask(cTMap *Mask, QString name)
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
//---------------------------------------------------------------------------
TMMap *TWorld::InitMask(QString name)
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

    //msleep(100);
    //emit debug(_M->PathName);

    return(_M);

}
//---------------------------------------------------------------------------
void TWorld::GetInputData(void)
{
//  LDD = ReadMap(getvaluename("ldd"));
//  InitMask(LDD);
  LDD = InitMask(getvaluename("ldd"));
  Grad = ReadMap(getvaluename("grad"));
  Outlet = ReadMap(getvaluename("outlet"));
  RainZone = ReadMap(getvaluename("id"));
  N = ReadMap(getvaluename("manning"));
  RR = ReadMap(getvaluename("RR"));
  PlantHeight = ReadMap(getvaluename("CH"));
  LAI = ReadMap(getvaluename("lai"));
  Cover = ReadMap(getvaluename("cover"));

  ThetaS1 = ReadMap(getvaluename("ThetaS1"));
  ThetaI1 = ReadMap(getvaluename("ThetaI1"));
  Psi1 = ReadMap(getvaluename("Psi1"));
  Ksat1 = ReadMap(getvaluename("Ksat1"));
  SoilDepth1 = ReadMap(getvaluename("SoilDep1"));
  if (SwitchInfilCrust)
  {
    CrustFraction = ReadMap(getvaluename("crustfrc"));
    KsatCrust = ReadMap(getvaluename("ksatcrst"));
  }
  else
    CrustFraction = NewMap(0);
  if (SwitchInfilCompact)
  {
    CompactFraction = ReadMap(getvaluename("compfrc"));
    KsatCompact = ReadMap(getvaluename("ksatcomp"));
  }
  else
    CompactFraction = NewMap(0);

  GrassPresent = NewMap(0);
  if (SwitchInfilGrass)
  {
      KsatGrass = ReadMap(getvaluename("ksatgras"));
      GrassWidthDX = ReadMap(getvaluename("grasswidth"));
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
      ThetaS2 = ReadMap(getvaluename("ThetaS2"));
      ThetaI2 = ReadMap(getvaluename("ThetaI2"));
      Psi2 = ReadMap(getvaluename("Psi2"));
      Ksat2 = ReadMap(getvaluename("Ksat2"));
      SoilDepth2 = ReadMap(getvaluename("SoilDep2"));
  }
  StoneFraction  = ReadMap(getvaluename("stonefrc"));
 // WheelWidth  = ReadMap(getvaluename("wheelwidth"));
  RoadWidthDX  = ReadMap(getvaluename("road"));

  if (SwitchErosion)
  {
    Cohesion = ReadMap(getvaluename("coh"));
    RootCohesion = ReadMap(getvaluename("cohadd"));
    AggrStab = ReadMap(getvaluename("AggrStab"));
    D50 = ReadMap(getvaluename("D50"));
  }

  if (SwitchIncludeChannel)
  {
    LDDChannel = InitMaskChannel(getvaluename("lddchan"));
    //ReadMapMask(getvaluename("lddchan"));
    ChannelWidth = ReadMapMask(LDDChannel, getvaluename("chanwidth"));
    ChannelSide = ReadMapMask(LDDChannel, getvaluename("chanside"));
    ChannelGrad = ReadMapMask(LDDChannel, getvaluename("changrad"));
    ChannelGrad->calcV(0.001, ADD);
    ChannelN = ReadMapMask(LDDChannel, getvaluename("chanman"));
    ChannelCohesion = ReadMapMask(LDDChannel, getvaluename("chancoh"));

    ChannelGrad->cover(0);
    ChannelSide->cover(0);
    ChannelWidth->cover(0);
    ChannelN->cover(0);
  }

  PointMap = ReadMap(getvaluename("outpoint"));

}
//---------------------------------------------------------------------------
void TWorld::IntializeData(void)
{

  //TO DO add units and descriptions

      //totals for mass balance
      Qtot = 0;
      MB = 0;
      InfilTot = 0;
      InfilKWTot = 0;
      IntercTot = 0;
      WaterVolTot = 0;
      RainTot = 0;
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

      tm = NewMap(0); // temp map for aux calculations

      //terrain maps
      Grad->calcV(0.001, ADD);
      DX = NewMap(0);
      FOR_ROW_COL_MV
      {
          DX->Drc = _dx/cos(Grad->Drc);
      }
      WheelWidthDX = NewMap(0);
      SoilWidthDX = NewMap(0);
      GullyWidthDX = NewMap(0);
      StoneWidthDX = NewMap(0);
      MDS = NewMap(0);
      FOR_ROW_COL_MV
      {
         double RRmm = 10* RR->Drc;
         MDS->Drc = max(0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*Grad->Drc*100);
         MDS->Drc /= 1000; // convert to m
      }

      // rainfall and interception maps
      Rain = NewMap(0);
      RainCum = NewMap(0);
      RainNet = NewMap(0);
      RainIntensity = NewMap(0);
      RainM3 = NewMap(0);
      CStor = NewMap(0);
      Interc = NewMap(0);

      // infiltration maps
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
      Soilwater->calc2(ThetaI1, SoilDepth1, MUL);
      if (SwitchTwoLayer)
          Soilwater2->calc2(ThetaI2, SoilDepth2, MUL);
      InfilMethod = getvalueint("Infil Method");

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
      Qoutflow = NewMap(0);
      q = NewMap(0);
      R = NewMap(0);
      Perim = NewMap(0);
      WaterVol = NewMap(0);
      WaterVolRunoff = NewMap(0);

      // calibration
      ksatCalibration = getvaluedouble("Ksat calibration");
      nCalibration = getvaluedouble("N calibration");
      ChnCalibration = getvaluedouble("Channel Ksat calibration");
      ChKsatCalibration = getvaluedouble("Channel N calibration");
      SplashDelivery = getvaluedouble("Splash Delivery Ratio");

      CanopyStorage = NewMap(0); //m
      FOR_ROW_COL_MV
      {
         CanopyStorage->Drc = 0.935+0.498*LAI->Drc-0.00575*pow(LAI->Drc, 2);
         CanopyStorage->Drc *= 0.001; // to m
      }

      // erosion maps
      Qs = NewMap(0);
      Qsn = NewMap(0);
      Qsoutflow = NewMap(0);
      DETSplash = NewMap(0);
      DETFlow = NewMap(0);
      DEP = NewMap(0);
      SedVol = NewMap(0);

      if (SwitchErosion)
      {
        TC = NewMap(0);
        Conc = NewMap(0);
        CG = NewMap(0);
        DG = NewMap(0);
        SettlingVelocity = NewMap(0);
        CohesionSoil = NewMap(0);
        Y = NewMap(0);

        FOR_ROW_COL_MV
        {
           CG->Drc = pow((D50->Drc+5)/0.32, -0.6);
           DG->Drc = pow((D50->Drc+5)/300, 0.25);
           SettlingVelocity->Drc = 2*(2650-1000)*9.80*pow(D50->Drc/2000000, 2)/(9*0.001);
           CohesionSoil->Drc = (1.0-Cover->Drc)*Cohesion->Drc + Cover->Drc*RootCohesion->Drc;
  //         CohesionSoil->Drc = Cohesion->Drc + RootCohesion->Drc;
           Y->Drc = min(1.0, 1.0/(0.89+0.56*CohesionSoil->Drc));
        }
      }

      // channel maps
      ChannelWidthUpDX = NewMap(0);
      ChannelWaterVol = NewMap(0);
      ChannelQoutflow = NewMap(0);
      RunoffVolinToChannel = NewMap(0);
      SedToChannel = NewMap(0);
      ChannelQsoutflow = NewMap(0);

      if (SwitchIncludeChannel)
      {
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
