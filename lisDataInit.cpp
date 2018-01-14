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
 \file lisDatainit.cpp
 \brief All data handling functions, read, parse, initialize, make and delete the map structures etc.

 functions: \n
 - cTMap TWorld::NewMap(double value) Make a map filled with a value and LDD as mask \n
 - cTMap TWorld::ReadMap(cTMap *Mask, QString name) Make a map and fill it with data from a file\n
 - void TWorld::DestroyData(void) Destroy the list with all maps and pointers, no need to do this one by one.\n
 - cTMap TWorld::InitMask(QString name) Make and read the reference map which is the LDD \n
 - cTMap TWorld::InitMaskChannel(QString name) Make and read the channel map (channel LDD)\n
 - cTMap TWorld::InitMaskTiledrain(QString name) Make and read the tiell drain map (tile LDD)\n
 - void TWorld::InitTiledrains(void) read and intiialize all Tile drain variables and maps \n
 - void TWorld::InitBuffers(void)  read and intiialize all buffer and sediment fence variables and maps \n
 - void TWorld::InitChannel(void)  read and intiialize all channel variables and maps \n
 - void TWorld::GetInputData(void) Read all input maps \n
 - void TWorld::IntializeData(void) Initialize all auxillary maps and variables\n
 - void TWorld::IntializeOptions(void) Initialize all options (Switch...) before the run file.\n
*/

#include <algorithm>
#include <qstring.h>
#include "io.h"
#include "model.h"
#include "operation.h"


#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )



//---------------------------------------------------------------------------
/** \n void TWorld::InitMapList(void)
*  blabla
*/
void TWorld::InitMapList(void)
{

    maplistnr = 0;
    for (int i = 0; i < NUMNAMES; i++)
    {
        maplistCTMap[i].m = NULL;
    }
}
//---------------------------------------------------------------------------
cTMap *TWorld::NewMap(double value)
{
    cTMap *_M = new cTMap();

    _M->MakeMap(LDD, value);
    // changed to LDD instead of Mask

    if (_M)
    {
        maplistCTMap[maplistnr].m = _M;
        maplistnr++;
    }

    return(_M);
}
//---------------------------------------------------------------------------
cTMap *TWorld::ReadMap(cTMap *Mask, QString name)
{

    cTMap *_M = new cTMap(readRaster(/*inputdir + */name));

    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
            if (!pcr::isMV(Mask->Drc) && pcr::isMV(_M->Drc))
            {
                QString sr, sc;
                sr.setNum(r); sc.setNum(c);
                ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+name+".\n \
                        This is within the flow domain (the LDD has a value here).\n \
                        This iusually happens when you have maps of different origin";
                        throw 1;
            }

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    return(_M);

}
//---------------------------------------------------------------------------
void TWorld::DestroyData(void)
{
    DEBUG("clear all maps");
    for (int i = 0; i < maplistnr; i++)
    {
        if (maplistCTMap[i].m != NULL)
        {
            delete maplistCTMap[i].m;
            maplistCTMap[i].m = NULL;
        }
    }
    DEBUG("clear rainfall structure");
    if (nrRainfallseries > 1)
    {
        for (int r=0; r < nrRainfallseries; r++)
            RainfallSeriesM[r].intensity.clear();
        RainfallSeriesM.clear();
    }
    if (nrSnowmeltseries > 1)
    {
        for (int r=0; r < nrSnowmeltseries; r++)
        {
            SnowmeltSeriesM[r].intensity.clear();
        }
        SnowmeltSeriesM.clear();
    }
    DEBUG("kill swatre structure");

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
    DEBUG("kill display data");
    ClearComboMaps();
    ClearHydrographData();


}
//---------------------------------------------------------------------------
cTMap *TWorld::InitMask(QString name)
{
    // read map and make a mask map

    cTMap *_M = new cTMap(readRaster(/*inputdir + */name));

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    _dx = _M->cellSize()*1.0000000;
    _nrRows = _M->nrRows();
    _nrCols = _M->nrCols();

    return(_M);

}
//---------------------------------------------------------------------------
cTMap *TWorld::InitMaskChannel(QString name)
{

    cTMap *_M = new cTMap(readRaster(/*inputdir + */name));

    maplistCTMap[maplistnr].m = _M;
    maplistnr++;

    return(_M);

}
//---------------------------------------------------------------------------
// Make a shaded relief map from the DEM for map display
//WARNING barriers are added here!!!!
void TWorld::InitShade(void)
{
    Shade = NewMap(0);
    //shade=cos?(I)sin?(S)cos(A-D)+sin?(I)cos(S)
    //double MaxDem = DEM->mapMaximum();
    // double MinDem = DEM->mapMinimum();
    double PI=3.14159265;
    FOR_ROW_COL_MV
    {
        double Incl = 20.0/180.0*PI;
        double Decl = 300/180.0*PI;
        double mat[9];
        double dx, dy, aspect;
        for (int i = 0; i < 9; i++)
            mat[i] = DEM->Drc;
        if (r > 0 && r < _nrRows-1 && c > 0 && c < _nrCols-1)
        {
            if(!pcr::isMV(LDD->data[r-1][c-1]))
                mat[0] = DEM->data[r-1][c-1];
            if(!pcr::isMV(LDD->data[r-1][c  ]))
                mat[1] = DEM->data[r-1][c  ];
            if(!pcr::isMV(LDD->data[r-1][c+1]))
                mat[2] = DEM->data[r-1][c+1];
            if(!pcr::isMV(LDD->data[r  ][c-1]))
                mat[3] = DEM->data[r  ][c-1];

            if(!pcr::isMV(LDD->data[r  ][c+1]))
                mat[5] = DEM->data[r  ][c+1];
            if(!pcr::isMV(LDD->data[r+1][c-1]))
                mat[6] = DEM->data[r+1][c-1];
            if(!pcr::isMV(LDD->data[r+1][c  ]))
                mat[7] = DEM->data[r+1][c  ];
            if(!pcr::isMV(LDD->data[r+1][c+1]))
                mat[8] = DEM->data[r+1][c+1];
        }
        dx = (mat[2] + 2*mat[5] + mat[8] - mat[0] -2*mat[3] - mat[6])/(8*_dx);
        dy = (mat[0] + 2*mat[1] + mat[2] - mat[6] -2*mat[7] - mat[8])/(8*_dx);
        if (dy < 0)
            aspect = atan(dx/dy)+2*PI;
        else
            if (dy > 0)
                aspect = atan(dx/dy)+PI;
            else
                aspect = 0;
        //qDebug() << r << c << aspect;
        Shade->Drc = cos(Incl)*Grad->Drc*cos(aspect-Decl) + sin(Incl)*cos(asin(Grad->Drc));
    }
    double MaxV = mapMaximum(*Shade);
    double MinV = mapMinimum(*Shade);

    FOR_ROW_COL_MV
    {
        Shade->Drc = (Shade->Drc-MinV)/(MaxV-MinV);
    }

}
//---------------------------------------------------------------------------
// read and Intiialize all channel variables and maps
void TWorld::InitChannel(void)
{
    // channel vars and maps that must be there even if channel is switched off
    ChannelVolTot = 0;
    ChannelSedTot = 0;
    ChannelDepTot = 0;
    ChannelDetTot = 0;
    FloodSedTot = 0;
    FloodDepTot = 0;
    FloodDetTot = 0;

    SedToChannel = NewMap(0);
    ChannelWidthUpDX = NewMap(0);
    //ChannelQoutflow = NewMap(0);
    RunoffVolinToChannel = NewMap(0);
    //ChannelQsoutflow = NewMap(0);
    ChannelAlpha = NewMap(0);
    ChannelPerimeter = NewMap(0); //VJ 110109 added for channel infil
    ChannelDX = NewMap(0);


    hmx = NewMap(0);
    UVflood = NewMap(0);
    floodTimeStart = NewMap(0);
    floodHmxMax = NewMap(0);
    floodVMax = NewMap(0);
    floodTime = NewMap(0);
    floodWaterVol= NewMap(0);

    dfhmx = NewMap(0);
    dfUV = NewMap(0);
    dfTimeStart = NewMap(0);
    dfHmxMax = NewMap(0);
    dfVMax = NewMap(0);
    dfTime = NewMap(0);

    maxChannelflow = NewMap(0);
    maxChannelWH = NewMap(0);
    ChannelAdj = NewMap(_dx);

    if (SwitchIncludeChannel)
    {
        //## channel maps
        LDDChannel = InitMaskChannel(getvaluename("lddchan"));
        // must be first" LDDChannel is the mask for channels
        ChannelWidth = ReadMap(LDDChannel, getvaluename("chanwidth"));
        //     ChannelWidth->checkMap(LARGER, _dx, "Channel width must be smaller than cell size");
        //ChannelWidth->checkMap(SMALLEREQUAL, 0, "Channel width must be larger than 0 in channel cells");
        //      ChannelWidth->calcValue(0.9*_dx, MIN);
        FOR_ROW_COL_MV_CH
        {
            if (ChannelWidth->Drc <= 0)
            {
                ErrorString = QString("Map %1 contains channel cells with width = 0").arg(getvaluename("chanwidth"));
                throw 1;
            }
        }

        ChannelSide = ReadMap(LDDChannel, getvaluename("chanside"));
        ChannelGrad = ReadMap(LDDChannel, getvaluename("changrad"));
        checkMap(*ChannelGrad, LARGER, 1.0, "Channel Gradient must be SINE of slope angle (not tangent)");
        calcValue(*ChannelGrad, 0.001, MAX);
        ChannelN = ReadMap(LDDChannel, getvaluename("chanman"));
        ChannelCohesion = ReadMap(LDDChannel, getvaluename("chancoh"));
        cover(*ChannelGrad, *LDD, 0);
        cover(*ChannelSide, *LDD, 0);
        cover(*ChannelWidth, *LDD, 0);
        cover(*ChannelN, *LDD, 0);

        calcValue(*ChannelN, ChnCalibration, MUL);
        if (SwitchChannelInfil)
        {
            ChannelKsat = ReadMap(LDDChannel, getvaluename("chanksat"));
            cover(*ChannelKsat, *LDD, 0);
            calcValue(*ChannelKsat, ChKsatCalibration, MUL);
            // store not used?
        }
        copy(*ChannelWidthUpDX, *ChannelWidth);
        cover(*ChannelWidthUpDX, *LDD, 0);
        //        double v = 0.9*_dx;
        //        ChannelWidthUpDX->calcValue(v, MIN);
        FOR_ROW_COL_MV
        {
            ChannelAdj->Drc = std::max(0.01*_dx, _dx - ChannelWidthUpDX->Drc);
            ChannelWidthUpDX->Drc = _dx - ChannelAdj->Drc;
        }

        FOR_ROW_COL_MV_CH
        {
            ChannelDX->Drc = _dx/cos(asin(ChannelGrad->Drc));
            if (SwitchEfficiencyDET == 1)
                ChannelY->Drc = std::min(1.0, 1.0/(0.89+0.56*CohesionSoil->Drc));
            else
                if (SwitchEfficiencyDET == 2)
                    ChannelY->Drc = std::min(1.0, 0.79*exp(-0.85*CohesionSoil->Drc));
                else
                    if (SwitchEfficiencyDET == 3)
                        ChannelY->Drc = std::min(1.0, 1.0/(2.0*CohesionSoil->Drc));
        }



        if(SwitchChannelFlood)
        {
            ChannelDepth = ReadMap(LDDChannel, getvaluename("chandepth"));
            cover(*ChannelDepth, *LDD,0);
            /* obsolete
            ChannelLevee = NewMap(0);

              if (SwitchLevees)
                ChannelLevee = ReadMap(LDD, getvaluename("chanlevee"));
            if (!SwitchLevees)
                fill(*ChannelLevee, 0.0);
            */
        }

    }



}


double TWorld::LogNormalDist(double d50,double s, double d)
{

    double dev = log(1.0 + s/d50);
    double dev2 = (log(d)  - log(d50));
    return (1.0/(d *sqrt(2.0*3.14159) * log(1.0 + s/d50)))*exp(-dev2*dev2)/(4*dev*dev);

}
//---------------------------------------------------------------------------
// NOT USED FOR NOW
void TWorld::InitMulticlass(void)
{

    unity = NewMap(1.0);

    Qs = NewMap(0);
    Qsn = NewMap(0);

    //### erosion maps
    DetSplashTot = 0;
    DetFlowTot = 0;
    DepTot = 0;
    DetTot = 0;
    DepTot = 0;
    SoilLossTot = 0;
    SoilLossTotT= 0;
    SoilLossTotOutlet = 0;
    SoilLossTotSub = 0;
    SedTot = 0;

    //Qsoutflow = NewMap(0);
    DETFlow = NewMap(0);
    DETSplash = NewMap(0);
    DETSplashCum = NewMap(0);
    DETFlowCum = NewMap(0);
    DEP = NewMap(0);
    Sed = NewMap(0);
    TC = NewMap(0);
    Conc = NewMap(0);
    CG = NewMap(0);
    DG = NewMap(0);
    SettlingVelocity = NewMap(0);
    CohesionSoil = NewMap(0);
    Y = NewMap(0);

    GWDepth = NewMap(0);


    TotalDetMap = NewMap(0);
    TotalDepMap = NewMap(0);
    TotalChanDetMap = NewMap(0);
    TotalChanDepMap = NewMap(0);
    TotalSoillossMap = NewMap(0);
    TotalSed = NewMap(0);
    TotalConc = NewMap(0);


    if(SwitchErosion)
    {
        FOR_ROW_COL_MV
        {

            CohesionSoil->Drc = Cohesion->Drc + Cover->Drc*RootCohesion->Drc;
            // soil cohesion everywhere, plantcohesion only where plants

            if (SwitchEfficiencyDET == 1)
                Y->Drc = std::min(1.0, 1.0/(0.89+0.56*CohesionSoil->Drc));
            else
                if (SwitchEfficiencyDET == 2)
                    Y->Drc = std::min(1.0, 0.79*exp(-0.85*CohesionSoil->Drc));
                else
                    if (SwitchEfficiencyDET == 3)
                        Y->Drc = std::min(1.0, 1.0/(2.0*CohesionSoil->Drc));

            //            if (StoneFraction->Drc > 0)
            //                Y->Drc = 0.84*exp(-6*StoneFraction->Drc);
            // GOED IDEE
        }

    }

    if(SwitchIncludeChannel)
    {
        ChannelY = NewMap(0);

    }


    if(SwitchErosion)
    {
        graindiameters.clear();
        settlingvelocities.clear();

        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV
            {
                CG->Drc = pow((D50->Drc+5)/0.32, -0.6);
                DG->Drc = pow((D50->Drc+5)/300, 0.25);
                SettlingVelocity->Drc = GetSV(D50->Drc);
            }
        }else if(SwitchUseGrainSizeDistribution)
        {

            if(SwitchEstimateGrainSizeDistribution)
            {
                if(numgrainclasses == 0)
                {
                    ErrorString = "Could not simulate 0 grain classes" +QString("\n")
                            + "Please provide a positive number";
                            throw 1;

                }


                distD50 = 0;
                distD90 = 0;
                int count = 0;
                FOR_ROW_COL_MV
                {
                    distD50 += D50->Drc;
                    distD90 += D90->Drc;
                    count++;
                }
                distD50 = distD50/count;
                distD90 = distD90/count;

                double s = distD90- distD50;
                double s2l = std::max(distD50 - 2*s,distD50);
//VJ leads always to dist D50 of course
                double s2r = 2 * s;

                int classesleft = numgrainclasses;
                int mod2 = classesleft % 2;
                if(mod2 == 1)
                {
                    classesleft -= 1;
                }

                for(int i = 1; i < classesleft/2 + 1 ; i++)
                {
                    double d = (distD50 - s2l) + ((double)i) * s2l/(1.0 + double(classesleft/2.0) );
                    graindiameters.append(d);
                    W_D.append(NewMap(s2l/(1.0 + double(classesleft/2.0) )));
                }
                if(mod2 == 1)
                {
                    graindiameters.append(distD50);
                    W_D.append(NewMap(0.5 *s2l/(1.0 + double(classesleft/2.0) ) + 0.5 * s2r/(1.0 + double(classesleft/2.0))));
                }

                for(int i = 1; i < classesleft/2 + 1; i++)
                {
                    double d = (distD50) + ((double)i) *s2r/(1.0 + double(classesleft/2.0) );
                    graindiameters.append(d);
                    W_D.append(NewMap(s2r/(1.0 + double(classesleft/2.0))));
                }

                FOR_GRAIN_CLASSES
                {

                    settlingvelocities.append(GetSV(graindiameters.at(d)));

                    FOR_ROW_COL_MV
                    {
                        W_D.Drcd = W_D.Drcd*LogNormalDist(D50->Drc,D90->Drc -D50->Drc,graindiameters.at(d));
                    }

                    StorageDep_D.append(NewMap(0.0));
                    Storage_D.append(NewMap(0.0));
                    RStorageDep_D.append(NewMap(0.0));
                    RStorage_D.append(NewMap(0.0));
                }

                FOR_ROW_COL_MV
                {
                    double wtotal = 0;
                    FOR_GRAIN_CLASSES
                    {
                        wtotal += (W_D).Drcd;
                    }

                    if(wtotal != 0)
                    {
                        FOR_GRAIN_CLASSES
                        {
                            (W_D).Drcd = (W_D).Drcd/wtotal;
                        }
                    }
                }

            }

            if(SwitchUseGrainSizeDistribution && SwitchReadGrainSizeDistribution)
            {

                numgrainclasses = 0;
                QStringList diamlist = getvaluename("Grain size class maps").split(",", QString::SkipEmptyParts);

                for(int i = 0; i < diamlist.count(); i++)
                {
                    double diam = diamlist.at(i).toDouble();
                    if( diam > 0.0)
                    {
                        numgrainclasses++;
                        graindiameters.append(diam);

                        settlingvelocities.append(GetSV(diam));

                        W_D.append(ReadMap(LDD,"GSD_"+diamlist.at(i)));

                        graindiameters.clear();

                        StorageDep_D.append(NewMap(0.0));
                        Storage_D.append(NewMap(0.0));
                        RStorageDep_D.append(NewMap(0.0));
                        RStorage_D.append(NewMap(0.0));
                    }
                }

                if(numgrainclasses == 0)
                {
                    ErrorString = "Could not interpret grain classes from the string: \n"
                                  +  getvaluename("Grain size class maps") + "\n"
                            + "Please provide positive values seperated by commas.";
                            throw 1;
                }
            }
        }


    }


}
//---------------------------------------------------------------------------
void TWorld::GetInputData(void)
{


    //## calibration factors
    ksatCalibration = getvaluedouble("Ksat calibration");
    nCalibration = getvaluedouble("N calibration");
    if (nCalibration == 0)
    {
        ErrorString = QString("Calibration: the calibration factor for Mannings n for slopes cannot be zero.");
        throw 1;
    }

    thetaCalibration = getvaluedouble("Theta calibration");
    psiCalibration = getvaluedouble("Psi calibration");
    ChnCalibration = getvaluedouble("Channel N calibration");
    if (ChnCalibration == 0)
    {
        ErrorString = QString("Calibration: the calibration factor for Mannings n for channels cannot be zero.");
        throw 1;
    }

    ChKsatCalibration = getvaluedouble("Channel Ksat calibration");
    SplashDelivery = getvaluedouble("Splash Delivery Ratio");
    if (SplashDelivery == 0)
    {
        ErrorString = QString("Calibration: the splash delivery factor cannot be zero.");
        throw 1;
    }
    DepositedCohesion = getvaluedouble("Particle Cohesion of Deposited Layer");

    st_scCalibration = getvaluedouble("Soil Cohesion Calibration");
    st_sifaCalibration = getvaluedouble("Soil Internal Friction Angle Calibration");
    st_sdCalibration = getvaluedouble("Soil Depth Calibration");
    st_csdCalibration = getvaluedouble("Create Stable Initial Safety Factor");
    st_csdsfCalibration = getvaluedouble("Minimum Safety Factor Calibration");


    Calibrate_EP = getvaluedouble("Erosive Power Calibration");
    Calibrate_TC = getvaluedouble("Transport Capacity Calibration");
    Calibrate_SV = getvaluedouble("Settling Velocity Calibration");
    Calibrate_YS = getvaluedouble("Internal Friction Angle");
    Calibrate_DV = getvaluedouble("Dynamic Viscosity Calibration");
    Calibrate_DF = getvaluedouble("Drag Force Calibration");
    Calibrate_SPF = getvaluedouble("Solid Phase Friction Calibration");
    Calibrate_DC = getvaluedouble("Deposition Criteria Calibration");
    Calibrate_ESC = getvaluedouble("Erosion Cohesion Calibration");
    Calibrate_EGC = getvaluedouble("Erosion Grain Size Calibration");


    StemflowFraction = getvaluedouble("Stemflow fraction");
    CanopyOpeness = getvaluedouble("Canopy Openess");
    //  maxFloodLevel = getvaluedouble("Max flood level");
    //  minFloodDt = getvaluedouble("Min flood dt");

    //VJ 110829 water repellency
    waterRep_a = getvaluedouble("Water Repellency A");
    waterRep_b = getvaluedouble("Water Repellency B");
    waterRep_c = getvaluedouble("Water Repellency C");
    waterRep_d = getvaluedouble("Water Repellency D");

    //## catchment data
    LDD = InitMask(getvaluename("ldd"));
    // THIS SHOULD BE THE FIRST MAP
    // LDD is also mask and reference file, everthing has to fit LDD
    // channels use channel LDD as mask

    tm = NewMap(0); // temp map for aux calculations
    tma = NewMap(0); // temp map for aux calculations
    tmb = NewMap(0); // temp map for aux calculations
    tmc = NewMap(0); // temp map for aux calculations
    tmd = NewMap(0); // temp map for aux calculations
    tme = NewMap(0); // temp map for aux calculations
    tmf = NewMap(0); // temp map for aux calculations
    difkin =  NewMap(0); // temp map for aux calculations

    Qmax = NewMap(0.0);


    for (int i = 0; i < 32; i++)
        SubsMaps[i].m = NULL;  // initialize substance structures

    // flood maps
    DEM = ReadMap(LDD, getvaluename("dem"));
    Grad = ReadMap(LDD, getvaluename("grad"));  // must be SINE of the slope angle !!!
    checkMap(*Grad, LARGER, 1.0, "Gradient must be SINE of slope angle (not tangent)");
    calcValue(*Grad, 0.001, MAX);

    Outlet = ReadMap(LDD, getvaluename("outlet"));
    cover(*Outlet, *LDD, 0);

    bool needset =true;
    // fill outlet with zero, some users have MV where no outlet
    FOR_ROW_COL_MV
    {
        if(needset && Outlet->Drc > 0)
        {
            needset = false;
            c_outlet = c;
            r_outlet = r;
            r_plot = r_outlet;
            c_plot = c_outlet;
        }
        /*if (Outlet->Drc == 1)
        {
            if (LDD->Drc != 5)
            {
                ErrorString = "Main outlet gridcell does not coincide with pit in LDD";
                //throw 1;
            }
            else
            {
                c_outlet = c;
                r_outlet = r;
                r_plot = r_outlet;
                c_plot = c_outlet;
            }
        }*/
    }

    PointMap = ReadMap(LDD,getvaluename("outpoint"));
    //map with points for output data

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

    //## landuse and surface data
    N = ReadMap(LDD,getvaluename("manning"));
    calcValue(*N, nCalibration, MUL); //VJ 110112 moved

    RR = ReadMap(LDD,getvaluename("RR"));
    PlantHeight = ReadMap(LDD,getvaluename("CH"));
    LAI = ReadMap(LDD,getvaluename("lai"));
    Cover = ReadMap(LDD,getvaluename("cover"));
    if (SwitchLitter)
    {
        Litter = ReadMap(LDD,getvaluename("litter"));

        checkMap(*Litter, LARGER, 1.0, "vegetation litter/herb cover fraction cannot be more than 1");
        checkMap(*Litter, SMALLER, 0.0, "Litter cover fraction must be >= 0");
        checkMap(*Litter, LARGER, 1.0, "Litter cover fraction must be <= 1.0");
    }
    else
    {
        Litter = NewMap(0);
    }

    LitterSmax = getvaluedouble("Litter interception storage");


    checkMap(*RR, SMALLER, 0.0, "Raindom roughness RR must be >= 0");
    checkMap(*N, SMALLER, 1e-6, "Manning's N must be > 0.000001");
    checkMap(*LAI, SMALLER, 0.0, "LAI must be >= 0");
    checkMap(*Cover, SMALLER, 0.0, "Cover fraction must be >= 0");
    checkMap(*Cover, LARGER, 1.0, "Cover fraction must be <= 1.0");
    checkMap(*PlantHeight, SMALLER, 0.0, "Cover fraction must be >= 0");

    LandUnit = ReadMap(LDD,getvaluename("landunit"));  //VJ 110107 added
    GrassFraction = NewMap(0);

    if (SwitchGrassStrip)
    {
        KsatGrass = ReadMap(LDD,getvaluename("ksatgras"));
        GrassWidthDX = ReadMap(LDD,getvaluename("grasswidth"));
        copy(*GrassFraction, *GrassWidthDX);
        calcValue(*GrassFraction, _dx, DIV);
        StripN = getvaluedouble("Grassstrip Mannings n");
        FOR_ROW_COL_MV
        {
            if (GrassWidthDX != 0)
            {
                N->Drc = N->Drc*(1-GrassFraction->Drc)+StripN*GrassFraction->Drc;
                Cover->Drc = Cover->Drc*(1-GrassFraction->Drc) + 0.95*GrassFraction->Drc;
                LAI->Drc = LAI->Drc*(1-GrassFraction->Drc) + 5.0*LAI->Drc;
            }
            //adjust mann N Cover and height
        }
    }
    else
        KsatGrass = NewMap(0);

    StoneFraction  = ReadMap(LDD,getvaluename("stonefrc"));
    // WheelWidth  = ReadMap(LDD,getvaluename("wheelwidth"));

    if (SwitchRoadsystem)
    {
        RoadWidthDX  = ReadMap(LDD,getvaluename("road"));
        checkMap(*RoadWidthDX, LARGER, _dx, "road width cannot be larger than gridcell size");
        FOR_ROW_COL_MV
        {
            N->Drc = N->Drc * (1-RoadWidthDX->Drc/_dx) + 0.001*RoadWidthDX->Drc/_dx;
        }
    }
    else
        RoadWidthDX = NewMap(0);

    if (SwitchHardsurface)
    {
        HardSurface = ReadMap(LDD,getvaluename("hardsurf"));
        calcValue(*HardSurface, 1.0, MIN);
        calcValue(*HardSurface, 0.0, MAX);
        FOR_ROW_COL_MV
        {
            N->Drc = N->Drc * (1-HardSurface->Drc) + 0.001*HardSurface->Drc;
        }
    }
    else
        HardSurface = NewMap(0);

    firstssreport = true;
    sfset = false;

    //## infiltration data
    if(InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE)
    {
        Ksat1 = ReadMap(LDD,getvaluename("ksat1"));
        SoilDepth1 = ReadMap(LDD,getvaluename("soildep1"));
        calcValue(*SoilDepth1, 1000, DIV);
        //VJ 101213 fixed bug: convert from mm to m
        // can be zero for outcrops
        FOR_ROW_COL_MV
        {
            SoilDepth1->Drc *= st_sdCalibration;
            if (SoilDepth1->Drc < 0)
            {
                ErrorString = QString("SoilDepth1 values < 0 at row %1, col %2").arg(r).arg(c);
                throw 1;
            }
        }

        ThetaS1 = ReadMap(LDD,getvaluename("thetas1"));
        ThetaI1 = ReadMap(LDD,getvaluename("thetai1"));

        calcValue(*ThetaI1, thetaCalibration, MUL); //VJ 110712 calibration of theta
        calcMap(*ThetaI1, *ThetaS1, MIN); //VJ 110712 cannot be more than porosity

        //VJ 101221 all infil maps are needed except psi
        if(InfilMethod != INFIL_KSAT)
        {
            Psi1 = ReadMap(LDD,getvaluename("psi1"));
            calcValue(*Psi1, psiCalibration, MUL); //VJ 110712 calibration of psi
            calcValue(*Psi1, 0.01, MUL); // convert to meter
        }

        if (SwitchTwoLayer)
        {
            ThetaS2 = ReadMap(LDD,getvaluename("thetaS2"));
            ThetaI2 = ReadMap(LDD,getvaluename("thetaI2"));

            calcValue(*ThetaI2, thetaCalibration, MUL); //VJ 110712 calibration of theta
            calcMap(*ThetaI2, *ThetaS2, MIN); //VJ 110712 cannot be more than porosity

            //VJ 101221 all infil maps are needed except psi
            if(InfilMethod != INFIL_KSAT)
            {
                Psi2 = ReadMap(LDD,getvaluename("psi2"));
                calcValue(*Psi2, psiCalibration, MUL); //VJ 110712 calibration of psi
                calcValue(*Psi2, 0.01, MUL);
            }

            Ksat2 = ReadMap(LDD,getvaluename("ksat2"));
            SoilDepth2 = ReadMap(LDD,getvaluename("soilDep2"));
            calcValue(*SoilDepth2, 1000, DIV);
            //VJ 101213 fixed bug: convert from mm to m

            FOR_ROW_COL_MV
            {
                SoilDepth2->Drc *= st_sdCalibration;
                if (SoilDepth2->Drc < 0)
                {
                    ErrorString = QString("SoilDepth2 values < 0 at row %1, col %2").arg(r).arg(c);
                    throw 1;
                }
            }
        }
        if (SwitchInfilCrust)
        {
            CrustFraction = ReadMap(LDD,getvaluename("crustfrc"));
            checkMap(*CrustFraction, LARGER, 1.0, "crust fraction cannot be more than 1");
            KsatCrust = ReadMap(LDD,getvaluename("ksatcrst"));
        }
        else
        {
            CrustFraction = NewMap(0);
            KsatCrust = NewMap(0);
        }

        if (SwitchInfilCompact)
        {
            CompactFraction = ReadMap(LDD,getvaluename("compfrc"));
            checkMap(*CompactFraction, LARGER, 1.0, "compacted area fraction cannot be more than 1");
            KsatCompact = ReadMap(LDD,getvaluename("ksatcomp"));
        }
        else
        {
            CompactFraction = NewMap(0);
            KsatCompact = NewMap(0);
        }
        FOR_ROW_COL_MV
        {
            if (CrustFraction->Drc +  CompactFraction->Drc > 1.0)
            {
                CrustFraction->Drc = 1.0-CompactFraction->Drc;
            }
        }
    }

    // SWATRE infiltration read maps and structures
    if (InfilMethod == INFIL_SWATRE)
    {
        // read all Swatre profiles
        ProfileID = ReadMap(LDD,getvaluename("profmap"));

        if (SwitchGrassStrip)
            ProfileIDGrass = ReadMap(LDD,getvaluename("profgrass"));

        if (SwitchInfilCrust)
        {
            CrustFraction = ReadMap(LDD,getvaluename("crustfrc"));
            ProfileIDCrust = ReadMap(LDD,getvaluename("profcrst"));
        }
        else
            CrustFraction = NewMap(0);

        RepellencyFraction = NewMap(1.0);
        if (SwitchWaterRepellency)
        {
            RepellencyCell = ReadMap(LDD,getvaluename("repelcell"));
            // values of 1 calculate repellency
        }
        else
            RepellencyCell = NewMap(0); //no repellency anywhere


        // repellency to 1, no effect

        if (SwitchInfilCompact)
        {
            CompactFraction = ReadMap(LDD,getvaluename("compfrc"));
            ProfileIDCompact = ReadMap(LDD,getvaluename("profcomp"));
        }
        else
            CompactFraction = NewMap(0);

        // read the swatre tables and make the information structure ZONE etc
        ReadSwatreInputNew();

        // obsolete
        //      int res = ReadSwatreInput(SwatreTableName, SwatreTableDir);
        //      if (res)
        //         throw res;
    }

    if(SwitchErosion)
    {
        Cohesion = ReadMap(LDD,getvaluename("coh"));
        RootCohesion = ReadMap(LDD,getvaluename("cohadd"));
        AggrStab = ReadMap(LDD,getvaluename("AggrStab"));

        D50 = ReadMap(LDD,getvaluename("D50"));
        D90 = ReadMap(LDD,getvaluename("D90"));

        FOR_ROW_COL_MV
        {
            D50->Drc = D50->Drc *Calibrate_EGC;
            D50->Drc = D90->Drc *Calibrate_EGC;
            Cohesion->Drc = Cohesion->Drc *Calibrate_ESC;
        }
    }


    EntrainmentDet = NewMap(0.0);
    EntrainmentDep = NewMap(0.0);
    TotalEntrainmentDet = NewMap(0.0);
    TotalEntrainmentDep = NewMap(0.0);

    EntrainmentTC = NewMap(0.0);
    Entrainmentshearstressc = NewMap(0.0);
    Entrainmentshearstress = NewMap(0.0);
    DepositionT = NewMap(0.0);

    ChannelEntrainmentDet = NewMap(0.0);
    ChannelEntrainmentDep = NewMap(0.0);
    ChannelTotalEntrainmentDet = NewMap(0.0);
    ChannelTotalEntrainmentDep = NewMap(0.0);

    if(SwitchSolidPhase && SwitchEntrainment)
    {

        SoilRockMaterial = ReadMap(LDD,getvaluename("debrismaterial"));
        SoilRockWater = NewMap(0.0);
        SoilRockSize = ReadMap(LDD,getvaluename("rocksize"));
        SoilRockDensity = ReadMap(LDD,getvaluename("rockdensity"));
        SoilRockIFA = ReadMap(LDD,getvaluename("rockifa"));

        FOR_ROW_COL_MV
        {
            SoilRockIFA->Drc = SoilRockIFA->Drc * Calibrate_YS;
            SoilRockMaterial->Drc = SoilRockMaterial->Drc* _dx*_dx;

            if(InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE)
            {
                SoilRockWater->Drc = SoilRockMaterial->Drc * ThetaI1->Drc;
            }else
            {
                SoilRockWater->Drc = 0;
            }

        }
        if(SwitchIncludeChannel)
        {
            RSoilRockMaterial = NewMap(0.0);//ReadMap(LDD,getvaluename("channeldebrismaterial"));
            RSoilRockWater =NewMap(0.0);
            RSoilRockSize = NewMap(0.0);//ReadMap(LDD,getvaluename("channelrocksize"));
            RSoilRockDensity = NewMap(2500.0);//ReadMap(LDD,getvaluename("channelrockdensity"));
            RSoilRockIFA = NewMap(0.0);//ReadMap(LDD,getvaluename("channelrockifa"));

        }

    }else
    {
        SoilRockMaterial = NewMap(0.0);
        SoilRockWater = NewMap(0.0);
        SoilRockSize = NewMap(0.0);
        SoilRockDensity = NewMap(2500.0);
        SoilRockIFA = NewMap(0.0);

        if(SwitchIncludeChannel)
        {
            RSoilRockMaterial = NewMap(0.0);
            RSoilRockWater =NewMap(0.0);
            RSoilRockSize = NewMap(0.0);
            RSoilRockDensity = NewMap(2500.0);
            RSoilRockIFA = NewMap(0.0);

        }

    }
    LDDChange = NewMap(0.0);
    UF2D_DEMOriginal = NewMap(0.0);
    DEMChange = NewMap(0.0);
    DFTotalInitiationHeight = NewMap(0.0);
    DFSoilDensity = NewMap(2000.0);
    if(SwitchSlopeStability)
    {
        DFSafetyFactor = NewMap(0.0);
        DEMOriginal = NewMap(0.0);
        DFUnstable = NewMap(0.0);
        DFSFCalibration = NewMap(0.0);
        DFInitiationHeight = NewMap(0.0);
        DFSFIterations = NewMap(0.0);
        DEMIterate = NewMap(0.0);
        FailLimit = NewMap(0.0);
        DFSoilDepth = NewMap(0.0);
        DFSurfaceWaterHeight = NewMap(0.0);
        DFSoilCohesion = NewMap(0.0);
        DFWaterHeight = NewMap(0.0);
        DFWaterSuction = NewMap(0.0);
        DFPlantCohesion = NewMap(0.0);
        DFPlantPressure = NewMap(0.0);
        DFThreshold = NewMap(0.0);
        DFThreshold1 = NewMap(0.0);
        DFSlope = NewMap(0.0);
        DFSlopeX = NewMap(0.0);
        DFSlopeY = NewMap(0.0);

        DFSlopeXf = NewMap(0.0);
        DFSlopeYf = NewMap(0.0);
        DFFailureMask = NewMap(1.0);

        DFForcingDemand = NewMap(0.0);
        DFForcingCapacity = NewMap(0.0);
        DFForcing = NewMap(0.0);
        DFForcingAdded = NewMap(0.0);
        DFForcingUp = NewMap(0.0);
        DFForcingUpAdded = NewMap(0.0);
        DFForcing2 = NewMap(0.0);
        DFForcingUp2 = NewMap(0.0);

        if(SwitchUpslopeForcing)
        {
            MaximumUpslopeForcing = NewMap(0.0);
        }

        if(SwitchDownslopeForcing)
        {
            MinimumDownslopeForcing = NewMap(1e30);
        }

        SF_Calibrate_Initial = (getvalueint("Create Stable Initial Safety Factor") == 1);
        SF_Calibrate_Margin = getvaluedouble("Minimum Safety Factor Calibration");
        SF_Calibrate_First = true;
        SF_BedRock_First = true;
        if(SF_Calibrate_LF)
        {
            DFFailureMask = ReadMap(LDD,getvaluename("FailureMask"));
        }


        //calibration values for slope failure aspect
        st_scCalibration = getvaluedouble("Soil Cohesion Calibration");
        st_sifaCalibration = getvaluedouble("Soil Internal Friction Angle Calibration");
        st_sdCalibration = getvaluedouble("Soil Depth Calibration");
        st_csdCalibration = getvaluedouble("Create Stable Initial Safety Factor");
        st_csdsfCalibration = getvaluedouble("Minimum Safety Factor Calibration");

        //SF_Calibrate_Mask = true;

        DFSoilInternalFrictionAngleOrg = ReadMap(LDD,getvaluename("soilifa"));
        DFSoilInternalFrictionAngle = NewMap(0.0);
        DFSoilDensity = ReadMap(LDD,getvaluename("soildensity"));
        DFSoilRockSize = ReadMap(LDD,getvaluename("soilrocksize"));
        DFZERO = NewMap(0.0);

        if(SwitchBedrock)
        {

            DFUnstable2 = NewMap(0.0);
            DFAddedPressure = NewMap(0.0);
            DFInitiationHeight2 = NewMap(0.0);
            DFSFIterations2 = NewMap(0.0);
            DFSoilDepth2 = NewMap(0.0);
            DFSafetyFactor2 = NewMap(0.0);
            DFSoilInternalFrictionAngle2Org = ReadMap(LDD,getvaluename("soilifa2"));
            DFSoilInternalFrictionAngle2 = NewMap(0.0);
            DFSoilCohesion2Org = ReadMap(LDD,getvaluename("soilcohesion2"));
            DFSoilCohesion2 = NewMap(0.0);
            DFSoilDensity2 = ReadMap(LDD,getvaluename("soildensity2"));
            DFSoilRockSize2 = ReadMap(LDD,getvaluename("soilrocksize2"));
            DFSFCalibration2 = NewMap(0.0);

            FOR_ROW_COL_MV
            {
                DFSoilInternalFrictionAngle2Org->Drc = DFSoilInternalFrictionAngle2Org->Drc * st_sifaCalibration;
                DFSoilCohesion2Org->Drc = DFSoilCohesion2Org->Drc * st_scCalibration;
            }
        }

        FOR_ROW_COL_MV
        {
            DEMOriginal->Drc = DEM->Drc;
            DFSoilInternalFrictionAngleOrg->Drc *= st_sifaCalibration;
        }

        if(SwitchSeismic)
        {
            PGA = ReadMap(LDD,getvaluename("pga"));
            PGATiming = ReadMap(LDD,getvaluename("pgatiming"));
            StrengthLoss = NewMap(0.0);
            StrengthLoss2 = NewMap(0.0);
            PGAInitiated = NewMap(0.0);
            PGACurrent = NewMap(0.0);
            if(SwitchBedrock)
            {
                StrengthLoss2 = NewMap(0.0);
            }
        }

    }

    if(SwitchSlopeFailure || SwitchSlopeStability )
    {
        SFMIN = getvaluedouble("Minimum Safety Factor");
        SFMAX = getvaluedouble("Maximum Safety Factor");
        SFDISPLAYMAX = getvaluedouble("Maximum safety factor for display");
        SFFAILMIN = getvaluedouble("Minimum Failure Height");

    }

    InitMulticlass();

    //## read and initialize all channel maps and variables
    InitChannel();

    //## make shaded relief map for display.
    InitShade();



}
//---------------------------------------------------------------------------
/// called after get input data, initializes non-input maps and variables
void TWorld::IntializeData(void)
{
    //TO DO add units and descriptions --> TMmapVariables.h


    TSList_point.clear();
    TSList_rainav.clear();
    TSList_snowav.clear();
    TSList_q.clear();
    TSList_h.clear();
    TSList_qs.clear();
    TSList_c.clear();

    //totals for mass balance
    MB = 0;
    MBs = 0;
    nrCells = 0;
    FOR_ROW_COL_MV
    {
        nrCells+=1;
    }

    DX = NewMap(0);
    CellArea = NewMap(0);
    FOR_ROW_COL_MV
    {
        DX->Drc = _dx/cos(asin(Grad->Drc));
//        if (SwitchIncludeChannel)
//            if (ChannelDX->Drc > 0)
//            {
//                DX->Drc = ChannelDX->Drc;
//            }
        CellArea->Drc = DX->Drc * _dx;
    }
    CatchmentArea = mapTotal(*CellArea);

    SoilWidthDX = NewMap(0);
    // not implemented
    // WheelWidth = NewMap(0);
    // WheelWidthDX = NewMap(0);
    // GullyWidthDX = NewMap(0);

    // surface storage
    MDS = NewMap(0);
    FOR_ROW_COL_MV
    {
        double RRmm = 10 * RR->Drc;
        MDS->Drc = std::max(0.0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*tan(asin(Grad->Drc))*100);
        MDS->Drc /= 1000; // convert to m
    }

    //combination display
    COMBO_QOFCH = NewMap(0);
    COMBO_SS = NewMap(0);

    //### rainfall and interception maps
    BaseFlow = 0;
    RainTot = 0;
    RainTotmm = 0;
    Rainpeak = 0;
    RainpeakTime = 0;
    RainstartTime = -1;
    rainStarted = false;
    RainAvgmm = 0;
    SnowAvgmm = 0;
    SnowTot = 0;
    SnowTotmm = 0;
    Snowpeak = 0;
    SnowpeakTime = 0;
    Rain = NewMap(0);
    Rainc = NewMap(0);
    RainCum = NewMap(0);
    RainCumFlat = NewMap(0);
    RainNet = NewMap(0);
    LeafDrain = NewMap(0);
    //not used RainIntensity = NewMap(0);
    //not used RainM3 = NewMap(0);
    CStor = NewMap(0);
    Interc = NewMap(0);
    // litter
    LCStor = NewMap(0);
    LInterc = NewMap(0);
    LRainCum = NewMap(0);
    //houses
    HStor = NewMap(0);
    IntercHouse = NewMap(0);
    DStor = NewMap(0);

    Snowmelt = NewMap(0);
    Snowmeltc = NewMap(0);
    SnowmeltCum = NewMap(0);

    InterceptionLAIType = getvalueint("Canopy storage equation");
    if (InterceptionLAIType == 8)
    {
        SwitchInterceptionLAI = false;
        CanopyStorage = ReadMap(LDD,getvaluename("smax"));
    }
    if (SwitchInterceptionLAI)
    {
        CanopyStorage = NewMap(0); //in m !!!
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
    calcValue(*CanopyStorage, 0.001, MUL); // from mm to m
    //NOTE: LAI is still needed for canopy openness, can be circumvented with cover
    if (SwitchHouses)
    {
        //houses info:
        //housecover.map;Fraction of hard roof surface per cell (-);housecover");
        //roofstore.map;Size of interception storage of rainwater on roofs (mm);roofstore");
        //drumstore.map;Size of storage of rainwater drums (m3);drumstore");
        HouseCover = ReadMap(LDD,getvaluename("housecover"));

        RoofStore = ReadMap(LDD,getvaluename("roofstore"));
        calcValue(*RoofStore, 0.001, MUL);
        // from mm to m
        DrumStore = ReadMap(LDD,getvaluename("drumstore"));
    }
    else
        HouseCover = NewMap(0);

    HouseWidthDX = NewMap(0);
    FOR_ROW_COL_MV
    {
        //HouseCover->Drc = std::min(0.9, HouseCover->Drc);
        HouseWidthDX->Drc = std::min(_dx,  HouseCover->Drc *_dx);
        // assume there is always space next to house
        N->Drc = N->Drc * (1-HouseCover->Drc) + 0.25*HouseCover->Drc;
    }

    //### infiltration maps
    InfilTot = 0;
    InfilTotmm = 0;
    InfilKWTot = 0;
    IntercTot = 0;
    IntercTotmm = 0;
    //houses
    IntercHouseTot = 0;
    IntercHouseTotmm = 0;
    WaterVolTot = 0;
    WaterVolSoilTot = 0;
    WaterVolTotmm = 0;
    WaterVolRunoffmm = 0;

    floodTotmm= 0;
    floodVolTot = 0;
    floodVolTotInit = 0;
    floodVolTotMax = 0;
    floodAreaMax = 0;
    floodBoundaryTot = 0;

    InfilVolFlood = NewMap(0);
    InfilVolKinWave = NewMap(0);
    InfilVol = NewMap(0);
    InfilmmCum = NewMap(0);
    InfilVolCum = NewMap(0);
    fact = NewMap(0);
    fpot = NewMap(0);
    factgr = NewMap(0);
    fpotgr = NewMap(0);
    Ksateff = NewMap(0);
    FSurplus = NewMap(0);
    hesinfil = NewMap(0);
    FFull = NewMap(0);
    runoffFractionCell = NewMap(0);
    runoffTotalCell = NewMap(0);

    L1 = NewMap(0);
    L2 = NewMap(0);

    if (InfilMethod != INFIL_SWATRE && InfilMethod != INFIL_NONE)
    {
        Fcum = NewMap(0);

        Fcumgr = NewMap(1e-10);
        L1gr = NewMap(1e-10);
        L2gr = NewMap(1e-10);
        if (InfilMethod != INFIL_KSAT)
        {
            Soilwater = NewMap(0);
            Soilwater2 = NewMap(0);
            calc2Maps(*Soilwater, *ThetaI1, *SoilDepth1, MUL);
            if (SwitchTwoLayer)
            {
                calc2Maps(*Soilwater2, *ThetaI2, *SoilDepth2, MUL);
            }
        }
    }

    //### runoff maps
    Qtot = 0;
    QtotT = 0;
    QtotOutlet = 0;
    Qtotmm = 0;
    Qpeak = 0;
    QpeakTime = 0;
    WH = NewMap(0);
    WHbef = NewMap(0);
    WHrunoff = NewMap(0);
    WHmax = NewMap(0);
    WHstore = NewMap(0);
    WHroad = NewMap(0);
    WHGrass = NewMap(0);
    FlowWidth = NewMap(0);
    fpa = NewMap(0);
    V = NewMap(0);
    Vx = NewMap(0);
    Vy = NewMap(0);
    Alpha = NewMap(0);

    //    AlphaF = NewMap(0);
    //    QF = NewMap(0);
    //    QnF = NewMap(0);

    Q = NewMap(0);
    Qn = NewMap(0);

    QinKW = NewMap(0);
    QoutKW = NewMap(0);
    Qoutput = NewMap(0);
    Houtput = NewMap(0);
    Qsoutput = NewMap(0);
    //Qoutflow = NewMap(0); // value of Qn*dt in pits only
    q = NewMap(0);
    R = NewMap(0);
    Perim = NewMap(0);
    WaterVolin = NewMap(0);
    WaterVolall = NewMap(0);

    SwatreSoilModel = NULL;
    SwatreSoilModelCrust = NULL;
    SwatreSoilModelCompact = NULL;
    SwatreSoilModelGrass = NULL;
    // swatre get input data is called before, ReadSwatreInput
    if (InfilMethod == INFIL_SWATRE)
    {
        thetaTop = NewMap(0);

        precision = 5.0;
        // note "5" is a precision factor dewtermining next timestep, set to 5 in old lisem

        // VJ 110420 added tiledrain depth for all profiles, is all used in infiltration
        SwatreSoilModel = InitSwatre(ProfileID);//, initheadName, TileDepth, swatreDT);
        if (SwatreSoilModel == NULL)
            throw 3;

        if (SwitchInfilCrust)// || SwitchWaterRepellency)
        {
            SwatreSoilModelCrust = InitSwatre(ProfileIDCrust);//, initheadName, TileDepth, swatreDT);
            if (SwatreSoilModelCrust == NULL)
                throw 3;
        }
        if (SwitchInfilCompact)
        {
            SwatreSoilModelCompact = InitSwatre(ProfileIDCompact);//, initheadName, TileDepth, swatreDT);
            if (SwatreSoilModelCompact == NULL)
                throw 3;
        }
        if (SwitchGrassStrip)
        {
            SwatreSoilModelGrass = InitSwatre(ProfileIDGrass);//, initheadName, TileDepth, swatreDT);
            if (SwatreSoilModelGrass == NULL)
                throw 3;
        }
        initSwatreStructure = true;
        // flag: structure is created and can be destroyed in function destroydata
    }

    if(SwitchErosion && SwitchUseMaterialDepth)
    {
        Storage = ReadMap(LDD, getvaluename("detmat"));

        FOR_ROW_COL_MV
        {
            if(Storage->Drc != -1)
            {
                Storage->Drc = Storage->Drc * ChannelAdj->Drc * DX->Drc;
            }else
            {
                Storage->Drc = -999999;
            }
        }
        StorageDep = NewMap(0.0);
        SedimentMixingDepth = ReadMap(LDD, getvaluename("sedmixdepth"));

        if(SwitchIncludeChannel)
        {
            if(SwitchErosion && SwitchUseMaterialDepth)
            {
                RStorage = ReadMap(LDD, getvaluename("chandetmat"));

                FOR_ROW_COL_MV
                {
                    if(RStorage->Drc != -1)
                    {
                        RStorage->Drc = RStorage->Drc * ChannelWidth->Drc * DX->Drc;
                    }else
                    {
                        Storage->Drc = -999999;
                    }
                }
                RStorageDep = NewMap(0.0);
                RSedimentMixingDepth = ReadMap(LDD, getvaluename("chansedmixdepth"));
            }
        }

    }

    if(SwitchErosion && SwitchUseGrainSizeDistribution)
    {
        IW_D.clear();
        RW_D.clear();

        FOR_GRAIN_CLASSES
        {
            if(SwitchIncludeChannel)
            {
                RW_D.append(NewMap(0.0));
            }
            IW_D.append(NewMap(0.0));

            FOR_ROW_COL_MV
            {
                IW_D.Drcd = W_D.Drcd;
                if(SwitchUseMaterialDepth)
                {
                    if(!(Storage->Drc > 0))
                    {
                            Storage_D.Drcd = W_D.Drcd * Storage->Drc;
                    }else
                    {
                        Storage_D.Drcd = -999999;
                    }
                }
                if(SwitchIncludeChannel)
                {
                    RW_D.Drcd = W_D.Drcd;
                    if(SwitchUseMaterialDepth)
                    {
                        if(!(RStorage->Drc > 0))
                        {
                                RStorage_D.Drcd = RW_D.Drcd * RStorage->Drc;
                        }else
                        {
                            RStorage_D.Drcd = -999999;
                        }
                    }
                }
            }
        }
    }

    MaximumDebrisFlowHeight = NewMap(0.0);
    MaximumDebrisFlowVelocity = NewMap(0.0);
    TotalSlopeFailure = NewMap(0.0);
    MinimumSafetyFactor = NewMap(0.0);

    //VJ 110113 all channel and buffer initialization moved to separate functions
    //calculate slope, outlets and pitches for kinematic 2D
    //K2Dslope also used for transport capacity of overland flow!

   // MakeWatersheds();
    FindBaseFlow();

    //Initialize all variables related to surface flow
    UF_Init();

    //read and initialize the flow barrier table and map
    InitFlowBarriers();

    //read and initialize the flow barrier table and map
    InitInflow();

    //create a function object referring to the cellprocesses wrapper
    fcompute = std::bind((&TWorld::CellProcesses),this,std::placeholders::_1);

    //creation callable function object by binding TWorld oject to member function
    //this leaves only the thread id as a parameter for the compute function
    flowcompute = std::bind((&TWorld::UF_Compute),this,std::placeholders::_1);
}
//---------------------------------------------------------------------------
void TWorld::IntializeOptions(void)
{
    nrRainfallseries = 0;
    nrSnowmeltseries = 0;

    //dirs and names
    resultDir.clear();
    inputDir.clear();
    outflowFileName = QString("totals.txt");//.clear();
    totalErosionFileName = QString("erosion.map");//.clear();
    totalDepositionFileName = QString("deposition.map");//.clear();
    totalChanErosionFileName = QString("chandet.map");//.clear();
    totalChanDepositionFileName = QString("chandep.map");//.clear();
    totalSoillossFileName = QString("soilloss.map");//.clear();
    totalLandunitFileName = QString("totlandunit.txt");//.clear();
    outflowFileName = QString("hydrohgraph.csv");//.clear();

    FinalFluidPhaseFileName = QString("finalfv.map");//.clear();
    FinalSolidPhaseFileName = QString("finalsv.map");//.clear();

    floodLevelFileName = QString("floodmaxH.map");//.clear();
    floodTimeFileName = QString("floodtime.map");//.clear();
    floodFEWFileName = QString("floodstart.map");//.clear();
    floodMaxVFileName = QString("floodmaxV.map");//.clear();

    floodStatsFileName = QString("floodstats.csv");//.clear();

    rainfallMapFileName = QString("rainfall.map");
    interceptionMapFileName = QString("interception.map");
    infiltrationMapFileName = QString("infiltration.map");
    runoffMapFileName = QString("runoff.map");
    runoffFractionMapFileName = QString("rofraction.map");
    channelDischargeMapFileName = QString("chandism3.map");

    FileName_MaxDebrisFlowHeight = QString("maximumdebrisflowheight.map");
    FileName_MaxDebrisFlowVelocity = QString("maximumdebrisflowvelocity.map");
    FileName_DebrisFlowStart = QString("debrisflowstart.map");
    FileName_Entrainment = QString("entrainment.map");
    FileName_SlopeFailure = QString("slopefailure.map");
    FileName_MinimumSafetyFactor = QString("minimumsafetyfactor.map");
    FileName_EntrainmentDet = QString("SolidEnt.map");
    FileName_EntrainmentDep = QString("SolidDep.map");
    FileName_MaximumUpSlopeForcing = QString("MaxForcingUpslope.map");
    FileName_MinimumDownSlopeForcing = QString("MinForcingDownslope.map");

    rainFileName.clear();
    SwitchLimitTC = false;
    rainFileDir.clear();
    snowmeltFileName.clear();
    snowmeltFileDir.clear();
    SwatreTableDir.clear();
    SwatreTableName = QString("profile.inp");//.clear();
    resultFileName.clear();

    SwitchUse2Layer = false;
    SwitchUseGrainSizeDistribution = false;
    SwitchReadGrainSizeDistribution = false;
    SwitchHardsurface = false;
    SwitchLimitDepTC = false;
    SwatreInitialized = false;
    SwitchInfilGA2 = false;
    SwitchWheelPresent = false;
    SwitchCompactPresent = false;
    SwitchIncludeChannel = false;
    SwitchChannelFlood = false;
    SwitchChannelBaseflow = false;
    startbaseflowincrease = false;
    SwitchChannelInfil = false;
    SwitchAllinChannel = false;
    SwitchErosion = false;
    SwitchAltErosion = false;
    SwitchSimpleDepression = false;
    SwitchBuffers = false;
    SwitchSedtrap = false;
    SwitchRainfall = true; //VL 110103 add rainfall default true
    SwitchSnowmelt = false;
    SwitchRunoffPerM = false;
    SwitchInfilCompact = false;
    SwitchInfilCrust = false;
    SwitchGrassStrip = false;
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
    SwitchEfficiencyDET = 1;
    SwitchStoninessDET = false;
    SwitchWritePCRnames = false;
    SwitchWriteCommaDelimited = true;
    SwitchWritePCRtimeplot = false;
    SwitchNoErosionOutlet = false;
    SwitchDrainage = false;
    SwitchPestout = false;
    SwitchSeparateOutput = false;
    SwitchInterceptionLAI = false;
    SwitchTwoLayer = false;
    SwitchSimpleSedKinWave = false;
    SwitchSOBEKoutput = false;
    SwitchPCRoutput = false;
    SwitchGeometric = true;
    SwitchPercolation = true;

    SwitchWriteHeaders = true; // write headers in output files in first timestep

    initSwatreStructure = false;
    // check to flag when swatre 3D structure is created, needed to clean up data
    SwitchFlowBarriers = true;
    SwitchInflow = false;
    SwitchBarriers = true;
    SwitchMaxVolume = true;
    SwitchChannelMaxVolume = true;
    SwitchUFInitial = true;
    SwitchUFForced = true;


}
//---------------------------------------------------------------------------
void TWorld::FindBaseFlow()
{

    if(SwitchChannelBaseflow)
    {
        BaseFlowDischarges = ReadMap(LDD, getvaluename("baseflow"));
        BaseFlowInflow = NewMap(0.0);
        BaseFlowInitialVolume = NewMap(0.0);
        FOR_ROW_COL_MV_CH
        {
            pcr::setMV(tma->Drc);
        }
        FOR_ROW_COL_MV_CH
        {
            pcr::setMV(tmb->Drc);
        }
        FOR_ROW_COL_MV_CH
        {
            pcr::setMV(tma->Drc);
            pcr::setMV(tmb->Drc);
            tmc->Drc = 0;
            tmd->Drc = 0;
        }

        for (int  ro = 0; ro < _nrRows; ro++){
        for (int  co = 0; co < _nrCols; co++){
        if(!pcr::isMV(LDDChannel->data[ro][co]))
        {
            if(LDDChannel->data[ro][co] == 5)
            {

                int ncells = 0;
                double infiltration = 0;
                double inflow = 0;
                double baseflow = BaseFlowDischarges->data[ro][co];

                int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
                int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

                /// Linked list of cells in order of LDD flow network, ordered from pit upwards
                LDD_LINKEDLIST *list = NULL, *temp = NULL;
                list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                list->prev = NULL;
                /// start gridcell: outflow point of area
                list->rowNr = ro;
                list->colNr = co;

                while (list != NULL)
                {
                    int i = 0;
                    bool  subCachDone = true; // are sub-catchment cells done ?
                    int rowNr = list->rowNr;
                    int colNr = list->colNr;

                    /** put all points that have to be calculated to calculate the current point in the list,
                     before the current point */
                    for (i=1; i<=9; i++)
                    {
                        int r, c;
                        int ldd = 0;

                        // this is the current cell
                        if (i==5)
                            continue;

                        r = rowNr+dy[i];
                        c = colNr+dx[i];

                        if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                            ldd = (int) LDDChannel->Drc;
                        else
                            continue;

                        // check if there are more cells upstream, if not subCatchDone remains true
                        if (pcr::isMV(tma->Drc) &&
                                FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                                INSIDE(r, c))
                        {
                            temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                            temp->prev = list;
                            list = temp;
                            list->rowNr = r;
                            list->colNr = c;
                            subCachDone = false;
                        }
                    }

                    // all cells above a cell are linked in a "sub-catchment or branch
                    // continue with water and sed calculations
                    // rowNr and colNr are the last upstreM cell linked
                    if (subCachDone)
                    {
                        int r = rowNr;
                        int c = colNr;
                        tma->Drc = 0;
                        ncells ++;

                        if(InfilMethod != INFIL_NONE && InfilMethod != INFIL_SWATRE)
                        {
                            double ksat = 0;//Ksat1->Drc;

                            if(SwitchChannelInfil)
                            {
                                ksat = ChannelKsat->Drc;
                            }

                            infiltration += (ChannelWidth->Drc) * DX->Drc * ksat *1.0/3600000.0;
                            tmd->Drc = (ChannelWidth->Drc) * DX->Drc * ksat *1.0/3600000.0;
                        }

                        temp=list;
                        list=list->prev;
                        free(temp);
                        // go to the previous cell in the list

                    }/* eof subcatchment done */
                } /* eowhile list != NULL */


                inflow = (baseflow + infiltration)/ ncells;

                qDebug() << baseflow << infiltration << ncells;

                list = NULL;
                temp = NULL;
                list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                list->prev = NULL;
                /// start gridcell: outflow point of area
                list->rowNr = ro;
                list->colNr = co;

                while (list != NULL)
                {
                    int i = 0;
                    bool  subCachDone = true; // are sub-catchment cells done ?
                    int rowNr = list->rowNr;
                    int colNr = list->colNr;

                    /** put all points that have to be calculated to calculate the current point in the list,
                     before the current point */
                    for (i=1; i<=9; i++)
                    {
                        int r, c;
                        int ldd = 0;

                        // this is the current cell
                        if (i==5)
                            continue;

                        r = rowNr+dy[i];
                        c = colNr+dx[i];

                        if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                            ldd = (int) LDDChannel->Drc;
                        else
                            continue;

                        // check if there are more cells upstream, if not subCatchDone remains true
                        if (pcr::isMV(tmb->Drc) &&
                                FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                                INSIDE(r, c))
                        {
                            temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                            temp->prev = list;
                            list = temp;
                            list->rowNr = r;
                            list->colNr = c;
                            subCachDone = false;
                        }
                    }

                    // all cells above a cell are linked in a "sub-catchment or branch
                    // continue with water and sed calculations
                    // rowNr and colNr are the last upstreM cell linked
                    if (subCachDone)
                    {
                        int r = list->rowNr;
                        int c = list->colNr;
                        tmb->Drc = 0;
                        BaseFlowInflow->Drc = inflow;

                        tmc->Drc += 1;

                        for (i=1;i<=9;i++)
                        {

                            int r, c, ldd = 0;

                            if (i==5)  // Skip current cell
                                continue;

                            r = rowNr+dy[i];
                            c = colNr+dx[i];

                            if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                                ldd = (int) LDDChannel->Drc;
                            else
                                continue;

                            if (INSIDE(r, c) &&
                                    FLOWS_TO(ldd, r,c,rowNr, colNr) &&
                                    !pcr::isMV(LDDChannel->Drc) )
                            {
                                tmc->data[list->rowNr][list->colNr] += tmc->Drc;
                                tmd->data[list->rowNr][list->colNr] += tmd->Drc;
                            }
                        }

                        r = list->rowNr;
                        c = list->colNr;

                        double q = (tmc->Drc * inflow - tmd->Drc);

                        double h, h1;
                        double w = ChannelWidth->Drc;
                        h = 1;
                        // first guess new h with old alpha
                        h1 = h;
                        double A;

                        if (q > 0)
                        {
                            double _23 = 2.0/3.0;
                            double F, dF;
                            int count = 0;

                            do{
                                h = h1;
                                if (h < 1e-10)
                                    break;
                                //double P = w+2*h;
                                //double A = h*w;

                                double P,R;
                                const double beta = 0.6;
                                const double _23 = 2.0/3.0;
                                double wh = h;
                                double FW = ChannelWidth->Drc;
                                double dw = /*0.5* */(ChannelWidthUpDX->Drc - FW); // extra width when non-rectamgular

                                if (dw > 0)
                                {
                                    //Perim = FW + 2*sqrt(wh*wh + dw*dw);
                                    P = FW + 2.0*wh/cos(atan(ChannelSide->Drc));
                                    A = FW*wh + wh*dw;
                                }
                                else
                                {
                                    P = FW + 2.0*wh/cos(atan(ChannelSide->Drc));
                                    A = FW*wh;
                                }

                                R = A/P;
                                F = std::max(0.0, 1.0 - q/(sqrt(ChannelGrad->Drc)/ChannelN->Drc*A*pow(R,_23)));
                                dF = (5.0*w+6.0*h)/(3.0*h*P);
                                h1 = h - F/dF;
                                // function divided by derivative
                                count++;
                            }while(fabs(h1-h) > 1e-10 && count < 20);
                        }

                        BaseFlowInitialVolume->data[list->rowNr][list->colNr] = A*DX->Drc;

                        temp=list;
                        list=list->prev;
                        free(temp);
                        // go to the previous cell in the list

                    }/* eof subcatchment done */
                } /* eowhile list != NULL */

            }
        }

        }}

        FOR_ROW_COL_MV_CH
        {
            tmc->Drc = 0;
        }
        FOR_ROW_COL_MV_CH
        {
            tmd->Drc = 0;
        }
    }



}
