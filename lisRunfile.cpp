/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

#include "lisemqt.h"
#include "global.h"

/*!
  \file lisRunfile.cpp
  \brief Read and parse the temporary runfile produce by the interface.

  Functions in here are doubled from the interface. The idea is to keep interface and model\n
  completely separate. In principle the model could be called directly with a runfile (not implemented). \n

functions: \n
- QString TWorld::getvaluename(QString vname) \n
- double TWorld::getvaluedouble(QString vname) \n
- int TWorld::getvalueint(QString vname) \n
- QString TWorld::CheckDir(QString p, bool makeit) \n
- QString TWorld::GetName(QString p) \n
- void TWorld::ParseRunfileData() \n
- void TWorld::GetRunFile() \n
*/

//---------------------------------------------------------------------------
QString TWorld::getvaluename(QString vname)
{
    for (int i = 0; i < nrrunnamelist; i++)
    {
        if(vname.toUpper() == runnamelist[i].name.toUpper())
        {
            // VJ 110420 special case
            if (InfilMethod == INFIL_SWATRE && runnamelist[i].name.toUpper() == QString("INITHEAD"))
            {
                QFileInfo info(inputDir + runnamelist[i].value + QString(".001"));
                if (!info.exists())
                {
                    ErrorString = "Filename not found for map \"<I>"+runnamelist[i].name + "\" - " + info.fileName();
                    throw 1;
                }
                else
                {
                    return inputDir + info.baseName();
                }
            }
            else
            {
                QFileInfo info(inputDir + runnamelist[i].value);
                if (!info.exists())
                {
                    ErrorString = "Filename not found for map \"<I>"+runnamelist[i].name + "\" - " + info.fileName();
                    throw 1;
                }
                else
                {
                    return inputDir + info.fileName();
                }
            }
        }
    }

    ErrorString = QString("Map ID: \"%1\" not found! You could be using an old runfile,\nor a map is not present.").arg(vname);
    throw 3;
}
//---------------------------------------------------------------------------
double TWorld::getvaluedouble(QString vname)
{
    for (int i = 0; i < nrrunnamelist; i++)
        if(vname.toUpper() == runnamelist[i].name.toUpper())
        {
//            return runnamelist[i].value.replace(",",".").toDouble();
            return runnamelist[i].value.toDouble();
        }

    ErrorString = QString("Variable ID: \"%1\" not found! You could be using an old runfile,\nor a variable has been added that is not present.").arg(vname);
    throw 3;
}
//---------------------------------------------------------------------------
int TWorld::getvalueint(QString vname)
{
    for (int i = 0; i < nrrunnamelist; i++)
        if(vname.toUpper() == runnamelist[i].name.toUpper())
        {
            return runnamelist[i].value.toInt();
        }

    ErrorString = QString("Variable ID: \"%1\" not found! You could be using an old runfile,\nor a variable has been added that is not present.").arg(vname);
    throw 3;
}
//---------------------------------------------------------------------------
QString TWorld::CheckDir(QString p, bool makeit)
{
    QString path;
    path = QDir(p).fromNativeSeparators(p);
    path = QDir(path).absoluteFilePath(path);

    if (!path.endsWith("/") || !path.endsWith("\\"))
        path = path + '/';

    if (!QDir(path).exists())
    {
        if (makeit)
        {
            QDir(path).mkpath(path);
            DEBUG("NOTE: Result dir created !");
            qDebug() << "NOTE: Result dir created !";
        }
        else
            path.clear();
    }

    return path;
}
//---------------------------------------------------------------------------
QString TWorld::GetName(QString p)
{
    QFileInfo fi(p);
    QStringList ss = fi.filePath().split("/");
    int n = ss.count();
    return(ss[n-1]);
}
//---------------------------------------------------------------------------
QString TWorld::checkOutputMapName(QString p, QString S, int i)
{
    if (p.isEmpty())
    {
        ErrorString = "Please give a name for the "+S+".";
        throw 1;
    }
    if (!p.contains("."))
    {
        if (i == 0)
            ErrorString = "Please give a name with an extention, such as \".map\".";
        else
            ErrorString = "Please give a name with an extention, such as \".txt\" or \".csv\".";
        throw 1;
    }
    return p;
}
//---------------------------------------------------------------------------
void TWorld::ParseRunfileData(void)
{
    int j=0;
   // int dummyAdvancedSed= 0;
    SwitchRainfall = true;
    SwitchSV = 1;

    // do all switches (checkbox options) first
    for (j = 0; j < nrrunnamelist; j++)
    {
        int iii = runnamelist[j].value.toInt();
        QString p1 = runnamelist[j].name;
        QString p = runnamelist[j].value;

    //    if (p1.compare("Nr user Cores") == 0) userCores = iii;
        //options in the main code, order is not important
        if (p1.compare("Include Erosion simulation")==0)     SwitchErosion =          iii == 1;
        if (p1.compare("Include main channels")==0)          SwitchIncludeChannel =   iii == 1;
        if (p1.compare("Include channel infil")==0)          SwitchChannelInfil     = iii == 1;
        if (p1.compare("Include channel baseflow")==0)       SwitchChannelBaseflow  = iii == 1;
        if (p1.compare("Include channel culverts")==0)       SwitchCulverts  = iii == 1;

        if (p1.compare("Variable Timestep")==0) SwitchVariableTimestep = iii == 1;
     //   if (p1.compare("Minimum Timestep Method")==0) SwitchMinDTfloodMethod = iii == 1;
        if (p1.compare("Use SWOF 2.0")==0) SwitchSWOFopen = iii == 1;
        if (p1.compare("Use MUSCL")==0) SwitchMUSCL = iii == 1;
        if (p1.compare("Use time avg V")==0) SwitchTimeavgV = iii == 1;
        if (p1.compare("Flow Boundary 2D")==0)     FlowBoundaryType = iii;
        if (p1.compare("Advanced Options")==0)     SwitchAdvancedOptions = iii == 1;

        if (p1.compare("Detachment efficiency")==0)           SwitchEfficiencyDET = iii;
        if (p1.compare("SettlingVelocity")==0)                SwitchSV = iii;
        if (p1.compare("Use material depth")==0)              SwitchUseMaterialDepth  = iii == 1;
        if (p1.compare("No detachment boundary")==0)          SwitchNoBoundarySed  = iii == 1;
   //     if (p1.compare("Advanced sediment")==0) SwitchAdvancedSed = iii;
        if (p1.compare("Use 2 phase flow")==0)                SwitchUse2Phase = iii;
        if (p1.compare("Include River diffusion")==0)         SwitchIncludeRiverDiffusion = iii == 1;
        if (p1.compare("Include diffusion")==0)               SwitchIncludeDiffusion = iii == 1;
        if (p1.compare("Use grain size distribution")==0)     SwitchMulticlass = iii == 1;
        if (p1.compare("Estimate grain size distribution")==0)SwitchEstimateGrainSizeDistribution = iii == 1;
        if (p1.compare("Read grain distribution maps")==0)    SwitchReadGrainSizeDistribution    = iii == 1;
        if (p1.compare("Number of grain size classes (simulated)")==0)  numgrainclasses    = iii ;
        if (p1.compare("Grain size class maps")==0)     GrainMaps  = p;

        if (p1.compare("Flood initial level map")==0)         SwitchFloodInitial     = iii == 1;
        if (p1.compare("Include house storage")==0)            SwitchHouses    =   iii == 1;
        if (p1.compare("Include raindrum storage")==0)         SwitchRaindrum  =   iii == 1;

        if (p1.compare("Include Snowmelt")==0)               SwitchSnowmelt =         iii == 1;
        if (p1.compare("Include Satellite Image")==0)        SwitchImage =            iii == 1;
        if (p1.compare("Hard Surfaces")==0)                  SwitchHardsurface      = iii == 1;
        if (p1.compare("Include road system")==0)            SwitchRoadsystem     = iii == 1;
        if (p1.compare("Include tile drains")==0)            SwitchIncludeTile      = iii == 1;
        if (p1.compare("Include storm drains")==0)           SwitchIncludeStormDrains      = iii == 1;

        if (p1.compare("Include litter interception")==0)    SwitchLitter =           iii == 1;

        if (p1.compare("Include flow barriers")==0)          SwitchFlowBarriers = iii == 1;
        if (p1.compare("Flow barrier table filename")==0)    FlowBarriersFileName = p;
        if (p1.compare("Include buffers")==0)          SwitchBuffers = iii == 1;

        if (p1.compare("Include Sediment traps")==0)         SwitchSedtrap =          iii == 1;
        if (p1.compare("Include compacted")==0)              SwitchInfilCompact =     iii == 1;
        if (p1.compare("Include grass strips")==0)           SwitchGrassStrip =       iii == 1;
        if (p1.compare("Include crusts")==0)                 SwitchInfilCrust =       iii == 1;
        if (p1.compare("Impermeable sublayer")==0)           SwitchImpermeable =      iii == 1;
        if (p1.compare("Matric head files")==0)              SwitchDumphead =         iii == 1;
        if (p1.compare("Geometric mean Ksat")==0)            SwitchGeometric =    	  iii == 1;
        if (p1.compare("Use Water Repellency")==0)           SwitchWaterRepellency  = iii == 1;
        if (p1.compare("Timeplot as PCRaster")==0) {
            SwitchWritePCRtimeplot = iii == 1;
            SwitchWriteCommaDelimited = iii < 1;
        }
        if (p1.compare("Regular runoff output")==0)          SwitchOutputTimeStep =   iii == 1;
        if (p1.compare("User defined output")==0)            SwitchOutputTimeUser =   iii == 1;
        if (p1.compare("Output interval")==0)				 printinterval = iii;
        if (p1.compare("Subsoil drainage")==0)               SwitchDrainage =         iii == 1;
        if (p1.compare("Report point output separate")==0)   SwitchSeparateOutput =   iii == 1;
        if (p1.compare("Report digits out")==0)   ReportDigitsOut = iii;
        if (p1.compare("Report end run")==0)   SwitchEndRun = iii == 1;

        if (p1.compare("KE parameters EQ1")==0)
        {
            QStringList param;
            param = p.split(";",Qt::SkipEmptyParts);
            if (param[0].toInt() == 1)
                KEequationType = KE_EXPFUNCTION;
            KEParamater_a1 = param[1].toDouble();
            KEParamater_b1 = param[2].toDouble();
            KEParamater_c1 = param[3].toDouble();
        }
        if (p1.compare("KE parameters EQ2")==0)
        {
            QStringList param;
            param = p.split(";",Qt::SkipEmptyParts);
            if (param[0].toInt() == 1)
                KEequationType = KE_LOGFUNCTION;
            KEParamater_a2 = param[1].toDouble();
            KEParamater_b2 = param[2].toDouble();
        }
        if (p1.compare("KE parameters EQ3")==0)
        {
            QStringList param;
            param = p.split(";",Qt::SkipEmptyParts);
            if (param[0].toInt() == 1)
                KEequationType = KE_POWERFUNCTION;
            KEParamater_a3 = param[1].toDouble();
            KEParamater_b3 = param[2].toDouble();
        }
        if (p1.compare("KE time based")==0)   SwitchKETimebased = iii == 1;

        if (p1.compare("OutRunoff")==0)         SwitchOutrunoff = iii == 1;
        if (p1.compare("OutWH")==0)             SwitchOutwh = iii == 1;
        if (p1.compare("OutV")==0)              SwitchOutvelo = iii == 1;
        if (p1.compare("OutInterception")==0)   SwitchOutInt = iii == 1;
        if (p1.compare("OutSurfStor")==0)       SwitchOutss = iii == 1;
        if (p1.compare("OutInf")==0)            SwitchOutinf = iii == 1;
        if (p1.compare("OutTileDrain")==0)      SwitchOutTiledrain = iii == 1;
        if (p1.compare("OutTileVolume")==0)     SwitchOutTileVol = iii == 1;
        if (p1.compare("OutDet")==0)            SwitchOutDet = iii == 1;
        if (p1.compare("OutDep")==0)            SwitchOutDep = iii == 1;
        if (p1.compare("OutTC")==0)             SwitchOutTC = iii == 1;
        if (p1.compare("OutConc")==0)           SwitchOutConc = iii == 1;
        if (p1.compare("OutSed")==0)            SwitchOutSed = iii == 1;
        if (p1.compare("OutSL")==0)             SwitchOutSL = iii == 1;
        if (p1.compare("OutSedSS")==0)          SwitchOutSedSS = iii == 1;
        if (p1.compare("OutSedBL")==0)          SwitchOutSedBL = iii == 1;

        if (p1.compare("Erosion map units (0/1/2)")==0)  ErosionUnits = iii;

        InfilMethod = getvalueint("Infil Method");

    }// first loop of runnamelist

    //##########################

    //SwitchUserCores = userCores > 0;
    //qDebug() << userCores;

    // check a few things
    if (InfilMethod == INFIL_GREENAMPT2 || InfilMethod == INFIL_SMITH2)
        SwitchTwoLayer = true;
    else
        SwitchTwoLayer = false;
    if (InfilMethod == INFIL_SWATRE)
    {
        swatreDT = _dt/10; //getvaluedouble("SWATRE internal minimum timestep");

        SwitchGeometric = (getvalueint("Geometric mean Ksat") == 1);
      //  initheadName = getvaluename("inithead");
        // only map name is needed, data is read in swatre lib
        //profileName = getname("profile");//?????????????????????
        // profile map name
    }
    if (SwitchImpermeable)
        SwitchPercolation = false;
    // cannot have both

    //SwitchUse2Phase = SwitchAdvancedSed;
    SwitchUseGrainSizeDistribution = (getvalueint("Use grain size distribution") == 1);
    //qDebug() << SwitchUse2Phase << SwitchAdvancedSed;

    SwitchChannelFlood = true; // always true
    if (!SwitchIncludeChannel)
    {
        SwitchChannelBaseflow = false;
        //SwitchChannelFlood = false;
        SwitchChannelInfil = false;
    }

    // next get the main input directory
    for (j = 0; j < nrrunnamelist; j++)
    {
        QString p1 = runnamelist[j].name;
        QString p = runnamelist[j].value;

        // input ourput dirs and file names
        if (p1.compare("Map Directory")==0)
            inputDir=CheckDir(p);
    }
    // start again and do the rest of the variables, map names etc.
    // choice of options in first loop determines what happens in this loop
    for (j = 0; j < nrrunnamelist; j++)
    {
        QString p1 = runnamelist[j].name;
        QString p = runnamelist[j].value;

        if (InfilMethod == INFIL_SWATRE)
        {
            if (p1.compare("Table Directory")==0)
                SwatreTableDir = CheckDir(p);
            if (p1.compare("Table File")==0)
                SwatreTableName = p;
            initheadName = getvaluename("inithead");
        }

        if (SwitchRainfall)
        {
            if (p1.compare("Rainfall Directory")==0) rainFileDir = CheckDir(p);
            if (p1.compare("Rainfall file")==0) rainFileName = rainFileDir + "/" + p;
        }
//        if (SwitchSnowmelt)
//        {
//            if (p1.compare("Snowmelt Directory")==0) snowmeltFileDir = CheckDir(p);
//            if (p1.compare("Snowmelt file")==0) snowmeltFileName = snowmeltFileDir+ "/" + p;
//        }
        if (SwitchImage)
        {
            if (p1.compare("satImage Directory")==0) satImageFileDir = CheckDir(p);
            if (p1.compare("satImage file")==0) satImageFileName = satImageFileDir + "/" + p;
            //qDebug() << satImageFileName;
        }

        // OUTPUT FILES
        if (p1.compare("Result Directory")==0) {
            resultDir = CheckDir(p, true);
        }

        if (p1.compare("Main results file")==0)
            resultFileName = checkOutputMapName(p, "main results file", 1);
        if (p1.compare("Total Series file")==0)
            totalSeriesFileName = checkOutputMapName(p, "Total Series file", 1);
        if (p1.compare("Filename point output")==0)
            outflowFileName = checkOutputMapName(p, "hydrograph file(s)", 1);
        if (p1.compare("Filename landunit output")==0)
            totalLandunitFileName = checkOutputMapName(p, "Landunit stats output file",1);

        if (p1.compare("Rainfall map")==0)
            rainfallMapFileName = checkOutputMapName(p, "rainfall map", 0);
        if (p1.compare("Interception map")==0)
            interceptionMapFileName = checkOutputMapName(p, "interception map", 0);
        if (p1.compare("Infiltration map")==0)
            infiltrationMapFileName = checkOutputMapName(p, "infiltration map", 0);
        if (p1.compare("Runoff map")==0)
            runoffMapFileName = checkOutputMapName(p, "runoff map", 0);
        if (p1.compare("WH max level map")==0)
            floodWHmaxFileName = checkOutputMapName(p, "WH max level map",0);
        if (p1.compare("Max Velocity")==0)
            floodMaxVFileName = checkOutputMapName(p, "Max Velocity",0);
        if (p1.compare("Max Momentum")==0)
            floodMaxVHFileName = checkOutputMapName(p, "Max Momentum",0);

        if (SwitchIncludeChannel)
        {
            if (p1.compare("Channel discharge map")==0)
                channelDischargeMapFileName = checkOutputMapName(p, "Channel discharge map", 0);
            if (p1.compare("Channel Max Q")==0)
                floodMaxQFileName =  p = checkOutputMapName(p, "channel max discharge",0);
            if (p1.compare("Channel Max WH")==0)
                floodMaxChanWHFileName =  p = checkOutputMapName(p, "channel max water height",0);
        }
        if(SwitchErosion)
        {
            if (p1.compare("Erosion map")==0)
                totalErosionFileName = checkOutputMapName(p, "detachment map",0);
            if (p1.compare("Deposition map")==0)
                totalDepositionFileName = checkOutputMapName(p, "deposition map",0);
            if (p1.compare("Soilloss map")==0)
                totalSoillossFileName = checkOutputMapName(p, "soil loss map",0);
            if (p1.compare("Channel detachment map")==0)
                totalChanErosionFileName = checkOutputMapName(p, "Channel detachment map",0);
            if (p1.compare("Channel deposition map")==0)
                totalChanDepositionFileName = checkOutputMapName(p, "Channel deposition map",0);
        }

        if(SwitchChannelFlood)
        {
            if (p1.compare("Flood time map")==0)
                floodTimeFileName = checkOutputMapName(p, "flood time map",0);
            if (p1.compare("Flood stats")==0)
                floodStatsFileName =  p = checkOutputMapName(p, "flood statistics file",1);
            if (p1.compare("Flood start time")==0)
                floodFEWFileName =  p = checkOutputMapName(p, "flood start time",0);
        }

    }

    if(SwitchSnowmelt) {
        SwitchRainfall = false;
        snowmeltFileName = rainFileName;
    }

    SwitchResultDatetime = getvalueint("Result datetime") == 1;
    if (SwitchResultDatetime) {
        QDir(resultDir).mkpath(QString("res"+op.timeStartRun+"/"));
        resultDir = resultDir + QString("res"+op.timeStartRun+"/");
    }


    Outrunoff = "ro";
    Outconc   = "conc";
    Outwh     = "wh";
    Outrwh    = "";
    OutInt    = "int";
    Outtc     = "tc";
    Outeros   = "det";
    Outdepo   = "dep";
    Outvelo   = "v";
    Outinf    = "inf";
    Outss     = "sstor";
    Outchvol  = "";
    OutTiledrain = "Qtile";
    OutTileVol = "Voltile";
    OutTileV = "Vtile";
    OutHmx  = "";
    OutQf  = "Qf";
    OutVf  = "";
    OutHmxWH  = "";
    OutSL  = "sloss";
    OutSed  = "sed";
    OutSedSS  = "sedSS";
    OutSedBL  = "sedBL";

}
//------------------------------------------------------------------------------
void TWorld::GetRunFile(void)
{
    QFile fin(temprunname);

    if (!fin.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        ErrorString = "Cannot open runfile: " + temprunname;
        throw 2;
    }

    for (int i = 0; i < NUMNAMES; i++)
    {
        runnamelist[i].name.clear();
        runnamelist[i].value.clear();
    }
    nrrunnamelist = 0;

    while (!fin.atEnd())
    {
        QString S = fin.readLine();
        if (!S.trimmed().isEmpty())
        {
            if (S.contains("="))
            {
                QStringList SL = S.split(QRegExp("="));
                runnamelist[nrrunnamelist].name = SL[0].trimmed();
                runnamelist[nrrunnamelist].value = SL[1].trimmed();

                nrrunnamelist++;
            }
        }
    }
}
//------------------------------------------------------------------------------
