
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

#include "model.h"


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

    if (!path.endsWith("/"))
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

    // do all switches (checkbox options) first
    for (j = 0; j < nrrunnamelist; j++)
    {
        int iii = runnamelist[j].value.toInt();
        QString p1 = runnamelist[j].name;
        QString p = runnamelist[j].value;

        //options in the main code, order is not important
        if (p1.compare("Include Erosion simulation")==0)     SwitchErosion =          iii == 1;
        if (p1.compare("Include main channels")==0)          SwitchIncludeChannel =   iii == 1;
        if (p1.compare("Include channel infil")==0)          SwitchChannelInfil     = iii == 1;
        if (p1.compare("Include channel baseflow")==0)       SwitchChannelBaseflow  = iii == 1;
        if (p1.compare("Include channel flooding")==0)       SwitchChannelFlood     = iii == 1;
        if (p1.compare("Include rainfall flooding")==0)      SwitchRainfallFlood     = iii == 1;
        if (p1.compare("Include road system")==0)            SwitchRoadsystem     = iii == 1;
        if (p1.compare("Include tile drains")==0)            SwitchIncludeTile      = iii == 1;
     //   if (p1.compare("Flood method explicit")==0)          SwitchFloodExplicit    = iii == 1;
        if (p1.compare("Flood method SWOF2D order 1")==0)    SwitchFloodSWOForder1  = iii == 1;
        if (p1.compare("Flood method SWOF2D order 2")==0)    SwitchFloodSWOForder2  = iii == 1;

        if (p1.compare("Close boundary Kin Wave 2D")==0)     SwitchClosedBoundaryOF = iii == 1;

//        if (p1.compare("D90 for distribution")==0)          distD90 = p.toDouble();
//        if (p1.compare("D50 for distribution")==0)          distD50 = p.toDouble();
    //    if (p1.compare("OF method")==0)                       OF_Method     = iii;

        if (p1.compare("River BL method")==0)                 R_BL_Method     = iii;
        if (p1.compare("River SS method")==0)                 R_SS_Method     = iii;
      //  if (p1.compare("Estimate d90")==0)                    SwithEstimated90     = iii == 1;
        if (p1.compare("Use material depth")==0)              SwitchUseMaterialDepth  = iii == 1;

        if (p1.compare("Estimate grain size distribution")==0)SwitchEstimateGrainSizeDistribution = iii == 1;

        if (p1.compare("Read grain distribution maps")==0)    SwitchReadGrainSizeDistribution    = iii == 1;

        if (p1.compare("Advanced sediment configuration")==0)
        {
            if(iii == 2)
            {
                SwitchUse2Layer = true;
                SwitchUseGrainSizeDistribution = true;

            }else if(iii == 1)
            {
                SwitchUse2Layer = true;
                SwitchUseGrainSizeDistribution = false;
            }else
            {
                SwitchUse2Layer = false;
                SwitchUseGrainSizeDistribution = false;
            }
        }

        if(SwitchEstimateGrainSizeDistribution)
        {
                if (p1.compare("Number of grain size classes (simulated)")==0)  numgrainclasses    = iii ;
                //if (p1.compare("Grain size distribution type")==0)    GrainSizeDistributionType     = iii;
        }else
            if(SwitchReadGrainSizeDistribution)
        {
               // if (p1.compare("Number of grain size classes (maps)")==0)  numgrainclasses     = iii;
                if (p1.compare("Grain size class maps")==0)    GrainMaps    = p;
        }

        if (p1.compare("Detachment efficiency")==0)          SwitchEfficiencyDET = iii;
    //    if (p1.compare("Detachment stoniness")==0)           SwitchStoninessDET   = iii == 1;

        if (p1.compare("Flood initial level map")==0)          SwitchFloodInitial     = iii == 1;
   //     if (p1.compare("Flood calc as watershed")==0)          SwitchWatershed     = iii == 1;
        if (p1.compare("Flood sediment transport method")==0)  SwitchFloodSedimentMethod     = iii == 1;

        //        if (p1.compare("Minimum stats flood height")==0)     SwitchLevees     = iii == 1;
        //houses
        if (p1.compare("Include house storage")==0)            SwitchHouses    =   iii == 1;
        if (p1.compare("Include raindrum storage")==0)         SwitchRaindrum  =   iii == 1;
        if (p1.compare("All water and sediment to outlet")==0) SwitchAllinChannel  =  iii == 1;
        SwitchAllinChannel = true;
        //VJ 100526 always true in old LISEM

        if (p1.compare("Include Rainfall")==0)               SwitchRainfall =         iii == 1;
        if (p1.compare("Include Snowmelt")==0)               SwitchSnowmelt =         iii == 1;
        //  if (p1.compare("Alternative flow detachment")==0)    SwitchAltErosion =       iii == 1;
        if (p1.compare("Simple depression storage")==0)      SwitchSimpleDepression = iii == 1;
        if (p1.compare("Hard Surfaces")==0)                  SwitchHardsurface      = iii == 1;
        if (p1.compare("Limit TC")==0)                       SwitchLimitTC =          iii == 1;
        //if (p1.compare("Limit Deposition TC")==0)            SwitchLimitDepTC =       iii == 1;
        if (p1.compare("Include litter interception")==0)    SwitchLitter =          iii == 1;

        if (p1.compare("Include flow barriers")==0)          SwitchFlowBarriers = iii == 1;
        if (p1.compare("Flow barrier table filename")==0)    FlowBarriersFileName = p;

        if (p1.compare("Include buffers")==0)                SwitchBuffers =          iii == 1;
        if (p1.compare("Include Sediment traps")==0)         SwitchSedtrap =          iii == 1;
        if (p1.compare("Include wheeltracks")==0)            SwitchInfilCompact =     iii == 1;
        if (p1.compare("Include compacted")==0)              SwitchInfilCompact =     iii == 1;
        if (p1.compare("Include grass strips")==0)           SwitchGrassStrip =       iii == 1;
        if (p1.compare("Include crusts")==0)                 SwitchInfilCrust =       iii == 1;
        if (p1.compare("Impermeable sublayer")==0)           SwitchImpermeable =      iii == 1;
        if (p1.compare("Include percolation")==0)            SwitchPercolation =      iii == 1;
        if (p1.compare("Matric head files")==0)              SwitchDumphead =         iii == 1;
        if (p1.compare("Geometric mean Ksat")==0)            SwitchGeometric =    		iii == 1;
        if (p1.compare("Use Water Repellency")==0)           SwitchWaterRepellency  = iii == 1;
        if (p1.compare("Runoff maps in l/s/m")==0)           SwitchRunoffPerM =       iii == 1;
        if (p1.compare("Timeseries as PCRaster")==0)         SwitchWritePCRnames =    iii == 1;
        if (p1.compare("Timeseries as CSV")==0)              SwitchWriteCommaDelimited =    iii == 1;
        if (p1.compare("Timeplot as PCRaster")==0)           SwitchWritePCRtimeplot = iii == 1;
        if (p1.compare("Regular runoff output")==0)          SwitchOutputTimeStep =   iii == 1;
        if (p1.compare("User defined output")==0)            SwitchOutputTimeUser =   iii == 1;
        if (p1.compare("Output interval")==0)				 printinterval = iii;
        // if (p1.compare("No erosion at outlet")==0)           SwitchNoErosionOutlet =  iii == 1;
        if (p1.compare("Subsoil drainage")==0)               SwitchDrainage =         iii == 1;
        if (p1.compare("Gully infiltration")==0)             SwitchGullyInfil =       iii == 1;
        if (p1.compare("Use initial gully dimensions")==0)   SwitchGullyInit =        iii == 1;
        if (p1.compare("Report point output separate")==0)   SwitchSeparateOutput =   iii == 1;
        if (p1.compare("Report point output for SOBEK")==0)  SwitchSOBEKoutput = iii == 1;
        if (p1.compare("SOBEK date string")==0)
        {
            SOBEKdatestring = p;
            SOBEKdatestring.remove(10,100);
        }
        if (p1.compare("Report digits out")==0)   ReportDigitsOut = iii;


        if (p1.compare("Use canopy storage map")==0)   	   SwitchInterceptionLAI =  iii == 0;

        if (p1.compare("KE parameters EQ1")==0)
        {
            QStringList param;
            param = p.split(",",QString::SkipEmptyParts);
            if (param[0].toInt() == 1)
                KEequationType = KE_EXPFUNCTION;
            KEParamater_a1 = param[1].toDouble();
            KEParamater_b1 = param[2].toDouble();
            KEParamater_c1 = param[3].toDouble();
        }
        if (p1.compare("KE parameters EQ2")==0)
        {
            QStringList param;
            param = p.split(",",QString::SkipEmptyParts);
            if (param[0].toInt() == 1)
                KEequationType = KE_LOGFUNCTION;
            KEParamater_a2 = param[1].toDouble();
            KEParamater_b2 = param[2].toDouble();
        }
        if (p1.compare("KE parameters EQ3")==0)
        {
            QStringList param;
            param = p.split(",",QString::SkipEmptyParts);
            if (param[0].toInt() == 1)
                KEequationType = KE_POWERFUNCTION;
            KEParamater_a3 = param[1].toDouble();
            KEParamater_b3 = param[2].toDouble();
        }
        if (p1.compare("KE time based")==0)   SwitchKETimebased = iii == 1;

        if (p1.compare("CheckOutputMaps")==0)   outputcheck = p.split(",");
        // outputcheck is a string with 0,1,0,1,... etc

        if (p1.compare("Erosion map units (0/1/2)")==0)  ErosionUnits = iii;

        InfilMethod = getvalueint("Infil Method");

    }// first loop of runnamelist
qDebug() << SwitchClosedBoundaryOF;
    if (SwitchFloodSWOForder2)
    {
        SwitchFloodSWOForder1 = false;
    }
    if (SwitchFloodSWOForder1)
    {
        SwitchFloodSWOForder2 = false;
    }

    SwitchFlood1D2DCoupling = getvalueint("Flooding 1D2D coupling");
    SwitchKinematic2D = std::max(getvalueint("Routing Kin Wave 2D"), 1);
    CourantKin = getvaluedouble("Courant Kin Wave 2D");

    TimestepKinMin = getvaluedouble("Timestep Kin Wave 2D");

    OF_Method = (SwitchUseGrainSizeDistribution? OFHAIRSINEROSE : OFGOVERS);

    if (SwitchChannelFlood && !SwitchFloodSWOForder1 && !SwitchFloodSWOForder2)
        SwitchFloodSWOForder1 = true;

    // check a few things
    if (InfilMethod == INFIL_GREENAMPT2 || InfilMethod == INFIL_SMITH2)
        SwitchTwoLayer = true;
    else
        SwitchTwoLayer = false;
    if (InfilMethod == INFIL_SWATRE)
    {
        swatreDT = getvaluedouble("SWATRE internal minimum timestep");
        SwitchGeometric = (getvalueint("Geometric mean Ksat") == 1);
      //  initheadName = getvaluename("inithead");
        // only map name is needed, data is read in swatre lib
        //profileName = getname("profile");//?????????????????????
        // profile map name
    }

    if (SwitchImpermeable)
        SwitchPercolation = false;
    // cannot have both


    // fill up outputcheck for older runfiles
    if (outputcheck.count() < 20)
        for (int k = outputcheck.count(); k < 20; k++)
            outputcheck << "0";

//    if (!SwitchIncludeTile && outputcheck.count() > 11)
//        outputcheck[11] = "0";  //????????????

    if (!SwitchIncludeChannel)
    {
        SwitchChannelBaseflow = false;
        SwitchChannelFlood = false;
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

        // input ourput dirs and file names
    //    if (p1.compare("Map Directory")==0)
    //        inputDir=CheckDir(p);

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
        if (SwitchSnowmelt)
        {
            if (p1.compare("Snowmelt Directory")==0) snowmeltFileDir = CheckDir(p);
            if (p1.compare("Snowmelt file")==0) snowmeltFileName = snowmeltFileDir + p;
        }

        // OUTPUT FILES
        if (p1.compare("Result Directory")==0)
            resultDir = CheckDir(p, true);
           // QFileInfo("resultDir");
           // WHAT IS THIS ?????????????????

        if (p1.compare("Main results file")==0)
            resultFileName = checkOutputMapName(p, "main results file", 1);
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
        if (p1.compare("Runoff fraction map")==0)
            runoffFractionMapFileName = checkOutputMapName(p, "runoff fraction map", 0);
        if (p1.compare("Channel discharge map")==0)
            channelDischargeMapFileName = checkOutputMapName(p, "Channel discharge map", 0);
        if (p1.compare("WH max level map")==0)
            floodWHmaxFileName = checkOutputMapName(p, "WH max level map",0);

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
            if (p1.compare("Flood level map")==0)
                floodLevelFileName = checkOutputMapName(p, "flood level map",0);
            if (p1.compare("Flood time map")==0)
                floodTimeFileName = checkOutputMapName(p, "flood time map",0);
            if (p1.compare("Flood stats")==0)
                floodStatsFileName =  p = checkOutputMapName(p, "flood statistics file",1); ;
            if (p1.compare("Channel Max Q")==0)
                floodMaxQFileName =  p = checkOutputMapName(p, "channel max discharge",0); ;
            if (p1.compare("Channel Max WH")==0)
                floodMaxChanWHFileName =  p = checkOutputMapName(p, "channel max water height",0); ;
            if (p1.compare("Flood start time")==0)
                floodFEWFileName =  p = checkOutputMapName(p, "flood start time",0); ;
            if (p1.compare("Flood Max V")==0)
                floodMaxVFileName = checkOutputMapName(p, "flood max V",0);
        }

        // output map timeseries, standard names, to avoid unreadable pcraster names
        if (p1.compare("OUTRUNOFF")==0)  Outrunoff = GetName(p);
        if (p1.compare("OUTCONC"  )==0)  Outconc   = GetName(p);
        if (p1.compare("OUTWH"    )==0)  Outwh     = GetName(p);
        if (p1.compare("OUTRWH"   )==0)  Outrwh    = GetName(p);
        if (p1.compare("OUTTC"    )==0)  Outtc     = GetName(p);
        if (p1.compare("OUTEROS"  )==0)  Outeros   = GetName(p);
        if (p1.compare("OUTDEPO"  )==0)  Outdepo   = GetName(p);
        if (p1.compare("OUTVELO"  )==0)  Outvelo   = GetName(p);
        if (p1.compare("OUTINF"   )==0)  Outinf    = GetName(p);
        if (p1.compare("OUTSS"    )==0)  Outss     = GetName(p);
        if (p1.compare("OUTCHVOL" )==0)  Outchvol  = GetName(p);
        if (p1.compare("OUTTILED" )==0)  OutTiledrain  = GetName(p);
        if (p1.compare("OUTHMX"   )==0)    OutHmx  = GetName(p);
        if (p1.compare("OUTQF"    )==0)     OutQf  = GetName(p);
        if (p1.compare("OUTVF"    )==0)     OutVf  = GetName(p);
        if (p1.compare(""
                       ""
                       "OUTHMXWH" )==0)  OutHmxWH  = GetName(p);
        if (p1.compare("OUTSOILLOSS" )==0)  OutSL  = GetName(p);
        if (p1.compare("OUTSED" )==0)    OutSed  = GetName(p);
    }


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
