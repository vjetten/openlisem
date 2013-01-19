
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
- QString TWorld::CheckDir(QString p) \n
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
QString TWorld::CheckDir(QString p)
{
   /*
 if (QDir(p).exists())
 {
  if (!p.endsWith("/"))
   p = p + '/';
  return p;
 }
 else
 {
  if (!QDir(p).exists())
   ErrorString = p1 +": "+p+ " is not an existing directory!";
  throw 1;
 }
 return "";
 */
   //	p.replace("/","\\");
   //	//	if (!p.endsWith("/"))
   //	if (!p.endsWith("\\"))
   //		p = p + "\\";
   //	if (QDir(p).exists())
   //		return p;
   //	else
   //		return "";

   QString path;
   path = QDir(p).fromNativeSeparators(p);
   path = QDir(path).absoluteFilePath(path);

   if (!path.endsWith("/"))
      path = path + '/';

   if (!QDir(path).exists())
      path.clear();

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
void TWorld::ParseRunfileData(void)
{
   int j=0;

   // do all switches first
   for (j = 0; j < nrrunnamelist; j++)
   {
      int iii = runnamelist[j].value.toInt();
      QString p1 = runnamelist[j].name;
      QString p = runnamelist[j].value;

      //fprintf(fout,"%s=%s\n",(const char *)p1.toLatin1(),(const char *)p.toLatin1());
      /*
  // main lisem types
  if (p1.compare("LISEM Type")==0)
  {
   SwitchWheelAsChannel = iii == LISEMWHEELTRACKS;
   SwitchMulticlass = iii == LISEMMULTICLASS;
   SwitchNutrients = iii == LISEMNUTRIENTS;
   SwitchGullies = iii == LISEMGULLIES;
  }
*/
      //options in the main code, order is not important
      if (p1.compare("No Erosion simulation")==0)          SwitchErosion =          iii == 0;
      if (p1.compare("Include main channels")==0)          SwitchIncludeChannel =   iii == 1;
      if (p1.compare("Include channel infil")==0)          SwitchChannelInfil =     iii == 1;
      if (p1.compare("Include channel baseflow")==0)       SwitchChannelBaseflow =  iii == 1;
      if (p1.compare("Include channel flooding")==0)       SwitchChannelFlood    =   iii == 1;
      if (p1.compare("Include tile drains")==0)            SwitchIncludeTile    =   iii == 1;
      //houses
      if (p1.compare("Include house storage")==0)            SwitchHouses    =   iii == 1;
      if (p1.compare("Include raindrum storage")==0)         SwitchRaindrum  =   iii == 1;
      if (p1.compare("All water and sediment to outlet")==0) SwitchAllinChannel  =  iii == 1;
      SwitchAllinChannel = true;
      //VJ 100526 always true in old LISEM

      if (p1.compare("Include Rainfall")==0)               SwitchRainfall =         iii == 1;
      if (p1.compare("Include Snowmelt")==0)               SwitchSnowmelt =         iii == 1;
      if (p1.compare("Alternative flow detachment")==0)    SwitchAltErosion =       iii == 1;
      if (p1.compare("Simple depression storage")==0)      SwitchSimpleDepression = iii == 1;
      if (p1.compare("Hard Surfaces")==0)                  SwitchHardsurface      = iii == 1;
      if (p1.compare("Limit TC")==0)                       SwitchLimitTC =          iii == 1;
      //if (p1.compare("Limit Deposition TC")==0)            SwitchLimitDepTC =       iii == 1;
      if (p1.compare("Include buffers")==0)                SwitchBuffers =          iii == 1;
      if (p1.compare("Include Sediment traps")==0)         SwitchSedtrap =          iii == 1;
      if (p1.compare("Include wheeltracks")==0)            SwitchInfilCompact =     iii == 1;
      if (p1.compare("Include grass strips")==0)           SwitchGrassStrip =       iii == 1;
      if (p1.compare("Include crusts")==0)                 SwitchInfilCrust =       iii == 1;
      if (p1.compare("Impermeable sublayer")==0)           SwitchImpermeable =      iii == 1;
      if (p1.compare("Matric head files")==0)              SwitchDumphead =         iii == 1;
      if (p1.compare("Geometric mean Ksat")==0)            SwitchGeometric =    		iii == 1;
      if (p1.compare("Use Water Repellency")==0)           SwitchWaterRepellency  = iii == 1;
      if (p1.compare("Runoff maps in l/s/m")==0)           SwitchRunoffPerM =       iii == 1;
      if (p1.compare("Timeseries as PCRaster")==0)         SwitchWritePCRnames =    iii == 1;
      if (p1.compare("Timeseries as CSV")==0)              SwitchWriteCommaDelimited =    iii == 1;
      if (p1.compare("Timeplot as PCRaster")==0)           SwitchWritePCRtimeplot = iii == 1;
      if (p1.compare("Regular runoff output")==0)          SwitchOutputTimeStep =   iii == 1;
      if (p1.compare("User defined output")==0)            SwitchOutputTimeUser =   iii == 1;
      if (p1.compare("Output interval")==0)					  printinterval = iii;
      if (p1.compare("No erosion at outlet")==0)           SwitchNoErosionOutlet =  iii == 1;
      if (p1.compare("Subsoil drainage")==0)               SwitchDrainage =         iii == 1;
      if (p1.compare("Gully infiltration")==0)             SwitchGullyInfil =       iii == 1;
      if (p1.compare("Use initial gully dimensions")==0)   SwitchGullyInit =        iii == 1;
      if (p1.compare("Report point output separate")==0)   SwitchSeparateOutput =   iii == 1;
      if (p1.compare("Report point output for SOBEK")==0)
         SwitchSOBEKoutput =      iii == 1;
      if (p1.compare("SOBEK date string")==0)              SOBEKdatestring = p;
      SOBEKdatestring.remove(10,100);
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

      if (p1.compare("Erosion map units (0/1/2)")==0)      ErosionUnits = iii;
      // VJ 111001
   }

   InfilMethod = getvalueint("Infil Method");
   if (InfilMethod == INFIL_GREENAMPT2 || InfilMethod == INFIL_SMITH2)
      SwitchTwoLayer = true;

   for (j = 0; j < nrrunnamelist; j++)
   {
      QString p1 = runnamelist[j].name;
      QString p = runnamelist[j].value;

      // input ourput dirs and file names
      if (p1.compare("Map Directory")==0)
         inputDir=CheckDir(p);
      if (p1.compare("Result Directory")==0)
         resultDir = CheckDir(p);

      if (InfilMethod == INFIL_SWATRE)
      {
         if (p1.compare("Table Directory")==0)
         {
            SwatreTableDir = CheckDir(p);
         }
         if (p1.compare("Table File")==0)
         {
            SwatreTableName = p;
         }
      }
      if (p1.compare("Main results file")==0)
      {
         resultFileName =  p;
         if (p.isEmpty())
         {
            ErrorString = "Please give a name for the main results file";
            throw 1;
         }
      }
      if (p1.compare("Filename point output")==0)
      {
         outflowFileName =  p;
         if (p.isEmpty())
         {
            ErrorString = "Please give a name for the hydrograph file";
            throw 1;
         }
      }
      if (p1.compare("Filename landunit output")==0)
      {
         totalLandunitFileName =  p;
         if (p.isEmpty())
         {
            ErrorString = "Please give a name for the Landunit stats output file";
            throw 1;
         }
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
      if (SwitchErosion)
      {
         if (p1.compare("Erosion map")==0)
         {
            totalErosionFileName = p;
            if (p.isEmpty())
            {
               ErrorString = "Please give a name for the detachment map";
               throw 1;
            }
            if (!p.contains("."))
            {
               ErrorString = "Please give a name with an extention (such as \".map\")";
               throw 1;
            }
         }
         if (p1.compare("Deposition map")==0)
         {
            totalDepositionFileName =  p;
            if (p.isEmpty())
            {
               ErrorString = "Please give a name for the deposition map";
               throw 1;
            }
            if (!p.contains("."))
            {
               ErrorString = "Please give a name with an extention (such as \".map\")";
               throw 1;
            }
         }
         if (p1.compare("Soilloss map")==0)
         {
            totalSoillossFileName =  p;
            if (p.isEmpty())
            {
               ErrorString = "Please give a name for the soil loss map";
               throw 1;
            }
            if (!p.contains("."))
            {
               ErrorString = "Please give a name with an extention (such as \".map\")";
               throw 1;
            }
         }
         // resultDir is added in report operation
      }

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
   }

   if (InfilMethod == INFIL_SWATRE)
   {
      swatreDT = getvaluedouble("SWATRE internal minimum timestep");
      SwitchGeometric = (getvalueint("Geometric mean Ksat") == 1);
      initheadName = getvaluename("inithead");
      // only map name is needed, data is read in swatre lib
      //profileName = getname("profile");//?????????????????????
      // profile map name
   }

   if (!SwitchIncludeTile && outputcheck.count() > 11)
      outputcheck[11] = "0";

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
