/*
 * runfile.cpp
 *
 *  Created on: Feb 23, 2010
 *      Author: jetten
 */

#include "model.h"

#define addslash(_s) (!_s.endsWith("\\")? _s += "\\" : _s)

//---------------------------------------------------------------------------
QString TWorld::getvaluename(const char *vname)
{
   QString vvname(vname);

   for (int i = 0; i < nrnamelist; i++)
     if(vvname.toUpper().compare(namelist[i].name.toUpper()) == 0)
       {
         if (namelist[i].name.trimmed().isEmpty())
         {
             ErrorString = "Map filename not found for : " + vvname;
             throw 1;
         }
         else
             return namelist[i].value;
       }
   ErrorString = "wrong internal map ID "+vvname;
   throw 3;
}
//---------------------------------------------------------------------------
double TWorld::getvaluedouble(const char *vname)
{
   QString vvname = vname;

   for (int i = 0; i < nrnamelist; i++)
     if(vvname.toUpper().compare(namelist[i].name.toUpper()) == 0)
     {
        return namelist[i].value.toDouble();
     }

   ErrorString = "wrong internal map ID "+vvname;
   throw 3;
}
//---------------------------------------------------------------------------
int TWorld::getvalueint(const char *vname)
{
   QString vvname = vname;

   for (int i = 0; i < nrnamelist; i++)
     if(vvname.toUpper().compare(namelist[i].name.toUpper()) == 0)
     {
        return namelist[i].value.toInt();
     }

   ErrorString = "wrong internal map ID "+vvname;
   throw 3;
}
//------------------------------------------------------------------------------
void TWorld::GetRunFile()
{
    QFile fin(temprunname);

    if (!fin.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        ErrorString = "Cannot open runfile: " + temprunname;
        throw 2;
    }

    for (int i = 0; i < NUMNAMES; i++)
    {
        namelist[i].name.clear();
        namelist[i].value.clear();
    }
    nrnamelist = 0;

    while (!fin.atEnd())
    {
        QString S = fin.readLine();
        if (!S.trimmed().isEmpty())
        {
          if (S.contains("="))
            {
              QStringList SL = S.split(QRegExp("="));
              namelist[nrnamelist].name = SL[0].trimmed();
              namelist[nrnamelist].value = SL[1].trimmed();
              nrnamelist++;
            }
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::ParseInputData()
{
    int j=0;

    for (j = 0; j < nrnamelist; j++)
    {
         int iii = namelist[j].value.toInt();
         QString p1 = namelist[j].name;
         QString p = namelist[j].value;

          // main lisem types
          if (p1.compare("LISEM Type")==0)
          {
              SwitchWheelAsChannel = iii == LISEMWHEELTRACKS;
              SwitchMulticlass = iii == LISEMMULTICLASS;
              SwitchNutrients = iii == LISEMNUTRIENTS;
              SwitchGullies = iii == LISEMGULLIES;
          }

          // input ourput dirs and file names
          if (p1.compare("Map Directory")==0) inputDir=addslash(p);
          if (p1.compare("Result Directory")==0) resultDir = addslash(p);
          if (p1.compare("Table Directory")==0) tableDir = addslash(p);
          if (p1.compare("Main results file")==0) resultFileName = p + resultDir;
          if (p1.compare("Erosion map")==0) totalErosionFileName =  p + resultDir;
          if (p1.compare("Deposition map")==0) totalDepositionFileName =  p + resultDir;
          if (p1.compare("Soilloss map")==0) totalSoillossFileName =  p + resultDir;
          if (p1.compare("Filename point output")==0) outflowFileName =  p + resultDir;
          if (p1.compare("Rainfall Directory")==0) rainFileDir = addslash(p);
          if (p1.compare("Rainfall file")==0) rainFileName = rainFileDir + p;
          if (p1.compare("Snowmelt Directory")==0) snowmeltFileDir = addslash(p);
          if (p1.compare("Snowmelt file")==0) snowmeltFileName = snowmeltFileDir + p;

          //options in the main code, order is not important
          if (p1.compare("No Erosion simulation")==0)          SwitchErosion =          iii == 0;
          if (p1.compare("Include main channels")==0)          SwitchIncludeChannel =   iii == 1;
          if (p1.compare("Include channel infil")==0)          SwitchChannelInfil =     iii == 1;
          if (p1.compare("Include channel baseflow")==0)       SwitchChannelBaseflow =  iii == 1;
          if (p1.compare("Include snowmelt")==0)               SwitchSnowmelt =         iii == 1;
          if (p1.compare("Alternative flow detachment")==0)    SwitchAltErosion =       iii == 1;
          if (p1.compare("Simple depression storage")==0)      SwitchSimpleDepression = iii == 1;
          if (p1.compare("Hard Surfaces")==0)                  SwitchHardsurface      = iii == 1;
          if (p1.compare("Include buffers")==0)                SwitchBuffers =          iii == 1;
          if (p1.compare("Include wheeltracks")==0)            SwitchInfilCompact =     iii == 1;
          if (p1.compare("Include grass strips")==0)           SwitchInfilGrass =       iii == 1;
          if (p1.compare("Include crusts")==0)                 SwitchInfilCrust =       iii == 1;
          if (p1.compare("Impermeable sublayer")==0)           SwitchImpermeable =      iii == 1;
          if (p1.compare("Matric head files")==0)              SwitchDumphead =         iii == 1;
          if (p1.compare("Geometric mean Ksat")==0)            SwitchGeometricMean =    iii == 1;
          if (p1.compare("Runoff maps in l/s/m")==0)           SwitchRunoffPerM =       iii == 1;
          if (p1.compare("Timeseries as PCRaster")==0)         SwitchWritePCRnames =    iii == 1;
          if (p1.compare("Timeplot as PCRaster")==0)           SwitchWritePCRtimeplot = iii == 1;
          if (p1.compare("Regular runoff output")==0)          SwitchOutputTimeStep =   iii == 1;
          if (p1.compare("User defined output")==0)            SwitchOutputTimeUser =   iii == 1;
          if (p1.compare("No erosion at outlet")==0)           SwitchNoErosionOutlet =  iii == 1;
          if (p1.compare("Subsoil drainage")==0)               SwitchDrainage =         iii == 1;
          if (p1.compare("Gully infiltration")==0)             SwitchGullyInfil =       iii == 1;
          if (p1.compare("Use initial gully dimensions")==0)   SwitchGullyInit =        iii == 1;
          if (p1.compare("Report point output separate")==0)   SwitchSeparateOutput =   iii == 1;
          if (p1.compare("Report point output for SOBEK")==0)  SwitchSOBEKOutput =      iii == 1;
          if (p1.compare("SOBEK date string")==0)              SOBEKdatestring = p;
          SOBEKdatestring.remove(10,100);
          if (p1.compare("Use canopy storage map")==0)   SwitchInterceptionLAI =        iii == 0;

    }

}
//------------------------------------------------------------------------------
