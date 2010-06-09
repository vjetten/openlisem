/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

#include "model.h"


//---------------------------------------------------------------------------
QString TWorld::getvaluename(QString vname)
{
	for (int i = 0; i < nrnamelist; i++)
	{
		if(vname.toUpper() == namelist[i].name.toUpper())
		{
			if (namelist[i].value.trimmed().isEmpty())
			{
//				DEBUG(vname.toUpper() + " " + namelist[i].name.toUpper()+"<==");
				ErrorString = "Filename not found for : " + vname;
				throw 1;
			}
			else
			{
				QFileInfo info(namelist[i].value);
				//namelist[i].value = inputDir + info.fileName();
				return inputDir + info.fileName();//namelist[i].value;
			}
		}
	}
	ErrorString = "wrong internal map ID "+vname;
	throw 3;
}
//---------------------------------------------------------------------------
double TWorld::getvaluedouble(QString vname)
{
	for (int i = 0; i < nrnamelist; i++)
		if(vname.toUpper() == namelist[i].name.toUpper())
		{
			return namelist[i].value.toDouble();
		}

	ErrorString = "wrong internal map ID "+vname;
	throw 3;
}
//---------------------------------------------------------------------------
int TWorld::getvalueint(QString vname)
{
	for (int i = 0; i < nrnamelist; i++)
		if(vname.toUpper() == namelist[i].name.toUpper())
		{
			return namelist[i].value.toInt();
		}

	ErrorString = "wrong internal map ID "+vname;
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
	QString sss;
	sss.setNum(nrnamelist);
}
//---------------------------------------------------------------------------
QString TWorld::CheckDir(QString p, QString p1)
{
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
void TWorld::ParseInputData()
{
	int j=0;
	//FILE *fout = fopen("c:\\try.txt","w");

	// do all switches first
	for (j = 0; j < nrnamelist; j++)
	{
		int iii = namelist[j].value.toInt();
		QString p1 = namelist[j].name;
		QString p = namelist[j].value;

		//fprintf(fout,"%s=%s\n",(const char *)p1.toLatin1(),(const char *)p.toLatin1());

		// main lisem types
		if (p1.compare("LISEM Type")==0)
		{
			SwitchWheelAsChannel = iii == LISEMWHEELTRACKS;
			SwitchMulticlass = iii == LISEMMULTICLASS;
			SwitchNutrients = iii == LISEMNUTRIENTS;
			SwitchGullies = iii == LISEMGULLIES;
		}

		//options in the main code, order is not important
		if (p1.compare("No Erosion simulation")==0)          SwitchErosion =          iii == 0;
		if (p1.compare("Include main channels")==0)          SwitchIncludeChannel =   iii == 1;
		if (p1.compare("Include channel infil")==0)          SwitchChannelInfil =     iii == 1;
		if (p1.compare("Include channel baseflow")==0)       SwitchChannelBaseflow =  iii == 1;
		if (p1.compare("All water and sediment to outlet")==0) SwitchAllinChannel    =  iii == 1;
		SwitchAllinChannel = true;
		//VJ 100526 always true in old LISEM

		if (p1.compare("Include snowmelt")==0)               SwitchSnowmelt =         iii == 1;
		if (p1.compare("Alternative flow detachment")==0)    SwitchAltErosion =       iii == 1;
		if (p1.compare("Simple depression storage")==0)      SwitchSimpleDepression = iii == 1;
		if (p1.compare("Hard Surfaces")==0)                  SwitchHardsurface      = iii == 1;
		if (p1.compare("Include buffers")==0)                SwitchBuffers =          iii == 1;
		if (p1.compare("Include Sediment traps")==0)         SwitchSedtrap =          iii == 1;
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
		if (p1.compare("Use canopy storage map")==0)   	   SwitchInterceptionLAI =  iii == 0;

		if (p1.compare("CheckOutputMaps")==0)   outputcheck = p.split(",");
		if (p1.compare("Output interval")==0)   printinterval = iii;
	}

	for (j = 0; j < nrnamelist; j++)
	{
		QString p1 = namelist[j].name;
		QString p = namelist[j].value;

		// input ourput dirs and file names
		if (p1.compare("Map Directory")==0) inputDir=CheckDir(p, p1);
		if (p1.compare("Result Directory")==0) p = resultDir = CheckDir(p, p1);

		//   if (p1.compare("Table Directory")==0) tableDir = CheckDir(p, p1);
		// move to swatre later when infiltration method is known!
		if (p1.compare("Main results file")==0) resultFileName = p;
		if (p1.compare("Filename point output")==0) outflowFileName =  p;
		// resultDir is added in report operation

		if (p1.compare("Rainfall Directory")==0) rainFileDir = CheckDir(p, p1);
		if (p1.compare("Rainfall file")==0) rainFileName = rainFileDir + p;
		if (SwitchErosion)
		{
			if (p1.compare("Erosion map")==0) totalErosionFileName =  p;
			if (p1.compare("Deposition map")==0) totalDepositionFileName =  p;
			if (p1.compare("Soilloss map")==0) totalSoillossFileName =  p;
			// resultDir is added in report operation
		}
		if (SwitchSnowmelt)
		{
			if (p1.compare("Snowmelt Directory")==0) snowmeltFileDir = CheckDir(p, p1);
			if (p1.compare("Snowmelt file")==0) snowmeltFileName = snowmeltFileDir + p;
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
	}
	//fclose(fout);
//	if (SwitchSedtrap)
	//	SwitchBuffers = true;

}
//------------------------------------------------------------------------------
