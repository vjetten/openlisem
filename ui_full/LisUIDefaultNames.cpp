/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * Default runfile variable names
 */

#include "lisemqt.h"

//---------------------------------------------------------------------------
// change runfile strings with current interface options
void lisemqt::DefaultRunFile()
{
   int i;
   for (i = 0; i < NUMNAMES; i++)
   {
      defnamelist[i].name.clear();
      defnamelist[i].value.clear();
   }
   // clear first

   // runfile has structure:
   // name=value
   i = 0;
	defnamelist[i++].name = QString("[openLISEM runfile version 4]");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[LISEM main type]");
	defnamelist[i++].name = QString("LISEM Type");
	defnamelist[i++].name = QString("");
//	defnamelist[i++].name = QString("[Work Directory]");
//	defnamelist[i++].name = QString("WorkDir");
//	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Input]");
	defnamelist[i++].name = QString("Map Directory");
	defnamelist[i++].name = QString("Include Rainfall");
	defnamelist[i++].name = QString("Rainfall Directory");
	defnamelist[i++].name = QString("Rainfall file");
	defnamelist[i++].name = QString("Include Snowmelt");
	defnamelist[i++].name = QString("Snowmelt Directory");
	defnamelist[i++].name = QString("Snowmelt file");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Output]");
	defnamelist[i++].name = QString("Result Directory");
   defnamelist[i].value = QString("totals.txt");
   defnamelist[i++].name = QString("Main results file");
   defnamelist[i].value = QString("hydrograph.csv");
   defnamelist[i++].name = QString("Filename point output");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Report point output separate");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Report point output for SOBEK");
	defnamelist[i++].name = QString("SOBEK date string");
   defnamelist[i].value = QString("eros.map");
   defnamelist[i++].name = QString("Erosion map");
   defnamelist[i].value = QString("depo.map");
   defnamelist[i++].name = QString("Deposition map");
   defnamelist[i].value = QString("soilloss.map");
   defnamelist[i++].name = QString("Soilloss map");
   defnamelist[i].value = QString("totlandunit.txt");
   defnamelist[i++].name = QString("Filename landunit output");
   defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Simulation times]");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Begin time");
   defnamelist[i].value = QString("100");
   defnamelist[i++].name = QString("End time");
   defnamelist[i].value = QString("0.15");
   defnamelist[i++].name = QString("Timestep");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[General options]");
   defnamelist[i].value = QString("1");
	defnamelist[i++].name = QString("Include Rainfall");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include snowmelt");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("No Erosion simulation");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include main channels");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include tile drains");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include channel infil");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include channel baseflow");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("All water and sediment to outlet");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("No erosion at outlet");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Alternative flow detachment");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Simple depression storage");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Hard Surfaces");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Interception]");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Use canopy storage map");
   defnamelist[i].value = QString("1");
   defnamelist[i++].name = QString("Canopy storage equation");
   defnamelist[i].value = QString("0.05");
   defnamelist[i++].name = QString("Stemflow fraction");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Conservation]");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include grass strips");
	defnamelist[i++].name = QString("Grassstrip Mannings n");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include buffers");
   defnamelist[i++].name = QString("Sediment bulk density");
   defnamelist[i].value = QString("0");
   defnamelist[i++].name = QString("Include Sediment traps");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Calibration]");
   defnamelist[i].value = QString("1.0");
   defnamelist[i++].name = QString("Ksat calibration");
   defnamelist[i].value = QString("1.0");
   defnamelist[i++].name = QString("N calibration");
   defnamelist[i].value = QString("1.0");
   defnamelist[i++].name = QString("Channel Ksat calibration");
   defnamelist[i].value = QString("1.0");
   defnamelist[i++].name = QString("Channel N calibration");
   defnamelist[i].value = QString("0.1");
   defnamelist[i++].name = QString("Splash Delivery Ratio");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Gully options]");
	defnamelist[i++].name = QString("Fcrit relation");
	defnamelist[i++].name = QString("Threshold gradient");
	defnamelist[i++].name = QString("QW relation");
	defnamelist[i++].name = QString("QW param A");
	defnamelist[i++].name = QString("QW param B");
	defnamelist[i++].name = QString("Gully infiltration");
	defnamelist[i++].name = QString("Use initial gully dimensions");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Infiltration]");
	defnamelist[i++].name = QString("Infil Method");
	defnamelist[i++].name = QString("Include wheeltracks");
	defnamelist[i++].name = QString("Include crusts");
	defnamelist[i++].name = QString("Impermeable sublayer");
	defnamelist[i++].name = QString("Subsoil drainage");
	defnamelist[i++].name = QString("Table Directory");
	defnamelist[i++].name = QString("Table File");
	defnamelist[i++].name = QString("SWATRE internal minimum timestep");
	defnamelist[i++].name = QString("Matric head files");
	defnamelist[i++].name = QString("Geometric mean Ksat");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Output maps]");
	defnamelist[i++].name = QString("Runoff maps in l/s/m");
	defnamelist[i++].name = QString("Timeseries as PCRaster");
	defnamelist[i++].name = QString("Timeplot as PCRaster");
	defnamelist[i++].name = QString("Erosion map units (0/1/2)");
	defnamelist[i++].name = QString("Regular runoff output");
	defnamelist[i++].name = QString("Output interval");
	defnamelist[i++].name = QString("User defined output");
	defnamelist[i++].name = QString("Output times");
	defnamelist[i++].name = QString("CheckOutputMaps");
	defnamelist[i++].name = QString("CheckOutputMapsNUT");
	defnamelist[i++].name = QString("CheckOutputMapsMC");
	defnamelist[i++].name = QString("CheckOutputMapsGUL");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Texture classes]");
	defnamelist[i++].name = QString("ClassMu");
	defnamelist[i++].name = QString("[map names]");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputBASIC]");
	defnamelist[i++].name = QString("OUTRUNOFF");
	defnamelist[i++].name = QString("OUTCONC");
	defnamelist[i++].name = QString("OUTWH");
	defnamelist[i++].name = QString("OUTRWH");
	defnamelist[i++].name = QString("OUTTC");
	defnamelist[i++].name = QString("OUTEROS");
	defnamelist[i++].name = QString("OUTDEPO");
	defnamelist[i++].name = QString("OUTVELO");
	defnamelist[i++].name = QString("OUTINF");
	defnamelist[i++].name = QString("OUTSS");
	defnamelist[i++].name = QString("OUTCHVOL");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputMC]");
	defnamelist[i++].name = QString("OUTMU0");
	defnamelist[i++].name = QString("OUTMU1");
	defnamelist[i++].name = QString("OUTMU2");
	defnamelist[i++].name = QString("OUTMU3");
	defnamelist[i++].name = QString("OUTMU4");
	defnamelist[i++].name = QString("OUTMU5");
	defnamelist[i++].name = QString("OUTD50SUSP");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputNut]");
	defnamelist[i++].name = QString("OUTPSOLUT");
	defnamelist[i++].name = QString("OUTPSUS");
	defnamelist[i++].name = QString("OUTPINF");
	defnamelist[i++].name = QString("OUTNH4SOLUT");
	defnamelist[i++].name = QString("OUTNH4SUS");
	defnamelist[i++].name = QString("OUTNH4INF");
	defnamelist[i++].name = QString("OUTNO3SOLUT");
	defnamelist[i++].name = QString("OUTNO3SUS");
	defnamelist[i++].name = QString("OUTNO3INF");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputNutErosDep]");
	defnamelist[i++].name = QString("OUTPDEP");
	defnamelist[i++].name = QString("OUTNH4DEP");
	defnamelist[i++].name = QString("OUTNO3DEP");
	defnamelist[i++].name = QString("OUTPDET");
	defnamelist[i++].name = QString("OUTNH4DET");
	defnamelist[i++].name = QString("OUTNO3DET");
	defnamelist[i++].name = QString("");
   defnamelist[i++].name = QString("[OutputGul]");
	defnamelist[i++].name = QString("OUTGULD");
	defnamelist[i++].name = QString("OUTGULW");
	defnamelist[i++].name = QString("OUTGULA");
	defnamelist[i++].name = QString("OUTGULF");
	defnamelist[i++].name = QString("OUTGULDEM");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Catchment]");
	defnamelist[i++].name = QString("grad");
	defnamelist[i++].name = QString("ldd");
	defnamelist[i++].name = QString("outlet");
	defnamelist[i++].name = QString("ID");
	defnamelist[i++].name = QString("outpoint");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Landuse]");
   defnamelist[i++].name = QString("landunit");
   defnamelist[i++].name = QString("cover");
   defnamelist[i++].name = QString("lai");
	defnamelist[i++].name = QString("ch");
	defnamelist[i++].name = QString("smax");
	defnamelist[i++].name = QString("road");
	defnamelist[i++].name = QString("grasswidth");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Buffers]");
	defnamelist[i++].name = QString("bufferID");
	defnamelist[i++].name = QString("bufferVolume");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Snowmelt]");
	defnamelist[i++].name = QString("SnowID");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Erosion]");
	defnamelist[i++].name = QString("coh");
	defnamelist[i++].name = QString("cohadd");
	defnamelist[i++].name = QString("aggrstab");
	defnamelist[i++].name = QString("d50");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Surface]");
	defnamelist[i++].name = QString("rr");
	defnamelist[i++].name = QString("manning");
	defnamelist[i++].name = QString("crustfrc");
	defnamelist[i++].name = QString("compfrc");
	defnamelist[i++].name = QString("stonefrc");
	defnamelist[i++].name = QString("hardsurf");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilSwatre]");
	defnamelist[i++].name = QString("profmap");
	defnamelist[i++].name = QString("profcrst");
	defnamelist[i++].name = QString("profwltr");
	defnamelist[i++].name = QString("profgras");
	defnamelist[i++].name = QString("inithead");
	defnamelist[i++].name = QString("headout");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilExtra]");
	defnamelist[i++].name = QString("ksatcrst");
	defnamelist[i++].name = QString("ksatcomp");
	defnamelist[i++].name = QString("ksatgras");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilGA1]");
	defnamelist[i++].name = QString("ksat1");
	defnamelist[i++].name = QString("psi1");
	defnamelist[i++].name = QString("thetas1");
	defnamelist[i++].name = QString("thetai1");
	defnamelist[i++].name = QString("soildep1");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilGA2]");
	defnamelist[i++].name = QString("ksat2");
	defnamelist[i++].name = QString("psi2");
	defnamelist[i++].name = QString("thetas2");
	defnamelist[i++].name = QString("thetai2");
	defnamelist[i++].name = QString("soildep2");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Channelinfil]");
	defnamelist[i++].name = QString("chanksat");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Channels]");
	defnamelist[i++].name = QString("lddchan");
	defnamelist[i++].name = QString("chanwidth");
	defnamelist[i++].name = QString("chanside");
	defnamelist[i++].name = QString("changrad");
	defnamelist[i++].name = QString("chanman");
	defnamelist[i++].name = QString("chancoh");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[ChannelBaseflow]");
	defnamelist[i++].name = QString("chanbaseflux");
	defnamelist[i++].name = QString("chanincrease");
	defnamelist[i++].name = QString("chanvolini");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Wheeltrack]");
	defnamelist[i++].name = QString("lddwheel");
	defnamelist[i++].name = QString("wheelnbr");
	defnamelist[i++].name = QString("wheelwidth");
	defnamelist[i++].name = QString("wheeldepth");
	defnamelist[i++].name = QString("wheelgradient");
	defnamelist[i++].name = QString("wheelman");
	defnamelist[i++].name = QString("wheelcohesion");
	defnamelist[i++].name = QString("ksatwt");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Texture]");
	defnamelist[i++].name = QString("fractionmu0");
	defnamelist[i++].name = QString("fractionmu1");
	defnamelist[i++].name = QString("fractionmu2");
	defnamelist[i++].name = QString("fractionmu3");
	defnamelist[i++].name = QString("fractionmu4");
	defnamelist[i++].name = QString("fractionmu5");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsP]");
	defnamelist[i++].name = QString("pcont");
	defnamelist[i++].name = QString("psolute");
	defnamelist[i++].name = QString("pefficiency");
	defnamelist[i++].name = QString("psorp");
	defnamelist[i++].name = QString("pconv");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsNO3]");
	defnamelist[i++].name = QString("no3cont");
	defnamelist[i++].name = QString("no3solute");
	defnamelist[i++].name = QString("no3efficiency");
	defnamelist[i++].name = QString("no3sorp");
	defnamelist[i++].name = QString("no3conv");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsNH4]");
	defnamelist[i++].name = QString("nh4cont");
	defnamelist[i++].name = QString("nh4solute");
	defnamelist[i++].name = QString("nh4efficiency");
	defnamelist[i++].name = QString("nh4sorp");
	defnamelist[i++].name = QString("nh4conv");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsBD]");
	defnamelist[i++].name = QString("bulk");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Gully]");
	defnamelist[i++].name = QString("dem");
	defnamelist[i++].name = QString("gullyn");
	defnamelist[i++].name = QString("bulkdens1");
	defnamelist[i++].name = QString("gulksat1");
	defnamelist[i++].name = QString("gullydep");
	defnamelist[i++].name = QString("gullycoh");
	defnamelist[i++].name = QString("bulkdens2");
	defnamelist[i++].name = QString("gulksat2");
	defnamelist[i++].name = QString("nonfcrit");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[GullyInit]");
	defnamelist[i++].name = QString("gulwinit");
	defnamelist[i++].name = QString("guldinit");
	nrdefnamelist = i;

   // fill with map variables with default mapnames
   for (int j = 0; j < nrdefnamelist; j++)
      for (int k = 0; k < nrmaplist; k++)
      {
         if (mapList[k].id.toUpper() == defnamelist[j].name.toUpper())
         {
            defnamelist[j].value = mapList[k].name;
            //qDebug() << "default" << k << mapList[k].name << mapList[k].id;
         }
      }

}
//---------------------------------------------------------------------------

