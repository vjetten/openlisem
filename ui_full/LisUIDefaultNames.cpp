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
      namelist[i].name.clear();
      namelist[i].value.clear();
   }
   // clear first

   // runfile has structure:
   // name=value
   i = 0;
   namelist[i++].name = QString("[openLISEM runfile version 4]");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[LISEM main type]");
   namelist[i++].name = QString("LISEM Type");
   namelist[i++].name = QString("");
//	namelist[i++].name = QString("[Work Directory]");
//	namelist[i++].name = QString("WorkDir");
//	namelist[i++].name = QString("");
   namelist[i++].name = QString("[Input]");
   namelist[i++].name = QString("Map Directory");
   namelist[i++].name = QString("Include Rainfall");
   namelist[i++].name = QString("Rainfall Directory");
   namelist[i++].name = QString("Rainfall file");
   namelist[i++].name = QString("Include Snowmelt");
   namelist[i++].name = QString("Snowmelt Directory");
   namelist[i++].name = QString("Snowmelt file");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Output]");
   namelist[i++].name = QString("Result Directory");
   namelist[i].value = QString("totals.txt");
   namelist[i++].name = QString("Main results file");
   namelist[i].value = QString("hydrograph.csv");
   namelist[i++].name = QString("Filename point output");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Report point output separate");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Report point output for SOBEK");
   namelist[i++].name = QString("SOBEK date string");
   namelist[i].value = QString("eros.map");
   namelist[i++].name = QString("Erosion map");
   namelist[i].value = QString("depo.map");
   namelist[i++].name = QString("Deposition map");
   namelist[i].value = QString("soilloss.map");
   namelist[i++].name = QString("Soilloss map");
   namelist[i].value = QString("totlandunit.txt");
   namelist[i++].name = QString("Filename landunit output");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Simulation times]");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Begin time");
   namelist[i].value = QString("100");
   namelist[i++].name = QString("End time");
   namelist[i].value = QString("0.15");
   namelist[i++].name = QString("Timestep");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[General options]");
   namelist[i].value = QString("1");
   namelist[i++].name = QString("Include Rainfall");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include snowmelt");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("No Erosion simulation");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include main channels");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include channel infil");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include channel baseflow");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include tile drains");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("All water and sediment to outlet");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("No erosion at outlet");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Alternative flow detachment");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Simple depression storage");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Hard Surfaces");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Interception]");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Use canopy storage map");
   namelist[i].value = QString("1");
   namelist[i++].name = QString("Canopy storage equation");
   namelist[i].value = QString("0.05");
   namelist[i++].name = QString("Stemflow fraction");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Conservation]");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include grass strips");
   namelist[i++].name = QString("Grassstrip Mannings n");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include buffers");
   namelist[i++].name = QString("Sediment bulk density");
   namelist[i].value = QString("0");
   namelist[i++].name = QString("Include Sediment traps");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Calibration]");
   namelist[i].value = QString("1.0");
   namelist[i++].name = QString("Ksat calibration");
   namelist[i].value = QString("1.0");
   namelist[i++].name = QString("N calibration");
   namelist[i].value = QString("1.0");
   namelist[i++].name = QString("Channel Ksat calibration");
   namelist[i].value = QString("1.0");
   namelist[i++].name = QString("Channel N calibration");
   namelist[i].value = QString("0.1");
   namelist[i++].name = QString("Splash Delivery Ratio");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Gully options]");
   namelist[i++].name = QString("Fcrit relation");
   namelist[i++].name = QString("Threshold gradient");
   namelist[i++].name = QString("QW relation");
   namelist[i++].name = QString("QW param A");
   namelist[i++].name = QString("QW param B");
   namelist[i++].name = QString("Gully infiltration");
   namelist[i++].name = QString("Use initial gully dimensions");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Infiltration]");
   namelist[i++].name = QString("Infil Method");
   namelist[i++].name = QString("Include wheeltracks");
   namelist[i++].name = QString("Include crusts");
   namelist[i++].name = QString("Impermeable sublayer");
   namelist[i++].name = QString("Subsoil drainage");
   namelist[i++].name = QString("Table Directory");
   namelist[i++].name = QString("Table File");
   namelist[i++].name = QString("SWATRE internal minimum timestep");
   namelist[i++].name = QString("Matric head files");
   namelist[i++].name = QString("Geometric mean Ksat");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Output maps]");
   namelist[i++].name = QString("Runoff maps in l/s/m");
   namelist[i++].name = QString("Timeseries as PCRaster");
   namelist[i++].name = QString("Timeplot as PCRaster");
   namelist[i++].name = QString("Erosion map units (0/1/2)");
   namelist[i++].name = QString("Regular runoff output");
   namelist[i++].name = QString("Output interval");
   namelist[i++].name = QString("User defined output");
   namelist[i++].name = QString("Output times");
   namelist[i++].name = QString("CheckOutputMaps");
   namelist[i++].name = QString("CheckOutputMapsNUT");
   namelist[i++].name = QString("CheckOutputMapsMC");
   namelist[i++].name = QString("CheckOutputMapsGUL");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Texture classes]");
   namelist[i++].name = QString("ClassMu");
   namelist[i++].name = QString("[map names]");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[OutputBASIC]");
   namelist[i++].name = QString("OUTRUNOFF");
   namelist[i++].name = QString("OUTCONC");
   namelist[i++].name = QString("OUTWH");
   namelist[i++].name = QString("OUTRWH");
   namelist[i++].name = QString("OUTTC");
   namelist[i++].name = QString("OUTEROS");
   namelist[i++].name = QString("OUTDEPO");
   namelist[i++].name = QString("OUTVELO");
   namelist[i++].name = QString("OUTINF");
   namelist[i++].name = QString("OUTSS");
   namelist[i++].name = QString("OUTCHVOL");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[OutputMC]");
   namelist[i++].name = QString("OUTMU0");
   namelist[i++].name = QString("OUTMU1");
   namelist[i++].name = QString("OUTMU2");
   namelist[i++].name = QString("OUTMU3");
   namelist[i++].name = QString("OUTMU4");
   namelist[i++].name = QString("OUTMU5");
   namelist[i++].name = QString("OUTD50SUSP");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[OutputNut]");
   namelist[i++].name = QString("OUTPSOLUT");
   namelist[i++].name = QString("OUTPSUS");
   namelist[i++].name = QString("OUTPINF");
   namelist[i++].name = QString("OUTNH4SOLUT");
   namelist[i++].name = QString("OUTNH4SUS");
   namelist[i++].name = QString("OUTNH4INF");
   namelist[i++].name = QString("OUTNO3SOLUT");
   namelist[i++].name = QString("OUTNO3SUS");
   namelist[i++].name = QString("OUTNO3INF");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[OutputNutErosDep]");
   namelist[i++].name = QString("OUTPDEP");
   namelist[i++].name = QString("OUTNH4DEP");
   namelist[i++].name = QString("OUTNO3DEP");
   namelist[i++].name = QString("OUTPDET");
   namelist[i++].name = QString("OUTNH4DET");
   namelist[i++].name = QString("OUTNO3DET");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[OutputGul]");
   namelist[i++].name = QString("OUTGULD");
   namelist[i++].name = QString("OUTGULW");
   namelist[i++].name = QString("OUTGULA");
   namelist[i++].name = QString("OUTGULF");
   namelist[i++].name = QString("OUTGULDEM");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Catchment]");
   namelist[i++].name = QString("grad");
   namelist[i++].name = QString("ldd");
   namelist[i++].name = QString("outlet");
   namelist[i++].name = QString("ID");
   namelist[i++].name = QString("outpoint");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Landuse]");
   namelist[i++].name = QString("landunit");
   namelist[i++].name = QString("cover");
   namelist[i++].name = QString("lai");
   namelist[i++].name = QString("ch");
   namelist[i++].name = QString("smax");
   namelist[i++].name = QString("road");
   namelist[i++].name = QString("grasswidth");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Buffers]");
   namelist[i++].name = QString("bufferID");
   namelist[i++].name = QString("bufferVolume");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Snowmelt]");
   namelist[i++].name = QString("SnowID");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Erosion]");
   namelist[i++].name = QString("coh");
   namelist[i++].name = QString("cohadd");
   namelist[i++].name = QString("aggrstab");
   namelist[i++].name = QString("d50");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Surface]");
   namelist[i++].name = QString("rr");
   namelist[i++].name = QString("manning");
   namelist[i++].name = QString("crustfrc");
   namelist[i++].name = QString("compfrc");
   namelist[i++].name = QString("stonefrc");
   namelist[i++].name = QString("hardsurf");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[InfilSwatre]");
   namelist[i++].name = QString("profmap");
   namelist[i++].name = QString("profcrst");
   namelist[i++].name = QString("profwltr");
   namelist[i++].name = QString("profgras");
   namelist[i++].name = QString("inithead");
   namelist[i++].name = QString("headout");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[InfilExtra]");
   namelist[i++].name = QString("ksatcrst");
   namelist[i++].name = QString("ksatcomp");
   namelist[i++].name = QString("ksatgras");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[InfilGA1]");
   namelist[i++].name = QString("ksat1");
   namelist[i++].name = QString("psi1");
   namelist[i++].name = QString("thetas1");
   namelist[i++].name = QString("thetai1");
   namelist[i++].name = QString("soildep1");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[InfilGA2]");
   namelist[i++].name = QString("ksat2");
   namelist[i++].name = QString("psi2");
   namelist[i++].name = QString("thetas2");
   namelist[i++].name = QString("thetai2");
   namelist[i++].name = QString("soildep2");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Channelinfil]");
   namelist[i++].name = QString("chanksat");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Channels]");
   namelist[i++].name = QString("lddchan");
   namelist[i++].name = QString("chanwidth");
   namelist[i++].name = QString("chanside");
   namelist[i++].name = QString("changrad");
   namelist[i++].name = QString("chanman");
   namelist[i++].name = QString("chancoh");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[ChannelBaseflow]");
   namelist[i++].name = QString("chanbaseflux");
   namelist[i++].name = QString("chanincrease");
   namelist[i++].name = QString("chanvolini");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Wheeltrack]");
   namelist[i++].name = QString("lddwheel");
   namelist[i++].name = QString("wheelnbr");
   namelist[i++].name = QString("wheelwidth");
   namelist[i++].name = QString("wheeldepth");
   namelist[i++].name = QString("wheelgradient");
   namelist[i++].name = QString("wheelman");
   namelist[i++].name = QString("wheelcohesion");
   namelist[i++].name = QString("ksatwt");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Texture]");
   namelist[i++].name = QString("fractionmu0");
   namelist[i++].name = QString("fractionmu1");
   namelist[i++].name = QString("fractionmu2");
   namelist[i++].name = QString("fractionmu3");
   namelist[i++].name = QString("fractionmu4");
   namelist[i++].name = QString("fractionmu5");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[NutsP]");
   namelist[i++].name = QString("pcont");
   namelist[i++].name = QString("psolute");
   namelist[i++].name = QString("pefficiency");
   namelist[i++].name = QString("psorp");
   namelist[i++].name = QString("pconv");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[NutsNO3]");
   namelist[i++].name = QString("no3cont");
   namelist[i++].name = QString("no3solute");
   namelist[i++].name = QString("no3efficiency");
   namelist[i++].name = QString("no3sorp");
   namelist[i++].name = QString("no3conv");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[NutsNH4]");
   namelist[i++].name = QString("nh4cont");
   namelist[i++].name = QString("nh4solute");
   namelist[i++].name = QString("nh4efficiency");
   namelist[i++].name = QString("nh4sorp");
   namelist[i++].name = QString("nh4conv");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[NutsBD]");
   namelist[i++].name = QString("bulk");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[Gully]");
   namelist[i++].name = QString("dem");
   namelist[i++].name = QString("gullyn");
   namelist[i++].name = QString("bulkdens1");
   namelist[i++].name = QString("gulksat1");
   namelist[i++].name = QString("gullydep");
   namelist[i++].name = QString("gullycoh");
   namelist[i++].name = QString("bulkdens2");
   namelist[i++].name = QString("gulksat2");
   namelist[i++].name = QString("nonfcrit");
   namelist[i++].name = QString("");
   namelist[i++].name = QString("[GullyInit]");
   namelist[i++].name = QString("gulwinit");
   namelist[i++].name = QString("guldinit");
   nrnamelist = i;

   // fill with map variables with default mapnames
   for (int j = 0; j < nrnamelist; j++)
      for (int k = 0; k < nrmaplist; k++)
      {
         if (mapList[k].name.toUpper() == namelist[j].name.toUpper())
         {
            namelist[j].value = mapList[k].value;
            //qDebug() << "default" << k << mapList[k].name << mapList[k].id;
         }
      }

}
//---------------------------------------------------------------------------

