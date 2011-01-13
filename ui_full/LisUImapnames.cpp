/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/
/*
 * fill the map list structure with names and basic descriptions
 */


#include "lisemqt.h"

//--------------------------------------------------------------------
// DEFmaps has default var names, filenames and descriptions.
// first number 0/1/2 is flag if title treenode or subnode
// fill the mapList structure
void lisemqt::fillMapnames()
{
   int subbranch = 0, branch = -1, nr = 0, dec = 0;
   QStringList SL;

   DefaultMapnames();
   //VJ 110113 take the default map list

   for (int i = 0; i < DEFmaps.count(); i++)\
   {

      if (DEFmaps[i].startsWith("0"))
      {
         branch++;
         subbranch = 0;
         dec = 0;
      }
      if (DEFmaps[i].startsWith("1"))
      {
         dec += 10;
         subbranch = 0;
      }
      if (DEFmaps[i].startsWith("2"))
      {
         nr++;
         subbranch++;

         SL = DEFmaps[i].split(";",QString::SkipEmptyParts);
         mapList[nr].groupnr=branch;
         mapList[nr].varnr=subbranch+dec;
         mapList[nr].value=SL[4];
         mapList[nr].name=SL[3];
         mapList[nr].dir="";
      }
   }
   nrmaplist = nr;
}
//--------------------------------------------------------------------
// enables or disables a branch and expands or contracts it
void lisemqt::change_MapNameModel(int parentrow, int selrow, bool setit)
{
	if (MapNameModel)
	{
		QModelIndex indexParent = MapNameModel->index(parentrow, 0);
		// select and expand chosen main and sub-branch
		if (selrow >= 0)
		{
			MapNameModel->setFlag(setit, parentrow);
			// set top level node
			if (selrow == 0)
				for (int k = 0; k < MapNameModel->rowCount(indexParent); k++)
					MapNameModel->setFlag(setit, k, indexParent);
			else
				if (selrow >= 10)
				{
					selrow -= 10;
					MapNameModel->setFlag(setit, selrow, indexParent);

					QModelIndex indexChild = MapNameModel->index(selrow, 0, indexParent);

					for (int k = 0; k < MapNameModel->rowCount(indexChild); k++)
						MapNameModel->setFlag(setit, k, indexChild);
				}

			if (setit)
			{
				treeView->expand(MapNameModel->index(parentrow,0));
				if (selrow >= 0)
					treeView->expand(MapNameModel->index(selrow, 0, indexParent));
			}
			else
			{
				if (selrow >= 0)
					treeView->collapse(MapNameModel->index(selrow, 0, indexParent));
				if (selrow == 0)
					treeView->collapse(MapNameModel->index(parentrow,0));
			}
		}
	}
}
//--------------------------------------------------------------------
// this is called twice: first to initialize the interface
//also after each call of a runfile so that the runfile mapnames are loaded
void lisemqt::FillMapTree()
{
	if (MapNameModel)
	{
		delete MapNameModel;
		MapNameModel = NULL;
	}

   MapNameModel = new TreeModel(DEFmaps);

	treeView->setModel(MapNameModel);

	treeView->setColumnWidth(0,196);
	treeView->setColumnWidth(1,196);
	treeView->setColumnWidth(2,400);
	treeView->setColumnWidth(3,0);
	treeView->setColumnWidth(4,0);
	treeView->setAlternatingRowColors(true);

	change_MapNameModel(RAINFALLMAPS,0, true);
	change_MapNameModel(CATCHMENTMAPS,0, true);
	change_MapNameModel(LANDUSEMAPS,0, true);
	change_MapNameModel(SURFACEMAPS,0, true);
	change_MapNameModel(EROSIONMAPS,0, true);
	// enable basic maps, tree nodes 0-4 = first 5 branches
   
	treeView->collapseAll();
}
//--------------------------------------------------------------------
//VJ 101222 this function is not used, replaced by the next function
// get user edited name from tre and put in listMaps structure for saving
void lisemqt::doChangeMapname(QModelIndex topLeft, QModelIndex bottomRight )
{
   int groupnr = topLeft.parent().parent().row(); // assume 3 level structure
   int varnr = topLeft.row();

   // if not 3 level structure then use 2nd level
   if (topLeft.parent().parent().row() < 0)
      groupnr = topLeft.parent().row();

   if (groupnr == INFILTRATIONMAPS || groupnr == CHANNELSMAPS || groupnr == NUTRIENTSMAPS)
      varnr = (topLeft.parent().row()+1)*10 + topLeft.row();

   for (int k = 0; k < nrmaplist; k++)
   {
      if (mapList[k].groupnr == groupnr && mapList[k].varnr == varnr)
      {
         QVariant d = MapNameModel->data(topLeft,Qt::DisplayRole);;//MapNameModel->data(MapNameModel->index(j, k, indexParent),0);
         //qDebug() << d;
         mapList[k].name = d.toString();
      }
   }
}
//--------------------------------------------------------------------
// double click on treeView mapname opens this function to get the map name
// it looks for the correct map in mapList and copies the file name into mapList
void lisemqt::doOpenMapname(QModelIndex topLeft)
{
   if (topLeft.column() != 1 || (topLeft.row() == 0 && topLeft.parent().row() < 0))
      return;
   // if not a mapname return
   // else find out which mapname it is

   int groupnr = topLeft.parent().parent().row();
   // assume 3 level structure: e.g. infil->G%A1-> map x
   int k, varnr = topLeft.row();

   // if not 3 level structure then use 2nd level
   if (topLeft.parent().parent().row() < 0)
      groupnr = topLeft.parent().row();

   if (groupnr == INFILTRATIONMAPS || groupnr == CHANNELSMAPS || groupnr == NUTRIENTSMAPS)
      varnr = (topLeft.parent().row()+1)*10 + topLeft.row();
   // correct for 3 level structures

   // find the map in mapList
   for (k = 0; k < nrmaplist; k++)
   {
      if (mapList[k].groupnr == groupnr && mapList[k].varnr == varnr)
      {
         QVariant d = MapNameModel->data(topLeft,Qt::DisplayRole);
         break;
      }
   }

   QString path = QFileDialog::getOpenFileName(this,	QString("Select the map: %1;").arg(mapList[k].name),E_MapDir->text(),QString("*.map *.csf;;*.*"));
   // open file dialog


   if (!path.isEmpty())// && QFileInfo(path).exists())
   {
      MAP *m = Mopen(QFileInfo(path).absoluteFilePath().toAscii().constData(), M_READ);
      if (m == NULL)
      {
         QMessageBox::critical(this, "openLISEM",
                              QString("File \"%1\" is not a PCRaster map.")
                              .arg(path));
         return;
      }

      mapList[k].value = QFileInfo(path).fileName();
      mapList[k].dir = QFileInfo(path).dir().path();
      // put the name and path into he mapList structure
      //qDebug() << "mapname edit" <<  mapList[k].name << mapList[k].id << k;
      QVariant d(mapList[k].value);
      MapNameModel->setData(topLeft, d, Qt::EditRole);
      // put the name into the treeView and MapNameModel
   }
}
//--------------------------------------------------------------------
void lisemqt::DefaultMapnames()
{

//# interface maplist, DO NOT CHANGE if you don't know what you are doing
//# syntax: branch level; keyword; default mapname; description; variable name
   DEFmaps.append("0;Rainfall");
   DEFmaps.append("2;ID;ID.map;Raingauge zone ID numbers;ID");
   DEFmaps.append("0;Catchment");
   DEFmaps.append("2;Gradient;grad.map;Sine of slope gradient in direction of flow;grad");
   DEFmaps.append("2;LDD;ldd.map;Local surface Drainage Direction network;ldd");
   DEFmaps.append("2;Outlet;outlet.map;Main catchment outlet corresponding to LDD map;outlet");
   DEFmaps.append("2;Points;outpoint.map;Reporting points for hydrograph/sedigraph (1 to nr);outpoint");
   DEFmaps.append("0;Landuse");
   DEFmaps.append("2;Units;landunit.map;Classified land unit map (integers 0-n) for output of erosion values;landunit");
   DEFmaps.append("2;Cover;per.map;Fraction surface cover by vegetation and residue;cover");
   DEFmaps.append("2;LAI;lai.map;Leaf area index of the plant cover in a gridcell (m2/m2);lai");
   DEFmaps.append("2;Height;ch.map;Plant height (m);ch");
   DEFmaps.append("2;Road width;roadwidt.map;Width of impermeable roads (m);road");
   DEFmaps.append("2;Grass strips;grasswid.map;Width of grass strips (m);grasswidth");
   DEFmaps.append("2;Canopy storage;smax.map;Maximum canopy storage (mm);smax");
   DEFmaps.append("0;Surface");
   DEFmaps.append("2;RR;rr.map;Random Roughness (here standard deviation of heights) (cm);rr");
   DEFmaps.append("2;n;n.map;Manning's n (-);manning");
   DEFmaps.append("2;Stoniness;stonefrc.map;Fraction covered by stones (affects only splash det.) (-);stonefrc");
   DEFmaps.append("2;Crust;crustfrc.map;Fraction of gridcell covered with Crust (-) (see also ksat crust);crustfrc");
   DEFmaps.append("2;Compacted;compfrc.map;Fraction of gridcell compacted (e.g. wheeltracks)(-) (see also ksat compacted);compfrc");
   DEFmaps.append("2;Hard Surface;hardsurf.map;No interception/infiltration/detachment (value 1);hardsurf");
   DEFmaps.append("0;Erosion");
   DEFmaps.append("2;Cohesion;coh.map;Cohesion (kPa);coh");
   DEFmaps.append("2;Cohesion;cohadd.map;Extra cohesion factor by e.g. plant roots (kPa);cohadd");
   DEFmaps.append("2;Aggregates;aggrstab.map;Aggregate stability for splash erosion (-);aggrstab");
   DEFmaps.append("2;D50;d50.map;Median of the texture of the suspendeed matter (mu);d50");
   DEFmaps.append("0;Infiltration");
   DEFmaps.append("1;Swatre");
   DEFmaps.append("2;Profile soil;profile.map;ID numbers corresponding to land units in profile table;profmap");
   DEFmaps.append("2;Prof. Crust;profcrst.map;ID numbers of crusted soils (using also profile table);profcrst");
   DEFmaps.append("2;Prof. Wheel;profwltr.map;ID numbers of compacted wheel tracks (using also profile table);profwltr");
   DEFmaps.append("2;Prof. Grass;profgras.map;ID numbers of grasstrips (using also profile table);profgras");
   DEFmaps.append("2;Initial suction;inithead;initial matrix potential (cm) of layers 001 to nnn (filename witout extension);inithead");
   DEFmaps.append("1;1st layer Green&Ampt/Smith&Parlange");
   DEFmaps.append("2;Ksat1;ksat1.map;Layer 1: Saturated Hydraulic Conductivity (mm/h);ksat1");
   DEFmaps.append("2;Psi1;psi1.map;Layer 1: Average suction at the wetting front (cm);psi1");
   DEFmaps.append("2;Thetas1;thetas1.map;Layer 1: Porosity (-);thetas1");
   DEFmaps.append("2;Thetai1;thetai1.map;Layer 1: Initial moisture content (-);thetai1");
   DEFmaps.append("2;Depth1;soildep1.map;Layer 1: Depth (mm) to bottom of layer 1;soildep1");
   DEFmaps.append("1;2nd layer Green&Ampt/Smith&Parlange");
   DEFmaps.append("2;Ksat2;ksat2.map;Layer 2: Saturated Hydraulic Conductivity (mm/h);ksat2");
   DEFmaps.append("2;Psi2;psi2.map;Layer 2: Average suction at the wetting front (cm);psi2");
   DEFmaps.append("2;Thetas2;thetas2.map;Layer 2: Porosity (-);thetas2");
   DEFmaps.append("2;Thetai2;thetai2.map;Layer 2: Initial moisture content (-);thetai2");
   DEFmaps.append("2;Depth2;soildep2.map;Layer 2: Depth (mm) to bottom of layer 2;soildep2");
   DEFmaps.append("1;Ksat subtraction");
   DEFmaps.append("2;Ksat1;ksat1.map;Saturated Hydraulic Conductivity (mm/h);ksat1");
   DEFmaps.append("1;Special surfaces");
   DEFmaps.append("2;Ksat Crust;ksatcrst.map;Ksat of crusts (all models except SWATRE) (mm/h);ksatcrst");
   DEFmaps.append("2;Ksat Compact;ksatcomp.map;Ksat of compacted areas (all models except SWATRE) (mm/h);ksatcomp");
   DEFmaps.append("2;Ksat Grass;ksatgras.map;Ksat of grassstrips (all models except SWATRE) (mm/h);ksatgras");
   DEFmaps.append("0;Channels");
   DEFmaps.append("1;Channel properties");
   DEFmaps.append("2;LDD;lddchan.map;LDD of main channel (must be 1 branch connected to the outlet);lddchan");
   DEFmaps.append("2;Width;chanwidt.map;Channel width (m);chanwidth");
   DEFmaps.append("2;Side angle;chanside.map;Channel side angle (tan angle  channel side and surface: 0 is rectangular);chanside");
   DEFmaps.append("2;Gradient;changrad.map;Slope gradient of channel bed (-);changrad");
   DEFmaps.append("2;N;chanman.map;Mannings n of channel bed (-);chanman");
   DEFmaps.append("2;Cohesion;chancoh.map;Cohesion of channel bed (kPa);chancoh");
   DEFmaps.append("1;Channelinfil");
   DEFmaps.append("2;Ksat;chanksat.map;Infiltration rate of channel bed (mm/h);chanksat");
   DEFmaps.append("1;ChannelBaseflow");
   DEFmaps.append("2;Inflow flux;chanbaseflux.map;Incoming flux into channel from the two sides (m3/s);chanbaseflux");
   DEFmaps.append("2;Increase in baseflow;chanincrease.map;Increase in basevolume during rainstorm (-);chanincrease");
   DEFmaps.append("2;Initial volume;chanvini.map;Initial baseflow water volume in channel (m3);chanvolini");
   DEFmaps.append("0;Buffers");
   DEFmaps.append("2;Buffer ID nr;bufferid.map;ID number for each buffer starting with 1 (0 is outside area);bufferID");
   DEFmaps.append("2;Buffer volume;buffervol.map;Buffer volumes at the locations of the buffers (m3);bufferVolume");
   DEFmaps.append("0;Snowmelt");
   DEFmaps.append("2;Snowmelt ID;snowid.map;Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area);SnowID");
   DEFmaps.append("0;Tile drains");
   DEFmaps.append("2;LDD;lddtile.map;LDD of tile drain system (must be 1 system connected to the outlet);lddtile");
   DEFmaps.append("2;Width;tilewidth.map;Tile drain pipe width (m);tilewidth");
   DEFmaps.append("2;Height;tileheight.map;Tile drain height for max volume (m);tileheight");
   DEFmaps.append("2;Gradient;tilegrad.map;Slope gradient of tile drain bed (-);tilegrad");
   DEFmaps.append("2;N;tileman.map;Mannings n of tile drain bed (-);tileman");
   DEFmaps.append("0;Wheeltracks");
   DEFmaps.append("2;LDD;lddwheel.map;LDD of wheeltrack network (can be separate branches with pits);lddwheel");
   DEFmaps.append("2;Number;wheelnbr.map;Number of wheeltrack channels in a gridcell (-);wheelnbr");
   DEFmaps.append("2;Width;wheelwid.map;Sum of widths of wheeltracks in a gridcell (m);wheelwidth");
   DEFmaps.append("2;Depth;wheeldep.map;Wheel track overflow depth (cm);wheeldepth");
   DEFmaps.append("2;Gradient;wheelgrd.map;DEFmapsope gradient of wheel tracks (-);wheelgradient");
   DEFmaps.append("2;N;wheelman.map;Mannings n of Wheel tracks (-);wheelman");
   DEFmaps.append("2;Cohesion;wheelcoh.map;Cohesion of wheel tracks (kPa);wheelcohesion");
   DEFmaps.append("2;Ksat;ksatwt.map;Saturated hydraulic conductivity of wheel tracks (mm/h);ksatwt");
   DEFmaps.append("0;Texture classes");
   DEFmaps.append("2;Class 0;mu0.map;Clay fraction (MUST BE CLAY < 2mu);fractionmu0");
   DEFmaps.append("2;Class 1;mu1.map;Soil texture fraction for class 1 (-);fractionmu1");
   DEFmaps.append("2;Class 2;mu2.map;Soil texture fraction for class 2 (-);fractionmu2");
   DEFmaps.append("2;Class 3;mu3.map;Soil texture fraction for class 3 (-);fractionmu3");
   DEFmaps.append("2;Class 4;mu4.map;Soil texture fraction for class 4 (-);fractionmu4");
   DEFmaps.append("2;Class 5;mu5.map;Soil texture fraction for class 5 (-);fractionmu5");
   DEFmaps.append("0;Nutrients");
   DEFmaps.append("1;Phosphorus");
   DEFmaps.append("2;Bulk Dens.;bulkdens.map;Bulk density of the topsoil (kg/m3);bulk");
   DEFmaps.append("2;P Content;pcont.map;Phosphate (P) content of the soil (kg/kg);pcont");
   DEFmaps.append("2;P Solute;psolut.map;Initial solute store P in surface layer (kg/m2);psolute");
   DEFmaps.append("2;P Efficiency;peff.map;Extraction efficiency (s-1);pefficiency");
   DEFmaps.append("2;P Sorption;Psorp.map;Sorption isotherm kd (m3/kg);psorp");
   DEFmaps.append("2;P Conversion;Pconv.map;Conversion P from soil content to clay content(-);pconv");
   DEFmaps.append("1;NH4");
   DEFmaps.append("2;NH4 Content;nh4cont.map;Ammonium (NH4+) content of the soil (kg/kg);nh4cont");
   DEFmaps.append("2;NH4 Solute;nh4solut.map;Initial solute store NH4 in surface layer (kg/m2);nh4solute");
   DEFmaps.append("2;NH4 Efficiency;nh4eff.map;Extraction efficiency (s-1);nh4efficiency");
   DEFmaps.append("2;NH4 Sorption;NH4sorp.map;Sorption isotherm kd (m3/kg);nh4sorp");
   DEFmaps.append("2;NH4 Conversion;NH4conv.map;Conversion NH4 from soil content to clay content(-);nh4conv");
   DEFmaps.append("1;NO3");
   DEFmaps.append("2;NO3 Content;NO3cont.map;Nitrate (NO3-) content of the soil (kg/kg);no3cont");
   DEFmaps.append("2;NO3 Solute;NO3solut.map;Initial solute store NO3 in surface layer (kg/m2);no3solute");
   DEFmaps.append("2;NO3 Efficiency;NO3eff.map;Extraction efficiency (s-1);no3efficiency");
   DEFmaps.append("2;NO3 Sorption;NO3sorp.map;Sorption isotherm kd (m3/kg);no3sorp");
   DEFmaps.append("2;NO3 Conversion;NO3conv.map;Conversion NO3 from soil content to clay content(-);no3conv");
   DEFmaps.append("0;Gullies");
   DEFmaps.append("1;General");
   DEFmaps.append("2;DEM;dem.map;Digital elevation model (m);dem");
   DEFmaps.append("2;mannings N;gullyman.map;manning's n gully bottom (-);gullyn");
   DEFmaps.append("2;Excluded areas;noncrit.map;areas to be excluded (1) and rest (0);nonfcrit");
   DEFmaps.append("2;Gully initial Width;gulwinit.map; initial gully width (m);gulwinit");
   DEFmaps.append("2;Gully initial Depth;guldinit.map; initial gully depth (m);guldinit");
   DEFmaps.append("1;Soil Layer 1");
   DEFmaps.append("2;BulkDensity;bulkdens.map;Bulkdensity of topsoil (kg/m3);bulkdens1");
   DEFmaps.append("2;Ksat;ksat1.map;Ksat of topsoil for gully infil (mm/h);gulksat1");
   DEFmaps.append("1;Soil Layer 2");
   DEFmaps.append("2;Depth layer 2;soilDep2.map;Depth to subsoil (cm);gullydep");
   DEFmaps.append("2;Cohesion layer 2;coh2.map;Cohesion of subsoil (kPa);gullycoh");
   DEFmaps.append("2;BulkDensity 2;bulkden2.map;Bulkdensity of subsoil (kg/m3);bulkdens2");
   DEFmaps.append("2;Ksat 2;gulksat2.map;Ksat of subsoil for gully infil (mm/h);gulksat2");

}
