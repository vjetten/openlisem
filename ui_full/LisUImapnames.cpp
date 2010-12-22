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
#define putmap(k,j,i) SL=DEFmaps[DEFmaps.count()-1].split(";",QString::SkipEmptyParts);\
   mapList[k].groupnr=j;mapList[k].varnr=i;mapList[k].name=SL[4];mapList[k].dir="";mapList[k].id=SL[6]
// this define fills a list of mapnames to store user changes


#include "lisemqt.h"

//--------------------------------------------------------------------
// stringlist with default var names, filenames and descriptions.
// first number 0/1/2 is flag if title treenode or subnode
// putmap is a define that fills the mapList structure
void lisemqt::DefaultMapnames()
{
   int i, j, k;
   QStringList SL;
   // example: SL[0]=tree node level, SL[1] and SL[2] are variable numbers, SL[3] is first visible field
   //SL[4] = map filename, SL[5] = description, SL[6] is internal variable name, SL[7] = pathname, empty here

   j = 0;i = 0;k = 0; // j and i are needed to find the data back in the tree structure, k is mapcount
   DEFmaps.append("0;Rainfall");
   DEFmaps.append(QString("2;%1;%2;ID;ID.map;Raingauge zone ID numbers;id;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Catchment");
   DEFmaps.append(QString("2;%1;%2;Gradient;grad.map;Sine of slope gradient in direction of flow;grad;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;LDD;ldd.map;Local surface Drainage Direction network;ldd;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Outlet;outlet.map;Main catchment outlet corresponding to LDD map;outlet;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Points;outpoint.map;Reporting points for hydrograph/sedigraph (1 to nr);outpoint;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Landuse");
   DEFmaps.append(QString("2;%1;%2;Cover;per.map;Fraction surface cover by vegetation and residue;cover;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;LAI;lai.map;Leaf area index of the plant cover in a gridcell (m2/m2);lai;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Height;ch.map;Plant height (m);ch;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Road width;roadwidt.map;Width of impermeable roads (m);road;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Grass strips;grasswid.map;Width of grass strips (m);grasswidth;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Surface");
   DEFmaps.append(QString("2;%1;%2;RR;rr.map;Random Roughness (here standard deviation of heights) (cm);rr;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;n;n.map;Manning's n (-);manning;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Stoniness;stonefrc.map;Fraction covered by stones (affects only splash det.) (-);stonefrc;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Crust;crustfrc.map;Fraction of gridcell covered with Crust (-) (see also ksat crust);crustfrc;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Compacted;compfrc.map;Fraction of gridcell compacted (e.g. wheeltracks)(-) (see also ksat compacted);compfrc;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Hard Surface;hardsurf.map;No interception/infiltration/detachment (value 1);hardsurf;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Erosion");
   DEFmaps.append(QString("2;%1;%2;Cohesion;coh.map;Cohesion (kPa);coh;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Cohesion;cohadd.map;Extra cohesion factor by e.g. plant roots (kPa);cohadd;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Aggregates;aggrstab.map;Aggregate stability for splash erosion (-);aggrstab;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;D50;d50.map;Median of the texture of the suspendeed matter (mu);d50;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 10;
   DEFmaps.append("0;Infiltration");

   DEFmaps.append("1;Swatre");
   DEFmaps.append(QString("2;%1;%2;Profile soil;profile.map;ID numbers corresponding to land units in profile table;profmap;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Prof. Crust;profcrst.map;ID numbers of crusted soils (using also profile table);profcrst;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Prof. Wheel;profwltr.map;ID numbers of compacted wheel tracks (using also profile table);profwltr;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Prof. Grass;profgras.map;ID numbers of grasstrips (using also profile table);profgras;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Initial suction;inithead;initial matrix potential (cm) of layers 001 to nnn (filename witout extension);inithead;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   //DEFmaps.append(QString("2;%1;%2;Output head;headout.map;Locations to write tables of the matrix potential;headout;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 20;
   DEFmaps.append("1;1st layer Green&Ampt/Smith&Parlange");
   DEFmaps.append(QString("2;%1;%2;Ksat1;ksat1.map;Layer 1: Saturated Hydraulic Conductivity (mm/h);ksat1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Psi1;psi1.map;Layer 1: Average suction at the wetting front (cm);psi1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Thetas1;thetas1.map;Layer 1: Porosity (-);thetas1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Thetai1;thetai1.map;Layer 1: Initial moisture content (-);thetai1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Depth1;soildep1.map;Layer 1: Depth (mm) to bottom of layer 1;soildep1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 30;
   DEFmaps.append("1;2nd layer Green&Ampt/Smith&Parlange");
   DEFmaps.append(QString("2;%1;%2;Ksat2;ksat2.map;Layer 2: Saturated Hydraulic Conductivity (mm/h);ksat2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Psi2;psi2.map;Layer 2: Average suction at the wetting front (cm);psi2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Thetas2;thetas2.map;Layer 2: Porosity (-);thetas2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Thetai2;thetai2.map;Layer 2: Initial moisture content (-);thetai2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Depth2;soildep2.map;Layer 2: Depth (mm) to bottom of layer 2;soildep2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 40;
   DEFmaps.append("1;Ksat subtraction");
   DEFmaps.append(QString("2;%1;%2;Ksat1;ksat1.map;Saturated Hydraulic Conductivity (mm/h);ksat1);").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 50;
   DEFmaps.append("1;Special surfaces");
   DEFmaps.append(QString("2;%1;%2;Ksat Crust;ksatcrst.map;Ksat of crusts (all models except SWATRE) (mm/h);ksatcrst;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Ksat Compact;ksatcomp.map;Ksat of compacted areas (all models except SWATRE) (mm/h);ksatcomp;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Ksat Grass;ksatgras.map;Ksat of grassstrips (all models except SWATRE) (mm/h);ksatgras;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Channels");
   i = 10;
   DEFmaps.append("1;Channel properties");
   DEFmaps.append(QString("2;%1;%2;LDD;lddchan.map;LDD of main channel (must be 1 branch connected to the outlet);lddchan;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Width;chanwidt.map;Channel width (m);chanwidth;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Side angle;chanside.map;Channel side angle (tan angle  channel side and surface: 0 is rectangular);chanside;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Gradient;changrad.map;Slope gradient of channel bed (-);changrad;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;N;chanman.map;Mannings n of channel bed (-);chanman;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Cohesion;chancoh.map;Cohesion of channel bed (kPa);chancoh;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 20;
   DEFmaps.append("1;Channelinfil");
   DEFmaps.append(QString("2;%1;%2;Ksat;chanksat.map;Infiltration rate of channel bed (mm/h);chanksat;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 30;
   DEFmaps.append("1;ChannelBaseflow");
   DEFmaps.append(QString("2;%1;%2;Inflow flux;chanbaseflux.map;Incoming flux into channel from the two sides (m3/s);chanbaseflux;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Increase in baseflow;chanincrease.map;Increase in basevolume during rainstorm (-);chanincrease;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Initial volume;chanvini.map;Initial baseflow water volume in channel (m3);chanvolini;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Buffers");
   DEFmaps.append(QString("2;%1;%2;Buffer ID nr;bufferid.map;ID number for each buffer starting with 1 (0 is outside area);bufferID;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Buffer volume;buffervol.map;Buffer volumes at the locations of the buffers (m3);bufferVolume;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   //DEFmaps.append(QString("2;%1;%2;Buffer area;bufferarea.map;Buffer area at locations of the buffers (m2);bufferarea;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Snowmelt");
   DEFmaps.append(QString("2;%1;%2;Snowmelt ID;snowid.map;Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area);SnowID;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Wheeltracks");
   DEFmaps.append(QString("2;%1;%2;LDD;lddwheel.map;LDD of wheeltrack network (can be separate branches with pits);lddwheel;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Number;wheelnbr.map;Number of wheeltrack channels in a gridcell (-);wheelnbr;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Width;wheelwid.map;Sum of widths of wheeltracks in a gridcell (m);wheelwidth;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Depth;wheeldep.map;Wheel track overflow depth (cm);wheeldepth;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Gradient;wheelgrd.map;DEFmapsope gradient of wheel tracks (-);wheelgradient;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;N;wheelman.map;Mannings n of Wheel tracks (-);wheelman;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Cohesion;wheelcoh.map;Cohesion of wheel tracks (kPa);wheelcohesion;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Ksat;ksatwt.map;Saturated hydraulic conductivity of wheel tracks (mm/h);ksatwt;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Texture classes");
   DEFmaps.append(QString("2;%1;%2;Class 0;mu0.map;Clay fraction (MUST BE CLAY < 2mu);fractionmu0;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Class 1;mu1.map;Soil texture fraction for class 1 (-);fractionmu1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Class 2;mu2.map;Soil texture fraction for class 2 (-);fractionmu2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Class 3;mu3.map;Soil texture fraction for class 3 (-);fractionmu3;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Class 4;mu4.map;Soil texture fraction for class 4 (-);fractionmu4;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Class 5;mu5.map;Soil texture fraction for class 5 (-);fractionmu5;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   j++; i = 0;
   DEFmaps.append("0;Nutrients");
   i = 10;
   DEFmaps.append("1;Phosphorus");
   DEFmaps.append(QString("2;%1;%2;Bulk Dens.;bulkdens.map;Bulk density of the topsoil (kg/m3);bulk;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;P Content;pcont.map;Phosphate (P) content of the soil (kg/kg);pcont;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;P Solute;psolut.map;Initial solute store P in surface layer (kg/m2);psolute;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;P Efficiency;peff.map;Extraction efficiency (s-1);pefficiency;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;P Sorption;Psorp.map;Sorption isotherm kd (m3/kg);psorp;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;P Conversion;Pconv.map;Conversion P from soil content to clay content(-);pconv;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 20;
   DEFmaps.append("1;NH4");
   DEFmaps.append(QString("2;%1;%2;NH4 Content;nh4cont.map;Ammonium (NH4+) content of the soil (kg/kg);nh4cont;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NH4 Solute;nh4solut.map;Initial solute store NH4 in surface layer (kg/m2);nh4solute;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NH4 Efficiency;nh4eff.map;Extraction efficiency (s-1);nh4efficiency;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NH4 Sorption;NH4sorp.map;Sorption isotherm kd (m3/kg);nh4sorp;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NH4 Conversion;NH4conv.map;Conversion NH4 from soil content to clay content(-);nh4conv;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   i = 30;
   DEFmaps.append("1;NO3");
   DEFmaps.append(QString("2;%1;%2;NO3 Content;NO3cont.map;Nitrate (NO3-) content of the soil (kg/kg);no3cont;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NO3 Solute;NO3solut.map;Initial solute store NO3 in surface layer (kg/m2);no3solute;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NO3 Efficiency;NO3eff.map;Extraction efficiency (s-1);no3efficiency;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NO3 Sorption;NO3sorp.map;Sorption isotherm kd (m3/kg);no3sorp;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;NO3 Conversion;NO3conv.map;Conversion NO3 from soil content to clay content(-);no3conv;").arg(j).arg(i)); putmap(k,j,i);k++;i++;

   j++; i = 0;
   DEFmaps.append("0;Gullies");
   DEFmaps.append(QString("2;%1;%2;DEM;dem.map;Digital elevation model (m);dem;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;mannings N;gullyman.map;manning's n gully bottom (-);gullyn;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;BulkDensity;bulkdens.map;Bulkdensity of topsoil (kg/m3);bulkdens1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Ksat;ksat1.map;Ksat of topsoil for gully infil (mm/h);gulksat1;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Depth layer 2;soilDep2.map;Depth to subsoil (cm);gullydep;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Cohesion layer 2;coh2.map;Cohesion of subsoil (kPa);gullycoh;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;BulkDensity 2;bulkden2.map;Bulkdensity of subsoil (kg/m3);bulkdens2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Ksat 2;gulksat2.map;Ksat of subsoil for gully infil (mm/h);gulksat2;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Excluded areas;noncrit.map;areas to be excluded (1) and rest (0);nonfcrit;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Gully initial Width;gulwinit.map; initial gully width (m);gulwinit;").arg(j).arg(i)); putmap(k,j,i);k++;i++;
   DEFmaps.append(QString("2;%1;%2;Gully initial Depth;guldinit.map; initial gully depth (m);guldinit;").arg(j).arg(i)); putmap(k,j,i);k++;i++;

   nrmaplist = k;
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
//also after each call of a runfile so that the runfile mapnames are loadd
void lisemqt::FillMapList()
{
	QStringList headers;

	if (MapNameModel)
	{
		delete MapNameModel;
		MapNameModel = NULL;
	}

	headers  <<"Variable name"<< "Map name"<<"Description"<<" ";
	MapNameModel = new TreeModel(headers, DEFmaps);

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
         qDebug() << d;
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

      mapList[k].name = QFileInfo(path).fileName();
      mapList[k].dir = QFileInfo(path).dir().path();
      // put the name and path into he mapList structure

      QVariant d(mapList[k].name);
      MapNameModel->setData(topLeft, d, Qt::EditRole);
      // put the name into the treeView and MapNameModel
   }
}
//--------------------------------------------------------------------
