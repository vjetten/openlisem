
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file lisReportfile.cpp
  \brief reporting maps, hydrographs, outlet/area totals and land unit stats

functions: \n
- void TWorld::OutputUI() fill output structure 'op' with results to talk to the interface.\n
- void TWorld::reportAll() \n
- void TWorld::ReportTimeseriesNew() report outlet data to textfile\n
- void TWorld::ReportTotalsNew() report totals to text file\n
- void TWorld::ReportMaps() report maps and mapseries\n
- void TWorld::CountLandunits() make a list of landunit numbers\n
- void TWorld::ReportLandunits() report text data per landunit\n
 */

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

//---------------------------------------------------------------------------
// Make the maps to bedrawn in the interface as a copy in the op starcture
// reason is that all pointers are destroyed after the run so when lisem finishes
// the information on the output screen points to an empty pointer
// by copying the info remains available
// initialize maps for output to screen
// must be done after Initialize Data because then we know how large the map is
void TWorld::setupDisplayMaps()
{
    if (op.baseMap != 0)
    {
        delete op.baseMap;
        delete op.baseMapDEM;
        delete op.channelMap;
        delete op.outletMap;
        delete op.roadMap;
        delete op.houseMap;
        delete op.hardsurfaceMap;
    }

    op.baseMap = new cTMap();
    op.baseMapDEM = new cTMap();
    op.channelMap = new cTMap();
    op.outletMap = new cTMap();
    op.roadMap = new cTMap();
    op.houseMap = new cTMap();
    op.hardsurfaceMap = new cTMap();

    op.baseMap->MakeMap(LDD, 0);
    op.baseMapDEM->MakeMap(LDD, 0);
    op.channelMap->MakeMap(LDD, 0);
    op.outletMap->MakeMap(LDD, 0);
    op.roadMap->MakeMap(LDD, 0);
    op.houseMap->MakeMap(LDD, 0);
    op.hardsurfaceMap->MakeMap(LDD, 0);
}
//---------------------------------------------------------------------------
void TWorld::setLegendColors()
{
    Legend.clear();
    LegendMap.clear();

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.01);
    Colormap.append(0.1);
    Colormap.append(0.3);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#9eccee");
    Colors.append("#427dc6");
    Colors.append("#204ab5");
    Colors.append("#072a9c");
    Colors.append("#df2a36");

    Legend<<Colors; //0
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.25);
    Colormap.append(0.5);
    Colormap.append(0.75);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#fdecc6");//ffffb2");
    Colors.append("#fecc5c");
    Colors.append("#fd8d3c");
    Colors.append("#f03b20");
    Colors.append("#bd0026");

    Legend<<Colors; //1
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.25);
    Colormap.append(0.5);
    Colormap.append(0.75);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#9eccee");
    Colors.append("#427dc6");
    Colors.append("#204ab5");
    Colors.append("#072a9c");
    Colors.append("#07215e");

    Legend<<Colors; //2
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.15);
    Colormap.append(0.3);
    Colormap.append(0.5);
    Colormap.append(0.75);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#ffff99");
    Colors.append("#ffff51");
    Colors.append("#c7e55a");
    Colors.append("#32b1df");
    Colors.append("#3271ca");
    Colors.append("#2c3898");

    Legend<<Colors; //3
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.001);
    Colormap.append(0.5);
    Colormap.append(1.0);
    Colors.clear();

    Colors.append("#77aa66");
    Colors.append("#32b1df");
    Colors.append("#3271ca");

    Legend<<Colors; //4
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.25);
    Colormap.append(0.5);
    Colormap.append(0.75);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#d7191c");
    Colors.append("#fdae61");
    Colors.append("#fdfd7e");
    Colors.append("#abdda4");
    Colors.append("#2b83ba");

    Legend<<Colors; //5
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
  //  Colormap.append(0.25);
    Colormap.append(0.49);
    Colormap.append(0.5);
    Colormap.append(0.51);
    Colormap.append(0.75);
    Colormap.append(1.0);

    Colors.clear();
    Colors.append("#2b83ba");
   // Colors.append("#a4ddd9");
    Colors.append("#fffff0");
    Colors.append("#dddddd");
    Colors.append("#fffff0");
    Colors.append("#d3b03e");
    Colors.append("#d7191c");

    Legend<<Colors; //6
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.3);
    Colormap.append(0.5);
    Colormap.append(0.70);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#616ca2");
    Colors.append("#50B547");
    Colors.append("#FFFFFF");
    Colors.append("#ffff88");
    Colors.append("#FF0000");

    Legend<<Colors; //7
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
    Colormap.append(0.25);
    Colormap.append(0.5);
    Colormap.append(0.75);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#fffff0");
    Colors.append("#fed98e");
    Colors.append("#fe9929");
    Colors.append("#d95f0e");
    Colors.append("#a94400");


    Legend<<Colors; //8
    LegendMap << Colormap;

    Colormap.clear();
    Colormap.append(0.0);
 //   Colormap.append(0.25);
  //  Colormap.append(0.5);
  //  Colormap.append(0.75);
    Colormap.append(1.0);
    Colors.clear();
    Colors.append("#fffff0");
  //  Colors.append("#bae4bc");
   // Colors.append("#7bccc4");
  //  Colors.append("#43a2ca");
    Colors.append("#2b83ba");

    Legend<<Colors; //9
    LegendMap << Colormap;

}
//---------------------------------------------------------------------------

// setup display list of maps in combo boxes
//create combo box maps, once before run
void TWorld::GetComboMaps()
{
    ClearComboMaps();

    setLegendColors();

    int cl = 0;
    if (QUnits == 0)
        AddComboMap(0,"Total Discharge","l/s",Qoutput,LegendMap[cl],Legend[cl],true,false,1.0, 1.0);
    else
        AddComboMap(0,"Total Discharge","m3/s",Qoutput,LegendMap[cl],Legend[cl],true,false,1.0, 0.001);
    //factor is already done in Qoutput, so that reportfile is also done, not only screen
    // the only thing that needs to change here is the text "m3/s"
    //AddComboMap(0,"Added Discharge","m3/s",Qbase,LegendMap[cl],Legend[cl],true,false,1.0, 1.0);//0.001);

    cl = 2;
    AddComboMap(0,"Water Height","m",hmxWH,LegendMap[cl],Legend[cl],false,false,1.0,0.001);
 //   AddComboMap(0,"Water inflow","m3",ChannelQSide,LegendMap[cl],Legend[cl],true,false,1.0,1.0);
//    if (Switch2DDiagonalFlow)
//       AddComboMap(0,"Diagonal Discharge","l/s",Qdiag,LegendMap[cl],Legend[cl],false,false,1.0, 0.01);
    cl = 1;
    AddComboMap(0,"Overland flow Velocity","m/s",COMBO_V,LegendMap[cl],Legend[cl],false,false,1.0, 0.001);
 //   AddComboMap(0,"Flow Velocity","m/s",Uflood,LegendMap[cl],Legend[cl],false,false,1.0, 0.01);
 //   AddComboMap(0,"Flow Velocity","m/s",Vflood,LegendMap[cl],Legend[cl],false,false,1.0, 0.01);
 //   AddComboMap(0,"Flow Velocity","m/s",K2DOutlets,LegendMap[cl],Legend[cl],false,false,1.0, 0.01);
    AddComboMap(0,"Overland flow Momentum","m2/s",VH,LegendMap[cl],Legend[cl],false,false,1.0, 0.001); //VH
    //AddComboMap(0,"boundary","-",K2DOutlets,LegendMap[cl],Legend[cl],false,false,1.0, 0.01);

    AddComboMap(0,"Cumulative overland flow","m3",Qm3total,LegendMap[0],Legend[0],false,false,1.0, 1.0);//0.001);
    if (SwitchKinematic2D == K2D_METHOD_DYN || SwitchKinematic2D == K2D_METHOD_KINDYN) {
        cl = 2;
        QString txt = QString("Max flood Height (h>%1m)").arg(minReportFloodHeight);
        AddComboMap(0,txt,"m",floodHmxMax,LegendMap[cl],Legend[cl],false,false,1.0,0.01);
        cl = 5;
        AddComboMap(0,"Flood Start Time","min",floodTimeStart,LegendMap[cl],Legend[cl],false,false,1.0,1.0);
        cl = 6;
        AddComboMap(0,"Flood duration","min",floodTime,LegendMap[cl],Legend[cl],false,false,1.0,1.0);
        cl = 0;
        AddComboMap(0,"Flood Hazard Index [WH(V+0.5)]","-",FHI,LegendMap[0],Legend[0],true,false,1.0, 0.001);
    }


    if(SwitchIncludeChannel)
    {

        cl = 0;
        if (QUnits == 0)
            AddComboMap(0,"Channel Discharge","l/s",ChannelQn,LegendMap[cl],Legend[cl],true,false,1000.0, 1.0);
        else
            AddComboMap(0,"Channel Discharge","m3/s",ChannelQn,LegendMap[cl],Legend[cl],true,false,1.0, 0.001);
        cl = 2;
        AddComboMap(0,"Channel Water Height","m",ChannelWH,LegendMap[cl],Legend[cl],false,false,1.0,0.001);
//        if (SwitchChannelBaseflow)
//            AddComboMap(0,"Baseflow inflow","m3/s",Qbase,LegendMap[cl],Legend[cl],false,false,1.0,0.01);

        cl = 1;
        AddComboMap(0,"Channel Velocity","m/s",ChannelV,LegendMap[cl],Legend[cl],false,false,1.0,0.001);
    }

    if(SwitchIncludeTile || SwitchIncludeStormDrains) {
        cl = 0;
        AddComboMap(0,"Storm Drain Volume","m3",TileWaterVol,LegendMap[cl],Legend[cl],false,false,1.0,0.001);
        AddComboMap(0,"Storm Drain Discharge","m3/s",TileQn,LegendMap[cl],Legend[cl],true,false,1.0,0.001);
    }

    cl = 3;
    AddComboMap(0,"Interception","mm",InterceptionmmCum,LegendMap[cl],Legend[cl],false,false,1.0,1.0);

    if(SwitchInfiltration)
    {
        AddComboMap(0,"Infiltration","mm",InfilmmCum,LegendMap[cl],Legend[cl],false,false,1.0,1.0);
        if (InfilMethod > 1) {
            // Show a weighed avregae of cell wetting front else it seems partly impermeable surface infiltrate very deep
            AddComboMap(0,"Average depth wetting front","mm",Lwmm,LegendMap[cl],Legend[cl],false,false,1.0,1.0);  // swatre?
        }

        if (SwitchGWflow) {
            AddComboMap(0,"Groundwater level","m",GWWH,LegendMap[cl],Legend[cl],false,false,1.0,0.001);
            AddComboMap(0,"Groundwater level max","m",GWWHmax,LegendMap[cl],Legend[cl],false,false,1.0,0.001);
            //AddComboMap(0,"SD2","m",SoilDepth2,LegendMap[cl],Legend[cl],false,false,1.0,0.001);
        }
        cl = 6;
        if (SwitchSlopeStability)
            AddComboMap(0,"Slope Stability","m",FSlope,LegendMap[cl],Legend[cl],false,false,1.0,0.001);


        if (InfilMethod != INFIL_SWATRE) {
            cl = 3;
            //AddComboMap(0,"Avg Moisture content layer 1","-",Thetaeff,LegendMap[cl],Legend[cl],false,false,1.0,1.0);
            AddComboMap(0,"Avg Moisture content layer 1","-",ThetaI1a,LegendMap[cl],Legend[cl],false,false,1.0,1.0);
            if (SwitchTwoLayer)
                AddComboMap(0,"Avg Moisture content layer 2","-",ThetaI2a,LegendMap[cl],Legend[cl],false,false,1.0,1.0);
            if (!SwitchImpermeable || SwitchGWflow)
                AddComboMap(0,"Percolation","mm",Perc,LegendMap[cl],Legend[cl],false,false,1000,1.0);
        }
    }

    cl = 4;
    double factor = 3600000.0/_dt; //from m to mm/h

    AddComboMap(0,"Rainfall Cumulative","mm",RainCumFlat,LegendMap[cl],Legend[cl],false,false,1000.0,0.1);
    AddComboMap(0,"Rainfall Intensity","mm/h",Rain,LegendMap[cl],Legend[cl],false,false,factor,0.1);
    //  AddComboMap(0,"ETa cumulative","mm",ETa,LegendMap[cl],Legend[cl],false,false,1000.0,0.1);
    if (SwitchIncludeET) {
        AddComboMap(0,"ETa Cumulative","mm",ETaCum,LegendMap[3],Legend[3],false,false,1000.0,0.1);
       // AddComboMap(0,"ETp Cumulative","mm",ETpCum,LegendMap[3],Legend[3],false,false,1000.0,0.1);
       }

    if(SwitchErosion)
    {
        double step = 0.01;
        cl = 6;

        QString unit = "kg/cell";
        double factor = 1.0;
        if(ErosionUnits == 2)
        {
            factor = 1.0/(_dx*_dx);
            unit = "kg/m2";
        }else if (ErosionUnits == 0)
        {
            factor = 10.0/(_dx*_dx);
            unit = "t/ha";
        }
        AddComboMap(1,"Total Soil Loss",unit,TotalSoillossMap,LegendMap[cl],Legend[cl],false,true,factor, step);

        cl = 8;
        AddComboMap(1,"Splash detachment",unit,DETSplashCum,LegendMap[cl],Legend[cl],false,false,factor, step);
        AddComboMap(1,"Flow detachment",unit,DETFlowCum,LegendMap[cl],Legend[cl],false,false,factor, step);
        cl = 9;
        AddComboMap(1,"Deposition",unit,DEPCum,LegendMap[cl],Legend[cl],false,false,-factor, step);

        cl = 8;
        AddComboMap(1,"Sed. Concentration","kg/m3",TotalConc,LegendMap[cl],Legend[cl],false,false,1.0, step);
        if (SwitchSedtrap)
            AddComboMap(1,"Sed trap","kg/m3",SedMaxVolume,LegendMap[cl],Legend[cl],false,false,1.0, step);

        double factor_g = 1000/(_dx*_dx);
        QString unit_g = "g/m2";
        AddComboMap(1,"Suspended sed.",unit_g,COMBO_SS,LegendMap[cl],Legend[cl],false,false,factor_g, step);

        AddComboMap(1,"TC suspended","kg/m3",COMBO_TC,LegendMap[cl],Legend[cl],false,false,1.0, step);
        if(SwitchUse2Phase) {
            AddComboMap(1,"Bedload sed.",unit_g,COMBO_BL,LegendMap[cl],Legend[cl],false,false,factor_g, step);
         //   AddComboMap(1,"TC bedload","kg/m3",BLTCFlood,LegendMap[cl],Legend[cl],false,false,1.0, step);
         //   AddComboMap(1,"SS depth","m",SSDepthFlood,LegendMap[cl],Legend[cl],false,false,1.0, step);
         //   AddComboMap(1,"BL depth","m",BLDepthFlood,LegendMap[cl],Legend[cl],false,false,1.0, step);
        }

        // cl = 9;
        // if(SwitchUseMaterialDepth) {
        //     AddComboMap(1,"Storage",unit,Storage,LegendMap[cl],Legend[cl],false,false,-factor, step);
        //     AddComboMap(1,"Storage",unit,StorageDep,LegendMap[cl],Legend[cl],false,false,-factor, step);
        // }
    }
}
//---------------------------------------------------------------------------
void TWorld::ClearComboMaps()
{

//    for(int i =op.ComboMapsSafe.length() - 1; i >-1 ; i--)
//    {
//        delete op.ComboMapsSafe.at(i);
//    }
//    op.ComboMapsSafe.clear();

    op.ComboLists.clear();
    op.ComboMaps.clear();
    op.ComboColorMap.clear();
    op.ComboColors.clear();
    op.ComboLogaritmic.clear();
    op.ComboSymColor.clear();
    op.ComboMapNames.clear();
    op.ComboUnits.clear();
    op.ComboScaling.clear();
    op.userMinV.clear();
    op.userMaxV.clear();
    op.comboStep.clear();

    op.comboboxset = false;
}
//---------------------------------------------------------------------------
//because this is moved outside the timeloop and the pointer of a map is copied to this structure (not the map itself),
// you cannot display this after the original maps are destroyed
void TWorld::AddComboMap(int listn, QString name, QString unit,cTMap * map,QList<double> ColorMap, QList<QString> Colors,
                         bool log,bool symcol, double scale, double step)
{
    op.ComboLists.append(listn);
    op.ComboMaps.append(map);
    // copy pointer or make a map and copy content
    //op.ComboMapsSafe.append(new cTMap());
    //op.ComboMapsSafe.at(op.ComboMapsSafe.length()-1)->MakeMap(LDD,0.0);

    op.ComboColorMap.append(ColorMap);
    op.ComboColors.append(Colors);
    op.ComboLogaritmic.append(log);
    op.ComboSymColor.append(symcol);
    op.ComboMapNames.append(name);
    op.ComboUnits.append(unit);
    op.ComboScaling.append(scale);  // multiplier for display
    op.userMinV.append(0);  //initialize to 0, used to save users choice
    op.userMaxV.append(0);
    op.comboStep.append(step);

    op.comboboxset = false;
}

void TWorld::CopyComboMap(int i, cTMap * map)
{
    op.ComboMaps.append(new cTMap);
    op.ComboMaps.at(i)->MakeMap(LDD,0.0);
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        op.ComboMaps.at(i)->Drc = map->Drc;
    }}
    // copy the map content
}
