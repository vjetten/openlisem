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

#ifndef GL3DWORLDCREATOR_H
#define GL3DWORLDCREATOR_H

#include <ui_full/3D/GL3DWidget.h>
#include <3D/World/GL3DSurface.h>
#include <3D/World/GL3DFlowSurface.h>
#include <3D/World/GL3DSkyBox.h>
#include <3D/World/GL3DCameraController.h>
#include <3D/Objects/GL3DTree.h>
#include <3D/Objects/GL3DPPRain.h>
#include <3D/Objects/GL3DBuilding.h>
#include <3D/Objects/GL3DRoads.h>
#include <3D/Objects/GL3DInstancedTree.h>
#include <3D/Objects/GL3DGrass.h>
#include <3D/Objects/GL3DClouds.h>

struct Output3D
{
    cTMap * MASK;
    cTMap * DEM;
    cTMap * DEM_Filled;
    cTMap * DEM_Change;
    cTMap * SlopeX;
    cTMap * SlopeY;
    cTMap * FlowSlopeX;
    cTMap * FlowSlopeY;
    cTMap * ImageR;
    cTMap * ImageG;
    cTMap * ImageB;
    cTMap * FlowH;
    cTMap * FlowU;
    cTMap * FlowV;
    cTMap * FlowS;
    cTMap * temp;
    cTMap * VegCover;
    cTMap * VegHeight;
    cTMap * SoilCover;
    cTMap * RandomRoughness;
    cTMap * Buildings;
    cTMap * Roads;
    cTMap * ChannelLDD;
    cTMap * ChannelDepth;
    cTMap * ChannelWidth;
    cTMap * ChannelDepthFilled;
    cTMap * ChannelFlowH;
    cTMap * ChannelFlowU;
    cTMap * ChannelFlowS;
    float rainfall;
};

struct Settings3D
{
    QVector4D Light_Ambient;
    QVector4D Light_Directional;
    QVector3D Light_Directional_Direction;
    bool Surface_Draw;
    double Surface_Micro_Elevation_Scale;
    double Surface_Mipmap_Distance_1;
    double Surface_Mipmap_Distance_2;
    QVector3D Surface_Vegetated_Small_Color;
    QVector3D Surface_Vegetated_Large_Color;
    QVector3D Surface_Bare_Color;
    QVector3D Surface_Roads_Color;
    QVector3D Surface_Buildings_Color;
    QVector4D Surface_Erosion_Color;
    QVector4D Surface_Deposition_Color;
    bool Water_Draw;
    double Water_Reflectivity;
    double Water_Refractivity;
    double Water_Velocity_Scale;
    double Water_Micro_Elevation_Scale;
    double Water_Transparancy;
    QVector4D Water_Deep_Color;
    QVector4D Water_Shallow_Color;
    QVector4D Water_Sediment_Color;
    bool Roads_Draw;
    double Roads_Distance;
    bool Buildings_Draw;
    double Buildings_Distance;
    bool Trees_Draw;
    double Trees_Distance;
    double Trees_Instances;
    double Trees_Increment;
    bool Grass_Draw;
    double Grass_Distance;
    double Grass_Instances;
    double Grass_Increment;
    double Grass_Vertical_Scale;
    bool Rain_Draw;
    bool Clouds_Draw;
};


class GL3DWorldCreator;

class GL3DWorldLoader : public QThread
{
public:

    GL3DWorldCreator * creator;


private:
    void run();

};


class GL3DWorldCreator
{


    bool LisemWorldCreated = false;

public:

    GL3DWorldCreator()
    {
        loader = new GL3DWorldLoader();
    }

    GL3DWorldLoader * loader;
    QMutex * m_Mutex;

   inline void LinkToWidget(GL3DWidget * w)
   {
       m_Widget = w;
   }

   bool DoneLoading = false;

   GL3DWidget * m_Widget;

   bool done_creating = false;

   void CreateWorldFromLisem();
   void CreateWorldFromLisemThread();
   void UpdateWorldFromLisem();
   void UpdateWorldSettings(Settings3D s);
   void DestroyWorldFromLisem();

   //maps from lisem model
   Output3D input;

   //Settings from user interface
   Settings3D settings;

   //all used objects
   GL3DSurface * surface;
   GL3DFlowSurface * fsurface;
   GL3DCameraController * controller;
   GL3DSkyBox * skybox;
   GL3DPPRain * rain;
   GL3DInstancedTree * treesobjecti;
   GL3DGrass * grassobjecti;
   GL3DBuildings * buildingobject;
   GL3DRoads * roadsobject;
   GL3DClouds * cloudsobject;
};


#endif // GL3DWORLDCREATOR_H
