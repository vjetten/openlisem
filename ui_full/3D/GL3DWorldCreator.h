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
#include <3D/Objects/GL3DSurface.h>
#include <3D/Objects/GL3DFlowSurface.h>
#include <3D/Objects/GL3DSkyBox.h>
#include <3D/World/GL3DCameraController.h>

struct Output3D
{
    cTMap * MASK;
    cTMap * DEM;
    cTMap * DEM_Filled;
    cTMap * SlopeX;
    cTMap * SlopeY;
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

};


class GL3DWorldCreator
{


    bool LisemWorldCreated = false;

public:

    GL3DWorldCreator()
    {

    }

   inline void LinkToWidget(GL3DWidget * w)
   {
       m_Widget = w;
   }

   GL3DWidget * m_Widget;

   void CreateWorldFromLisem();
   void UpdateWorldFromLisem();
   void DestroyWorldFromLisem();

   //maps from lisem model
   Output3D input;

   //all used objects
   GL3DSurface * surface;
   GL3DFlowSurface * fsurface;
   GL3DCameraController * controller;
   GL3DSkyBox * skybox;

};



#endif // GL3DWORLDCREATOR_H
