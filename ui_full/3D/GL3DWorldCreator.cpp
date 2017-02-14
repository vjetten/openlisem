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

#include <ui_full/3D/GL3DWidget.h>
#include <3D/Objects/GL3DSurface.h>
#include <3D/Objects/GL3DSkyBox.h>
#include <3D/World/GL3DCameraController.h>
#include <3D/GL3DWorldCreator.h>

void GL3DWorldCreator::CreateWorldFromLisem()
{

    qDebug() <<"create lisem world 3d";
    if(LisemWorldCreated)
    {
        DestroyWorldFromLisem();
    }
    LisemWorldCreated = true;

    input.MASK = new cTMap();
    input.DEM = new cTMap();
    input.DEM_Filled = new cTMap();
    input.SlopeX = new cTMap();
    input.SlopeY = new cTMap();
    input.ImageR = new cTMap();
    input.ImageG = new cTMap();
    input.ImageB = new cTMap();
    input.FlowH = new cTMap();
    input.FlowU = new cTMap();
    input.FlowV = new cTMap();
    input.FlowS = new cTMap();
    input.temp = new cTMap();
    input.VegCover = new cTMap();
    input.VegHeight = new cTMap();
    input.SoilCover = new cTMap();
    input.RandomRoughness = new cTMap();
    input.Buildings = new cTMap();
    input.Roads = new cTMap();

    input.MASK->MakeMap(op.baseMapDEM,0.0);
    input.DEM->MakeMap(op.baseMapDEM,0.0);
    input.DEM_Filled->MakeMap(op.baseMapDEM,0.0);
    input.SlopeX->MakeMap(op.baseMapDEM,0.0);
    input.SlopeY->MakeMap(op.baseMapDEM,0.0);
    input.ImageR->MakeMap(op.baseMapDEM,0.0);
    input.ImageG->MakeMap(op.baseMapDEM,0.0);
    input.ImageB->MakeMap(op.baseMapDEM,0.0);
    input.FlowH->MakeMap(op.baseMapDEM,0.0);
    input.FlowU->MakeMap(op.baseMapDEM,0.0);
    input.FlowV->MakeMap(op.baseMapDEM,0.0);
    input.FlowS->MakeMap(op.baseMapDEM,0.0);
    input.temp->MakeMap(op.baseMapDEM,0.0);
    input.VegCover->MakeMap(op.baseMapDEM,0.0);
    input.VegHeight->MakeMap(op.baseMapDEM,0.0);
    input.SoilCover->MakeMap(op.baseMapDEM,0.0);
    input.RandomRoughness->MakeMap(op.baseMapDEM,0.0);
    input.Buildings->MakeMap(op.baseMapDEM,0.0);
    input.Roads->MakeMap(op.baseMapDEM,0.0);

    FOR_ROW_COL_MV(input.DEM,input.MASK)
    {
        input.MASK->Drc = op.baseMapDEM->Drc;
        input.DEM->Drc = op.baseMapDEM->Drc;
        input.VegCover->Drc = op.vegcover->Drc;
        input.VegHeight->Drc = op.vegheight->Drc;
        input.SoilCover->Drc = (1.0-op.vegcover->Drc);
        input.RandomRoughness->Drc = op.randomroughness->Drc;
        input.Buildings->Drc = op.houseMap->Drc;
        input.Roads->Drc = op.roadMap->Drc;
        input.FlowH->Drc = op.gl_flow_height->Drc;
        input.FlowU->Drc = op.gl_flow_u->Drc;
        input.FlowV->Drc = op.gl_flow_v->Drc;
        input.FlowS->Drc = op.gl_flow_c->Drc;
    }

    MapMath::FillDem(input.DEM,input.DEM_Filled,input.temp,50);
    MapMath::SlopeMap(input.DEM,input.SlopeX,input.SlopeY);

    m_Widget->m_World->RemoveAllObjects();

    m_Widget->m_World->AddObject(new GL3DSkyBox());

    surface=  new GL3DSurface();
    surface->SetSurfaceMap(input.DEM,input.DEM_Filled,input.SlopeX,input.SlopeY);
    surface->SetTerrainProperties(input.VegCover,input.VegHeight,input.RandomRoughness,input.Buildings,input.Roads);

    m_Widget->m_World->AddObject(surface);

    fsurface = new GL3DFlowSurface();
    fsurface->SetSurface(surface);
    fsurface->SetFlowProperties(input.FlowH,input.FlowU,input.FlowV,input.FlowS);

    m_Widget->m_World->AddObject(fsurface);

    controller = new GL3DCameraController();
    controller->SetSurface(surface);
    controller->SetCamera(m_Widget->m_Camera);

    controller->SetStartPosition();

    m_Widget->m_World->AddObject(controller);
}

void GL3DWorldCreator::UpdateWorldFromLisem()
{

    FOR_ROW_COL_MV(input.DEM,input.MASK)
    {
        input.FlowH->Drc = op.gl_flow_height->Drc;
        input.FlowU->Drc = op.gl_flow_u->Drc;
        input.FlowV->Drc = op.gl_flow_v->Drc;
        input.FlowS->Drc = op.gl_flow_c->Drc;
    }

    fsurface->SetFlowProperties(input.FlowH,input.FlowU,input.FlowV,input.FlowS);

}

void GL3DWorldCreator::DestroyWorldFromLisem()
{
    qDebug() <<"destroy lisem world 3d";

    m_Widget->m_World->RemoveAndDestroyAllObjects();
    LisemWorldCreated = false;
}

