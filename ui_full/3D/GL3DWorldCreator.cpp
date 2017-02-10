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

void GL3DWidget::CreateWorldFromLisem()
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
    input.WaterHeight = new cTMap();
    input.temp = new cTMap();
    input.VegCover = new cTMap();
    input.VegHeight = new cTMap();
    input.SoilCover = new cTMap();
    input.RandomRoughness = new cTMap();

    input.MASK->MakeMap(op.baseMapDEM,0.0);
    input.DEM->MakeMap(op.baseMapDEM,0.0);
    input.DEM_Filled->MakeMap(op.baseMapDEM,0.0);
    input.SlopeX->MakeMap(op.baseMapDEM,0.0);
    input.SlopeY->MakeMap(op.baseMapDEM,0.0);
    input.ImageR->MakeMap(op.baseMapDEM,0.0);
    input.ImageG->MakeMap(op.baseMapDEM,0.0);
    input.ImageB->MakeMap(op.baseMapDEM,0.0);
    input.WaterHeight->MakeMap(op.baseMapDEM,0.0);
    input.temp->MakeMap(op.baseMapDEM,0.0);
    input.VegCover->MakeMap(op.baseMapDEM,0.0);
    input.VegHeight->MakeMap(op.baseMapDEM,0.0);
    input.SoilCover->MakeMap(op.baseMapDEM,0.0);
    input.RandomRoughness->MakeMap(op.baseMapDEM,0.0);

    FOR_ROW_COL_MV(input.DEM,input.MASK)
    {
        input.MASK->Drc = op.baseMapDEM->Drc;
        input.DEM->Drc = op.baseMapDEM->Drc;

    }

    MapMath::FillDem(input.DEM,input.DEM_Filled,input.temp,50);
    MapMath::SlopeMap(input.DEM,input.SlopeX,input.SlopeY);

    this->m_World->RemoveAllObjects();

    m_World->AddObject(new GL3DSkyBox());

    GL3DSurface * surface=  new GL3DSurface();
    surface->SetSurfaceMap(input.DEM,input.DEM_Filled,input.SlopeX,input.SlopeY);

    m_World->AddObject(surface);

    GL3DCameraController * controller = new GL3DCameraController();
    controller->SetSurface(surface);
    controller->SetCamera(m_Camera);

    controller->SetStartPosition();

    m_World->AddObject(controller);
}

void GL3DWidget::UpdateWorldFromLisem()
{



}

void GL3DWidget::DestroyWorldFromLisem()
{
    qDebug() <<"destroy lisem world 3d";

    this->m_World->RemoveAndDestroyAllObjects();
    LisemWorldCreated = false;
}

