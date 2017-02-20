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
#include <3D/Objects/GL3DTree.h>
#include <3D/Objects/GL3DBuilding.h>
#include <3D/Objects/GL3DRoads.h>

void GL3DWorldLoader::run()
{
    qDebug() << "start loading thread";

    creator->m_Mutex->lock();


    creator->CreateWorldFromLisemThread();

    creator->m_Mutex->unlock();

    qDebug() << "End loading thread";
}


void GL3DWorldCreator::CreateWorldFromLisemThread()
{
    qDebug() <<"create lisem world 3d thread";
    this->m_Widget->context()->makeCurrent();

    MapMath::FillDem(input.DEM,input.DEM_Filled,input.temp,50);
    MapMath::FillDem(input.ChannelDepth,input.ChannelDepthFilled,input.temp,10);
    MapMath::SlopeMap(input.DEM,input.SlopeX,input.SlopeY);

    m_Widget->m_World->RemoveAllObjects();

    GL3DSkyBox * sb = new GL3DSkyBox();
    m_Widget->m_World->AddObject(sb,false,true,false);

    surface=  new GL3DSurface();
    surface->SetSurfaceMap(input.DEM,input.DEM_Filled,input.SlopeX,input.SlopeY);
    surface->SetTerrainProperties(input.VegCover,input.VegHeight,input.RandomRoughness,input.Buildings,input.Roads);

    if(op.has_channel)
    {
        surface->SetChannel(m_Widget,input.ChannelLDD,input.ChannelWidth,input.ChannelDepthFilled,input.temp,op.has_channelflooding);
    }
    m_Widget->m_World->AddObject(surface,true);

    GL3DRoads * roadsobject = new GL3DRoads();
    roadsobject->SetRoadDistribution(m_Widget,surface,input.Roads);
    m_Widget->m_World->AddObject(roadsobject);

    fsurface = new GL3DFlowSurface();
    fsurface->SetSurface(surface);
    fsurface->SetFlowProperties(input.FlowH,input.FlowU,input.FlowV,input.FlowS);
    fsurface->SetFlowPropertiesChannel(input.ChannelFlowH,input.ChannelFlowU,input.ChannelFlowS);
    fsurface->SetSkyBox(sb);
    m_Widget->m_World->AddObject(fsurface,false,false,true);

    controller = new GL3DCameraController();
    controller->SetSurface(surface);
    controller->SetCamera(m_Widget->m_Camera);
    controller->SetStartPosition();
    m_Widget->m_World->AddObject(controller);

    GL3DTrees * treesobject = new GL3DTrees();
    treesobject->SetTreeDistribution(surface,input.VegCover,input.VegHeight);
    m_Widget->m_World->AddObject(treesobject);

    GL3DBuildings * buildingobject = new GL3DBuildings();
    buildingobject->SetBuildingsDistribution(surface,input.Buildings,input.temp);
    m_Widget->m_World->AddObject(buildingobject);

    this->rain = new GL3DPPRain();
    m_Widget->m_World->AddObject(rain);

    m_Widget->ReadyToDraw = true;
    done_creating = true;

    m_Widget->context()->doneCurrent();
    m_Widget->context()->moveToThread(m_Widget->thread());

    m_Widget->gl_context_control->lock();
    m_Widget->gl_context_try = true;
    m_Widget->setUpdatesEnabled(true);
    m_Widget->gl_context_control->unlock();


}

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
    input.DEM_Change = new cTMap();
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
    input.ChannelLDD = new cTMap();
    input.ChannelDepth = new cTMap();
    input.ChannelDepthFilled = new cTMap();
    input.ChannelWidth = new cTMap();
    input.ChannelFlowH = new cTMap();
    input.ChannelFlowU = new cTMap();
    input.ChannelFlowS = new cTMap();

    input.MASK->MakeMap(op.baseMapDEM,0.0);
    input.DEM->MakeMap(op.baseMapDEM,0.0);
    input.DEM_Filled->MakeMap(op.baseMapDEM,0.0);
    input.DEM_Change->MakeMap(op.baseMapDEM,0.0);
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
    input.ChannelLDD->MakeMap(op.baseMapDEM,0.0);
    input.ChannelDepth->MakeMap(op.baseMapDEM,0.0);
    input.ChannelDepthFilled->MakeMap(op.baseMapDEM,0.0);
    input.ChannelWidth->MakeMap(op.baseMapDEM,0.0);
    input.ChannelFlowH->MakeMap(op.baseMapDEM,0.0);
    input.ChannelFlowU->MakeMap(op.baseMapDEM,0.0);
    input.ChannelFlowS->MakeMap(op.baseMapDEM,0.0);

    input.DEM_Filled->setAllMV();
    input.ChannelDepth->setAllMV();
    input.ChannelDepthFilled->setAllMV();

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
        input.DEM_Change->Drc = op.gl_dem_change->Drc;
        if(op.has_channel)
        {
            if(!pcr::isMV(op.ch_ldd->Drc))
            {
               input.ChannelDepth->Drc = op.ch_d->Drc;
            }
            input.ChannelLDD->Drc = op.ch_ldd->Drc;
            input.ChannelWidth->Drc = op.ch_w->Drc;
            input.ChannelFlowH->Drc = op.gl_ch_flow_height->Drc;
            input.ChannelFlowU->Drc = op.gl_ch_flow_v->Drc;
            input.ChannelFlowS->Drc = op.gl_ch_flow_c->Drc;
        }

    }

    input.rainfall = op.rain_average;

    m_Widget->gl_context_control->lock();
    m_Widget->gl_context_try = false;
    m_Widget->setUpdatesEnabled(false);
    m_Widget->gl_context_control->unlock();

    m_Mutex = new QMutex();
    m_Mutex->lock();

    loader->creator = this;
    loader->start();

    m_Widget->context()->doneCurrent();
    m_Widget->context()->moveToThread(loader);

    m_Mutex->unlock();

    //loader->wait();
    //this->m_Widget->context()->makeCurrent();

}

void GL3DWorldCreator::UpdateWorldFromLisem()
{
    if(done_creating)
    {

        FOR_ROW_COL_MV(input.DEM,input.MASK)
        {
            input.FlowH->Drc = op.gl_flow_height->Drc;
            input.FlowU->Drc = op.gl_flow_u->Drc;
            input.FlowV->Drc = op.gl_flow_v->Drc;
            input.FlowS->Drc = op.gl_flow_c->Drc;
            input.DEM_Change->Drc = op.gl_dem_change->Drc;

            if(op.has_channel)
            {
                input.ChannelFlowH->Drc = op.gl_ch_flow_height->Drc;
                input.ChannelFlowU->Drc = op.gl_ch_flow_v->Drc;
                input.ChannelFlowS->Drc = op.gl_ch_flow_c->Drc;
            }
        }

        FOR_ROW_COL_MV(input.DEM,input.MASK)
        {
            if(pcr::isMV(op.ch_ldd->Drc))
            {
                pcr::setMV(op.gl_ch_flow_height->Drc);
                pcr::setMV(op.gl_ch_flow_v->Drc);
                pcr::setMV(op.gl_ch_flow_c->Drc);
            }
        }

        MapMath::FillDem(op.gl_ch_flow_height,input.ChannelFlowH,input.temp,10);
        MapMath::FillDem(op.gl_ch_flow_v,input.ChannelFlowU,input.temp,10);
        MapMath::FillDem(op.gl_ch_flow_c,input.ChannelFlowS,input.temp,10);

        input.rainfall = op.rain_average;

        fsurface->SetFlowProperties(input.FlowH,input.FlowU,input.FlowV,input.FlowS);
        fsurface->SetFlowPropertiesChannel(input.ChannelFlowH,input.ChannelFlowU,input.ChannelFlowS);
        rain->SetRainfall(input.rainfall);
    }
}

void GL3DWorldCreator::DestroyWorldFromLisem()
{
    qDebug() <<"destroy lisem world 3d";

    m_Widget->m_World->RemoveAndDestroyAllObjects();
    LisemWorldCreated = false;
}

