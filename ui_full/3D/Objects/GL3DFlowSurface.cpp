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

#include <3D/Graphics/GL3DTextures.h>
#include <3D/Objects/GL3DFlowSurface.h>

void GL3DFlowSurface::SetSurface(GL3DSurface * surface)
{
    m_Surface = surface;
}

void GL3DFlowSurface::SetFlowProperties(cTMap * h, cTMap * u, cTMap * v, cTMap * s)
{
    m_FlowH = h;
    m_FlowU = u;
    m_FlowV = v;
    m_FlowS = s;

    text_updated = true;

}

void GL3DFlowSurface::CreateTextures(GL3DWidget * widget)
{
    if(text_updated || !text_created)
    {
        if(!text_created)
        {
            if(m_FlowH != 0)
            {
                m_Texture_FlowH = widget->m_Textures->LoadTextureFromMap(false,m_FlowH,0,0,true);
                m_Texture_FlowU = widget->m_Textures->LoadTextureFromMap(false,m_FlowH,0,0,true);
                m_Texture_FlowV = widget->m_Textures->LoadTextureFromMap(false,m_FlowH,0,0,true);
                m_Texture_FlowS = widget->m_Textures->LoadTextureFromMap(false,m_FlowH,0,0,true);

                text_created = true;
                text_updated = false;
            }else
            {
                text_created = false;

            }

        }else if(m_FlowH != 0)
        {

            m_Texture_FlowH->UpdateTextureFromMap(widget, m_FlowH,0,0,true);
            m_Texture_FlowU->UpdateTextureFromMap(widget, m_FlowH,0,0,true);
            m_Texture_FlowV->UpdateTextureFromMap(widget, m_FlowH,0,0,true);
            m_Texture_FlowS->UpdateTextureFromMap(widget, m_FlowH,0,0,true);

            text_created = true;
            text_updated = false;
        }else
        {
            text_created = false;

        }
    }
}

void GL3DFlowSurface::OnCreate(GL3DWidget * widget)
{

    this->recieve_render = true;
    this->recieve_render_late = true;
    m_Shader_Flow = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SURFACE_FLOW);
}

void GL3DFlowSurface::OnRenderLate(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{
    widget->gl->glDisable(GL_BLEND);

    this->CreateTextures(widget);

    QMatrix4x4 matrix;
    matrix.setToIdentity();

    m_Shader_Flow->m_program->bind();
    m_Shader_Flow->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    m_Shader_Flow->m_program->setUniformValue("Mmatrix",matrix);
    m_Shader_Flow->m_program->setUniformValue("Cposition",camera->m_Position);

    m_Shader_Flow->m_program->setUniformValue("microElevationScaleX",(float) m_Surface->m_CellSize);
    m_Shader_Flow->m_program->setUniformValue("microElevationScaleY",(float) 0.25f);
    m_Shader_Flow->m_program->setUniformValue("microElevationScaleZ",(float) m_Surface->m_CellSize);

    m_Shader_Flow->m_program->setUniformValue("SurfaceExtentX",(float)m_Surface->m_XExtent);
    m_Shader_Flow->m_program->setUniformValue("SurfaceExtentZ",(float)m_Surface->m_ZExtent);
    m_Shader_Flow->m_program->setUniformValue("ElevationMin",(float)m_Surface->m_ElevationMin);
    m_Shader_Flow->m_program->setUniformValue("ElevationMax",(float)m_Surface->m_ElevationMax);
    m_Shader_Flow->m_program->setUniformValue("viewportSize",camera->m_viewportSize);

    m_Shader_Flow->m_program->setUniformValue("CellSize",(float)m_Surface->m_CellSize);

    m_Surface->m_GLObject_Tesselated.bind();

    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture,"heightMap",0);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_MicroElevation,"microElevation",1);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_MicroElevation_Normal,"microElevation_normal",2);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_Mask,"mask",3);

    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_SlopeX,"slopeX",4);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_SlopeY,"slopeY",5);

    m_Shader_Flow->ActivateTextureOn(widget,m_Texture_FlowH,"FlowH",6);
    m_Shader_Flow->ActivateTextureOn(widget,m_Texture_FlowU,"FlowU",7);
    m_Shader_Flow->ActivateTextureOn(widget,m_Texture_FlowV,"FlowV",8);
    m_Shader_Flow->ActivateTextureOn(widget,m_Texture_FlowS,"FlowS",9);

    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_VegCover,"VegCover",10);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_VegHeight,"VegHeight",11);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_RandomRoughness,"RandomRoughness",12);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_Buildings,"Buildings",13);
    m_Shader_Flow->ActivateTextureOn(widget,m_Surface->m_Texture_Roads,"Roads",14);

    m_Shader_Flow->m_program->setUniformValue("TextureSizeX",(float)2.5);
    m_Shader_Flow->m_program->setUniformValue("TextureSizeY",(float)2.5);

    m_Shader_Flow->m_program->setPatchVertexCount(3);
    widget->gl->glDrawElements(GL_PATCHES, m_Surface->m_Geometry_Tesselated->m_IndexCount, GL_UNSIGNED_INT, 0);


    int e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() <<e;
    }

    m_Surface->m_GLObject_Tesselated.release();
    m_Shader_Flow->m_program->release();

    widget->gl->glEnable(GL_BLEND);

}
void GL3DFlowSurface::OnDestroy(GL3DWidget * widget)
{


}
