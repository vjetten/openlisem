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

#ifndef FLOWSURFACE3D
#define FLOWSURFACE3D

#include <3D/GL3DWidget.h>
#include <3D/Objects/GL3DSurface.h>

class GL3DFlowSurface : public GL3DObject
{

public:
    GL3DFlowSurface() : GL3DObject() {

    }


    GL3DSurface * m_Surface;

    cTMap * m_FlowH =0;
    cTMap * m_FlowU =0;
    cTMap * m_FlowV =0;
    cTMap * m_FlowS =0;

    GL3DTexture * m_Texture_FlowH;
    GL3DTexture * m_Texture_FlowU;
    GL3DTexture * m_Texture_FlowV;
    GL3DTexture * m_Texture_FlowS;

    GL3DShader * m_Shader_Flow;

    bool text_created = false;
    bool text_updated = false;

    void SetSurface(GL3DSurface * surface);
    void SetFlowProperties(cTMap * h, cTMap * u, cTMap * v, cTMap * s);

    void CreateTextures(GL3DWidget * widget);

    void OnCreate(GL3DWidget * widget);
    void OnRenderLate(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget * widget);


};


#endif
