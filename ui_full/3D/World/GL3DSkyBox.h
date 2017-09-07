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

#ifndef SKYBOX3D
#define SKYBOX3D

#include "ui_full/3D/GL3DWidget.h"
#include "ui_full/3D/Graphics/GL3DShaders.h"
#include "ui_full/3D/Graphics/GL3DTextures.h"

class GL3DSkyBox
{

public:
    GL3DSkyBox() {
    }


    GL3DGeometry * m_Geometry;
    GL3DShader * m_Shader;
    GL3DTexture * m_Texture;

    QOpenGLVertexArrayObject m_GLObject;

    void OnCreate(GL3DWidget *widget);
    void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

};


#endif