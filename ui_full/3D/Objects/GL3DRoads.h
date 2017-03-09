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

#ifndef GL3DROADS_H
#define GL3DROADS_H

#include <3D/GL3DWidget.h>
#include <3D/Graphics/GL3DModels.h>
#include <3D/World/GL3DSurface.h>

#include <3D/Graphics/GL3DTextures.h>
#include <3D/Graphics/GL3DGeometry.h>
#include <3D/Graphics/GL3DShaders.h>

class GL3DRoads : public GL3DObject
{
public:
    GL3DRoads() : GL3DObject()
    {
    }

    bool draw = true;

    inline void SetDraw(bool in_draw)
    {
        draw = in_draw;
    }

    GL3DTexture * m_Texture_Color;
    GL3DTexture * m_Texture_Normal;
    GL3DTexture * m_Texture_Bump;
    GL3DTexture * m_Texture_Specular;

    GL3DGeometry * m_Geometry;

    QOpenGLVertexArrayObject m_GLObject;

    GL3DModel * m_Model;

    GL3DShader * m_Shader;

    GL3DSurface * m_Surface;

    void OnCreate(GL3DWidget *widget);
    void OnRenderBefore(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

    void SetRoadDistribution(GL3DWidget * w,GL3DSurface * s,cTMap * roadwidth);





};



#endif // GL3DROADS_H
