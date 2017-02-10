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

#ifndef Surface3D
#define Surface3D

#include "ui_full/3D/Objects/GL3DObject.h"
#include "ui_full/3D/GL3DWidget.h"
#include "ui_full/3D/Graphics/GL3DShaders.h"
#include "ui_full/3D/Graphics/GL3DTextures.h"

class GL3DSurface : public GL3DObject
{

public:
    GL3DSurface() : GL3DObject() {
    }

    double m_XResolution;
    double m_ZResolution;
    double m_XExtent;
    double m_ZExtent;
    double m_CellSize;
    double m_ElevationMin;
    double m_ElevationMax;


    cTMap * m_Elevation = 0;
    cTMap * m_ElevationFilled = 0;
    cTMap * m_SlopeX = 0;
    cTMap * m_SlopeY = 0;

    GL3DShader * m_Shader;
    GL3DShader * m_Shader_Tesselated;

    GL3DTexture * m_Texture;
    GL3DTexture * m_Texture_SlopeX;
    GL3DTexture * m_Texture_SlopeY;

    GL3DTexture * m_Texture_Mask;
    GL3DTexture * m_Texture_MicroElevation;
    GL3DTexture * m_Texture_MicroElevation_Normal;

    GL3DTexture * m_Texture_Grass_Color;
    GL3DTexture * m_Texture_Grass_Bump;
    GL3DTexture * m_Texture_Grass_Normal;
    GL3DTexture * m_Texture_Bare_Color;
    GL3DTexture * m_Texture_Bare_Bump;
    GL3DTexture * m_Texture_Bare_Normal;

    GL3DGeometry * m_Geometry;
    GL3DGeometry * m_Geometry_Tesselated;

    QOpenGLVertexArrayObject m_GLObject;
    QOpenGLVertexArrayObject m_GLObject_Tesselated;

    void SetSurfaceMap(cTMap * Elevation,cTMap * Elevationf,cTMap * sx,cTMap * sy);
    void OnCreate(GL3DWidget *widget);
    void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

    double GetElevation(double x, double z);

};

#endif
