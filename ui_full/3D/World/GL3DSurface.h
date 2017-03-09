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

#include "ui_full/3D/GL3DWidget.h"
#include "ui_full/3D/Graphics/GL3DShaders.h"
#include "ui_full/3D/Graphics/GL3DTextures.h"


typedef struct GLLDD_LINKEDLIST {
    int rowNr;
    int colNr;
    struct GLLDD_LINKEDLIST *prev;
}  GLLDD_LINKEDLIST;

class GL3DSurface
{

public:
    GL3DSurface()
    {
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

    GL3DTexture * m_Texture_ChannelDepth;

    GL3DTexture * m_Texture_Mask;
    GL3DTexture * m_Texture_MicroElevation;
    GL3DTexture * m_Texture_MicroElevation_Normal;

    GL3DTexture * m_Texture_DemChange;
    bool dem_updated = false;

    GL3DTexture * m_Texture_grass_slope0_Color;
    GL3DTexture * m_Texture_grass_slope0_Bump;
    GL3DTexture * m_Texture_grass_slope0_Normal;
    GL3DTexture * m_Texture_grass_slope0_Spec;

    GL3DTexture * m_Texture_grass_slope1_Color;
    GL3DTexture * m_Texture_grass_slope1_Bump;
    GL3DTexture * m_Texture_grass_slope1_Normal;
    GL3DTexture * m_Texture_grass_slope1_Spec;

    GL3DTexture * m_Texture_bare_slope0_Color;
    GL3DTexture * m_Texture_bare_slope0_Bump;
    GL3DTexture * m_Texture_bare_slope0_Normal;
    GL3DTexture * m_Texture_bare_slope0_Spec;

    GL3DTexture * m_Texture_bare_slope1_Color;
    GL3DTexture * m_Texture_bare_slope1_Bump;
    GL3DTexture * m_Texture_bare_slope1_Normal;
    GL3DTexture * m_Texture_bare_slope1_Spec;

    GL3DTexture * m_Texture_erosion_Color;
    GL3DTexture * m_Texture_erosion_Bump;
    GL3DTexture * m_Texture_erosion_Normal;
    GL3DTexture * m_Texture_erosion_Spec;

    GL3DTexture * m_Texture_deposition_Color;
    GL3DTexture * m_Texture_deposition_Bump;
    GL3DTexture * m_Texture_deposition_Normal;
    GL3DTexture * m_Texture_deposition_Spec;

    GL3DTexture * m_Texture_Channel;

    double m_Geometry_Micro_Width = 25;
    int m_Geometry_Micro_Cells = 3;
    int m_Geometry_Multiplier = 45.0;
    GL3DGeometry * m_Geometry;
    GL3DGeometry * m_Geometry_Tesselated;
    GL3DGeometry * m_Geometry_Micro;

    QOpenGLVertexArrayObject m_GLObject;
    QOpenGLVertexArrayObject m_GLObject_Tesselated;

    void SetSurfaceMap(cTMap * Elevation,cTMap * Elevationf,cTMap * sx,cTMap * sy);


    cTMap * m_VegCover = 0;
    cTMap * m_VegHeight = 0;
    cTMap * m_RandomRoughness = 0;
    cTMap * m_Buildings = 0;
    cTMap * m_Roads = 0;
    cTMap * m_DemChange = 0;

    cTMap * ChannelLDD;
    cTMap * ChannelWidth;
    cTMap * ChannelDepth;

    bool has_channel = false;
    bool has_channeldepth = false;

    GL3DGeometry * m_Geometry_Channel;
    QOpenGLVertexArrayObject * m_Object_Channel;
    GL3DShader * m_Shader_Channel;

    GL3DTexture * m_Texture_VegCover;
    GL3DTexture * m_Texture_VegHeight;
    GL3DTexture * m_Texture_RandomRoughness;
    GL3DTexture * m_Texture_Buildings;
    GL3DTexture * m_Texture_Roads;

    void SetTerrainProperties(cTMap * VegCover,cTMap * VegHeight,cTMap * RandomRoughness,cTMap * Buildings,cTMap * Roads);
    void SetDemChange(cTMap*DemChange);
    void RotateRight(int &dx, int &dy);
    void RotateLeft(int &dx, int &dy);
    int GetRotationType(int dx1,int dy1,int dx2,int dy2);

    void SetChannel(GL3DWidget * w,cTMap * LDD, cTMap * width, cTMap * depth, cTMap * temp, bool flooding);

    void OnCreate(GL3DWidget *widget);
    void OnRenderBefore(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

    double GetElevation(double x, double z);

};

#endif
