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

#include <3D/Objects/GL3DSurface.h>
#include <QtOpenGL>
#include <QColor>
#include <QGLWidget>
#include <QGLWidget>
#include <QOpenGLFunctions>
#include <QBasicTimer>
#include <QtGui>
#include <QKeyEvent>
#include <QMetaEnum>
#include <Qt>

void GL3DSurface::SetSurfaceMap(cTMap * Elevation,cTMap * Elevationf,cTMap * sx,cTMap * sy)
{
    this->m_Elevation = Elevation;
    this->m_ElevationFilled = Elevationf;
    this->m_SlopeX = sx;
    this->m_SlopeY = sy;


    this->m_XResolution = m_Elevation->nrCols();
    this->m_ZResolution = m_Elevation->nrRows();
    this->m_XExtent = m_Elevation->nrCols() * m_Elevation->cellSize();
    this->m_ZExtent = m_Elevation->nrRows() * m_Elevation->cellSize();
    this->m_CellSize = m_Elevation->cellSize();

    bool first = true;
    this->m_ElevationMin = 0;
    this->m_ElevationMax = 0;

    FOR_ROW_COL_MV(Elevation,Elevation)
    {
        if(first)
        {
            this->m_ElevationMin = Elevation->Drc;
            this->m_ElevationMax = Elevation->Drc;
            first = false;

        }else
        {
            this->m_ElevationMin = std::min(float(Elevation->Drc),float(this->m_ElevationMin));
            this->m_ElevationMax = std::max(float(Elevation->Drc),float(this->m_ElevationMax));
        }
    }
}

void GL3DSurface::OnCreate(GL3DWidget *widget)
{
    //this->m_Shader = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SIMPLE);
    this->recieve_render = true;
    this->m_Shader = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SIMPLE);
    this->m_Shader_Tesselated = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SURFACE_TESSELATED);


    this->m_Texture = widget->m_Textures->LoadTextureFromMap(false, this->m_ElevationFilled,0,0,true);
    this->m_Texture_Mask = widget->m_Textures->LoadTextureFromMap(false, this->m_Elevation,0,0,true,true);
    this->m_Texture_SlopeX = widget->m_Textures->LoadTextureFromMap(true, this->m_SlopeX,0,0,true);
    this->m_Texture_SlopeY = widget->m_Textures->LoadTextureFromMap(true, this->m_SlopeY,0,0,true);

    this->m_Texture_MicroElevation = widget->m_Textures->LoadTextureFromFile("micro-elevation.bmp");
    this->m_Texture_MicroElevation_Normal = widget->m_Textures->LoadTextureFromFile("micro-elevation_normal.bmp");

    this->m_Texture_Grass_Color = widget->m_Textures->LoadTextureFromFile("grass_color.bmp");
    this->m_Texture_Grass_Bump = widget->m_Textures->LoadTextureFromFile("grass_bmp.bmp");
    this->m_Texture_Grass_Normal = widget->m_Textures->LoadTextureFromFile("grass_normal.bmp");
    this->m_Texture_Bare_Color = widget->m_Textures->LoadTextureFromFile("bare_color.bmp",true,true);
    this->m_Texture_Bare_Bump = widget->m_Textures->LoadTextureFromFile("bare_bmp.bmp",true,true);
    this->m_Texture_Bare_Normal = widget->m_Textures->LoadTextureFromFile("bare_normal.bmp",true,true);

    if(this->m_Elevation != 0)
    {

        this->m_Geometry = widget->m_Geometries->LoadGeometryFromMap(this->m_Elevation,0,0);
        this->m_GLObject.create();
        widget->BindGeometry(m_GLObject,m_Shader,m_Geometry);

        this->m_Geometry_Tesselated = widget->m_Geometries->LoadGeometryFromMap(this->m_Elevation,0,0);
        this->m_GLObject_Tesselated.create();
        widget->BindGeometry(m_GLObject_Tesselated,m_Shader_Tesselated,m_Geometry_Tesselated);
    }
}

void GL3DSurface::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{

    QMatrix4x4 matrix;
    matrix.setToIdentity();
    m_Shader_Tesselated->m_program->bind();
    m_Shader_Tesselated->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    m_Shader_Tesselated->m_program->setUniformValue("Mmatrix",matrix);
    m_Shader_Tesselated->m_program->setUniformValue("Cposition",camera->m_Position);

    m_Shader_Tesselated->m_program->setUniformValue("microElevationScaleX",(float) this->m_CellSize);
    m_Shader_Tesselated->m_program->setUniformValue("microElevationScaleY",(float) 0.25f);
    m_Shader_Tesselated->m_program->setUniformValue("microElevationScaleZ",(float) this->m_CellSize);

    m_Shader_Tesselated->m_program->setUniformValue("SurfaceExtentX",(float)this->m_XExtent);
    m_Shader_Tesselated->m_program->setUniformValue("SurfaceExtentZ",(float)this->m_ZExtent);
    m_Shader_Tesselated->m_program->setUniformValue("ElevationMin",(float)this->m_ElevationMin);
    m_Shader_Tesselated->m_program->setUniformValue("ElevationMax",(float)this->m_ElevationMax);
    m_Shader_Tesselated->m_program->setUniformValue("viewportSize",camera->m_viewportSize);

    m_Shader_Tesselated->m_program->setUniformValue("CellSize",(float)this->m_CellSize);

    this->m_GLObject_Tesselated.bind();

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture,"heightMap",0);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_MicroElevation,"microElevation",1);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_MicroElevation_Normal,"microElevation_normal",2);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Mask,"mask",3);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_SlopeX,"slopeX",4);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_SlopeY,"slopeY",5);

    m_Shader_Tesselated->m_program->setUniformValue("TextureSizeX",(float)2.5);
    m_Shader_Tesselated->m_program->setUniformValue("TextureSizeY",(float)2.5);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Grass_Color,"grass_color",6);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Grass_Bump,"grass_bump",7);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Grass_Normal,"grass_normal",8);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Bare_Color,"gravel_color",9);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Bare_Bump,"gravel_bump",10);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Bare_Normal,"gravel_normal",11);

    //widget->gl->glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    m_Shader_Tesselated->m_program->setPatchVertexCount(3);
    widget->gl->glDrawElements(GL_PATCHES, this->m_Geometry_Tesselated->m_IndexCount, GL_UNSIGNED_INT, 0);

    //widget->gl->glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    int e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() <<e;
    }
    this->m_GLObject_Tesselated.release();
    m_Shader_Tesselated->m_program->release();

}

void GL3DSurface::OnDestroy(GL3DWidget *widget)
{


}


double GL3DSurface::GetElevation(double x, double z)
{
        double cs = m_ElevationFilled->cellSize();
        int cn = std::floor(x/cs);
        int rn = std::floor(z/cs);

        double wx = 1.0-(x - (cn * cs))/cs;
        double wy = 1.0-(z - (rn * cs))/cs;

        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
        double v4 = 0;

        double w1 = 0;
        double w2 = 0;
        double w3 = 0;
        double w4 = 0;

        if(!OUTOFMAP(m_ElevationFilled,rn,cn))
        {
            if(!MV(m_ElevationFilled,rn,cn))
            {
                w1 = (1.0-wx) * (1.0-wy);
                v1 =m_ElevationFilled->data[rn][cn];
            }
        }
        if(!OUTOFMAP(m_ElevationFilled,rn-1,cn))
        {
            if(!MV(m_ElevationFilled,rn-1,cn))
            {
                w2 = (1.0-wx) * (wy);
                v2 =m_ElevationFilled->data[rn-1][cn];
            }
        }
        if(!OUTOFMAP(m_ElevationFilled,rn,cn-1))
        {
            if(!MV(m_ElevationFilled,rn,cn-1))
            {
                w3 = (wx) * (1.0-wy);
                v3 =m_ElevationFilled->data[rn][cn-1];
            }
        }
        if(!OUTOFMAP(m_ElevationFilled,rn-1,cn-1))
        {
            if(!MV(m_ElevationFilled,rn-1,cn-1))
            {
                w4 = (wx) * (wy);
                v4 =m_ElevationFilled->data[rn-1][cn-1];
            }
        }

        return (w1+w2+w3+w4) > 0? ((w1*v1 + w2*v2 + w3*v3 + w4 * v4)/(w1+w2+w3+w4)): 0.0;



}
