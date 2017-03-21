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

#include <3D/World/GL3DSurface.h>
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

void GL3DSurface::SetTerrainProperties(cTMap * VegCover,cTMap * VegHeight,cTMap * RandomRoughness,cTMap * Buildings,cTMap * Roads)
{

    m_VegCover = VegCover;
    m_VegHeight = VegHeight;
    m_RandomRoughness = RandomRoughness;
    m_Buildings = Buildings;
    m_Roads = Roads;

}

void GL3DSurface::SetDemChange(cTMap*DemChange)
{

    dem_updated = true;
    m_DemChange = DemChange;

}

void GL3DSurface::RotateRight(int &dx, int &dy)
{
    int dxo = dx;
    int dyo = dy;

    if(dxo == 1 && dyo == 1)
    {
        dx = 1;
        dy = 0;

    }else if(dxo ==0 && dyo == 1)
    {
        dx = 1;
        dy = 1;

    }else if(dxo == 1 && dyo == 0)
    {
        dx = 1;
        dy = -1;

    }else if(dxo == -1 && dyo == -1)
    {
        dx = -1;
        dy = 0;

    }else if(dxo == -1 && dyo == 0)
    {
        dx = -1;
        dy = 1;

    }else if(dxo == 0 && dyo == -1)
    {
        dx = -1;
        dy = -1;

    }else if(dxo == -1 && dyo == 0)
    {
        dx = 0;
        dy = 1;

    }else if(dxo == 1 && dyo == -1)
    {
        dx = 0;
        dy = -1;
    }
}

void GL3DSurface::RotateLeft(int &dx, int &dy)
{
    int dxo = dx;
    int dyo = dy;

    if(dxo == 1 && dyo == 0)
    {
        dx = 1;
        dy = 1;

    }else if(dxo == 1 && dyo == 1)
    {
        dx = 0;
        dy = 1;

    }else if(dxo == 1 && dyo == -1)
    {
        dx = 1;
        dy = 0;

    }else if(dxo == -1 && dyo == 0)
    {
        dx = -1;
        dy = -1;

    }else if(dxo == -1 && dyo == 1)
    {
        dx = -1;
        dy = 0;

    }else if(dxo == -1 && dyo == -1)
    {
        dx = 0;
        dy = -1;

    }else if(dxo == 0 && dyo == 1)
    {
        dx = -1;
        dy = 0;

    }else if(dxo == 0 && dyo == -1)
    {
        dx = 1;
        dy = -1;
    }

}

int GL3DSurface::GetRotationType(int dx1, int dy1, int dx2, int dy2)
{
    int z_cross = dx1 * dy2 - dy1 * dx2;
    if(z_cross == 0)
    {
        return 0;
    }else if(z_cross > 0)
    {
        return 1;
    }else if(z_cross < 0)
    {
        return -1;
    }
    return 0;
}

#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

/// shortcut to check if r,c is inside map boundaries, used in kinematic and flooding
#define INSIDE(r, c) (r>=0 && r<_nrRows && c>=0 && c<_nrCols)

void GL3DSurface::SetChannel(GL3DWidget * widget,cTMap * LDD, cTMap * width, cTMap * depth, cTMap * temp, bool flooding)
{
    has_channel = true;
    has_channeldepth = flooding;

    ChannelLDD = LDD;
    ChannelWidth = width;
    ChannelDepth = depth;

    cTMap * tmb= temp;

    int _nrRows = LDD->nrRows();
    int _nrCols = LDD->nrCols();

    double cs = LDD->cellSize();


    //set mask for cells that are already done
    int n_channel = 0;
    tmb->setAllMV();

    //count total channel cells
    for (int  ro = 0; ro < _nrRows; ro++){
    for (int  co = 0; co < _nrCols; co++){
    if(!pcr::isMV(LDD->data[ro][co]))
    {
        if(LDD->data[ro][co] != 5)
        {
            n_channel++;
        }

    }}}


    //allocate data
    Vertex * vertices = (Vertex *)malloc(sizeof(Vertex) * (n_channel * 12));
    GLuint * indices = (GLuint *)malloc(sizeof(GLuint)*(  n_channel * 18));

    int n_current = 0;

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    for (int  ro = 0; ro < _nrRows; ro++){
    for (int  co = 0; co < _nrCols; co++){
    if(!pcr::isMV(LDD->data[ro][co]))
    {
        if(LDD->data[ro][co] != 5)
        {
                    int r = ro;
                    int c = co;
                    tmb->Drc = 0;

                    int rp = r;
                    int cp = c;
                    int np = 0;

                    int i = 0;
                    for (i=1; i<=9; i++)
                    {
                        int rt, ct;
                        int ldd = 0;

                        // this is the current cell
                        if (i==5)
                            continue;

                        rt = r+dy[i];
                        ct = c+dx[i];

                        if (INSIDE(rt, ct) && !pcr::isMV(LDD->data[rt][ct]))
                            ldd = (int) LDD->data[rt][ct];
                        else
                            continue;

                        // check if there are more cells upstream, if not subCatchDone remains true
                        if (FLOWS_TO(ldd, rt, ct, r, c) &&
                                INSIDE(rt, ct))
                        {
                            rp = rt;
                            cp = ct;
                            np ++;
                        }
                    }

                    if(np > 1)
                    {
                        rp = r;
                        cp = c;
                    }

                    //subcache done

                    int ldd_this = LDD->Drc;
                    if(ldd_this != 5)
                    {
                        int rn = r+dy[ldd_this];
                        int cn = c+dx[ldd_this];

                        if(INSIDE(rn,cn))
                        {
                            int ldd_next = LDD->data[rn][cn];
                            int ldd_prev = LDD->data[rp][cp];
                            float depth_this = abs( depth->Drc);
                            float depth_next = abs( depth->data[rn][cn]);
                            float depth_prev = abs( depth->data[rp][cp]);
                            float width_this = width->Drc;
                            float width_next = width->data[rn][cn];
                            float width_prev = width->data[rp][cp];

                            int dx1 = dx[ldd_this];
                            int dy1 = dy[ldd_this];

                            int dx2 = dx[ldd_next];
                            int dy2 = dy[ldd_next];

                            int rot_type1 = 0;//GetRotationType(dx[ldd_this],dy[ldd_this],dx[ldd_prev],dy[ldd_prev]);
                            int rot_type2 = 0;//GetRotationType(dx[ldd_next],dy[ldd_next],dx[ldd_this],dy[ldd_this]);

                            //if rot_type = 0, oppsite of direction is with vector

                            int dxw1 = dx1;
                            int dyw1 = dy1;
                            int dxw2 = dx2;
                            int dyw2 = dy2;

                            if(rot_type1 > 0)
                            {
                                RotateRight(dxw1,dyw1);
                            }else
                            {
                                RotateLeft(dxw1,dyw1);
                            }
                            if(rot_type2 > 0)
                            {
                                RotateRight(dxw2,dyw2);
                            }else
                            {
                                RotateLeft(dxw2,dyw2);
                            }

                            float xl1 = ((double)(c) + 0.5) * cs + 0.5 * width_this * (rot_type1 == 0? (-dy[ldd_this]) : dyw1);
                            float yl1 = ((double)(r) + 0.5) * cs + 0.5 * width_this * (rot_type1 == 0? (dx[ldd_this]) : dxw1);
                            float zl1 = this->GetElevation(xl1,yl1);
                            float xr1 = ((double)(c) + 0.5) * cs - 0.5 * width_this * (rot_type1 == 0? (-dy[ldd_this]) : dyw1);
                            float yr1 = ((double)(r) + 0.5) * cs - 0.5 * width_this * (rot_type1 == 0? (dx[ldd_this]) : dxw1);
                            float zr1 = this->GetElevation(xr1,yr1);

                            float xl2 = ((double)(cn) + 0.5) * cs + 0.5 * width_next * (rot_type2 == 0? (-dy[ldd_next]) : dyw2);
                            float yl2 = ((double)(rn) + 0.5) * cs + 0.5 * width_next * (rot_type2 == 0? (dx[ldd_next]) : dxw2);
                            float zl2 = this->GetElevation(xl2,yl2);
                            float xr2 = ((double)(cn) + 0.5) * cs - 0.5 * width_next * (rot_type2 == 0? (-dy[ldd_next]) : dyw2);
                            float yr2 = ((double)(rn) + 0.5) * cs - 0.5 * width_next * (rot_type2 == 0? (dx[ldd_next]) : dxw2);
                            float zr2 = this->GetElevation(xr2,yr2);

                            vertices[n_current*12 + 0] = Vertex(QVector3D(xl1,zl1,yl1));
                            vertices[n_current*12 + 1] = Vertex(QVector3D(xl2,zl2,yl2));
                            vertices[n_current*12 + 2] = Vertex(QVector3D(xl1,zl1-depth_this,yl1));
                            vertices[n_current*12 + 3] = Vertex(QVector3D(xl2,zl2-depth_next,yl2));
                            vertices[n_current*12 + 4] = Vertex(QVector3D(xl1,zl1-depth_this,yl1));
                            vertices[n_current*12 + 5] = Vertex(QVector3D(xr1,zr1-depth_this,yr1));
                            vertices[n_current*12 + 6] = Vertex(QVector3D(xl2,zl2-depth_next,yl2));
                            vertices[n_current*12 + 7] = Vertex(QVector3D(xr2,zr2-depth_next,yr2));
                            vertices[n_current*12 + 8] = Vertex(QVector3D(xr1,zr1,yr1));
                            vertices[n_current*12 + 9] = Vertex(QVector3D(xr2,zr2,yr2));
                            vertices[n_current*12 + 10] = Vertex(QVector3D(xr1,zr1-depth_this,yr1));
                            vertices[n_current*12 + 11] = Vertex(QVector3D(xr2,zr2-depth_next,yr2));

                            indices[n_current*18 + 0] = n_current*12 + 0;
                            indices[n_current*18 + 1] = n_current*12 + 2;
                            indices[n_current*18 + 2] = n_current*12 + 1;
                            indices[n_current*18 + 3] = n_current*12 + 1;
                            indices[n_current*18 + 4] = n_current*12 + 2;
                            indices[n_current*18 + 5] = n_current*12 + 3;
                            indices[n_current*18 + 6] = n_current*12 + 4;
                            indices[n_current*18 + 7] = n_current*12 + 5;
                            indices[n_current*18 + 8] = n_current*12 + 6;
                            indices[n_current*18 + 9] = n_current*12 + 5;
                            indices[n_current*18 + 10] = n_current*12 + 7;
                            indices[n_current*18 + 11] = n_current*12 + 6;
                            indices[n_current*18 + 12] = n_current*12 + 8;
                            indices[n_current*18 + 13] = n_current*12 + 9;
                            indices[n_current*18 + 14] = n_current*12 + 10;
                            indices[n_current*18 + 15] = n_current*12 + 9;
                            indices[n_current*18 + 16] = n_current*12 + 11;
                            indices[n_current*18 + 17] = n_current*12 + 10;

                            //draw these sides transparantly


                        }

                        n_current ++;
                    }

        }
    }}}

    qDebug() << "channel created Nr: " << n_channel << n_current;


    m_Geometry_Channel = new GL3DGeometry();
    m_Geometry_Channel->CreateGeometry(widget,vertices,12 * n_channel,indices,n_channel*18);
    m_Object_Channel = new QOpenGLVertexArrayObject();
    m_Object_Channel->create();

    m_Shader_Channel = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_CHANNEL);
    GL3DDrawFunctions::BindGeometry(widget,*m_Object_Channel,m_Shader_Channel,m_Geometry_Channel);

}

void GL3DSurface::CreateMicroSurface(GL3DWidget *widget)
{

}


void GL3DSurface::OnCreate(GL3DWidget *widget)
{
    //this->m_Shader = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SIMPLE);
    this->m_Shader = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SIMPLE);
    this->m_Shader_Tesselated = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SURFACE_TESSELATED);

    if(has_channel)
    {
        m_Texture_ChannelDepth = widget->m_Textures->LoadTextureFromMap(false, this->ChannelDepth,0,0,true);
    }

    this->m_Texture = widget->m_Textures->LoadTextureFromMap(false, this->m_ElevationFilled,0,0,true);
    this->m_Texture_Mask = widget->m_Textures->LoadTextureFromMap(false, this->m_Elevation,0,0,true,true);
    this->m_Texture_SlopeX = widget->m_Textures->LoadTextureFromMap(false, this->m_SlopeX,0,0,true);
    this->m_Texture_SlopeY = widget->m_Textures->LoadTextureFromMap(false, this->m_SlopeY,0,0,true);

    this->m_Texture_MicroElevation = widget->m_Textures->LoadTextureFromFile("micro-elevation.bmp");
    this->m_Texture_MicroElevation_Normal = widget->m_Textures->LoadTextureFromFile("micro-elevation_normal.bmp");

    this->m_Texture_DemChange = widget->m_Textures->LoadTextureFromMap(false, this->m_DemChange,0,0,true);

    m_Texture_grass_slope0_Color = widget->m_Textures->LoadTextureFromFile("grass_slope0_color.bmp");
    m_Texture_grass_slope0_Bump = widget->m_Textures->LoadTextureFromFile("grass_slope0_bmp.bmp");
    m_Texture_grass_slope0_Normal = widget->m_Textures->LoadTextureFromFile("grass_slope0_normal.bmp");
    m_Texture_grass_slope0_Spec = widget->m_Textures->LoadTextureFromFile("grass_slope0_spec.bmp");

    m_Texture_grass_slope1_Color = widget->m_Textures->LoadTextureFromFile("grass_slope3_color.bmp");
    m_Texture_grass_slope1_Bump = widget->m_Textures->LoadTextureFromFile("grass_slope3_bmp.bmp");
    m_Texture_grass_slope1_Normal = widget->m_Textures->LoadTextureFromFile("grass_slope3_normal.bmp");
    m_Texture_grass_slope1_Spec = widget->m_Textures->LoadTextureFromFile("grass_slope3_spec.bmp");

    m_Texture_bare_slope0_Color = widget->m_Textures->LoadTextureFromFile("bare_slope0_color.bmp");
    m_Texture_bare_slope0_Bump = widget->m_Textures->LoadTextureFromFile("bare_slope0_bmp.bmp");
    m_Texture_bare_slope0_Normal = widget->m_Textures->LoadTextureFromFile("bare_slope0_normal.bmp");
    m_Texture_bare_slope0_Spec = widget->m_Textures->LoadTextureFromFile("bare_slope0_spec.bmp");

    m_Texture_bare_slope1_Color = widget->m_Textures->LoadTextureFromFile("bare_slope1_color.bmp");
    m_Texture_bare_slope1_Bump = widget->m_Textures->LoadTextureFromFile("bare_slope1_bmp.bmp");
    m_Texture_bare_slope1_Normal = widget->m_Textures->LoadTextureFromFile("bare_slope1_normal.bmp");
    m_Texture_bare_slope1_Spec = widget->m_Textures->LoadTextureFromFile("bare_slope1_spec.bmp");

    m_Texture_erosion_Color = widget->m_Textures->LoadTextureFromFile("dep2_color.bmp");
    m_Texture_erosion_Bump = widget->m_Textures->LoadTextureFromFile("dep2_color.bmp");
    m_Texture_erosion_Normal = widget->m_Textures->LoadTextureFromFile("dep2_color.bmp");
    m_Texture_erosion_Spec = widget->m_Textures->LoadTextureFromFile("dep2_color.bmp");

    m_Texture_deposition_Color = widget->m_Textures->LoadTextureFromFile("dep1_color.bmp");
    m_Texture_deposition_Bump = widget->m_Textures->LoadTextureFromFile("dep1_color.bmp");
    m_Texture_deposition_Normal = widget->m_Textures->LoadTextureFromFile("dep1_color.bmp");
    m_Texture_deposition_Spec = widget->m_Textures->LoadTextureFromFile("dep1_color.bmp");

    m_Texture_Channel = widget->m_Textures->LoadTextureFromFile("channel.jpg");

    this->m_Texture_VegCover = widget->m_Textures->LoadTextureFromMap(false, this->m_VegCover,0,0,true);
    this->m_Texture_VegHeight = widget->m_Textures->LoadTextureFromMap(false, this->m_VegHeight,0,0,true);
    this->m_Texture_RandomRoughness = widget->m_Textures->LoadTextureFromMap(false, this->m_RandomRoughness,0,0,true);
    this->m_Texture_Buildings = widget->m_Textures->LoadTextureFromMap(false, this->m_Buildings,0,0,true);
    this->m_Texture_Roads = widget->m_Textures->LoadTextureFromMap(false, this->m_Roads,0,0,true);


    if(this->m_Elevation != 0)
    {
        this->has_microsurface = true;

        if(this->m_Elevation->cellSize() < 0.2 || m_XExtent * m_ZExtent < 100)
        {
            this->has_microsurface = false;
            return;
        }

        this->m_Geometry_MicroTesselated = widget->m_Geometries->LoadGeometryRaster(this->m_MicroPatches,this->m_MicroPatches,this->m_MicroPatchSize);
        this->m_GLObject_MicroTesselated.create();
        GL3DDrawFunctions::BindGeometry(widget,m_GLObject_MicroTesselated,m_Shader_Tesselated,m_Geometry_MicroTesselated);

        /*this->m_Geometry = widget->m_Geometries->LoadGeometryFromMap(this->m_Elevation,30);
        this->m_GLObject.create();
        GL3DDrawFunctions::BindGeometry(widget,m_GLObject,m_Shader,m_Geometry);*/

        this->m_Geometry_Tesselated = widget->m_Geometries->LoadGeometryFromMap(this->m_ElevationFilled,(float) m_Geometry_Multiplier);
        this->m_GLObject_Tesselated.create();
        GL3DDrawFunctions::BindGeometry(widget,m_GLObject_Tesselated,m_Shader_Tesselated,m_Geometry_Tesselated);
    }
}

void GL3DSurface::OnRenderBefore(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{

    if(dem_updated)
    {
        qDebug() << "GL update dem";
        FOR_ROW_COL_MV(m_Elevation,m_Elevation)
        {
            m_ElevationFilled->Drc = m_Elevation->Drc + m_DemChange->Drc;
        }
        this->m_Texture = widget->m_Textures->LoadTextureFromMap(false, this->m_ElevationFilled,0,0,true);
        this->m_Texture_DemChange = widget->m_Textures->LoadTextureFromMap(false, this->m_DemChange,0,0,true);
        dem_updated = false;
    }

    if(has_channel)
    {

        QMatrix4x4 matrix;
        matrix.setToIdentity();

        //Draw channel over surface

        //widget->gl->glDisable(GL_DEPTH_TEST);
        widget->gl->glDisable(GL_CULL_FACE);

        matrix.setToIdentity();

        //widget->DrawModelObject(m_Model,camera,matrix);

        m_Shader_Channel->m_program->bind();
        m_Shader_Channel->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
        m_Shader_Channel->m_program->setUniformValue("Mmatrix",matrix);
        m_Shader_Channel->m_program->setUniformValue("CPosition",camera->m_Position);

        m_Shader_Channel->m_program->setUniformValue("SurfaceExtentX",(float)m_XExtent);
        m_Shader_Channel->m_program->setUniformValue("SurfaceExtentZ",(float)m_ZExtent);

        m_Shader_Channel->m_program->setUniformValue("CellSize",(float)m_CellSize *(float) m_Geometry_Multiplier);

        m_Shader_Channel->m_program->setUniformValue("Light_Ambient",world->Light_Ambient);
        m_Shader_Channel->m_program->setUniformValue("Light_Directional",world->Light_Directional);
        m_Shader_Channel->m_program->setUniformValue("Light_Directional_Direction",world->Light_Directional_Direction);

        //m_Model_Channel->m_Shader->ActivateTextureOn(widget,m_Texture_Color,"Texture_Color",4);
        //m_Model_Channel->m_Shader->ActivateTextureOn(widget,m_Texture_Normal,"Texture_Normal",5);
        //m_Model_Channel->m_Shader->ActivateTextureOn(widget,m_Texture_Bump,"Texture_Bump",6);
        //m_Model_Channel->m_Shader->ActivateTextureOn(widget,m_Texture_Specular,"Texture_Specular",7);

        m_Shader_Channel->ActivateTextureOn(widget,m_Texture,"heightMap",0);
        m_Shader_Channel->ActivateTextureOn(widget,m_Texture_Mask,"mask",1);

        m_Shader_Channel->ActivateTextureOn(widget,m_Texture_SlopeX,"slopeX",2);
        m_Shader_Channel->ActivateTextureOn(widget,m_Texture_SlopeY,"slopeY",3);

        m_Shader_Channel->ActivateTextureOn(widget,m_Texture_Channel,"Texture_Color",4);

        m_Object_Channel->bind();

        widget->gl->glDrawElements(GL_TRIANGLES, m_Geometry_Channel->m_IndexCount, GL_UNSIGNED_INT, 0);
        m_Object_Channel->release();

        int e = widget->gl->glGetError();
        if( e != GL_NO_ERROR)
        {
            qDebug() << "opengl error in drawing roads " <<e;
        }


        m_Shader_Channel->m_program->release();

        //widget->gl->glEnable(GL_DEPTH_TEST);
        widget->gl->glEnable(GL_CULL_FACE);
    }


}

void GL3DSurface::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{

    if(dem_updated)
    {

        this->m_Texture_DemChange = widget->m_Textures->LoadTextureFromMap(false, this->m_DemChange,0,0,true);

        dem_updated= false;
    }

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

    m_Shader_Tesselated->m_program->setUniformValue("CellSize",(float)this->m_CellSize*(float) m_Geometry_Multiplier);

    m_Shader_Tesselated->m_program->setUniformValue("Light_Ambient",world->Light_Ambient);
    m_Shader_Tesselated->m_program->setUniformValue("Light_Directional",world->Light_Directional);
    m_Shader_Tesselated->m_program->setUniformValue("Light_Directional_Direction",world->Light_Directional_Direction);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture,"heightMap",0);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_MicroElevation,"microElevation",1);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_MicroElevation_Normal,"microElevation_normal",2);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Mask,"mask",3);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_SlopeX,"slopeX",4);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_SlopeY,"slopeY",5);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_VegCover,"VegCover",6);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_VegHeight,"VegHeight",7);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_RandomRoughness,"RandomRoughness",8);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Buildings,"Buildings",9);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Roads,"Roads",10);

    m_Shader_Tesselated->m_program->setUniformValue("TextureSizeX",(float)2.5);
    m_Shader_Tesselated->m_program->setUniformValue("TextureSizeY",(float)2.5);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope0_Color,"grass_slope0_Color",11);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope0_Bump,"grass_slope0_Bump",12);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope0_Normal,"grass_slope0_Normal",13);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope0_Spec,"grass_slope0_Spec",14);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope0_Color,"bare_slope0_Color",15);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope0_Bump,"bare_slope0_Bump",16);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope0_Normal,"bare_slope0_Normal",17);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope0_Spec,"bare_slope0_Spec",18);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Color,"grass_slope1_Color",19);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Bump,"grass_slope1_Bump",20);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Normal,"grass_slope1_Normal",21);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Spec,"grass_slope1_Spec",22);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Color,"bare_slope1_Color",23);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Bump,"bare_slope1_Bump",24);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Normal,"bare_slope1_Normal",25);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Spec,"bare_slope1_Spec",26);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_erosion_Color,"erosion_Color",19);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_erosion_Bump,"erosion_Bump",20);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_erosion_Normal,"erosion_Normal",21);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_erosion_Spec,"erosion_Spec",22);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_deposition_Color,"deposition_Color",23);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_deposition_Bump,"deposition_Bump",24);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_deposition_Normal,"deposition_Normal",25);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_deposition_Spec,"deposition_Spec",26);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_DemChange,"demChange",27);




    //draw total surface
    m_Shader_Tesselated->m_program->setUniformValue("fix_tesselation",false);
    m_Shader_Tesselated->m_program->setUniformValue("fix_tesselationlevel",1.0f);

    m_Shader_Tesselated->m_program->setUniformValue("TerrainOffset",QVector3D(0.0,0.0,0.0));


    this->m_GLObject_Tesselated.bind();

    m_Shader_Tesselated->m_program->setPatchVertexCount(3);
    widget->gl->glDrawElements(GL_PATCHES, this->m_Geometry_Tesselated->m_IndexCount, GL_UNSIGNED_INT, 0);

    this->m_GLObject_Tesselated.release();

    //draw microsurface
    m_Shader_Tesselated->m_program->setUniformValue("fix_tesselation",true);
    m_Shader_Tesselated->m_program->setUniformValue("fix_tesselationlevel",64.0f);

    m_Shader_Tesselated->m_program->setUniformValue("TerrainOffset",camera->m_Position - QVector3D(0.5 * float(m_MicroPatches) * m_MicroPatchSize,0.0,0.5 * float(m_MicroPatches) * m_MicroPatchSize));

    this->m_GLObject_MicroTesselated.bind();

    m_Shader_Tesselated->m_program->setPatchVertexCount(3);
    widget->gl->glDrawElements(GL_PATCHES, this->m_Geometry_MicroTesselated->m_IndexCount, GL_UNSIGNED_INT, 0);

    this->m_GLObject_MicroTesselated.release();


    //widget->gl->glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    int e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() <<e;
    }

    m_Shader_Tesselated->m_program->release();




}

void GL3DSurface::OnDestroy(GL3DWidget *widget)
{


}


double GL3DSurface::GetElevation(double x, double z)
{
    return MapMath::GetValueAt(m_ElevationFilled,x,z);

}
