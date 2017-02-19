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

void GL3DSurface::SetTerrainProperties(cTMap * VegCover,cTMap * VegHeight,cTMap * RandomRoughness,cTMap * Buildings,cTMap * Roads)
{

    m_VegCover = VegCover;
    m_VegHeight = VegHeight;
    m_RandomRoughness = RandomRoughness;
    m_Buildings = Buildings;
    m_Roads = Roads;

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
    widget->BindGeometry(*m_Object_Channel,m_Shader_Channel,m_Geometry_Channel);

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

    m_Texture_grass_bare0_Color = widget->m_Textures->LoadTextureFromFile("grass_slope0_color.bmp");
    m_Texture_grass_bare0_Bump = widget->m_Textures->LoadTextureFromFile("grass_slope0_bmp.bmp");
    m_Texture_grass_bare0_Normal = widget->m_Textures->LoadTextureFromFile("grass_slope0_normal.bmp");
    m_Texture_grass_bare0_Spec = widget->m_Textures->LoadTextureFromFile("grass_slope0_spec.bmp");

    m_Texture_grass_bare1_Color = widget->m_Textures->LoadTextureFromFile("grass_bare1_color.bmp");
    m_Texture_grass_bare1_Bump = widget->m_Textures->LoadTextureFromFile("grass_bare1_bmp.bmp");
    m_Texture_grass_bare1_Normal = widget->m_Textures->LoadTextureFromFile("grass_bare1_normal.bmp");
    m_Texture_grass_bare1_Spec = widget->m_Textures->LoadTextureFromFile("grass_bare1_spec.bmp");

    m_Texture_grass_bare2_Color = widget->m_Textures->LoadTextureFromFile("grass_bare2_color.bmp");
    m_Texture_grass_bare2_Bump = widget->m_Textures->LoadTextureFromFile("grass_bare2_bmp.bmp");
    m_Texture_grass_bare2_Normal = widget->m_Textures->LoadTextureFromFile("grass_bare2_normal.bmp");
    m_Texture_grass_bare2_Spec = widget->m_Textures->LoadTextureFromFile("grass_bare2_spec.bmp");

    m_Texture_grass_bare3_Color = widget->m_Textures->LoadTextureFromFile("grass_bare3_color.bmp");
    m_Texture_grass_bare3_Bump = widget->m_Textures->LoadTextureFromFile("grass_bare3_bmp.bmp");
    m_Texture_grass_bare3_Normal = widget->m_Textures->LoadTextureFromFile("grass_bare3_normal.bmp");
    m_Texture_grass_bare3_Spec = widget->m_Textures->LoadTextureFromFile("grass_bare3_spec.bmp");

    m_Texture_bare_slope1_Color = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope1_color.bmp");
    m_Texture_bare_slope1_Bump = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope1_bmp.bmp");
    m_Texture_bare_slope1_Normal = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope1_normal.bmp");
    m_Texture_bare_slope1_Spec = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope1_spec.bmp");

    m_Texture_bare_slope4_Color = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope4_color.bmp");
    m_Texture_bare_slope4_Bump = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope4_bmp.bmp");
    m_Texture_bare_slope4_Normal = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope4_normal.bmp");
    m_Texture_bare_slope4_Spec = widget->m_Textures->LoadTextureFromFile("grass_bare4_slope4_spec.bmp");

    m_Texture_grass_slope1_Color = widget->m_Textures->LoadTextureFromFile("grass_slope1_color.bmp");
    m_Texture_grass_slope1_Bump = widget->m_Textures->LoadTextureFromFile("grass_slope1_bmp.bmp");
    m_Texture_grass_slope1_Normal = widget->m_Textures->LoadTextureFromFile("grass_slope1_normal.bmp");
    m_Texture_grass_slope1_Spec = widget->m_Textures->LoadTextureFromFile("grass_slope1_spec.bmp");

    m_Texture_grass_slope3_Color = widget->m_Textures->LoadTextureFromFile("grass_slope3_color.bmp");
    m_Texture_grass_slope3_Bump = widget->m_Textures->LoadTextureFromFile("grass_slope3_bmp.bmp");
    m_Texture_grass_slope3_Normal = widget->m_Textures->LoadTextureFromFile("grass_slope3_normal.bmp");
    m_Texture_grass_slope3_Spec = widget->m_Textures->LoadTextureFromFile("grass_slope3_spec.bmp");

    m_Texture_grass_slope4_Color = widget->m_Textures->LoadTextureFromFile("grass_slope4_color.bmp");
    m_Texture_grass_slope4_Bump = widget->m_Textures->LoadTextureFromFile("grass_slope4_bmp.bmp");
    m_Texture_grass_slope4_Normal = widget->m_Textures->LoadTextureFromFile("grass_slope4_normal.bmp");
    m_Texture_grass_slope4_Spec = widget->m_Textures->LoadTextureFromFile("grass_slope4_spec.bmp");

    m_Texture_Channel = widget->m_Textures->LoadTextureFromFile("channel.jpg");

    this->m_Texture_VegCover = widget->m_Textures->LoadTextureFromMap(false, this->m_VegCover,0,0,true);
    this->m_Texture_VegHeight = widget->m_Textures->LoadTextureFromMap(false, this->m_VegHeight,0,0,true);
    this->m_Texture_RandomRoughness = widget->m_Textures->LoadTextureFromMap(false, this->m_RandomRoughness,0,0,true);
    this->m_Texture_Buildings = widget->m_Textures->LoadTextureFromMap(false, this->m_Buildings,0,0,true);
    this->m_Texture_Roads = widget->m_Textures->LoadTextureFromMap(false, this->m_Roads,0,0,true);

    if(this->m_Elevation != 0)
    {

        this->m_Geometry = widget->m_Geometries->LoadGeometryFromMap(this->m_Elevation,0,0);
        this->m_GLObject.create();
        widget->BindGeometry(m_GLObject,m_Shader,m_Geometry);

        this->m_Geometry_Tesselated = widget->m_Geometries->LoadGeometryFromMap(this->m_ElevationFilled,0,0);
        this->m_GLObject_Tesselated.create();
        widget->BindGeometry(m_GLObject_Tesselated,m_Shader_Tesselated,m_Geometry_Tesselated);
    }
}

void GL3DSurface::OnRenderBefore(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{

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
        m_Shader_Channel->m_program->setUniformValue("Cposition",camera->m_Position);

        m_Shader_Channel->m_program->setUniformValue("SurfaceExtentX",(float)m_XExtent);
        m_Shader_Channel->m_program->setUniformValue("SurfaceExtentZ",(float)m_ZExtent);

        m_Shader_Channel->m_program->setUniformValue("CellSize",(float)m_CellSize);

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

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_VegCover,"VegCover",6);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_VegHeight,"VegHeight",7);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_RandomRoughness,"RandomRoughness",8);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Buildings,"Buildings",9);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_Roads,"Roads",10);

    m_Shader_Tesselated->m_program->setUniformValue("TextureSizeX",(float)2.5);
    m_Shader_Tesselated->m_program->setUniformValue("TextureSizeY",(float)2.5);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare0_Color,"grass_bare0_Color",11);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare0_Bump,"grass_bare0_Bump",12);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare0_Normal,"grass_bare0_Normal",13);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare0_Spec,"grass_bare0_Spec",14);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare1_Color,"grass_bare1_Color",15);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare1_Bump,"grass_bare1_Bump",16);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare1_Normal,"grass_bare1_Normal",17);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare1_Spec,"grass_bare1_Spec",18);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare2_Color,"grass_bare2_Color",19);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare2_Bump,"grass_bare2_Bump",20);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare2_Normal,"grass_bare2_Normal",21);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare2_Spec,"grass_bare2_Spec",22);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare3_Color,"grass_bare3_Color",23);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare3_Bump,"grass_bare3_Bump",24);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare3_Normal,"grass_bare3_Normal",25);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_bare3_Spec,"grass_bare3_Spec",26);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Color,"bare_slope1_Color",27);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Bump,"bare_slope1_Bump",28);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Normal,"bare_slope1_Normal",29);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope1_Spec,"bare_slope1_Spec",30);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope4_Color,"bare_slope4_Color",27);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope4_Bump,"bare_slope4_Bump",28);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope4_Normal,"bare_slope4_Normal",29);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_bare_slope4_Spec,"bare_slope4_Spec",30);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Color,"grass_slope1_Color",35);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Bump,"grass_slope1_Bump",36);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Normal,"grass_slope1_Normal",37);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope1_Spec,"grass_slope1_Spec",38);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope3_Color,"grass_slope3_Color",39);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope3_Bump,"grass_slope3_Bump",40);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope3_Normal,"grass_slope3_Normal",41);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope3_Spec,"grass_slope3_Spec",42);

    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope4_Color,"grass_slope4_Color",43);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope4_Bump,"grass_slope4_Bump",44);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope4_Normal,"grass_slope4_Normal",45);
    m_Shader_Tesselated->ActivateTextureOn(widget,m_Texture_grass_slope4_Spec,"grass_slope4_Spec",46);

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

        int sgnx = -1;
        int sgny = -1;

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
        if(!OUTOFMAP(m_ElevationFilled,rn + sgnx *1,cn))
        {
            if(!MV(m_ElevationFilled,rn + sgnx *1,cn))
            {
                w2 = (1.0-wx) * (wy);
                v2 =m_ElevationFilled->data[rn + sgnx *1][cn];
            }
        }
        if(!OUTOFMAP(m_ElevationFilled,rn,cn + sgny *1))
        {
            if(!MV(m_ElevationFilled,rn,cn + sgny *1))
            {
                w3 = (wx) * (1.0-wy);
                v3 =m_ElevationFilled->data[rn][cn + sgny *1];
            }
        }
        if(!OUTOFMAP(m_ElevationFilled,rn + sgnx *1,cn + sgny *1))
        {
            if(!MV(m_ElevationFilled,rn + sgnx *1,cn + sgny *1))
            {
                w4 = (wx) * (wy);
                v4 =m_ElevationFilled->data[rn + sgnx *1][cn + sgny *1];
            }
        }

        return (w1+w2+w3+w4) > 0? ((w1*v1 + w2*v2 + w3*v3 + w4 * v4)/(w1+w2+w3+w4)): 0.0;



}
