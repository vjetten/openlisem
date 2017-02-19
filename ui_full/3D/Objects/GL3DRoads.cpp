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

#include "3D/Objects/GL3DRoads.h"

void GL3DRoads::OnCreate(GL3DWidget *widget)
{
    this->recieve_render = true;

    m_Texture_Color = widget->m_Textures->LoadTextureFromFile("asphalt.jpg");
    m_Texture_Normal = widget->m_Textures->LoadTextureFromFile("asphalt_normal.jpg");
    m_Texture_Bump = widget->m_Textures->LoadTextureFromFile("asphalt_bump.jpg");
    m_Texture_Specular = widget->m_Textures->LoadTextureFromFile("asphalt_spec.jpg");

    m_Model = new GL3DModel();
    m_Model->Create(widget);
    //m_Model->LoadObjectFile(widget,"house_highp_larger/Houses_Pack.obj",11);
    m_Model->BindCustomShader(widget, widget->m_Shaders->GetDefaultShader(GL3D_SHADER_ROADS));

    m_Model->AddCustomGeometry(widget,m_Geometry,0);
}

void GL3DRoads::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{

    widget->gl->glDisable(GL_DEPTH_TEST);
    widget->gl->glDisable(GL_CULL_FACE);

    QMatrix4x4 matrix;
    matrix.setToIdentity();

    //widget->DrawModelObject(m_Model,camera,matrix);

    widget->gl->glDisable(GL_CULL_FACE);

    m_Model->m_Shader->m_program->bind();
    m_Model->m_Shader->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    m_Model->m_Shader->m_program->setUniformValue("Mmatrix",matrix);
    m_Model->m_Shader->m_program->setUniformValue("Cposition",camera->m_Position);

    m_Model->m_Shader->m_program->setUniformValue("SurfaceExtentX",(float)m_Surface->m_XExtent);
    m_Model->m_Shader->m_program->setUniformValue("SurfaceExtentZ",(float)m_Surface->m_ZExtent);

    m_Model->m_Shader->m_program->setUniformValue("CellSize",(float)m_Surface->m_CellSize);

    m_Model->m_Shader->ActivateTextureOn(widget,m_Texture_Color,"Texture_Color",4);
    m_Model->m_Shader->ActivateTextureOn(widget,m_Texture_Normal,"Texture_Normal",5);
    m_Model->m_Shader->ActivateTextureOn(widget,m_Texture_Bump,"Texture_Bump",6);
    m_Model->m_Shader->ActivateTextureOn(widget,m_Texture_Specular,"Texture_Specular",7);

    m_Model->m_Shader->ActivateTextureOn(widget,m_Surface->m_Texture,"heightMap",0);
    m_Model->m_Shader->ActivateTextureOn(widget,m_Surface->m_Texture_Mask,"mask",1);

    m_Model->m_Shader->ActivateTextureOn(widget,m_Surface->m_Texture_SlopeX,"slopeX",2);
    m_Model->m_Shader->ActivateTextureOn(widget,m_Surface->m_Texture_SlopeY,"slopeY",3);

    for(int i = 0; i < m_Model->GLVAO_List.length(); i++)
    {
        m_Model->GLVAO_List.at(i)->bind();

        widget->gl->glDrawElements(GL_TRIANGLES, m_Model->Geometry_List.at(i)->m_IndexCount, GL_UNSIGNED_INT, 0);
        m_Model->GLVAO_List.at(i)->release();
    }

    int e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in drawing roads " <<e;
    }


    m_Model->m_Shader->m_program->release();

    widget->gl->glEnable(GL_DEPTH_TEST);
    widget->gl->glEnable(GL_CULL_FACE);

}

void GL3DRoads::OnDestroy(GL3DWidget *widget)
{



}

void GL3DRoads::SetRoadDistribution(GL3DWidget * w,GL3DSurface * s,cTMap * rw)
{
    qDebug() << "initialize roads" << rw;

    m_Surface = s;

    int n = 0;
    FOR_ROW_COL_MV(rw,rw)
    {

        if(rw->Drc > 0.0)
        {
            n++;
        }
    }

    int nrrows = rw->nrRows();
    int nrcols = rw->nrCols();

    double cs = rw->cellSize();

    Vertex * vertices = (Vertex *)malloc(sizeof(Vertex) * (n * 4));
    GLuint * indices = (GLuint *)malloc(sizeof(GLuint)*(  n * 6));

    int i = 0;

    qDebug() << "create for cells" << n;
    FOR_ROW_COL_MV(rw,rw)
    {

        if(rw->Drc > 0.0)
        {
            bool top = false;
            bool bottom = false;
            bool left = false;
            bool right = false;

            bool nw = false;
            bool ne = false;
            bool sw = false;
            bool se = false;
            if(r + 1 < nrrows)
            {
                if(!pcr::isMV(rw->data[r+1][c]))
                {
                    if(rw->data[r+1][c] > 0.0)
                    {
                        bottom = true;
                    }
                }
            }
            if(r - 1 > 0)
            {
                if(!pcr::isMV(rw->data[r-1][c]))
                {
                    if(rw->data[r-1][c] > 0.0)
                    {
                        top = true;
                    }
                }
            }
            if(c + 1 < nrcols)
            {
                if(!pcr::isMV(rw->data[r][c+1]))
                {
                    if(rw->data[r][c+1] > 0.0)
                    {
                        right = true;
                    }
                }
            }
            if(c - 1 > 0)
            {
                if(!pcr::isMV(rw->data[r][c-1]))
                {
                    if(rw->data[r][c-1] > 0.0)
                    {
                        left = true;
                    }
                }
            }



            if(r + 1 < nrrows && c - 1 > 0)
            {
                if(!pcr::isMV(rw->data[r+1][c-1]))
                {
                    if(rw->data[r+1][c-1] > 0.0)
                    {
                        sw = true;
                    }
                }
            }
            if(r + 1< nrrows && c + 1 < nrcols)
            {
                if(!pcr::isMV(rw->data[r+1][c-1]))
                {
                    if(rw->data[r+1][c-1] > 0.0)
                    {
                        se = true;
                    }
                }
            }
            if(r - 1 > 0 && c - 1 > 0)
            {
                if(!pcr::isMV(rw->data[r-1][c-1]))
                {
                    if(rw->data[r-1][c-1] > 0.0)
                    {
                        nw = true;
                    }
                }
            }
            if(r - 1 > 0 && c + 1 < nrcols)
            {
                if(!pcr::isMV(rw->data[r-1][c+1]))
                {
                    if(rw->data[r-1][c+1] > 0.0)
                    {
                        ne = true;
                    }
                }
            }

            bool both_y = top && bottom;
            bool both_x = right && left;

            double x1,x2,x3,x4,y1,y2,y3,y4;

            double cx1 = c * cs;
            double cx2 = c * cs + cs;
            double cy1 = r * cs;
            double cy2 = r * cs + cs;

            double crw = rw->Drc;
            double crwf = crw/cs;
            double hcrwf = 0.5 * crwf;
            double dist = (cs - crw)/2.0;//hcrwf * cs;
            double distd = (rw->Drc*0.5)/sqrt(2.0);

            if(both_y && !both_x)
            {
                //road is aligned in y-direction
                x1 = cx1 + dist;
                x2 = cx2 - dist;
                x3 = cx1 + dist;
                x4 = cx2 - dist;
                y1 = cy1;
                y2 = cy1;
                y3 = cy2;
                y4 = cy2;

            }else if(both_x && !both_y)
            {
                //road is aligned in x-direction
                x1 = cx1;
                x2 = cx1;
                x3 = cx2;
                x4 = cx2;
                y1 = cy1 + dist;
                y2 = cy2 - dist;
                y3 = cy1 + dist;
                y4 = cy2 - dist;

            }else if(left && bottom)
            {
                x1 = cx1;
                x2 = cx1;
                x3 = cx1 + dist;
                x4 = cx2 - dist;
                y1 = cy1 + dist;
                y2 = cy2 - dist;
                y3 = cy2;
                y4 = cy2;

            }else if(right && bottom)
            {
                x1 = cx2;
                x2 = cx2;
                x3 = cx1 + dist;
                x4 = cx2 - dist;
                y1 = cy1 + dist;
                y2 = cy2 - dist;
                y3 = cy2;
                y4 = cy2;

            }else if(top && left)
            {
                x1 = cx1 + dist;
                x2 = cx2 - dist;
                x3 = cx1;
                x4 = cx1;
                y1 = cy1;
                y2 = cy1;
                y3 = cy1 + dist;
                y4 = cy2 - dist;

            }else if(top && right)
            {
                x1 = cx1 + dist;
                x2 = cx2 - dist;
                x3 = cx2;
                x4 = cx2;
                y1 = cy1;
                y2 = cy1;
                y3 = cy1 + dist;
                y4 = cy2 - dist;

            }else if(left && ne)
            {

                x1 = cx1;
                x2 = cx1;
                x3 = cx2 + distd;
                x4 = cx2 - distd;
                y1 = cy1 + dist;
                y2 = cy2 - dist;
                y3 = cy1 + distd;
                y4 = cy1 - distd;
            }else if(left && se)
            {
                x1 = cx1;
                x2 = cx1;
                x3 = cx2 + distd;
                x4 = cx2 - distd;
                y1 = cy1 + dist;
                y2 = cy2 - dist;
                y3 = cy2 + distd;
                y4 = cy2 - distd;
            }else if(right && nw)
            {

                x1 = cx1 - distd;
                x2 = cx1 + distd;
                x3 = cx2;
                x4 = cx2;
                y1 = cy1 - distd;
                y2 = cy1 + distd;
                y3 = cy1 + dist;
                y4 = cy2 - dist;

            }else if(right && sw)
            {
                x1 = cx1 - distd;
                x2 = cx1 + distd;
                x3 = cx2;
                x4 = cx2;
                y1 = cy2 - distd;
                y2 = cy2 + distd;
                y3 = cy1 + dist;
                y4 = cy2 - dist;

            }else if(bottom && nw)
            {
                x1 = cx1 - distd;
                x2 = cx1 + distd;
                x3 = cx1 - dist;
                x4 = cx2 + dist;
                y1 = cy1 - distd;
                y2 = cy1 + distd;
                y3 = cy2;
                y4 = cy2;

            }else if(bottom && ne)
            {
                x1 = cx2 - distd;
                x2 = cx2 + distd;
                x3 = cx1 + dist;
                x4 = cx2 - dist;
                y1 = cy1 - distd;
                y2 = cy1 + distd;
                y3 = cy2;
                y4 = cy2;
            }else if(top && sw)
            {
                x1 = cx1 + dist;
                x2 = cx2 - dist;
                x3 = cx1 + distd;
                x4 = cx1 - distd;
                y1 = cy1;
                y2 = cy1;
                y3 = cy2 + distd;
                y4 = cy2 - distd;
            }else if(top && se)
            {
                x1 = cx1 + dist;
                x2 = cx2 - dist;
                x3 = cx2 + distd;
                x4 = cx2 - distd;
                y1 = cy1;
                y2 = cy1;
                y3 = cy2 + distd;
                y4 = cy2 - distd;
            }else if(nw && se)
            {
                x1 = cx1 + distd;
                x2 = cx1 - distd;
                x3 = cx2 + distd;
                x4 = cx2 - distd;
                y1 = cy2 + distd;
                y2 = cy2 - distd;
                y3 = cy1 + distd;
                y4 = cy1 - distd;
            }else if(ne && sw)
            {
                x1 = cx2 + distd;
                x2 = cx2 - distd;
                x3 = cx1 + distd;
                x4 = cx1 - distd;
                y1 = cy2 + distd;
                y2 = cy2 - distd;
                y3 = cy1 + distd;
                y4 = cy1 - distd;
            }else if(bottom)
            {
                x1 = cx1 + dist;
                x2 = cx2 - dist;
                x3 = cx1 + dist;
                x4 = cx2 - dist;
                y1 = cy1 + 0.5 * cs;
                y2 = cy1 + 0.5 * cs;
                y3 = cy2;
                y4 = cy2;

            }else if(top)
            {
                x1 = cx1 + dist;
                x2 = cx2 - dist;
                x3 = cx1 + dist;
                x4 = cx2 - dist;
                y1 = cy1;
                y2 = cy1;
                y3 = cy2 - 0.5 * cs;
                y4 = cy2 - 0.5 * cs;

            }else if(right)
            {
                x1 = cx1+0.5 *cs;
                x2 = cx1+0.5 *cs;
                x3 = cx2;
                x4 = cx2;
                y1 = cy1 + dist;
                y2 = cy2 - dist;
                y3 = cy1 + dist;
                y4 = cy2 - dist;
            }else if(left)
            {
                x1 = cx1;
                x2 = cx1;
                x3 = cx2 - 0.5 * cs;
                x4 = cx2 - 0.5 * cs;
                y1 = cy1 + dist;
                y2 = cy2 - dist;
                y3 = cy1 + dist;
                y4 = cy2 - dist;
            }else

            double z1 = s->GetElevation(x1,y1);
            double z2 = s->GetElevation(x2,y2);
            double z3 = s->GetElevation(x3,y3);
            double z4 = s->GetElevation(x4,y4);

            Vertex v1 = Vertex(QVector3D(x1,0.0,y1));
            Vertex v2 = Vertex(QVector3D(x2,0.0,y2));
            Vertex v3 = Vertex(QVector3D(x3,0.0,y3));
            Vertex v4 = Vertex(QVector3D(x4,0.0,y4));

            if(i < n)
            {
                vertices[i * 4 + 0] = v1;
                vertices[i * 4 + 1] = v2;
                vertices[i * 4 + 2] = v3;
                vertices[i * 4 + 3] = v4;

                indices[i*6 + 0] = (GLuint)(i * 4 + 0);
                indices[i*6 + 1] = (GLuint)(i * 4 + 2);
                indices[i*6 + 2] = (GLuint)(i * 4 + 1);
                indices[i*6 + 3] = (GLuint)(i * 4 + 2);
                indices[i*6 + 4] = (GLuint)(i * 4 + 1);
                indices[i*6 + 5] = (GLuint)(i * 4 + 3);
            }

            i++;

        }

    }


    m_Geometry = new GL3DGeometry();
    m_Geometry->CreateGeometry(w,vertices,n * 4,indices,n * 6);
    m_GLObject.create();


    m_Shader = w->m_Shaders->GetDefaultShader(GL3D_SHADER_ROADS);
    w->BindGeometry(m_GLObject,m_Shader,m_Geometry);

    qDebug() << " road placed in cells Nr: " << n;
}
