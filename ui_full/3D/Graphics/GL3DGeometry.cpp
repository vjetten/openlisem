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

#include "ui_full/3D/Graphics/GL3DGeometry.h"

void GL3DGeometries::Create(GL3DWidget * widget)
{

    m_Widget = widget;
    QuadGeometry = new GL3DGeometry();
    QuadGeometry->CreateGeometry(widget,g_quad_vertex_buffer_data,6);
}

GL3DGeometry * GL3DGeometries::LoadGeometryFromMap(cTMap * elevation, int m, bool data2d)
{

    GL3DGeometry *g = new GL3DGeometry();
    g->CreateGeometry(m_Widget,elevation,m,data2d);
    this->m_GeometryList.append(g);
    return g;
}

GL3DGeometry * GL3DGeometries::LoadGeometryFromArray(const Vertex * data,int lv,const GLuint * indices, int li)
{
    GL3DGeometry *g = new GL3DGeometry();
    g->CreateGeometry(m_Widget,data,lv,indices,li);
    this->m_GeometryList.append(g);
    return g;
}

GL3DGeometry * GL3DGeometries::LoadGeometryFromFile(QString file)
{



    return 0;
}


void GL3DGeometries::ClearUnused()
{


}

void GL3DGeometries::Destroy()
{



}

double GL3DGeometry::GetMapValue(cTMap * map,double y, double x)
{

    double cs = map->cellSize();
    int cn = std::floor(x/cs);
    int rn = std::floor(y/cs);

    double v1 = 0;
    double v2 = 0;
    double v3 = 0;
    double v4 = 0;

    double w1 = 0;
    double w2 = 0;
    double w3 = 0;
    double w4 = 0;

    if(!OUTOFMAP(map,rn,cn))
    {
        if(!MV(map,rn,cn))
        {
            w1 = 1;
            v1 =map->data[rn][cn];
        }
    }
    if(!OUTOFMAP(map,rn+1,cn))
    {
        if(!MV(map,rn+1,cn))
        {
            w2 = 1;
            v2 =map->data[rn+1][cn];
        }
    }
    if(!OUTOFMAP(map,rn,cn+1))
    {
        if(!MV(map,rn,cn+1))
        {
            w3 = 1;
            v3 =map->data[rn][cn+1];
        }
    }
    if(!OUTOFMAP(map,rn+1,cn+1))
    {
        if(!MV(map,rn+1,cn+1))
        {
            w4 = 1;
            v4 =map->data[rn+1][cn+1];
        }
    }


    return (w1+w2+w3+w4) > 0? ((w1*v1 + w2*v2 + w3*v3 + w4 * v4)/(w1+w2+w3+w4)): 0.0;
}

void GL3DGeometry::CreateGeometry(QGLWidget * widget,cTMap * map,int m, bool data2d)
{
    uses_index = true;
    this->is_2d_data = data2d;

    int n = 0;
    for(int r = 0; r < map->nrRows(); r = r + m){
    for (int c = 0; c < map->nrCols(); c= c + m){
                n++;
    }}

    int oresx = map->nrCols();
    int oresy = map->nrRows();

    int xmin = oresx;
    int xmax = 0;
    int ymin = oresy;
    int ymax = 0;
    bool first = true;
    double emax = 0;
    double emin = 0;

    FOR_ROW_COL_MV(map,map)
    {
        if(first)
        {
            first = false;
            emax = map->Drc;
            emin = map->Drc;
        }

        xmin = std::min(c,xmin);
        xmax = std::max(c,xmax);
        ymin = std::min(r,ymin);
        ymax = std::max(r,ymax);
        emax = std::max(emax,map->Drc);
        emin = std::min(emin,map->Drc);
    }

    int needed_resx = xmax-xmin;
    int needed_resy = ymax-ymin;

    double erange = emax-emin;

    int nr_triangles;
    int nr_vertex;
    int nr_index;

    if(!data2d)
    {
        nr_triangles = n*2;
        nr_vertex = n*2 * 3;
        nr_index = n*2 * 3;
    }else
    {
        nr_triangles = n;
        nr_vertex = n*1;
        nr_index = n*1;
    }

    m_IndexCount = nr_index;
    m_PatchCount = nr_index;

    qDebug() << "allocate vertex memory for surface";
    Vertex * vertexdata;
    Vertex2D * vertexdata2d;
    if(!data2d)
    {
        vertexdata =(Vertex*) malloc(sizeof(Vertex) * nr_triangles * 3);
    }else
    {
        vertexdata2d =(Vertex2D*) malloc(sizeof(Vertex) * nr_triangles * 1);
    }

    GLuint * indexdata;

    if(!data2d)
    {
        indexdata = (GLuint*)malloc(sizeof(GLuint) * nr_triangles * 3);
    }

    qDebug() << "succesfully allocated memory";

    int nt =0;
    double cs =  map->cellSize();
    double csm = m * map->cellSize();

    qDebug() << "create vertices";


    for(int r = 0; r < map->nrRows(); r = r + m){
    for (int c = 0; c < map->nrCols(); c= c + m){

        double c2 = (double(c))-0.5;
        double r2 = (double(r))-0.5;

        double e1 = GetMapValue(map,r2*cs,c2*cs);
        double e2 = GetMapValue(map,r2*cs,(c2+1)*cs);
        double e3 = GetMapValue(map,(r2+1)*cs,c2*cs);

        if(!data2d)
        {
            vertexdata[nt*6 + 0] = Vertex(QVector3D((double)c2*cs,e1,(double)r2*cs),QVector3D((e1-emin)/(erange),(e1-emin)/(erange),(e1-emin)/(erange)));
            vertexdata[nt*6 + 1] = Vertex(QVector3D((double)(c2+m)*cs,e2,(double)r2*cs),QVector3D((e2-emin)/(erange),(e2-emin)/(erange),(e2-emin)/(erange)));
            vertexdata[nt*6 + 2] = Vertex(QVector3D((double)c2*cs,e3,(double)(r2+m)*cs),QVector3D((e3-emin)/(erange),(e3-emin)/(erange),(e3-emin)/(erange)));

            e1 = GetMapValue(map,(r2+1)*cs,(c2+1)*cs);
            vertexdata[nt*6 + 3] = Vertex(QVector3D((double)(c2+m)*cs,e1,(r2+m)*cs),QVector3D((e1-emin)/(erange),(e1-emin)/(erange),(e1-emin)/(erange)));
            vertexdata[nt*6 + 4] = Vertex(QVector3D((double)c2*cs,e3,(double)(r2+m)*cs),QVector3D((e3-emin)/(erange),(e3-emin)/(erange),(e3-emin)/(erange)));
            vertexdata[nt*6 + 5] = Vertex(QVector3D((double)(c2+m)*cs,e2,(double)r2*cs),QVector3D((e2-emin)/(erange),(e2-emin)/(erange),(e2-emin)/(erange)));

            indexdata[nt*6 + 0] = (GLuint)(nt*6 + 0);
            indexdata[nt*6 + 1] = (GLuint)(nt*6 + 1);
            indexdata[nt*6 + 2] = (GLuint)(nt*6 + 2);
            indexdata[nt*6 + 3] = (GLuint)(nt*6 + 3);
            indexdata[nt*6 + 4] = (GLuint)(nt*6 + 4);
            indexdata[nt*6 + 5] = (GLuint)(nt*6 + 5);

        }else
        {
            vertexdata2d[nt*4 + 0] = Vertex2D(QVector2D((double)c2*cs,(double)r2*cs));
            //vertexdata2d[nt*4 + 1] = Vertex2D(QVector2D((double)(c2+1.0)*cs,(double)r2*cs));
            //vertexdata2d[nt*4 + 2] = Vertex2D(QVector2D((double)c2*cs,(double)(r2+1.0)*cs));

            //vertexdata2d[nt*4 + 3] = Vertex2D(QVector2D((double)(c2+1.0)*cs,(r2+1.0)*cs));
            //vertexdata2d[nt*4 + 4] = Vertex2D(QVector2D((double)c2*cs,(double)(r2+1.0)*cs));
            //vertexdata2d[nt*6 + 5] = Vertex2D(QVector2D((double)(c2+1.0)*cs,(double)r2*cs));

        }

        nt ++;
    }}
    qDebug() << "succesfully created vertices";

    qDebug() << "create and fill buffers";
    QOpenGLBuffer temp1(QOpenGLBuffer::VertexBuffer);
    m_vertex = temp1;
    m_vertex.create();
    m_vertex.bind();
    m_vertex.setUsagePattern(QOpenGLBuffer::StaticDraw);
    if(!data2d)
    {
        m_vertex.allocate(vertexdata, sizeof(Vertex) * nr_vertex);
    }else
    {
        m_vertex.allocate(vertexdata2d, sizeof(Vertex2D) * nr_vertex);
    }
    m_vertex.release();

    if(!data2d)
    {
        QOpenGLBuffer temp2(QOpenGLBuffer::IndexBuffer);
        m_index = temp2;
        m_index.create();
        m_index.bind();
        m_index.setUsagePattern(QOpenGLBuffer::StaticDraw);
        m_index.allocate(indexdata, sizeof(GLuint) * nr_index);
        m_index.release();
    }

    qDebug() << "succesfully created buffers";

    if(!data2d)
    {
        free(vertexdata);
        free(indexdata);
    }else
    {
        free(vertexdata2d);
    }

}

void GL3DGeometry::CreateGeometry(QGLWidget * widget,const Vertex * data,int lv,const GLuint * indices, int li)
{
    uses_index = true;
    is_2d_data = false;

    QOpenGLBuffer temp1(QOpenGLBuffer::VertexBuffer);
    m_vertex = temp1;
    m_vertex.create();
    m_vertex.bind();
    m_vertex.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_vertex.allocate(data, sizeof(Vertex) * lv);
    m_vertex.release();

    QOpenGLBuffer temp2(QOpenGLBuffer::IndexBuffer);
    m_index = temp2;
    m_index.create();
    m_index.bind();
    m_index.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_index.allocate(indices, sizeof(GLuint) * li);
    m_index.release();

    m_IndexCount = li;
}

void GL3DGeometry::CreateGeometry(QGLWidget * widget,const Vertex2D * data,int lv,const GLuint * indices, int li)
{
    uses_index = true;
    is_2d_data = true;

    QOpenGLBuffer temp1(QOpenGLBuffer::VertexBuffer);
    m_vertex = temp1;
    m_vertex.create();
    m_vertex.bind();
    m_vertex.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_vertex.allocate(data, sizeof(Vertex2D) * lv);
    m_vertex.release();

    QOpenGLBuffer temp2(QOpenGLBuffer::IndexBuffer);
    m_index = temp2;
    m_index.create();
    m_index.bind();
    m_index.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_index.allocate(indices, sizeof(GLuint) * li);
    m_index.release();

    m_IndexCount = li;
}


void GL3DGeometry::CreateGeometry(QGLWidget * widget,const Vertex * data,int lv)
{
    uses_index = false;
    is_2d_data = false;

    QOpenGLBuffer temp1(QOpenGLBuffer::VertexBuffer);
    m_vertex = temp1;
    m_vertex.create();
    m_vertex.bind();
    m_vertex.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_vertex.allocate(data, sizeof(Vertex) * lv);
    m_vertex.release();

    m_IndexCount = lv;
}

void GL3DGeometry::CreateGeometry(QGLWidget * widget,QString file)
{



}

void GL3DGeometry::ClearGeometry()
{



}
