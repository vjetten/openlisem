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

#ifndef Geometries3D
#define Geometries3D

#include <QtOpenGL/QGLWidget>
#include "QDebug"
#include <QVector>
#include <vector>
#include <QtOpenGL/QtOpenGL>
#include "csfMap.h"
#include "LisUIoutput.h"
#include "global.h"
#include <QString>
#include "ui_full/3D/GL3DWidget.h"
#include "ui_full/3D/Graphics/GL3DMath.h"

class GL3DWidget;
class GL3DGeometry;

//Colored triangle
static const Vertex sg_vertexes[] = {
  Vertex( QVector3D( 0.0f,  1.0f, 0.0f), QVector3D(1.0f, 0.0f, 0.0f) ),
  Vertex( QVector3D( 0.0f, 0.0f, 1.0f), QVector3D(0.0f, 1.0f, 0.0f) ),
  Vertex( QVector3D( 0.0f, 1.0f, 1.0f), QVector3D(0.0f, 0.0f, 1.0f) )
};
static const GLuint sg_indices[] = {0,1,2};

#define SQRT_3_3 10
static const Vertex vertices[8] = {
    Vertex( QVector3D( -SQRT_3_3,	-SQRT_3_3,	-SQRT_3_3),QVector3D(1.0f, 0.0f, 0.0f)),
    Vertex( QVector3D( SQRT_3_3,	-SQRT_3_3,	-SQRT_3_3),QVector3D(0.0f, 1.0f, 0.0f)),
    Vertex( QVector3D( -SQRT_3_3,	SQRT_3_3,	-SQRT_3_3),QVector3D(1.0f, 0.0f, 1.0f)),
    Vertex( QVector3D( SQRT_3_3,	SQRT_3_3,	-SQRT_3_3),QVector3D(1.0f, 1.0f, 0.0f)),
    Vertex( QVector3D( -SQRT_3_3,	-SQRT_3_3,	SQRT_3_3),QVector3D(0.0f, 0.0f, 1.0f)),
    Vertex( QVector3D( SQRT_3_3,	-SQRT_3_3,	SQRT_3_3),QVector3D(1.0f, 1.0f, 1.0f)),
    Vertex( QVector3D( -SQRT_3_3,	SQRT_3_3,	SQRT_3_3),QVector3D(0.0f, 1.0f, 1.0f)),
    Vertex( QVector3D( SQRT_3_3,	SQRT_3_3,	SQRT_3_3),QVector3D(0.0f, 0.0f, 0.0f))
};
static const GLuint indices[24] = {
    1,			5,			7,			3,	// positive x
    2,			0,			4,			6,	// negative x
    4,			5,			7,			6,	// positive y
    0,			1,			3,			2,	// negative y
    0,			1,			5,			4,	// positive z
    3,			2,			6,			7	// negative z
};

static const Vertex vertices_SKYBOX2[36] =
{
    Vertex( QVector3D(-10.0f,  10.0f, -10.0f)),
    Vertex( QVector3D(-10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f, -10.0f)),
    Vertex( QVector3D(-10.0f,  10.0f, -10.0f)),

    Vertex( QVector3D(-10.0f, -10.0f,  10.0f)),
    Vertex( QVector3D(-10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D(-10.0f,  10.0f, -10.0f)),
    Vertex( QVector3D(-10.0f,  10.0f, -10.0f)),
    Vertex( QVector3D(-10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D(-10.0f, -10.0f,  10.0f)),

    Vertex( QVector3D( 10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f, -10.0f)),

    Vertex( QVector3D(-10.0f, -10.0f,  10.0f)),
    Vertex( QVector3D(-10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f,  10.0f)),
    Vertex( QVector3D(-10.0f, -10.0f,  10.0f)),

    Vertex( QVector3D(-10.0f,  10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D(-10.0f,  10.0f,  10.0f)),
    Vertex( QVector3D(-10.0f,  10.0f, -10.0f)),

    Vertex( QVector3D(-10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D(-10.0f, -10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f, -10.0f)),
    Vertex( QVector3D(-10.0f, -10.0f,  10.0f)),
    Vertex( QVector3D( 10.0f, -10.0f,  10.0f))
};
static const Vertex vertices_SKYBOX[36] =
{
    Vertex( QVector3D(-1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f, -1.0f)),

    Vertex( QVector3D(-1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D(-1.0f, -1.0f,  1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f, -1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f, -1.0f)),
    Vertex( QVector3D(-1.0f, -1.0f,  1.0f)),

    Vertex( QVector3D( 1.0f, -1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f, -1.0f)),

    Vertex( QVector3D(-1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D(-1.0f, -1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D(-1.0f, -1.0f,  1.0f)),

    Vertex( QVector3D( 1.0f,  1.0f, -1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f,  1.0f,  1.0f)),
    Vertex( QVector3D(-1.0f,  1.0f, -1.0f)),

    Vertex( QVector3D(-1.0f, -1.0f,  1.0f)),
    Vertex( QVector3D(-1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D(-1.0f, -1.0f,  1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f, -1.0f)),
    Vertex( QVector3D( 1.0f, -1.0f,  1.0f))
};

static const GLuint indices_SKYBOX[36] = {
    0,			1,			2,			3,
    4,			5,			6,			7,
    8,			9,			10,			11,
    12,			13,			14,			15,
    16,			17,			18,			19,
    20,			21,			22,			23,
    24,			25,			26,			27,
    28,			29,			30,			31,
    32,			33,			34,			35
};

static const Vertex g_quad_vertex_buffer_data[] = {
    Vertex( QVector3D(-1.0f, -1.0f, 0.0f)),
    Vertex( QVector3D(-1.0f, 1.0f, 0.0f)),
    Vertex( QVector3D(1.0f, -1.0f, 0.0f)),
    Vertex( QVector3D(1.0f,  -1.0f, 0.0f)),
    Vertex( QVector3D(-1.0f, 1.0f, 0.0f)),
    Vertex( QVector3D(1.0f,  1.0f, 0.0f)),
};

class GL3DGeometries
{

public:
    QList<GL3DGeometry *> m_GeometryList;

    GL3DGeometry * LoadGeometryRaster(int r, int c, double size);
    GL3DGeometry * LoadGeometryFromMap(cTMap * elevation, int m = 1, bool data2d = false);
    GL3DGeometry * LoadGeometryFromArray(const Vertex * data, int lv,const GLuint * indices, int li);
    GL3DGeometry * LoadGeometryFromFile(QString file);

    GL3DWidget * m_Widget;

    GL3DGeometry * QuadGeometry;

    GL3DGeometries()
    {

    };

    void Create(GL3DWidget * widget);
    void ClearUnused();
    void Destroy();

};

class GL3DGeometry
{


public:

    bool is_2d_data = false;
    bool is_created = false;
    bool from_matrix = false;
    QString Source = "";

    int m_IndexCount = 0;
    int m_PatchCount = 0;

    int nr_using_objects;

    bool uses_index;

    QOpenGLBuffer m_vertex;
    QOpenGLBuffer m_index;

    GL3DGeometry()
    {

    };
    double GetMapValue(cTMap * map,double x, double y);

    void CreateGeometryRaster(QGLWidget * widget,int r, int c, double size);
    void CreateGeometry(QGLWidget * widget,cTMap * map, int m = 1, bool data2d = false);
    void CreateGeometry(QGLWidget * widget,const Vertex * data, int lv,const GLuint * indices, int li);
    void CreateGeometry(QGLWidget * widget,const Vertex2D * data, int lv,const GLuint * indices, int li);
    void CreateGeometry(QGLWidget * widget,const Vertex * data, int lv);
    void CreateGeometry(QGLWidget * widget,QString file);
    void ClearGeometry();

};

#endif
