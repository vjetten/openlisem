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

#ifndef GL3DMODELS_H
#define GL3DMODELS_H

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

class GL3DWidget;
class GL3DModel;
class GL3DMaterial;

class GL3DModels
{
public:
    GL3DModels()
    {

    };

    GL3DWidget * m_Widget;

    void Create(GL3DWidget * widget);
    void ClearUnused();
    void Destroy();


};

class GL3DModel
{

public:

    bool Is_Loaded = false;

    GL3DModel()
    {

    };

    GL3DShader * m_Shader;

    QList<GL3DMaterial *> Material_List;
    QList<GL3DGeometry *> Geometry_List;
    QList<int> Material_Pointer;
    QList<QOpenGLVertexArrayObject *> GLVAO_List;

    void LoadObjectFile(GL3DWidget * w, QString file, float rescale = 1.0);
    QList<GL3DMaterial *> LoadMaterialsFile(GL3DWidget * widget,QString file,QString path);

    void AddCustomGeometry(GL3DWidget * w, GL3DGeometry * g, GL3DMaterial * m);
    void BindCustomShader(GL3DWidget * w, GL3DShader * s);

    void Create(GL3DWidget * widget);
    void Destroy(GL3DWidget * widget);

    struct vec3d
    {
        float x;
        float y;
        float z;
    };

};

class GL3DMaterial
{
public:
    QString name;

    GL3DTexture * Texture = 0;
    GL3DTexture * Texture_ka = 0;
    GL3DTexture * Texture_kd = 0;
    GL3DTexture * Texture_ks = 0;
    GL3DTexture * Texture_ns = 0;
    GL3DTexture * Texture_alpha = 0;
    GL3DTexture * Texture_disp = 0;
    GL3DTexture * Texture_bump = 0;
    GL3DTexture * Texture_normal = 0;

    QVector3D color_ka = QVector3D(1.0,1.0,1.0);
    QVector3D color_kd = QVector3D(1.0,1.0,1.0);
    QVector3D color_ks = QVector3D(1.0,1.0,1.0);
    float spec_power = 2.0;
    double alpha = 1.0;

    int illum = 2;


};

#endif // GL3DMODELS_H
