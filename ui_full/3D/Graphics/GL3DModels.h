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

    QList<GL3DMaterial *> Material_List;
    QList<GL3DGeometry *> Geometry_List;
    QList<int> Material_Pointer;

    void LoadObjectFile(GL3DWidget * w, QString file);
    QList<GL3DMaterial *> LoadMaterialsFile(GL3DWidget * widget,QString file,QString path);

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

    GL3DTexture * Texture;
    GL3DTexture * Texture_ka;
    GL3DTexture * Texture_kd;
    GL3DTexture * Texture_ks;
    GL3DTexture * Texture_ns;
    GL3DTexture * Texture_alpha;
    GL3DTexture * Texture_disp;
    GL3DTexture * Texture_bump;
    GL3DTexture * Texture_normal;

    QVector3D color_ka;
    QVector3D color_kd;
    QVector3D color_ks;
    int spec_power;
    double alpha;

    int illum;

    GL3DShader * m_Shader;

};

#endif // GL3DMODELS_H
