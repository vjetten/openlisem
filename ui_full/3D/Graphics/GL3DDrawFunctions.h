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


#ifndef GL3DDRAWFUNCTIONS_H
#define GL3DDRAWFUNCTIONS_H

#include "3D/World/GL3DWorld.h"
#include "3D/World/GL3DCamera.h"
#include "3D/Graphics/GL3DShaders.h"
#include "3D/Graphics/GL3DTextures.h"
#include "3D/Graphics/GL3DGeometry.h"
#include "ui_full/3D/Objects/GL3DObject.h"
#include "3D/Graphics/GL3DMapMath.h"
#include "3D/Graphics/GL3DModels.h"
#include "3D/World/GL3DSurface.h"

class GL3DSurface;

class GL3DDrawFunctions
{

public:



    static void BindGeometryInstanced(GL3DWidget * gl, QOpenGLVertexArrayObject &object,GL3DShader * s, GL3DGeometry * g, QOpenGLBuffer &MatrixBuffer);
    static void BindGeometry(GL3DWidget * gl, QOpenGLVertexArrayObject &object,GL3DShader * s, GL3DGeometry * g);
    static void DrawModelObject(GL3DWidget * gl,GL3DModel * m, GL3DCamera * camera, QMatrix4x4 ModelMatrix);
    static void DrawModelGeometryWithMaterial(GL3DWidget * gl,GL3DGeometry * g, GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, QMatrix4x4 modelmatrix);
    static void DrawModelGeometryWithMaterialMultipleStart(GL3DWidget * gl,GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera);
    static void DrawModelGeometryWithMaterialMultiple(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, QMatrix4x4 ModelMatrix);
    static void DrawModelGeometryWithMaterialMultipleEnd(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera);
    static void DrawModelInstanced(GL3DWidget * gl, GL3DModel * m, GL3DCamera * camera);
    static void DrawModelGeometryWithMaterialInstanced(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, int count);
    static void DrawModelGLInstanced(GL3DWidget * gl, GL3DModel * m, GL3DCamera * camera, GL3DSurface * surface,double dist_min, double dist_max, double dist_fade, double increment, double rand_location, double rand_rotation, double rand_scale, float * rand);
    static void DrawModelGeometryWithMaterialGLInstanced(GL3DWidget * gl, GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, int count, GL3DSurface * surface,double dist_min, double dist_max, double dist_fade, double increment, double rand_location, double rand_rotation, double rand_scale, float * rand);




};

#endif // GL3DDRAWFUNCTIONS_H
