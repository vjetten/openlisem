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

#ifndef Shaders3D
#define Shaders3D

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

#define SHADER_NULL 0
class GL3DWidget;
class GL3DShader;


class GL3DShaders
{

public:

    QList<GL3DShader *> m_DefaultShaderList;
    QList<GL3DShader *> m_LoadedShaderList;

    GL3DShader * LoadShaderFromText(const char * text_ps,const char * text_vs, bool add = true,const char * text_tcs = SHADER_NULL,const char * text_tes = SHADER_NULL,const char * text_gs = SHADER_NULL);
    GL3DShader * LoadShaderFromFile(QString file_ps,QString file_vs, bool add = true,QString file_tcs="",QString file_tes = "",QString file_gs = "");

    GL3DShader * GetDefaultShader(int index);

    GL3DShader * rainshader;
    GL3DShader * channelshader;
    GL3DShader * copyshader;
    GL3DShader * blendshader;

    GL3DShaders()
    {

    };

    GL3DWidget * m_Widget;

    void Create(GL3DWidget * widget);
    void ClearUnused();
    void Destroy();

};

class GL3DShader
{


public:


    QOpenGLShaderProgram *m_program;

    bool is_created;
    bool from_text;
    QString source;

    int nr_using_ojects;

    GL3DShader()
    {

    };

    void LoadShaderFromText(GL3DWidget * widget,const char * text_ps,const char * text_vs,const char * text_tcs = SHADER_NULL,const char * text_tes = SHADER_NULL,const char * text_gs = SHADER_NULL);
    void LoadShaderFromFile(GL3DWidget * widget,QString file_ps,QString file_vs,QString file_tcs = "",QString file_tes = "",QString file_gs = "");

    void ActivateTextureOn(GL3DWidget * widget, GL3DTexture *t, const char * name, int slot);
    void ActivateTextureOn(GL3DWidget * widget, GLuint t, const char * name, int slot);

    void ActivateCubeTextureOn(GL3DWidget * widget, GLuint t, const char * name, int slot);

    void ClearShader();

};

#define GL3D_SHADER_SIMPLE 0
#define GL3D_SHADER_SKYBOX 1
#define GL3D_SHADER_SURFACE_TESSELATED 2
#define GL3D_SHADER_SURFACE_FLOW 3
#define GL3D_SHADER_MODEL 4
#define GL3D_SHADER_ROADS 5
#define GL3D_SHADER_CHANNEL 6

/////////////////////////////////////////////////////////////////////////////////////
///Simple shader
/////////////////////////////////////////////////////////////////////////////////////
static const char *vertexShaderSource =
    "in highp vec3 posAttr;\n"
    "in lowp vec3 colAttr;\n"
    "uniform highp mat4 Cmatrix;\n"
    "uniform highp mat4 Mmatrix;\n"
    "out lowp vec4 col;\n"
    "void main() {\n"
    "   gl_Position = Cmatrix *vec4(posAttr, 1.0);\n"
    "   col = vec4(colAttr,1.0);\n"
    "}\n";
//;
static const char *fragmentShaderSource =
    "in lowp vec4 col;\n"
    "void main() {\n"
    "   gl_FragColor = col;\n"
    "}\n";

#endif
