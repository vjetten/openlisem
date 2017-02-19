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

#ifndef Textures3D
#define Textures3D

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
#include "3D/Graphics/GL3DColorRamp.h"

class GL3DWidget;
class GL3DTexture;

class GL3DTextures
{

public:
    QList<GL3DTexture *> m_DefaultTextureList;
    QList<GL3DTexture *> m_LoadedTextureList;

    GL3DTexture * LoadTextureFromMap(bool debug, cTMap * elevation, int res_x = 0,int res_y = 0, bool data = false,bool mask = false,bool fill = false,GL3DColorRamp * color_ramp = 0);
    GL3DTexture * LoadTextureFromMatrix(float * ,int lx, int ly);
    GL3DTexture * LoadTextureFromFile(QString file, bool add = true, bool is_repeat = false);
    GL3DTexture * LoadCubeTextureFromFile(QString file, bool add = true);

    void CreateFrameBuffers(GL3DWidget * widget,int w, int h);

    bool bufferscreated = false;

    // The depth buffer
    GLuint depthrenderbuffer;
    GLuint channeldepthrenderbuffer;
    GLuint RenderTexture;
    GLuint RenderTextureCopy;
    GLuint RenderTextureWater;
    GLuint LocationTexture;
    GLuint NormalTexture;
    GLuint InfoTexture;
    GLuint ChannelInfoTexture;
    GLuint ChannelTexture;
    GLuint ChannelFramebuffer;
    GLuint Framebuffer;
    GLuint FramebufferCopy;
    GLuint FramebufferWater;

    GL3DWidget * m_Widget;

    GL3DTextures()
    {

    };

    void Create(GL3DWidget * widget);
    void ClearUnused();
    void Destroy();

};

class GL3DTexture
{


public:

    QImage m_QImage;

    QImage m_QGLImage;
    QImage m_QGLImage_pos_x;
    QImage m_QGLImage_neg_x;
    QImage m_QGLImage_pos_y;
    QImage m_QGLImage_neg_y;
    QImage m_QGLImage_pos_z;
    QImage m_QGLImage_neg_z;
    QImage i_pos_x;
    QImage i_neg_x;
    QImage i_pos_y;
    QImage i_neg_y;
    QImage i_pos_z;
    QImage i_neg_z;

    int width;
    int height;

    GLuint m_GLTexture;

    bool is_created;
    bool from_matrix;
    QString Source;

    QOpenGLTexture * t;

    int nr_using_objects;

    GL3DTexture()
    {

    };

    void CreateTexture(bool debug, GL3DWidget * widget,cTMap * elevation, int res_x = 0,int res_y = 0, bool data = false, bool mask = false, bool fill = false, GL3DColorRamp * color_ramp = 0);
    void CreateTexture(GL3DWidget * widget,float * ,int lx, int ly);
    void CreateTexture(GL3DWidget * widget,QString file, bool is_repeat);
    void CreateTextureDirectPath(GL3DWidget * widget,QString file, bool is_repeat);

    void CreateCubeTexture(GL3DWidget * widget,QString file);
    QImage createSubImage(QImage* image, const QRect & rect);

    void UpdateTextureFromMap(GL3DWidget * widget, cTMap * elevation, int res_x = 0,int res_y = 0, bool data = false, bool mask = false, bool fill = false, GL3DColorRamp * color_ramp = 0);

    void ClearTexture(GL3DWidget * widget);

};

#endif
