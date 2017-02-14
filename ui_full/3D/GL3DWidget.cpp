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

#include "ui_full/3D/GL3DWidget.h"

void GL3DWidget::initializeGL()
{

    initializeOpenGLFunctions();

    gl = this;
    if ( !gl )
    {
        qFatal("Requires OpenGL >= 4.0");
        exit( 1 );
    }
    gl->initializeOpenGLFunctions();


    m_Directory = QDir::currentPath();


    m_Geometries = new GL3DGeometries();
    m_Textures = new GL3DTextures();
    m_Shaders = new GL3DShaders();
    m_Models = new GL3DModels();

    m_Geometries->Create(this);
    m_Textures->Create(this);
    m_Shaders->Create(this);
    m_Models->Create(this);

    setMouseTracking(true);
    setMinimumSize(100, 100);

    for(int i = 0; i < GL3D_INPUT_NRKEYS; i++)
    {
        KeyPressed[i] = false;
        KeyPressedT[i] = false;
    }
    for(int i = 0; i < GL3D_INPUT_NRBUTTONS; i++)
    {
        MousePressed[i] = false;
        MousePressedT[i] = false;
    }
    m_Camera = new GL3DCamera();
    m_World = new GL3DWorld();

    m_Camera->Create(this);
    m_World->Create(this);

    int error = gl->glGetError();
    if(error != GL_NO_ERROR)
    {
        qDebug() << "initialisation error 1 " <<  error;
    }

    // Set up the rendering context, define display lists etc.:
    gl->glClearColor(0.0, 0.0, 0.0, 1.0);
    gl->glEnable(GL_DEPTH_TEST);
    gl->glEnable(GL_CULL_FACE);
    gl->glCullFace(GL_FRONT);
    gl->glEnable(GL_TEXTURE_CUBE_MAP);
    gl->glEnable(GL_TEXTURE_2D);
    gl->glEnable(GL_BLEND);
    gl->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    error = gl->glGetError();
    if(error != GL_NO_ERROR)
    {
        qDebug() << "initialisation error 2 " <<  error;
    }

    this->m_Time = QDateTime::currentMSecsSinceEpoch();


    //quad for post-processing
    m_GLQuadObject.create();
    BindGeometry(m_GLQuadObject,m_Shaders->copyshader,m_Geometries->QuadGeometry);


}

void GL3DWidget::timerEvent(QTimerEvent *e)
{

    //Drawing function of parent class
    this->update();
}

void GL3DWidget::showEvent(QShowEvent *e)
{

    m_Timer.start(12, this);
}

void GL3DWidget::hideEvent(QHideEvent *e)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    m_Timer.stop();
}

void GL3DWidget::resizeGL(int w, int h)
{
    m_Camera->ResizeViewPort(w,h);

    m_Textures->CreateFrameBuffers(this,w,h);

    //for now, widget resolution is render resolution
    this->Width = w;
    this->Height = h;
}

void GL3DWidget::paintGL()
{
    //Custom function
    this->Update();

    this->Onrender();
}

void GL3DWidget::Onrender()
{

    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //gl->glBindFramebuffer(GL_FRAMEBUFFER, 0);

    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->Framebuffer);

    GLenum DrawBuffers[3] = {GL_COLOR_ATTACHMENT0,GL_COLOR_ATTACHMENT2,GL_COLOR_ATTACHMENT3};
    gl->glDrawBuffers(3, DrawBuffers);

    if(gl->glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        qDebug() << "error with framebuffer";
    }

    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //gl->glDepthMask(GL_TRUE);

    this->m_World->OnRender(this,m_Camera,m_DT);

    //gl->glDepthMask(GL_FALSE);
    GLenum DrawBuffersCopy[1] = {GL_COLOR_ATTACHMENT4};
    gl->glDrawBuffers(1, DrawBuffersCopy);

    this->m_Shaders->copyshader->m_program->bind();
    this->m_GLQuadObject.bind();
    m_Shaders->copyshader->ActivateTextureOn(this,m_Textures->RenderTexture,"tex_input",0);

    this->gl->glDrawArrays(GL_TRIANGLES,0,m_Geometries->QuadGeometry->m_IndexCount);

    GLenum DrawBuffersfinal[3] = {GL_COLOR_ATTACHMENT5};
    gl->glDrawBuffers(3, DrawBuffersfinal);

    gl->glClear(GL_DEPTH_BUFFER_BIT);

    this->m_World->OnRenderLate(this,m_Camera,m_DT);

    gl->glBindFramebuffer(GL_FRAMEBUFFER, 0);

    this->m_Shaders->blendshader->m_program->bind();

    this->m_GLQuadObject.bind();

    m_Shaders->blendshader->ActivateTextureOn(this,m_Textures->RenderTexture,"tex_input",0);
    m_Shaders->blendshader->ActivateTextureOn(this,m_Textures->RenderTextureWater,"tex_input1",1);

    this->gl->glDrawArrays(GL_TRIANGLES,0,m_Geometries->QuadGeometry->m_IndexCount);
}

void GL3DWidget::Update()
{
    qint64 time_old = m_Time;
    m_Time = QDateTime::currentMSecsSinceEpoch();
    m_DT = float(m_Time - time_old)/1000.0f;


    m_Camera->SetCurrentMatrix();

    this->m_World->OnUpdate(this,m_DT);

    UseInput();


}

void GL3DWidget::Close()
{
}

