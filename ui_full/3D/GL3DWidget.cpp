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
        return;
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
        return;
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
        return;
    }

    this->m_Time = QDateTime::currentMSecsSinceEpoch();
    this->m_Time_Start = QDateTime::currentMSecsSinceEpoch();

    //quad for post-processing
    m_GLQuadObject.create();
    GL3DDrawFunctions::BindGeometry(this,m_GLQuadObject,m_Shaders->copyshader,m_Geometries->QuadGeometry);

    m_GLQuadObjectChannel.create();
    GL3DDrawFunctions::BindGeometry(this,m_GLQuadObjectChannel,m_Shaders->channelshader,m_Geometries->QuadGeometry);

    is_created = true;

    return;
}

void GL3DWidget::timerEvent(QTimerEvent *e)
{
    qDebug() << "timerevent";

    if(gl_context_try)
    {
        //if(gl_context_control->tryLock())
        {
            //this->setUpdatesEnabled(true);
            //Drawing function of parent class
            this->update();
        }

    }



}

void GL3DWidget::showEvent(QShowEvent *e)
{

        qDebug() << "show";

    m_Timer.start(1, this);
}

void GL3DWidget::hideEvent(QHideEvent *e)
{
        qDebug() << "hide";

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    m_Timer.stop();
}

void GL3DWidget::resizeGL(int w, int h)
{
        qDebug() << "resize";

    m_Camera->ResizeViewPort(w,h);

    m_Textures->CreateFrameBuffers(this,w,h);

    //for now, widget resolution is render resolution
    this->Width = w;
    this->Height = h;
}

void GL3DWidget::paintGL()
{
        qDebug() << "paint";

    //Custom function
    this->Update();


   if(ReadyToDraw)
   {
       this->Onrender();
   }else
   {
       gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   }
}

void GL3DWidget::Onrender()
{

    if(m_World->m_Surface == 0 || m_World->m_SkyBox == 0 || m_World->m_WaterSurface == 0)
    {
        return;
    }
    gl->glBindFramebuffer(GL_FRAMEBUFFER, 0);
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->ChannelFramebuffer);
    GLenum chDrawBuffers[2] = {GL_COLOR_ATTACHMENT0,GL_COLOR_ATTACHMENT1};
    gl->glDrawBuffers(2, chDrawBuffers);
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    this->m_World->OnRenderBefore(this,m_Camera,m_DT);
    this->m_World->m_Surface->OnRenderBefore(this,m_World,m_Camera,m_DT);


    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->Framebuffer);
    //draw color, position/depth and normal/surfacemask
    GLenum DrawBuffers[4] = {GL_COLOR_ATTACHMENT0,GL_COLOR_ATTACHMENT1,GL_COLOR_ATTACHMENT2,GL_COLOR_ATTACHMENT3};
    gl->glDrawBuffers(4, DrawBuffers);
    if(gl->glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        qDebug() << "error with framebuffer";
    }
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    //skybox,surface and objects
    this->m_World->m_SkyBox->OnRender(this,m_World,m_Camera,m_DT);
    this->m_World->m_Surface->OnRender(this,m_World,m_Camera,m_DT);

    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->FramebufferWater);
    //copy color texture
    GLenum DrawBuffersCopy1[1] = {GL_COLOR_ATTACHMENT0};
    gl->glDrawBuffers(1, DrawBuffersCopy1);
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    this->m_Shaders->channelshader->m_program->bind();
    this->m_GLQuadObjectChannel.bind();
    m_Shaders->channelshader->ActivateTextureOn(this,m_Textures->RenderTexture,"tex_input",0);
    m_Shaders->channelshader->ActivateTextureOn(this,m_Textures->ChannelTexture,"tex_input1",1);
    m_Shaders->channelshader->ActivateTextureOn(this,m_Textures->ChannelInfoTexture,"tex_input2",2);
    m_Shaders->channelshader->ActivateTextureOn(this,m_Textures->InfoTexture,"tex_input3",3);
    this->gl->glDrawArrays(GL_TRIANGLES,0,m_Geometries->QuadGeometry->m_IndexCount);
    this->m_Shaders->channelshader->m_program->release();

    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->Framebuffer);
    //copy color texture
    GLenum DrawBuffersCopy2[1] = {GL_COLOR_ATTACHMENT0};
    gl->glDrawBuffers(1, DrawBuffersCopy2);
    //gl->glClear(GL_COLOR_BUFFER_BIT);
    this->m_Shaders->copyshader->m_program->bind();
    this->m_GLQuadObject.bind();
    gl->glDisable(GL_DEPTH_TEST);
    m_Shaders->copyshader->ActivateTextureOn(this,m_Textures->RenderTextureWater,"tex_input",0);
    this->gl->glDrawArrays(GL_TRIANGLES,0,m_Geometries->QuadGeometry->m_IndexCount);
    gl->glEnable(GL_DEPTH_TEST);

    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->Framebuffer);
    gl->glDrawBuffers(4, DrawBuffers);

    this->m_World->OnRender(this,m_Camera,m_DT);


    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->FramebufferCopy);
    //copy color texture
    GLenum DrawBuffersCopy3[1] = {GL_COLOR_ATTACHMENT0};
    gl->glDrawBuffers(1, DrawBuffersCopy3);
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    this->m_Shaders->copyshader->m_program->bind();
    this->m_GLQuadObject.bind();
    m_Shaders->copyshader->ActivateTextureOn(this,m_Textures->RenderTexture,"tex_input",0);
    this->gl->glDrawArrays(GL_TRIANGLES,0,m_Geometries->QuadGeometry->m_IndexCount);

    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->Framebuffer);
    //only add to the color buffer now
    GLenum DrawBuffersfinal[1] = {GL_COLOR_ATTACHMENT0};
    gl->glDrawBuffers(1, DrawBuffersfinal);
    gl->glClear(GL_DEPTH_BUFFER_BIT);


    //draw water
    this->m_World->m_WaterSurface->OnRenderLate(this,m_World,m_Camera,m_DT);
    this->m_World->OnRenderLate(this,m_Camera,m_DT);

    gl->glBindFramebuffer(GL_FRAMEBUFFER, m_Textures->FramebufferCopy);
    //only add to the color buffer now
    GLenum DrawBuffersfinal1[1] = {GL_COLOR_ATTACHMENT0};
    gl->glDrawBuffers(1, DrawBuffersfinal1);
    gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    this->m_Shaders->copyshader->m_program->bind();
    this->m_GLQuadObject.bind();
    m_Shaders->copyshader->ActivateTextureOn(this,m_Textures->RenderTexture,"tex_input",0);
    this->gl->glDrawArrays(GL_TRIANGLES,0,m_Geometries->QuadGeometry->m_IndexCount);

    GLuint source_buffer = m_Textures->FramebufferCopy;
    GLuint source_texture = m_Textures->RenderTextureCopy;
    GLuint target_buffer = m_Textures->Framebuffer;
    GLuint target_texture = m_Textures->RenderTexture;

    for(int i = 0; i < m_World->m_RenderPostObjectList.length();i++)
    {

        gl->glBindFramebuffer(GL_FRAMEBUFFER, target_buffer);
        //only add to the color buffer now
        GLenum DrawBuffersfinal2[1] = {GL_COLOR_ATTACHMENT0};
        gl->glDrawBuffers(1, DrawBuffersfinal2);

        m_World->m_RenderPostObjectList.at(i)->OnRenderPost(this,source_texture,target_texture,this->m_World,m_Camera,m_DT);

        GLuint buffer_temp = source_buffer;
        GLuint texture_temp = source_texture;
        source_buffer = target_buffer;
        source_texture = target_texture;
        target_buffer = buffer_temp;
        target_texture = texture_temp;

    }


    //bind original window framebuffer
    gl->glBindFramebuffer(GL_FRAMEBUFFER, 0);
    this->m_Shaders->copyshader->m_program->bind();
    this->m_GLQuadObject.bind();
    m_Shaders->copyshader->ActivateTextureOn(this,source_texture,"tex_input",0);//
    this->gl->glDrawArrays(GL_TRIANGLES,0,m_Geometries->QuadGeometry->m_IndexCount);
    this->m_Shaders->copyshader->m_program->release();




}

void GL3DWidget::Update()
{
    qint64 time_old = m_Time;
    m_Time = QDateTime::currentMSecsSinceEpoch();
    m_DT = float(m_Time - time_old)/1000.0f;
    m_Time_s = float(m_Time - m_Time_Start)/1000.0f;

    m_Camera->SetCurrentMatrix();

    this->m_World->OnUpdate(this,m_DT);

    UseInput();


}

void GL3DWidget::Close()
{
}

