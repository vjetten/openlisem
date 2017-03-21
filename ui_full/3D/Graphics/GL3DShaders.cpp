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

#include <3D/Graphics/GL3DShaders.h>



GL3DShader * GL3DShaders::LoadShaderFromText(const char * text_ps,const char * text_vs, bool add,const char * text_tcs,const char * text_tes,const char * text_gs)
{
    GL3DShader * sh = new GL3DShader();

    sh->LoadShaderFromText(m_Widget,text_ps,text_vs,text_tcs,text_tes,text_gs);
    if(add)
    {
        this->m_LoadedShaderList.append(sh);
    }
    return sh;
}

GL3DShader * GL3DShaders::LoadShaderFromFile(QString file_ps, QString file_vs, bool add, QString file_tcs, QString file_tes,QString file_gs)
{
    GL3DShader * sh = new GL3DShader();

    sh->LoadShaderFromFile(m_Widget,file_ps,file_vs,file_tcs,file_tes,file_gs);
    if(add)
    {
        this->m_LoadedShaderList.append(sh);
    }
    return sh;


}


GL3DShader * GL3DShaders::GetDefaultShader(int index)
{
    if((index > -1) && (index < m_DefaultShaderList.length()))
    {

        return m_DefaultShaderList.at(index);
    }else
    {
        return 0;
    }

}

void GL3DShaders::Create(GL3DWidget * widget)
{

    m_Widget = widget;
    m_DefaultShaderList.clear();
    m_LoadedShaderList.clear();

    rainshader = LoadShaderFromFile("postprocess/rain_f.glsl","postprocess/rain_v.glsl",false);

    channelshader = LoadShaderFromFile("postprocess/channel_f.glsl","postprocess/channel_v.glsl",false);

    copyshader = LoadShaderFromFile("postprocess/copy_f.glsl","postprocess/copy_v.glsl",false);

    blendshader = LoadShaderFromFile("postprocess/blend_f.glsl","postprocess/blend_v.glsl",false);

    //GL3D_SHADER_SIMPLE = 0
    m_DefaultShaderList.append(LoadShaderFromText(fragmentShaderSource,vertexShaderSource,false));
    //GL3D_SHADER_SKYBOX = 1
    m_DefaultShaderList.append(LoadShaderFromFile("sky/sky_f.glsl","sky/sky_v.glsl",false));
    //GL3D_SHADER_SURFACE_TESSELATED = 2
    m_DefaultShaderList.append(LoadShaderFromFile("surface/f.glsl","surface/v.glsl",false,"surface/tc.glsl","surface/te.glsl","surface/g.glsl"));
    //GL3D_SHADER_SURFACE_FLOW = 3
    m_DefaultShaderList.append(LoadShaderFromFile("flow/f.glsl","flow/v.glsl",false,"flow/tc.glsl","flow/te.glsl","flow/g.glsl"));
    //GL3D_SHADER_MODEL = 4
    m_DefaultShaderList.append(LoadShaderFromFile("objects/object_f.glsl","objects/object_v.glsl",false));
    //GL3D_SHADER_MODEL_INSTANCED = 5
    m_DefaultShaderList.append(LoadShaderFromFile("objects/object_i_f.glsl","objects/object_i_v.glsl",false));
    //GL3D_SHADER_MODEL_INSTANCED = 6
    m_DefaultShaderList.append(LoadShaderFromFile("objects/object_gli_f.glsl","objects/object_gli_v.glsl",false));
    //GL3D_SHADER_ROADS = 7
    m_DefaultShaderList.append(LoadShaderFromFile("roads/roads_f.glsl","roads/roads_v.glsl",false));
    //GL3D_SHADER_CHANNEL = 8
    m_DefaultShaderList.append(LoadShaderFromFile("channel/channel_f.glsl","channel/channel_v.glsl",false,"","","channel/channel_g.glsl"));
    //GL3D_SHADER_CLOUDS_GLINSTANCED = 9
    m_DefaultShaderList.append(LoadShaderFromFile("clouds/clouds_gli_f.glsl","clouds/clouds_gli_v.glsl",false));

}

void GL3DShaders::ClearUnused()
{


}

void GL3DShaders::Destroy()
{



}


void GL3DShader::LoadShaderFromText(GL3DWidget * widget,const char * text_ps,const char * text_vs,const char * text_tcs ,const char * text_tes,const char * text_gs)
{

    m_program = new QOpenGLShaderProgram(widget);
    m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, text_vs);
    m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, text_ps);

    if(text_tcs != SHADER_NULL)
    {
        m_program->addShaderFromSourceCode(QOpenGLShader::TessellationControl, text_tcs);
    }
    if(text_tes != SHADER_NULL)
    {
        m_program->addShaderFromSourceCode(QOpenGLShader::TessellationEvaluation, text_tes);
    }
    if(text_gs != SHADER_NULL)
    {
        m_program->addShaderFromSourceCode(QOpenGLShader::Geometry, text_gs);
    }



    bool suc = m_program->link();

    if(!suc)
    {
        qDebug() << "shader not compiled ";
    }else
    {
        qDebug() << "shader succesfully compiled ";
    }
}

void GL3DShader::LoadShaderFromFile(GL3DWidget * widget,QString file_ps,QString file_vs,QString file_tcs ,QString file_tes,QString file_gs)
{
    m_program = new QOpenGLShaderProgram(widget);

    m_program->addShaderFromSourceFile(QOpenGLShader::Vertex, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_vs);
    m_program->addShaderFromSourceFile(QOpenGLShader::Fragment, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_ps);

    if(file_tcs.length() != 0)
    {
        m_program->addShaderFromSourceFile(QOpenGLShader::TessellationControl, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_tcs);
    }
    if(file_tes.length() != 0)
    {
        m_program->addShaderFromSourceFile(QOpenGLShader::TessellationEvaluation, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_tes);
    }
    if(file_gs.length() != 0)
    {
        m_program->addShaderFromSourceFile(QOpenGLShader::Geometry, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_gs);
    }

    bool suc = m_program->link();

    if(!suc)
    {
        qDebug() << "shader not compiled " << file_vs << " " << file_ps;
    }else
    {
        qDebug() << "shader succesfully compiled";
    }
}

void GL3DShader::ActivateTextureOn(GL3DWidget * widget, GL3DTexture *t, const char * name, int slot)
{
    m_program->setUniformValue(name, slot);
    widget->gl->glActiveTexture(GL_TEXTURE0 + slot);
    widget->gl->glBindTexture(GL_TEXTURE_2D, t->m_GLTexture);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

}

void GL3DShader::ActivateTextureOn(GL3DWidget * widget, GLuint t, const char * name, int slot)
{
    m_program->setUniformValue(name, slot);
    widget->gl->glActiveTexture(GL_TEXTURE0 + slot);
    widget->gl->glBindTexture(GL_TEXTURE_2D, t);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    widget->gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

}

void GL3DShader::ActivateCubeTextureOn(GL3DWidget * widget, GLuint t, const char * name, int slot)
{
    m_program->setUniformValue(name, slot);
    widget->gl->glActiveTexture(GL_TEXTURE0 + slot);
    widget->gl->glBindTexture(GL_TEXTURE_CUBE_MAP, t);
    widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_REPEAT);
    widget->gl->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_REPEAT);

}


void GL3DShader::ClearShader()
{



}

