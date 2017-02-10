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

    //GL3D_SHADER_SIMPLE = 0
    m_DefaultShaderList.append(LoadShaderFromText(fragmentShaderSource,vertexShaderSource,false));
    //GL3D_SHADER_SKYBOX = 1
    m_DefaultShaderList.append(LoadShaderFromText(fragmentShaderSource_SKYBOX,vertexShaderSource_SKYBOX,false));
    //GL3D_SHADER_SURFACE_TESSELATED = 2
    //m_DefaultShaderList.append(LoadShaderFromText(fragmentShaderSource_SURFACE_TESSELATED,vertexShaderSource_SURFACE_TESSELATED,false,tesselationcontrolShaderSource_SURFACE_TESSELATED,tesselationevaluationShaderSource_SURFACE_TESSELATED,geometryShaderSource_SURFACE_TESSELATED));
    m_DefaultShaderList.append(LoadShaderFromFile("surface/f.glsl","surface/v.glsl",true,"surface/tc.glsl","surface/te.glsl","surface/g.glsl"));

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
        qCritical( "shader not compiled");
        qDebug() << "shader not compiled";
    }else
    {
        qDebug() << "shader succesfully compiled";
    }
}

void GL3DShader::LoadShaderFromFile(GL3DWidget * widget,QString file_ps,QString file_vs,QString file_tcs ,QString file_tes,QString file_gs)
{

    m_program = new QOpenGLShaderProgram(widget);


    m_program->addShaderFromSourceFile(QOpenGLShader::Vertex, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_vs);
    m_program->addShaderFromSourceFile(QOpenGLShader::Fragment, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_ps);

    if(file_tcs != SHADER_NULL)
    {
        m_program->addShaderFromSourceFile(QOpenGLShader::TessellationControl, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_tcs);
    }
    if(file_tes != SHADER_NULL)
    {
        m_program->addShaderFromSourceFile(QOpenGLShader::TessellationEvaluation, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_tes);
    }
    if(file_gs != SHADER_NULL)
    {
        m_program->addShaderFromSourceFile(QOpenGLShader::Geometry, widget->m_Directory + "/" + GL3D_DIR_SHADERS + file_gs);
    }

    bool suc = m_program->link();

    if(!suc)
    {
        qCritical( "shader not compiled");
        qDebug() << "shader not compiled";
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

void GL3DShader::ClearShader()
{



}

