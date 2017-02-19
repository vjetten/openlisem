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

#include <3D/Objects/GL3DSkyBox.h>


void GL3DSkyBox::OnCreate(GL3DWidget *widget)
{
    //this->m_Shader = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SIMPLE);
    this->m_Geometry = widget->m_Geometries->LoadGeometryFromArray(vertices_SKYBOX,36,indices_SKYBOX,36);
    this->m_Shader = widget->m_Shaders->GetDefaultShader(GL3D_SHADER_SKYBOX);
    this->m_GLObject.create();
    this->m_Texture = widget->m_Textures->LoadCubeTextureFromFile("skybox_1.bmp");

    widget->BindGeometry(m_GLObject,m_Shader,m_Geometry);

}

void GL3DSkyBox::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{

    widget->gl->glDepthMask(GL_FALSE);
    widget->gl->glEnable(GL_TEXTURE_CUBE_MAP);

    QMatrix4x4 matrix;
    matrix.setToIdentity();
    m_Shader->m_program->bind();


    m_Shader->m_program->setUniformValue("Cmatrix",camera->m_CLookAtNoTranslation);
    m_Shader->m_program->setUniformValue("Mmatrix",matrix);

    int e = widget->gl->glGetError();

    widget->gl->glActiveTexture(GL_TEXTURE0);
    widget->gl->glBindTexture(GL_TEXTURE_CUBE_MAP, this->m_Texture->m_GLTexture);
    e = widget->gl->glGetError();

    
    m_Shader->m_program->setUniformValue("cubetexture",GL_TEXTURE0+ this->m_Texture->m_GLTexture);
    e = widget->gl->glGetError();


    this->m_GLObject.bind();

    widget->gl->glDrawElements(GL_TRIANGLES, this->m_Geometry->m_IndexCount, GL_UNSIGNED_INT, 0);

    e = widget->gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "3 opengl error" <<e;
    }

    this->m_GLObject.release();

    m_Shader->m_program->release();

    widget->gl->glDepthMask(GL_TRUE);
    widget->gl->glDisable(GL_TEXTURE_CUBE_MAP);
}

void GL3DSkyBox::OnDestroy(GL3DWidget *widget)
{


}
