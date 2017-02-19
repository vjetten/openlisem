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
#include <QtOpenGL>
#include <QColor>
#include <QGLWidget>


void GL3DWidget::BindGeometry(QOpenGLVertexArrayObject &object,GL3DShader * s, GL3DGeometry * g)
{
    if(!g->is_2d_data)
    {
        if(g->uses_index)
        {
            object.bind();
            g->m_vertex.bind();
            g->m_index.bind();

            s->m_program->bind();

            int vertexlocation = s->m_program->attributeLocation("posAttr");
            int colorlocation = s->m_program->attributeLocation("colAttr");
            int uvlocation = s->m_program->attributeLocation("uvcAttr");
            int normallocation = s->m_program->attributeLocation("norAttr");
            int tangentlocation = s->m_program->attributeLocation("tanAttr");
            int bitangentlocation = s->m_program->attributeLocation("btaAttr");

            s->m_program->enableAttributeArray(vertexlocation);
            s->m_program->enableAttributeArray(colorlocation);

            s->m_program->enableAttributeArray(uvlocation);
            s->m_program->enableAttributeArray(normallocation);

            s->m_program->enableAttributeArray(tangentlocation);
            s->m_program->enableAttributeArray(bitangentlocation);

            s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(colorlocation, GL_FLOAT, Vertex::colorOffset(), Vertex::ColorTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(normallocation, GL_FLOAT, Vertex::normalOffset(), Vertex::NormalTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(uvlocation, GL_FLOAT, Vertex::uvOffset(), Vertex::UVTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(tangentlocation, GL_FLOAT, Vertex::tangentOffset(), Vertex::TangentTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(bitangentlocation, GL_FLOAT, Vertex::bitangentOffset(), Vertex::BiTangentTupleSize, Vertex::stride());

            object.release();
            g->m_vertex.release();
            g->m_index.release();
            s->m_program->release();
        }else
        {

            object.bind();
            g->m_vertex.bind();

            s->m_program->bind();

            int vertexlocation = s->m_program->attributeLocation("posAttr");
            int colorlocation = s->m_program->attributeLocation("colAttr");
            int texturelocation = 0;

            s->m_program->enableAttributeArray(vertexlocation);
            s->m_program->enableAttributeArray(colorlocation);

            s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(colorlocation, GL_FLOAT, Vertex::colorOffset(), Vertex::ColorTupleSize, Vertex::stride());

            object.release();
            g->m_vertex.release();
            s->m_program->release();

        }
    }else
    {
        object.bind();
        g->m_vertex.bind();

        s->m_program->bind();

        int vertexlocation = s->m_program->attributeLocation("posAttr");
        int texturelocation = 0;

        s->m_program->enableAttributeArray(vertexlocation);

        s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex2D::positionOffset(), Vertex2D::PositionTupleSize, Vertex2D::stride());

        object.release();
        g->m_vertex.release();
        s->m_program->release();


    }
}


void GL3DWidget::DrawModelObject(GL3DModel * m, GL3DCamera * camera, QMatrix4x4 ModelMatrix)
{
    /*gl->glDisable(GL_CULL_FACE);

    QMatrix4x4 matrix;
    matrix.setToIdentity();
    m->m_Shader->m_program->bind();


    m->m_Shader->m_program->setUniformValue("Light_Ambient",QVector3D(0.3,0.3,0.3));
    m->m_Shader->m_program->setUniformValue("Light_Directional",QVector3D(-1.0,-1.0,0.0));
    m->m_Shader->m_program->setUniformValue("Light_Directional_Color",QVector3D(0.7,0.7,0.7));
    m->m_Shader->m_program->setUniformValue("CPosition",camera->m_Position);
    m->m_Shader->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    m->m_Shader->m_program->setUniformValue("Mmatrix",matrix);

    int e = gl->glGetError();

    //set textures


    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {

        QOpenGLVertexArrayObject * vao = m->GLVAO_List.at(i);
        vao->bind();

        GL3DMaterial *mat = m->Material_List.at(m->Material_Pointer.at(i));

        m->m_Shader->m_program->setUniformValue("has_Texture_ka",!(mat->Texture_ka == 0));
        if(!(mat->Texture_ka == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_ka->m_GLTexture,"Texture_ka",0);
        }

        m->m_Shader->m_program->setUniformValue("has_Texture_kd",!(mat->Texture_kd == 0));
        if(!(mat->Texture_kd == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_kd->m_GLTexture,"Texture_kd",1);
        }

        m->m_Shader->m_program->setUniformValue("has_Texture_ks",!(mat->Texture_ks == 0));
        if(!(mat->Texture_ks == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_ks->m_GLTexture,"Texture_ks",2);
        }

        m->m_Shader->m_program->setUniformValue("has_Texture_ns",!(mat->Texture_ns == 0));
        if(!(mat->Texture_ns == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_ns->m_GLTexture,"Texture_ns",3);
        }

        m->m_Shader->m_program->setUniformValue("has_Texture_alpha",!(mat->Texture_alpha == 0));
        if(!(mat->Texture_alpha == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_alpha->m_GLTexture,"Texture_alpha",4);
        }

        m->m_Shader->m_program->setUniformValue("has_Texture_disp",!(mat->Texture_disp == 0));
        if(!(mat->Texture_disp == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_disp->m_GLTexture,"Texture_disp",5);
        }

        m->m_Shader->m_program->setUniformValue("has_Texture_bump",!(mat->Texture_bump == 0));
        if(!(mat->Texture_bump == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_bump->m_GLTexture,"Texture_bump",6);
        }

        m->m_Shader->m_program->setUniformValue("has_Texture_normal",!(mat->Texture_normal == 0));
        if(!(mat->Texture_normal == 0))
        {
            m->m_Shader->ActivateTextureOn(this,mat->Texture_normal->m_GLTexture,"Texture_normal",7);
        }

        m->m_Shader->m_program->setUniformValue("color_ka",mat->color_ka);
        m->m_Shader->m_program->setUniformValue("color_kd",mat->color_kd);
        m->m_Shader->m_program->setUniformValue("color_ks",mat->color_ks);
        m->m_Shader->m_program->setUniformValue("spec_power",mat->spec_power);
        m->m_Shader->m_program->setUniformValue("alpha",(float)mat->alpha);
        m->m_Shader->m_program->setUniformValue("illum",mat->illum);

        gl->glDrawElements(GL_TRIANGLES, m->Geometry_List.at(i)->m_IndexCount, GL_UNSIGNED_INT, 0);

        e = gl->glGetError();
        if( e != GL_NO_ERROR)
        {
            qDebug() << "opengl error in drawing model object" <<e;
        }

        vao->release();
    }

    m->m_Shader->m_program->release();

    gl->glEnable(GL_CULL_FACE);
    gl->glCullFace(GL_FRONT);*/

    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {
        DrawModelGeometryWithMaterialMultipleStart(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);
        DrawModelGeometryWithMaterialMultiple(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);
        DrawModelGeometryWithMaterialMultipleEnd(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);


    }



}

void GL3DWidget::DrawModelGeometryWithMaterial(GL3DGeometry * g, GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, QMatrix4x4 modelmatrix)
{


}

void GL3DWidget::DrawModelGeometryWithMaterialMultipleStart(GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera)
{
    gl->glDisable(GL_CULL_FACE);

    Shader->m_program->bind();

    Shader->m_program->setUniformValue("Light_Ambient",QVector3D(0.3,0.3,0.3));
    Shader->m_program->setUniformValue("Light_Directional",QVector3D(-1.0,-1.0,0.0));
    Shader->m_program->setUniformValue("Light_Directional_Color",QVector3D(0.7,0.7,0.7));


    int e = gl->glGetError();

    //set textures
    vao->bind();


    if((mat == 0))
    {
        mat = new GL3DMaterial();
    }

    Shader->m_program->setUniformValue("has_Texture_ka",!(mat->Texture_ka == 0));
    if(!(mat->Texture_ka == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_ka->m_GLTexture,"Texture_ka",0);
    }

    Shader->m_program->setUniformValue("has_Texture_kd",!(mat->Texture_kd == 0));
    if(!(mat->Texture_kd == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_kd->m_GLTexture,"Texture_kd",1);
    }

    Shader->m_program->setUniformValue("has_Texture_ks",!(mat->Texture_ks == 0));
    if(!(mat->Texture_ks == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_ks->m_GLTexture,"Texture_ks",2);
    }

    Shader->m_program->setUniformValue("has_Texture_ns",!(mat->Texture_ns == 0));
    if(!(mat->Texture_ns == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_ns->m_GLTexture,"Texture_ns",3);
    }

    Shader->m_program->setUniformValue("has_Texture_alpha",!(mat->Texture_alpha == 0));
    if(!(mat->Texture_alpha == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_alpha->m_GLTexture,"Texture_alpha",4);
    }

    Shader->m_program->setUniformValue("has_Texture_disp",!(mat->Texture_disp == 0));
    if(!(mat->Texture_disp == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_disp->m_GLTexture,"Texture_disp",5);
    }

    Shader->m_program->setUniformValue("has_Texture_bump",!(mat->Texture_bump == 0));
    if(!(mat->Texture_bump == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_bump->m_GLTexture,"Texture_bump",6);
    }

    Shader->m_program->setUniformValue("has_Texture_normal",!(mat->Texture_normal == 0));
    if(!(mat->Texture_normal == 0))
    {
        Shader->ActivateTextureOn(this,mat->Texture_normal->m_GLTexture,"Texture_normal",7);
    }

    Shader->m_program->setUniformValue("color_ka",mat->color_ka);
    Shader->m_program->setUniformValue("color_kd",mat->color_kd);
    Shader->m_program->setUniformValue("color_ks",mat->color_ks);
    Shader->m_program->setUniformValue("spec_power",mat->spec_power);
    Shader->m_program->setUniformValue("alpha",(float)mat->alpha);
    Shader->m_program->setUniformValue("illum",mat->illum);

}

void GL3DWidget::DrawModelGeometryWithMaterialMultiple(GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera, QMatrix4x4 ModelMatrix)
{


    Shader->m_program->setUniformValue("CPosition",camera->m_Position);
    Shader->m_program->setUniformValue("Cmatrix",camera->m_CameraMatrix);
    Shader->m_program->setUniformValue("Mmatrix",ModelMatrix);

    gl->glDrawElements(GL_TRIANGLES, g->m_IndexCount, GL_UNSIGNED_INT, 0);

    int e = gl->glGetError();
    if( e != GL_NO_ERROR)
    {
        qDebug() << "opengl error in drawing model object" <<e;
    }


}

void GL3DWidget::DrawModelGeometryWithMaterialMultipleEnd(GL3DGeometry * g,GL3DShader * Shader, QOpenGLVertexArrayObject * vao,GL3DMaterial * mat, GL3DCamera * camera)
{

    vao->release();

    Shader->m_program->release();

    gl->glEnable(GL_CULL_FACE);
    gl->glCullFace(GL_FRONT);

}

