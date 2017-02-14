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
            int texturelocation = 0;

            s->m_program->enableAttributeArray(vertexlocation);
            s->m_program->enableAttributeArray(colorlocation);

            s->m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
            s->m_program->setAttributeBuffer(colorlocation, GL_FLOAT, Vertex::colorOffset(), Vertex::ColorTupleSize, Vertex::stride());

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

void GL3DWidget::BindTexture(QOpenGLVertexArrayObject &object,GL3DShader * s, GL3DTexture * t,int index)
{
    object.bind();

}

