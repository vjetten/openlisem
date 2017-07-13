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

#ifndef GL3DPPRAIN_H
#define GL3DPPRAIN_H

#include "ui_full/3D/Objects/GL3DObject.h"
#include "ui_full/3D/GL3DWidget.h"
#include "ui_full/3D/Graphics/GL3DShaders.h"
#include "ui_full/3D/Graphics/GL3DTextures.h"




class GL3DPPRain : public GL3DObject
{

public:
    GL3DPPRain() : GL3DObject() {
    }

    GL3DTexture * m_RainTexture;
    GL3DShader * m_RainShader;
    QOpenGLVertexArrayObject * m_Object;

    bool draw = false;

    inline void SetDraw(bool in_draw)
    {
        draw = in_draw;
    }

    float rainfall = 0;

    inline void OnCreate(GL3DWidget * widget)
    {
        this->recieve_render_post = true;
        m_RainShader = widget->m_Shaders->rainshader;
        m_RainTexture = widget->m_Textures->LoadTextureFromFile("rainfall1.png",true,true);
        m_Object = new QOpenGLVertexArrayObject();
        m_Object->create();
        GL3DDrawFunctions::BindGeometry(widget,*m_Object,m_RainShader,widget->m_Geometries->QuadGeometry);

    }

    inline void SetRainfall(float r)
    {
        rainfall = r;
    }

    inline void OnRenderPost(GL3DWidget * widget,GLuint Texture_Source, GLuint Texture_target, GL3DWorld * world, GL3DCamera* camera, double dt)
    {
        if(draw)
        {

            m_RainShader->m_program->bind();
            this->m_Object->bind();
            m_RainShader->ActivateTextureOn(widget,Texture_Source,"tex_input",0);
            m_RainShader->ActivateTextureOn(widget,m_RainTexture->m_GLTexture,"tex_rain",1);
            m_RainShader->m_program->setUniformValue("time", (float) widget->m_Time_s);
            m_RainShader->m_program->setUniformValue("rainfall", (float) rainfall);
            widget->gl->glDrawArrays(GL_TRIANGLES,0,widget->m_Geometries->QuadGeometry->m_IndexCount);
            this->m_Object->release();
        }
    }



};
#endif // GL3DPPRAIN_H
