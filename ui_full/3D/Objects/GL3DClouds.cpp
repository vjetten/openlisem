
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

#include <3D/Objects/GL3DClouds.h>




void GL3DClouds::SetDraw(bool draw)
{



}

void GL3DClouds::CreateCloudData()
{


}

void GL3DClouds::OnCreate(GL3DWidget *widget)
{
    m_Model_highp = new GL3DModel();
    m_Model_highp->Create(widget,false, true);
    m_Model_highp->LoadObjectFile(widget,"cloud_highp/cloud_highp.obj",1.0);
    m_Model_highp->BindCustomShader(widget,widget->m_Shaders->GetDefaultShader(GL3D_SHADER_CLOUDS_GLINSTANCED));
    m_Model_highp->CreateVAOs(widget);

    m_Model_lowp = new GL3DModel();
    m_Model_lowp->Create(widget,false, true);
    m_Model_lowp->LoadObjectFile(widget,"cloud_highp/cloud_highp.obj",10.0);
    m_Model_lowp->CreateVAOs(widget);


    //seed random number generator
    srand(12345678);

    randomdata = (float *) malloc(sizeof(float) * 4 * 256);

    for(int i = 0; i < 4 * 256; i++)
    {
        randomdata[i] = ((double) rand() / double(RAND_MAX));

    }

    QList<float> wave= {20,10,5,2};
    QList<float> ampl= {1.0,0.5,0.25,0.1};

    m_Texture_Perlin1 = widget->m_Textures->LoadTextureFromPerlin(256,256,wave,ampl);
    m_Texture_Perlin2 = widget->m_Textures->LoadTextureFromPerlin(256,256,wave,ampl);

}
void GL3DClouds::OnCreatSurfaceBasedObjects(GL3DWidget * widget,GL3DWorld * world, GL3DSurface * surface, QVector3D current_pos)
{

    m_Surface = surface;

}

void GL3DClouds::OnUpdate(GL3DWidget * widget,GL3DWorld * world, QVector3D Position, double dt)
{



}

void GL3DClouds::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{
    widget->gl->glDepthMask( GL_FALSE);
    widget->gl->glBlendFunc(GL_SRC_ALPHA,GL_ONE);

    //GL3DDrawFunctions::DrawModelGLInstancedCubic(widget,m_Model_highp,camera,m_Surface,64,64,256,256,4,30.0,10.0,15.0,3.14159*2.0,0.97,this->randomdata,m_Texture_Perlin1,m_Texture_Perlin2);

    GL3DDrawFunctions::DrawModelGLInstancedCubic(widget,m_Model_highp,camera,m_Surface,64,64,256,256,2,30.0,10.0,15.0,3.14159*2.0,0.85,this->randomdata,m_Texture_Perlin1,m_Texture_Perlin2,4);
    GL3DDrawFunctions::DrawModelGLInstancedCubic(widget,m_Model_highp,camera,m_Surface,64,64,256,256,2,30.0,10.0,15.0,3.14159*2.0,0.85,this->randomdata,m_Texture_Perlin1,m_Texture_Perlin2,32);

    widget->gl->glDepthMask( GL_TRUE);
    widget->gl->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void GL3DClouds::OnDestroy(GL3DWidget *widget)
{



}
