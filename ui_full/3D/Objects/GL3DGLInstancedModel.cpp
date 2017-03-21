
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

#include <3D/Objects/GL3DGLInstancedModel.h>

void GL3DGLInstancedModel::SetModelHighp(QString model, double max_dist, double vert_scale, QVector3D offset)
{
    smodel_highp = model;
    max_distance_highp = max_dist;
    has_highp = true;
    highp_vert_scale = vert_scale;
}

void GL3DGLInstancedModel::SetModelLowp(QString model, double max_dist, double vert_scale, QVector3D offset)
{
    smodel_lowp = model;
    max_distance_lowp = max_dist;
    has_lowp = true;
    lowp_vert_scale = vert_scale;
}

void GL3DGLInstancedModel::SetSmoothing(bool in_smooth, double smooth_lowp, double smooth_highp)
{
    has_smooth = in_smooth;
    distance_smooth_highp = smooth_lowp;
    distance_smooth_lowp = 200;

}

void GL3DGLInstancedModel::SetMaxInstances(int max_lowp,int max_highp)
{

    max_nr_instances_lowp = max_lowp;
    max_nr_instances_highp = max_highp;

}

void GL3DGLInstancedModel::SetIncrementDistance(double increment)
{
    length_increment = increment;

}

void GL3DGLInstancedModel::SetRandomParameters(double rand_increment, double rand_scale)
{
    length_random = rand_increment;
    scale_random = rand_scale;
}

void GL3DGLInstancedModel::SetRandomRotate(bool rotate)
{
    rotate_random = rotate;

}

void GL3DGLInstancedModel::SetDraw(bool in_draw)
{
    this->draw = in_draw;
}

void GL3DGLInstancedModel::OnCreate(GL3DWidget *widget)
{

    needs_reset;

    if(has_highp)
    {
        m_Model_highp = new GL3DModel();
        m_Model_highp->Create(widget,false, true);
        m_Model_highp->LoadObjectFile(widget,smodel_highp,highp_vert_scale);
        m_Model_highp->CreateVAOs(widget);
    }

    if(has_lowp)
    {
        m_Model_lowp = new GL3DModel();
        m_Model_lowp->Create(widget,false, true);
        m_Model_lowp->LoadObjectFile(widget,smodel_lowp,lowp_vert_scale);
        m_Model_lowp->CreateVAOs(widget);
    }


}

void GL3DGLInstancedModel::OnCreatSurfaceBasedObjects(GL3DWidget * widget,GL3DWorld * world, GL3DSurface * surface, QVector3D current_pos)
{
    m_Surface = surface;

    CreateRandomData();
}


void GL3DGLInstancedModel::CreateRandomData()
{
    //seed random number generator
    srand(12345678);

    randomdata = (float *) malloc(sizeof(float) * 4 * 256);

    for(int i = 0; i < 4 * 256; i++)
    {
        randomdata[i] = ((double) rand() / double(RAND_MAX));

    }

}

void GL3DGLInstancedModel::OnUpdate(GL3DWidget * widget,GL3DWorld * world, QVector3D position,double dt)
{

}


void GL3DGLInstancedModel::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{
    if(draw)
    {

        if(has_highp && !has_lowp)
        {
            GL3DDrawFunctions::DrawModelGLInstanced(widget,m_Model_highp,camera,m_Surface,0.0,max_distance_highp,has_smooth? distance_smooth_highp:0.0,length_increment,length_random,rotate_random? 3.14159 * 2.0 :0.0,scale_random,this->randomdata);
        }else if(has_highp && has_lowp)
        {
            GL3DDrawFunctions::DrawModelGLInstanced(widget,m_Model_highp,camera,m_Surface,0.0,max_distance_highp,has_smooth? distance_smooth_highp:0.0,length_increment,length_random,rotate_random? 3.14159 * 2.0 :0.0,scale_random,this->randomdata);
        }
        if(has_lowp && !has_highp)
        {
            GL3DDrawFunctions::DrawModelGLInstanced(widget,m_Model_lowp,camera,m_Surface,0,max_distance_lowp,has_smooth? distance_smooth_lowp:0.0,length_increment,length_random,rotate_random? 3.14159 * 2.0 :0.0,scale_random,this->randomdata);

        }else if(has_lowp && has_highp)
        {
            GL3DDrawFunctions::DrawModelGLInstanced(widget,m_Model_lowp,camera,m_Surface,0.0,max_distance_lowp,has_smooth? distance_smooth_lowp:0.0,length_increment,length_random,rotate_random? 3.14159 * 2.0 :0.0,scale_random,this->randomdata);

        }
    }

}

void GL3DGLInstancedModel::OnDestroy(GL3DWidget *widget)
{


}

