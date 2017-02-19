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

#include <3D/Objects/GL3DTree.h>
#include <3D/Graphics/GL3DModels.h>
#include <GL3DMath.h>

void GL3DTrees::OnCreate(GL3DWidget *widget)
{
    recieve_render = true;
    m_Tree_highp = new GL3DModel();
    m_Tree_highp->Create(widget);
    m_Tree_highp->LoadObjectFile(widget,"tree_highp/Tree1.obj",0.5);

    m_Tree_medp = new GL3DModel();
    m_Tree_medp->Create(widget);
    m_Tree_medp->LoadObjectFile(widget,"tree_medp/tree_mid.obj");

    m_Tree_lowp = new GL3DModel();
    m_Tree_lowp->Create(widget);
    m_Tree_lowp->LoadObjectFile(widget,"tree_lowp/tree_low.obj");


}

void GL3DTrees::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{
    QMatrix4x4 ModelMatrix;
    ModelMatrix.setToIdentity();

    GL3DModel * m = this->m_Tree_highp;

    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {
        widget->DrawModelGeometryWithMaterialMultipleStart(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

        for(int j = 0;j < this->Tree_Positions.length(); j++)
        {
            ModelMatrix.setToIdentity();

            float dist = (Tree_Positions.at(j) - camera->m_Position).lengthSquared();
            if(dist < 500)
            {
                ModelMatrix.translate(Tree_Positions.at(j));
                ModelMatrix.rotate(Tree_Rotation.at(j),0.0,1.0,0.0);
                widget->DrawModelGeometryWithMaterialMultiple(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);

            }
        }

        widget->DrawModelGeometryWithMaterialMultipleEnd(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);


    }

    m = this->m_Tree_medp;

    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {
        widget->DrawModelGeometryWithMaterialMultipleStart(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

        for(int j = 0;j < this->Tree_Positions.length(); j++)
        {
            ModelMatrix.setToIdentity();

            float dist = (Tree_Positions.at(j) - camera->m_Position).lengthSquared();
            if(dist > 500 && dist < 10000)
            {
                ModelMatrix.translate(Tree_Positions.at(j));
                ModelMatrix.rotate(Tree_Rotation.at(j),0.0,1.0,0.0);
                widget->DrawModelGeometryWithMaterialMultiple(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);

            }
        }

        widget->DrawModelGeometryWithMaterialMultipleEnd(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);


    }

    m = this->m_Tree_lowp;

    for(int i = 0; i < m->GLVAO_List.length(); i++)
    {
        widget->DrawModelGeometryWithMaterialMultipleStart(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

        for(int j = 0;j < this->Tree_Positions.length(); j++)
        {
            ModelMatrix.setToIdentity();

            float dist = (Tree_Positions.at(j) - camera->m_Position).lengthSquared();
            if(dist > 10000)
            {
                ModelMatrix.translate(Tree_Positions.at(j));
                ModelMatrix.rotate(Tree_Rotation.at(j),0.0,1.0,0.0);
                widget->DrawModelGeometryWithMaterialMultiple(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);

            }
        }

        widget->DrawModelGeometryWithMaterialMultipleEnd(m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

    }


}

void GL3DTrees::OnDestroy(GL3DWidget *widget)
{
    Tree_Positions.clear();
    Tree_Rotation.clear();
    Tree_Height.clear();



}

void GL3DTrees::SetTreeDistribution(GL3DSurface * s, cTMap * veg_cover, cTMap * veg_h)
{
    //seed random number generator
    srand(12345678);

    double x_increment  = 10.0;
    double y_increment  = 10.0;
    double x_rand = 1.0 * x_increment;
    double y_rand = 1.0 * y_increment;
    double cellsize = veg_cover->cellSize();
    double sizex = cellsize * veg_cover->nrCols();
    double sizey = cellsize * veg_cover->nrRows();

    int optionsx = sizex/x_increment;
    int optionsy = sizey/y_increment;

    int trees_placed = 0;

    for(int i = 0; i < optionsx; i++)
    {
        for(int j = 0; j < optionsy; j++)
        {
            int r = j * y_increment/cellsize;
            int c = i * x_increment/cellsize;

            if(r < veg_cover->nrRows() && c < veg_cover->nrCols())
            {
                if(!pcr::isMV(veg_cover->Drc))
                {

                    float rand_x= ((double) rand() / (RAND_MAX));
                    float rand_y= ((double) rand() / (RAND_MAX));
                    float rand_r= 2.0 * 3.14159 *((double) rand() / (RAND_MAX));
                    float rand_b= ((double) rand() / (RAND_MAX));

                    if((rand_b) > (1.0 - veg_cover->Drc))
                    {
                        double x = i * x_increment + rand_x * x_rand;
                        double z = j * y_increment + rand_y * y_rand;
                        Tree_Positions.append(QVector3D(x,s->GetElevation(x,z) - 0.1,z));
                        Tree_Rotation.append(rand_r);
                        Tree_Height.append(2.0);
                        trees_placed ++;
                    }
                }
            }
        }
    }

    qDebug() << "trees placed nr:" << trees_placed;
}
