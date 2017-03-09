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

#include <3D/Objects/GL3DBuilding.h>
#include <3D/Graphics/GL3DModels.h>
#include <GL3DMath.h>

void GL3DBuildings::OnCreate(GL3DWidget *widget)
{
    recieve_render = true;
    m_Building_h_larger = new GL3DModel();
    m_Building_h_larger->Create(widget);
    m_Building_h_larger->LoadObjectFile(widget,"house_highp_larger/Houses_Pack.obj",11);
    m_Building_h_larger->CreateVAOs(widget);

    m_Building_h_large = new GL3DModel();
    m_Building_h_large->Create(widget);
    m_Building_h_large->LoadObjectFile(widget,"house_highp_large/Final_House.obj",2.0);
    m_Building_h_large->CreateVAOs(widget);

    m_Building_h_med = new GL3DModel();
    m_Building_h_med->Create(widget);
    m_Building_h_med->LoadObjectFile(widget,"house_highp_medium/house_highp.obj",0.45);
    m_Building_h_med->CreateVAOs(widget);

    m_Building_h_small = new GL3DModel();
    m_Building_h_small->Create(widget);
    m_Building_h_small->LoadObjectFile(widget,"house_highp_small/house.obj",1.2);
    m_Building_h_small->CreateVAOs(widget);

    /*m_Building_l_larger;
    m_Building_l_large;
    m_Building_l_med;
    m_Building_l_small;*/
}

void GL3DBuildings::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{
    if(draw)
    {
        QMatrix4x4 ModelMatrix;
        ModelMatrix.setToIdentity();

        GL3DModel * m = this->m_Building_h_larger;

        for(int i = 0; i < m->GLVAO_List.length(); i++)
        {
            GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleStart(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

            for(int j = 0;j < this->Building_h_larger_Positions.length(); j++)
            {
                ModelMatrix.setToIdentity();

                    ModelMatrix.translate(this->Building_h_larger_Positions.at(j));
                    ModelMatrix.rotate(this->Building_h_larger_Rotation.at(j),0.0,1.0,0.0);
                    GL3DDrawFunctions::DrawModelGeometryWithMaterialMultiple(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);


            }

            GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleEnd(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

        }

        m = this->m_Building_h_large;

        for(int i = 0; i < m->GLVAO_List.length(); i++)
        {
            GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleStart(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

            for(int j = 0;j < this->Building_h_large_Positions.length(); j++)
            {
                ModelMatrix.setToIdentity();

                    ModelMatrix.translate(this->Building_h_large_Positions.at(j));
                    ModelMatrix.rotate(this->Building_h_large_Rotation.at(j),0.0,1.0,0.0);
                    GL3DDrawFunctions::DrawModelGeometryWithMaterialMultiple(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);


            }

            GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleEnd(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

        }

        m = this->m_Building_h_med;

        for(int i = 0; i < m->GLVAO_List.length(); i++)
        {
            GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleStart(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

             for(int j = 0;j < this->Building_h_med_Positions.length(); j++)
            {
                ModelMatrix.setToIdentity();

                    ModelMatrix.translate(this->Building_h_med_Positions.at(j));
                    ModelMatrix.rotate(this->Building_h_med_Rotation.at(j),0.0,1.0,0.0);
                    GL3DDrawFunctions::DrawModelGeometryWithMaterialMultiple(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);


            }

            GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleEnd(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

         }

          m = this->m_Building_h_small;

          for(int i = 0; i < m->GLVAO_List.length(); i++)
          {
              GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleStart(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

              for(int j = 0;j < this->Building_h_small_Positions.length(); j++)
              {
                  ModelMatrix.setToIdentity();

                     ModelMatrix.translate(this->Building_h_small_Positions.at(j));
                     ModelMatrix.rotate(this->Building_h_small_Rotation.at(j),0.0,1.0,0.0);
                     GL3DDrawFunctions::DrawModelGeometryWithMaterialMultiple(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera, ModelMatrix);


              }

              GL3DDrawFunctions::DrawModelGeometryWithMaterialMultipleEnd(widget,m->Geometry_List.at(i),m->m_Shader,m->GLVAO_List.at(i),m->Material_List.at(m->Material_Pointer.at(i)),camera);

          }

    }

}

void GL3DBuildings::OnDestroy(GL3DWidget *widget)
{



}

void GL3DBuildings::SetBuildingsDistribution(GL3DSurface * s,cTMap * bc, cTMap * temp)
{
    FOR_ROW_COL_MV(bc,bc)
    {
        temp->Drc = 0;
    }

    Building_h_larger_Positions.clear();
    Building_h_larger_Scale.clear();
    Building_h_larger_Rotation.clear();

    Building_h_large_Positions.clear();
    Building_h_large_Scale.clear();
    Building_h_large_Rotation.clear();

    Building_h_med_Positions.clear();
    Building_h_med_Scale.clear();
    Building_h_med_Rotation.clear();

    Building_h_small_Positions.clear();
    Building_h_small_Scale.clear();
    Building_h_small_Rotation.clear();

    //seed random number generator
    srand(12345678);

    m_BuildingCover = bc;

    double cs = m_BuildingCover->cellSize();

    int nrcols = m_BuildingCover->nrCols();
    int nrrows = m_BuildingCover->nrRows();

    double length_larger = 60;
    double length_large = 30;
    double length_med = 15;
    double length_small = 5;

    double minimum_spacing = 5;
    double minimum_spacing_cells = std::ceil(5.0/cs);


    double area_larger = length_larger*length_larger;
    double area_large = length_large*length_large;
    double area_med = length_med*length_med;
    double area_small = length_small*length_small;

    int max_iterations_search = std::max(1.0,area_large/cs);

    double cover_threshold = 0.5;

    int n_buildingsll = 0;
    int n_buildingsl = 0;
    int n_buildingsm = 0;
    int n_buildingss = 0;

    FOR_ROW_COL_MV(bc,bc)
    {

        //building cover here
        //start looking for an area where a building might be
        if(bc->Drc > cover_threshold && temp->Drc == 0)
        {
            int i = 1;
            bool bbreak = false;

            //untill max iterations are reached, increase square size and check connnected building area
            for(i = 2; i < max_iterations_search; i++)
            {

                for(int j = 0; j < i; j++)
                {
                    int r2 = r + (j+1);
                    int c2 = c + (i-1);

                    if(r2 < nrrows && c2 < nrcols)
                    {
                        if(temp->data[r2][c2] == 1 || bc->data[r2][c2] < cover_threshold || pcr::isMV(bc->data[r2][c2]))
                        {
                            bbreak = true;
                            i--;
                            break;

                        }
                    }
                }
                for(int j = 0; j < i-1; j++)
                {
                    int r2 = r + (i-1);
                    int c2 = c + (j+1);

                    if(r2 < nrrows && c2 < nrcols)
                    {
                        if(temp->data[r2][c2] == 1 || bc->data[r2][c2] < cover_threshold || pcr::isMV(bc->data[r2][c2]))
                        {
                            bbreak = true;
                            i--;
                            break;
                        }
                    }
                }
                if(bbreak)
                {
                    break;
                }

            }

            //now check the end resulting value of i, to get connected size of building cover

            double area = (double(i)) *(double(i)) * cs*cs;

            float rand_x= ((double) rand() / (RAND_MAX));
            float rand_y= ((double) rand() / (RAND_MAX));
            float rand_r= 2.0 * 3.14159 *((double) rand() / (RAND_MAX));

            bool placed = false;

            if(area >= area_larger )
            {
                placed = true;
                double x1 = ((double)(c)) * cs ;
                double x2 = ((double)(c)) * cs + length_larger;
                double z1 = ((double)(r)) * cs ;
                double z2 = ((double)(r)) * cs + length_larger;

                double y = std::min(std::min(std::min(s->GetElevation(x1,z1),s->GetElevation(x1,z2)),s->GetElevation(x2,z1)),s->GetElevation(x2,z2));

                this->Building_h_larger_Positions.append(QVector3D((x1+x2)/2.0,y,(z1 + z2)/2.0));
                this->Building_h_larger_Rotation.append(0.0);
                this->Building_h_larger_Scale.append(QVector3D(1.0,1.0,1.0));

                n_buildingsll ++;

            }else if(area >= area_large )
            {
                placed = true;
                double x1 = ((double)(c)) * cs ;
                double x2 = ((double)(c)) * cs + length_large;
                double z1 = ((double)(r)) * cs ;
                double z2 = ((double)(r)) * cs + length_large;

                double y = std::min(std::min(std::min(s->GetElevation(x1,z1),s->GetElevation(x1,z2)),s->GetElevation(x2,z1)),s->GetElevation(x2,z2));

                this->Building_h_large_Positions.append(QVector3D((x1+x2)/2.0,y,(z1 + z2)/2.0));
                this->Building_h_large_Rotation.append(0.0);
                this->Building_h_large_Scale.append(QVector3D(1.0,1.0,1.0));

                n_buildingsl ++;

            }else if(area >= area_med)
            {
                placed = true;
                double x1 = ((double)(c)) * cs ;
                double x2 = ((double)(c)) * cs + length_med;
                double z1 = ((double)(r)) * cs ;
                double z2 = ((double)(r)) * cs + length_med;

                double y = std::min(std::min(std::min(s->GetElevation(x1,z1),s->GetElevation(x1,z2)),s->GetElevation(x2,z1)),s->GetElevation(x2,z2));

                this->Building_h_med_Positions.append(QVector3D((x1+x2)/2.0,y,(z1 + z2)/2.0));
                this->Building_h_med_Rotation.append(0.0);
                this->Building_h_med_Scale.append(QVector3D(1.0,1.0,1.0));

                n_buildingsm ++;

            }else if(area >= area_small)
            {
                placed = true;
                double x1 = ((double)(c)) * cs ;
                double x2 = ((double)(c)) * cs + length_small;
                double z1 = ((double)(r)) * cs ;
                double z2 = ((double)(r)) * cs + length_small;

                double y = std::min(std::min(std::min(s->GetElevation(x1,z1),s->GetElevation(x1,z2)),s->GetElevation(x2,z1)),s->GetElevation(x2,z2));

                this->Building_h_small_Positions.append(QVector3D((x1+x2)/2.0,y,(z1 + z2)/2.0));
                this->Building_h_small_Rotation.append(0.0);
                this->Building_h_small_Scale.append(QVector3D(1.0,1.0,1.0));

                n_buildingss ++;

            }

            //depending on the required spacing, increase i when marking building presence

            int k = i + minimum_spacing_cells;
            if(placed)
            {
                for(i = 1; i < k; i++)
                {
                    for(int j = 0; j < i; j++)
                    {
                        int r2 = r + (j+1);
                        int c2 = c + (i-1);

                        if(r2 < nrrows && c2 < nrcols)
                        {
                            temp->data[r2][c2] = 1;
                        }
                    }
                    for(int j = 0; j < i-1; j++)
                    {
                        int r2 = r + (i-1);
                        int c2 = c + (j+1);

                        if(r2 < nrrows && c2 < nrcols)
                        {
                            temp->data[r2][c2] = 1;
                        }
                    }
                }
            }

        }

    }

    qDebug() << "buildings place Nr: " << n_buildingsll << n_buildingsl << n_buildingsm << n_buildingss;

}
