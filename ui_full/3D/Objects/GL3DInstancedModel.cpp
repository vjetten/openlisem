
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

#include <3D/Objects/GL3DInstancedModel.h>

void GL3DInstancedModel::SetModelHighp(QString model, double max_dist, double vert_scale, QVector3D offset)
{
    smodel_highp = model;
    max_distance_highp = max_dist;
    has_highp = true;
    highp_vert_scale = vert_scale;
}

void GL3DInstancedModel::SetModelLowp(QString model, double max_dist, double vert_scale, QVector3D offset)
{
    smodel_lowp = model;
    max_distance_lowp = max_dist;
    has_lowp = true;
    lowp_vert_scale = vert_scale;
}

void GL3DInstancedModel::SetSmoothing(bool in_smooth, double smooth_lowp, double smooth_highp)
{
    has_smooth = in_smooth;
    distance_smooth_highp = smooth_lowp;
    distance_smooth_lowp = 200;

}

void GL3DInstancedModel::SetMaxInstances(int max_lowp,int max_highp)
{

    max_nr_instances_lowp = max_lowp;
    max_nr_instances_highp = max_highp;

}

void GL3DInstancedModel::SetLocations(QList<QVector3D> positions)
{
    Actual_Positions = positions;

}

void GL3DInstancedModel::SetIncrementDistance(double increment)
{
    length_increment = increment;

}
void GL3DInstancedModel::SetProbability(cTMap * prob_map)
{

    m_ProbMap = prob_map;
}


void GL3DInstancedModel::SetRandomParameters(double rand_increment, double rand_scale)
{
    length_random = rand_increment;
    scale_random = rand_scale;
}

void GL3DInstancedModel::SetRandomRotate(bool rotate)
{
    rotate_random = rotate;

}

void GL3DInstancedModel::SetDraw(bool in_draw)
{
    this->draw = in_draw;
}

void GL3DInstancedModel::OnCreate(GL3DWidget *widget)
{

    needs_reset;
    old_row = 0;
    old_col = 0;

    if(has_highp)
    {
        m_Model_highp = new GL3DModel();
        m_Model_highp->Create(widget,true);
        m_Model_highp->LoadObjectFile(widget,smodel_highp,highp_vert_scale);
        m_Model_highp->CreateMatrixBuffer(this->max_nr_instances_highp);
        m_Model_highp->CreateVAOs(widget);
    }

    if(has_lowp)
    {
        m_Model_lowp = new GL3DModel();
        m_Model_lowp->Create(widget,true);
        m_Model_lowp->LoadObjectFile(widget,smodel_lowp,lowp_vert_scale);
        m_Model_lowp->CreateMatrixBuffer(this->max_nr_instances_lowp);
        m_Model_lowp->CreateVAOs(widget);
    }


}

void GL3DInstancedModel::OnCreatSurfaceBasedObjects(GL3DWidget * widget,GL3DWorld * world, GL3DSurface * surface, QVector3D current_pos)
{
    m_Surface = surface;

    nr_Cols = ((int) floor(m_Surface->m_XExtent/length_increment));
    nr_Rows = ((int) floor(m_Surface->m_ZExtent/length_increment));

    CreateRandomData();

    UpdateBuffer(widget,world,current_pos);
}

int GL3DInstancedModel::GetRandomDataIndexAt(int r, int c)
{
    if(r < 0 || c < 0)
    {
        return 0;
    }
    return ((r * nr_Cols + c) % max_nr_random);

}

void GL3DInstancedModel::CreateRandomData()
{
    //seed random number generator
    srand(12345678);

    Random_Positions.clear();
    Random_Rotation.clear();
    Random_Height.clear();
    Random_Alpha.clear();
    Random_Random.clear();

    for(int i = 0; i < max_nr_random; i++)
    {
        float rand_x= ((double) rand() / (RAND_MAX));
        float rand_y= ((double) rand() / (RAND_MAX));
        float rand_s= ((double) rand() / (RAND_MAX));
        float rand_r= 2.0 * 3.14159 *((double) rand() / (RAND_MAX));
        float rand_b= ((double) rand() / (RAND_MAX));

        Random_Positions.append(QVector3D(rand_x * length_random,0.0,rand_y * length_random));
        Random_Rotation.append(rand_r);
        Random_Height.append(1.0 + (rand_s *scale_random * 2.0 - scale_random));
        Random_Alpha.append(rand_b);
        Random_Random.append(rand_b);

    }



}

void GL3DInstancedModel::CreateModelDistribution(GL3DWidget * widget,GL3DSurface * s,cTMap * veg_cover, cTMap * veg_h)
{


}


void GL3DInstancedModel::UpdateBuffer(GL3DWidget * widget,GL3DWorld * world, QVector3D pos)
{

    bool redraw_buffer =false;

    if(has_smooth)
    {
        current_cols = floor(std::max(max_distance_lowp + distance_smooth_lowp*0.5,max_distance_highp + distance_smooth_highp*0.5)/length_increment);
    }else
    {
        current_cols = floor(std::max(max_distance_lowp,max_distance_highp)/length_increment);
    }
    if(current_cols % 2 == 1)
    {
        current_cols += 1;
    }
    current_cols = current_cols + 1;
    current_rows = current_cols;

    current_col = (int) floor(pos.x() / length_increment);
    current_row = (int) floor(pos.z()/  length_increment);

    int dr = current_row - old_row;
    int dc = current_col - old_col;

    //if too much movement, it is pointless to update tree list from old data
    //then just do a hard reset
    if(abs(dr) + abs(dc) > 15)
    {
        needs_reset = true;
    }

    if(needs_reset)
    {
        Current_Positions.clear();
        Current_Rotation.clear();
        Current_Height.clear();
        Current_Alpha.clear();
        Current_Random.clear();
        Current_b_highp.clear();
        Current_b_lowp.clear();

        for(int c = current_col - 0.5 * (current_cols-1); c < 1+current_col + 0.5 * (current_cols-1); c++)
        {
            for(int r = current_row - 0.5 * (current_rows-1); r < 1+current_row + 0.5 * (current_rows-1); r++)
            {


                int rand_index = GetRandomDataIndexAt(r,c);
                double treecover = MapMath::GetValueAt(m_ProbMap,c * length_increment,r * length_increment);

                if(Random_Random.at(rand_index) < treecover)
                {
                    QVector3D pos_rand = Random_Positions.at(rand_index);
                    QVector3D location = QVector3D(c * length_increment + length_increment * 2.0*(pos_rand.x()-0.5),0.0,r * length_increment + length_increment * 2.0*(pos_rand.z()-0.5));
                    location.setY(m_Surface->GetElevation(location.x(),location.z()));
                    Current_Positions.append(highp_offset+location);
                    Current_RowCol.append(QVector2D(r,c));

                    double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                    if(dist < max_distance_highp*max_distance_highp)
                    {
                        Current_b_highp.append(true);
                        Current_b_lowp.append(false);
                    }else if(dist < max_distance_lowp*max_distance_lowp)
                    {
                        Current_b_highp.append(false);
                        Current_b_lowp.append(true);
                    }else
                    {
                        Current_b_highp.append(false);
                        Current_b_lowp.append(false);
                    }

                }else
                {
                    Current_Positions.append(QVector3D(0.0,0.0,0.0));
                    Current_RowCol.append(QVector2D(-1,-1));

                    Current_b_highp.append(false);
                    Current_b_lowp.append(false);

                }


            }

        }


        old_col = current_col;
        old_row = current_row;


        needs_reset = false;
        redraw_buffer = true;

    }else if(abs(dr) + abs(dc) > 0)
    {
        //in this method of rewriting the list of currently drawn trees, several things could be improved
        //when seperately moving x and y directions of movement, the corners are not done efficiently
        //secondly, this method would be much faster if implemented using a normal dynamicly allocated array

        for(int i = 0; i < abs(dr); i++)
        {
            //move everything top or bottom direction
            if(dr > 0)
            {

                //delete first
                Current_RowCol.removeFirst();
                Current_Positions.removeFirst();
                Current_b_highp.removeFirst();
                Current_b_lowp.removeFirst();

                //now change the bottom locations to the top locations
                for(int j = 1; j < current_cols+1; j++)
                {
                    int index = j * current_rows;
                    int r = old_row + 0.5 * (current_rows - 1)+1;
                    int c = old_col - 0.5 * (current_cols - 1) + (j-1);

                    int rand_index = GetRandomDataIndexAt(r,c);
                    double treecover = MapMath::GetValueAt(m_ProbMap,c * length_increment,r * length_increment);

                    QVector3D pos_rand;
                    QVector3D location;
                    QVector2D row_col_loc;

                    if(Random_Random.at(rand_index) < treecover)
                    {
                        pos_rand = Random_Positions.at(rand_index);
                        location = QVector3D(c * length_increment + length_increment * 2.0*(pos_rand.x()-0.5),0.0,r * length_increment + length_increment * 2.0*(pos_rand.z()-0.5));
                        location.setY(m_Surface->GetElevation(location.x(),location.z()));
                        row_col_loc = QVector2D(r,c);
                    }else
                    {
                        location = QVector3D(0.0,0.0,0.0);
                        row_col_loc = QVector2D(-1,-1);

                    }

                    if(index < Current_Positions.length())
                    {
                        Current_Positions.replace(index,highp_offset+location);
                        Current_RowCol.replace(index,row_col_loc);
                        double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                        if(dist < max_distance_highp)
                        {
                            Current_b_highp.replace(index,true);
                            Current_b_lowp.replace(index,false);
                        }else if(dist < max_distance_lowp)
                        {
                            Current_b_highp.replace(index,false);
                            Current_b_lowp.replace(index,true);
                        }else
                        {
                            Current_b_highp.replace(index,false);
                            Current_b_lowp.replace(index,false);
                        }
                    }else
                    {
                        Current_Positions.append(highp_offset+location);
                        Current_RowCol.append(row_col_loc);
                        double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                        if(dist < max_distance_highp)
                        {
                            Current_b_highp.append(true);
                            Current_b_lowp.append(false);
                        }else if(dist < max_distance_lowp)
                        {
                            Current_b_highp.append(false);
                            Current_b_lowp.append(true);
                        }else
                        {
                            Current_b_highp.append(false);
                            Current_b_lowp.append(false);
                        }
                    }


                }

                old_row = old_row+1;

            }
            //move everything top or bottom direction
            if(dr < 0)
            {

                //delete last
                Current_RowCol.removeLast();
                Current_Positions.removeLast();
                Current_b_highp.removeLast();
                Current_b_lowp.removeLast();

                //now change the bottom locations to the top locations
                for(int j = 1; j < current_cols+1; j++)
                {

                    int index = (j * current_rows);

                    int r = old_row - 0.5 * (current_rows - 1)-1;
                    int c = old_col - 0.5 * (current_cols - 1) + (j-1);

                    int rand_index = GetRandomDataIndexAt(r,c);
                    double treecover = MapMath::GetValueAt(m_ProbMap,c * length_increment,r * length_increment);

                    QVector3D pos_rand;
                    QVector3D location;
                    QVector2D row_col_loc;

                    if(Random_Random.at(rand_index) < treecover)
                    {
                        pos_rand = Random_Positions.at(rand_index);
                        location = QVector3D(c * length_increment + length_increment * 2.0*(pos_rand.x()-0.5),0.0,r * length_increment + length_increment * 2.0*(pos_rand.z()-0.5));
                        location.setY(m_Surface->GetElevation(location.x(),location.z()));
                        row_col_loc = QVector2D(r,c);
                    }else
                    {
                        location = QVector3D(0.0,0.0,0.0);
                        row_col_loc = QVector2D(-1,-1);

                    }

                    if(index > 0 && index < Current_Positions.length())
                    {
                        Current_Positions.replace(index,highp_offset+location);
                        Current_RowCol.replace(index,row_col_loc);
                        double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                        if(dist < max_distance_highp)
                        {
                            Current_b_highp.replace(index,true);
                            Current_b_lowp.replace(index,false);
                        }else if(dist < max_distance_lowp)
                        {
                            Current_b_highp.replace(index,false);
                            Current_b_lowp.replace(index,true);
                        }else
                        {
                            Current_b_highp.replace(index,false);
                            Current_b_lowp.replace(index,false);
                        }
                    }else
                    {
                        Current_Positions.prepend(highp_offset+location);
                        Current_RowCol.prepend(row_col_loc);
                        double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                        if(dist < max_distance_highp)
                        {
                            Current_b_highp.prepend(true);
                            Current_b_lowp.prepend(false);
                        }else if(dist < max_distance_lowp)
                        {
                            Current_b_highp.prepend(false);
                            Current_b_lowp.prepend(true);
                        }else
                        {
                            Current_b_highp.prepend(false);
                            Current_b_lowp.prepend(false);
                        }
                    }


                }

                int r = old_row - 1 - 0.5 * (current_rows - 1);
                int c = old_col - 1 - 0.5 * (current_cols - 1);

                int rand_index = GetRandomDataIndexAt(r,c);
                double treecover = MapMath::GetValueAt(m_ProbMap,r * length_increment,c * length_increment);

                QVector3D pos_rand;
                QVector3D location;
                QVector2D row_col_loc;

                if(Random_Random.at(rand_index) < treecover)
                {
                    pos_rand = Random_Positions.at(rand_index);
                    location = QVector3D(c * length_increment + length_increment * 2.0*(pos_rand.x()-0.5),0.0,r * length_increment + length_increment * 2.0*(pos_rand.z()-0.5));
                    location.setY(m_Surface->GetElevation(location.x(),location.z()));
                    row_col_loc = QVector2D(r,c);
                }else
                {
                    location = QVector3D(0.0,0.0,0.0);
                    row_col_loc = QVector2D(-1,-1);

                }

                Current_Positions.prepend(highp_offset+location);
                Current_RowCol.prepend(row_col_loc);

                double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                if(dist < max_distance_highp)
                {
                    Current_b_highp.prepend(true);
                    Current_b_lowp.prepend(false);
                }else if(dist < max_distance_lowp)
                {
                    Current_b_highp.prepend(false);
                    Current_b_lowp.prepend(true);
                }else
                {
                    Current_b_highp.prepend(false);
                    Current_b_lowp.prepend(false);
                }

                old_row = old_row-1;
            }


        }


        for(int i = 0; i < abs(dc); i++)
        {
            //move everything top or bottom direction
            if(dc > 0)
            {

                //delete first colum
                for(int j = 0; j < current_rows; j++)
                {
                    Current_RowCol.removeFirst();
                    Current_Positions.removeFirst();
                    Current_b_highp.removeFirst();
                    Current_b_lowp.removeFirst();
                }

                //now change the bottom locations to the top locations
                for(int j = 0; j < current_rows; j++)
                {
                    int r = old_row - 0.5 * (current_rows - 1)+ (j);
                    int c = old_col + 0.5 * (current_cols - 1) +1;

                    int rand_index = GetRandomDataIndexAt(r,c);
                    double treecover = MapMath::GetValueAt(m_ProbMap,c * length_increment,r * length_increment);

                    QVector3D pos_rand;
                    QVector3D location;
                    QVector2D row_col_loc;

                    if(Random_Random.at(rand_index) < treecover)
                    {
                        pos_rand = Random_Positions.at(rand_index);
                        location = QVector3D(c * length_increment + length_increment * 2.0*(pos_rand.x()-0.5),0.0,r * length_increment + length_increment * 2.0*(pos_rand.z()-0.5));
                        location.setY(m_Surface->GetElevation(location.x(),location.z()));
                        row_col_loc = QVector2D(r,c);
                    }else
                    {
                        location = QVector3D(0.0,0.0,0.0);
                        row_col_loc = QVector2D(-1,-1);

                    }

                    Current_Positions.append(highp_offset+location);
                    Current_RowCol.append(row_col_loc);

                    double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                    if(dist < max_distance_highp)
                    {
                        Current_b_highp.append(true);
                        Current_b_lowp.append(false);
                    }else if(dist < max_distance_lowp)
                    {
                        Current_b_highp.append(false);
                        Current_b_lowp.append(true);
                    }else
                    {
                        Current_b_highp.append(false);
                        Current_b_lowp.append(false);
                    }

                }
                old_col = old_col+1;
            }
            //move everything top or bottom direction
            if(dc < 0)
            {

                //delete last colum
                for(int j = 0; j < current_rows; j++)
                {
                    Current_RowCol.removeLast();
                    Current_Positions.removeLast();
                    Current_b_highp.removeLast();
                    Current_b_lowp.removeLast();
                }

                //now change the bottom locations to the top locations
                for(int j = 0; j < current_rows; j++)
                {
                    int r = old_row - 0.5 * (current_rows - 1)+ (j);
                    int c = old_col - 0.5 * (current_cols - 1) - 1;

                    int rand_index = GetRandomDataIndexAt(r,c);
                    double treecover = MapMath::GetValueAt(m_ProbMap,c * length_increment,r * length_increment);

                    QVector3D pos_rand;
                    QVector3D location;
                    QVector2D row_col_loc;

                    if(Random_Random.at(rand_index) < treecover)
                    {
                        pos_rand = Random_Positions.at(rand_index);
                        location = QVector3D(c * length_increment + length_increment * 2.0*(pos_rand.x()-0.5),0.0,r * length_increment + length_increment * 2.0*(pos_rand.z()-0.5));
                        location.setY(m_Surface->GetElevation(location.x(),location.z()));
                        row_col_loc = QVector2D(r,c);
                    }else
                    {
                        location = QVector3D(0.0,0.0,0.0);
                        row_col_loc = QVector2D(-1,-1);
                    }

                    Current_Positions.prepend(highp_offset+location);
                    Current_RowCol.prepend(row_col_loc);

                    double dist = QVector3D(location.x() - pos.x(),location.y()-pos.y(),location.z()-pos.z()).lengthSquared();
                    if(dist < max_distance_highp)
                    {
                        Current_b_highp.prepend(true);
                        Current_b_lowp.prepend(false);
                    }else if(dist < max_distance_lowp)
                    {
                        Current_b_highp.prepend(false);
                        Current_b_lowp.prepend(true);
                    }else
                    {
                        Current_b_highp.prepend(false);
                        Current_b_lowp.prepend(false);
                    }

                }

                old_col = old_col-1;
            }

        }

        redraw_buffer = true;
    }

    if(redraw_buffer)
    {
        this->m_Model_highp->BindMatrixBuffer();
        this->m_Model_highp->WriteToMatrixBuffer(&Current_b_highp,&Current_Positions,0,0);
        this->m_Model_highp->ReleaseMatrixBuffer();

        this->m_Model_lowp->BindMatrixBuffer();
        this->m_Model_lowp->WriteToMatrixBuffer(&Current_b_lowp,&Current_Positions,0,0);
        this->m_Model_lowp->ReleaseMatrixBuffer();
    }

}


void GL3DInstancedModel::OnUpdate(GL3DWidget * widget,GL3DWorld * world, QVector3D position,double dt)
{

    if(draw)
    {
        UpdateBuffer(widget,world,position);
    }

}


void GL3DInstancedModel::OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)
{
    if(draw)
    {
        GL3DDrawFunctions::DrawModelInstanced(widget,m_Model_highp,camera);
        GL3DDrawFunctions::DrawModelInstanced(widget,m_Model_lowp,camera);
    }

}

void GL3DInstancedModel::OnDestroy(GL3DWidget *widget)
{


}

