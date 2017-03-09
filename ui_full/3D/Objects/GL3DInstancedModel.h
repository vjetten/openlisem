

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
#ifndef GL3DINSTANCEDMODEL_H
#define GL3DINSTANCEDMODEL_H

#include <3D/Objects/GL3DObject.h>
#include <3D/World/GL3DSurface.h>


class GL3DInstancedModel : public GL3DObject
{

public:
    GL3DInstancedModel() : GL3DObject()
    {
            recieve_render = true;
    }

    bool needs_reset;

    bool draw = true;

    double max_distance_highp = 75;
    double max_distance_lowp = 500;

    bool has_smooth = true;
    double distance_smooth_highp = 50;
    double distance_smooth_lowp = 200;

    double length_increment;
    double length_random;
    double scale_random;
    bool rotate_random;

    bool has_highp;
    bool has_lowp;
    QString smodel_lowp;
    QString smodel_highp;

    int current_nr_instances = 0;
    int max_nr_instances_lowp = 5000;
    int max_nr_instances_highp = 100;

    int max_nr_random = 1000;

    bool has_actual;

    QList<QVector3D> Actual_Positions;

    QList<QVector3D> Random_Positions;
    QList<float> Random_Rotation;
    QList<float> Random_Height;
    QList<float> Random_Alpha;
    QList<float> Random_Random;

    int current_cols;
    int current_rows;

    int current_col;
    int current_row;

    int old_col;
    int old_row;

    int nr_Cols;
    int nr_Rows;

   QList<bool> Current_b_highp;
   QList<bool> Current_b_lowp;

    QList<QVector2D> Original_RowCol;
    QList<QVector2D> Current_RowCol;
    QList<QVector3D> Current_Positions;
    QList<float> Current_Rotation;
    QList<float> Current_Height;
    QList<float> Current_Alpha;
    QList<float> Current_Random;

    cTMap * m_ProbMap;

    GL3DModel * m_Model_highp;
    GL3DModel * m_Model_lowp;

    double highp_vert_scale = 1.0;
    double lowp_vert_scale = 1.0;

    QVector3D highp_offset = QVector3D(0.0,0.0,0.0);
    QVector3D lowp_offset = QVector3D(0.0,0.0,0.0);

    GL3DSurface * m_Surface;

    void SetModelHighp(QString model, double max_dist, double vert_scale = 1.0, QVector3D offset = QVector3D(0.0,0.0,0.0));
    void SetModelLowp(QString model, double max_dist, double vert_scale = 1.0, QVector3D offset = QVector3D(0.0,0.0,0.0));

    void SetMaxInstances(int max_lowp,int max_highp);

    void SetSmoothing(bool smooth, double smooth_lowp, double smooth_highp);

    void SetProbability(cTMap * pob_map);

    void SetLocations(QList<QVector3D> positions);
    void SetIncrementDistance(double increment);
    void SetRandomParameters(double rand_increment, double rand_scale);
    void SetRandomRotate(bool rotate);

    void SetDraw(bool draw);

    void OnCreate(GL3DWidget *widget);
    void OnCreatSurfaceBasedObjects(GL3DWidget * widget,GL3DWorld * world, GL3DSurface * surface, QVector3D current_pos);

    int GetRandomDataIndexAt(int r, int c);
    void CreateRandomData();
    void CreateModelDistribution(GL3DWidget * widget,GL3DSurface * s,cTMap * veg_cover, cTMap * veg_h);

    void UpdateBuffer(GL3DWidget * widget,GL3DWorld * world, QVector3D pos);

    void OnUpdate(GL3DWidget * widget,GL3DWorld * world, QVector3D Position, double dt);

    void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

};


#endif
