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

#ifndef Buildings3D
#define Buildings3D

#include <3D/GL3DWidget.h>
#include <3D/Graphics/GL3DModels.h>
#include <3D/World/GL3DSurface.h>

class GL3DBuildings : public GL3DObject
{

public:
    GL3DBuildings() : GL3DObject()
    {
    }

    QList<QVector3D> Building_h_larger_Positions;
    QList<QVector3D> Building_h_larger_Scale;
    QList<float> Building_h_larger_Rotation;

    QList<QVector3D> Building_h_large_Positions;
    QList<QVector3D> Building_h_large_Scale;
    QList<float> Building_h_large_Rotation;

    QList<QVector3D> Building_h_med_Positions;
    QList<QVector3D> Building_h_med_Scale;
    QList<float> Building_h_med_Rotation;

    QList<QVector3D> Building_h_small_Positions;
    QList<QVector3D> Building_h_small_Scale;
    QList<float> Building_h_small_Rotation;


    GL3DModel * m_Building_h_larger;
    GL3DModel * m_Building_h_large;
    GL3DModel * m_Building_h_med;
    GL3DModel * m_Building_h_small;

    GL3DModel * m_Building_l_larger;
    GL3DModel * m_Building_l_large;
    GL3DModel * m_Building_l_med;
    GL3DModel * m_Building_l_small;


    cTMap * m_BuildingCover;

    bool draw = true;

    inline void SetDraw(bool in_draw)
    {
        draw = in_draw;
    }

    void OnCreate(GL3DWidget *widget);
    void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

    void SetBuildingsDistribution(GL3DSurface * s,cTMap * building_cover,cTMap * temp);
};


#endif
