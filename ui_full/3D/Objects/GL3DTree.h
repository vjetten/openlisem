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

#ifndef GL3DTREE
#define GL3DTREE

#include <3D/GL3DWidget.h>
#include <3D/Graphics/GL3DModels.h>
#include <3D/World/GL3DSurface.h>

class GL3DTrees : public GL3DObject
{

public:
    GL3DTrees() : GL3DObject()
    {
    }

    QList<QVector3D> Tree_Positions;
    QList<float> Tree_Rotation;
    QList<float> Tree_Height;

    GL3DModel * m_Tree_highp;
    GL3DModel * m_Tree_medp;
    GL3DModel * m_Tree_lowp;
    GL3DModel * m_Tree_highp_i;
    GL3DModel * m_Tree_medp_i;
    GL3DModel * m_Tree_lowp_i;

    void OnCreate(GL3DWidget *widget);
    void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

    void SetTreeDistribution(GL3DWidget * widget,GL3DSurface * s,cTMap * veg_cover, cTMap * veg_h);
};


#endif
