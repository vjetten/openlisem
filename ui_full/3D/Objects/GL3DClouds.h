
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

#ifndef GL3DCLOUDS_H
#define GL3DCLOUDS_H

#include <3D/Objects/GL3DObject.h>
#include <3D/World/GL3DSurface.h>
#include <3D/Graphics/GL3DMath.h>


class GL3DClouds : public GL3DObject
{

public:
    GL3DClouds() : GL3DObject()
    {
            recieve_render = true;
            m_Noise = new PerlinNoise();
    }

    PerlinNoise *m_Noise;

    float * randomdata;

    GL3DModel * m_Model_highp;
    GL3DModel * m_Model_lowp;

    GL3DTexture * m_Texture_Perlin1;
    GL3DTexture * m_Texture_Perlin2;

    GL3DSurface * m_Surface;

    QList<QVector3D> m_ParticleLocations;
    QList<QVector3D> m_ParticleScale;
    QList<QVector3D> m_ParticleRotations;
    QList<float> m_ParticleAlpha;

    void SetDraw(bool draw);

    void CreateCloudData();

    void OnCreate(GL3DWidget *widget);
    void OnCreatSurfaceBasedObjects(GL3DWidget * widget,GL3DWorld * world, GL3DSurface * surface, QVector3D current_pos);

    void OnUpdate(GL3DWidget * widget,GL3DWorld * world, QVector3D Position, double dt);

    void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt);
    void OnDestroy(GL3DWidget *widget);

};



#endif // GL3DCLOUDS_H
