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

#ifndef OBJECT3D
#define OBJECT3D

#include "ui_full/3D/GL3DWidget.h"
#include "ui_full/3D/World/GL3DWorld.h"
#include "ui_full/3D/World/GL3DCamera.h"
#include "3D/World/GL3DSurface.h"

class GL3DSurface;

class GL3DObject
{

public:

    GL3DObject()
    {
        added = false;
        created = false;
    };

    int draw_order = 0;

    //bounding box size (for determining if needs to be rendered)
    //dual-directional length, so x width from -0.5sx and 0.5sx
    double sx;
    double sy;
    double sz;

    //position of center of object
    double x;
    double y;
    double z;

    //rotation around the 3 axes
    double rx;
    double ry;
    double rz;

    //is created (usually called when added to world)?
    bool created;
    //is added to world?
    bool added;

    //if current_active is true, OnUpdate is called()
    bool current_active;
    //if current_render is true, OnRender/OnRenderPost is called()
    bool current_render;

    //based on these values, the object is registered in the GL3DWorld
    bool recieve_render_before;
    bool recieve_render;
    bool recieve_render_late;
    bool recieve_render_post;
    bool recieve_keys;
    bool recieve_button;
    bool recieve_mousemove;

    //Overwrite these functions in the specific classes that inherit the object classs
    virtual inline void OnCreate(GL3DWidget * widget)  {return;};
    virtual inline void OnAddToWorld(GL3DWorld * world)  {return;};
    virtual inline void OnCreatSurfaceBasedObjects(GL3DWidget * widget,GL3DWorld * world, GL3DSurface * surface, QVector3D current_pos) {return;};
    virtual inline void OnUpdate(GL3DWidget * widget,GL3DWorld * world, QVector3D position, double dt)  {return;};
    virtual inline void OnRenderBefore(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)  {return;};
    virtual inline void OnRender(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)  {return;};
    virtual inline void OnRenderLate(GL3DWidget * widget,GL3DWorld * world, GL3DCamera* camera, double dt)  {return;};
    virtual inline void OnRenderPost(GL3DWidget * widget,GLuint Texture_source, GLuint Texture_target,GL3DWorld * world, GL3DCamera* camera, double dt)  {return;};
    virtual inline void OnRemoveFromWorld(GL3DWorld * world)  {return;};
    virtual inline void OnDestroy(GL3DWidget * widget)  {return;};
    virtual inline void OnKey(int key)  {return;};
    virtual inline void OnButton(int button)  {return;};
    virtual inline void OnMouseMove(int dx, int dy, double dt)  {return;};

};

#endif
