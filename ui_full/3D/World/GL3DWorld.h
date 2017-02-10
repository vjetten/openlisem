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
#ifndef World3D_H
#define World3D_H

#include "ui_full/3D/GL3DWidget.h"
#include "ui_full/3D/Objects/GL3DObject.h"

class GL3DObject;
class GL3DWidget;

class GL3DWorld
{

    GL3DWidget * m_Widget;

public:

    GL3DWorld()
    {
        is_created = false;
    };

    QList<GL3DObject*> m_ObjectList;
    QList<GL3DObject*> m_KeyObjectList;
    QList<GL3DObject*> m_ButtonObjectList;
    QList<GL3DObject*> m_MouseMoveObjectList;
    QList<GL3DObject*> m_RenderObjectList;
    QList<GL3DObject*> m_RenderPostObjectList;

    double time_multiplier;
    double time;
    double timestep;

    bool is_created;

    void Create(GL3DWidget * widget);
    void AddObject(GL3DObject * object);
    void RemoveObject(GL3DObject * object);
    void RemoveAllObjects();
    void RemoveAndDestroyAllObjects();
    void Destroy();

    void OnRender(GL3DWidget * widget, GL3DCamera* camera, double dt);
    void OnRenderPost(GL3DWidget * widget, GL3DCamera* camera, double dt);
    void OnUpdate(GL3DWidget * widget, double dt);
    void OnKey(int key);
    void OnButton(int button);
    void OnMouseMove(int dx, int dy, double dt);

};

#endif
