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

#include <3D/World/GL3DWorld.h>

void GL3DWorld::Create(GL3DWidget * widget)
{
    m_Widget = widget;
    is_created = true;

}

void GL3DWorld::AddObject(GL3DObject * object)
{

    if(object->added)
    {

        return;
    }
    if(!object->created)
    {
        object->OnCreate(this->m_Widget);
    }
    object->OnAddToWorld(this);

    if(!this->m_ObjectList.contains(object))
    {
        m_ObjectList.append(object);

        if(object->recieve_button)
        {
            m_ButtonObjectList.append(object);
        }
        if(object->recieve_keys)
        {
            m_KeyObjectList.append(object);
        }
        if(object->recieve_mousemove)
        {
            m_MouseMoveObjectList.append(object);
        }
        if(object->recieve_render)
        {
            m_RenderObjectList.append(object);
        }
        if(object->recieve_render_late)
        {
            m_RenderLateObjectList.append(object);
        }
        if(object->recieve_render_post)
        {
            m_RenderPostObjectList.append(object);
        }

    }

}

void GL3DWorld::RemoveObject(GL3DObject * object)
{
    if(this->m_ObjectList.contains(object))
    {
        m_ObjectList.removeAll(object);
        m_KeyObjectList.removeAll(object);
        m_ButtonObjectList.removeAll(object);
        m_MouseMoveObjectList.removeAll(object);
        m_RenderObjectList.removeAll(object);
        m_RenderLateObjectList.removeAll(object);
        m_RenderPostObjectList.removeAll(object);
    }

    object->OnRemoveFromWorld(this);
}
void GL3DWorld::RemoveAllObjects()
{
    for(int i = 0; i < m_ObjectList.length(); i++)
    {
        m_ObjectList.at(i)->OnRemoveFromWorld(this);
    }

    m_ObjectList.clear();
    m_KeyObjectList.clear();
    m_ButtonObjectList.clear();
    m_MouseMoveObjectList.clear();
    m_RenderObjectList.clear();
    m_RenderPostObjectList.clear();

}
void GL3DWorld::RemoveAndDestroyAllObjects()
{
    for(int i = 0; i < m_ObjectList.length(); i++)
    {
        m_ObjectList.at(i)->OnRemoveFromWorld(this);
        m_ObjectList.at(i)->OnDestroy(this->m_Widget);
    }

    m_ObjectList.clear();
    m_KeyObjectList.clear();
    m_ButtonObjectList.clear();
    m_MouseMoveObjectList.clear();
    m_RenderObjectList.clear();
    m_RenderLateObjectList.clear();
    m_RenderPostObjectList.clear();
}

void GL3DWorld::Destroy()
{
    RemoveAndDestroyAllObjects();

    is_created = false;


}


void GL3DWorld::OnRender(GL3DWidget * widget, GL3DCamera* camera, double dt)
{

    for(int i = 0; i < m_RenderObjectList.length();i++)
    {

        m_RenderObjectList.at(i)->OnRender(widget,this,camera,dt);
    }
}

void GL3DWorld::OnRenderLate(GL3DWidget * widget, GL3DCamera* camera, double dt)
{

    for(int i = 0; i < m_RenderLateObjectList.length();i++)
    {

        m_RenderLateObjectList.at(i)->OnRenderLate(widget,this,camera,dt);
    }
}

void GL3DWorld::OnRenderPost(GL3DWidget * widget, GL3DCamera* camera, double dt)
{
    for(int i = 0; i < m_RenderPostObjectList.length();i++)
    {
        m_RenderPostObjectList.at(i)->OnRenderPost(widget,this,camera,dt);
    }
}


void GL3DWorld::OnUpdate(GL3DWidget * widget, double dt)
{
    for(int i = 0; i < m_ObjectList.length();i++)
    {
        m_ObjectList.at(i)->OnUpdate(widget,this,dt);
    }
}

void GL3DWorld::OnKey(int key)
{
    for(int i = 0; i < m_KeyObjectList.length();i++)
    {
        m_KeyObjectList.at(i)->OnKey(key);
    }
}

void GL3DWorld::OnButton(int button)
{
    for(int i = 0; i < m_ButtonObjectList.length();i++)
    {
        m_ButtonObjectList.at(i)->OnButton(button);
    }
}

void GL3DWorld::OnMouseMove(int dx, int dy, double dt)
{
    for(int i = 0; i < m_MouseMoveObjectList.length();i++)
    {
        m_MouseMoveObjectList.at(i)->OnMouseMove(dx,dy,dt);
    }
}

