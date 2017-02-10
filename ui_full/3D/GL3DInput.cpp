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

#include "ui_full/3D/GL3DWidget.h"
#include <QtGui>
#include <QKeyEvent>
#include <QMetaEnum>

void GL3DWidget::keyPressEvent( QKeyEvent *keyEvent )
{
    int key = keyEvent->key();
    {
        if(Qt::Key_Up == key)
        {
            KeyPressed[GL3D_KEY_UP] = true;
        }
        if(Qt::Key_Down == key)
        {
            KeyPressed[GL3D_KEY_UP] = true;
        }
        if(Qt::Key_Left == key)
        {
           KeyPressed[GL3D_KEY_LEFT] = true;
        }
        if(Qt::Key_Right == key)
        {
           KeyPressed[GL3D_KEY_RIGHT] = true;
        }
        if(Qt::Key_W == key)
        {
           KeyPressed[GL3D_KEY_W] = true;
        }
        if(Qt::Key_A == key)
        {
           KeyPressed[GL3D_KEY_A] = true;
        }
        if(Qt::Key_S == key)
        {
           KeyPressed[GL3D_KEY_S] = true;
        }
        if(Qt::Key_D == key)
        {
           KeyPressed[GL3D_KEY_D] = true;
        }
        if(Qt::Key_Shift == key)
        {
           KeyPressed[GL3D_KEY_SHIFT] = true;
        }
        if(Qt::Key_Control == key)
        {
           KeyPressed[GL3D_KEY_CTRL] = true;
        }
        if(Qt::Key_Q == key)
        {
           KeyPressed[GL3D_KEY_Q] = true;
        }
        if(Qt::Key_E == key)
        {
           KeyPressed[GL3D_KEY_E] = true;
        }
    }

    {

        if(Qt::Key_Up == key)
        {
            KeyPressedT[GL3D_KEY_UP] = true;
        }
        if(Qt::Key_Down == key)
        {
            KeyPressedT[GL3D_KEY_UP] = true;
        }
        if(Qt::Key_Left == key)
        {
            KeyPressedT[GL3D_KEY_LEFT] = true;
        }
        if(Qt::Key_Right == key)
        {
            KeyPressedT[GL3D_KEY_RIGHT] = true;
        }
        if(Qt::Key_W == key)
        {
            KeyPressedT[GL3D_KEY_W] = true;
        }
        if(Qt::Key_A == key)
        {
            KeyPressedT[GL3D_KEY_A] = true;
        }
        if(Qt::Key_S == key)
        {
            KeyPressedT[GL3D_KEY_S] = true;
        }
        if(Qt::Key_D == key)
        {
            KeyPressedT[GL3D_KEY_D] = true;
        }
        if(Qt::Key_Shift == key)
        {
            KeyPressedT[GL3D_KEY_SHIFT] = true;
        }
        if(Qt::Key_Control == key)
        {
            KeyPressedT[GL3D_KEY_CTRL] = true;
        }
        if(Qt::Key_Q == key)
        {
           KeyPressedT[GL3D_KEY_Q] = true;
        }
        if(Qt::Key_E == key)
        {
           KeyPressedT[GL3D_KEY_E] = true;
        }
    }

}

void GL3DWidget::keyReleaseEvent( QKeyEvent *keyEvent )
{
    int key = keyEvent->key();
    {
        if(Qt::Key_Up == key)
        {
            KeyPressed[GL3D_KEY_UP] = false;
        }
        if(Qt::Key_Down == key)
        {
            KeyPressed[GL3D_KEY_UP] = false;
        }
        if(Qt::Key_Left == key)
        {
           KeyPressed[GL3D_KEY_LEFT] = false;
        }
        if(Qt::Key_Right == key)
        {
           KeyPressed[GL3D_KEY_RIGHT] = false;
        }
        if(Qt::Key_W == key)
        {
           KeyPressed[GL3D_KEY_W] = false;
        }
        if(Qt::Key_A == key)
        {
           KeyPressed[GL3D_KEY_A] = false;
        }
        if(Qt::Key_S == key)
        {
           KeyPressed[GL3D_KEY_S] = false;
        }
        if(Qt::Key_D == key)
        {
           KeyPressed[GL3D_KEY_D] = false;
        }
        if(Qt::Key_Shift == key)
        {
           KeyPressed[GL3D_KEY_SHIFT] = false;
        }
        if(Qt::Key_Control == key)
        {
           KeyPressed[GL3D_KEY_CTRL] = false;
        }
        if(Qt::Key_Q == key)
        {
           KeyPressed[GL3D_KEY_Q] = false;
        }
        if(Qt::Key_E == key)
        {
           KeyPressed[GL3D_KEY_E] = false;
        }
    }

}

void GL3DWidget::mouseMoveEvent(QMouseEvent *event)
{
    QTransform t;
    t.scale(1, -1);
    t.translate(0, -Height+1);
    QPoint pos = event->pos() * t;

    int tempx = MouseOldX;
    int tempy = MouseOldY;

    MouseOldX = MouseX;
    MouseOldY = MouseY;

    MouseX = pos.x();
    MouseY = pos.y();

    MouseDX = MouseX - tempx;
    MouseDY = MouseY - tempy;

    if(KeepMouseMiddle)
    {
        QCursor c = cursor();
        c.setPos(mapToGlobal(QPoint(Width / 2, Height / 2)));
        c.setShape(Qt::BlankCursor);
        setCursor(c);

        QPoint pos = this->mapFromGlobal(c.pos());
        MouseX = pos.x();
        MouseY = pos.y();
    }else
    {
        QCursor c = cursor();
        c.setShape(Qt::ArrowCursor);
    }

}
void GL3DWidget::mousePressEvent(QMouseEvent *event)
{
    switch (event->button())
    {
        case Qt::LeftButton:
            MousePressed[GL3D_MOUSE_LEFT] = true;
        case Qt::RightButton:
            MousePressed[GL3D_MOUSE_RIGHT] = true;
        case Qt::MiddleButton:
            MousePressed[GL3D_MOUSE_MIDDLE] = true;
    }
    switch (event->button())
    {
        case Qt::LeftButton:
            MousePressedT[GL3D_MOUSE_LEFT] = true;
        case Qt::RightButton:
            MousePressedT[GL3D_MOUSE_RIGHT] = true;
        case Qt::MiddleButton:
            MousePressedT[GL3D_MOUSE_MIDDLE] = true;
    }
}

void GL3DWidget::mouseReleaseEvent(QMouseEvent *event)
{
    switch (event->button())
    {
        case Qt::LeftButton:
            MousePressed[GL3D_MOUSE_LEFT] = false;
        case Qt::RightButton:
            MousePressed[GL3D_MOUSE_RIGHT] = false;
        case Qt::MiddleButton:
            MousePressed[GL3D_MOUSE_MIDDLE] = false;
    }

}

void GL3DWidget::UseInput()
{
    if(KeyPressedT[GL3D_KEY_SHIFT] == true)
    {
        KeepMouseMiddle = !KeepMouseMiddle;
    }

    for(int i = 0; i < GL3D_INPUT_NRKEYS; i++)
    {
        KeyPressedT[i] = false;
    }
    for(int i = 0; i < GL3D_INPUT_NRBUTTONS; i++)
    {
        MousePressedT[i] = false;
    }
    MouseDX = 0;
    MouseDY = 0;
    MouseScroll = 0;

}

void GL3DWidget::wheelEvent(QWheelEvent *event)
{
    MouseScroll += event->angleDelta().y();

}
