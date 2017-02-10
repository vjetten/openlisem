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

#include <3D/World/GL3DCamera.h>

void GL3DCamera::Create(GL3DWidget * widget)
{

    m_Widget = widget;

    m_Projection.setToIdentity();

    m_CTranslate.setToIdentity();
    m_CRotate.setToIdentity();

    m_Position.setX(-50);
    m_Position.setY(0);
    m_Position.setZ(0);
    m_Rotation.setX(0);
    m_Rotation.setY(0);
    m_Rotation.setY(0);

    m_CameraMatrix.setToIdentity();

    m_ZNear = 0.00001;
    m_ZFar = 99999999999;

    m_FoV = 55;


}

void GL3DCamera::ResizeViewPort(int w, int h)
{
    // setup viewport, projection etc.:
    glViewport(0, 0, (GLint)w, (GLint)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    aspect = qreal(w) / qreal(h ? h : 1);

    m_Projection.setToIdentity();
    m_Projection.perspective(m_FoV, aspect,m_ZNear,m_ZFar);
}

void GL3DCamera::SetCurrentMatrix()
{
    this->m_CTranslate.setToIdentity();
    this->m_CRotate.setToIdentity();
    this->m_CameraMatrix.setToIdentity();

    // setup world translation/rotation
    this->m_CTranslate.translate(this->m_Position);
    this->m_CRotate.rotate(this->m_Rotation);
    m_CameraMatrix = ( m_Projection*m_CRotate*m_CTranslate);

    QVector3D x;
    QVector3D y;
    QVector3D z;
    this->m_Rotation.getAxes(&x,&y,&z);
    x.normalize();
    y.normalize();
    z.normalize();

    m_CameraMatrix.setToIdentity();
    m_CameraMatrix.lookAt(m_Position,m_Position +QVector3D(1,0,0),QVector3D(0,0,1));
    m_CameraMatrix.perspective(1.0, aspect,m_ZNear,m_ZFar);
#define degreesToRadians(x) x*(3.141592f/180.0f)
}

void GL3DCamera::SetRotation(double a, double b, double c)
{
    this->m_Rotation.fromAxisAndAngle(1.0,0.0,0.0,a);
    this->m_Rotation.fromAxisAndAngle(0.0,1.0,0.0,b);
    this->m_Rotation.fromAxisAndAngle(0.0,0.0,1.0,c);

}

void GL3DCamera::SetPosition(double a, double b, double c)
{
    this->m_Position.setX(a);
    this->m_Position.setY(b);
    this->m_Position.setZ(c);

}

void GL3DCamera::Rotate(double a, double b, double c)
{
    this->rx += a;
    this->ry += b;
    this->rz += c;
    this->m_Rotation *= QQuaternion::fromAxisAndAngle(1.0,0.0,0.0,a);
    this->m_Rotation *= QQuaternion::fromAxisAndAngle(0.0,1.0,0.0,b);
    this->m_Rotation *= QQuaternion::fromAxisAndAngle(0.0,0.0,1.0,c);
}

void GL3DCamera::Move(double a, double b, double c)
{
    this->m_Position.setX(m_Position.x() +a);
    this->m_Position.setY(m_Position.y() +b);
    this->m_Position.setZ(m_Position.z() +c);

}

void GL3DCamera::RotateRelative(double a, double b, double c)
{
    QVector3D x;
    QVector3D y;
    QVector3D z;
    this->m_Rotation.getAxes(&x,&y,&z);
    x.normalize();
    y.normalize();
    z.normalize();

    this->m_Rotation *= QQuaternion::fromAxisAndAngle(x,a);
    this->m_Rotation *= QQuaternion::fromAxisAndAngle(y,b);
    this->m_Rotation *= QQuaternion::fromAxisAndAngle(z,c);
}

void GL3DCamera::MoveRelative(double a, double b, double c)
{
    QVector3D x;
    QVector3D y;
    QVector3D z;
    this->m_Rotation.getAxes(&x,&y,&z);
    x.normalize();
    y.normalize();
    z.normalize();

    this->m_Position = this->m_Position + (x * a);
    this->m_Position = this->m_Position + (y * b);
    this->m_Position = this->m_Position + (z * c);

}
