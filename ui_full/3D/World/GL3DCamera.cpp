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
#define degreesToRadians(x) x*(3.141592f/180.0f)

void GL3DCamera::Create(GL3DWidget * widget)
{

    m_Widget = widget;

    m_Projection.setToIdentity();
    m_CLookAtNoTranslation.setToIdentity();
    m_CameraMatrix.setToIdentity();

    m_Position.setX(200);
    m_Position.setY(200);
    m_Position.setZ(200);
    m_ViewDir.setX(1);
    m_ViewDir.setY(0);
    m_ViewDir.setZ(0);
    m_Up.setX(0);
    m_Up.setY(1);
    m_Up.setZ(0);
    m_Right.setX(0);
    m_Right.setY(0);
    m_Right.setZ(1);

    RotatedX = 0;
    RotatedY = 0;
    RotatedZ = 0;

    m_ZNear = 0.1;
    m_ZFar = 999999;

    m_FoV = 55;

}

void GL3DCamera::ResizeViewPort(int w, int h)
{
    // setup viewport, projection etc.:
    glViewport(0, 0, (GLint)w, (GLint)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    aspect = qreal(w) / qreal(h ? h : 1);

    m_viewportSize = QVector2D(float(w),float(h));

    m_Projection.setToIdentity();
    m_Projection.perspective(m_FoV, aspect,m_ZNear,m_ZFar);
}

void GL3DCamera::SetCurrentMatrix()
{
    m_Projection.setToIdentity();
    m_Projection.perspective(m_FoV, aspect,m_ZNear,m_ZFar);

    QMatrix4x4 matrix1;
    QMatrix4x4 matrix2;
    matrix1.setToIdentity();
    matrix2.setToIdentity();

    this->m_CLookAtNoTranslation.setToIdentity();
    this->m_CameraMatrix.setToIdentity();

    m_CameraMatrix = m_Projection;
    m_CLookAtNoTranslation = m_Projection;

    matrix2.lookAt(m_Position,m_Position+m_ViewDir,m_Up);
    matrix1.lookAt(QVector3D(0,0,0),m_ViewDir,m_Up);

    m_CLookAtNoTranslation = m_CLookAtNoTranslation*matrix1;
    m_CameraMatrix = m_CameraMatrix*matrix2;


}

void GL3DCamera::SetRotation(double a, double b, double c)
{

}

void GL3DCamera::SetPosition(double a, double b, double c)
{
    this->m_Position.setX(a);
    this->m_Position.setY(b);
    this->m_Position.setZ(c);

}

void GL3DCamera::RotateX(double Angle)
{
    RotatedX += Angle;

    //Rotate viewdir around the right vector:
    m_ViewDir = QVector3D(m_ViewDir*cos(Angle*PIdiv180)
                                + m_Up*sin(Angle*PIdiv180)).normalized();

    //now compute the new UpVector (by cross product)
    m_Up = QVector3D(0,1,0);//QVector3D::crossProduct(m_ViewDir, m_Right)*-1;

}

void GL3DCamera::RotateY(double Angle)
{
    RotatedY += Angle;

    //Rotate viewdir around the up vector:
    m_ViewDir =  QVector3D(m_ViewDir*cos(Angle*PIdiv180)
                                - m_Right*sin(Angle*PIdiv180)).normalized();

    //now compute the new RightVector (by cross product)
    m_Right = QVector3D::crossProduct(m_ViewDir, m_Up);
}

void GL3DCamera::RotateZ(double Angle)
{
    RotatedZ += Angle;

    //Rotate viewdir around the right vector:
    m_Right =  QVector3D(m_Right*cos(Angle*PIdiv180)
                                + m_Up*sin(Angle*PIdiv180)).normalized();

    //now compute the new UpVector (by cross product)
    m_Up = QVector3D(0,1,0);//QVector3D::crossProduct(m_ViewDir, m_Right)*-1;
}

void GL3DCamera::Move(double a, double b, double c)
{
    this->m_Position.setX(m_Position.x() +a);
    this->m_Position.setY(m_Position.y() +b);
    this->m_Position.setZ(m_Position.z() +c);

}

void GL3DCamera::MoveForward(double Distance )
{
    m_Position += (m_ViewDir*(-Distance));
}
void GL3DCamera::StrafeRight ( double Distance )
{
    m_Position += (m_Right*Distance);
}

void GL3DCamera::MoveUpward( double Distance )
{
    m_Position += (m_Up*Distance);
}

void GL3DCamera::LookAt(float x, float y, float z)
{
    QVector3D look= QVector3D(x,y,z);
    m_ViewDir = (look-m_Position);
    m_ViewDir.normalize();

}

void GL3DCamera::Zoom( double dangle )
{

    this->m_FoV += dangle;
}

