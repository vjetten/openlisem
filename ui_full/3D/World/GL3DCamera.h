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

#ifndef Camera3D_H
#define Camera3D_H

class GL3DWidget;

#include "ui_full/3D/GL3DWidget.h"
#include <QMatrix4x4>
#include <QVector2D>
#include <QVector3D>

#define PI 3.1415926535897932384626433832795
#define PIdiv180 (PI/180.0)

class GL3DCamera
{

public:

    double m_ZNear;
    double m_ZFar;
    double m_FoV;
    GL3DWidget * m_Widget;

    QMatrix4x4 m_Projection;

    QMatrix4x4 m_CLookAtNoTranslation;

    QVector3D m_Position;
    QVector3D m_ViewDir;
    QVector3D m_Up;
    QVector3D m_Right;

    QMatrix4x4 m_CameraMatrix;

    qreal aspect;

    QVector2D m_viewportSize;

    double RotatedX;
    double RotatedY;
    double RotatedZ;

    GL3DCamera()
    {
    };

    void Create(GL3DWidget * widget);

    void ResizeViewPort(int w, int h);

    void SetCurrentMatrix();

    void SetRotation(double a, double b, double c);
    void SetPosition(double a, double b, double c);
    void Rotate(double a, double b, double c);
    void Move(double a, double b, double c);
    void RotateX(double Angle);
    void RotateY(double Angle);
    void RotateZ(double Angle);
    void MoveForward(double Distance );
    void StrafeRight ( double Distance );
    void MoveUpward( double Distance );

    void LookAt(float x, float y, float z);
    void Zoom(double dAngle);
};

#endif
