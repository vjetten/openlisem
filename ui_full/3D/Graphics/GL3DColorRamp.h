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

#ifndef GL3DCOLORRAMP_H
#define GL3DCOLORRAMP_H

#include "3D/Graphics/GL3DMath.h"
#include <Qt>
#include <QVector4D>
#include <QList>

#define GL3D_COLOR_RAMP_LINEAR 0
#define GL3D_COLOR_RAMP_LOGARITMIC 0
#define GL3D_COLOR_RAMP_SHIFTED_LOGARITMIC 0

class GL3DColorRamp
{

public:

    QList<QVector4D> m_ColorList;
    QList<double> m_ColorStopList;
    int mode = GL3D_COLOR_RAMP_LINEAR;

    GL3DColorRamp()
    {

    };

    GL3DColorRamp(QVector4D color1)
    {

    };

    GL3DColorRamp(QVector4D color1, QVector4D color2)
    {

    };

    GL3DColorRamp(QVector4D color1, double stop2,QVector4D color2,QVector4D color3)
    {

    };

    GL3DColorRamp(QVector4D color1, double stop2,QVector4D color2,double stop3,QVector4D color3,QVector4D color4)
    {

    };

    void Clear();
    void AddColorStop(double stop, QVector4D color);
    void SetMode();
    QVector4D GetColor(double min, double max, double val);
    void Sort();
    bool is_sorted = true;
};

#endif
