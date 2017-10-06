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

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "global.h"
#include "ui_full/3D/GL3DWidget.h"

void lisemqt::Setup3DPlot()
{
#ifndef COMPILE_WITHOUT_3D

    //disable 3d tab
    tabWidget->setTabEnabled(3,false);

    glwidget = new GL3DWidget();
    layout_3D->addWidget(glwidget);

    creator = new GL3DWorldCreator();
    creator->LinkToWidget(glwidget);

    QPalette pal = Button3D->palette();
    pal.setColor(QPalette::Button, QColor(Qt::gray));
    Button3D->setAutoFillBackground(true);
    Button3D->setPalette(pal);
    Button3D->update();

    glwidget->setFocus();

#endif

}

void lisemqt::doSet3DFocus()
{
#ifndef COMPILE_WITHOUT_3D
    glwidget->setFocus();
#endif
}

void lisemqt::doUpdateGLSettings()
{
#ifndef COMPILE_WITHOUT_3D

    glwidget->setFocus();

    //get the 3d setting from interface and give them to the lisem world creator
    Settings3D s;
    s.Light_Ambient.setX(this->GL_Light_Ambient_R->value());
    s.Light_Ambient.setY(this->GL_Light_Ambient_G->value());
    s.Light_Ambient.setZ(this->GL_Light_Ambient_B->value());
    s.Light_Ambient.setW(this->GL_Light_Ambient_A->value());
    s.Light_Directional.setX(this->GL_Light_Directional_R->value());
    s.Light_Directional.setY(this->GL_Light_Directional_G->value());
    s.Light_Directional.setZ(this->GL_Light_Directional_B->value());
    s.Light_Directional.setW(this->GL_Light_Directional_A->value());
    s.Light_Directional_Direction.setX(this->GL_Light_Directional_X->value());
    s.Light_Directional_Direction.setY(this->GL_Light_Directional_Y->value());
    s.Light_Directional_Direction.setZ(this->GL_Light_Directional_Z->value());

    /*s.Surface_Draw = this->GL_Surface_Draw->isChecked();
    s.Surface_Micro_Elevation_Scale = this->GL_Surface_Micro_Elevation_Scale->value();
    s.Surface_Mipmap_Distance_1 = this->GL_Surface_Mipmap_Distance_1->value();
    s.Surface_Mipmap_Distance_2 = this->GL_Surface_Mipmap_Distance_2->value();
    s.Surface_Vegetated_Small_Color.setX(this->GL_Surface_Vegetated_Small_R->value());
    s.Surface_Vegetated_Small_Color.setY(this->GL_Surface_Vegetated_Small_G->value());
    s.Surface_Vegetated_Small_Color.setZ(this->GL_Surface_Vegetated_Small_B->value());
    s.Surface_Vegetated_Large_Color.setX(this->GL_Surface_Vegetated_Large_R->value());
    s.Surface_Vegetated_Large_Color.setY(this->GL_Surface_Vegetated_Large_G->value());
    s.Surface_Vegetated_Large_Color.setZ(this->GL_Surface_Vegetated_Large_B->value());
    s.Surface_Bare_Color.setX(this->GL_Surface_Bare_R->value());
    s.Surface_Bare_Color.setY(this->GL_Surface_Bare_G->value());
    s.Surface_Bare_Color.setZ(this->GL_Surface_Bare_B->value());
    s.Surface_Roads_Color.setX(this->GL_Surface_Road_R->value());
    s.Surface_Roads_Color.setY(this->GL_Surface_Road_G->value());
    s.Surface_Roads_Color.setZ(this->GL_Surface_Road_B->value());
    s.Surface_Buildings_Color.setX(this->GL_Surface_Building_R->value());
    s.Surface_Buildings_Color.setY(this->GL_Surface_Building_G->value());
    s.Surface_Buildings_Color.setZ(this->GL_Surface_Building_B->value());
    s.Surface_Erosion_Color.setX(this->GL_Surface_Erosion_Color_R->value());
    s.Surface_Erosion_Color.setY(this->GL_Surface_Erosion_Color_G->value());
    s.Surface_Erosion_Color.setZ(this->GL_Surface_Erosion_Color_B->value());
    s.Surface_Erosion_Color.setW(this->GL_Surface_Erosion_Color_R->value());
    s.Surface_Deposition_Color.setX(this->GL_Surface_Deposition_Color_G->value());
    s.Surface_Deposition_Color.setY(this->GL_Surface_Deposition_Color_B->value());
    s.Surface_Deposition_Color.setZ(this->GL_Surface_Deposition_Color_R->value());
    s.Surface_Deposition_Color.setW(this->GL_Surface_Deposition_Color_G->value());

    s.Water_Draw = this->GL_Water_Draw->isChecked();
    s.Water_Reflectivity = this->GL_Water_Reflectivity->value();
    s.Water_Refractivity = this->GL_Water_Refractivity->value();
    s.Water_Velocity_Scale = this->GL_Water_Velocity_Scale->value();
    s.Water_Micro_Elevation_Scale = this->GL_Water_Micro_Elevation_Scale->value();
    s.Water_Transparancy = this->GL_Water_Transparancy->value();
    s.Water_Deep_Color.setX(this->GL_Water_Deep_R->value());
    s.Water_Deep_Color.setY(this->GL_Water_Deep_G->value());
    s.Water_Deep_Color.setZ(this->GL_Water_Deep_B->value());
    s.Water_Deep_Color.setW(this->GL_Water_Deep_A->value());
    s.Water_Shallow_Color.setX(this->GL_Water_Shallow_R->value());
    s.Water_Shallow_Color.setY(this->GL_Water_Shallow_G->value());
    s.Water_Shallow_Color.setZ(this->GL_Water_Shallow_B->value());
    s.Water_Shallow_Color.setW(this->GL_Water_Shallow_A->value());
    s.Water_Sediment_Color.setX(this->GL_Water_Sediment_R->value());
    s.Water_Sediment_Color.setY(this->GL_Water_Sediment_G->value());
    s.Water_Sediment_Color.setZ(this->GL_Water_Sediment_B->value());
    s.Water_Sediment_Color.setW(this->GL_Water_Sediment_A->value());*/

    s.Clouds_Draw = this->GL_Objects_Clouds_Draw->isChecked();
    s.Rain_Draw = this->GL_Objects_Rain_Draw->isChecked();
    s.Roads_Draw = this->GL_Objects_Roads_Draw->isChecked();
    //s.Roads_Distance = this->GL_Objects_Roads_Distance->value();

    s.Buildings_Draw = this->GL_Objects_Buildings_Draw->isChecked();
    s.Buildings_Distance = this->GL_Objects_Buildings_Distance->value();

    s.Trees_Draw = this->GL_Objects_Trees_Draw->isChecked();
    s.Trees_Distance = this->GL_Objects_Trees_Distance->value();
    //s.Trees_Instances = this->GL_Objects_Trees_Instances->value();
    s.Trees_Increment = this->GL_Objects_Trees_Increment->value();

    s.Grass_Draw = this->GL_Objects_Grass_Draw->isChecked();
    s.Grass_Distance = this->GL_Objects_Grass_Distance->value();
    //s.Grass_Instances = this->GL_Objects_Grass_Instances->value();
    s.Grass_Increment = this->GL_Objects_Grass_Increment->value();
    //s.Grass_Vertical_Scale = this->GL_Objects_Grass_Vertical_Scale->value();

    creator->UpdateWorldSettings(s);

#endif

}
