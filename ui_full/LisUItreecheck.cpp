/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/


/*!
  \file LisUItreecheck.cpp
  \brief map tree interaction

 * function that determine reactions of the map tree structure when
 * the user checks main options in the interface: disable or enbale braches
 * e.g. channel maps become enabled when this option is chosen
 */


#include "lisemqt.h"
#include "global.h"
/*
//#define RAINFALLMAPS 0
//#define CATCHMENTMAPS 1
//#define LANDUSEMAPS 2
//#define SURFACEMAPS 3
//#define INFILTRATIONMAPS 4
//#define CHANNELMAPS 5
//#define HOUSESMAPS 6
//#define EROSIONMAPS 7
//#define CONSERVATIONMAPS 7
//#define TILEDRAINMAPS 8
//#define PESTMAPS 11

*/
//!!!  all previous code was obsolete and interfered with changes !!!


//--------------------------------------------------------------------
void lisemqt::on_checkExpandActive_clicked()
{
    if (!checkExpandActive->isChecked())
        treeView->collapseAll();
    else
        for (int i = 0; i < MapNameModel->rowCount(); i++)
        {
            if (MapNameModel->getFlag(i))
                treeView->expand(MapNameModel->index(i,0));
        }
}
//--------------------------------------------------------------------
void lisemqt::RunAllChecks()
{
    for (int i = 0; i < 12; i++)
        checkMapNameModel(i, 0, false);

    // infiltration has a second level
    checkMapNameModel(INFILTRATIONMAPS, 10, false);
    checkMapNameModel(INFILTRATIONMAPS, 11, false);
    checkMapNameModel(INFILTRATIONMAPS, 12, false);
    checkMapNameModel(INFILTRATIONMAPS, 13, false);
    checkMapNameModel(INFILTRATIONMAPS, 14, false);
    checkMapNameModel(INFILTRATIONMAPS, 15, false);

    checkExpandActive->setChecked(false);
    treeView->collapseAll();
}
