/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
\file swatinit.cpp
\brief SWATRE: initialize soil profile with inithead maps data and clean up after run

functions:
- SOIL_MODEL * TWorld::InitSwatre(cTMap *profileMap); \n
- void TWorld::CloseSwatre(SOIL_MODEL *s); \n
*/

#include "lerror.h"
#include "model.h"

//--------------------------------------------------------------------------------
// make the 3d structure PIXEL_INFO, based on profile numbers in map
// zone exists is done before
SOIL_MODEL *TWorld::InitSwatre(cTMap *profileMap)
{
    SOIL_MODEL *s = (SOIL_MODEL *)malloc(sizeof(SOIL_MODEL));
    /* TODO check if this needs freeing when error */

    // long nrCel = _nrCols*_nrRows;

    s->minDt = swatreDT;
    s->pixel = new PIXEL_INFO[(long)nrCells];

    // set initial values
    for (long i = 0; i < (long)nrCells; i++) {
        s->pixel[i].profile = nullptr;
        s->pixel[i].dumpHid = 0;  //set to 1 for output of a pixel
        s->pixel[i].tiledrain = 0;
        s->pixel[i].wh = 0;
        s->pixel[i].percolation = 0;
        s->pixel[i].tilenode = -1;      // set tiledrain to 0, and tiledepth to -1 (above surface)
    }

    // give each pixel a profile
    FOR_ROW_COL_MV_L {
        int profnr = swatreProfileNr.indexOf((int)profileMap->Drc);
        if (profnr > 0)
            s->pixel[i_].profile = profileList[profnr];  // pointer to profile
        else
            Error(QString("SWATRE: Profile number %1 in profile.map does not exist in the defenitions in profile.inp").arg((int)profileMap->Drc));

        if(SwitchDumpH || SwitchDumpTheta || SwitchDumpK) {
            s->pixel[i_].dumpHid = SwatreOutput->Drc;
        }
    }}

// fill the inithead structure of each pixel and set tiledrain depth if any
for (int k = 0; k < zone->nrNodes; k++) {
    QString fname = QString("%1.%2").arg(initheadName).arg(k+1, 3, 10, QLatin1Char('0'));
    // make inithead.001 to .00n name

    inith = ReadMap(LDD,fname);
    // get inithead information

    FOR_ROW_COL_MV_L {
        s->pixel[i_].h.append(inith->Drc);//*psiCalibration;
        // find depth of tilenode
        if (SwitchIncludeTile) {
            if (!pcr::isMV(TileDepth->Drc) && TileDepth->Drc > 0) {
                // NOTE depth is in m while node info is in cm, so *100
                // endComp is the depth at the bottom of the compartment, so the tile is <= endcomp
                if (s->pixel[i_].profile->zone->endComp[k] > TileDepth->Drc*100)
                    s->pixel[i_].tilenode = k-1;
            }
        }

        if (SHOWDEBUG) {
            qDebug() << fname <<i_ << s->pixel[i_].h.size() << s->pixel[i_].h[k] << inith->Drc;
        }

    }}
}
qDebug() << "DONE InitSwatre";
return(s);
}
//--------------------------------------------------------------------------------
/// soil model instance to be freed
void TWorld::CloseSwatre(SOIL_MODEL *s)
{
    if (s == nullptr)
        return;

    swatreProfileDef.clear();
    swatreProfileNr.clear();

    FOR_ROW_COL_MV_L {
        s->pixel[i_].h.clear();
    }}
    delete[] s->pixel;

    free(s);
    s = nullptr;
    qDebug() << "closed swatre";
}
//--------------------------------------------------------------------------------
