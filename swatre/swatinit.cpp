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
\file swatinit.cpp
\brief SWATRE: initialize soil profile with inithead maps data and clean up after run

functions:
- SOIL_MODEL *TWorld::InitSwatre(cTMap *profileMap); \n
- void TWorld::FreeSwatreInfo(void);\n
- void TWorld::CloseSwatre(SOIL_MODEL *s); \n
*/

#include "lerror.h"
#include "model.h"

//--------------------------------------------------------------------------------
// make the 3D structure PIXEL_INFO, based on profile numbers in map
// needs zone info which needs to be done before in readswatreinput
// read optional Hinit maps
SOIL_MODEL *TWorld::InitSwatre(cTMap *profileMap)
{
   //  int numThreads = omp_get_max_threads();
   // // QVector<NODES> threadBuffers;
   //  threadBuffers.reserve(numThreads);

   //  for (int i = 0; i < numThreads; ++i) {
   //      threadBuffers.append(NODES(MAX_NODES+3));
   //  }

    SOIL_MODEL *s = new SOIL_MODEL;

    s->pixel = new PIXEL_INFO[nrValidCells];

    // set initial values
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        s->pixel[i_].r = r;
        s->pixel[i_].c = c;
        s->pixel[i_].profile = nullptr;
        s->pixel[i_].tiledrain = 0;
        s->pixel[i_].wh = 0;
        s->pixel[i_].percolation = 0;
        s->pixel[i_].tilenode = -1;      // set tiledrain to 0, and tiledepth to -1 (above surface)
        s->pixel[i_].currDt = swatreDT;
    }}

    // give each pixel a profile
   // #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        int profilenr = static_cast <int>(profileMap->Drc);
        int profindex = swatreProfileNr.indexOf(profilenr);

        if (profilenr > 0)
            s->pixel[i_].profile = profileList[profindex];  // pointer to profile

    }}


    // fill the inithead structure of each pixel and set tiledrain depth if any
    for (int k = 0; k < zone->nrNodes; k++) {

        if (!SwitchHinit4all) {
            QString name = QString("%1.%2").arg(initheadName).arg(k+1, 3, 10, QLatin1Char('0'));
            // make inithead.001 to .00n name
            cTMap* map = ReadMap(LDD, name);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                map->Drc *= psiCalibration;
            }}
            inith->append(map);
        } else {
            cTMap* map = NewMap(HinitValue);
            // confusing to have a user defined value and a calibration on it
            // #pragma omp parallel for num_threads(userCores)
            // FOR_ROW_COL_MV_L {
            //     map->Drc *= psiCalibration;
            // }}
            inith->append(map);
        }

        // get inithead information
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            cTMap *map = inith->at(k);
            s->pixel[i_].h.append(map->Drc);
           // s->pixel[i_].theta.append(0.5);

            // find depth of tilenode
            if (SwitchIncludeTile) {
                if (!pcr::isMV(TileDepth->Drc) && TileDepth->Drc > 0) {
                    // NOTE depth is in m while node info is in cm, so *100
                    // endComp is the depth at the bottom of the compartment, so the tile is <= endcomp
                    if (s->pixel[i_].profile->zone->endComp[k] > TileDepth->Drc*100)
                        s->pixel[i_].tilenode = k-1;
                }
            }
        }}
    }

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

    delete[] s->pixel;
    //free(s);
    delete s;
    s = nullptr;

    //qDebug() << "closed swatre";
}
//--------------------------------------------------------------------------------
// free the zone, luts and profiles, these are only pointers in PIXEL_INFO
void  TWorld::FreeSwatreInfo(void)
{
    if (zone == nullptr)
       return;

    if (zone != nullptr) {
        zone->dz.clear();
        zone->z.clear();
        zone->endComp.clear();
        zone->disnod.clear();
        zone->rootz.clear();
        delete zone;
        zone = nullptr;
    }

    if (profileList != nullptr) {
        if (profileList[0] != nullptr) {
            for(int i=0; i < sizeProfileList; i++)
                if (profileList[i] != nullptr)
                    free(profileList[i]);
        }
        //free(profileList);
        delete profileList;
        profileList = nullptr;
    }

    if (horizonList != nullptr) {
        for(int i=0; i < nrHorizonList; i++)
        {
            for(int k = 0; k < 5; k++)
                horizonList[i]->lut->hydro[k].clear();
            //free(horizonList[i]->lut);
            delete horizonList[i]->lut;
            free(horizonList[i]);
        }
        free(horizonList);
        horizonList = nullptr;
    }

    nrHorizonList = 0;
    sizeHorizonList = 0;

    // free pixel_info
    if (SwatreSoilModel != nullptr)
        CloseSwatre(SwatreSoilModel);
    if (SwatreSoilModelCrust != nullptr)
        CloseSwatre(SwatreSoilModelCrust);
    if (SwatreSoilModelCompact != nullptr)
        CloseSwatre(SwatreSoilModelCompact);
    if (SwatreSoilModelGrass != nullptr)
        CloseSwatre(SwatreSoilModelGrass);

    initSwatreStructure = false;

    DEBUG("SWATRE mem freed");
}
//--------------------------------------------------------------------------------
