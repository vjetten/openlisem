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
    //SOIL_MODEL *s = (SOIL_MODEL *)malloc(sizeof(SOIL_MODEL));
    SOIL_MODEL *s = new SOIL_MODEL;

    s->minDt = swatreDT;
    s->pixel = new PIXEL_INFO[(long)nrCells];

    // set initial values
    for (long i = 0; i < (long)nrCells; i++) {
        s->pixel[i].profile = nullptr;
        //s->pixel[i].dumpHid = 0;  //set to 1 for output of a pixel
        s->pixel[i].tiledrain = 0;
        s->pixel[i].wh = 0;
        s->pixel[i].percolation = 0;
        s->pixel[i].tilenode = -1;      // set tiledrain to 0, and tiledepth to -1 (above surface)        
        s->pixel[i].impfrac = 0;        // fraction roads, houses etc, for first node

        s->pixel[i].corrKsOA = 1.0;
        s->pixel[i].corrKsOB = 0.0;
        s->pixel[i].corrKsDA = 1.0;
        s->pixel[i].corrKsDB = 0.0;
        s->pixel[i].corrPOA = 1.0;
        s->pixel[i].corrPOB = 0.0;
        s->pixel[i].corrPDA = 1.0;
        s->pixel[i].corrPDB = 0.0;
    }

    // give each pixel a profile
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        int profnr = swatreProfileNr.indexOf((int)profileMap->Drc);

        if (profnr > 0)
            s->pixel[i_].profile = profileList[profnr];  // pointer to profile
        // profile = <= 0 now set to impermeable
        s->pixel[i_].impfrac = fractionImperm->Drc;

        if (SwitchOMCorrection) {
            // these correction come from calculations based on Saxton and rawls
            // Ks in mm/h convert to cm/s, affects factor B of the regression
            double OM2 = OMcorr->Drc*OMcorr->Drc;
            s->pixel[i_].corrKsOA = 0.0026*OM2 + 0.0359*OMcorr->Drc + 1;
            s->pixel[i_].corrKsOB = 0.1/3600.0*(0.253*OM2 + 2.9368*OMcorr->Drc + 0.0007);
            s->pixel[i_].corrPOA  = -0.001*OM2 + 0.1014*OMcorr->Drc + 1.0;
            s->pixel[i_].corrPOB  = 0.0006*OM2 - 0.0282*OMcorr->Drc;
            // pore A -0.001x2 + 0.1014x + 1
            // pore B 0.0006x2 - 0.0282x
        }
        if (SwitchDensCorrection) {
            double D2 = DensFact->Drc*DensFact->Drc;
            // the regression is made with ks in cm/s, this affects B, not A: mm/h cm/s = *0.1/3600.0
            s->pixel[i_].corrKsDA = 3.1429*D2 - 9.5657*DensFact->Drc + 7.4229;
            s->pixel[i_].corrKsDB = 0.1/3600.0*(135.4*D2 - 311.07*DensFact->Drc + 175.67);
            s->pixel[i_].corrPDA  = DensFact->Drc;
            s->pixel[i_].corrPDB   = -1.0 * DensFact->Drc + 1.0;
        }
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
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                map->Drc *= psiCalibration;
            }}
            inith->append(map);
        }

        // get inithead information
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            cTMap *map = inith->at(k);
            s->pixel[i_].h.append(map->Drc);

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
        delete zone;
        zone = nullptr;
    }

    if (profileList != nullptr) {
        if (profileList[0] != nullptr) {
            for(int i=0; i < sizeProfileList; i++)
                if (profileList[i] != nullptr)
                    free(profileList[i]);
        }
        free(profileList);
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
