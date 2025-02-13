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



#ifndef SWATRE_P_H
#define SWATRE_P_H

#include <QtCore>

#define MAX_NODES            20
//#define MAX_NODES_P          (MAX_NODES+3)

// maximum amount of ponding that is regarded as no ponding (0)
#define POND_EPS             (1.0E-8)
// maximum amount of time for which it is not worth doing an iteration
#define TIME_EPS             (1.0E-6)

#define THETA_COL	    0 // moisture content theta (-)
#define H_COL           1 // suction (-cm)
#define K_COL           2 // hydraulic conductivity (cm/day)
#define DMCH_COL        3 // h for differential moisture capacity (cm)
#define DMCC_COL        4 // theta for differential moisture capacity (-)

//-------------------------------------------------------------
/// SWATRE structure geometry of profile: node distances etc.
/**	profile node setup:
-		 endComp is what is in the profile.inp file, the bottom of the layer
-		 dz = (endComp[i-1] - endComp[i]) is negative layer thickness
-		 z = 0.5*(dz[i-1]+dz[i]) is negative centre of compartment, nodes
-		 disnod = z[i]-z[i-1] is negative distance between centres, nodes
- schematics:

          -------   surface    -       - z[0]-
              o                  |dz[0] -      | disnod[0]
          -------   endComp[0] -        |z[1]-
              o                  |dz[1] -      | disnod[1]
          -------   endcomp[1] -        |z[2]-
              o                  |dz[2] -      | disnod[2]
          -------   endcomp[2] -
         etc.

*/
typedef struct ZONE   {
    int  nrNodes;
    QVector <double> dz;
    QVector <double> z;
    QVector <double> endComp;
    QVector <double> disnod;   
} ZONE;
//---------------------------------------------------------------------------
/// SWATRE Land use tables, Rows is nr of lines in table, always 3+2 cols
typedef struct LUT {
    int   Rows;
    QVector<double> hydro[5]; // theta,h,k,dmch,dmcc
} LUT;
//---------------------------------------------------------------------------
typedef struct HORIZON {
    QString name;
    LUT  *lut;     // lut of theta, h, k, dmch, dmcc
} HORIZON;
//---------------------------------------------------------------------------
typedef struct PROFILE {
    int            profileId; 	/** number identifying this profile  >= 0 */
    const ZONE     *zone; 		/** array with zone.nrNodes elements: containing z, dz, node distance etc*/
    const HORIZON  **horizon; 	/** ptr to horizon information this node belongs to */
    QVector <double> KsatCal;
} PROFILE;
//---------------------------------------------------------------------------
typedef double NODE_ARRAY[MAX_NODES+3];
//---------------------------------------------------------------------------
typedef struct PIXEL_INFO {
    const PROFILE *profile;    /** profile this pixel belongs to */
    QVector <double> h;
    double wh;
    double impfrac;
    double percolation;
    double theta; // for pesticides?
    double tiledrain;   /** drainage into tiledrin system at a given depth */
    double corrKsOA;
    double corrKsOB;
    double corrKsDA;
    double corrKsDB;
    double corrPOA;
    double corrPOB;
    double corrPDA;
    double corrPDB;
    int tilenode;    /** nearest node that has the tiledrain */
   // int dumpHid;     /** if 0 then no head output else write to file amed Hx where x is dumpH value */
} PIXEL_INFO;
//---------------------------------------------------------------------------
typedef struct SOIL_MODEL {
    struct PIXEL_INFO  *pixel;
    double minDt;
} SOIL_MODEL;
//---------------------------------------------------------------------------

#endif // SWATRE_P_H
