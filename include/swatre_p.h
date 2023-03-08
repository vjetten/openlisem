/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file swatre_p.h
  \brief SWATRE private declarations and structures
*/

#ifndef SWATRE_P_H
#define SWATRE_P_H

#include "swatresoillut.h"


#define MAX_NODES            20
#define MAX_NODES_P          (MAX_NODES+3)

/* maximum amount of ponding that is regarded as no ponding (0) */
#define POND_EPS             (1.0E-6)
/* maximum amount of time for which it is not worth doing an iteration */
#define TIME_EPS             (1.0E-6)

#define NrNodes(profile)        (zone->nrNodes)  //profile->zone->nrNodes)
#define Dz(profile)             (profile->zone->dz)
#define DistNode(profile)       (profile->zone->disnod)
#define Horizon(profile, node)  (profile->horizon[node])
/* sizeof intermediate arrays is fixed to optimize computation: */

//-------------------------------------------------------------
/// SWATRE structure geometry of profile: node distances etc.
/**	profile node setup:
-		 endComp is what is in the profile.inp file, the bottom of the layer
-		 dz = (endComp[i-1] - endComp[i]) is negative layer thickness
-		 z = 0.5*(dz[i-1]+dz[i]) is negative centre of compartment, nodes
-		 disnod = z[i]-z[i-1] is negative distance between centres, nodes
- schematics:
\code
		  -------   surface    -       - z[0]-
			  o                  |dz[0] -      | disnod[0]
		  -------   endComp[0] -        |z[1]-
			  o                  |dz[1] -      | disnod[1]
		  -------   endcomp[1] -        |z[2]-
			  o                  |dz[2] -      | disnod[2]
		  -------   endcomp[2] -
		 etc.
\endcode
*/
typedef struct ZONE   {
	/* structure explaining how the soil subdivided */
		int  nrNodes;     	/*!< nr. of nodes (equals nr. of compartments), arrays with nrNodes elements: */
		double *dz;       	/*!< compartment size (cm.) used as negative in [SWATRE]*/
		double *z;        	/*!< position of nodal point relative to top soil (Cm) e.g -21 means 21 cm below top of soil */
		double *endComp;   /*!< end of compartment i only used when reading the profiles, arrays with nrNodes+1 elements: */
		double *disnod;    /*!< distance between nodal points, 0 is between top-profile and first nodal point last is between bottom and last nodal point */

} ZONE;
//---------------------------------------------------------------------------
/* change this structure if we add VanGenughten eqs. */
/// SWATRE structure with names and pointers to land use tables
typedef struct HORIZON {
   char *name;  	/** name of horizon: filename of lut-table */
   LUT  *lut;     /** lut of theta, h, k, dmch, dmcc */
} HORIZON;
//---------------------------------------------------------------------------
/// SWATRE structure with horizon and node info
typedef struct PROFILE {
   int            profileId; 	/** number identifying this profile  >= 0 */
   const ZONE     *zone; 		/** array with zone.nrNodes elements: containing z, dz, node distance etc*/
   const HORIZON  **horizon; 	/** ptr to horizon information this node belongs to */
   int nrNodes;
} PROFILE;
//---------------------------------------------------------------------------
typedef double NODE_ARRAY[MAX_NODES_P];
//---------------------------------------------------------------------------
/// SWATRE structure for actual matrix head and profile information
/** SWATRE structure for actual matrix head and profile information
   Each pixel in the map profile.map has this info. PROFILE gives the profiel layout
   The h array contains the matrix potentials ported to the next timestep.
   Specific can be put here such as tile drain flux, drip irrigation flux etc
  */
typedef struct PIXEL_INFO {
   const PROFILE *profile;    /** profile this pixel belongs to */
   double        *h;          /** array of MAX_NODES nodes with matrix head */
   double        *theta;          /** array of MAX_NODES nodes with matrix head */
   double        *k;          /** array of MAX_NODES nodes with matrix head */
   double        currDt;      /** current size of SWATRE timestep */
   double        tiledrain;   /** drainage into tiledrin system at a given depth */
   int           tilenode;    /** nearest node that has the tiledrain */
   int           repellency;  /** water repellency will be calculated if 1 */
   int           dumpHid;     /** if 0 then no head output else write to file amed Hx where x is dumpH value */
   int nrNodes;
   double var; // variable can be used for anything
} PIXEL_INFO;
//---------------------------------------------------------------------------


#endif // SWATRE_P_H
