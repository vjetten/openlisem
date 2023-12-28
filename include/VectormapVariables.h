/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
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

#ifndef VECTORMAPVARIABLES_H
#define VECTORMAPVARIABLES_H

#endif // VECTORMAPVARIABLES_H

double *vtma;
double *vtmb;
double *vtmc;
double *vtmd;

// infiltration, percolation, redistribution
double *vKsat1;
double *vThetaS1;
double *vThetaI1;
double *vThetaI1a;
double *vThetaR1;
double *vPsi1;
double *vThetaFC1;
double *vlambda1;
double *vSoilDepth1;
double *vSoilDepth1init;

double *vKsat2;
double *vThetaS2;
double *vThetaI2;
double *vThetaI2a;
double *vlambda2;
double *vThetaFC2;
double *vPsi2;
double *vThetaR2;
double *vSoilDepth2;
double *vSoilDepth2init;

double *vKsat3;
double *vThetaS3;
double *vThetaI3;
double *vThetaI3a;
double *vlambda3;
double *vThetaFC3;
double *vPsi3;
double *vThetaR3;
double *vSoilDepth3;
double *vSoilDepth3init;

double *vCrustFraction;
double *vKsatCrust;
double *vPoreCrust;
double *vCompactFraction;
double *vKsatCompact;
double *vPoreCompact;
double *vGrassFraction;
double *vKsatGrass;
double *vPoreGrass;

double *vLw;
double *vInfilVol;
double *vKsateff;
double *vThetaeff;
double *vPoreeff;
double *vPerc;

// interception
double *vInterc;
double *vCStor;
double *vCanopyStorage;
double *vCover;
double *vLai;
double *vkLAI;
double *vLeafDrain;
double *vLitter;
double *vLCStor;
double *vLInterc;
double *vHStor;
double *vRoofStore;
double *vIntercHouse;
double *vDrumStore;
double *vDStor;
double *vIntercETa;

//double *vWH;
//double *vWHrunoff;
//double *vWHroad;
//double *vFloodDomain;
//double *vhmx;
//double *vDX;
//double *vSoilWidthDX;




