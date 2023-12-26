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

QVector <double> vtma;
QVector <double> vtmb;
QVector <double> vtmc;
QVector <double> vtmd;

// infiltration, percolation, redistribution
QVector <double> vKsat1;
QVector <double> vThetaS1;
QVector <double> vThetaI1;
QVector <double> vThetaI1a;
QVector <double> vThetaR1;
QVector <double> vPsi1;
QVector <double> vThetaFC1;
QVector <double> vlambda1;
QVector <double> vSoilDepth1;
QVector <double> vSoilDepth1init;

QVector <double> vKsat2;
QVector <double> vThetaS2;
QVector <double> vThetaI2;
QVector <double> vThetaI2a;
QVector <double> vlambda2;
QVector <double> vThetaFC2;
QVector <double> vPsi2;
QVector <double> vThetaR2;
QVector <double> vSoilDepth2;
QVector <double> vSoilDepth2init;

QVector <double> vKsat3;
QVector <double> vThetaS3;
QVector <double> vThetaI3;
QVector <double> vThetaI3a;
QVector <double> vlambda3;
QVector <double> vThetaFC3;
QVector <double> vPsi3;
QVector <double> vThetaR3;
QVector <double> vSoilDepth3;
QVector <double> vSoilDepth3init;

QVector <double> vCrustFraction;
QVector <double> vKsatCrust;
QVector <double> vPoreCrust;
QVector <double> vCompactFraction;
QVector <double> vKsatCompact;
QVector <double> vPoreCompact;
QVector <double> vGrassFraction;
QVector <double> vKsatGrass;
QVector <double> vPoreGrass;

QVector <double> vLw;
QVector <double> vInfilVol;
QVector <double> vKsateff;
QVector <double> vThetaeff;
QVector <double> vPoreeff;
QVector <double> vPerc;

// interception
QVector <double> vInterc;
QVector <double> vCStor;
QVector <double> vCanopyStorage;
QVector <double> vCover;
QVector <double> vkLAI;
QVector <double> vLeafDrain;
QVector <double> vLitter;
QVector <double> vLCStor;
QVector <double> vLInterc;
QVector <double> vHStor;
QVector <double> vRoofStore;
QVector <double> vIntercHouse;
QVector <double> vDrumStore;
QVector <double> vDStor;

//QVector <double> vWH;
//QVector <double> vWHrunoff;
//QVector <double> vWHroad;
//QVector <double> vFloodDomain;
//QVector <double> vhmx;
//QVector <double> vDX;
//QVector <double> vSoilWidthDX;




