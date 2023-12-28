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

#include <algorithm>
#include <qstring.h>
#include "io.h"
#include "lisemqt.h"
#include "global.h"

#include "model.h"
#include "operation.h"
#include "CsfRGBMap.h"

//---------------------------------------------------------------------------
double TWorld::MapTotal1D(QVector <double> &V)
{
    double total = 0.;
    #pragma omp parallel for reduction(+:total) num_threads(userCores)
    FOR_ROW_COL_MV_V {
        if (!pcr::isMV(V[i_]))
            total = total + V[i_];
    }
    return (total);
}
//---------------------------------------------------------------------------
void TWorld::NewMap1D(QVector <double> &V, double value)
{
    V.clear();
    for (int r = 0; r < LDD->nrRows(); r++)
        for (int c = 0; c < LDD->nrCols(); c++) {
            if (!pcr::isMV(LDD->Drc))
                V << 0.0;
        }
}
//---------------------------------------------------------------------------
void TWorld::ReadMap1D(cTMap *Mask, QVector <double> &V, QString name)
{
    cTMap *_M = new cTMap(readRaster(name));

    for (int r = 0; r < Mask->nrRows(); r++)
        for (int c = 0; c < Mask->nrCols(); c++) {
            if (!pcr::isMV(Mask->Drc) && pcr::isMV(_M->Drc)) {
                QString sr, sc;
                sr.setNum(r); sc.setNum(c);
                ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+name+".\n \
                                                                                                     This is a cell with missing values where a flow network esists (either LDD, Channel LDD, tile drain LDD).";
                                                                                                     throw 1;
            }
        }


    for (int r = 0; r < Mask->nrRows(); r++)
        for (int c = 0; c < Mask->nrCols(); c++) {
            if (!pcr::isMV(Mask->Drc)) {
                if (pcr::isMV(_M->Drc))
                    V << 0.0;
                else
                    V << _M->Drc;
            }
        }

    delete _M;

    if (V.size() != nrCells) {
        ErrorString = "Length of Vector in map: "+name+" is wrong.";
        throw 1;
    }

}
//---------------------------------------------------------------------------
void TWorld::checkMap1D(QVector <double> &V,int oper,double value, QString mapName, QString SS)
{
    FOR_ROW_COL_MV_L {
        if (oper == LARGER && V[i_] > value)
        {
            ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger than %4.\n").arg(r).arg(c).arg(mapName).arg(value) + SS;
            throw 1;
        }
        else
            if (oper == SMALLER && V[i_] < value)
            {
                ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller than %4.\n").arg(r).arg(c).arg(mapName).arg(value) + SS;
                throw 1;
            }
            else
                if (oper == LARGEREQUAL && V[i_] >= value)
                {
                    ErrorString = QString("Value at row=%1 and col=%2 in %3 is larger or equal than %4.\n").arg(r).arg(c).arg(mapName).arg(value) + SS;
                    throw 1;
                }
                else
                    if (oper == SMALLEREQUAL && V[i_] <= value)
                    {
                        ErrorString = QString("Value at row=%1 and col=%2 in %3 is smaller or equal than %4.\n").arg(r).arg(c).arg(mapName).arg(value) + SS;
                        throw 1;
                    }
    }}
}
//---------------------------------------------------------------------------

void TWorld::InitSoilInput1D(void)
{
    if (!Switch1Darrays)
        return;

    NewMap1D(vtma,0);
    NewMap1D(vtmb,0);
    NewMap1D(vtmc,0);
    NewMap1D(vtmd,0);

    if(InfilMethod != INFIL_SWATRE) {

        ReadMap1D(LDD,vSoilDepth1,getvaluename("soildep1"));
        NewMap1D(vSoilDepth1init,0);
        FOR_ROW_COL_MV_V {
            vSoilDepth1[i_] /= 1000.0;
            vSoilDepth1[i_] *= SD1Calibration;
            vSoilDepth1init[i_] = vSoilDepth1[i_];
        }

        ReadMap1D(LDD,vThetaS1,getvaluename("thetas1"));
        ReadMap1D(LDD,vThetaI1,getvaluename("thetai1"));
        NewMap1D(vThetaI1a,0); // used for screen output
        FOR_ROW_COL_MV_V {
            vThetaI1[i_] *= thetaCalibration;
            vThetaI1[i_] = std::min(vThetaI1[i_], vThetaS1[i_]);
            vThetaI1a[i_] = vThetaI1[i_];
        }

        ReadMap1D(LDD,vKsat1,getvaluename("ksat1"));
        NewMap1D(vThetaR1, 0);
        NewMap1D(vlambda1, 0);
        NewMap1D(vThetaFC1, 0);

        FOR_ROW_COL_MV_V {
            double ks = std::max(0.5,std::min(1000.0,log(vKsat1[i_])));
            vlambda1[i_] = 0.0849*ks+0.159;
            vlambda1[i_] = std::min(std::max(0.1,vlambda1[i_]),0.7);
            vtma[i_] = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m   //psi ae
            vThetaR1[i_] = 0.0673*exp(-0.238*log(ks));
            vThetaFC1[i_] = -0.0519*log(ks) + 0.3714;
        }

        if (SwitchPsiUser) {
            ReadMap1D(LDD, vPsi1, getvaluename("psi1"));
            FOR_ROW_COL_MV_V {
                vPsi1[i_] *= 0.01; // to m
            }
        } else {
            NewMap1D(vPsi1, 0);
            FOR_ROW_COL_MV_V {
                vPsi1[i_] = exp(-0.3382*log(vKsat1[i_]) + 3.3425)*0.01;
                vPsi1[i_] = std::max(vPsi1[i_],vtma[i_]);//psi1ae[i_]);
            }
        }

        FOR_ROW_COL_MV_V {
            vKsat1[i_] *= ksatCalibration;
            vKsat1[i_] *= _dt/3600000.0;
        }

        if (SwitchTwoLayer)
        {
            ReadMap1D(LDD,vSoilDepth2,getvaluename("soildep2"));
            NewMap1D(vSoilDepth2init,0);
            FOR_ROW_COL_MV_V {
                vSoilDepth2[i_] /= 1000.0;
                vSoilDepth2[i_] *= SD2Calibration;
                vSoilDepth2init[i_] = vSoilDepth2[i_];
            }


            ReadMap1D(LDD,vThetaS2,getvaluename("thetas2"));
            ReadMap1D(LDD,vThetaI2,getvaluename("thetai2"));
            NewMap1D(vThetaI2a,0); // used for screen output
            FOR_ROW_COL_MV_V {
                vThetaI2[i_] *= thetaCalibration;
                vThetaI2[i_] = std::min(vThetaI2[i_], vThetaS2[i_]);
                vThetaI2a[i_] = vThetaI2[i_];
            }

            ReadMap1D(LDD,vKsat2,getvaluename("ksat2"));
            NewMap1D(vThetaR2, 0);
            NewMap1D(vlambda2, 0);
            NewMap1D(vThetaFC2, 0);

            FOR_ROW_COL_MV_V {
                double ks = std::max(0.5,std::min(1000.0,log(vKsat2[i_])));
                vlambda2[i_] = 0.0849*ks+0.159;
                vlambda2[i_] = std::min(std::max(0.1,vlambda2[i_]),0.7);
                vtma[i_] = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m   //psi ae
                vThetaR2[i_] = 0.0673*exp(-0.238*log(ks));
                vThetaFC2[i_] = -0.0519*log(ks) + 0.3714;
            }

            if (SwitchPsiUser) {
                ReadMap1D(LDD, vPsi2, getvaluename("psi2"));
                FOR_ROW_COL_MV_V {
                    vPsi2[i_] *= 0.01; // to m
                }
            } else {
                NewMap1D(vPsi2, 0);
                FOR_ROW_COL_MV_V {
                    vPsi2[i_] = exp(-0.3382*log(vKsat2[i_]) + 3.3425)*0.01;
                    vPsi2[i_] = std::max(vPsi2[i_],vtma[i_]);//psi1ae[i_]);
                }
            }

            FOR_ROW_COL_MV_V {
                vKsat2[i_] *= ksat2Calibration;
                vKsat2[i_] *= _dt/3600000.0;

            }
        } // 2 layer

        if (SwitchThreeLayer)
        {
            ReadMap1D(LDD,vSoilDepth3,getvaluename("soildep3"));
            NewMap1D(vSoilDepth3init,0);
            FOR_ROW_COL_MV_V {
                vSoilDepth3[i_] /= 1000.0;
                vSoilDepth3[i_] *= SD2Calibration;
                vSoilDepth3init[i_] = vSoilDepth3[i_];
            }


            ReadMap1D(LDD,vThetaS3,getvaluename("thetas3"));
            ReadMap1D(LDD,vThetaI3,getvaluename("thetai3"));
            NewMap1D(vThetaI3a,0); // used for screen output
            FOR_ROW_COL_MV_V {
                vThetaI3[i_] *= thetaCalibration;
                vThetaI3[i_] = std::min(vThetaI3[i_], vThetaS3[i_]);
                vThetaI3a[i_] = vThetaI3[i_];
            }

            ReadMap1D(LDD,vKsat3,getvaluename("ksat3"));
            NewMap1D(vThetaR3, 0);
            NewMap1D(vlambda3, 0);
            NewMap1D(vThetaFC3, 0);

            FOR_ROW_COL_MV_V {
                double ks = std::max(0.5,std::min(1000.0,log(vKsat3[i_])));
                vlambda3[i_] = 0.0849*ks+0.159;
                vlambda3[i_] = std::min(std::max(0.1,vlambda3[i_]),0.7);
                vtma[i_] = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m   //psi ae
                vThetaR3[i_] = 0.0673*exp(-0.238*log(ks));
                vThetaFC3[i_] = -0.0519*log(ks) + 0.3714;
            }

            if (SwitchPsiUser) {
                ReadMap1D(LDD, vPsi3, getvaluename("psi3"));
                FOR_ROW_COL_MV_V {
                    vPsi3[i_] *= 0.01; // to m
                }
            } else {
                NewMap1D(vPsi3, 0);
                FOR_ROW_COL_MV_V {
                    vPsi3[i_] = exp(-0.3382*log(vKsat3[i_]) + 3.3425)*0.01;
                    vPsi3[i_] = std::max(vPsi3[i_],vtma[i_]);//psi1ae[i_]);
                }
            }

            FOR_ROW_COL_MV_V {
                vKsat3[i_] *= ksat3Calibration;
                vKsat3[i_] *= _dt/3600000.0;

            }
        } // 3 layer

        if (SwitchInfilCrust) {
            ReadMap1D(LDD, vCrustFraction, getvaluename("crustfrc"));
            ReadMap1D(LDD, vKsatCrust, getvaluename("ksatcrst"));
            ReadMap1D(LDD, vPoreCrust, getvaluename("porecrst"));
            FOR_ROW_COL_MV_V {
                vCrustFraction[i_] = std::min(1.0, vCrustFraction[i_]);
                vKsatCrust[i_] *= _dt/3600000.0;
            }
        }

        if (SwitchInfilCompact)
        {
            ReadMap1D(LDD,vKsatCompact,getvaluename("ksatcomp"));
            ReadMap1D(LDD,vPoreCompact,getvaluename("porecomp"));
            ReadMap1D(LDD,vCompactFraction,getvaluename("compfrc"));
            FOR_ROW_COL_MV_V {
                vCompactFraction[i_] = std::min(1.0, vCompactFraction[i_]);
                vKsatCompact[i_] *= _dt/3600000.0;
            }
        }

        if (SwitchInfilCompact && SwitchInfilCrust) {
            FOR_ROW_COL_MV_V {
                if (vCrustFraction[i_] + vCompactFraction[i_] > 1.0) {
                    vCrustFraction[i_] = 1.0-vCompactFraction[i_];
                }
            }
        }

        if (SwitchGrassStrip)
        {
            ReadMap1D(LDD,vKsatGrass,getvaluename("ksatgras"));
            ReadMap1D(LDD,vPoreGrass,getvaluename("poregras"));
            ReadMap1D(LDD,vGrassFraction,getvaluename("grasfrc"));
            FOR_ROW_COL_MV_V {
                vGrassFraction[i_] = std::min(1.0, vGrassFraction[i_]);
                vKsatGrass[i_] *= _dt/3600000.0;
            }
        }
    }

    NewMap1D(vLw, 0);
    NewMap1D(vInfilVol, 0);

    NewMap1D(vPerc, 0);
    NewMap1D(vPoreeff, 0);
    NewMap1D(vThetaeff, 0);
    NewMap1D(vKsateff, 0);

}

void TWorld::InitLULCInput1D(void)
{
    if (!Switch1Darrays)
        return;

    ReadMap1D(LDD,vtma,getvaluename("lai"));
    ReadMap1D(LDD,vCover,getvaluename("cover"));

    checkMap1D(vtma, SMALLER, 0.0, getvaluename("lai"),"LAI must be >= 0");
    checkMap1D(vCover, SMALLER, 0.0, getvaluename("cover"),"Cover fraction must be >= 0");
    checkMap1D(vCover, LARGER, 1.0, getvaluename("cover"),"Cover fraction must be <= 1.0");

    if (SwitchGrassStrip) {
        FOR_ROW_COL_MV_V {
            if (vGrassFraction[i_] > 0)
            {
                vCover[i_] = vCover[i_]*(1-vGrassFraction[i_]) + 0.95*vGrassFraction[i_];
                vtma[i_] = vtma[i_]*(1-vGrassFraction[i_]) + 5.0*vGrassFraction[i_];
            }
        }
    }


    if (SwitchInterceptionLAI)
    {
        NewMap1D(vCanopyStorage, 0); //in m !!!
        FOR_ROW_COL_MV_V {
            switch (InterceptionLAIType) {
                case 0: vCanopyStorage[i_] = 0.4376 * vtma[i_] + 1.0356;break; // gives identical results
                case 1: vCanopyStorage[i_] = 0.2331 * vtma[i_]; break;
                case 2: vCanopyStorage[i_] = 0.3165 * vtma[i_]; break;
                case 3: vCanopyStorage[i_] = 1.46 * pow(vtma[i_],0.56); break;
                case 4: vCanopyStorage[i_] = 0.0918 * pow(vtma[i_],1.04); break;
                case 5: vCanopyStorage[i_] = 0.2856 * vtma[i_]; break;
                case 6: vCanopyStorage[i_] = 0.1713 * vtma[i_]; break;
                case 7: vCanopyStorage[i_] = 0.59 * pow(vtma[i_],0.88); break;
            }
        }
    } else {
        ReadMap1D(LDD, vCanopyStorage, getvaluename("smax"));
    }

    // openness coefficient k
    NewMap1D(vkLAI, 0);
    FOR_ROW_COL_MV_V {
        vCanopyStorage[i_] *= SmaxCalibration;
        vCanopyStorage[i_] *= 0.001; // mm to m
        vkLAI[i_] = 1-exp(-CanopyOpeness*vtma[i_]);
    }
    NewMap1D(vInterc,0.0);
    NewMap1D(vCStor,0.0);
    NewMap1D(vCanopyStorage,0.0);
    NewMap1D(vLeafDrain,0.0);

    if (SwitchLitter)
    {
        NewMap1D(vLInterc,0.0);
        NewMap1D(vLCStor,0.0);
        ReadMap1D(LDD,vLitter,getvaluename("litter"));
        checkMap1D(vLitter, SMALLER, 0.0,getvaluename("litter"),"Litter cover fraction must be >= 0");
        checkMap1D(vLitter, LARGER, 1.0, getvaluename("litter"),"Litter cover fraction must be <= 1.0");
        LitterSmax = getvaluedouble("Litter interception storage");
    }

    if (SwitchHouses)
    {
        NewMap1D(vIntercHouse,0.0);
        NewMap1D(vHStor,0.0);
        ReadMap1D(LDD,vRoofStore,getvaluename("roofstore"));
        FOR_ROW_COL_MV_V {
            vRoofStore[i_] *= 0.001; // mm to m
        }

        if (SwitchRaindrum) {
            ReadMap1D(LDD,vDrumStore,getvaluename("drumstore"));
            NewMap1D(vDStor,0.0);
        }
    }

    if (SwitchIncludeET)
        NewMap1D(vIntercETa,0.0);

}
