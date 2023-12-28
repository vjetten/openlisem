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
double TWorld::MapTotal1D(double *V)
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
double* TWorld::NewMap1D(double value)
{
    double* V = new double[nrValidCells];

    long i = 0;
    for (int r = 0; r < LDD->nrRows(); r++)
        for (int c = 0; c < LDD->nrCols(); c++) {
            if (!pcr::isMV(LDD->Drc))
                V[i] = value;
                i++;
        }
    return V;
   // qDebug() << (double)V.size()/(double)(_nrRows*_nrCols);
}
//---------------------------------------------------------------------------
double* TWorld::ReadMap1D(cTMap *Mask, QString name)
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

    double* V = new double[nrValidCells];
    long i = 0;
    for (int r = 0; r < Mask->nrRows(); r++)
        for (int c = 0; c < Mask->nrCols(); c++) {
                V[i] = _M->Drc;
                i++;
            }


    delete _M;

    return V;
}
//---------------------------------------------------------------------------
void TWorld::checkMap1D(double *V,int oper,double value, QString mapName, QString SS)
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

    vtma = NewMap1D(0);
    vtmb = NewMap1D(0);
    vtmc = NewMap1D(0);

    if(InfilMethod != INFIL_SWATRE) {

        vSoilDepth1 = ReadMap1D(LDD,getvaluename("soildep1"));
        vSoilDepth1init = NewMap1D(0);
        FOR_ROW_COL_MV_V {
            vSoilDepth1[i_] /= 1000.0;
            vSoilDepth1[i_] *= SD1Calibration;
            vSoilDepth1init[i_] = vSoilDepth1[i_];
        }

        vThetaS1 = ReadMap1D(LDD,getvaluename("thetas1"));
        vThetaI1 = ReadMap1D(LDD,getvaluename("thetai1"));
        vThetaI1a = NewMap1D(0); // used for screen output
        FOR_ROW_COL_MV_V {
            vThetaI1[i_] *= thetaCalibration;
            vThetaI1[i_] = std::min(vThetaI1[i_], vThetaS1[i_]);
            vThetaI1a[i_] = vThetaI1[i_];
        }

        vKsat1 = ReadMap1D(LDD,getvaluename("ksat1"));
        vThetaR1 = NewMap1D(0);
        vlambda1 = NewMap1D(0);
        vThetaFC1 = NewMap1D(0);

        FOR_ROW_COL_MV_V {
            double ks = std::max(0.5,std::min(1000.0,log(vKsat1[i_])));
            vlambda1[i_] = 0.0849*ks+0.159;
            vlambda1[i_] = std::min(std::max(0.1,vlambda1[i_]),0.7);
            vtma[i_] = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m   //psi ae
            vThetaR1[i_] = 0.0673*exp(-0.238*log(ks));
            vThetaFC1[i_] = -0.0519*log(ks) + 0.3714;
        }

        if (SwitchPsiUser) {
            vPsi1 = ReadMap1D(LDD, getvaluename("psi1"));
            FOR_ROW_COL_MV_V {
                vPsi1[i_] *= 0.01; // to m
            }
        } else {
            vPsi1 = NewMap1D(0);
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
            vSoilDepth2 = ReadMap1D(LDD,getvaluename("soildep2"));
            vSoilDepth2init = NewMap1D(0);
            FOR_ROW_COL_MV_V {
                vSoilDepth2[i_] /= 1000.0;
                vSoilDepth2[i_] *= SD2Calibration;
                vSoilDepth2init[i_] = vSoilDepth2[i_];
            }


            vThetaS2 = ReadMap1D(LDD,getvaluename("thetas2"));
            vThetaI2 = ReadMap1D(LDD,getvaluename("thetai2"));
            vThetaI2a = NewMap1D(0); // used for screen output
            FOR_ROW_COL_MV_V {
                vThetaI2[i_] *= thetaCalibration;
                vThetaI2[i_] = std::min(vThetaI2[i_], vThetaS2[i_]);
                vThetaI2a[i_] = vThetaI2[i_];
            }

            vKsat2 = ReadMap1D(LDD,getvaluename("ksat2"));
            vThetaR2 = NewMap1D(0);
            vlambda2 = NewMap1D(0);
            vThetaFC2 = NewMap1D(0);

            FOR_ROW_COL_MV_V {
                double ks = std::max(0.5,std::min(1000.0,log(vKsat2[i_])));
                vlambda2[i_] = 0.0849*ks+0.159;
                vlambda2[i_] = std::min(std::max(0.1,vlambda2[i_]),0.7);
                vtma[i_] = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m   //psi ae
                vThetaR2[i_] = 0.0673*exp(-0.238*log(ks));
                vThetaFC2[i_] = -0.0519*log(ks) + 0.3714;
            }

            if (SwitchPsiUser) {
                vPsi2 = ReadMap1D(LDD, getvaluename("psi2"));
                FOR_ROW_COL_MV_V {
                    vPsi2[i_] *= 0.01; // to m
                }
            } else {
                vPsi2 =
                    NewMap1D(0);
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
            vSoilDepth3 = ReadMap1D(LDD,getvaluename("soildep3"));
            vSoilDepth3init = NewMap1D(0);
            FOR_ROW_COL_MV_V {
                vSoilDepth3[i_] /= 1000.0;
                vSoilDepth3[i_] *= SD2Calibration;
                vSoilDepth3init[i_] = vSoilDepth3[i_];
            }


            vThetaS3 = ReadMap1D(LDD,getvaluename("thetas3"));
            vThetaI3 = ReadMap1D(LDD,getvaluename("thetai3"));
            vThetaI3a = NewMap1D(0); // used for screen output
            FOR_ROW_COL_MV_V {
                vThetaI3[i_] *= thetaCalibration;
                vThetaI3[i_] = std::min(vThetaI3[i_], vThetaS3[i_]);
                vThetaI3a[i_] = vThetaI3[i_];
            }

            vKsat3 = ReadMap1D(LDD,getvaluename("ksat3"));
            vThetaR3 = NewMap1D(0);
            vlambda3 = NewMap1D(0);
            vThetaFC3 = NewMap1D(0);

            FOR_ROW_COL_MV_V {
                double ks = std::max(0.5,std::min(1000.0,log(vKsat3[i_])));
                vlambda3[i_] = 0.0849*ks+0.159;
                vlambda3[i_] = std::min(std::max(0.1,vlambda3[i_]),0.7);
                vtma[i_] = exp( -0.3012*ks + 3.5164) * 0.01; // 0.01 to convert to m   //psi ae
                vThetaR3[i_] = 0.0673*exp(-0.238*log(ks));
                vThetaFC3[i_] = -0.0519*log(ks) + 0.3714;
            }

            if (SwitchPsiUser) {
                vPsi3 = ReadMap1D(LDD, getvaluename("psi3"));
                FOR_ROW_COL_MV_V {
                    vPsi3[i_] *= 0.01; // to m
                }
            } else {
                vPsi3 = NewMap1D(0);
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
            vCrustFraction = ReadMap1D(LDD, getvaluename("crustfrc"));
            vKsatCrust = ReadMap1D(LDD, getvaluename("ksatcrst"));
            vPoreCrust = ReadMap1D(LDD, getvaluename("porecrst"));
            FOR_ROW_COL_MV_V {
                vCrustFraction[i_] = std::min(1.0, vCrustFraction[i_]);
                vKsatCrust[i_] *= _dt/3600000.0;
            }
        }

        if (SwitchInfilCompact)
        {
            vKsatCompact = ReadMap1D(LDD,getvaluename("ksatcomp"));
            vPoreCompact = ReadMap1D(LDD,getvaluename("porecomp"));
            vCompactFraction = ReadMap1D(LDD,getvaluename("compfrc"));
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
            vKsatGrass = ReadMap1D(LDD,getvaluename("ksatgras"));
            vPoreGrass = ReadMap1D(LDD,getvaluename("poregras"));
            vGrassFraction = ReadMap1D(LDD,getvaluename("grasfrc"));
            FOR_ROW_COL_MV_V {
                vGrassFraction[i_] = std::min(1.0, vGrassFraction[i_]);
                vKsatGrass[i_] *= _dt/3600000.0;
            }
        }
    }

    vLw = NewMap1D(0);
    //NewMap1D(vInfilVol, 0);

    vPerc =NewMap1D( 0);
    vPoreeff = NewMap1D( 0);
    vThetaeff = NewMap1D(0);
    vKsateff = NewMap1D(0);

}

void TWorld::InitLULCInput1D(void)
{
    if (!Switch1Darrays)
        return;

    vLai = ReadMap1D(LDD,getvaluename("lai"));
    vCover = ReadMap1D(LDD,getvaluename("cover"));

    checkMap1D(vLai, SMALLER, 0.0, getvaluename("lai"),"LAI must be >= 0");
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
        vCanopyStorage = NewMap1D(0); //in m !!!
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
        vCanopyStorage = ReadMap1D(LDD, getvaluename("smax"));
    }

    // openness coefficient k
    vkLAI = NewMap1D(0);
    FOR_ROW_COL_MV_V {
        vCanopyStorage[i_] *= SmaxCalibration;
        vCanopyStorage[i_] *= 0.001; // mm to m
        vkLAI[i_] = 1-exp(-CanopyOpeness*vtma[i_]);
    }
    vInterc = NewMap1D(0.0);
    vCStor = NewMap1D(0.0);
    vCanopyStorage = NewMap1D(0.0);
    vLeafDrain = NewMap1D(0.0);

    if (SwitchLitter)
    {
        vLInterc = NewMap1D(0.0);
        vLCStor = NewMap1D(0.0);
        vLitter = ReadMap1D(LDD,getvaluename("litter"));
        checkMap1D(vLitter, SMALLER, 0.0,getvaluename("litter"),"Litter cover fraction must be >= 0");
        checkMap1D(vLitter, LARGER, 1.0, getvaluename("litter"),"Litter cover fraction must be <= 1.0");
        LitterSmax = getvaluedouble("Litter interception storage");
    }

    if (SwitchHouses)
    {
        vIntercHouse = NewMap1D(0.0);
        vHStor = NewMap1D(0.0);
        vRoofStore = ReadMap1D(LDD,getvaluename("roofstore"));
        FOR_ROW_COL_MV_V {
            vRoofStore[i_] *= 0.001; // mm to m
        }

        if (SwitchRaindrum) {
            vDrumStore = ReadMap1D(LDD,getvaluename("drumstore"));
            vDStor = NewMap1D(0.0);
        }
    }

    if (SwitchIncludeET)
        vIntercETa = NewMap1D(0.0);

}
