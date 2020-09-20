/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
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
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisFlowBarriers.cpp
  \brief Get flowbarrier data, and make sure they are used at the correct times

functions: \n
 */

#include <memory>
#include "io.h"
#include "model.h"
#include "operation.h"
#include <QtCore/qstring.h>


void TWorld::GetFlowBarrierData(QString name)
{
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList BarrierTypes;

    if (!fi.exists())
    {
        ErrorString = "FlowBarrier file not found: " + name;
        throw 1;
    }

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.count()-1,1);
        // readLine also reads \n as a character on an empty line!
        if (!S.trimmed().isEmpty())
            BarrierTypes << S.trimmed();
    }
    fff.close();


    FBid.clear();

    FBHeightN.clear();
    FBHeightS.clear();
    FBHeightE.clear();
    FBHeightW.clear();

    FBTimeN.clear();
    FBTimeS.clear();
    FBTimeE.clear();
    FBTimeW.clear();

    for(int i = BarrierTypes.length() -1;  i > 0; i--)
    {
        if(BarrierTypes.at(i).at(0) == '/')
        {
            BarrierTypes.removeAt(i);
        }
    }


    for(int i = 0; i < BarrierTypes.length(); i++)
    {
        QStringList list = BarrierTypes.at(i).split(QRegExp("\\s+"),Qt::SkipEmptyParts );

        for(int j = 0; j < 9; j++ )
        {
            bool ok = true;
            if(j < list.length())
            {
                if(j == 0)
                {
                    FBid.append(list.at(j).toInt(&ok));
                    if(ok == false){ErrorString = "FlowBarrier file has a non-number value "; throw 1;}
                }else if( j == 1)
                {
                    FBHeightN.append(list.at(j).toDouble(&ok));
                    if(ok == false){ErrorString = "FlowBarrier file has a non-number value "; throw 1;}
                }else if( j == 2)
                {
                    FBHeightE.append(list.at(j).toDouble(&ok));
                    if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 3)
                {
                    FBHeightS.append(list.at(j).toDouble(&ok));
                    if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}

                }else if( j == 4)
                {
                    FBHeightW.append(list.at(j).toDouble(&ok));
                    if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 5)
                {
                        FBTimeN.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 6)
                {
                        FBTimeE.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 7)
                {
                        FBTimeS.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 8)
                {
                        FBTimeW.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }
            }else
            {
                if( j == 0)
                {
                        FBid.append(0);
                }else if( j == 1)
                {
                        FBHeightN.append(0.0);
                }else if( j == 2)
                {
                        FBHeightW.append(0.0);
                }else if( j == 3)
                {
                        FBHeightE.append(0.0);
                }else if( j == 4)
                {
                        FBHeightS.append(0.0);
                }else if( j == 5)
                {
                        FBTimeN.append(-1.0);
                }else if( j == 6)
                {
                        FBTimeS.append(-1.0);
                }else if( j == 7)
                {
                        FBTimeE.append(-1.0);
                }else if( j == 8)
                {
                        FBTimeW.append(-1.0);

                }

            }

        }

    }



}

void TWorld::SetFlowBarriers()
{
    if (!SwitchFlowBarriers)
        return;
    FOR_ROW_COL_MV_L {
        if(this->time > FlowBarrierNT->Drc && !(FlowBarrierNT->Drc < 0))
        {
            FlowBarrierN->Drc = 0;
        }
        if(this->time > FlowBarrierET->Drc && !(FlowBarrierET->Drc < 0))
        {
            FlowBarrierE->Drc = 0;
        }
        if(this->time > FlowBarrierST->Drc && !(FlowBarrierST->Drc < 0))
        {
            FlowBarrierS->Drc = 0;
        }
        if(this->time > FlowBarrierWT->Drc && !(FlowBarrierWT->Drc < 0))
        {
            FlowBarrierW->Drc = 0;
        }

    }}
}

double TWorld::FBW(double h, int r, int c, int dr, int dc)
{

    if(dr == 0 && dc == 0)
    {
        return 1.0;
    }

    if(OUTORMV(r+dr,c+dc))
    {
        return 0.0;
    }

    if((h + DEM->Drc) > (FB(r,c,dr,dc) + DEM->data[r + dr][c+dc]))
    {
        return 1.0;
    }else
    {
        return 0.0;
    }
}

double TWorld::FB(int r, int c, int rd, int cd)
{
    double dem = 0;
    if(rd == 0 && cd == 0)
    {
        return dem;
    }

    if(OUTORMV(r+rd,c+cd))
    {
        return 0.0;
    }

    if(rd == 1 && cd == 0)
    {
        return dem + (FlowBarrierS->Drc);
    }
    else if(rd == 1 && cd == 1)
    {
        return dem + (std::max(std::max(FlowBarrierS->Drc,FlowBarrierE->Drc),std::max(FB(r,c +cd,0,rd),FB(r+rd,c,cd,0))));
    }
    else if(rd == 0 && cd == 1)
    {
        return dem + (FlowBarrierE->Drc);
    }
    else if(rd == -1 && cd == 1)
    {
        return dem + (std::max(std::max(FlowBarrierN->Drc,FlowBarrierE->Drc),std::max(FB(r,c  +cd,0,rd),FB(r+rd,c,cd,0))));
    }
    else if(rd == -1 && cd == 0)
    {
        return dem + (FlowBarrierN->Drc);
    }
    else if(rd == -1 && cd == -1)
    {
        return dem + (std::max(std::max(FlowBarrierN->Drc,FlowBarrierW->Drc),std::max(FB(r,c  +cd,0,rd),FB(r+rd,c,cd,0))));
    }
    else if(rd == 0 && cd == -1)
    {
        return dem + (FlowBarrierW->Drc);
    }else if(rd == 1 && cd == -1)
    {
        return dem + (std::max(std::max(FlowBarrierS->Drc,FlowBarrierW->Drc),std::max(FB(r,c +cd,0,rd),FB(r+rd,c,cd,0))));
    }

    return 0;
}

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::DEMFB()
 * @brief Returns the digital elevation model height, with the addition of flow barriers
 *
 * @param r : row number
 * @param c : column number
 * @param rd : row direction (-1 for top, 1 for bottom)
 * @param cd : column direction (-1 for left, 1 for right)
 * @param addwh : include water height for overland flow
 * @return digital elevation model height, with the addition of flow barriers
 * @see K2DDEMA
 */
double TWorld::DEMFB(int r, int c, int rd, int cd, bool addwh)
{

    double wh = 0;
    double dem = 0;
    if(INSIDE(r+rd,c+cd))
    {
        if(!pcr::isMV(LDD->data[r+rd][c+cd]))
        {
            if(addwh)
            {
                wh = WHrunoff->data[r + rd][c + cd];
            }
            dem = DEM->data[r + rd][c + cd];
        }else
        {
            if(!pcr::isMV(LDD->data[r][c]))
            {
                wh = 0;
                dem = DEM->Drc;
            }else
            {
               return 0;
            }
        }

    }else if(INSIDE(r,c))
    {
        if(!pcr::isMV(LDD->data[r][c]))
        {
            wh = 0;
            dem = DEM->Drc;
        }else
        {
           return 0;
        }

    }else
    {
        return 0;
    }

    if(OUTORMV(r+rd,c+cd))
    {
        return dem;
    }


    if(rd == 0 && cd == 0)
    {
        return dem + wh;
    }

    if(rd == 1 && cd == 0)
    {
        return dem + std::max(wh,(FlowBarrierS->Drc));
    }
    else if(rd == 1 && cd == 1)
    {
        return dem + std::max(wh,(std::max(std::max(FlowBarrierS->Drc,FlowBarrierE->Drc),std::max(FB(r,c +cd,0,rd),FB(r+rd,c,cd,0)))));
    }
    else if(rd == 0 && cd == 1)
    {
        return dem + std::max(wh,(FlowBarrierE->Drc));
    }
    else if(rd == -1 && cd == 1)
    {
        return dem + std::max(wh,(std::max(std::max(FlowBarrierN->Drc,FlowBarrierE->Drc),std::max(FB(r,c  +cd,0,rd),FB(r+rd,c,cd,0)))));
    }
    else if(rd == -1 && cd == 0)
    {
        return dem + std::max(wh,(FlowBarrierN->Drc));
    }
    else if(rd == -1 && cd == -1)
    {
        return dem + std::max(wh,(std::max(std::max(FlowBarrierN->Drc,FlowBarrierW->Drc),std::max(FB(r,c  +cd,0,rd),FB(r+rd,c,cd,0)))));
    }
    else if(rd == 0 && cd == -1)
    {
        return dem + std::max(wh,(FlowBarrierW->Drc));
    }else if(rd == 1 && cd == -1)
    {
        return dem + std::max(wh,(std::max(std::max(FlowBarrierS->Drc,FlowBarrierW->Drc),std::max(FB(r,c +cd,0,rd),FB(r+rd,c,cd,0)))));
    }
    return 0;
}

//---------------------------------------------------------------------------
void TWorld::InitFlowBarriers(void)
{

    FlowBarrier = NewMap(0);

    FlowBarrierN = NewMap(0);
    FlowBarrierW = NewMap(0);
    FlowBarrierS = NewMap(0);
    FlowBarrierE = NewMap(0);

    FlowBarrierNT = NewMap(-1);
    FlowBarrierWT = NewMap(-1);
    FlowBarrierST = NewMap(-1);
    FlowBarrierET = NewMap(-1);

    if(SwitchFlowBarriers)
    {
        QString filename = getvaluename("Flow barrier table filename");

        GetFlowBarrierData(filename);

        FlowBarrier = ReadMap(LDD,getvaluename("flowbarrierindex"));


        for(int i = 0; i < FBid.length(); i++)
        {

            int index = FBid.at(i);
            if(index <= 0)
            {
                continue;
            }
            FOR_ROW_COL_MV
            {

                if(FlowBarrier->Drc == index)
                {
                    FlowBarrierN->Drc = FBHeightN.at(i);
                    FlowBarrierS->Drc = FBHeightS.at(i);
                    FlowBarrierE->Drc = FBHeightE.at(i);
                    FlowBarrierW->Drc = FBHeightW.at(i);

                    FlowBarrierNT->Drc = FBTimeN.at(i);
                    FlowBarrierST->Drc = FBTimeS.at(i);
                    FlowBarrierET->Drc = FBTimeE.at(i);
                    FlowBarrierWT->Drc = FBTimeW.at(i);
                }

            }
        }
    }
}
