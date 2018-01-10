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
#include <cassert>
#include "CsfMap.h"
#include "error.h"
#include "model.h"

void TWorld::InitInflow(void)
{
    InflowID = NewMap(0);

    if(SwitchInflow)
    {
        QString filename = getvaluename("Inflow table filename");

        GetInflowData(filename);

        InflowID = ReadMap(LDD,getvaluename("inflowid"));

    }

}

void TWorld::GetInflowData(QString name)
{
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList InflowTypes;

    if (!fi.exists())
    {
        ErrorString = "Inflow file not found: " + name;
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
            InflowTypes << S.trimmed();
    }
    fff.close();


    for(int i = InflowTypes.length() -1;  i > 0; i--)
    {
        if(InflowTypes.at(i).at(0) == '/')
        {
            InflowTypes.removeAt(i);
        }
    }


    IFTime.clear();
    IFQ.clear();

    int nr_id = 0;
    bool  ok = true;

    for(int i = 0; i < InflowTypes.length(); i++)
    {
        QStringList list = InflowTypes.at(i).split(QRegExp("\\s+"),QString::SplitBehavior::SkipEmptyParts );

        qDebug() << list;

        IFTime.append(list.at(0).toDouble(&ok));
        if(!ok) { ErrorString = "Not a Number: " + list.at(0);}

        if( i == 0)
        {
           nr_id = std::max(0,list.length() - 1);
        }else if( list.length() - 1 != nr_id)
        {
            ErrorString = "wrong length in file: " + name + " row " + QString::number(i);
            throw 1;
        }

        for(int j = 1; j < list.length(); j++ )
        {
            if(i == 0)
            {
                QList<double> * nlist = new QList<double>();
                IFQ.append(nlist);
            }

            IFQ.at(j-1)->append(list.at(j).toDouble(&ok));
            if(!ok) { ErrorString = "Not a Number: " + list.at(0);}
        }


    }

}


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
        QStringList list = BarrierTypes.at(i).split(QRegExp("\\s+"),QString::SplitBehavior::SkipEmptyParts );

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
                        FBHeightW.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 3)
                {
                        FBHeightE.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 4)
                {
                        FBHeightS.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 5)
                {
                        FBTimeN.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 6)
                {
                        FBTimeS.append(list.at(j).toDouble(&ok));
                        if(ok == false){ErrorString = "FlowBarrier file has a non-number value"; throw 1;}
                }else if( j == 7)
                {
                        FBTimeE.append(list.at(j).toDouble(&ok));
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
    FOR_ROW_COL_MV
    {
        if(this->time > FlowBarrierNT->Drc && !(FlowBarrierNT->Drc < 0))
        {
            FlowBarrierN->Drc = 0;
        }
        if(this->time > FlowBarrierST->Drc && !(FlowBarrierST->Drc < 0))
        {
            FlowBarrierS->Drc = 0;
        }
        if(this->time > FlowBarrierWT->Drc && !(FlowBarrierWT->Drc < 0))
        {
            FlowBarrierW->Drc = 0;
        }
        if(this->time > FlowBarrierET->Drc && !(FlowBarrierET->Drc < 0))
        {
            FlowBarrierE->Drc = 0;
        }
    }


}

double TWorld::GetFlowBarrierHeight(int r, int c, int rd, int cd)
{
    if(INSIDE(r+rd,c+cd))
        {
            if(!pcr::isMV(LDD->data[r+rd][c+cd]))
            {
            }else
            {
                if(!pcr::isMV(LDD->data[r][c]))
                {

                }else
                {
                   return 0.0;
                }
            }

        }else if(INSIDE(r,c))
        {
            if(!pcr::isMV(LDD->data[r][c]))
            {

            }else
            {
               return 0.0;
            }

        }else
        {
            return 0;
        }

        if(rd == 0 && cd == 0)
        {
            return 0.0;
        }

        if(rd == 1 && cd == 0)
        {
            return FlowBarrierS->Drc;
        }
        else if(rd == 1 && cd == 1)
        {
            return std::max(FlowBarrierS->Drc,FlowBarrierE->Drc);
        }
        else if(rd == 0 && cd == 1)
        {
            return FlowBarrierE->Drc;
        }
        else if(rd == -1 && cd == 1)
        {
            return std::max(FlowBarrierN->Drc,FlowBarrierE->Drc);
        }
        else if(rd == -1 && cd == 0)
        {
            return FlowBarrierN->Drc;
        }
        else if(rd == -1 && cd == -1)
        {
            return std::max(FlowBarrierN->Drc,FlowBarrierW->Drc);
        }
        else if(rd == 0 && cd == -1)
        {
            return FlowBarrierW->Drc;
        }else if(rd == 1 && cd == -1)
        {
            return std::max(FlowBarrierS->Drc,FlowBarrierW->Drc);
        }



}
