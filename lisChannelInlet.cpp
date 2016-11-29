
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
  \file lisChannelInlet.cpp
  \brief read and process files related to channel inlets/outlets with specified discharge capacity

 */

#include <memory>
#include "io.h"
#include "model.h"
#include "operation.h"


void TWorld::InitSubInlets()
{
    DichargeInlets.clear();

    if((!SwitchIncludeChannel) || (!SwitchChannelSubInlets))
    {
        return;
    }

    qDebug() << "read sub inlets";

    //read map with id values
    SubInletID = ReadMap(LDDChannel,getvaluename("subinlets"));
    SubInletQtransfer = NewMap(0.0);
    SubInletVchange = NewMap(0.0);

    //read all relevant discharge files, and fill DichargeInlets data structure
    ReadChannelInletFiles();

    //finish up preperations

}

void TWorld::DestroySubInlets()
{
    if((!SwitchIncludeChannel) || (!SwitchChannelSubInlets))
    {
        return;
    }

    for(int i = 0; i < DichargeInlets.length(); i++)
    {
        DichargeInlets.at(i)->values.clear();
        delete DichargeInlets.at(i);
    }

    DichargeInlets.clear();

}


void TWorld::ReadChannelInletFiles()
{
    QList<double> found_id;
    QList<double> found_link;
    found_id.clear();
    found_link.clear();

    //for every found id value, we need to find a subinlet_id.txt file to read
    FOR_ROW_COL_MV_CH
    {
        if(SubInletID->Drc == 0)
        {
            continue;
        }

        bool found  = false;
        for(int i = 0; i < found_id.length(); i++)
        {
            if(std::fabs(found_id.at(i)) == std::fabs(SubInletID->Drc))
            {
                if(found_id.at(i) == SubInletID->Drc)
                {
                    ErrorString = "Duplicate index in subinlet map" + QString::number(SubInletID->Drc);
                    throw 1;
                }

                found = true;
                for(int j = 0; j < DichargeInlets.length();j++)
                {
                    if(DichargeInlets.at(j)->id == std::fabs(SubInletID->Drc))
                    {
                        DichargeInlets.at(j)->is_linked = true;
                        if(SubInletID->Drc < 0)
                        {
                            DichargeInlets.at(j)->r2 = r;
                            DichargeInlets.at(j)->c2 = c;
                        }else
                        {
                            DichargeInlets.at(j)->r = r;
                            DichargeInlets.at(j)->c = c;
                        }

                    }
                }

                break;
            }
        }

        if(!found)
        {
            qDebug() << " found inlet id " << SubInletID->Drc;
            //so we found a new id, add to list
            found_id.append(SubInletID->Drc);

            //read discharge file
            int pointid = std::fabs(SubInletID->Drc);
            InletData * data = new InletData();
            data->id = pointid;
            data->is_linked = false;
            if(SubInletID->Drc < 0)
            {
                data->r2 = r;
                data->c2 = c;
            }else
            {
                data->r = r;
                data->c = c;
            }
            ReadInletTimeSeriesFile(inputDir + "subinlet_" + QString::number(pointid,10) + ".txt",data);
            DichargeInlets.append(data);

        }
    }

}

void TWorld::ReadDischargeFile(int pointid)
{


}

void TWorld::ReadInletTimeSeriesFile(QString name,InletData * data)
{

    data->times.clear();
    data->values.clear();

    QFile fff(name);
    QFileInfo fi(name);

    if (!fi.exists())
    {
        ErrorString = "hahafile not found: " + name;
        throw 1;
    }

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    QString S;
    QStringList qRecs;
    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.count()-1,1);
        // readLine also reads \n as a character on an empty line!
        if (!S.trimmed().isEmpty())
            qRecs << S.trimmed();
        //qDebug() << S;
    }
    fff.close();
    bool ok = true;
    int nrcols = 0;
    if(qRecs.length() > 4)
    {
        nrcols = qRecs[1].toInt(&ok, 10);

        if(!ok)
        {
            ErrorString = "Cannot read number of columns in" + name;
            throw 1;
        }
    }else
    {
        ErrorString = "Cannot read/interpret file " + name;
        throw 1;
    }
    qDebug () << "file id subinlet" << data->id << "nrcols " << nrcols;
    for(int i = 1; i < nrcols; i++)
    {
        QList<double> * l  = new QList<double>();
        data->values.append(l);
    }

    int nrskip = 2 + nrcols;
    for(int r = nrskip; r < (qRecs.length()); r++)
    {
        QStringList SL = qRecs[r].split(QRegExp("\\s+"));
        if(SL.length() != nrcols)
        {
            ErrorString = "Wrong numbr of columns in line " + QString::number(r) + " in file " +name;
            throw 1;
        }

        data->times.append(SL.at(0).toDouble(&ok));
        if(!ok)
        {
            ErrorString = "Cannot read/interpret number in file " + name;
        }
        for(int i = 1; i < nrcols; i++)
        {
            qDebug() << "point " << data->id << "time " << SL.at(0) << "value" << SL.at(1);

            data->values.at(i-1)->append(SL.at(i).toDouble(&ok));
            if(!ok)
            {
                ErrorString = "Cannot read/interpret number in file " + name;
            }
        }
    }
}

double TWorld::GetDischargeAtTime(int index, double _time, double novalue)
{
    _time = _time/60;
    if(!(index < DichargeInlets.length()))
    {
        return novalue;
    }


    if(DichargeInlets.at(index)->values.at(0)->length() > 0)
    {
        double rval = novalue;

        for(int i = 0; i < DichargeInlets.at(index)->times.length(); i++)
        {
            if(DichargeInlets.at(index)->times.at(i) <= _time)
            {
                if(i != DichargeInlets.at(index)->times.length())
                {
                    if(DichargeInlets.at(index)->times.at(i+1) > _time)
                    {
                        rval = DichargeInlets.at(index)->values.at(0)->at(i);
                    }
                }else
                {
                    rval = DichargeInlets.at(index)->values.at(0)->at(DichargeInlets.at(index)->values.at(0)->length()-1);
                }
            }
        }

        return rval;
    }

    return novalue;

}

double TWorld::GetDischargeOutlet(double q_old, int r, int c)
{
    if((!SwitchIncludeChannel) || (!SwitchChannelSubInlets))
    {
        return q_old;
    }

    if(SubInletID->Drc  == 0)
    {
        return q_old;
    }

    for(int i = 0; i < DichargeInlets.length() ; i++)
    {
        if(std::fabs(SubInletID->Drc) == DichargeInlets.at(i)->id)
        {
            double qval = GetDischargeAtTime((i),time,0);

            if(DichargeInlets.at(i)->is_linked)
            {
                if(qval > 0)
                {
                    if(r == DichargeInlets.at(i)->r && c == DichargeInlets.at(i)->c)
                    {

                         double temp = q_old;
                         q_old = std::max(0.0,q_old - std::fabs(qval));
                         //SubInletVchange->data[r][c] -= _dt * (temp - q_old);
                         SubInletQtransfer->data[DichargeInlets.at(i)->r2][DichargeInlets.at(i)->c2] += (temp - q_old);

                    }else if(r == DichargeInlets.at(i)->r2 && c == DichargeInlets.at(i)->c2)
                    {



                    }
                }else
                {
                    if(r == DichargeInlets.at(i)->r2 && c == DichargeInlets.at(i)->c2)
                    {
                        double temp = q_old;
                        q_old = std::max(0.0,q_old - std::fabs(qval));
                        //SubInletVchange->data[r][c] -= _dt * (temp-q_old);
                        SubInletQtransfer->data[DichargeInlets.at(i)->r][DichargeInlets.at(i)->c] += (temp-q_old);


                    }else if(r == DichargeInlets.at(i)->r && c == DichargeInlets.at(i)->c)
                    {


                    }
                }

            }else
            {
                if(qval > 0)
                {

                }else
                {
                    //SubInletVchange->data[r][c] -= _dt *std::min(q_old,std::fabs(qval));
                    q_old = std::max(0.0,q_old - std::fabs(qval));

                }
            }

        }

    }

    return std::max(0.0,q_old);
}

double TWorld::GetDischargeInlet(double q_old, int r, int c)
{
    if((!SwitchIncludeChannel) || (!SwitchChannelSubInlets))
    {
        return q_old;
    }
    if(SubInletID->Drc  == 0)
    {
        return q_old;
    }

    for(int i = 0; i < DichargeInlets.length() ; i++)
    {
        if(std::fabs(SubInletID->Drc) == DichargeInlets.at(i)->id)
        {
            double qval = GetDischargeAtTime((i),time,0);

            if(DichargeInlets.at(i)->is_linked)
            {
                if(qval > 0)
                {
                    if(r == DichargeInlets.at(i)->r && c == DichargeInlets.at(i)->c)
                    {

                    }else if(r == DichargeInlets.at(i)->r2 && c == DichargeInlets.at(i)->c2)
                    {
                        if(SubInletQtransfer->Drc > 0)
                        {
                            q_old = q_old + SubInletQtransfer->Drc;
                            //SubInletVchange->data[r][c] += _dt * SubInletQtransfer->Drc;
                        }

                        SubInletQtransfer->Drc = 0;


                    }
                }else
                {
                    if(r == DichargeInlets.at(i)->r2 && c == DichargeInlets.at(i)->c2)
                    {

                    }else if(r == DichargeInlets.at(i)->r && c == DichargeInlets.at(i)->c)
                    {
                        if(SubInletQtransfer->Drc > 0)
                        {
                            q_old = q_old + SubInletQtransfer->Drc;
                            //SubInletVchange->data[r][c] += _dt * SubInletQtransfer->Drc;
                        }

                        SubInletQtransfer->Drc = 0;

                    }
                }

            }else
            {
                if(qval > 0)
                {
                    //SubInletVchange->data[r][c] += _dt * (qval-q_old);
                    q_old = qval;
                }else
                {

                }
            }

        }

    }

    return std::max(0.0,q_old);
}


double TWorld::GetVolumeInlet(double v_old, int r, int c)
{

    //does not do anything anymore!!!
    //used to add subinlet water through volume change



    if((!SwitchIncludeChannel) || (!SwitchChannelSubInlets))
    {
        return v_old;
    }

    /*if(SubInletQtransfer->Drc > 0)
    {
        v_old = v_old + SubInletQtransfer->Drc;
        SubInletVchange->data[r][c] += _dt * SubInletQtransfer->Drc;

    }

    SubInletQtransfer->Drc = 0;*/

    return v_old;
}
