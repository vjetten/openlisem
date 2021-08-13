
/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/


#include "model.h"



//---------------------------------------------------------------------------
void TWorld::GetDischargeData(QString name)
{
    Q_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList QRecs;
    QStringList SL;
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    double time = 0.0;
    QVector <int> locationnnrsmap;
    QVector <int> locationnnrsrec;

    locationnnrsmap.clear();
    locationnnrsrec.clear();

    DischargeInSeries.clear();

    if (!fi.exists())
    {
        ErrorString = "Input discharge file not found: " + name;
        throw 1;
    }

    nrDischargeseries = 0;

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    // read the timeseries data
    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.count()-1,1);
        // readLine also reads \n as a character on an empty line!
        if (!S.trimmed().isEmpty())
            QRecs << S.trimmed();
    }
    fff.close();

    // count the nr of inflow points
    FOR_ROW_COL_MV_L {
        if (QinLocation->Drc > 0) {

            locationnnrsmap << (int) QinLocation->Drc;
            nrStations++;
        }
    }}

    int count = QRecs[1].toInt(&ok, 10); // nr of columns stated in the file

    if (nrStations != count - 1) {
        ErrorString = QString("Nr columns is in the discharge inflow file (%1) is not equal to nr locations in the discharge inflow point map").arg(name);
        throw 1;
    }

    nrSeries = QRecs.size() - nrStations - 3;  // nr of timesteps/values
    if (nrSeries <= 1) {
        ErrorString = QString("Records in the discharge inflow file <= 1: %1").arg(name);
        throw 1;
    }
    // check if all data rows have the correct nr of columns
    bool check = false;
    int err = 0;
    for (int i = count+2; i < QRecs.size(); i++) {
        SL = QRecs[i].split(QRegExp("\\s+"));
        int nrrecords = SL.size();
        if (nrrecords != count) {
            err = i;
            check = true;
            break;
        }
    }
    if (check) {
        ErrorString = QString("Nr columns is not equal to nr of records %1 in row in file: %2").arg(err).arg(name);
        throw 1;
    }


    for (int i = 3; i < nrStations+3; i++)
            locationnnrsrec << QRecs[i].toInt(&ok, 10);

    if(locationnnrsmap != locationnnrsrec) {
        ErrorString = QString("Discharge input location numbers in timeseries and map are not the same.").arg(err).arg(name);
        throw 1;
    }

    //done checking !

    FOR_ROW_COL_MV_L {
        if (QinLocation->Drc > 0) {
            LDD_COORloc cr;

            cr.loc = (int) QinLocation->Drc;
            cr.r = r;
            cr.c = c;

            // find col nr in timeseries for this location id
            for(int j=0; j < nrStations; j++)
                if(QRecs[j+3].toInt(&ok, 10) == cr.loc){
                    cr.nr = j+1;
                }
            crQin_ << cr;
        }
    }}

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.Qin.clear();

        QString dirname;
        dirname = rainFileDir;

        QStringList SL = QRecs[r+nrStations+3].split(QRegExp("\\s+"), Qt::SkipEmptyParts);
        // split rainfall record row with whitespace

        rl.time = SL[0].toDouble();
        if (r == 0)
            time = rl.time;

        if (r > 0 && rl.time <= time)
        {
            ErrorString = QString("Discharge input file at time %1 has unreadable value.").arg(rl.time);
            throw 1;
        }
        else
            time = rl.time;

        for (int i = 1; i <= nrStations; i++)
        {
            bool ok = false;
            rl.Qin << SL[i].toDouble(&ok);

            if (!ok)
            {
                ErrorString = QString("Discharge input file at time %1 has unreadable value.").arg(SL[0]);
                throw 1;
            }
        }
        DischargeInSeries << rl;
    }

    //add 0 to timestep of 1 year in munites
    rl.Qin.clear();
    rl.time = 1440*365;
    nrSeries++;
    for (int i = 1; i < nrStations; i++)
        rl.Qin << 0.0;
    DischargeInSeries << rl;

    nrDischargeseries = nrSeries;
}
//---------------------------------------------------------------------------

void TWorld::DischargeInflow(void)
{
    double timeminprev = (time-_dt) / 60; //prev time in minutes
    int  qplace;
    double tt = 3600000.0;
    Q_LIST rl;

    if (!SwitchChannelInflow)
        return;

    for (qplace = 0; qplace < nrDischargeseries; qplace++)
        if (timeminprev < DischargeInSeries[qplace].time)
            break;


    #pragma omp parallel for num_threads(userCores)
    for(int i_ = 0; i_ < crQin_.size(); i_++) {
        int r = crQin_[i_].r;
        int c = crQin_[i_].c;
        int loc = crQin_[i_].loc;
        int nr = crQin_[i_].nr;

//        for (int j = 0; j < DischargeInSeries[qplace].Qin.size(); j++) {
//            if(loc == DischargeInSeries.Qin[qplace])
//                Qinflow->Drc = DischargeInSeries[qplace].Qin[j]*_dt/tt;
//            break;
//        }
    }

    report(*Qinflow,"qinf");
}
