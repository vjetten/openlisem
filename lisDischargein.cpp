
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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/


#include "model.h"


//---------------------------------------------------------------------------
void TWorld::GetDischargeDataNew(QString name)
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

    if (!fi.exists())
    {
        ErrorString = "User defined input discharge file not found: " + name;
        throw 1;
    }

    nrDischargeseries = 0;
    DischargeSeries.clear();
    currentDischargerow = 0;

    // read rainfall text file
    fff.open(QIODevice::ReadOnly | QIODevice::Text);
    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\r\n"))
            S.remove(S.count()-2,2);
        if (S.contains("\n"))
            S.remove(S.count()-1,1);

        if (!S.trimmed().isEmpty())
            QRecs << S.trimmed();
    }
    fff.close();

    // check first if PCRaster graph format is present: header, number of vars, columns equal vars
    int count = QRecs[1].toInt(&ok, 10); // nr of cols in file
    // header
    // second line is only an integer
    if (ok)
    {
        SL = QRecs[count+2].split(QRegExp("\\s+"));
        // check nr of columns in file
        if (count != SL.count()) {
            ErrorString = "Error: The nr of columns/stations in the discharge file does not equal the number on the second row.";
            throw 1;
        }

        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as first column
    }

    // get station numbers from header, or fill in 1,2 ... n
    stationQID.clear();
    for (int i = 0; i < nrStations; i++) {
        SL = QRecs[i+3].split(QRegExp("\\s+"));
        int tmp = SL.last().toInt(&ok, 10);
        if (ok)
            stationQID << tmp;
        // if the header ends with a number, that number is the ID corresponding with the map
        else
            stationQID << i+1;
        // elese the position is the ID number

    }

    nrSeries = QRecs.size() - nrStations - 3;
    if (nrSeries <= 1)
    {
        ErrorString = "Discharge records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.Qin.clear();
        rl.stationnr.clear();
        int r_ = r+nrStations+3;

        // split rainfall record row with whitespace
        QStringList SL = QRecs[r_].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);
        time = rl.time;

        // check if time is increasing with next row
        if (r+1 < nrSeries) {
            QStringList SL1 = QRecs[r_+1].split(QRegExp("\\s+"), Qt::SkipEmptyParts);
            int time1 = getTimefromString(SL1[0]);
            if (time1 < time) {
                ErrorString = QString("Time in discharge records is not increasing from row %1 to %2. Check your file!").arg(r_).arg(r_+1);
                throw 1;
            }
        }

        // rainfall values in this row
        for (int i = 1; i <= nrStations; i++)
        {
            bool ok = false;

            rl.Qin << SL[i].toDouble(&ok);
            if (!ok)
            {
                ErrorString = QString("Discharge records at time %1 has an unreadable value: %2.").arg(SL[0]).arg(SL[i]);
                throw 1;
            }
            rl.stationnr << stationID.at(i-1);
        }

        DischargeSeries << rl;
    }

    nrDischargeseries = DischargeSeries.size();
}
//---------------------------------------------------------------------------
void TWorld::GetDischargeMapfromStations(void)
{
    double currenttime = (time)/60;
    bool sameQ = false;

    // from time t to t+1 the rain is the rain of t

    // if time is outside records then use map with zeros
    if (currenttime < RainfallSeries[0].time || currenttime > RainfallSeries[nrRainfallseries-1].time)
    {
        //DEBUG("run time outside rainfall records");
        #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                Rain->Drc = 0;
                Rainc->Drc = 0;
                RainNet->Drc = 0;
            }}
            return;
        }

    // where are we in the series
    int currentrow = 0;
    for (int j = 0; j < DischargeSeries.count(); j++) {
        if (currenttime >= DischargeSeries[j].time && currenttime < DischargeSeries[j+1].time) {
            currentrow = j;
            break;
        }
    }

    if (currentrow == currentDischargerow && currentrow > 0)
        sameQ = true;

    // get the next map from file
    if (!sameQ) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            QuserIn->Drc = 0;
            for (int k = 0; k < stationQID.size(); k++) {
                if ((int) DischargeUserPoints->Drc == DischargeSeries[currentrow].stationnr.at(k)) {
                    QuserIn->Drc = DischargeSeries[currentrow].Qin[k]; // assuming m3/s
                   // qDebug() <<  QuserIn->Drc << currentrow << k << DischargeSeries[currentrow].stationnr.at(k);
                }
            }
        }}
    }

    currentDischargerow = currentrow;

//    double hoi =MapTotal(*QuserIn);
//    qDebug() << hoi;

}

/*
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
        if (S.contains("\r\n"))
            S.remove(S.count()-2,2);
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

    int count = QRecs[1].toInt(&ok, 10);
    // nr of columns stated in the file, 1st col is time
    if (count-1 != nrStations) {
        ErrorString = "Dischsarge input file error: the nr stations in the file does not equal the number of stations in the map.";
        throw 1;
    }

    SL = QRecs[nrStations+3].split(QRegExp("\\s+"));
    // check nr of columns in file for the second data record
    if (count != SL.count()) {
        ErrorString = "Dischsarge input file error: the nr of columns in the file does not equal the number on the second row.";
        throw 1;
    }

    nrSeries = QRecs.size() - nrStations - 3;  // nr of timesteps/values
    if (nrSeries <= 1) {
        ErrorString = QString("Records in the discharge inflow file <= 1: %1").arg(name);
        throw 1;
    }

    // check if all data rows have the correct nr of columns
//    bool check = false;
//    int err = 0;
//    for (int i = count+2; i < QRecs.size(); i++) {
//        SL = QRecs[i].split(QRegExp("\\s+"));
//        int nrrecords = SL.size();
//        if (nrrecords != count) {
//            err = i;
//            check = true;
//            break;
//        }
//    }
//    if (check) {
//        ErrorString = QString("Nr columns is not equal to nr of records %1 in row in file: %2").arg(err).arg(name);
//        throw 1;
//    }


    for (int i = 2; i < nrStations+2; i++) {
            locationnnrsrec << QRecs[i].toInt(&ok, 10);
    }

    //
//    if(locationnnrsmap != locationnnrsrec) {
//        ErrorString = QString("Discharge input location numbers in timeseries and map are not the same.").arg(err).arg(name);
//        throw 1;
//    }

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




//---------------------------------------------------------------------------

*/
