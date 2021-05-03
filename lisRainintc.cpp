
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

/*!
  \file lisRainintc.cpp
  \brief Get rainfall, make a rainfall map and calculate interception
  \brief Get Discharge input at user defined points

functions: \n
- void TWorld::GetDischargeData(QString name) \n
- void TWorld::GetRainfallData(void) \n
- void TWorld::RainfallMap(void) \n
 */

#include <memory>
#include "io.h"
#include "model.h"
#include "operation.h"


//---------------------------------------------------------------------------
/// read rainfall files of different types and put data in RainfallSeries
/// reads also old RUU lisem rain files
/// can read rainfall maps in between intensity values
/// format: first line ends with integer that is nr of data columns excl time
void TWorld::GetRainfallData(QString name)
{
    RAIN_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList rainRecs;
    QStringList SL;
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    int skiprows = 0;
    double time = 0.0;
    bool oldformat = true;

    if (!fi.exists())
    {
        ErrorString = "Rainfall file not found: " + name;
        throw 1;
    }

    nrRainfallseries = 0;

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.count()-1,1);
        if (!S.trimmed().isEmpty())
            rainRecs << S.trimmed();
    }
    fff.close();

    // check first if PCRaster graph format is present: header, number of vars, columns equal vars
    int count = rainRecs[1].toInt(&ok, 10);
    // header
    // second line is only an integer
    if (ok)
    {
        SL = rainRecs[count+2].split(QRegExp("\\s+"));
        if (count == SL.count())
            oldformat = false;
        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as forst column
    }

    if (rainRecs[0].contains("RUU"))
        oldformat = true;

    if (oldformat)
    {
        QStringList SL = rainRecs[0].split(QRegExp("\\s+"));
        // get first line, white space character as split for header

        nrStations = SL[SL.size()-1].toInt(&ok, 10);
        // read nr stations from last value in old style header
        // failure gives 0
        SL = rainRecs[rainRecs.count()-1].split(QRegExp("\\s+"));
        oldformat = (nrStations == SL.count()-1);
    }

    //check if nr stations found equals nr columns-1, 1st column is time
    if (oldformat)
        skiprows = 1;
    else
        skiprows = 3;

    int nrmap = 0;
    nrmap = countUnits(*RainZone);

    if (nrmap > nrStations)
    {
        ErrorString = QString("Number of stations in rainfall file (%1) < nr of rainfall zones in ID map (%2)").arg(nrStations).arg(nrmap);
        throw 1;
    }

    nrSeries = rainRecs.size() - nrStations - skiprows;
    // count rainfall or snowmelt records

    if (nrSeries <= 1)
    {
        ErrorString = "Rainfall records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.intensity.clear();
        rl.isMap = false;
        rl.name = "";
        QString dirname;
        dirname = rainFileDir;

        // split rainfall record row with whitespace
        QStringList SL = rainRecs[r+nrStations+skiprows].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);

        qDebug() << rl.time << SL[0];

        if (r == 0)
            time = rl.time;

        if (r > 0 && rl.time <= time)
        {
            qDebug() << r << time << rl.time;
            ErrorString = QString("Rainfall records at time %1 has unreadable value.").arg(rl.time);
            throw 1;
        }
        else
            time = rl.time;

        // check if record has characters, then filename assumed
        if (SwitchRainfallSatellite)//  SL[1].contains(QRegExp("[A-Za-z]")))
        {
            QFileInfo fi(QDir(dirname), SL[1]);
            // asume second record is name
            if (!fi.exists())
            {
                ErrorString = QString("Rainfall map %1 not found. Rainfall maps must be in the rainfall directory.").arg(SL[1]);
                throw 1;
            }

            rl.name = fi.absoluteFilePath();
            rl.isMap = true;
            // a mapname if the file exists
        }
        else
        {
            // record is a assumed to be a double

            for (int i = 1; i <= nrStations; i++)
            {
                bool ok = false;
                rl.intensity << SL[i].toDouble(&ok);
                if (!ok)
                {
                    ErrorString = QString("Rainfall records at time %1 has unreadable value.").arg(SL[0]);
                    throw 1;
                }
            }
        }
        RainfallSeriesM << rl;
    }

    if(!SwitchETSatellite) {
        nrSeries++;
        rl.time = rl.time+1440;

        for (int i = 1; i < nrStations; i++)
            rl.intensity << 0.0;

        RainfallSeriesM << rl;
    }

    nrRainfallseries = nrSeries;
}
//---------------------------------------------------------------------------
/**
rainfall intensity read is that reported with the current line: example\n
0 0\n
5 2.3   ->from 0 to 5 minutes intensity is 0\n
7.5 4.5 ->from 5 to 7.5 minutes intensity is 2.3\n
etc. */
void TWorld::RainfallMap(void)
{
    double timeminprev = (time-_dt) / 60; //prev time in minutes
    int  rainplace;
    double tt = 3600000.0;

    if (!SwitchRainfall)
        return;

    for (rainplace = 0; rainplace < nrRainfallseries; rainplace++)
        if (timeminprev < RainfallSeriesM[rainplace].time)
            break;

    if (RainfallSeriesM[rainplace].isMap)
    {
        auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(
            RainfallSeriesM[rainplace].name)));

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            if (pcr::isMV(_M->Drc))
            {
                QString sr, sc;
                sr.setNum(r); sc.setNum(c);
                ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+RainfallSeriesM[rainplace].name;
                throw 1;
            }
            else
                Rain->Drc = _M->Drc *_dt/tt;
        }}
    }
    else
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            Rain->Drc = RainfallSeriesM[rainplace].intensity[(int) RainZone->Drc-1]*_dt/tt;
            //TODO: weighted average if dt larger than table dt
        }}
    }

    if (!rainStarted) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if(Rain->Drc > 0) rainStarted = true;
        }}
    }


    if (rainStarted && RainstartTime == -1)
        RainstartTime = time;

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        Rainc->Drc = Rain->Drc * _dx/DX->Drc;
        // correction for slope dx/DX, water spreads out over larger area
        RainCumFlat->Drc += Rain->Drc;
        // cumulative rainfall
        RainCum->Drc += Rainc->Drc;
        // cumulative rainfall corrected for slope, used in interception
        RainNet->Drc = Rainc->Drc;
        // net rainfall in case of interception
    }}
}
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
qDebug() << locationnnrsrec;
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
    //---------------------------------------------------------------------------
void TWorld::GetETData(QString name)
{
    RAIN_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList ETRecs;
    QStringList SL;
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    int skiprows = 0;
    double time = 0.0;

    if (!SwitchIncludeET)
        return;

    if (!fi.exists())
    {
        ErrorString = "ET file not found: " + name;
        throw 1;
    }

    nrETseries = 0;

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    // get all rainfall records without empty lines
    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.count()-1,1);
        if (!S.trimmed().isEmpty())
            ETRecs << S.trimmed();
    }
    fff.close();

    // check first if PCRaster graph format is present: header, number of vars, columns equal vars
    int cols = ETRecs[1].toInt(&ok, 10); // nr of columns specified in file

    // header
    // second line is only an integer
    bool oldformat = true;
    if (ok)
    {
        SL = ETRecs[cols+2].split(QRegExp("\\s+"));  // get first data row and split in cols
        if (cols == SL.count())
            oldformat = false;
        //if the number of columns equals the integer then new format
        nrStations = cols-1;
        // nr stations is cols-2 with day/hour/min as first column
    }
    if (ETRecs[0].contains("RUU"))
        oldformat = true;

    if(oldformat) {
        ErrorString = "old style rainfall/ET records no longer supported";
        throw 1;
    }

    //check if nr stations found equals nr columns-1, 1st column is time
    skiprows = 3;

    int nrmap = 0;
    nrmap = countUnits(*ETZone);

    if (nrmap > nrStations)
    {
        ErrorString = QString("Number of stations in ET file (%1) < nr of rainfall zones in ETID map (%2)").arg(nrStations).arg(nrmap);
        throw 1;
    }

    nrSeries = ETRecs.size() - nrStations - skiprows;
    // count nr data records

    if (nrSeries <= 1)
    {
        ErrorString = "ET data records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.intensity.clear();
        rl.isMap = false;
        rl.name = "";
        QString dirname;
        dirname = ETFileDir;

        // split rainfall record row with whitespace
        QStringList SL = ETRecs[r+nrStations+skiprows].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);

        if (r == 0)
            time = rl.time;

        if (r > 0 && rl.time <= time)
        {
            ErrorString = QString("ET records at time %1 has unreadable value.").arg(rl.time);
            throw 1;
        }
        else
            time = rl.time;

        // check if record has characters, then filename assumed

        if (SwitchETSatellite)
        {
            QFileInfo fi(QDir(dirname), SL[1]);
            // asume second record is name
            if (!fi.exists())
            {
                ErrorString = QString("ET maps %1 not found. ET maps must be in the ET directory.").arg(SL[1]);
                throw 1;
            }

            rl.name = fi.absoluteFilePath();
            rl.isMap = true;
            // a mapname if the file exists
        } else {
            // record is assumed to be a double

            for (int i = 1; i <= nrStations; i++)
            {
                bool ok = false;
                rl.intensity << SL[i].toDouble(&ok);
                if (!ok)
                {
                    ErrorString = QString("ET records at time %1 has unreadable value.").arg(SL[0]);
                    throw 1;
                }
            }
        }
        ETSeriesM << rl;
    }

    if(!SwitchETSatellite) {
        nrSeries++;
        rl.time = rl.time+1440;
        for (int i = 1; i < nrStations; i++)
            rl.intensity << 0.0;
        ETSeriesM << rl;
    }

    nrETseries = nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetSnowmeltData(QString name)
{
    RAIN_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList rainRecs;
    QStringList SL;
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    int skiprows = 0;
    double time = 0.0;
    bool oldformat = true;

    if (!SwitchSnowmelt)
        return;

    if (!fi.exists())
    {
        ErrorString = "SnowMelt file not found: " + name;
        throw 1;
    }

    nrSnowmeltseries = 0;

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.count()-1,1); // readLine also reads \n as a character on an empty line!
        if (!S.trimmed().isEmpty())
            rainRecs << S.trimmed();
    }
    fff.close();

    // check first if PCRaster graph format is present: header, number of vars, columns equal vars
    int count = rainRecs[1].toInt(&ok, 10);
    // header
    // second line is only an integer
    if (ok)
    {
        SL = rainRecs[count+2].split(QRegExp("\\s+"));
        if (count == SL.count())
            oldformat = false;
        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as forst column
    }

    if (rainRecs[0].contains("RUU"))
        oldformat = true;

    if (oldformat)
    {
        QStringList SL = rainRecs[0].split(QRegExp("\\s+"));
        // get first line, white space character as split for header

        nrStations = SL[SL.size()-1].toInt(&ok, 10);
        // read nr stations from last value in old style header
        // failure gives 0
        SL = rainRecs[rainRecs.count()-1].split(QRegExp("\\s+"));
        oldformat = (nrStations == SL.count()-1);
    }

    //check if nr stations found equals nr columns-1, 1st column is time
    if (oldformat)
        skiprows = 1;
    else
        skiprows = 3;

    int nrmap = 0;
    nrmap = countUnits(*SnowmeltZone);
    if (nrmap > nrStations)
    {
        ErrorString = QString("Number of stations in Snowmelt file (%1) < nr of rainfall zones in SNOWID map (%2)").arg(nrStations).arg(nrmap);
        throw 1;
    }
    nrSeries = rainRecs.size() - nrStations - skiprows;
    // count rainfall or snowmelt records

    if (nrSeries <= 1)
    {
        ErrorString = "Snowmelt records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize rainfall record structure
        rl.time = 0;
        rl.intensity.clear();
        rl.isMap = false;
        rl.name = "";
        QString dirname;
        dirname = snowmeltFileDir;

        // split rainfall record row with whitespace
        QStringList SL = rainRecs[r+nrStations+skiprows].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);

        if (r == 0)
            time = rl.time;


        if (r > 0 && rl.time <= time)
        {
            ErrorString = QString("Snowmelt records at time %1 has unreadable value.").arg(rl.time);
            throw 1;
        }
        else
            time = rl.time;

        // check if record has characters, then filename assumed
        if (SwitchRainfallSatellite)//  SL[1].contains(QRegExp("[A-Za-z]")))
        {
            QFileInfo fi(QDir(dirname), SL[1]);
            // asume second record is name
            if (!fi.exists())
            {
                ErrorString = QString("Snowmelt map %1 not found. Snowmelt maps must be in the rainfall directory.").arg(SL[1]);
                throw 1;
            }

            rl.name = fi.absoluteFilePath();
            rl.isMap = true;
            // a mapname if the file exists
        }
        else
        {
            // record is a assumed to be a double
            for (int i = 1; i <= nrStations; i++)
            {
                bool ok = false;
                rl.intensity << SL[i].toDouble(&ok);
                if (!ok)
                {
                    ErrorString = QString(" records at time %1 has unreadable value.").arg(SL[0]);
                    throw 1;
                }
            }
        }

        SnowmeltSeriesM << rl;
    }

    if(!SwitchETSatellite) {
        nrSeries++;
        rl.time = rl.time+1440;

        for (int i = 1; i < nrStations; i++)
            rl.intensity << 0.0;

        SnowmeltSeriesM << rl;
    }

    nrSnowmeltseries = nrSeries;
}
//---------------------------------------------------------------------------
double TWorld::getTimefromString(QString sss)
{
    // read date time string and convert to time in minutes

    if (!sss.contains(QRegExp("[-:/]"))) {
        return(sss.toDouble());

    } else {
        // DDD/HH/MM or DDD-HH-MM or DDD:HH:MM
        QStringList DHM = sss.split(QRegExp("[-:/]"));
        bool nohour = (DHM.count() == 2);
        double day = DHM.at(0).toDouble();
        double hour = 0;
        double min = 0;
        if (nohour)
            min = DHM.at(1).toDouble();
        else {
            hour = DHM.at(1).toDouble();
            min = DHM.at(2).toDouble();
        }

        return((day-1)*1440+hour*60+min);
    }
}
//---------------------------------------------------------------------------
void TWorld::ETMap(void)
{
    double timeminprev = (time-_dt) / 60; //prev time in minutes
    int  ETplace;
    double tt = 3600000.0;

    if (!SwitchRainfall)
        return;

    for (ETplace = 0; ETplace < nrETseries; ETplace++)
        if (timeminprev < ETSeriesM[ETplace].time)
            break;

    if (ETSeriesM[ETplace].isMap)
    {
        auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(ETSeriesM[ETplace].name)));

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            if (pcr::isMV(_M->Drc))
            {
                QString sr, sc;
                sr.setNum(r); sc.setNum(c);
                ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+ETSeriesM[ETplace].name;
                throw 1;
            }
            else
                ETp->Drc = _M->Drc *_dt/tt;
        }}
    }
    else
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            ETp->Drc = ETSeriesM[ETplace].intensity[(int) ETZone->Drc-1]*_dt/tt;
        }}
    }

    if (!ETStarted) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if(ETp->Drc > 0) ETStarted = true;
        }}
    }


    if (ETStarted && ETstartTime == -1)
        ETstartTime = time;
}
