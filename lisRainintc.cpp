
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
  \brief Get rainfall, make a rainfall map
  \brief Get ET, make a ETa map

functions: \n
- void TWorld::GetSpatialMeteoData(QString name, int type) \n
- void TWorld::GetRainfallData(void) \n
- void TWorld::GetETData(QString name) \n
- void TWorld::GetRainfallMap(void) \n
- void TWorld::GetETMap(void) \n
- double TWorld::getTimefromString(QString sss) \n
 */

#include <memory>
#include "io.h"
#include "model.h"
#include "operation.h"

// read the text file with list of maps, check the filenames and add to a record with time in minutes
void TWorld::GetSpatialMeteoData(QString name, int type)
{
    METEO_LIST rl;
    QFile fff(name);
    QFileInfo fi(name);
    QString S;
    QStringList rainRecs;
    QStringList SL;
    int skip = 2;
    int nrSeries = 0;
    double time = 0.0;

    if (!fi.exists())
    {
        if (type == 0)
            ErrorString = "Rainfall file not found: " + name;
        if (type == 1)
            ErrorString = "Et file not found: " + name;
        if (type == 2)
            ErrorString = "Snowmelt file not found: " + name;
        throw 1;
    }

    // read all lines in the text file
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
   // int count = rainRecs[1].toInt(&ok, 10);

    // format
    //header
    // 2 (variables)
    // // DDD/HH/MM or DDD-HH-MM or DDD:HH:MM
    // map name

    nrSeries = rainRecs.size() - skip;
    // count records

    if (nrSeries <= 1)
    {
        if (type == 0)
            ErrorString = "Rainfall records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        if (type == 1)
            ErrorString = "ET records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        if (type == 2)
            ErrorString = "Snowmelt records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    // initalize
    RainfallSeriesMaps.clear();
    ETSeriesMaps.clear();
    SnowmeltSeriesMaps.clear();

    currentRainfallrow = 0;
    currentETrow = 0;
    currentSnowmeltrow = 0;


    QString dirname;
    if (type == 0)
        dirname = rainFileDir;
    if (type == 1)
        dirname = ETFileDir;
    if (type == 2)
        dirname = snowmeltFileDir;

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize record structure
        rl.time = 0;
        rl.name = "";

        // split rainfall record row with whitespace
        QStringList SL = rainRecs[r+skip].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);

        qDebug() << rl.time << SL[0];

        if (r == 0)
            time = rl.time;

        //??????????????????? what does this check???
        if (r > 0 && rl.time <= time)
        {
            qDebug() << r << time << rl.time;
            ErrorString = QString("Rainfall records at time %1 has unreadable value.").arg(rl.time);
            throw 1;
        }
        else
            time = rl.time;

        // check if filename exists
        QFileInfo fi(QDir(dirname), SL[1]);
            // asume second record is name
        if (!fi.exists())
        {
            if (type == 0)
                ErrorString = QString("Rainfall map %1 not found. Rainfall maps must be in the rainfall directory.").arg(SL[1]);
            if (type == 1)
                ErrorString = QString("ET map %1 not found. Rainfall maps must be in the rainfall directory.").arg(SL[1]);
            if (type == 2)
                ErrorString = QString("Snowmelt map %1 not found. Rainfall maps must be in the rainfall directory.").arg(SL[1]);
            throw 1;
        }

        rl.name = fi.absoluteFilePath();

        // add the record to the list
        if (type == 0)
            RainfallSeriesMaps << rl;
        if (type == 1)
            ETSeriesMaps << rl;
        if (type == 2)
            SnowmeltSeriesMaps << rl;

    }

    if (type == 0)
        nrRainfallseries = nrSeries;
    if (type == 1)
        nrSnowmeltseries = nrSeries;
    if (type == 2)
        nrETseries = nrSeries;
}

//---------------------------------------------------------------------------
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
    RainfallSeriesM.clear();
    currentRainfallrow = 0;

    // read rainfall text file
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

    // count gauge areas in the ID.map
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

        // split rainfall record row with whitespace
        QStringList SL = rainRecs[r+nrStations+skiprows].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);

        qDebug() << rl.time << SL[0];

        if (r == 0)
            time = rl.time;

        if (r > 0 && rl.time <= time)
        {
            qDebug() <<"e"<< r << time << rl.time;
            ErrorString = QString("Rainfall records at time %1 has unreadable value.").arg(rl.time);
            throw 1;
        }
        else
            time = rl.time;

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

        RainfallSeriesM << rl;
    }

    // sometimes not an increasing timeseries
    for(int i = 1; i < nrSeries; i++){
        if (RainfallSeriesM [i].time <= RainfallSeriesM[i-1].time) {
            ErrorString = QString("Rainfall records time is not increasing at row %1.").arg(i);
            throw 1;
        }
    }

    //add a zero at the end
//    nrSeries++;
//    rl.time = RainfallSeriesM[nrSeries-1].time+1440;

//    for (int i = 1; i < nrStations; i++)
//        rl.intensity << 0.0;

  //  RainfallSeriesM << rl;

    nrRainfallseries = RainfallSeriesM.size();//nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetRainfallMap(void)
{
    double currenttime = (time)/60;
    int  rainplace;
    double tt = 3600000.0;
    bool norain = false;
    bool samerain = false;

    if (!SwitchRainfall)
        return;
    if (SwitchRainfallSatellite)
        return;

    // from time t to t+1 the rain is the rain of t

    // where are we in the series
    int currentrow = 0;
    // if time is outside records then use map with zeros
    if (currenttime < RainfallSeriesM[0].time || currenttime > RainfallSeriesM[nrRainfallseries-1].time)
        norain = true;
 //   qDebug() << norain << currenttime << RainfallSeriesM[0].time << RainfallSeriesM[nrRainfallseries-1].time;

    if (norain) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Rain->Drc = 0;
        }}
    } else {
        // find current record
        for (rainplace = currentRainfallrow; rainplace < nrRainfallseries-1; rainplace++) {
    //qDebug() << currenttime << RainfallSeriesM[rainplace].time << RainfallSeriesM[rainplace+1].time;
            if (currenttime >= RainfallSeriesM[rainplace].time && currenttime < RainfallSeriesM[rainplace+1].time) {
                currentrow = rainplace;                
                break;
            }
        }
        if (currentrow == currentRainfallrow && rainplace > 0)
            samerain = true;
        else
            currentRainfallrow = currentrow;

        // get the next map from file
        if (!samerain) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                Rain->Drc = RainfallSeriesM[rainplace].intensity[(int) RainZone->Drc-1]*_dt/tt;
            }}
        }
    }

    // find start time of rainfall, for flood peak and rain peak
    if (!rainStarted) {
        FOR_ROW_COL_MV {
            if(Rain->Drc > 0) {
                rainStarted = true;
                break;
            }
        }
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
void TWorld::GetRainfallSatMap(void)
{
    double currenttime = (time)/60;
    int  rainplace;
    double tt = 3600000.0;
    bool norain = false;
    bool samerain = false;

    if (!SwitchRainfall)
        return;
    if (!SwitchRainfallSatellite)
        return;

    // from time t to t+1 the rain is the rain of t

    // where are we in the series
    int currentrow = 0;
    // if time is outside records then use map with zeros
    if (currenttime < RainfallSeriesM[0].time || currenttime > RainfallSeriesM[nrRainfallseries-1].time)
        norain = true;

    if (norain) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Rain->Drc = 0;
        }}
    } else {
        // find current record
        for (rainplace = currentRainfallrow; rainplace < nrRainfallseries-1; rainplace++) {
            if (currenttime >= RainfallSeriesM[rainplace].time && currenttime < RainfallSeriesM[rainplace+1].time) {
                currentrow = rainplace;
                break;
            }
        }
        if (currentrow == currentRainfallrow && rainplace > 0)
            samerain = true;
        else
            currentRainfallrow = currentrow;

        // get the next map from file
        if (!samerain) {
            auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(RainfallSeriesMaps[rainplace].name)));

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (pcr::isMV(_M->Drc)) {
                        QString sr, sc;
                        sr.setNum(r); sc.setNum(c);
                        ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+RainfallSeriesMaps[rainplace].name;
                        throw 1;
                } else {
                  Rain->Drc = _M->Drc *_dt/tt;
                }
            }}
        }
    }

    // find start time of rainfall, for flood peak and rain peak
    if (!rainStarted) {
        FOR_ROW_COL_MV {
            if(Rain->Drc > 0) {
                rainStarted = true;
                break;
            }
        }
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

        return(day*1440+hour*60+min);
    }
}
//---------------------------------------------------------------------------

