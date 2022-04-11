
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
    int skip = 4;
    int nrSeries = 0;

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
     //   qDebug() << S;
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

    QString dirname;
    if (type == 0) {
        RainfallSeriesMaps.clear();
        dirname = rainSatFileDir;
        currentRainfallrow = 0;
    }
    if (type == 1) {
        ETSeriesMaps.clear();
        dirname = ETSatFileDir;
        currentETrow = 0;
    }
    if (type == 2) {
        SnowmeltSeriesMaps.clear();
        dirname = snowmeltFileDir;
        currentSnowmeltrow = 0;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize record structure
        rl.time = 0;
        rl.name = "";
        rl.calib = 1.0;

        // split rainfall record row with whitespace
        QStringList SL = rainRecs[r+skip].split(QRegExp("\\s+"), Qt::SkipEmptyParts);
//qDebug() << SL;
        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);

        // check if filename exists
        QFileInfo fi(QDir(dirname), SL[1]);
            // asume second record is name
        if (SL.count() > 2) {
            bool ok;
            double v = SL[2].toDouble(&ok);
            if (ok)
                rl.calib = v;
        }

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
      //  qDebug() << rl.time << rl.name;

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
        nrETseries = nrSeries;
    if (type == 2)
        nrSnowmeltseries = nrSeries;
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
    RainfallSeries.clear();
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

        //qDebug() << rl.time << SL[0];

        if (r == 0)
            time = rl.time;

        if (r > 0 && rl.time <= time)
        {
            //qDebug() <<"e"<< r << time << rl.time;
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

        RainfallSeries << rl;
    }

    // sometimes not an increasing timeseries
    for(int i = 1; i < nrSeries; i++){
        if (RainfallSeries[i].time <= RainfallSeries[i-1].time) {
            ErrorString = QString("Rainfall records time is not increasing at row %1.").arg(i);
            throw 1;
        }
    }
//    for(int i = 0; i < nrSeries; i++)
//        qDebug() << RainfallSeries[i].time << RainfallSeries[i].intensity.at(0);

    nrRainfallseries = RainfallSeries.size();//nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetRainfallMapfromStations(void)
{
    double currenttime = (time)/60;
    double tt = _dt/3600000.0;
    bool samerain = false;

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
    int currentrow = 0;// rainplace;
    // find current record
//    while (currenttime >= RainfallSeries[rainplace].time
//        && currenttime < RainfallSeries[rainplace+1].time)
//    {
//        currentrow = rainplace;
//        rainplace++;
//    }
//qDebug() << time/86400 << currenttime << rainplace << currentrow  << RainfallSeries[currentrow].time << RainfallSeries[currentrow].intensity[0];
//    if (currentrow == currentRainfallrow && currentrow > 0)
//        samerain = true;

    for (int j = 0; j < RainfallSeries.count(); j++) {
        if (currenttime >= RainfallSeries[j].time && currenttime < RainfallSeries[j+1].time) {
            currentrow = j;
            break;
        }
    }
    if (currentrow == currentRainfallrow && currentrow > 0)
        samerain = true;


    // get the next map from file
    if (!samerain) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Rain->Drc = RainfallSeries[currentrow].intensity[(int) RainZone->Drc-1]*tt;
            if (Rain->Drc > 0)
                rainStarted = true;
        }}
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Rainc->Drc = Rain->Drc * _dx/DX->Drc;
        // correction for slope dx/DX, water spreads out over larger area
        RainCumFlat->Drc += Rain->Drc;
        // cumulative rainfall
        RainCum->Drc += Rainc->Drc;
        // cumulative rainfall corrected for slope, used in interception
        RainNet->Drc = Rainc->Drc;
        // net rainfall in case of interception
    }}

    currentRainfallrow = currentrow;

    // find start time of rainfall, for flood peak and rain peak
    if (rainStarted && RainstartTime == -1)
        RainstartTime = time;

}
//---------------------------------------------------------------------------
void TWorld::GetRainfallMap(void)
{
    double currenttime = (time)/60;
    double tt = _dt/3600000.0 * PBiasCorrection; // mm/h to m -> mm/h = mm X/3600*_dt -> X*0.0001
    bool samerain = false;

    // from time t to t+1 the rain is the rain of t

    // where are we in the series
    // if time is outside records then use map with zeros
    if (currenttime < RainfallSeriesMaps[0].time || currenttime > RainfallSeriesMaps[nrRainfallseries-1].time) {
        DEBUG("run time outside rainfall records");
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Rain->Drc = 0;
            RainNet->Drc = 0;
            Rainc->Drc = 0;
        }}
        return;
    }

    int currentrow = rainplace;
    // find current record
//    while (currenttime >= RainfallSeriesMaps[rainplace].time
//        && currenttime < RainfallSeriesMaps[rainplace+1].time)
//    {
//        currentrow = rainplace;
//        rainplace++;
//    }




    for (int j = 0; j < RainfallSeriesMaps.count(); j++) {
        if (currenttime >= RainfallSeriesMaps[j].time && currenttime < RainfallSeriesMaps[j+1].time) {
            currentrow = j;
            break;
        }
    }

    if (currentrow == currentRainfallrow && currentrow > 0)
        samerain = true;
    //qDebug() << currentrow << currenttime << currentRainfallrow << samerain;

    // get the next map from file
    if (!samerain) {
      //  qDebug() << currentrow << RainfallSeriesMaps[currentrow].name;
        // read a map
        //auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(RainfallSeriesMaps[currentrow].name)));
        cTMap *_M = new cTMap(readRaster(RainfallSeriesMaps[currentrow].name));
        double calibration = RainfallSeriesMaps[currentrow].calib;
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            double rain_ = 0;
            if (pcr::isMV(_M->Drc)) {
                QString sr, sc;
                sr.setNum(r); sc.setNum(c);
                ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+RainfallSeriesMaps[rainplace].name;
            } else
                rain_ = _M->Drc * tt; // * RainfallSeriesMaps[currentrow].calib;

            if (rain_ < 0)
                rain_ = 0;
            if (rain_ > 0)
                rainStarted = true;
            Rain->Drc = rain_ * calibration;
        }}
        delete _M;
    } //samerain

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Rainc->Drc = Rain->Drc * _dx/DX->Drc;
        // correction for slope dx/DX, water spreads out over larger area
        RainCumFlat->Drc += Rain->Drc;
        // cumulative rainfall
        RainCum->Drc += Rainc->Drc;
        // cumulative rainfall corrected for slope, used in interception
        RainNet->Drc = Rainc->Drc;
        // net rainfall in case of interception
    }}

    currentRainfallrow = currentrow;

    if (rainStarted && RainstartTime == -1)
        RainstartTime = time;

}

//---------------------------------------------------------------------------
    // not used
double TWorld::getmaxRainfall()
{
    double maxv = 0;
    double avg = 0;
    double tt = _dt/3600000.0;
    if (SwitchRainfallSatellite) {
        for (int i = 0; i < nrRainfallseries-1; i++) {
            auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(RainfallSeriesMaps[i].name)));
            avg = MapTotal(*_M)/nrCells;
            maxv = std::max(avg, maxv);
        }
    } else {
        avg = 0;
        for (int i = 0; i < nrRainfallseries-1; i++) {
            for (int j = 0; j < RainfallSeries[i].intensity.size(); j++)
                avg = avg + RainfallSeries[i].intensity[j]*tt;
        }
        maxv = std::max(maxv, avg);
    }
    return (maxv);
}
//---------------------------------------------------------------------------
double TWorld::getTimefromString(QString sss)
{
    // read date time string and convert to time in minutes
    double day = 0;
    double hour = 0;
    double min = 0;
    if (SwitchEventbased) {
        min = sss.toDouble();
        return(min);
    }

    if (!sss.contains(QRegExp("[-:/]"))) {
        day = sss.toDouble();
    } else {
        // DDD/HH/MM or DDD-HH-MM or DDD:HH:MM
        QStringList DHM = sss.split(QRegExp("[-:/]"));

        if (DHM.count() == 2) {
            day = DHM.at(0).toDouble();
            min = DHM.at(1).toDouble();
        } else {
            day = DHM.at(0).toDouble();
            hour = DHM.at(1).toDouble();
            min = DHM.at(2).toDouble();
        }
    }
    return(day*1440+hour*60+min);
}
//---------------------------------------------------------------------------

