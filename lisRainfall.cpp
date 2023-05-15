
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
// get station data for ID interpolation or tiesen polyhon map
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
        if (S.contains("\r\n"))
            S.remove(S.count()-2,2);
        if (S.contains("\n"))
            S.remove(S.count()-1,1);

        if (!S.trimmed().isEmpty())
            rainRecs << S.trimmed();
    }
    fff.close();

    oldformat = (rainRecs[0].contains("RUU"));
    // original very old format
    if (oldformat) {
        ErrorString = "The old RUU rainfall file format is not longer supported.";
        throw 1;
    }

    // check first if PCRaster graph format is present: header, number of vars, columns equal vars
    int count = rainRecs[1].toInt(&ok, 10); // nr of cols in file
    // header
    // second line is only an integer
    if (ok)
    {
        SL = rainRecs[count+2].split(QRegExp("\\s+"));
        // check nr of columns in file
        if (count != SL.count()) {
            ErrorString = "Rainfall file error: The nr of columns in the rainfall file does not equal the number on the second row.";
            throw 1;
        }

        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as first column
    }

    // get station numbers from header, or fill in 1,2 ... n
    stationID.clear();
    for (int i = 0; i < nrStations; i++) {
        SL = rainRecs[i+3].split(QRegExp("\\s+"));
        int tmp = SL.last().toInt(&ok, 10);
        if (ok)
            stationID << tmp;
        else
            stationID << i+1;

    }
   // qDebug() << "stations" << stationID;

    if (stationID.count() == 1) {
        SwitchIDinterpolation = false;
        DEBUG("Only one rain station found, Inverse distance interpolation not done.");
    }

    if (SwitchIDinterpolation) {

        IDIpointsRC.clear();

        if (SwitchUseIDmap) {
            // read points in the IDgauge.map and check if a corresponding number exists in the rainfall file
            // allow points outside MV mask of LDD
            for(int r = 0; r < _nrRows; r++) {
                for (int c = 0; c < _nrCols; c++) {
                    if (!pcr::isMV(IDRainPoints->Drc) && IDRainPoints->Drc > 0) {
                        IDI_POINT p;
                        p.r = r;
                        p.c = c;
                        p.nr = IDRainPoints->Drc;
                        p.V = 0;
                        IDIpointsRC << p;
                    }
                }
            }

            // check if station nrs correspond to map nrs
            for (int i = 0; i < stationID.count(); i++) {
                bool found = false;

                for (int j = 0; j < IDIpointsRC.count(); j++) {
                    if (IDIpointsRC.at(j).nr == stationID.at(i))
                        found = true;
                }
                if (!found) {
                    ErrorString = "Gauge ID number(s) in IDgauge.map not present in the rainfall input file";
                    throw 1;
                }
            }
        } else {

            for (int i = 0; i < stationID.count(); i++) {
                SL = rainRecs[i+3].split(QRegExp("\\s+"));
                if (SL.count() < 3)
                    break;
                IDI_POINT p;
                p.r = SL[0].toInt();
                p.c = SL[1].toInt();
                p.nr = SL[2].toInt();
                p.V = 0;
                qDebug() << p.r << p.c << p.V;
                IDIpointsRC << p;
            }
        }

    }  else {
        // count gauge areas in the ID.map

        QList <int> tmp;
        tmp = countUnits(*RainZone);
        int nrmap = tmp.count();

        if (nrmap > nrStations)
        {
            ErrorString = QString("Number of stations in rainfall file (%1) < nr of rainfall zones in ID map (%2)").arg(nrStations).arg(nrmap);
            throw 1;
        }
    }

    nrSeries = rainRecs.size() - nrStations - 3;
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
        rl.stationnr.clear();
        int r_ = r+nrStations+3;

        // split rainfall record row with whitespace
        QStringList SL = rainRecs[r_].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);
        time = rl.time;

        // check if time is increasing with next row
        if (r+1 < nrSeries) {
            QStringList SL1 = rainRecs[r_+1].split(QRegExp("\\s+"), Qt::SkipEmptyParts);
            int time1 = getTimefromString(SL1[0]);
            if (time1 < time) {
                ErrorString = QString("Time in rainfall records is not increasing from row %1 to %2. Check your file!").arg(r_).arg(r_+1);
                throw 1;
            }
        }

        // rainfall values in this row
        for (int i = 1; i <= nrStations; i++)
        {
            bool ok = false;

            rl.intensity << SL[i].toDouble(&ok);
            if (!ok)
            {
                ErrorString = QString("Rainfall records at time %1 has unreadable value: %2.").arg(SL[0]).arg(SL[i]);
                throw 1;
            }
            rl.stationnr << stationID.at(i-1);
        }

        RainfallSeries << rl;
    }

    // sometimes not an increasing timeseries
//    for(int i = 1; i < nrSeries; i++){
//        if (RainfallSeries[i].time <= RainfallSeries[i-1].time) {
//            ErrorString = QString("Rainfall records time is not increasing at row %1.").arg(i);
//            throw 1;
//        }
//    }

    nrRainfallseries = RainfallSeries.size();//nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetRainfallMapfromStations(void)
{
    double currenttime = (time)/60;
    double tt = _dt/3600000.0* PBiasCorrection;
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
        if (SwitchIDinterpolation) {
//            for (int j = 0; j < IDIpointsRC.size(); j++) {
//                IDI_POINT p;
//                p.r = IDIpointsRC.at(j).c;
//                p.c = IDIpointsRC.at(j).r;
//                p.nr = IDIpointsRC.at(j).nr;
//                p.V = RainfallSeries[currentrow].intensity[j]*tt;
//                IDIpointsRC.replace(j,p);
//               // qDebug() << j << IDIpointsRC.at(j).nr << IDIpointsRC.at(j).V;
//            }

            bool found = false;
            for (int j = 0; j < IDIpointsRC.size(); j++) {
                // select the right station nr
                for (int k = 0; k < IDIpointsRC.size(); k++) {
                //    qDebug() << j << k << IDIpointsRC.at(j).nr << RainfallSeries[currentrow].stationnr.at(k);
                    if (IDIpointsRC.at(j).nr == RainfallSeries[currentrow].stationnr.at(k)) {
                        IDI_POINT p;
                        p.r = IDIpointsRC.at(j).c;
                        p.c = IDIpointsRC.at(j).r;
                        p.nr = IDIpointsRC.at(j).nr;
                        p.V = RainfallSeries[currentrow].intensity[k]*tt;
                        IDIpointsRC.replace(j,p);
                        found = true;
                    }
                }
            }
            if (!found) {
                ErrorString = "IDI value not found";
                throw 1;
            }

            IDInterpolation();

        } else {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                double value = -1;
                for (int k = 0; k < stationID.size(); k++) {
                    if ((int) RainZone->Drc == RainfallSeries[currentrow].stationnr.at(k))
                        value = RainfallSeries[currentrow].intensity[k]*tt;
                }
                Rain->Drc = value;

                if (Rain->Drc > 0)
                    rainStarted = true;
            }}
        }
    }
    //qDebug() <<  MapTotal(*Rain)<< currentrow << currentRainfallrow << currenttime << RainfallSeries[currentrow].time;
    //report(*Rain,"rain");
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
    for (int j = 0; j < RainfallSeriesMaps.count(); j++) {
        if (currenttime >= RainfallSeriesMaps[j].time && currenttime < RainfallSeriesMaps[j+1].time) {
            currentrow = j;
            break;
        }
    }

    if (currentrow == currentRainfallrow && currentrow > 0)
        samerain = true;
  //  qDebug() << currentrow << currenttime << currentRainfallrow << samerain;

    SwitchdoRrainAverage = false;
    // get the next map from file
    if (!samerain) {
        // create an empty map and read the file
        auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(RainfallSeriesMaps[currentrow].name)));
        //  cTMap *_M = new cTMap(readRaster(RainfallSeriesMaps[currentrow].name));
        double calibration = RainfallSeriesMaps[currentrow].calib;

        if (_M->nrCols() != _nrCols || _M->nrRows() != _nrRows) {
            ErrorString = "Nr of rows or Cols in the rainfall map does not match the database";
            throw 1;
        }

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            double rain_ = 0;
            //tma->Drc = 0;

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
            if (!SwitchdoRrainAverage)
                Rain->Drc= rain_ * calibration;
            else
                tma->Drc = rain_ * calibration;
        }}

        if (SwitchdoRrainAverage) {
            double avg = mapAverage(*tma);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                Rain->Drc = avg;
            }}
        }

      //  delete _M;
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
void TWorld::IDInterpolation()
{
    #pragma omp parallel for collapse(2) num_threads(userCores)
    for(int r = 0; r < _nrRows; r++) {
        for (int c = 0; c < _nrCols; c++) {
            double w_total = 0.0;
            double val_total = 0.0;

            for(int i = 0; i < IDIpointsRC.size(); i++)
            {
                int rr = IDIpointsRC.at(i).r;
                int rc = IDIpointsRC.at(i).c;

                if (rr >= 0 && rr < _nrRows && rc >= 0 && rc < _nrCols) {
                    //calc distance factor
                    double dx = (r-rr) * _dx;
                    double dy = (c-rc) * _dx;

                    if(dx == 0 && dy == 0) {
                        val_total = IDIpointsRC.at(i).V;
                        w_total = 1.0;
                    } else {

                        double distancew = qSqrt(dx*dx + dy*dy);
                        if (rainIDIfactor > 1)
                            distancew = std::pow(distancew,rainIDIfactor);

                        w_total += 1.0/distancew;
                        val_total += IDIpointsRC.at(i).V * 1.0/distancew;

                    }
                }
            }

            if (!pcr::isMV(Rain->Drc)) {
                //if(w_total > 0.0)
                    Rain->Drc = val_total/w_total;
                //else
                  //  Rain->Drc = 0.0;
            }
        }
    }
}

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
        ErrorString = "Discharge file not found: " + name;
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
   // qDebug() << "stations" << stationID;
    // count gauge areas in the ID.map

//    QList <int> tmp;
//    tmp = countUnits(*RainZone);
//    int nrmap = tmp.count();

//    if (nrmap > nrStations)
//    {
//        ErrorString = QString("Number of stations in rainfall file (%1) < nr of rainfall zones in ID map (%2)").arg(nrStations).arg(nrmap);
//        throw 1;
//    }


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
        QStringList SL = rainRecs[r_].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

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
