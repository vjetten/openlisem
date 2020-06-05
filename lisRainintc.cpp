
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
  \file lisRainintc.cpp
  \brief Get rainfall, make a rainfall map and calculate interception

functions: \n
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
void TWorld::GetRainfallDataM(QString name, bool israinfall)
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
    QString errorS = (israinfall ? "Rainfall" : "Snowmelt");
    bool oldformat = true;


    if (!fi.exists())
    {
        ErrorString = errorS + " file not found: " + name;
        throw 1;
    }

    nrRainfallseries = 0;
    nrSnowmeltseries = 0;

    fff.open(QIODevice::ReadOnly | QIODevice::Text);

    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.contains("\n"))
            S.remove(S.count()-1,1);
        // readLine also reads \n as a character on an empty line!
        if (!S.trimmed().isEmpty())
            rainRecs << S.trimmed();
        //qDebug() << S;
    }
    fff.close();
    // get all rainfall records without empty lines
    //rainRecs << QString("1000000 0");
    // DOES NOT WORK for more than one station
    // add a very large number so that the rainfall after the last timestep is 0


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

//TO DO old format check doesn't always work
//    if (!rainRecs[0].contains("TIMESERIE INTENSITY NORMAL"))
//        oldformat = false;
/*
    // if not, check if new PCRaster style rainfall rec
    if (!oldformat || !ok)
    {
        if (rainRecs[1].count() <= 2)   //VJ 140624 <= 2 in case there is a eol character that is missed, to interpret tss files
        {
            SL = rainRecs[1].split(QRegExp("\\s+"));
            nrStations = SL[0].toInt(&ok, 10);
            // new style gives one nrstations more (includes time column) old style only counts rain stations
            skiprows = 2;
        }
        if (SL.count() != 1 || !ok)
        {
            ErrorString = "Cannot read/interpret nr " + errorS + " stations in header file";
            throw 1;
        }
    }
*/
    if (oldformat)
        skiprows = 1;
    else
        skiprows = 3;

    int nrmap = 0;
    if (israinfall)
        nrmap = countUnits(*RainZone);
    else
        nrmap = countUnits(*SnowmeltZone);

//    if (nrmap != nrStations)
//    {
//        ErrorString = QString("Number of stations in rainfall file (%1) != nr of rainfall zones in ID map (%2)").arg(nrStations).arg(nrmap);
//        throw 1;
//    }
    nrSeries = rainRecs.size() - nrStations - skiprows;
    // count rainfall or snowmelt records

    if (nrSeries <= 1)
    {
        ErrorString = errorS + " records <= 1, must at least have one interval with 2 rows: a begin and end time.";
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
        if (israinfall)
            dirname = rainFileDir;
        else
            dirname = snowmeltFileDir;

        QStringList SL = rainRecs[r+nrStations+skiprows].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

          // split rainfall record row with whitespace
        rl.time = SL[0].toDouble();
        // time in min
   //     qDebug() << SL << rl.time;
        if (r == 0)
            time = rl.time;

        if (r > 0 && rl.time <= time)
        {
            ErrorString = errorS + QString(" records at time %1 has unreadable value.").arg(rl.time);
            throw 1;
        }
        else
            time = rl.time;

        // check if record has characters, then filename assumed

        if (SL[1].contains(QRegExp("[A-Za-z]")))
        {
            QFileInfo fi(QDir(dirname), SL[1]);
            // asume second record is name
            if (!fi.exists())
            {
                ErrorString = errorS + QString(" map %1 not found. Rainfall maps must be in the rainfall directory.").arg(SL[1]);
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
                    ErrorString = errorS + QString(" records at time %1 has unreadable value.").arg(SL[0]);
                    throw 1;
                }
            }
        }
        if (israinfall)
            RainfallSeriesM << rl;
        else
            SnowmeltSeriesM << rl;
    }

    nrSeries++;
    rl.time = 1440*365;
    for (int i = 1; i < nrStations; i++)
        rl.intensity << 0.0;

    if (israinfall)
        RainfallSeriesM << rl;
    else
        SnowmeltSeriesM << rl;

    if (israinfall)
        nrRainfallseries = nrSeries;
    else
        nrSnowmeltseries = nrSeries;

//      for (int i = 0; i < RainfallSeriesM.count(); i++)
//         qDebug() << "rain" << RainfallSeriesM[i].time <<RainfallSeriesM[i].intensity << RainfallSeriesM[i].name;

//      for (int i = 0; i < SnowmeltSeriesM.count(); i++)
//         qDebug() << " snow" << SnowmeltSeriesM[i].time <<SnowmeltSeriesM[i].intensity << SnowmeltSeriesM[i].name;
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

        FOR_ROW_COL_MV
                if (pcr::isMV(_M->Drc))
        {
            QString sr, sc;
            sr.setNum(r); sc.setNum(c);
            ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+RainfallSeriesM[rainplace].name;
            throw 1;
        }
        FOR_ROW_COL_MV
        {
            Rain->Drc = _M->Drc *_dt/tt;
            if (!rainStarted && Rain->Drc  > 0)
                rainStarted = true;
        }
    }
    else
    {
        FOR_ROW_COL_MV
        {
            if (RainZone->Drc-1 <  RainfallSeriesM[rainplace].intensity.count())
                Rain->Drc = RainfallSeriesM[rainplace].intensity[(int) RainZone->Drc-1]*_dt/tt;
            else
            {
                ErrorString = QString("No rainfall data for ID map zone %1").arg(RainZone->Drc);
                throw 1;
            }
            // Rain in m per timestep from mm/h, rtecord nr corresponds map ID value -1
            //TODO: weighted average if dt larger than table dt
            if (!rainStarted && Rain->Drc  > 0)
                rainStarted = true;
        }
    }
    if (rainStarted && RainstartTime == -1)
    {
        RainstartTime = time;
    }


    FOR_ROW_COL_MV
    {
        Rainc->Drc = Rain->Drc * _dx/DX->Drc;
        // correction for slope dx/DX, water spreads out over larger area
        RainCumFlat->Drc += Rain->Drc;
        // cumulative rainfall
        RainCum->Drc += Rainc->Drc;
        // cumulative rainfall corrected for slope, used in interception
        RainNet->Drc = Rainc->Drc;
        // net rainfall in case of interception

        if (Rain->Drc == 0)
            noRain->Drc += _dt;
        else
            noRain->Drc = 0;

//        if (noRain->Drc > 3.0*3600.0)
  //          RainCum->Drc = 0;
        // if dry spell for more than 23 hours

    }
}
//---------------------------------------------------------------------------

