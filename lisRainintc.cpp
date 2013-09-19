
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

#include "model.h"


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
    bool ok;
    int nrStations = 0;
    int nrSeries = 0;
    int skiprows = 0;
    double time;
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

    //   qDebug() << rainRecs.count();

    QStringList SL = rainRecs[0].split(QRegExp("\\s+"));
    // get first line, white space character as split for header

    nrStations = SL[SL.size()-1].toInt(&ok, 10);
    // read nr stations from last value in old style header
    // failure gives 0
    nrStations += 1;
    SL = rainRecs[rainRecs.count()-1].split(QRegExp("\\s+"));
    oldformat = (nrStations == SL.count());
    //check if nr stations found equals nr columns-1, 1st column is time

    // if not, check if new PCRaster style rainfall rec
    if (!oldformat || !ok)
    {
        if (rainRecs[1].count() == 1)
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

    int nr = RainZone->countUnits();
    if (nr != nrStations-1)
    {
        ErrorString = QString("Number of stations in rainfall file (%1) does not match number of units in ID map (%2)").arg(nrStations-1).arg(nr);
        throw 1;
    }

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

        QStringList SL = rainRecs[r+nrStations+skiprows].split(QRegExp("\\s+"), QString::SkipEmptyParts);
        // split rainfall record row with whitespace
        rl.time = SL[0].toDouble();
        // time in min

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

            for (int i = 1; i < nrStations; i++)
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

    if (israinfall)
        nrRainfallseries = nrSeries;
    else
        nrSnowmeltseries = nrSeries;

    //  for (int i = 0; i < RainfallSeriesM.count(); i++)
    //     qDebug() << RainfallSeriesM[i].time <<RainfallSeriesM[i].intensity << RainfallSeriesM[i].name;

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
    //double timemin = time / 60;  //time in minutes
    double timeminprev = (time-_dt) / 60; //prev time in minutes
    int  place;
    double tt = 3600000.0;

    if (!SwitchRainfall)
        return;

    for (place = 0; place < nrRainfallseries; place++)
        if (timeminprev < RainfallSeriesM[place].time)
            break;

    if (RainfallSeriesM[place].isMap)
    {
        TMMap *_M = new TMMap();
        _M->PathName = RainfallSeriesM[place].name;
        bool res = _M->LoadFromFile();
        if (!res)
        {
            ErrorString = "Cannot find map " +_M->PathName;
            throw 1;
        }

        for (int r = 0; r < _nrRows; r++)
            for (int c = 0; c < _nrCols; c++)
                if (!IS_MV_REAL8(&LDD->Drc) && IS_MV_REAL8(&_M->Drc))
                {
                    QString sr, sc;
                    sr.setNum(r); sc.setNum(c);
                    ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+RainfallSeriesM[place].name;
                    throw 1;
                }
        FOR_ROW_COL_MV
        {
            Rain->Drc = _M->Drc *_dt/tt;
            Rainc->Drc = Rain->Drc * _dx/DX->Drc;
            RainCum->Drc += Rainc->Drc;
            RainNet->Drc = Rainc->Drc;
        }

        _M->KillMap();
    }
    else
    {
        FOR_ROW_COL_MV
        {
            Rain->Drc = RainfallSeriesM[place].intensity[(int) RainZone->Drc-1]*_dt/tt;
            // Rain in m per timestep from mm/h, rtecord nr corresponds map nID value -1
            Rainc->Drc = Rain->Drc * _dx/DX->Drc;
            // correction for slope dx/DX, water spreads out over larger area

            //TODO: weighted average if dt larger than table dt

            RainCum->Drc += Rainc->Drc;
            // cumulative rainfall corrected for slope, used in interception
            RainNet->Drc = Rainc->Drc;
        }
    }

}
//---------------------------------------------------------------------------

/// read rainfall files of different types and put data in RainfallSeries
/// reads also old RUU lisem rain files
/// reads rainfall maps
//OBSOLETE
/*
void TWorld::GetRainfallData(void)
{
   QFile fff(rainFileName);
   QFileInfo fi(rainFileName);
   QString S;
   bool ok;
   int j = 0;

   if (!fi.exists())
   {
      ErrorString = "Rainfall file not found: " + rainFileName;
      throw 1;
   }

   nrstations = 0;
   nrrainfallseries = 0;

   fff.open(QIODevice::ReadOnly | QIODevice::Text);
   S = fff.readLine();
   // read the header
   while (S.isEmpty())
      S = fff.readLine();
   // skip empty lines

   QStringList SL = S.split(QRegExp("\\s+")); //<== white space character as split
   nrstations = SL[SL.size()-2].toInt(&ok, 10);
   if (!ok)
   {
      ErrorString = "Cannot read nr rainfall stations in header rainfall file";
      throw 1;
   }

   for(int r=0; r < nrstations+1; r++)
      S = fff.readLine();
   // read column headers

   while (!fff.atEnd())
   {
      S = fff.readLine();
      qDebug() << S << nrrainfallseries;
      if (!S.trimmed().isEmpty())
         nrrainfallseries++;
      qDebug() << S << nrrainfallseries;
   }
   // count rainfall records, skip empty lines


   if (nrrainfallseries <= 1)
   {
      ErrorString = "rainfall records <= 1, must at least have one interval with a begin and end time.";
      throw 1;
   }


   nrrainfallseries++;
   RainfallSeries = new double*[nrrainfallseries];
   for(int r=0; r < nrrainfallseries; r++)
      RainfallSeries[r] = new double[nrstations+1];
   // make structure to contain rainfall
   //RainfallSeries is matrix with rows is data and 1st column is time, other columns are stations

   fff.close();
   // close file and start again


   fff.open(QIODevice::ReadOnly | QIODevice::Text);
   S = fff.readLine();
   for (int i = 0; i < nrstations; i++)
      S = fff.readLine();
   // read header data
   while (!fff.atEnd())
   {
      S = fff.readLine();
      if (S.trimmed().isEmpty()) continue;

      QStringList SL = S.split(QRegExp("\\s+"), QString::SkipEmptyParts);
      if (SL.size()-1 < nrstations)
      {
         QString ss;
         ErrorString = "Rainfall: Nr stations specified in header = " + ss.setNum(nrstations);
         ErrorString += ", nr columns available = "+ ss.setNum(SL.size()-1);
         throw 1;
      }


      RainfallSeries[j][0] = SL[0].toDouble();
      // time in min
      for (int i = 1; i < nrstations+1; i++)
         RainfallSeries[j][i] = SL[i].toDouble();
      // rainfall intensities
      j++;
   }


   RainfallSeries[nrrainfallseries-1][0] = 1e20;
   for (int i = 1; i < nrstations+1; i++)
      RainfallSeries[nrrainfallseries-1][i] = 0;
   // end series with 0 value and extreme time
   fff.close();
}
*/
//---------------------------------------------------------------------------

