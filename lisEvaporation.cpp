
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
    bool oldformat = true;

    if (!fi.exists())
    {
        ErrorString = "ET file not found: " + name;
        throw 1;
    }

    nrETseries = 0;

    // read ET text file
    fff.open(QIODevice::ReadOnly | QIODevice::Text);
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
    int count = ETRecs[1].toInt(&ok, 10);
    // header
    // second line is only an integer
    if (ok)
    {
        SL = ETRecs[count+2].split(QRegExp("\\s+"));
        if (count == SL.count())
            oldformat = false;
        //if the number of columns equals the integer then new format
        nrStations = count-1;
        // nr stations is count-1 for time as forst column
    }

    if (ETRecs[0].contains("RUU"))
        oldformat = true;

    if (oldformat)
    {
        QStringList SL = ETRecs[0].split(QRegExp("\\s+"));
        // get first line, white space character as split for header

        nrStations = SL[SL.size()-1].toInt(&ok, 10);
        // read nr stations from last value in old style header
        // failure gives 0
        SL = ETRecs[ETRecs.count()-1].split(QRegExp("\\s+"));
        oldformat = (nrStations == SL.count()-1);
    }

    //check if nr stations found equals nr columns-1, 1st column is time
    if (oldformat)
        skiprows = 1;
    else
        skiprows = 3;

    // count gauge areas in the ID.map
    int nrmap = 0;
    nrmap = countUnits(*ETZone);

    if (nrmap > nrStations)
    {
        ErrorString = QString("Number of stations in ET file (%1) < nr of ET zones in ID map (%2)").arg(nrStations).arg(nrmap);
        throw 1;
    }

    nrSeries = ETRecs.size() - nrStations - skiprows;
    // count ET or snowmelt records

    if (nrSeries <= 1)
    {
        ErrorString = "ET records <= 1, must at least have one interval with 2 rows: a begin and end time.";
        throw 1;
    }

    for(int r = 0; r < nrSeries; r++)
    {
        // initialize ET record structure
        rl.time = 0;
        rl.intensity.clear();

        // split ET record row with whitespace
        QStringList SL = ETRecs[r+nrStations+skiprows].split(QRegExp("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);

        if (r == 0)
            time = rl.time;

        if (r > 0 && rl.time <= time)
        {
            //qDebug() << r << time << rl.time;
            ErrorString = QString("ET records at time %1 has unreadable value.").arg(rl.time);
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
                ErrorString = QString("ET records at time %1 has unreadable value.").arg(SL[0]);
                throw 1;
            }
        }

        ETSeriesM << rl;
    }

    // sometimes not an increasing timeseries
    for(int i = 1; i < nrSeries; i++){
        if (ETSeriesM [i].time <= ETSeriesM[i-1].time) {
            ErrorString = QString("ET records time is not increasing at row %1.").arg(i);
            throw 1;
        }
    }

    nrETseries = ETSeriesM.size();//nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetETMap(void)
{
    double currenttime = (time)/60;
    int  ETplace;
    double tt = 0.001; //mm/day to m/day
            //_dt/86400.0*0.001; // mm/day to m / timestep
    bool noET = false;
    bool sameET= false;

//    if (!SwitchIncludeET)
//        return;
//    if (SwitchETSatellite)
//        return;

    // from time t to t+1 the ET is the ET of t

    // where are we in the series
    int currentrow = 0;
    // if time is outside records then use map with zeros
    if (currenttime < ETSeriesM[0].time || currenttime > ETSeriesM[nrETseries-1].time)
        noET = true;

    if (noET) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ETp->Drc = 0;
        }}
    } else {
        // find current record
        for (ETplace = 0; ETplace < nrETseries-1; ETplace++) {
            if (currenttime >= ETSeriesM[ETplace].time && currenttime < ETSeriesM[ETplace+1].time) {
                currentrow = ETplace;
                break;
            }
        }

        if (currentrow == currentETrow && ETplace > 0)
            sameET = true;

        // get the next map from file
        if (!sameET) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                ETp->Drc = ETSeriesM[currentrow].intensity[(int) ETZone->Drc-1]*tt;
            }}
        } //sameET
    }
    currentETrow = currentrow;
}
//---------------------------------------------------------------------------
void TWorld::GetETSatMap(void)
{
    double currenttime = (time)/60;
    int  ETplace;
    double tt = 0.001; //mm/day to m/day
            //_dt/86400.0*0.001; // mm/day to m / timestep

    bool noET = false;
    bool sameET= false;

//    if (!SwitchIncludeET)
//        return;
//    if (!SwitchETSatellite)
//        return;

    // from time t to t+1 the ET is the ET of t

    // where are we in the series
    int currentrow = 0;
    // if time is outside records then use map with zeros
    if (currenttime < ETSeriesMaps[0].time || currenttime > ETSeriesMaps[nrETseries-1].time) {
        noET = true;
        DEBUG("run time outside ET records");
    }

    if (noET) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ETp->Drc = 0;
        }}
    } else {
        // find current record
        for (ETplace = 0; ETplace < nrETseries-1; ETplace++) {
            if (currenttime >= ETSeriesMaps[ETplace].time && currenttime < ETSeriesMaps[ETplace+1].time) {
                currentrow = ETplace;
                break;
            }
        }

        if (currentrow == currentETrow && ETplace > 0)
            sameET = true;
        // get the next map from file
        if (!sameET) {

            auto _M = std::unique_ptr<cTMap>(new cTMap(readRaster(ETSeriesMaps[ETplace].name)));

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (pcr::isMV(_M->Drc)) {
                    QString sr, sc;
                    sr.setNum(r); sc.setNum(c);
                    ErrorString = "Missing value at row="+sr+" and col="+sc+" in map: "+ETSeriesMaps[ETplace].name;
                    throw 1;
                } else {
                    ETp->Drc = _M->Drc *tt;
                }
            }}
        }
    }
    currentETrow = currentrow;
}
//---------------------------------------------------------------------------
void TWorld::doETa()
{
    double ETafactor = 1;
    double Ld = 12;

    SwitchDailyET = true;
    if (SwitchDailyET) {
        double day = trunc(time/(86400));
        double hour = time/3600.0-day*24.0;
        double declination = -23.45 * PI/180.0 * cos(2*M_PI*(day+10)/365.0);
        Ld = 24.0/M_PI*(acos(-tan(declination)*tan(latitude/180.0*M_PI)));  // daylength in hour
        ETafactor = std::max(0.,sin((-0.5-hour/Ld)*M_PI)) / Ld*_dt/3600.0*M_PI*0.5;
            //<= this ensures that the sum of all steps in a day amounts to the daily ET, regardless of _dt
    }
    // sum of ETafactor during Ld is always 1, so ETp is devided with a sine curve over daylength Ld

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        if (Rain->Drc > 0)
            ETp->Drc = 0;

        double ETp_ = ETp->Drc * ETafactor;
        // ETp->Drc is assumed to be in m/day, etfactor does the rest

        if (ETp_ > 0) {
            double tot = 0;
            double tmp = 0;
            double eta = 0;
            double etanet = 0;
            double Cover_ = Cover->Drc;

            // interception decrease, drying out canopy
            double CStor_  = CStor->Drc;
            if (CStor_ > 0) {
                double ETa_int = ETp_;

                tmp = CStor_;
                ETa_int = std::min(ETa_int, CStor_);
                CStor_ = CStor_- ETa_int;

                if (CStor_ == 0)
                    RainCum->Drc = 0;
                // restart the cumulative process when CStor is dried out

                eta = tmp - CStor_; // actual eta is store before-after
                tot = tot + Cover_*eta;
                etanet = Cover_*eta;

                Interc->Drc = Cover_ * CStor_ * SoilWidthDX->Drc * DX->Drc;
                IntercETa->Drc += Cover_ * eta * SoilWidthDX->Drc * DX->Drc;
                CStor->Drc = CStor_;


            }

            bool ponded = hmxWH->Drc > 0;
            if (!ponded) {
                double pore = Poreeff->Drc;
                double theta = Thetaeff->Drc;
                double thetar = 0.025*pore;
                double Lw_ = Lw->Drc;
                // soil moisture transpiration covered surface
                double theta_e = (theta-thetar)/(pore-thetar);
                double f = 1+qPow(theta_e/0.5,6.0);
                double ETa_soil1 = (1.0-1.0/f) * etanet;//ETp_  * Cover_;  //Transpiration
                double ETa_soil2 = theta_e * ETp_ * (1-Cover_);   //Evaporation

                // there is an infiltration front
                if (Lw_ > 0) {
                    double moist = Lw_ * (pore-thetar);
                    eta = std::min(moist, ETa_soil1-ETa_soil2);

                    moist = moist - eta;
                    Lw->Drc = moist/(pore-thetar);
                    tot = tot + eta;
                } else {
                    // soil moisture evaporation bare surface
                    double moist = (theta-thetar) * SoilDepth1->Drc;
                    eta = std::min(moist, ETa_soil1-ETa_soil2);
                    moist = moist - eta;
                    Thetaeff->Drc = moist/SoilDepth1->Drc + thetar;
                    tot = tot + eta;
                }
            }

            // ETa = ETp for any ponded surfaces
//            if (ponded) {
//                double ETa_pond = ETp_;
//                if (FloodDomain->Drc > 0) {
//                    ETa_pond = std::min(ETa_pond, hmx->Drc);
//                    tmp = hmx->Drc;
//                    hmx->Drc = hmx->Drc-ETa_pond;
//                    eta = tmp - hmx->Drc;
//                } else {
//                    double WHRunoff_ = WHrunoff->Drc;
//                    ETa_pond = std::min(ETa_pond, WHRunoff_);
//                    tmp = WHRunoff_;
//                    WHRunoff_ = WHRunoff_-ETa_pond;
//                    eta = tmp - WHRunoff_;
//                    WHroad->Drc = WHRunoff_;
//                    WH->Drc = WHRunoff_ + WHstore->Drc;
//                    WHrunoff->Drc = WHRunoff_;
//                }
//                tot = tot + eta;
//                WaterVolall->Drc = CHAdjDX->Drc * (WHrunoff->Drc + hmx->Drc) + WHstore->Drc*SoilWidthDX->Drc*DX->Drc;
//            }

            // put total Eta in Eta map
            ETa->Drc = tot;
            ETaCum->Drc += tot;
        }
   //     if(r==100&&c==200)
    //        qDebug()<< "e" << ETp->Drc << CStor->Drc << Rainc->Drc << RainCum->Drc << ETa->Drc << InterceptionmmCum->Drc;

    }}
//    report(*Thetaeff,"ti");
//    report(*tma,"ta");
//    report(*tmb,"tb");
//    report(*ETp,"ETp");
}

