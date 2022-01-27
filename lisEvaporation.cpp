
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
        nrStations = count-1;
        // nr stations is count-1 for time as first column
    }

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
    // count ET records

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
            ErrorString = QString("ET records at row %1 has unreadable value.").arg(r);
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

        ETSeries << rl;
    }

    // sometimes not an increasing timeseries
    for(int i = 1; i < nrSeries; i++){
        if (ETSeries [i].time <= ETSeries[i-1].time) {
            ErrorString = QString("ET records time is not increasing at row %1.").arg(i);
            throw 1;
        }
    }

    nrETseries = ETSeries.size();//nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetETMap(void)
{
    double currenttime = (time)/60; //time in min
    double tt = 0.001; //mm to m
    bool sameET= false;
    // from time t to t+1 the ET is the ET of t

    // if time is outside records then use map with zeros
    if (currenttime < ETSeries[0].time || currenttime > ETSeries[nrETseries-1].time) {
        DEBUG("run time outside ET records");
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ETp->Drc = 0;
        }}
        return;
    }

    // where are we in the series
    int currentrow = ETplace;
    // find current record
    while (currenttime >= ETSeries[ETplace].time
           && currenttime < ETSeries[ETplace+1].time)
    {
        currentrow = ETplace;
        ETplace++;
    }

    if (currentrow == currentETrow && currentrow > 0)
        sameET = true;

    // get the next map from file
    if (!sameET) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ETp->Drc = ETSeries[currentrow].intensity[(int) ETZone->Drc-1]*tt;

        }}
    } //sameET
    currentETrow = currentrow;
}
//---------------------------------------------------------------------------
void TWorld::GetETSatMap(void)
{
    double currenttime = (time)/60; //time in min
    double tt = 0.001*ETBiasCorrection; //mm/day to m/day
    bool noET = false;
    bool sameET= false;

    // from time t to t+1 the ET is the ET of t

    // if time is outside records then use map with zeros
    if (currenttime < ETSeriesMaps[0].time || currenttime > ETSeriesMaps[nrETseries-1].time) {
        DEBUG("run time outside ET records");
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ETp->Drc = 0;
        }}
        return;
    }

// where are we in the series
    int currentrow = ETplace;
    // find current record
    while (currenttime >= ETSeriesMaps[ETplace].time
           && currenttime < ETSeriesMaps[ETplace+1].time)
    {
        currentrow = ETplace;
        ETplace++;
    }

    if (currentrow == currentETrow && currentrow > 0)
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
                ETp->Drc = std::max(0.0,_M->Drc *tt);
            }
        }}
    }

    currentETrow = currentrow;
}
//---------------------------------------------------------------------------
void TWorld::doETa()
{
    double ETafactor = 1;
    double Ld = 12;


   // SwitchDailyET = true;
    if (SwitchDailyET) {
        double day = trunc(time/(86400));
        double hour = time/3600.0-day*24.0;
        double declination = -23.45 * M_PI/180.0 * cos(2*M_PI*(day+10)/365.0);
        Ld = 24.0/M_PI*(acos(-tan(declination)*tan(latitude/180.0*M_PI)));  // daylength in hour
        ETafactor = std::max(0.,sin((-0.5-hour/Ld)*M_PI)) / Ld*_dt/3600.0*M_PI*0.5;
            //<= this ensures that the sum of all steps in a day amounts to the daily ET, regardless of _dt
        //qDebug() << ETafactor << Ld << declination;
    }
    // sum of ETafactor during Ld is always 1, so ETp is devided with a sine curve over daylength Ld

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tma->Drc = 0;
        if (Rain->Drc* 3600000.0/_dt > rainfallETa_threshold) {
            ETa->Drc = 0;
            ETp->Drc = 0;
        }

      //  if (r==200 && c == 200)
        //   qDebug() << time/60 << ETp->Drc << ETafactor << Rain->Drc*3600000.0/_dt;

        if (ETp->Drc > 0) {
            double eta = 0;
            double AreaSoil = SoilWidthDX->Drc * DX->Drc;
            double tot = 0;
            double etanet = 0;
            double Cover_ = Cover->Drc;
            double ETp_ = ETp->Drc * ETafactor;

            ETpCum->Drc += ETp_;

           //  interception decrease, drying out canopy
            double CStor_  = CStor->Drc;
            if (CStor_ > 0) {
                double ETa_int = ETp_;

                ETa_int = std::min(ETa_int, CStor_);
                CStor_ = CStor_- ETa_int;

                RainCum->Drc = std::max(0.0, RainCum->Drc-ETa_int);
                if (CStor_ < 1e-5)
                   RainCum->Drc = 0;

                // restart the cumulative process when CStor is dried out

                etanet = std::max(0.0, ETp_ - ETa_int);

                Interc->Drc = Cover_ * CStor_ * AreaSoil;
                IntercETa->Drc += Cover_ * ETa_int * AreaSoil;
                CStor->Drc = CStor_;
            }

            if (SwitchLitter) {
                double CvL = Litter->Drc;
                double LCS = LCStor->Drc;

                double ETa_int = std::min(etanet, LCS);
                etanet = std::max(0.0, ETp_ - ETa_int);
                LCStor->Drc = LCS- ETa_int;
                IntercETa->Drc += CvL * ETa_int * AreaSoil;
                LInterc->Drc =  CvL * LCS * AreaSoil;
            }

            if (SwitchHouses)
            {
                double CvH = HouseCover->Drc;
                double HS = HStor->Drc;

                double ETa_int = std::min(etanet, HS);
                etanet = std::max(0.0, ETp_ - ETa_int);
                HStor->Drc = HS - ETa_int;
                IntercETa->Drc += CvH * ETa_int * AreaSoil;
                double roofsurface = (_dx * DX->Drc * CvH); // m2
                IntercHouse->Drc =  roofsurface * HS;
            }

//            if (r==96 && c == 164)
//                qDebug() << ETp_ << CStor_ << RainCum->Drc << Interc->Drc;
            bool ponded = hmxWH->Drc > 0;

            if (!ponded) {
                double pore = Poreeff->Drc;
                double theta = Thetaeff->Drc;
                double thetar = ThetaR1->Drc;
                double Lw_ = Lw->Drc;
                double theta_e = (theta-thetar)/(pore-thetar);
                double f = 1.0/(1.0+qPow(theta_e/0.5,6.0));
                //double ETa_soil = theta_e*ETp_;
                double ETa_soil = (1.0-f)*etanet*Cover_ + theta_e*ETp_*(1-Cover_);   //Transpiration + Evaporation

                // there is an infiltration front
                if (Lw_ > 0) {
                    if(!SwitchTwoLayer || Lw_ < SoilDepth1->Drc) {
                        double moist = Lw_ * (pore-thetar);
                        eta = std::min(moist, ETa_soil);
                        moist = moist - eta;
                        Lw->Drc = moist/(pore-thetar);
                        tot = tot + eta;
                        tma->Drc += eta;
                    } else {
                        if (SwitchTwoLayer){
                            double moist = (Lw_-SoilDepth1->Drc) * (ThetaS2->Drc-ThetaR2->Drc);
                            eta = std::min(moist, ETa_soil);
                            moist = moist - eta;
                            Lw->Drc = moist/(ThetaS2->Drc-ThetaR2->Drc)+SoilDepth1->Drc;
                            tot = tot + eta;
                            tma->Drc += eta;
                        }
                    }
                } else {
                    // soil moisture evaporation dry surface
                    double moist = (theta-thetar) * SoilDepth1->Drc;
                    eta = std::min(moist, ETa_soil);
                    moist = moist - eta;
                    Thetaeff->Drc = moist/SoilDepth1->Drc + thetar;
                    tot = tot + eta;
                    tma->Drc += eta;
                }
            }

            // ETa = ETp for any ponded surfaces
            if (WHrunoff->Drc > 0.01) {
                double ETa_pond = ETp_;

                // if kin wave + overflow is used
//                if (FloodDomain->Drc > 0) {
//                    ETa_pond = std::min(ETa_pond, hmx->Drc);
//                    hmx->Drc = hmx->Drc-ETa_pond;
//                    eta = ETa_pond;

//                    hmxflood->Drc = std::max(0.0, WHrunoff->Drc + hmx->Drc - minReportFloodHeight);
//                    FloodWaterVol->Drc = hmxflood->Drc * CHAdjDX->Drc;
//                    double WHrunoffOutput = std::min(WHrunoff->Drc + hmx->Drc, minReportFloodHeight);
//                    RunoffWaterVol->Drc = WHrunoffOutput * CHAdjDX->Drc;
//                }

                double WHRunoff_ = WHrunoff->Drc;
                ETa_pond = std::min(ETa_pond, WHRunoff_);
                WHRunoff_ = WHRunoff_ - ETa_pond;
                eta = ETa_pond;
                WHroad->Drc = WHRunoff_;
                WH->Drc = WHRunoff_ + WHstore->Drc;
                WHrunoff->Drc = WHRunoff_;

//                ETa_pond = std::min(ETa_pond, WH->Drc);
//                WH->Drc = WH->Drc - ETa_pond;
//                WHroad->Drc = std::max(0.0, WHroad->Drc -ETa_pond);

//                if (WH->Drc < WHstore->Drc) {
//                    WHrunoff->Drc = 0;
//                    WHstore->Drc = WH->Drc;
//                } else {
//                    WHrunoff->Drc = WH->Drc - WHstore->Drc;
//                }

                tot = tot + eta;
                WaterVolall->Drc = CHAdjDX->Drc * (WHrunoff->Drc + hmx->Drc) + WHstore->Drc*SoilWidthDX->Drc*DX->Drc;
            }

            // put total Eta in Eta map
            ETa->Drc = tot;
            ETaCum->Drc += tot;
        }
    }}
    SoilETMBcorrection += MapTotal(*tma); // ET water coming from the soil, in m
}

// calc average soil moisture content for output to screen and folder
void TWorld::avgTheta()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Lw_ = Lw->Drc;
        double SoilDep1 = SoilDepth1->Drc;
        ThetaI1a->Drc = Thetaeff->Drc;

        if (Lw_ > 0 && Lw_ < SoilDep1 - 1e-3) {
            double f = Lw_/SoilDep1;
            ThetaI1a->Drc = f * ThetaS1->Drc + (1-f) *Thetaeff->Drc;
        }
        if (Lw_ > SoilDep1 - 1e-3)
            ThetaI1a->Drc = ThetaS1->Drc;

        if (SwitchTwoLayer) {
            double SoilDep2 = SoilDepth2->Drc;
            ThetaI2a->Drc = ThetaI2->Drc;
            if (Lw_ > SoilDep1 && Lw_ < SoilDep2 - 1e-3) {
                double f = (Lw_-SoilDep1)/(SoilDep1-SoilDep1);
                ThetaI2a->Drc = f * ThetaS2->Drc + (1-f) *ThetaI2->Drc;
            }
            if (Lw_ > SoilDep2 - 1e-3)
                ThetaI2a->Drc = ThetaS2->Drc;
        }
    }}
}

