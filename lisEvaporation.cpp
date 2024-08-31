
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
            S.remove(S.size()-1,1);
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
        SL = ETRecs[count+2].split(QRegularExpression("\\s+"));
        nrStations = count-1;
        // nr stations is count-1 for time as first column
    }

    skiprows = 3;

    // count gauge areas in the ID.map
    QList < int> tmp;
    tmp = countUnits(*ETZone);
    int nrmap = tmp.count();

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
        int r_ = r+nrStations+skiprows;

        // split ET record row with whitespace
        QStringList SL = ETRecs[r_].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

        // read date time string and convert to time in minutes
        rl.time = getTimefromString(SL[0]);
        time = rl.time;

        // check is time is increasing with next row
        if (r+1 < nrSeries) {
            QStringList SL1 = ETRecs[r_+1].split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            int time1 = getTimefromString(SL1[0]);
            if (time1 < time) {
                ErrorString = QString("Time in evaporation records is not increasing from row %1 to %2. Check your file!").arg(r_).arg(r_+1);
                throw 1;
            }
        }

//        if (r == 0)
//            time = rl.time;

//        if (r > 0 && rl.time <= time)
//        {
//            ErrorString = QString("ET records at row %1 has unreadable value.").arg(r);
//            throw 1;
//        }
//        else
//            time = rl.time;

        // record is a assumed to be a double
        for (int i = 1; i <= nrStations; i++)
        {
            bool ok = false;
            rl.intensity << SL[i].toDouble(&ok);
            if (!ok)
            {
                ErrorString = QString("ET records at time %1 has unreadable value: %2.").arg(SL[0]).arg(SL[i]);
                throw 1;
            }
        }

        ETSeries << rl;
    }

    // sometimes not an increasing timeseries
//    for(int i = 1; i < nrSeries; i++){
//        if (ETSeries [i].time <= ETSeries[i-1].time) {
//            ErrorString = QString("ET records time is not increasing at row %1.").arg(i);
//            throw 1;
//        }
//    }

    nrETseries = ETSeries.size();//nrSeries;
}
//---------------------------------------------------------------------------
void TWorld::GetETMap(void)
{
    double currenttime = (time)/60; //time in min
    double tt = 0.001*ETBiasCorrection; //mm to m
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
      //  qDebug() << currentrow << ETSeries[currentrow].intensity[0];
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ETp->Drc = ETSeries[currentrow].intensity[(int) ETZone->Drc-1]*tt;

        }}
    } //sameET
   // report(*ETp,"etp");
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
double TWorld::getETaFactor()
{
    double day = floor(time/86400.0);
    double hour = std::min(24.0,std::max(0.0, time/3600.0-day*24.0));
    double Ld = (2.0*acos(-tan(latitude*0.01745329) * tan(asin(0.397789 * sin(0.017214*(day-1)))))) * 3.8197186;
    double ETafactor = std::max(0.0,sin((-0.5-hour/Ld)*PI)) / Ld*_dt/3600.0*PI*0.5;

    return ETafactor;
}
//---------------------------------------------------------------------------
void TWorld::cell_ETa(int r, int c)
{
    tma->Drc = 0;
    if (Rain->Drc* 3600000.0/_dt > rainfallETa_threshold) {
        ETa->Drc = 0;
        ETp->Drc = 0;
    }

  //      if (r==200 && c == 200)
    //       qDebug() << time/60 << ETp->Drc << ETafactor << Rain->Drc*3600000.0/_dt;

    if (ETp->Drc*ETafactor > 0) {
        double AreaSoil = SoilWidthDX->Drc * DX->Drc;
        double Cover_ = Cover->Drc;
        double ETp_ = ETp->Drc * ETafactor; // potential ETp
        double tot = 0;
        double etanet = ETp_;
        double eta = 0;

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

            etanet = ETp_;//std::max(0.0, ETp_ - ETa_int);

            Interc->Drc = Cover_ * CStor_ * AreaSoil;
            IntercETa->Drc += Cover_ * ETa_int * AreaSoil;
            CStor->Drc = CStor_;
        }

        if (SwitchHouses)
        {
            double CvH = HouseCover->Drc;
            double HS = HStor->Drc;

            double ETa_int = std::min(etanet, HS);
            etanet = ETp_;//std::max(0.0, ETp_ - ETa_int);
            HStor->Drc = HS - ETa_int;
            IntercETa->Drc += CvH * ETa_int * AreaSoil;
            double roofsurface = (_dx * DX->Drc * CvH); // m2
            IntercHouse->Drc =  roofsurface * HS;
        }

        // on ground level energy is shared
        if (SwitchLitter) {
            double CvL = Litter->Drc;
            double LCS = LCStor->Drc;

            double ETa_int = std::min(etanet, LCS);
            etanet = ETp_; //std::max(0.0, ETp_ - ETa_int);
            // we assume the canopy has a different input higher so the ETP for ground is not affected by the canopy
            LCStor->Drc = LCS- ETa_int;
            IntercETa->Drc += CvL * ETa_int * CHAdjDX->Drc;
            LInterc->Drc =  CvL * LCS * CHAdjDX->Drc;
        }

//            if (r==96 && c == 164)
//                qDebug() << ETp_ << CStor_ << RainCum->Drc << Interc->Drc;
        bool ponded = hmxWH->Drc > 0.01;

        if (!ponded) {
            double pore = Poreeff->Drc;
            double theta = Thetaeff->Drc;
            double thetar = ThetaR1->Drc;
            double thetafc = ThetaFC1->Drc;
            double Lw_ = Lw->Drc;
            double theta_e = (theta-thetar)/(pore-thetar);
            double f = 1.0/(1.0+qPow(theta_e/0.4,8.0));
            double ETa_soil = 0;

            if (SwitchLitter)
                ETa_soil = (1.0-f)*ETp_*Cover_ + theta_e*ETp_*(1-Litter->Drc);   //Transpiration + surface Evaporation
            else
                ETa_soil = (1.0-f)*ETp_*Cover_ + theta_e*ETp_*(1-Cover_);   //Transpiration + Evaporation

            // adjust soil moisture because of ETa
            // there is an infiltration front
            if (Lw_ > 0) {
                if(Lw_ < SoilDepth1->Drc) {
                    double moist = Lw_ * (pore-theta);
                    eta = std::min(moist, ETa_soil);
                    moist = moist - eta; // new mositure
                    Lw->Drc = moist/(pore-theta);  // new wetting front
                    tot = tot + eta;
                    tma->Drc += eta;
                    // adjust moisture content layer 1
                    double dL = std::max(0.0,Lw_-Lw->Drc);
                    double m1 = dL*thetafc; // av moist freed layer
                    double m2= (SoilDepth1->Drc-Lw_)*(Thetaeff->Drc-thetar);
                    // av moist below old wetting front
                    Thetaeff->Drc = (m1+m2)/(SoilDepth1->Drc - Lw->Drc)+thetar;
                    // new average moisture below new WF
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
          //  if(c == 200 && r == 200)
            //    qDebug() << Thetaeff->Drc << ThetaI2->Drc << Lw->Drc << SoilDepth1->Drc;

        }
        // ETa = ETp for any ponded surfaces
        if (WHrunoff->Drc > 0.01) {
            double ETa_pond = ETp_;
            double WHRunoff_ = WHrunoff->Drc;
            ETa_pond = std::min(ETa_pond, WHRunoff_);
            WHRunoff_ = WHRunoff_ - ETa_pond;
            eta = ETa_pond;
            //WHroad->Drc = WHRunoff_;
            WH->Drc = WHRunoff_ + WHstore->Drc;
            WHrunoff->Drc = WHRunoff_;
            tot = tot + eta;
            WaterVolall->Drc = CHAdjDX->Drc * (WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;
        }

        // put total Eta in Eta map
        ETa->Drc = tot;
        ETaCum->Drc += tot;
    }
    //SoilETMBcorrection += MapTotal(*tma); // ET water coming from the soil, in m
}

