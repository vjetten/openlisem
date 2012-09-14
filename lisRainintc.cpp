
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
- void TWorld::Interception(void) \n
 */

#include "model.h"

//---------------------------------------------------------------------------
/// read rainfall files of different types and put data in RainfallSeries
/// reads also old RUU lisem rain files
/// reads rainfall maps
//OBSOLETE
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
//---------------------------------------------------------------------------
/// read rainfall files of different types and put data in RainfallSeries
/// reads also old RUU lisem rain files
/// can read rainfall maps in between intensity values
/// format: first line ends with integer that is nr of data columns excl time
void TWorld::GetRainfallDataM(void)
{
   RAIN_LIST rl;
   QFile fff(rainFileName);
   QFileInfo fi(rainFileName);
   QString S;
   QStringList rainRecs;
   bool ok;
   //int j = 0;

   if (!fi.exists())
   {
      ErrorString = "Rainfall file not found: " + rainFileName;
      throw 1;
   }

   nrstations = 0;
   nrrainfallseries = 0;

   fff.open(QIODevice::ReadOnly | QIODevice::Text);

   while (!fff.atEnd())
   {
      S = fff.readLine();
      if (!S.trimmed().isEmpty())
         rainRecs << S;
      qDebug() << S << "hoi";
   }
   fff.close();
   // get all rainfall records

   QStringList SL = rainRecs[0].split(QRegExp("\\s+"));
   // white space character as split
   nrstations = SL[SL.size()-2].toInt(&ok, 10);

   if (!ok)
   {
      ErrorString = "Cannot read nr rainfall stations in header rainfall file";
      throw 1;
   }
   qDebug() << QString(" nr rainfall stations %1").arg(nrstations);
   nrrainfallseries = rainRecs.size() - nrstations - 1;
   // count rainfall records
   qDebug() << QString(" nr rainfall records %1").arg(nrrainfallseries);

   if (nrrainfallseries <= 1)
   {
      ErrorString = "rainfall records <= 1, must at least have 2 rows: one interval with a begin and end time.";
      throw 1;
   }

   //RainfallSeriesM = new RAIN_LIST[nrrainfallseries];
//   for(int r=0; r < nrrainfallseries; r++)
//   {
//      RainfallSeriesM[r].time = 0;
////      RainfallSeriesM[r].intensity = new double[nrstations+1];
//      for (int j = 0; j <= nrstations; j++)
////         RainfallSeriesM[r].intensity[j] = 0;
//         RainfallSeriesM[r].intensity.clear();
//      RainfallSeriesM[r].isMap = false;
//      RainfallSeriesM[r].name = "";
//   }
   for(int r = 0; r < nrrainfallseries; r++)
   {
      // initialize rainfall record structure
      rl.time = 0;
      rl.intensity.clear();
      rl.isMap = false;
      rl.name = "";

      QStringList SL = rainRecs[r+nrstations+1].split(QRegExp("\\s+"), QString::SkipEmptyParts);
      // split rainfall record row with whitespace

      rl.time = SL[0].toDouble();
      // time in min

      // check if record has characters, then filename assumed
      if (SL[1].contains(QRegExp("[A-Za-z]")))
      {
         QFileInfo fi(QDir(rainFileDir), SL[1]);
         // asume second record is name
         if (!fi.exists())
         {
            ErrorString = QString("rainfall map %1 not found.").arg(SL[1]);
            throw 1;
         }

         rl.name = fi.absoluteFilePath();
         rl.isMap = true;
         // a mapname if the file exists
      }
      else
      {
         // record is a assumed to be a double

         for (int i = 0; i < nrstations; i++)
         {
            bool ok = false;
//            RainfallSeriesM[r].intensity[i] = SL[i+1].toDouble(&ok);
            rl.intensity << SL[i+1].toDouble(&ok);
            if (!ok)
            {
               ErrorString = QString("rainfall records at time %1 has unreadable value.").arg(SL[0]);
               throw 1;
            }

         }
        RainfallSeriesM << rl;
      }

      //qDebug() << RainfallSeriesM[r].time << RainfallSeriesM[r].intensity << RainfallSeriesM[r].name;
   }
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

   for (place = 0; place < nrrainfallseries; place++)
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

   // find the interval in which we are now
//OBSOLETE
//   FOR_ROW_COL_MV
//   {
//      int col = (int) RainZone->Drc;
//      double tt = 3600000.0;

//      Rain->Drc = RainfallSeries[place][col]*_dt/tt;
//      // Rain in m per timestep from mm/h
//      Rainc->Drc = Rain->Drc * _dx/DX->Drc;
//      // correction for slope dx/DX, water spreads out over larger area


//      RainCum->Drc += Rainc->Drc;
//      // cumulative rainfall corrected for slope, used in interception
//      RainNet->Drc = Rainc->Drc;
//   }
}
//---------------------------------------------------------------------------
/// Interception()
/// - interception seen as rigid storage SMax filling up and overflowing\n
/// - overflow flux is identical to rainfall flux in intensity\n
/// - SMax is the storage of the plants inside the gridcell, not the average storage of the gridcell\n
/// - so if a single tree inside a cell has an SMax of 2mm even if it covers 10%, the Smax of that cell is 2\n
/// - therefore the same goes for LAI: the LAI of the plants inside the gridcell\n
/// - this is also easier to observe. The LAI from a satellite image is the average LAI of a cell, must be divided by Cover
void TWorld::Interception(void)
{
   // all variables are in m
   if (!SwitchRainfall)
      return;
   //VJ 110113 bug fix, no interception when no rainfall and only snowmelt

   FOR_ROW_COL_MV
   {
      double CS = CStor->Drc;
      //actual canopy storage in m
      double Smax = CanopyStorage->Drc;
      //max canopy storage in m
      double LAIv;
      if (SwitchInterceptionLAI)
         LAIv = LAI->Drc;
      else
         LAIv = (log(1-Cover->Drc)/-0.4)/max(0.9,Cover->Drc);
      //Smax is based on LAI and LAI is the average of a gridcell, already including the cover
      // a low cover means a low LAI means little interception
      // avoid division by 0

      if (SwitchBuffers && !SwitchSedtrap)
         if(BufferID->Drc > 0)
            Smax = 0;
      // no interception with buffers, but sedtrap can have interception

      if (SwitchHardsurface && HardSurface->Drc > 0)
         Smax =  0;
      //VJ 110111 no interception on hard surfaces

      if (Smax > 0)
      {
         double k = exp(-CanopyOpeness*LAIv);
         CS = Smax*(1-exp(-k*RainCum->Drc/Smax));
         //      CS = Smax*(1-exp(-0.0653*LAIv*RainCum->Drc/Smax));
         //VJ 110209 direct use of openess, astons value too open. A good guess is using the cover LAI relation
         //and interpreting cover as openess: k = exp(-0.45*LAI)
      }
      else
         CS = 0;
      // 0.0653 is canopy openess, based on Aston (1979), based on Merriam (1960/1973), De Jong & Jetten 2003
      // is not the same as canopy cover. it also deals with how easy rainfall drips through the canopy
      // possible to use equation from Ahston but for very open Eucalypt

      CS = max(0, CS * (1-StemflowFraction));
      //VJ 110206 decrease storage with stemflow fraction!

      LeafDrain->Drc = max(0, Cover->Drc*(Rainc->Drc - (CS - CStor->Drc)));
      // diff between new and old strage is subtracted from rainfall
      // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!
      // note: cover already implicit in LAI and Smax, part falling on LAI is cover*rainfall

      CStor->Drc = CS;
      // put new storage back in map
      Interc->Drc =  Cover->Drc * CS * SoilWidthDX->Drc * DX->Drc; //*
      // only on soil surface, not channels or roads, in m3
      // cover already implicit in CS, Smax

      RainNet->Drc = LeafDrain->Drc + (1-Cover->Drc)*Rainc->Drc;
      // net rainfall is direct rainfall + drainage
      // rainfall that falls on the soil, used in infiltration
   }
}
//---------------------------------------------------------------------------
void TWorld::InterceptionHouses(void)
{
   // all variables are in m
   if (!SwitchHouses)
      return;

   FOR_ROW_COL_MV
   {
      if (HouseCover->Drc > 0)
      {
         double HS, DS;
         //actual roof storage in m
         double Hmax = RoofStore->Drc;
         //max roof storage in m
         HouseCover->Drc = qMin(HouseCover->Drc, 0.95);
         double Dmax = DrumStore->Drc/(CellArea->Drc*HouseCover->Drc);;
         //max drum storage in m
         double housedrain = 0;
         //overflow in m

         if (Hmax > 0)
         {
            double k = 0.5;
            // speed of filling of roof storage, very quickly
            HS = Hmax*(1-exp(-k*RainCum->Drc/Hmax));
            //roof storage in m
         }
         else
            HS = 0;

         if (Dmax > 0)
         {
            double k = 0.05;
            // speed of filling of raindrum near house, slower
            DS = Dmax*(1-exp(-k*RainCum->Drc/Dmax));
            //drum storage in m

         }
         else
            DS = 0;

         housedrain = max(0, (RainNet->Drc - (HS - HStor->Drc) - (DS - DStor->Drc)));
         // diff between new and old strage is subtracted from rainfall
         // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!

         HStor->Drc = HS;
         DStor->Drc = DS;
         // put new storage back in maps in m and m

         IntercHouse->Drc = HouseCover->Drc  * (HS+DS) * SoilWidthDX->Drc * DX->Drc;//_dx*_dx
         // total interception in m3

         RainNet->Drc = HouseCover->Drc*housedrain + (1-HouseCover->Drc)*RainNet->Drc;
         // net rainfall is direct rainfall + drainage
         // rainfall that falls on the soil, used in infiltration
      }
      else
         IntercHouse->Drc = 0;
   }
}
//---------------------------------------------------------------------------



