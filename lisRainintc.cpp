/*---------------------------------------------------------------------------
project: openLISEM
name: lisRainintc.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisRainintc.cpp:
- read rainfall file
- transfer rainfall to a map
- calc interception
---------------------------------------------------------------------------*/

#include "model.h"

//---------------------------------------------------------------------------
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
    QStringList SL = S.split(QRegExp("\\s+"));
    nrstations = SL[SL.size()-2].toInt(&ok, 10);
    if (!ok)
    {
      ErrorString = "Cannot read nr rainfall stations in header rainfall file";
      throw 1;
    }
    for(int r=0; r < nrstations+1; r++)
        S = fff.readLine();

    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (!S.trimmed().isEmpty())
          nrrainfallseries++;
    }
    if (nrrainfallseries <= 1)
    {
      ErrorString = "rainfall records <= 1!";
      throw 1;
    }

    nrrainfallseries++;
    RainfallSeries = new double*[nrrainfallseries];
    for(int r=0; r < nrrainfallseries; r++)
        RainfallSeries[r] = new double[nrstations+1];
    fff.close();
    //RainfallSeries is matrix with rows is data and 1st column is time, other columns are stations

    fff.open(QIODevice::ReadOnly | QIODevice::Text);
    S = fff.readLine();
    for (int i; i < nrstations; i++)
        S = fff.readLine();
    // read header data
    while (!fff.atEnd())
    {
        S = fff.readLine();
        if (S.trimmed().isEmpty()) continue;

        QStringList SL = S.split(QRegExp("\\s+"));
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
void TWorld::Rainfall(void)
{
     double timemin = time / 60;
     double timeminp = (time-_dt) / 60;
     int placep, place;

     for (placep = 0; placep < nrrainfallseries; placep++)
         if (timeminp < RainfallSeries[placep][0])
            break;
     for (place = 0; place < nrrainfallseries; place++)
         if (timemin < RainfallSeries[place][0])
            break;

     FOR_ROW_COL_MV
     {
         int col = (int) RainZone->Drc;

         Rain->Drc = RainfallSeries[place][col]/3600000 * _dt * _dx/DX->Drc;
         // Rain in m per timestep
         //TO DO: weighted average if dt larger than table dt
         // correction for slope dx/DX

         RainCum->Drc += Rain->Drc;
         // cumulative rainfall
      }
}
//---------------------------------------------------------------------------
void TWorld::Interception(void)
{
    FOR_ROW_COL_MV
    {
        double CS = CStor->Drc;
        //actual canopy storage in m
        double rain = Rain->Drc;
        double Smax = CanopyStorage->Drc;
        //max canopy storage in m
        double drain = 0;

        if (Smax > 0)
           CS = max(0, Smax*(1-exp(-0.9*RainCum->Drc/Smax)));
        else
           CS = 0;
        // 0.9 is canopy openess, is not the same as canopy cover

        drain = max(0, rain - (CS - CStor->Drc));
        // diff between new and old strage is subtracted from rainfall
        // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!

        CStor->Drc = CS;
        // put new storage back in map
        Interc->Drc = CS * Cover->Drc * SoilWidthDX->Drc * DX->Drc;
        // only on soil surface, not channels or roads, in m3

        RainNet->Drc = Cover->Drc*drain + (1-Cover->Drc)*rain;
        // net rainfall is direct rainfall + drainage
        // rainfall that falls on the soil
    }
}
//---------------------------------------------------------------------------



