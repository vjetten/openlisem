/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file LisUIoutput.h
  \brief structure to pass variables form the model to the interface, visible by both
  */

#include <CsfMap.h>
#include <CsfRGBMap.h>
#include <QList>

#ifndef LISUIOUTPUT_H_
#define LISUIOUTPUT_H_

/// structure to pass variables form the model to the interface.
/// This tsructure is the link, visible by both



struct output{
    int runstep;
    int printstep;
    int maxstep;

    QList<int> OutletIndices;
    QList<int> OutletLocationX;
    QList<int> OutletLocationY;
    QList<QList<double>*> OutletQ;
    QList<QList<double>*> OutletQs;
    QList<QList<double>*> OutletC;
    QList<QList<double>*> OutletChannelWH;
    QList<double> OutletQpeak;
    QList<double> OutletQpeaktime;
    QList<double> OutletQtot;
    QList<double> OutletQstot;

    QVector <double> ChanDataX;
    QVector <double> ChanDataY;
    QVector<int> Chanbranch;
    QList <int> branches;
    QVector <double> CulvertX;
    QVector <double> CulvertY;
    QVector <double> EndPointX;
    QVector <double> EndPointY;

    double timestep, CatchmentArea, dx, t,time, maxtime, EndTime, BeginTime;

    double
    // water
    MB, Qtot,  Qtile, Qtiletot, RunoffFraction, RainpeakTime,
    Qtotmm,  IntercTotmm, IntercHouseTotmm, WaterVolTotmm,InfilTotmm,StormDrainTotmm,
    RainTotmm, SurfStormm, InfilKWTotmm, Pmm, BaseFlowtotmm,IntercLitterTotmm,WaterVolTotchannelmm,
    floodBoundaryTot, floodBoundarySedTot,
    // channel
    ChannelVolTotmm, ChannelSedTot, ChannelDepTot, ChannelDetTot, ChannelWH,
    // flood
    FloodTotMax, FloodAreaMax, FloodArea, WHflood, Qflood, volFloodmm,
    FloodDetTot, FloodDepTot, FloodSedTot,
    // sediment
    MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedTot;

    // map pointers for display
    cTMap *baseMap;
    cTMap *baseMapDEM;
    cTMap *channelMap;
    cTMap *outletMap;
    cTMap *roadMap;
    cTMap *houseMap;
    cTMap *flowbarriersMap;
    cTRGBMap *Image;

    QList<double> graindiameters;

    QList<int> ComboLists;
    QList<cTMap *> ComboMaps;
    QList<cTMap *> ComboMapsSafe;
    QList<QList<double>> ComboColorMap;
    QList<QList<QString>> ComboColors;
    QList<bool> ComboLogaritmic;
    QList<bool> ComboSymColor;
    QStringList ComboMapNames;
    QStringList ComboUnits;
    QList<double> ComboScaling;
    QList<double> userMinV;
    QList<double> userMaxV;
    QList<double> comboStep;

    bool comboboxset;
    bool has_image;

    QString runfilename;
    QString LisemDir;
    QString format;
    QString timeStartRun;
    QString datestamp;

    bool doBatchmode;
};



#endif /* LISUIOUTPUT_H_ */
