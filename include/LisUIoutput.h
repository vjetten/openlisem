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
  \file LisUIoutput.h
  \brief structure to pass variables form the model to the interface, visible by both
  */

#ifndef LISUIOUTPUT_H_
#define LISUIOUTPUT_H_

/// structure to pass variables form the model to the interface.
/// This tsructure is the link, visible by both

struct output{
    int runstep;
    int printstep;
    int maxstep;

    int outputpointnr;
    QString outputpointdata;

    int
    F_solution,
    F_scheme,
    F_SSMethod,
    F_BLMethod,
    F_SigmaDiffusion,
    F_fluxLimiter,
    F_replaceV,
    F_Maxiter;
    double
    F_maxVelocity,
    F_extremeHeight,
    F_extremeDiff,
    F_courant,
    F_courant_diffusive;

    double CatchmentArea, dx, t,time, maxtime, EndTime, BeginTime;

    double
    // water
    MB, Qtot, Q, Qtile, Qpeak, RunoffFraction, RainpeakTime, QpeakTime,
    Qtotmm,  IntercTotmm, IntercHouseTotmm, WaterVolTotmm, InfilTotmm,
    RainTotmm, SurfStormm, InfilKWTotmm, Pmm,
    // channel
    ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot, ChannelWH,
    // flood
    FloodTotMax, FloodAreaMax, WHflood, Qflood, volFloodmm,
    // sediment
    MBs, Qs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedTot, C,
    Cplot, Qsplot,
    // buffer
    BufferVolTot, BufferSedTot,
    // screen output
    QtotPlot, SoilLossTotPlot, QpeakPlot, QPlot;


    cTMap *baseMap;
    cTMap *baseMapDEM;
    cTMap *channelMap;
    cTMap *roadMap;
    cTMap *houseMap;

    QList<double> graindiameters;

    //combox selection of drawn map
    QList<cTMap *> ComboMaps;
    QList<cTMap *> ComboMapsSafe;
    QList<QList<double>> ComboColorMap;
    QList<QList<QString>> ComboColors;
    QList<bool> ComboLogoritmic;
    QList<bool> ComboSymColor;
    QStringList ComboMapNames;
    QStringList ComboUnits;
    QList<double> ComboScaling;

    bool comboboxset;

    QString runfilename;
    QString LisemDir;

    bool doBatchmode;
};



#endif /* LISUIOUTPUT_H_ */
