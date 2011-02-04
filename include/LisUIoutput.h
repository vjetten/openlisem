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
	double CatchmentArea, dx, t,time, maxtime, EndTime, BeginTime;

	double MB, Qtot, Qtotmm, Qpeak, IntercTotmm, WaterVolTotmm, InfilTotmm,
	RainTotmm, SurfStormm, InfilKWTotmm,
	MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedTot,
	ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot,
	RunoffFraction, RainpeakTime, QpeakTime, Q, Qs, C, P,
	BufferVolTot, BufferSedTot;

	TMMap *DrawMap;

	//bool SwitchErosion;
	//bool SwitchIncludeChannel;
	QString runfilename;
	QString LisemDir;
};



#endif /* LISUIOUTPUT_H_ */
