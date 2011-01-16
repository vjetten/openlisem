/*
 * LisUIoutput.h
 *
 *  Created on: Jun 5, 2010
 *      Author: jetten
 */

#ifndef LISUIOUTPUT_H_
#define LISUIOUTPUT_H_


/** structure with variables passed form model to interface each run step.
 this structure is visible both in the UI and in the model world */
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
