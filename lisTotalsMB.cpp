/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem

Functionality in lisTotalMB.cpp:
- calculating totals for output and mass balance
- mass balance calculations
---------------------------------------------------------------------------*/

#include "model.h"


//---------------------------------------------------------------------------
// totals for screen and file output and mass balance
void TWorld::Totals(void)
{
    double rainfall, snowmelt;
    double oldrainpeak, oldsnowpeak; 
    
	/***** WATER *****/
    
    if (SwitchRainfall)
    {
        RainAvgmm = Rain->MapTotal()*1000/nrCells;
        RainTotmm += RainAvgmm;        
        // avg area rainfall in mm
        
        tm->calc2V(Rain, (_dx*_dx), MUL); //in m3
        rainfall = tm->MapTotal();
        RainTot += rainfall; // in m3
        
        oldrainpeak = Rainpeak;
        Rainpeak = max(Rainpeak, rainfall);
        if (oldrainpeak  < Rainpeak)
            RainpeakTime = time;
    }
    
	if (SwitchSnowmelt)
	{
        SnowAvgmm += Snowmelt->MapTotal()*1000/nrCells;
        SnowTotmm += SnowAvgmm;
        
		tm->calc2V(Snowmelt, (_dx*_dx), MUL); //in m3
		snowmelt = tm->MapTotal();
		SnowTot += snowmelt; // in m3
        
        oldsnowpeak = Snowpeak;
		Snowpeak = max(Snowpeak, snowmelt);
		if (oldsnowpeak < Snowpeak)
			SnowpeakTime = time;
	}
    
	IntercTot = Interc->MapTotal();
	IntercTotmm = IntercTot*1000/(_dx*_dx*nrCells);//CatchmentArea;
	// interception in mm and m3
    
	InfilTot += InfilVol->MapTotal() + InfilVolKinWave->MapTotal(); //m3
	InfilKWTot += InfilVolKinWave->MapTotal(); // not really used, available for output when needed
	InfilTotmm = max(0,InfilTot*1000/(_dx*_dx*nrCells));//CatchmentArea;
	// infiltration mm and m3
    
	tm->calc2V(WHstore, 1000, MUL); //mm
	SurfStoremm = tm->MapTotal()/nrCells;
	// surface storage CHECK THIS
	// does not go to MB, is already in tot water vol
    
	WaterVolTot = WaterVolall->MapTotal();//m3
	WaterVolTotmm = WaterVolTot*1000/(_dx*_dx*nrCells);//CatchmentArea; //mm
	// water on the surface in runoff in m3 and mm
	//NOTE: surface storage is already in here so does not need to be accounted for in MB
    
	Qtot += Qoutflow->MapTotal();
	// sum outflow m3 for all timesteps for all pits, is already mult by dt
	// needed for mass balance
	Qtotmm = Qtot*1000/(_dx*_dx*nrCells);//CatchmentArea;
	// in mm for screen output
    
	QtotOutlet += Qoutflow->DrcOutlet;
	// for screen output, total main outlet in m3
	TotalWatervol->copy(WaterVolall);
	// for sed conc calc output
    
    
	if (SwitchIncludeChannel)
	{
		WaterVolTot += ChannelWaterVol->MapTotal(); //m3
		// add channel vol to total
		WaterVolTotmm = WaterVolTot*1000/(_dx*_dx*nrCells);//CatchmentArea; //mm
		// recalc in mm for screen output
        
		Qtot += ChannelQoutflow->MapTotal();
		// add channel outflow (in m3) to total for all pits
		Qtotmm = Qtot*1000/(_dx*_dx*nrCells);//CatchmentArea;
		// recalc in mm for screen output
        
		QtotOutlet += ChannelQoutflow->DrcOutlet;
		// add channel outflow (in m3) to total for main outlet
		TotalWatervol->calc(ChannelWaterVol,ADD);
		// add channel volume to total for sed conc calc
        
	}
    
	if (SwitchBuffers)
	{
		BufferVolTot = BufferVol->MapTotal(); // in m3
		if (SwitchIncludeChannel)
			BufferVolTot += ChannelBufferVol->MapTotal();
		//sum up all volume remaining in all buffers (so the non-water!)
		BufferVolTotInit = BufferVolInit->MapTotal() + ChannelBufferVolInit->MapTotal();
		BufferVolTot = BufferVolTotInit - BufferVolTot;
		//subtract this from the initial volume to get the water in the buffers
	}
    
	// output fluxes for reporting
	FOR_ROW_COL_MV
	{
		Qoutput->Drc = 1000*(Qn->Drc + ChannelQn->Drc); // in l/s
	}
    
	oldrainpeak = Qpeak;
	Qpeak = max(Qpeak, Qoutput->DrcOutlet);
	if (oldrainpeak < Qpeak)
		QpeakTime = time;
	// peak flow and peak time calculation, based on sum channel and runoff
    
	/***** SEDIMENT *****/
    
    if (SwitchErosion)
	{
		DetSplashTot += DETSplash->MapTotal();
		DetFlowTot += DETFlow->MapTotal();
		DepTot += DEP->MapTotal();
		DetTot += DETSplash->MapTotal() + DETFlow->MapTotal();
		SedTot = Sed->MapTotal();
		// all in kg/cell
        
		SoilLossTot += Qsoutflow->MapTotal();
		// sum all sed in all pits (in kg), needed for mass balance
        
		SoilLossTotOutlet += Qsoutflow->DrcOutlet;
		// for screen output, total main outlet sed loss
		TotalSed->copy(Sed);
		// for sed conc
        
		if (SwitchIncludeChannel)
		{
			ChannelDetTot += ChannelDetFlow->MapTotal();
			ChannelDepTot += ChannelDep->MapTotal();
			ChannelSedTot = ChannelSed->MapTotal();
            
			SoilLossTot += ChannelQsoutflow->MapTotal();
			// add sed outflow for all pits to total soil loss
            
			SoilLossTotOutlet += ChannelQsoutflow->DrcOutlet;
            // add channel outflow (in kg) to total for main outlet
            
			TotalSed->calc(ChannelSed, ADD);
            // for sed conc file output
		}
        
		FOR_ROW_COL_MV
		{
			TotalConc->Drc = min(850,(TotalWatervol->Drc > 1e-6? TotalSed->Drc/TotalWatervol->Drc : 0));
		}
		// for file output
        
		if (SwitchBuffers || SwitchSedtrap)
		{
			BufferSedTot = BufferSed->MapTotal();
			if (SwitchIncludeChannel)
				BufferSedTot += ChannelBufferSed->MapTotal();
			BufferSedTot = BufferSedTotInit - BufferSedTot;
		}
      /** TODO add gully, wheeltracks etc */
        
      // spatial totals for output all in kg/cell
		FOR_ROW_COL_MV
		{
			Qsoutput->Drc = Qsn->Drc + ChannelQsn->Drc;  // sum channel and OF sed output in kg/s
            
			TotalDetMap->Drc += DETSplash->Drc + DETFlow->Drc;
			TotalDepMap->Drc += DEP->Drc;
			if (SwitchIncludeChannel)
			{
				TotalDetMap->Drc += ChannelDetFlow->Drc;
				TotalDepMap->Drc += ChannelDep->Drc;
			}
			TotalSoillossMap->Drc = TotalDetMap->Drc + TotalDepMap->Drc;
		}
	}
}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
	// Mass Balance water
	if (RainTot + SnowTot > 0)
		MB = (RainTot + SnowTot - IntercTot - InfilTot - WaterVolTot
              - BufferVolTot - Qtot)/(RainTot + SnowTot)*100;
    
	// Mass Balance sediment
	if (SwitchErosion && DetTot > 0)
		MBs = (DetTot + ChannelDetTot - SoilLossTot - SedTot - ChannelSedTot +
               DepTot + ChannelDepTot - BufferSedTot)/DetTot*100;
}
//---------------------------------------------------------------------------
