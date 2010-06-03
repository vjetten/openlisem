// METHOD 3
// Holtan/Overton infiltration


//Report("using the Holtan equation\n");

// initialize static maps

_static_spatial_init(REAL4, FC ,"fc.map");
_static_spatial_init(REAL4, A  ,"a.map");
 //   FC and A in mm/h
_static_spatial_init(REAL4, TP ,"tp.map");
 //   TP as a fraction, e.g. 0.45
_static_spatial_init(REAL4, DF ,"df.map");
 //   DF in mm
_static_spatial_init(REAL4, P  ,"p.map");
 //   P , e.g. 0.65 for silt loam

_static_spatial(REAL4, PIV);
_static_spatial(REAL4, GWC);

if (stepNr == 0)
begin{
	 // initialize at first timestep
	_spatial_input(REAL4, ASM,"asm.map");
   celltest(LDD, ASM);
	_spatial_input(REAL4, FP ,"fp.map");
   celltest(LDD, FP);
        // ASM and FP as a fraction of TP, e.g. 0.65 ( = 0.65*0.45 = 0.2925 moist)
        // only used to init PIV and GWC

	 calc(" PIV = TP*DF * (1.0-ASM)");
	 calc(" GWC = TP*DF * (1.0-FP)");

	 // Holtan error messages
	 rangetest(FC ,R_GE_LE, 0    , 100, "Infiltration rate at saturation");
	 rangetest(A  ,R_GE_LE, 0    , 100, "Difference initial and saturation rate");
	 rangetest(TP ,R_GE_LE, 0.001, 0.60,"Soil Porosity");
	 rangetest(DF ,R_GE_LE, 0.001, 300, "Control Zone Depth (mm)");
	 rangetest(P  ,R_GE_LE, 0.4  , 0.8, "Infiltration coefficient P");
	 rangetest(ASM,R_GE_LE, 0.2  , 1,   "Initial Soil Moisture Content");
	 rangetest(FP ,R_GE_LE, 0.4  , 0.9, "Field Capacity Value");

}end     // init code

_spatial(REAL4, FMAX);
calc(" FMAX = FC + A *((PIV div (TP *DF))**P)");
 // FMAX is maximum infiltration

_spatial(REAL4, FILT);
calc(" FILT = mif(FMAX lt RainIntensity,FMAX,PondAreaFract*FMAX+(1.0-PondAreaFract)*RainIntensity)");
 // FILT is the water that infiltrates (mm/h)
 // if FMAX is less than the rain
 // rate, than FILT is equal to FMAX

_spatial(REAL4, FiltSize);
calc(" FiltSize = (FILT  * (DT div 3600.0))");
 // infiltration (mm) per time interval
calc(" FiltSize = min(WH,FiltSize)");
 // Maximum infiltration is equal to WH

calc(" WH -= FiltSize ");
 // WH is the water height after infiltration, in mm
calc(" WHWheelTrack = max(WHWheelTrack-FiltSize,0.0) ");
calc(" WHCompact = max(WHCompact-FiltSize,0.0) ");
calc(" WHGrass = max(WHGrass-FiltSize,0.0) ");
calc(" WHCrust = max(WHCrust-FiltSize,0.0) ");

// **** process drainage (percolation) ***********************

_spatial(REAL4, DR);
calc(" DR = mif(PIV lt GWC,FC * ((1.0-(PIV div GWC))),0)");
calc(" DR *= sqr(DR)"); /* this a optimization for DR = DR**3 */
 // drainage in mm/h

calc(" PIV = PIV-FiltSize+(DR * DT div 3600.0)");
 // the PIV change is calculated

calc(" PIV = max(PIV,0.0)");
 // PIV for next timeinterval
 // PIV cannot be less than zero
 // PIV cannot be less than zero
