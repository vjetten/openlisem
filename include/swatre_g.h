#ifndef SWATRE_G_H
#define SWATRE_G_H

#include "csf.h"

typedef struct SOIL_MODEL {
	struct PIXEL_INFO  *pixel;	 // array of PIXEL_INFO  structs indexed as in compressed maps
	double precision;
	double minDt;
	double calibrationfactor;
	bool geometric;	
	bool swatreBottomClosed;
	long nrCells;
} SOIL_MODEL;

void FreeSwatreInfo(void);

#endif // SWATRE_G_H
