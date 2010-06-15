#ifndef SWATRE_G_H
#define SWATRE_G_H

#include "csf.h"

//char *ErrorMessage;

typedef struct MEM_HANDLE {
	REAL8 **Data;
	int nrCols;
	int nrRows;
} MEM_HANDLE;


typedef struct SOIL_MODEL {
	struct PIXEL_INFO  *pixel;	 /* array of PIXEL_INFO  structs indexed as in compressed maps */
	double precision;
	double minDt;
	double *nodeData;
} SOIL_MODEL;

SOIL_MODEL *InitSwatre(
		MEM_HANDLE *profileMapHandle, /* r- profile id map */
		const  char *initHeadMaps,    /* init head maps path */
		double dtMin,   	      		/* minumum timestep, is also initial timestep */
		double precis,                /* precision factor to adapt timestep */
		double curtime);              /* current simulation time, same type as variable "index" in LISEM */

void CloseSwatre(
		SOIL_MODEL *s); /* soil model instance to be freed */

int ReadSwatreInput(
		const char *fileName,    /* profile definition file     */
		const char *tablePath);  /* directory containing tables */

void FreeSwatreInfo(void);

void SwatreStep(
		SOIL_MODEL *s,             /* rw soil models state */
		MEM_HANDLE *waterHeight, 	/* rw waterheight map  */
		double    lisemTimeStep,   /* size of time step used in LISEM */
		double     curtime);			/* current simulation time, same type as variable "index" in LISEM */

// average moisture content from surface to layernr, used in nutrients
void SwatreTheta(
		SOIL_MODEL *s,
		MEM_HANDLE *Theta,
		int layernr,
		int avg);

#endif /* SWATRE_G_H */
