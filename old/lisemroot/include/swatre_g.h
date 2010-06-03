#ifndef SWATRE_G_H
#define SWATRE_G_H



#ifdef __cplusplus
 extern "C" {
#endif 



typedef struct SOIL_MODEL {
	struct PIXEL_INFO  *pixel;	 /* array of PIXEL_INFO  structs
	                                  * indexed as in compressed maps
	                                  */
	double precision;
	double minDt;
	double *nodeData;
} SOIL_MODEL;

void GetHoutputLoc(MEM_HANDLE **HHandle);

SOIL_MODEL *InitSwatre(
		MEM_HANDLE *profileMapHandle, /* r- profile id map */
		const  char *initHeadMaps,    /* init head name */
		MEM_HANDLE *OutHeadMapHandle,
          	const char *respath,  /* VJ map for output heads */	    
		double dtMin,   	      /* minumum timestep, is also initial timestep */
		double precis,                /* precision factor to adapt timestep */
                float  curtime);              /* current simulation time, same type as
		                                 variable "index" in LISEM */


void CloseSwatre(
	SOIL_MODEL *s); /* soil model instance to be freed */

void ReadSwatreInput(
	const char *fileName,    /* profile definition file     */
	const char *tablePath);  /* directory containing tables */

void FreeSwatreInfo(void);

void SwatreStep(
		SOIL_MODEL *s,                  /* rw soil models state */
		MEM_HANDLE *waterHeight, 	/* rw waterheight map  */
		MEM_HANDLE *InfilPot, 	/* rw waterheight map  */
       	MEM_HANDLE *outHeadMap,  /* VJ map for output heads */
      	const char *respath,  /* VJ map for output heads */
		double    lisemTimeStep,        /* size of time step used in LISEM */
        double     curtime,		/* current simulation time, same type as
						   variable "index" in LISEM */
        double calibrationfactor,
        BOOL geom);

void SwatreTheta(
		SOIL_MODEL *s,
		MEM_HANDLE *Theta,
		int layernr,
        int avg);

#ifdef __cplusplus
 }
#endif

#endif /* SWATRE_G_H */
