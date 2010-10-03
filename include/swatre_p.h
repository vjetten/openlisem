/*************************************************************/
/* private include file of swatre implementation: swatre_p.h */
/*************************************************************/
#include "soillut.h"
//#define min(a,b) ((a > b)? b : a)
//#define max(a,b) ((a > b)? a : b)

typedef struct ZONE   {
	/* structure explaining how the soil subdivided */
	int  nrNodes;     	/* nr. of nodes (equals nr. of compartments), arrays with nrNodes elements: */
	double *dz;       	/* compartment size (cm.) used as negative in [SWATRE]*/
	double *z;        	/* position of nodal point relative to top soil (Cm) e.g -21 means 21 cm below top of soil */
	double *endComp;   	/* end of compartment i only used when reading the profiles, arrays with nrNodes+1 elements: */
	double *disnod;    	/* distance between nodal points, 0 is between top-profile and first nodal point
								  last is between bottom and last nodal point */
} ZONE;

/* change this structure if we add VanGenughten eqs. */
typedef struct HORIZON {
	char *name;  	/* name of horizon: filename of lut-table */
	LUT  *lut;     /* lut of theta, h, k, dmch, dmcc */
} HORIZON;

typedef struct PROFILE {
	int            profileId; 	/* number identifying this profile  >= 0 */
	const ZONE     *zone; 		/* array with zone.nrNodes elements: */
	const HORIZON  **horizon; 	/* ptr to horizon information this node belongs to */
} PROFILE;

#define MAX_NODES            20
#define MAX_NODES_P          (MAX_NODES+3)

/* maximum amount of ponding that is regarded as no ponding (0) */
#define POND_EPS             (1.0E-6)
/* maximum amount of time for which it is not worth doing an iteration */
#define TIME_EPS             (1.0E-6)

typedef double NODE_ARRAY[MAX_NODES_P];

typedef struct PIXEL_INFO {
	double        *h;          /* array of MAX_NODES nodes with matrix head */
	double        currDt;      /* current size of SWATRE timestep */
	const PROFILE *profile;    /* profile this pixel belongs to */
	//int           dumpHid;   /* if 0 then no head output else write to file amed Hx where x is dumpH value */
} PIXEL_INFO;

//#define PIXEL_IS_DEFINED(p)	(p.h != NULL)


#define NrNodes(profile)        (profile->zone->nrNodes)
#define Dz(profile)             (profile->zone->dz)
#define DistNode(profile)       (profile->zone->disnod)
#define Horizon(profile, node)  (profile->horizon[node])
/* sizeof intermediate arrays is fixed to optimize computation: */
