#ifndef __CPS
#define __CPS



#ifdef __cplusplus
 extern "C" {
#endif 



/* stdlibx.c */
void StrToUpper(char *string);
 /*converts given string to only uppercase letters*/
char *CatPath(const char *fileName,const char *path);
char *MakeFileName(char *fileName,float minutes);
 /*P fileName r- is the base filename (first 5 characters of filename) */
 /*P minutes r- is the rest of the filename (extention is fraction of min)*/
 /* function adds number to filename*/
 /* returned filename exists until!! next call to MakeFileName*/
/* stdlibx.c */
extern char *MakePCRFileName(char *iname, int atTime);

#include "csf.h"

#define CPS_SPATIAL TRUE
#define CPS_NONSPATIAL FALSE
#define CPS_STATIC TRUE
#define CPS_NONSTATIC FALSE

extern int nrNormalMapsAllocated;
extern int nrCompressedMapsAllocated;

typedef enum CPS_MAP_FORM {
  CPS_EMPTY=0,
  CPS_NORMAL=1,
  CPS_COMPRESSED=2
} CPS_MAP_FORM;

typedef union VAR_ADR { /* address of variable */
	/* non-spatial: address of variable to value if nonspatial*/
	void *voidNonSpatial;
	UINT1 *uint1NonSpatial;
	 INT2 *int2NonSpatial;
	REAL4 *real4NonSpatial;
	REAL8 *real8NonSpatial;
	/* spatial: address of variable */
	union MEM_HANDLE *memH;
	/* static out of scope */
	void *keepStaticData;
	/* generic */
	void *generic;
} VAR_ADR;


/*********************************************************************/
/*                   data structure of VARINFO-list                  */
/*********************************************************************/

#define CR_LDD       CR_UNDEFINED
typedef struct VARINFO {
	const char *name;       /* name of variable
	                         * ptr to string constant
	                         */
	CSF_CR type;           	/* type constant  CR_UINT1 or CR_REAL4 
	                         * CR_LDD also possible as intermediate state
	                         */
	BOOL spatial;        	/* spatial or nonspatial*/
	int  declDepth;		/* the depth at wich this variable is declared*/
	CPS_MAP_FORM state;		/* the state in which a spatial is*/
	BOOL cpsStatic;		/* variable is CPS_STATIC or CPS_NONSTATIC*/
        BOOL staticInScope;     /* TRUE is static is in active scope 
                                 * FALSE otherwise
                                 */
	char mapName[128];     /* empty string if intermediate variable 
	                        * name if initialized with a map
	                        */
	VAR_ADR addr;
	struct VARINFO *prev;   /* ptr to previous VARINFO or NULL if this is*/
				/* the first one*/
} VARINFO;

typedef union MEM_HANDLE {
	/* spatial normal : ptr to 2-dim array */
	void **voidSpatialNorm;
	UINT1 **uint1SpatialNorm;
	REAL4 **real4SpatialNorm;
	REAL8 **real8SpatialNorm;
	/* spatial compressed : ptr to 1-dim array */
	void *voidSpatialComp;
	UINT1 *uint1SpatialComp;
	REAL4 *real4SpatialComp;
	REAL8 *real8SpatialComp;
	/* if static and not in scope this is a varinfo */
	VARINFO *staticKept;
	/* the generic ptr: */
	void *generic;
} MEM_HANDLE;

/* typedef full const qual maps */
typedef const REAL4 * const *CONST_REAL4_MAP;
typedef const UINT1 * const *CONST_UINT1_MAP;

/* TYPE STUFF for timeseries
 */
typedef enum TSER_TYPE {
	TSER_QUANTITY  =01,
	TSER_INTENSITY =02} TSER_TYPE;

typedef enum TSER_CUMM {
	TSER_UNDEF     =00,
	TSER_NORMAL    =01,
        TSER_CUMULATIEF=02} TSER_CUMM;

typedef struct COORDINATES
	{
		size_t rowNr;
		size_t colNr;
	}   COORDINATES;

typedef struct TIMESERIE
	{
		size_t numberStations;
		TSER_TYPE type;
		TSER_CUMM cumulatief;
                size_t lines;    /* ==  the number of timeintervals */
		REAL4 *time;
		REAL4 **timeSerie;
		COORDINATES **station;

	} TIMESERIE;

/*********************************************************************/
/*                      VARINFO - list operations                    */ 
/*********************************************************************/

void CpsRangeTest( void *a, const char *ops, double low, double high,
	const char *paramDescr, const char *fileName, int  lineNr);

#define R_GE_LE "[]"
#define R_GE_LT "[>"
#define R_GT_LE "<]"
#define R_GT_LT "<>"
#define R_GT_LE "<]"
#define R_GT    "<_"
#define R_GE    "[_"
#define R_LT    "_>"
#define R_LE    "_]"
#define R_DUMMY 0

#define rangetest(m,ops,low,high,paramDescr) \
  CpsRangeTest(&m,ops,low,high,paramDescr, __FILE__, __LINE__)

void CpsDeclare(const char *name, void *memHandle, const CSF_CR type, const BOOL spatial, const BOOL cpsStatic);

void CpsDeclareAndRead(const char *name, void *memHndl, CSF_CR type, BOOL spatial, BOOL cpsStatic, 
			const char *fileName);


VARINFO *CpsFindByName(const char *name);
/*P name r- name of variable */
/*P RET     first structure with that name, NULL if none found */

VARINFO *CpsFindByAddress(const void *address);
/*P address r- address of variable */
/*P RET     first structure containing that address, NULL if none found */

VARINFO *CpsFindByMemoryHandle(const MEM_HANDLE *handle);
/*P handle r- memory handle */
/*function returns structure containing that handle, NULL if not found */

void CpsEnterBlock(void);   /* called from MACRO begin*/

void CpsLeaveBlock(void);   /* called form MACRO end*/

void CpsRead(MEM_HANDLE *variable,const char *fileName);
void CpsWrite(MEM_HANDLE *vMap,const char *filename);

/*******************************************************************/
/*                       DEFINE OF MACROS                          */
/*******************************************************************/



#define begin		{ CpsEnterBlock();
#define begin_model	{ CpsEnterBlock();
#define end		CpsLeaveBlock(); }
#define end_model	CpsLeaveBlock(); exit(0); return 0; }

#define _spatial(type, varName) \
  MEM_HANDLE varName={NULL}; \
  CpsDeclare(#varName, &varName, CR_##type, CPS_SPATIAL, CPS_NONSTATIC)

#define _spatial_input(type, varName, mapName) \
  MEM_HANDLE varName={NULL}; \
  CpsDeclare(#varName, &varName, CR_##type, CPS_SPATIAL, CPS_NONSTATIC); \
  readpath(varName, mapName, PATH)

#define _static_spatial(type, varName) \
 static MEM_HANDLE varName={NULL}; \
 CpsDeclare(#varName, &varName, CR_##type, CPS_SPATIAL, CPS_STATIC)

#define _static_spatial_init(type, varName, mapName)   \
 static MEM_HANDLE varName={NULL}; \
 CpsDeclare(#varName, &varName, CR_##type, CPS_SPATIAL, CPS_STATIC); \
 readpath(varName, mapName, PATH)


#define _nonspatial(type, varName) \
	  type varName;CpsDeclare(#varName, &varName, CR_##type, \
	                                      CPS_NONSPATIAL, CPS_NONSTATIC)

#define _static_nonspatial(type, varName) \
	  static type varName;CpsDeclare(#varName, &varName, CR_##type, \
	                                          CPS_NONSPATIAL, CPS_STATIC)
	  
#define _static_nonspatial_init(type, varName, value) \
	  static type varName=value;CpsDeclare(#varName, &varName, CR_##type, \
	                                      CPS_NONSPATIAL, CPS_NONSTATIC)


#define read(varName, fileName) CpsRead(&varName,fileName)
	  /*can't use only varName without & and then find the address in
	    the VARINFO list by FindByMemoryHandle because the contents can
	    be NULL or a variable which is nonspatial can be passed (can't
	    find memoryhandle in that case)
	  */

#define readpath(varName, fileName, path) CpsRead(&varName,CatPath(fileName,path))


#define write(varName, fileName) CpsWrite(&varName,fileName)

/* disperse is very old (written by Harm, to transport to multiple
 * cells on base of an aspect map
 */
REAL4 CpsDisperse(MEM_HANDLE *cMap,MEM_HANDLE *roMap,MEM_HANDLE *aspectMap);
#define disperse(changeMap,runOverMap,aspectMap) \
		CpsDisperse(&changeMap,&runOverMap,&aspectMap)


void CpsTimeIn(MEM_HANDLE *vResultMap,MEM_HANDLE *vIdMap,const char *serieFileName,double intervalEnd, double intervalLength);
#define timeinput(newMap,idMap,fileName,intervalEnd,intervalLenght) \
		CpsTimeIn(&newMap,&idMap,fileName,(double)intervalEnd,(double)intervalLenght)
	 /*idMap has to be UINT1 */

void CpsTimeOutput(const char *fileName,MEM_HANDLE *vIdMap,MEM_HANDLE *vDataMap,double time);
#define timeoutput(fileName,idMap,dataMap,time) \
		CpsTimeOutput(fileName,&idMap,&dataMap,time)
	 /*idMap has to be UINT1*/

/* should bedefined in cpslang.c: */
extern MEM_HANDLE NULL_MAP;

/* kinemati.c */
void Kinematic(
	const char *srcFile, int srcLineNr,
	REAL4 *Qout, REAL4 *Sout,
	MEM_HANDLE *Q1, MEM_HANDLE *Q0, MEM_HANDLE *q, 
	MEM_HANDLE *S1, MEM_HANDLE *S0, MEM_HANDLE *s, 
	MEM_HANDLE *ldd, MEM_HANDLE *alpha, double beta, double deltaT, double deltaX, double epsilon);
#define kinematic(Qout,Sout,Q1,Q0,q,S1,S0,s,ldd,alpha,beta,deltaT,deltaX,epsilon) \
		Kinematic( __FILE__, __LINE__, (REAL4 *)Qout, (REAL4 *)Sout, \
		           &Q1,&Q0, \
		             (((q).generic == NULL) ? &NULL_MAP : &q),\
		            &S1, &S0,\
		             (((s).generic == NULL) ? &NULL_MAP : &s),\
			    &ldd,&alpha, \
			      beta,deltaT,deltaX, epsilon)
	/*ldd map has to be UINT1*/

/* upstream.c */
extern void Upstream(const char *srcFile, int srcLineNr, MEM_HANDLE *amountOutMap, MEM_HANDLE *amountInMap, MEM_HANDLE *lddMap);
extern int Downstream(const char *srcFile, int srcLineNr, MEM_HANDLE *amountOutMap, MEM_HANDLE *amountInMap, MEM_HANDLE *lddMap);

#define upstream(amountOut, amountIn, ldd) \
	Upstream(__FILE__,__LINE__,&amountOut, &amountIn, &ldd)
#define downstream(amountOut, amountIn, ldd) \
	Downstream(__FILE__,__LINE__,&amountOut, &amountIn, &ldd)


void CalcHeight(MEM_HANDLE *Q,MEM_HANDLE *b, MEM_HANDLE *h, MEM_HANDLE *z,
		MEM_HANDLE *n, MEM_HANDLE *s, double epsilon);
/*P Q r- is the Q wherefor h has to be appoximated*/
/*P b r- is the cell width */
/*P h rw is the first approximated h */
/*P z r- see Applied Hydrology p. 162 */
/*P n r- is manning's n*/
/*P s r- is the slope*/
/*P epsilon r- is the difference allowed between the appoxed Q and real Q*/
#define calcheight(Q,b,h,z, n,s,epsilon) \
		CalcHeight(&Q,&b,&h,&z,&n,&s,epsilon)

/* sprdldd.c */
extern void LddBuffer(const char *srcFile, int srcLineNr, 
	MEM_HANDLE *h_idRes, MEM_HANDLE *h_points, MEM_HANDLE *h_ldd, double bufSize);
#define lddbuffer(outId,inId,ldd,bufSize) \
                LddBuffer(__FILE__,__LINE__,&outId,&inId,&ldd,bufSize) 

/* NOTE: initswatre changed by VJ */
#define initswatre(idMap, hMap, initheads, outheads, dtMin, precision, curTime)\
	 InitSwatre(&idMap, &hMap, initheads, outheads, dtMin, precision, curTime)
#define initsoil(soilDescrFile, tableDir)\
	 ReadSwatreInput(soilDescrFile, tableDir)

#define swatrestep(s, waterHeightMap,   outheads, timeStep, curTime)\
     SwatreStep(s, &waterHeightMap,   outheads, timeStep, curTime)
/****************************************************************************/
/*                         used in calc                                    */
/****************************************************************************/
void *CpsDataAddress(const VARINFO *vi);
BOOL CpsIsSpatial(const VARINFO *vi);
CSF_CR CpsCellRepr(const VARINFO *vi);
BOOL CpsInitialized(const VARINFO *vi);
void CpsInitWithDataCompressed(VARINFO *vi, MEM_HANDLE data);

/****************************************************************************/
/*                        rest of include files                             */
/****************************************************************************/

/* function prototypes */

/* cpsio.c */
extern void **CpsOpenAndReadAs2d(const char *fileName, CSF_CR inputType);

/* initmaps.c */
void InitMaps(MAP *idMapHandle);
void ResetInitMaps(void);
BOOL MapSpecsInit(void);
CSF_VS RgiveValueScale(void);
REAL8 RgiveX0(void);
REAL8 RgiveY0(void);
size_t RgiveNrRows(void);
size_t RgiveNrCells(void);
size_t RgiveNrCols(void);
REAL8 RgiveCellSizeX(void);
REAL8 RgiveCellSizeY(void);
CSF_PT MgiveProjection(void);
BOOL RtestSpecs(MAP *mapHandle);

/* compress.c */
extern void **CpsNormal(VARINFO *vi);
extern void *CpsNorm2Comp(void **norm, CSF_CR cr);
extern void *CpsCompressed(VARINFO *vi);
extern void *CpsCompressedHandle(MEM_HANDLE *m, CSF_CR cr);
extern void *CpsNormalHandle(MEM_HANDLE *m, CSF_CR cr);
#ifdef DEVELOP_DEBUG
 BOOL ChkCompPtr(const void *p);
#else
# define ChkCompPtr(x)	(1)
#endif

/* cpsmem.c */
void *CpsMallocCompressed(size_t elementSize);
void CpsFreeCompressed(void *hhMap);
void **CpsMalloc2d(CSF_CR cr);
void CpsFree2d(void **data);
void CpsFreeContents(VARINFO *v);


/* mask.c */
UINT1 *GiveMask(void);
 /*function returns the mask*/
size_t GiveNrDefinedCells(void);
size_t GiveNrMask(void);
void FreeMask(void);
 /*function free's all memory allocated to hold the mask-info
  called at exit*/
void InitMask(char *name);
 /*P name r- is the name of the file which will be read*/


/* cpswarn.c */
extern void LisemWarning(const char *fmt, ...);
extern void Report(const char *fmt, ...);
extern BOOL silentMode;

/* cpsparam.c */
extern void OpenParamFile(
BOOL *writeEachTime,
 int argc, char *argv[], const char *modelName, const char *version);
extern void ReadStringParamFile(char *str);
extern void ReadNumberParamFile(REAL4 *val);

#ifdef __cplusplus
 }
#endif 

#endif /*__CPS*/

