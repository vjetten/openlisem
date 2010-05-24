/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/


//---------------------------------------------------------------------------
#ifndef CsfMapH
#define CsfMapH
//---------------------------------------------------------------------------

#include "csf.h"



#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))


//---------------------------------------------------------------------------

class cTMap
{
protected:

public:
    CSF_RASTER_HEADER MH;
    UINT2 projection;
    REAL4 **Data;

    QString MapName;
    QString PathName;
    int nrRows, nrCols;
    bool Created;
    
    void KillMap();
    void GetMapHeader(QString Name);
    void CreateMap(QString Name);
    void _MakeMap(cTMap *dup, REAL4 value);
    void WriteMap(QString Name);
    void WriteMapSeries(QString Dir, QString Name, int count);
    bool LoadFromFile();
    void ResetMinMax(void);


    cTMap();
    ~cTMap();
};

#endif

/*
typedef struct CSF_RASTER_HEADER
{
	 UINT2    valueScale;
	 UINT2    cellRepr;
	 CSF_VAR_TYPE minVal;
	 CSF_VAR_TYPE maxVal;
	 REAL8    xUL;
	 REAL8    yUL;
	 UINT4    nrRows;
	 UINT4    nrCols;
	 REAL8    cellSizeX;
	 REAL8    cellSizeY;
	 REAL8    angle;
} CSF_RASTER_HEADER;
*/

