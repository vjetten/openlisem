#include "stddefx.h"

#include "cps.h"
#include "misc.h"
#include "string.h"

/* read a map
 */
void **CpsOpenAndReadAs2d(
	const char *fileName, /*  name of file which has to be read */
	CSF_CR inputType)     /*  is the type which should be read*/
{
	MAP *m;
	void *linBuf;
	size_t nrCols,nrRows,nrCells;
	CSF_CR type;

    	m = Mopen(fileName, M_READ);
    	if (m == NULL)
        {
	    		MperrorExit(fileName,1);
            return NULL;
        }
        InitMaps(m);  /* init Map specifications */

        if (!RtestSpecs(m))
		Error("location attributes of map: %s\n"
		      "differ from location attributes of earlier read map(s)\n"
		      ,fileName);

	type = inputType;
	if (type == CR_LDD)
		type = CR_UINT1;

    	if (RuseAs(m, type)) /* csf2 -> RuseAs */
		MperrorExit(fileName,1);

    	nrRows = RgiveNrRows();
    	nrCols = RgiveNrCols();
    	nrCells = nrRows *nrCols;

	linBuf=Rmalloc(m,nrCells);

        if ((RgetSomeCells(m,0,nrCells, linBuf)) != nrCells)
		MperrorExit(fileName,1);
        Mclose(m);

        if (inputType == CR_LDD)
        {
    	size_t i;
    	size_t p=0;
    	UINT1 *v = (UINT1 *)linBuf;
    	for(i=0; i < nrCells; i++)
    	{
    		if (v[i] != MV_UINT1 && (v[i] < 1 || v[i] > 9))
    			Error("LDD type map '%s' contains values outside [1,9] range",
    				fileName);
    		if (v[i] == 5)
    			p++;
    	}
    	if (p == 0)
    		Error("LDD type map '%s' has no pit (or outflow point, value 5)",
    			fileName);
//    	if (p > 1)
//    		Error("LDD type map '%s' has multiple (%u) pits (or outflow point, value 5)",
//    			fileName,p);
//VJ: is allowed now
        }
        linBuf = ChkRealloc(linBuf, nrRows*nrCols*CELLSIZE(type));
	IFDEBUG(nrNormalMapsAllocated++);
        return MallocIndex2d(nrRows,nrCols,CELLSIZE(type),linBuf);
}

void CpsRead(
	MEM_HANDLE *variable, /* address of variable in which the data will be stored*/
	const char *fileName)
{
	VARINFO *vi;
	void **vMap;

	vi=CpsFindByAddress(variable);

	PRECOND(vi->spatial);  /*variable must be spatial*/

	strcpy(vi->mapName,fileName);
	if (vi->state != CPS_EMPTY)  /*variable has already memory space*/
	    CpsFreeContents(vi);

	if (vi->type != CR_REAL4 && vi->type != CR_UINT1 && vi->type != CR_LDD && vi->type != CR_INT4)
		Error("Illegal type at definition of '%s' ('%s)\n"
		      " should be UINT1, REAL4 or LDD", vi->name, fileName);
	vMap=CpsOpenAndReadAs2d(fileName,vi->type);
	if (vi->type == CR_LDD)
		vi->type = CR_UINT1;

	if (vMap==NULL) 
		Error("failure at Read\n");
	else
	{
		vi->state=CPS_NORMAL;
		PRECOND(variable == vi->addr.memH);
		vi->addr.memH->voidSpatialNorm=vMap;
	}
}


void CpsWrite(
	MEM_HANDLE *vMap, /* the map which should be stored*/
	const char *filename) /* name off the file in which it must be stored*/
{
	VARINFO *vi;
	MAP *mapHandle;
	void *mapData; /* linear buffer */
	size_t nrCells = RgiveNrCells();
	size_t nrRows = RgiveNrRows();
	size_t nrCols = RgiveNrCols();
	CSF_VS vs;

	vi=CpsFindByMemoryHandle(vMap);
	PRECOND(vi->state != CPS_EMPTY); /* CW RUNTIME ERROR */
	PRECOND(vi->spatial); /* CW RUNTIME ERROR */

	if(vi->type == CR_UINT1) 
		vs = VS_NOMINAL;
	else {
	 PRECOND(vi->type == CR_REAL4);
	 vs = VS_SCALAR;
	}

	mapHandle=Rcreate(filename,nrRows,nrCols,vi->type, vs,
	 MgiveProjection(), RgiveX0(), RgiveY0(), 0.0 /* angle */,
	     RgiveCellSizeX());
	if (mapHandle==NULL)
		MperrorExit(filename,1);


	mapData = CpsNormal(vi)[0];
	if (RputSomeCells(mapHandle,0,nrCells,mapData) != nrCells)
		MperrorExit(filename,1);
	Mclose(mapHandle);
}
