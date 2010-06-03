#include "stddefx.h"
#include "misc.h"
#include "cps.h"

#include <stdlib.h>
#include <string.h>   /* strcmp */

#include "cps.h"

//#define DUMPNAMES (1)

MEM_HANDLE NULL_MAP;
static BOOL normalExit;

#ifdef COMPASDOS
  static int currDepth=0;
#else
  int currDepth;
#endif

static VARINFO *varList = NULL;
			/* points to the last record entered*/
			/* be aware that the way this list is traversed*/
			/* essential is to implement the local scope*/


typedef struct STATIC_LIST {
	VARINFO *address;
	struct STATIC_LIST *prev;
}  STATIC_LIST;


STATIC_LIST *staticList = NULL;
			/* every static spatial which is not in his*/
		        /* scope has in his variable a pointer to */
		 	/* his VARINFO with in it (in the address field)*/
			/* a MEMORYHANDLE*/
			/* these have to be freed at exit*/

/*function adds static to static spatial list*/
static void CpsAddStaticSpatial(VARINFO *vi)
{
	STATIC_LIST *p;
	p=(STATIC_LIST *)ChkMalloc(sizeof(STATIC_LIST));
	p->address=vi;
	p->prev=staticList;
	staticList=p;
}

/* add new element to varinfo-list */
void CpsDeclare(
	const char *name,    /*  name of the declared variable */
	      void *adr,     /* address of the declared variable
	                      * , which is a ptr to MEM_HANDLE union or
	                      *      to nonspatial value
	                      */
	const CSF_CR type, 
	const BOOL spatial, 
	const BOOL cpsStatic)
{
	VARINFO *vi;
	MEM_HANDLE *a = (MEM_HANDLE *)adr;
  FILE *dump;

	NULL_MAP.generic = NULL; /* just each time, who cares */


	if (cpsStatic && spatial && (a->staticKept != NULL))
	{
	    /*variable is static and spatial and it is not the first time used
	      in that case, memorylocation(address) contains a VARINFO * to
	      the VARINFO which contains data about this variable.
	      Now the VARINFO has to be put back into the VARINFO-list and
	      into memorylocation of (address) has to be put the memoryhandle
	      to the virtual location of the map. Also the current depth has
	      to be stored back
	    */

	    vi = a->staticKept;
	    /* add to varList
	     */
	    vi->prev = varList;
	    varList=vi;

	    /* in vi->adr the data ptr was stored
	     * 
	     */
	    a->generic = vi->addr.keepStaticData;
	    vi->addr.memH= a;
	    vi->declDepth = currDepth;
	    vi->staticInScope = TRUE;
	}
	else
	{ /* a true decl: */
	    vi =CHK_MALLOC_TYPE(VARINFO,1);
      vi->mapName[0] = '\0';
	    vi->prev      = varList;  /* add record at end of list*/
	    varList       = vi;

	    vi->type      = type;
	    vi->name      = name;
	    vi->declDepth = currDepth;
	    vi->spatial   = spatial;
	    vi->state     = CPS_EMPTY;
	    vi->cpsStatic = cpsStatic;
	    vi->staticInScope = TRUE;
	    vi->addr.generic   =  adr;

	    if (spatial) 
	    { 
	      /* put it's own address in the variable
	       * so it's known to defined (not NULL anymore) 
	       */
	      PRECOND(a->generic == NULL); /* first time called
	                                    * controlled by macro's in cps.h 
	                                    */
	      a->generic = a;
	    }

	    if (spatial && cpsStatic)
	    {
		      /*add to spatial to static list*/
          CpsAddStaticSpatial(vi);
	    }
	}
}

void CpsDeclareAndRead(
	const char *name,    /*  name of the declared variable */
	void *adr,     /* address of the declared variable
	                  which is a ptr to MEM_HANDLE union or
	                  to nonspatial value */
	CSF_CR type,   
	BOOL spatial,   
	BOOL cpsStatic, 
	const char *fileName) /* name of the file which will be read to initialize the param */
{
	MEM_HANDLE *a = (MEM_HANDLE *)adr;
	PRECOND(spatial);
	
	if ( a->generic == NULL)
	{
		CpsDeclare(name, adr, type, spatial, cpsStatic);
		/*after declare contents of address is changed*/
		/*not NULL anymore	*/

	  a = (MEM_HANDLE *)adr;
	  POSTCOND(a->generic != NULL);
		CpsRead(a,fileName);
	}
	else CpsDeclare(name, adr, type ,spatial,cpsStatic);
}



/*
	function frees all memory space used for static spatials
          etc.
 */

#ifdef COMPASDOS
  static void CpsCleanUp(void)
  {
	if (! normalExit)
		/* we should check the VARINFO list here
		 * beware that we mess up that list in staticList below
		 */
		return;
#else
  void CpsCleanUp(void)
  {
#endif
	while (staticList!=NULL)
	{
		VARINFO *v = staticList->address;
		STATIC_LIST *tmp;
		PRECOND(v->spatial);
		PRECOND(v->cpsStatic);
	  if (v->staticInScope)
	   		CpsFreeContents(v);
	  else 
	  	switch(v->state) {
	      case CPS_NORMAL : CpsFree2d((void **)v->addr.keepStaticData); break;
		    case CPS_COMPRESSED : CpsFreeCompressed((void *)v->addr.keepStaticData); break;
		    case CPS_EMPTY: break;
		  }
		Free(v);
		tmp=staticList;
		staticList=staticList->prev;
		Free(tmp);
	}

	FreeMask();
	POSTCOND(nrNormalMapsAllocated == 0);
	POSTCOND(nrCompressedMapsAllocated == 0);
}


void CpsEnterBlock(void)   /* called from MACRO begin*/
{
	if (currDepth==0)
	{
        	normalExit = FALSE;
        	exitOnError = 1;
		/* EFENCE CHOCKES ON IT
		 * Linux get inf. loop
		 */

/* VJ: initialize more here for windows version */
        staticList = NULL;
        varList = NULL;
        Merrno=NOERROR;
        ResetInitMaps();

       #ifdef COMPASDOS
   		  atexit(CpsCleanUp);
       #endif
	}
	currDepth++;	/* one nested-level deeper*/
}

void CpsLeaveBlock(void)   /* called form MACRO end */
{
	VARINFO *v = varList;
//for checking purposes	
#ifdef DUMPNAMES
  FILE *fout = fopen("names.txt","a");
  fprintf(fout,"Current depth = %i\n", currDepth);
#endif

	while( v != NULL && v->declDepth == currDepth)
	{	/* free all (but not statics) varinfo's from this block*/

		varList = varList->prev;
		if (v->spatial && v->cpsStatic)
		{
			/*put data ptr in VARINFO (in address) and
			  put pointer to VARINFO in variable where till
			  now the data ptr was stored
			*/
			MEM_HANDLE *tempAddress=v->addr.memH;
			v->addr.keepStaticData = tempAddress->generic;
	    	v->staticInScope = FALSE;
			tempAddress->staticKept=v;
		}
		else
		{
#ifdef DUMPNAMES
         fprintf(fout,"%s",v->name);
         if (v->spatial)
            fprintf(fout," spatial");
         else
            fprintf(fout," nonspatial");
         if (v->cpsStatic)
            fprintf(fout," static");
         fprintf(fout,"\n");
         fclose(fout);
         fout = fopen("names.txt","a");
#endif
			CpsFreeContents(v);
			Free(v);
      v = NULL; //!!!!!!!!!
		}
		v = varList;
	}
	currDepth--;  /* leave this level*/
   if (currDepth == 0)
    	normalExit = TRUE;

#ifdef DUMPNAMES
fclose(fout);
#endif
}

VARINFO *CpsFindByName(const char *name)
/*P name r- name of variable */
/*P RET     first structure with that name, NULL if none found */
{
	VARINFO *v = varList;

	while( v != NULL && strcmp(name, v->name))
		v = v->prev;
	return(v);
}


/* returns  first structure containing that address, aborts if none found */
VARINFO *CpsFindByAddress(
	const void *address) /* address  address of variable */
{
	VARINFO *vi = varList;

	while( vi != NULL && vi->addr.generic != address)
		vi = vi->prev;
	POSTCOND(vi != NULL);
	return(vi);
}

/* function returns structure containing that handle,
 * aborts if not found */
VARINFO *CpsFindByMemoryHandle(
	const MEM_HANDLE *handle)
/*P handle r- memory handle */
{
	VARINFO *v = varList;

	while (v != NULL && (v->addr.memH) != handle)
			v = v->prev;
	POSTCOND(v != NULL);
	return(v);
}


BOOL CpsIsSpatial(const VARINFO *vi)
/*P vi r- varinfo of variable which has to be tested if it is spatial*/
/*function returns TRUE is variable is spatial, FALSE otherwise*/
{
	return(vi->spatial);
}

BOOL CpsInitialized(const VARINFO *vi)
/*P vi r- varinfo of variable which has to be tested if it is initialized*/
/*function returns TRUE if variable was once initilized*/
/*FALSE if not*/
{
	PRECOND(vi->spatial);     /* must be spatial */

	return vi->state != CPS_EMPTY;
}


CSF_CR CpsCellRepr(const VARINFO *vi)
/*P vi r- varinfo of variable of which the cellrepresentation is wanted*/
{
	PRECOND(vi->type == CR_UINT1 || vi->type ==CR_INT2 || \
		vi->type==CR_REAL4 || vi->type==CR_REAL8);
	return(vi->type);
}

void *CpsDataAddress(const VARINFO *vi)
/*P vi r- varinfo of variable of which the address is wanted*/
/*function returns the address of the variable*/
{
	return(vi->addr.generic);
}
