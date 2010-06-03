typedef struct Liststruct {
	int rowNr;
	int colNr;
	struct Liststruct *prev;
	}  Liststruct;

/* static, put one in each module */
typedef struct LDD_DIR {
	int deltaX;
	int deltaY;
} LDD_DIR;

extern const LDD_DIR LddData[10];


typedef struct LDD_SORT {
	int r;
	int c;
   int d;
} LDD_SORT;


/*
  local drain direction maps have values for directions as following:
    7  8  9
     \ | /
   4 - 5 - 6
     / | \
    1  2  3
*/

/* determine if (rFrom,cFrom) flows to (rTo, cTo)  */
#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
		( ldd[rFrom][cFrom]!=MV_UINT1 && rFrom >= 0 && cFrom >= 0 &&\
		  rFrom+LddData[ldd[rFrom][cFrom]].deltaY==rTo &&\
		  cFrom+LddData[ldd[rFrom][cFrom]].deltaX==cTo )

/* lddutil.c */
extern void FindOutflowPoint(CONST_UINT1_MAP ldd, int *rowNr, int *colNr);
