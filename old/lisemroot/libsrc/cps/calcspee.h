#ifndef __CALCSPEE
#define __CALCSPEE

void CalcSpeed(REAL4 **speedMap,REAL4 **slopeMap,REAL4 **manningMap,REAL4 **detentionMap,REAL4 **widthMap);

/*P speedMap     -w map with speed for each element*/
/*P slopeMap     r- map with slope in each element */
/*P manningMap   r- map with manning's n for each element */
/*P detentionMap r- map with the detention in each map */
/*P widthMap	 r- map with width of elements (used for hydraulische r)*/
/*P		    NULL if map is not given, in that case RgiveCellSizeX*/
/*P		    will be used*/

#endif /*__CALCSPEE*/
