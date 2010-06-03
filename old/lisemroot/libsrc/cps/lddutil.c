#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "lddutil.h"

/* static, put one in each module */
const LDD_DIR LddData[10] = {{ 0,   0},
			{-1,   1},    /* 1 */
			{ 0,   1},    /* 2 */
			{ 1,   1},    /* 3 */
			{-1,   0},    /* 4 */
			{ 0,   0},    /* 5 */
			{ 1,   0},    /* 6 */
			{-1,  -1},    /* 7 */
			{ 0,  -1},    /* 8 */
			{ 1,  -1}     /* 9 */
};
