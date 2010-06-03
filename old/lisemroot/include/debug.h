
#ifndef __STDIO
#include <stdio.h>
#define __STDIO
#endif


#ifdef DEBUG

/* avoid statement has no effect in
 * a constant precondition 
 * a source file that use PRECOND_CONST
 * must insert USE_PRECOND_CONST above
 * the USE header
 */
#define USE_PRECOND_CONST static int constCond_t_=1;


# include <assert.h>


# define DEBUGCOND(cond)	assert(cond)
# define IFDEBUG(action)     	action

# define POSTCOND(cond)		assert(cond)
# define PRECOND(cond)		assert(cond)

# define PRECOND_CONST(cond)	assert(constCond_t_ && (cond))



#else /* DEBUG not defined */

#define USE_PRECOND_CONST

# define DEBUGCOND(cond)
# define IFDEBUG(action)

# define POSTCOND(cond)
# define PRECOND(cond)

# define PRECOND_CONST(cond)	

#endif 

# define ARG_IS_USED(x)		PRECOND( ((x) == (x)) ? 1 : 0)
