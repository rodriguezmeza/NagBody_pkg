
// CPUTIME: returns elapsed CPU in seconds.  This version can be called
// from fortran routines in double precision mode.


#ifdef UNIX
#include <sys/times.h>
#include <sys/param.h>
#else
#include <sys/types.h>
#include <time.h>
#endif

#ifndef HZ
#  include <time.h>
#  define HZ CLK_TCK
#endif

extern void my_cputime(double *t);				/* Interface with f90 programs */
double c_cputime(void);

extern void my_cputime(t)
double *t;
{
	*t = c_cputime();
}

#ifdef UNIX
double c_cputime_()
{
    struct tms buffer;
        
    if (times(&buffer) == -1) {
	printf("times() call failed\n");
	exit(1);
    }
    return( buffer.tms_utime / (60.0 * HZ) );	/* in minutes */
}
#else
double c_cputime(void)
{
	time_t ltime;								/* Gives the real time			*/
	time(&ltime);
	return(((double)ltime)/((double) 60.0));
}
#endif
