/*==============================================================================
	MODULE: cputime					[General_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: returns elapsed CPU in seconds.  This version can be 
			   called from fortran routines in double precision mode.
	Language: C
	Use: '#include "...."
	Use in routines and functions: model (main)
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 1999-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/


#include "machines.h"

#define CPUTIME_OPT1				/* CPU time option 1.			  */
#undef CPUTIME_OPT1					/* CPU time option 2.			  */
									/* Option 2 gives results		  */
									/* closer to cpu times given by   */
									/* cpu_time routine of Visual	  */
									/* f90.							  */
#define IN_MINUTES
#undef IN_MINUTES


#ifndef IN_MINUTES
#define IN_SECONDS
#endif

#ifdef UNIX
#include <sys/times.h>
#include <sys/param.h>
#if !defined(CPUTIME_OPT)
#include <time.h>					/* Here is CLOCKS_PER_SEC		  */
									/* defined						  */
#endif
#else
#include <sys/types.h>				/* En el original est'a llamado 
									   en UNIX (verificarlo). Ver 
									   rutina SECOND.C barnesf77	  */
#include <time.h>
#endif

#ifndef HZ
#  include <time.h>
#  define HZ CLK_TCK
#endif

#ifdef MEDUSA
extern void my_cputime_(double *t);	/* Interface with f90 programs	  */
#endif
#ifdef PC169
extern void _my_cputime(double *t);	/* Interface with f90 programs	  */
#endif
#ifdef DEMI1
extern void my_cputime(double *t);	/* Interface with f90 programs    */
#endif
#ifdef VISUALC
extern void my_cputime(double *t);	/* Interface with f90 programs    */
#endif

#ifdef MEDUSA
double c_cputime_(void);
#endif
#ifdef PC169
double c_cputime_(void);
#endif
#ifdef DEMI1
double c_cputime(void);
#endif
#ifdef VISUALC
double c_cputime(void);
#endif

#ifdef DEMI1
extern void my_cputime(t)
#endif
#ifdef VISUALC
extern void my_cputime(t)
#endif
#ifdef MEDUSA
extern void my_cputime_(t)
#endif
#ifdef PC169
extern void _my_cputime(t)
#endif
double *t;
{
#ifdef VISUALC
*t = c_cputime();
#else
#ifdef MEDUSA
*t = c_cputime_();
#endif
#ifdef DEMI1
*t = c_cputime();
#endif
#ifdef PC169
*t = c_cputime_();
#endif
#endif
}

#ifdef UNIX
#ifdef DEMI1
double c_cputime()
#else
double c_cputime_()
#endif
{
#ifdef CPUTIME_OPT1
    struct tms buffer;
        
    if (times(&buffer) == -1) {
	printf("times() call failed\n");
	exit(1);
    }
#ifdef IN_MINUTES
    return( buffer.tms_utime / (60.0 * HZ) );	/* in minutes		  */
#else
    return( buffer.tms_utime / (HZ) );			/* in seconds		  */
#endif

/* 
	Other option as found in barnesc13. 
	The difference is the additional term in the
	return line.
*/
/* 
    if (times(&buffer) == -1)
        error("cputime in %s: times() call failed\n", getargv0());
    return ((buffer.tms_utime + buffer.tms_stime) / (60.0 * HZ));
*/
#else											/* CPU time option 2  */
#ifdef IN_MINUTES
	return(((double)((unsigned int)clock()))/(60.0 * CLOCKS_PER_SEC) );
#else
	return(((double)((unsigned int)clock()))/(CLOCKS_PER_SEC) );  /* in seconds */
#endif
#endif
}
#else
#ifdef CPUTIME_OPT1
double c_cputime(void)
{
	time_t ltime;							/* Gives the real time	  */
	time(&ltime);
#ifdef IN_MINUTES
	return(((double)ltime)/((double) 60.0));		/* in minutes	  */
#else
	return((double)ltime);							/* in seconds	  */
#endif
}
#else
double c_cputime(void)				/* As defined in the original	  */
									/* gnbsph						  */
{									/* CPU time option 2			  */
#ifdef IN_MINUTES
  return(((double)((unsigned int)clock()))/(60.0 * CLOCKS_PER_SEC) );
#else
  return(((double)((unsigned int)clock()))/(CLOCKS_PER_SEC) );	/* in seconds */
#endif
}
#endif

#endif

