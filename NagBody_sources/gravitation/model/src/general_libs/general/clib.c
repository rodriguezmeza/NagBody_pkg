
#include "switchs.h"
#include "machines.h"

#include <math.h>
#include <time.h>

#include "stdinc.h"
#include "getparam.h"

#ifdef UNIX						        // 'UNIX' defined in machines.h
#include <sys/times.h>
#include <sys/param.h>
#else
#include <sys/types.h>
#endif

#include <sys/timeb.h>
#include <string.h>	

#ifdef VISUALC
#include <io.h>						    // dup is defined here
#endif
#include <unistd.h>                     // dup is defined here

#include <sys/stat.h>
#include <stdarg.h>


void *allocate(long int nb)
{
    void *mem;

    mem = calloc(nb, 1); 
    if (mem == NULL)
        error("allocate in %s: not enuf memory (%d bytes)\n",
              getargv0(), nb);
    return (mem);
}

real *AllocVecR(int size)				// Define array of size
{										// Index start in 1
	real *v;							// (as in fortran)
	v = (real *)malloc(size * sizeof(real));
	return (v - 1);
}

int *AllocVecI(int size)				// Define array of size
{										// Index start in 1
	int *v;								// (as in fortran)
	v = (int *)malloc(size * sizeof(int));
	return (v - 1);
}

int *AllocVecINormal(int size)			// Define array of size
{										// Index start in 0
	int *v;
	v = (int *)malloc(size * sizeof(int));
	return (v);
}

// Vector v have to be alloc with AllocVecR
void FreeVecR(real *v)
{
	free(v+1);
}

// Vector v have to be alloc with AllocVecI
void FreeVecI(int *v)
{
	free(v+1);
}

void FreeVecINormal(int *v)
{
	free(v);
}
 
#ifdef UNIX
double cputime(void)
{ 
    struct tms buffer;

    if (times(&buffer) == -1)
        error("cputime in %s: times() call failed\n", getargv0());
    return ((buffer.tms_utime + buffer.tms_stime) / (60.0 * HZ));
} 
#else
double cputime(void)
{
	time_t ltime;
	time(&ltime);
	return(((double)ltime)/((double) 60.0));
}
#endif


#ifdef MPI					            // MPI is defined in switchs.h

#ifdef WALLCLOCK
#include <mpi.h>
#endif

double second(void)
{
#ifdef WALLCLOCK
  return MPI_Wtime();
#else
  return ((double)clock())/CLOCKS_PER_SEC;
#endif  

}

#else

double second(void)
{
  return ((double)((unsigned int)clock()))/CLOCKS_PER_SEC;

}

#endif // ! MPI


double timediff(double t0, double t1)
{
  double dt;
  
  dt=t1-t0;

  if(dt<0)  // overflow has occured (for systems with 32bit tick counter)
    {
#ifdef WALLCLOCK
        dt = 0;
#else
      dt=t1 + pow(2,32)/CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}


void error(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);	
    fflush(stderr);
    va_end(ap);
    exit(1);
}

void endrun(int ierr)
{
  if(ierr)
    {
      fprintf(stdout,
		  "endrun called with an error level of %d\n\n\n", 
		  ierr);
      exit(1);
    }
  exit(0);
}


void eprintf(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    va_end(ap);
}


bool scanopt(string opt, string key)
{
    char *op, *kp;

    op = (char *) opt; 

// Debug:: Clean after solved
//    fprintf(stdout,"\n\n -> scanopt: %s %s %s \n\n", *opt, *key, *op);
//    fprintf(stdout,"\n\n -> scanopt: (%d) %d (%d) (%d)", opt, key, op, CHARNULL);
//    fprintf(stdout,"\n -> scanopt: (%s) %s (%s) (%s) (%s)", opt, key, op, CHARNULL, NULL);
    fprintf(stdout,"\n -> scanopt: (%s) %s (%s)", opt, key, op);
    fprintf(stdout,"\n -> scanopt: %d %d %d %d %d\n\n", *opt, *key, *op, NULL, CHARNULL);
    fflush(stdout);
//

// Original lines:
/*
    while (*op != (char)NULL) {
        kp = key;		
        while ((*op != ',' ? *op : (char) NULL) == *kp) {
            if (*kp++ == (char)NULL)
                return (TRUE);
            op++;
        }
        while (*op != (char)NULL && *op++ != ',')

	    continue;
    }
*/
//
    while (*op != CHARNULL) {
        kp = key;
        while ((*op != ',' ? *op : CHARNULL) == *kp) {
            if (*kp++ == CHARNULL)
                return (TRUE);
            op++;
        }
        while (*op != CHARNULL && *op++ != ',')

        continue;
    }
//
    return (FALSE);
}


stream stropen(string name, string mode)
{
    bool inflag;
    int fds;
    stream res;
    struct stat buf;

    inflag = streq(mode, "r");
    if (name[0] == '-') {                       
        if (streq(name, "-")) {
            fds = dup(fileno(inflag ? stdin : stdout));
            if (fds == -1)
                error("stropen in %s: cannot dup %s\n",
                      getargv0(), inflag ? "stdin" : "stdout");
        } else
            fds = atoi(&name[1]);
        res = fdopen(fds, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen in %s: cannot open f.d. %d for %s\n",
                  getargv0(), fds, inflag ? "input" : "output");
    } else {
        if (streq(mode, "w") && stat(name, &buf) == 0)
            error("stropen in %s: file \"%s\" already exists\n",
                  getargv0(), name);
        res = fopen(name, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen in %s: cannot open file \"%s\" for %s\n",
                  getargv0(), name, inflag ? "input" : "output");
    }
    return (res);
}

