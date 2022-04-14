/*
 * ERROR.C: routines to report errors, etc.
 */

#include "stdinc.h"
#include <stdarg.h>

/*
 * ERROR: scream and die quickly.
 */

void error(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);		/* invoke interface to printf       */
    fflush(stderr);			/* drain std error buffer 	    */
    va_end(ap);
    exit(1);				/* quit with error status	    */
}

/*
 * FATAL: scream and die a messy death.
 */

void fatal(string fmt, ...)
{
    va_list ap;
    stream errstr = fopen("fatal_error.log", "a");

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    if (errstr != NULL) {
	vfprintf(errstr, fmt, ap);
	fflush(errstr);
    }
    va_end(ap);
    abort();				/* quit, leave core image           */
}

/*
 * EPRINTF: scream, but don't die yet.
 */

void eprintf(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);		/* invoke interface to printf       */
    fflush(stderr);			/* drain std error buffer 	    */
    va_end(ap);
}

#ifdef TESTBED

main(int argc, string argv[])
{
    eprintf("warning: foo=%f\tbar=%d\tfum=\"%s\"\n", 3.1415, 32768, "waldo");
    error("error: foo=%f\tbar=%d\tfum=\"%s\"\n", 3.1415, 32768, "waldo");
}

#endif
