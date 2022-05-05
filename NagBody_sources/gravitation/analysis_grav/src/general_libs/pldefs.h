
#ifndef __PLCDEMOS_H__
#define __PLCDEMOS_H__

#include "plConfig.h"
#include "plplot.h"
#include <math.h>
#include <string.h>
#include <ctype.h>

// En la version plplot-5.9.10 _c_plwid aparece como simbolo no definido... debe cambiarse a c_plwidth
// /Users/mar/NagBody_pkg/local/plplot/include/plplot/plplot.h:755:0: note: this is the location of the previous definition
// #define    plwid                    c_plwid
//
#define plwid plwidth           // Version plplot-5.9.10 (NEW 2013_11_27)
//
// Ver el manual plplot-5.9.10.pdf para ver como se usa plwidth (Y no está plwid)
//

// define PI if not defined by math.h

// Actually M_PI seems to be more widely used so we deprecate PI.
#ifndef PI
#define PI 3.1415926535897932384
#endif

#ifndef M_PI
#define M_PI 3.1415926535897932384
#endif

/* various utility macros */

#ifndef MAX
#define MAX(a,b)    (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)    (((a) < (b)) ? (a) : (b))
#endif

#ifndef ROUND
#define ROUND(a)    (PLINT)((a)<0. ? ((a)-.5) : ((a)+.5))
#endif


// Declarations for save string functions

#ifdef PL_HAVE_SNPRINTF
// In case only _snprintf is declared (as for Visual C++ and
// Borland compiler toolset) we redefine the function names
#ifdef _PL_HAVE_SNPRINTF
#define snprintf    _snprintf
#define snscanf     _snscanf
#endif // _PL_HAVE_SNPRINTF
#else // !PL_HAVE_SNPRINTF
// declare dummy functions which just call the unsafe
// functions ignoring the size of the string
int plsnprintf( char *buffer, int n, const char *format, ... );
int plsnscanf( const char *buffer, int n, const char *format, ... );
#define snprintf    plsnprintf
#define snscanf     plsnscanf
#endif // PL_HAVE_SNPRINTF

// Add in missing isnan definition if required
#if defined ( PL__HAVE_ISNAN )
#  define isnan    _isnan
#  if defined ( _MSC_VER )
#    include <float.h>
#  endif
#endif

#if !defined ( PL_HAVE_ISNAN )
#  define isnan( x )    ( ( x ) != ( x ) )
#endif

#endif	/* __PLCDEMOS_H__ */
