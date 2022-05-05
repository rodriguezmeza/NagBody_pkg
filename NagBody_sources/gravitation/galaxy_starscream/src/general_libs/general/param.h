
#if !defined (HZ)
#define HZ        100
#endif


#ifndef __dj_include_sys_param_h_
#define __dj_include_sys_param_h_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __dj_ENFORCE_ANSI_FREESTANDING

#ifndef __STRICT_ANSI__

#ifndef _POSIX_SOURCE

#include <limits.h>

#define MAXPATHLEN	PATH_MAX
#define MAXGETHOSTNAME	128

#endif 
#endif 
#endif 

#ifndef __dj_ENFORCE_FUNCTION_CALLS
#endif 

#ifdef __cplusplus
}
#endif

#endif 
