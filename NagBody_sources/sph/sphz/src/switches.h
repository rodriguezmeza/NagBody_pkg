//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)

#ifndef _switches_h
#define _switches_h


#if defined(ENTROPY)
#  define THERMOPT		";Entropy formulation,",
#  if defined(ADIABATIC)
#    define GASOPT		";Adiabatic gas (constant entf),",
#  else
#    define GASOPT		";Non-adiabatic gas (shock heating),",
#  endif
#else
#  define THERMOPT		";Energy formulation,",
#  if defined(ISOTHERMAL)
#    define GASOPT		";Isothermal gas (constant uint),",
#  elif defined(RADIATING)
#    if defined(DIFFUSING)
#      define GASOPT		";Gas radiates (diffusion model),",
#    elif defined(OPAQUE)
#      define GASOPT		";Gas radiates (opaque model),",
#    else
#      define GASOPT		";Gas radiates (optically thin),",
#    endif
#  elif defined(CONDUCTING)
#    define GASOPT		";Gas conducts heat,",
#  else
#    define GASOPT		";Ideal gas (shock heating),",
#  endif
#endif

#if defined(STARFORM) && defined(COMPVISC)
#  define STAROPT		";Star formation (using udotvis),",
#elif defined(STARFORM)
#  define STAROPT		";Star formation (using udot),",
#else
#  define STAROPT
#endif

#if defined(GRAVITY)
#  define DYNOPT		";Self-consistent gravity.",
#elif defined(EXTGRAV)
#  define DYNOPT		";External gravitational field.",
#elif defined(NOACCEL)
#  define DYNOPT		";No accelerations computed.",
#endif

#endif  // ! _switches_h
