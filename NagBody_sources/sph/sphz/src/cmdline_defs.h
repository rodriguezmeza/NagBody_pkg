//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)

#ifndef _cmdline_defs_h
#define _cmdline_defs_h


#include "switches.h"

string defv[] = {		";SPH/N-body simulation code.",
				THERMOPT
				GASOPT
				STAROPT
				DYNOPT
    "in=",			";Input file for initial conditions",
    "out=",			";Output file for N-body frames",
    "gamma=1.666667",		";Ratio of specific heats",
#if defined(RADIATING)
    "uradmax=1.0",		";Uinternal at max. of cooling curve",
    "lambmax=1.0",		";Max. cooling rate at unit density",
#endif
#if defined(DIFFUSING)
    "sigmastar=1.0",		";Stefan-Boltzmann parameter.",
				";Units are flux / specific_energy^4.",
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    "opacity=1.0",		";Radiation opacity parameter",
#endif
#if defined(CONDUCTING)
    "conduct=1.0",		";Thermal conduction parameter",
#endif
#if defined(STARFORM)
    "starprob=0.01",		";Star form. prob. per unit time",
    "rhoindx=0.0",		";Power-law index for density",
    "udotindx=0.0",		";Power-law index for dissipation",
    "starseed=12345",		";Random seed for star formation",
#endif
    "alpha=1.0",		";Bulk viscosity parameter",
    "beta=2.0",			";vN-R viscosity parameter",
    "nsmooth=40",		";Bodies in smoothing volume",
    "nbucket=16",		";Bodies in leaves of KD tree",
    "slope=0.0",		";Kernel slope at origin",
    "courant=0.25",		";Courant condition parameter",
    "dtime=1/64",		";Basic integration timestep",
    "fdrag=0.0",		";Velocity damping factor (1/t)",
#if defined(GRAVITY)
    "eps=0.025",		";Gravitational smoothing length",
    "usequad=false",		";If true, use quad moments",
#if !defined(QUICKSCAN)
    "theta=1.0",		";Force accuracy parameter",
#endif
#elif defined(EXTGRAV)
    "gravgsp=",			";Input GSP for external gravity",
#endif
    "options=",			";Various control options.",
				";Choices include: corrfunc, levelhist,",
				";new-tout, nolimit, fixstep, lockstep,",
				";reset-time, micro-star, bh86, and sw94.",
    "outputs=",			";Optional particle data written to output.",
				";Key dynamical variables always included.",
    "tstop=2.0",		";Time to stop integration",
    "dtout=1/16",		";Data output timestep",
    "save=",			";Write state file as code runs",
    "restore=",			";Continue run from state file",
    "VERSION=0.2",		";Mar 2005-2010",
    NULL,
};

#endif  // ! _cmdline_defs_h
