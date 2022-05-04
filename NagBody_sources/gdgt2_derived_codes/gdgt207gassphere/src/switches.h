
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS


#ifndef switches_h
#define switches_h

#define PERIODIC						// lcdm_gas
#undef PERIODIC

#define UNEQUALSOFTENINGS				// galaxy, cluster
#undef UNEQUALSOFTENINGS

#define PEANOHILBERT					// galaxy, gassphere, lcdm_gas, cluster
//#undef PEANOHILBERT

#define WALLCLOCK						// galaxy, gassphere, lcdm_gas, cluster
//#undef WALLCLOCK

#define PMGRID	256						// lcdm_gas, cluster
#undef PMGRID

#define PLACEHIGHRESREGION	3
#undef PLACEHIGHRESREGION

#define ENLARGEREGION	1.2
#undef ENLARGEREGION

#define ASMTH	1.25
#undef ASMTH

#define RCUT	4.5
#undef RCUT

#define DOUBLEPRECISION
#undef DOUBLEPRECISION

#define DOUBLEPRECISION_FFTW
#undef DOUBLEPRECISION_FFTW

#define SYNCHRONIZATION					// galaxy, gassphere, lcdm_gas, cluster
//#undef SYNCHRONIZATION

#define FLEXSTEPS
#undef FLEXSTEPS

#define PSEUDOSYMMETRIC
#undef PSEUDOSYMMETRIC

#define NOSTOP_WHEN_BELOW_MINTIMESTEP
#undef NOSTOP_WHEN_BELOW_MINTIMESTEP

#define NOPMSTEPADJUSTMENT
#undef NOPMSTEPADJUSTMENT

#define HAVE_HDF5
#undef HAVE_HDF5

#define OUTPUTPOTENTIAL
#undef OUTPUTPOTENTIAL

#define OUTPUTACCELERATION
#undef OUTPUTACCELERATION

#define OUTPUTCHANGEOFENTROPY
#undef OUTPUTCHANGEOFENTROPY

#define OUTPUTTIMESTEP
#undef OUTPUTTIMESTEP

#define NOGRAVITY
#undef NOGRAVITY

#define NOTREERND
#undef NOTREERND

#define NOTYPEPREFIX_FFTW
#undef NOTYPEPREFIX_FFTW

#define LONG_X	60
#undef LONG_X

#define LONG_Y	5
#undef LONG_Y

#define LONG_Z	0.2
#undef LONG_Z

#define TWODIMS
#undef TWODIMS

#define SPH_BND_PARTICLES
#undef SPH_BND_PARTICLES

#define NOVISCOSITYLIMITER
#undef NOVISCOSITYLIMITER

#define COMPUTE_POTENTIAL_ENERGY
#undef COMPUTE_POTENTIAL_ENERGY

#define LONGIDS
#undef LONGIDS

#define ISOTHERMAL
#undef ISOTHERMAL

#define SELECTIVE_NO_GRAVITY	2+4+8+16
#undef SELECTIVE_NO_GRAVITY

#define FORCETEST	0.1
#undef FORCETEST

#define MAKEGLASS	262144
#undef MAKEGLASS


#endif		// ! switches_h

