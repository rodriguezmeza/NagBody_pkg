/*==============================================================================
	HEADER: Physical Constants				[General_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date:	November, 2001
	Purpose: Headers of utilities for input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes: Reference: Physics Vadem Mecum
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/


/*
	Constantes físicas en el sistema cgs.
*/

#ifndef _physconstants_h
#define _physconstants_h

#ifndef G_CONST
#define G_CONST		0.00000006672
// #define G_CONST		0.00000000006672
#endif

#ifndef	M_SUN
#define M_SUN	1.989E33
#endif

#define K_BOLTZMANN		1.380658E-16
#define BOLTZMANNCONSTANT	1.381E-23		// Joules/Kelvin
#define ANGSTROM			1.0E-10			// Metros
#define KILO			1.0E3
#define NANOSEC			1.0E-12

#define M_PROTON		1.6726231E-24

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10				// Problemas si se definen 
											// variables con el mismo
											// nombre (e.g., params.c)
#define  PLANCK      6.6262e-27
#define  PROTONMASS  1.6726e-24

											// MARCA DE REVISION COSMOS
#define  HUBBLE      3.2407789e-18			// in h/sec

// EQUATION OF STATE PARAMETERS .....
#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

#define	SOUNDSPEEDSQR		0.155502		// Added to this header file ... see eosparam.h ...


// CAMBIO DE UNIDADES ----------------------------------------------------------
#define  CM_PER_MPC  3.085678e24
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#endif // !_physconstants_h

