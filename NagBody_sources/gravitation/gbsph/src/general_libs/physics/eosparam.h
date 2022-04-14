/*==============================================================================
	HEADER: eosparam.h				[General_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date:	November, 2001
	Purpose: Headers of utilities for equations of state
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:
	Comments and notes: Parameters that define the equation of state
						and the thermodynamical processes
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2001-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/


/*
Definiciones para las ecuaciones de estado
y tipo de proceso a usar.
*/

#define ADIABATICIDEALEQUATIONOFSTATE			// Original in gadget
//#undef  ADIABATICIDEALEQUATIONOFSTATE

#if !defined(ADIABATICIDEALEQUATIONOFSTATE)

#define ISOTHERMALIDEALEQUATIONOFSTATE
//#undef  ISOTHERMALIDEALEQUATIONOFSTATE

#define BAROTROPICIDEALEQUATIONOFSTATE
#undef  BAROTROPICIDEALEQUATIONOFSTATE

#endif	// ! ADIABATICIDEALEQUATIONOFSTATE


#ifdef ADIABATICIDEALEQUATIONOFSTATE			// Original in gadget
#define  GAMMA         (5.0/3.0)
#define  GAMMA_MINUS1  (GAMMA-1)
#endif

#ifdef BAROTROPICIDEALEQUATIONOFSTATE
#define  GAMMA         (5.0/3.0)
#define  GAMMA_MINUS1  (GAMMA-1)
#endif

#ifdef ISOTHERMALIDEALEQUATIONOFSTATE

//	This three parameters are for calculations
//		with the perturbed isothermal gaussian gas sphere

//#define SOUNDSPEED_ISOTERM	0.369331351
//#define	SOUNDSPEEDSQR		0.136406
//#define RHO_0			1.04865				

												// Para una esfera uniforme de
												// gas (Units: [GMR]) que
												// corresponde a M=1M_sun y
												// R = 5.0 x 10^16 cm
 
//	This three parameters are for calculations
//		with the perturbed isothermal uniform gas sphere

// L.D.G. Sigalotti and J. Klapp,
//		AA 319, 547 (1997).

#define SOUNDSPEED_ISOTERM	0.31547			
#define	SOUNDSPEEDSQR		0.0995216		
#define RHO_0			0.23873241			

											// For a uniform gas sphere	(Units: [GMR])
											// which corresponds to a M=1M_sun and
											// R = 3.20 x 10^16 cm

 
//	This three parameters are for calculations
//		with the perturbed isothermal uniform gas sphere
//	Truelove et al.	ApJ 495, 821 (1998).


//#define SOUNDSPEED_ISOTERM	0.394338		// La condicion inicial se genera con
//#define	SOUNDSPEEDSQR		0.155502		// models model_type=pertisothermalsphere
//#define RHO_0			0.23873241
											// For a uniform gas sphere	(Units: [GMR])
											// which corresponds to a M=1M_sun and
											// R = 5 x 10^16 cm

// Los siguientes parametros fueron calculados en el notebook "Colapso.nb"
// y describen un gas en un disco galactico. 

//#define SOUNDSPEED_ISOTERM	0.0220162
//#define	SOUNDSPEEDSQR		0.000484712		
//#define RHO_0			0.23873241


#endif


#define RHOCRITICAL		3.14228e3
#define	KBAROTROPIC		6.3638e-4

