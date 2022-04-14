/*==============================================================================
    HEADER: constants_defs.h		[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date: January 2018
    Purpose: Definitions of global variables and parameters
    Language: C
    Use: '#include "global_defs.h"
    Use in routines and functions:
    External headers: stdinc.h, data_struc_defs.h
    Comments and notes:
    Info: M.A. Rodriguez-Meza
        Depto. de Fisica, ININ
        Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
        e-mail: marioalberto.rodriguez@inin.gob.mx
        https://github.com/rodriguezmeza

    Major revisions:
    Copyright: (c) 2005-2020 Mar.  All Rights Reserved
================================================================================
    Legal matters:
    The author does not warrant that the program and routines it contains
    listed below are free from error or suitable for particular applications,
    and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

 
#ifndef _constants_defs_h
#define _constants_defs_h

#define FOURPI      12.56637061435917
#define KILO        1.0E3
#define MEGA        1.0E6
#define PARSEC      3.0857E16                       // In meters
#define MPC         MEGA*PARSEC                     // In meters
#define R0          KILO*PARSEC                     // In meters
#define M0          4.62489E35                      // Unit of mass in kilograms
#define RHO0        M0/(FOURPI*rpow(R0,3.0))        // In kilograms meters^(-3)
#define GN          6.67259E-11                     // Gravitation constant
                                                    // in meters^3 kilogram^(-1)
                                                    // second^(-2)
#define MSUN        1.989E30                        // Sun's mass in kilograms
#define H0          100*(0.69)*KILO/(MEGA*PARSEC)   // In seconds^(-1) (h->0.69)
#define RHOCRITIC   3.0*rsqr(H0)/(8.0*PI*GN)        // In kilograms meters^(-3)
#define GDAGGER     1.2E-10                         // SPARC acceleration constant
#define RHOFACTOR   (M0/MSUN)/(FOURPI*1.0E9)        // Conversion factor for rho's
#define MFACTOR     (M0/MSUN)                       // Conversion factor for mass
#define DIGITS      4                               // Digit to save tables
#define DDIGITS     3                               // Digit to save error tables



#endif // ! _constants_defs_h

