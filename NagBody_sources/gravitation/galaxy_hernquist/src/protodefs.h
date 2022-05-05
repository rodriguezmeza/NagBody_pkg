/*==============================================================================
	HEADER: protodefs.h				[galaxy_hernquist]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: main, direct_gravcalc,
					nbody_n2_io, startrun, timestep
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:  November 2008;
	Copyright: (c) 2005-2014 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void MainLoop(void);
void StartRun(string, string, string, string);
void EndRun(void);

extern void create_hernquist(double *mtot_type, int *ntot_type, double *ScaleTable);
extern void create_disk(double *mtot_type, int *ntot_type, double *ScaleTable);

extern double bessi0(double), bessk0(double), bessi1(double), bessk1(double);
extern double bessi0(double), bessk0(double), bessi1(double), bessk1(double);
extern float qromb(float (*func)(float), float a, float b);

#endif // ! _protodefs_h
