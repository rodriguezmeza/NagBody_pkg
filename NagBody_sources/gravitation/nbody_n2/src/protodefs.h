/*==============================================================================
	HEADER: protodefs.h				[nbody_n2]
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
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.																		!
==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void direct_gravcalc(bodyptr, int);

void inputdata(void);
void output(void);
void restorestate(string);

void MainLoop(void);
void StartRun(string, string, string, string);
void StartOutput(void);
void checkstop(void);
void EndRun(void);

#endif // ! _protodefs_h
