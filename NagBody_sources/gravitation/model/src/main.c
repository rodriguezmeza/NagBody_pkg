/*==============================================================================
   NAME: main.c			[model]
	Written by: Mario Alberto Rodriguez-Meza
	Starting date: February 1, 2006
	Purpose: Generation of initial condition for Simulation of general body
			systems and smoothed particle hydrodynamics
	Language: C
	Info:	M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revision: July 23, 2007
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
==============================================================================

	Use: models -help
	Input:
	Output:
	Units:
	History:
	Comments and notes:
	References:
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/


#define global
#include "globaldefs.h"
#include "cmdline_defs.h"

int main(int argc, string argv[])
{
	gd.cpuinit = cputime();
    InitParam(argv, defv);
    StartRun(argv[0], HEAD1, HEAD2, HEAD3);
	output();
	EndRun();
}

