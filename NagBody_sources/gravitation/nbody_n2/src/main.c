/*==============================================================================
   NAME: main.c						[nbody_n2]
	Written by: Mario Alberto Rodriguez-Meza
	Starting date: January, 2005
	Purpose: Simulation (N-body) of a selfgravitating N-body system.
		N*N force computation is done.
	Language: C
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revision:	July 2007; November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================

	Use: nbody_n2 -help
	Input: 	Command line parameters, Parameters file and/or icfile
	Output: snap.dat ...
	Units:
	History:
	Comments and notes: ...
	References:	Barnes Treecode, NEMO project, Gadget, Rapaport's book
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
	MainLoop();
	EndRun();
    return 0;
}

