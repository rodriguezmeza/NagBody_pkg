/*==============================================================================
   NAME: main.c						[md_ic_model]
	Written by: Mario Alberto Rodriguez-Meza
	Starting date: January, 2012
	Purpose: Simulation (MD) of Lennard-Jones gas dynamics.
		Hierarchical tree force computation is done.
	Language: C
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revision:	March 2012
	Copyright: (c) 2005-2012 Mar.  All Rights Reserved
================================================================================

	Use: md_ic_model -help
	Input: Command line parameters, Parameters file and/or icfile
	Output: thermo.dat rdf.dat vel.dat snap.dat ...
	Units: eps=sigma=mass(p)=kB=1
	History:
	Comments and notes: ...
	References:	Barnes Treecode, NEMO project, Gadget, Rapaport's book
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "./general_libs/general/stdinc.h"
#include "./general_libs/math/vectmath.h"
#include "./general_libs/general/getparam.h"
#define global
#include "globaldefs.h"
#include "cmdline_defs.h"
#include "protodefs.h"

int main(int argc, string argv[])						// CHECK 2D --- OK!!!
{
    gd.cpuinit = cputime();
    InitParam(argv, defv);
    StartRun(argv[0], HEAD1, HEAD2, HEAD3);
	MainLoop();
	EndRun();
    return 0;
}

