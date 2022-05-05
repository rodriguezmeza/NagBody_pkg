/*==============================================================================
	NAME: main.c		[galaxy_models]
	Written by: Mario Alberto Rodriguez-Meza.
	Starting date: May 2006
	Purpose: Main routine
	Language: C
	Comments and notes:
	Info: M.A. Rodriguez-Meza,
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revision:	May 2006; June 6, 2008;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved.
================================================================================

	Use: ./galaxy_models -help (or man galaxy_models)
	Input: 	Command line parameters, Parameters file and/or icfile
	Output: ...
	Units:
	History:
	Comments and notes:
	References:
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their
	use.
==============================================================================*/


#define global

#include "globaldefs.h"
#include "cmdlinedefs.h"
#include "protodefs.h"

int main(int argc, string argv[])
{
	double cpustart;

    gd.cpuinit = cputime();

    InitParam(argv, defv);
    StartRun(argv[0], HEAD1, HEAD2, HEAD3);

	cpustart = cputime();

	InitDisk();
	InitBulge();
	InitHalo();

	fprintf(stdout,"\n\nInitDisk, InitBulge and InitHalo CPU time: %g\n",
			cputime()-cpustart);
	fflush(stdout);
	cpustart = cputime();

	DiskVel();
	DiskStat();

	fprintf(stdout,"\n\nDiskVel and DiskStat CPU time: %g\n",cputime()-cpustart);
	fflush(stdout);
	cpustart = cputime();

	if (cmd.usehalo && cmd.selfghal) HaloVel();
	if (cmd.usebulge && cmd.selfgbul) BulgeVel();

	fprintf(stdout,"\n\nHaloVel and BulgeVel CPU time: %g\n",
			cputime()-cpustart);
	fflush(stdout);
	cpustart = cputime();

	cmtv();
	AddSat();
	StackMod();

	fprintf(stdout,"\n\nCMTV, AddSat, StackMod CPU time: %g\n",
			cputime()-cpustart);
	fflush(stdout);

	output();
	EndRun();

    return 0;
}

