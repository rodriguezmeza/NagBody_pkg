/*==============================================================================
   NAME: main.c						[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 1, 2006
	Purpose: Analysis of data coming from molecular dynamics codes (md_lj_tree,
		md_blj_tree, md_lj_n2)
	Language: C
	Info:	Mario A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revision: November 2008;
	Copyright: (c) 2007-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================

	Use: analysis_md -help
	Input:
	Output:
	Units:
	History:
     		Version 0.2 Final


	Comments and notes:
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#define global
#include "globaldefs.h"
#include "cmdline_defs.h"

int main(int argc, string argv[])						// CHECK 2D --- OK!!!
{
	gd.cpuinit = cputime();

    InitParam(argv, defv);
	startrun(argv[0], HEAD1, HEAD2, HEAD3);
    startoutput();
    data_analysis();
	code_endrun(gd.outlog, gd.cpuinit);
}

