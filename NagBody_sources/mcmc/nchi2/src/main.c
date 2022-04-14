/*==============================================================================
    NAME: main.c				[nchi2]
    Written by: M.A. Rodriguez-Meza
    Starting date: January 2018
    Purpose: Main routine
    Language: C
    Comments and notes:
    Info: M.A. Rodriguez-Meza
    Depto. de Fisica, ININ
    Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
    e-mail: marioalberto.rodriguez@inin.gob.mx
    https://github.com/rodriguezmeza

    Major revision:
    Copyright: (c) 2005-2020 Mar.  All Rights Reserved
================================================================================
 
    Use: nchi2 -help
    Input: 	Command line parameters, Parameters file
    Output: ...
    Units:
    History:
    Comments and notes:
    References:
================================================================================
    Legal matters:
    The author does not warrant that the program and routines it contains
    listed below are free from error or suitable for particular applications,
    and he disclaims all liability from any consequences arising from their use.
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

