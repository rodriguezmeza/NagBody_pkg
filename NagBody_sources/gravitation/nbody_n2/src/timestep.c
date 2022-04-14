/*==============================================================================
	MODULE: timestep.c				[nbody_n2]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: compute a timestep evolution of the system
	Language: C
	Use: 'stepsystem();'
	Routines and functions: stepsystem
	Modules, routines and external headers: stdinc, vectmath
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 2007; November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.																		!
==============================================================================*/

#include "globaldefs.h"

local void forcecalc(void);  
local void stepsystem(void);

void MainLoop(void)
{

    if (gd.nstep == 0) {
        forcecalc();
        output();
    }
    if (gd.dtime != 0.0)
        while (cmd.tstop - gd.tnow > 0.01*gd.dtime) {
            stepsystem();
            output();
			checkstop();
			if (gd.stopflag) break;
        }
}

local void forcecalc(void)
{
    direct_gravcalc(bodytab, cmd.nbody);
}

local void stepsystem(void)
{
    bodyptr p;

	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);
        ADDMULVS(Pos(p), Vel(p), gd.dtime);
    }
    forcecalc();
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);
    }
    gd.nstep++;
    gd.tnow = gd.tnow + gd.dtime;
}

