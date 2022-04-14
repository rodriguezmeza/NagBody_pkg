/* =============================================================================
	MODULE: timestep.c					[gbsph]
	Written by: M.A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Evolve the N-body system
	Language: C
	Use: 'first_treeforce();', 'stepsystem();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/vectmath.h"
//#include "../../../General_libs/general/getparam.h"
//#include "../../../General_libs/math/numrec.h"
#include "globaldefs.h"
#include "protodefs.h"

local void zero_treeforce(void);
local void first_treeforce(void);
local void stepsystem(void);

void MainLoop(void)
{
    if (gd.nstep == 0) {                           
        first_treeforce();
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

local void zero_treeforce(void)
{
    bodyptr p;

    for (p = bodytab; p < bodytab+cmd.nbody; p++) 
        Update(p) = TRUE;                     
    maketree(bodytab, cmd.nbody);                  
    for (p = bodytab; p < bodytab+cmd.nbody; p++) 
        CLRV(Acc(p));                     
    forcereport(); 
    gd.nstep_grav++; 
}

local void first_treeforce(void)
{
    bodyptr p;

    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        Update(p) = TRUE;  
    maketree(bodytab, cmd.nbody);
    
    forcecalc(bodytab, cmd.nbody);

    if (scanopt(cmd.force_models, "external-potential"))
        external_forcecalc(bodytab, cmd.nbody); 	
    
    forcereport(); 
    gd.nstep_grav++;                                
}

local void stepsystem(void)
{
    bodyptr p;

    if (scanopt(cmd.forcecalc_method, "motion_without_leader")) {
        for (p = bodytab; p < bodytab+cmd.nbody; p++) {
            Vel(p)[0] = Acc(p)[0]*0.5*gd.dtime;
            Vel(p)[1] = Acc(p)[1]*0.5*gd.dtime;
            Vel(p)[2] = Acc(p)[2]*0.5*gd.dtime;
            ADDMULVS(Pos(p), Vel(p), gd.dtime);
        }
        first_treeforce();
        for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
            Vel(p)[0] = Acc(p)[0]*0.5*gd.dtime;
            Vel(p)[1] = Acc(p)[1]*0.5*gd.dtime;
            Vel(p)[2] = Acc(p)[2]*0.5*gd.dtime;
        }
        gd.nstep++;
        gd.tnow = gd.tnow + gd.dtime;
    } else {
        for (p = bodytab; p < bodytab+cmd.nbody; p++) {
            ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);
            ADDMULVS(Pos(p), Vel(p), gd.dtime);
        }
        first_treeforce();
        for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
            ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);
        }
        gd.nstep++;
        gd.tnow = gd.tnow + gd.dtime;
    }
}

