/*==============================================================================
	MODULE: direct_gravcalc.c			[nbody_n2]
	Written by: M.A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Direct force computation (Complexity N*N)
	Language: C
	Use: 'direct_gravcalc(bodyptr btab, int nbody);'
	Routines and functions: stepsystem
	Modules, routines and external headers: stdinc, vectmath, ...
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
        http://www.inin.gob.mx/

	Major revisions: July 2007; November 2008;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.																		!
==============================================================================*/
 
#include "globaldefs.h"

void direct_gravcalc(bodyptr btab, int nbody)
{
    bodyptr p, q;
    double cpustart;
    vector pos0, acc0;
    real phi0;
    real eps2, dr2, drab, phi_q, mr3i;
    vector dr;
 
    cpustart = cputime(); 

    eps2 = cmd.eps * cmd.eps;               

	DO_BODY(p, btab, btab+nbody) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    gd.nbbcalc = 0;

	DO_BODY(p, btab, btab+nbody-1) {
		DO_BODY(q, p+1, btab+nbody) {
			DOTPSUBV(dr2, dr, Pos(q), Pos(p));
			dr2 += eps2;
			drab = rsqrt(dr2);
			phi_q = Mass(q)/drab;
			Phi(p) -= phi_q;
			Phi(q) -= phi_q;
			mr3i = phi_q / dr2;
			ADDMULVS(Acc(p), dr, mr3i);
			ADDMULVS(Acc(q), dr, -mr3i);
			gd.nbbcalc += 2;						// nbbcalc = n(n-1)
		}
	}
    
    

	gd.cpuforce = cputime() - cpustart;
}

