/* ==============================================================================
!	MODULE: time_step.c															!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date:	February 2005												!
!	Purpose: compute a timestep evolution of the system							!
!	Language: C																	!
!	Use: tree_ljforce(); stepsystem()											!
!	Routines and functions:	tree_ljforce, stepsystem							!
!	External modules, routines and headers:	stdinc, vectmath, mathfns,			!
!		global_defs, proto_defs													!
!	Comments and notes:															!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.							!
!		e-mail: marioalberto.rodriguez@inin.gob.mx
!		http://www.astro.inin.mx/mar											!
!																				!
!	Major revisions:															!
!	Copyright: (c) 2005-2011 Mar.  All Rights Reserved.							!
!===============================================================================
!	Legal matters:																!
!	The author does not warrant that the program and routines it contains		!
!	listed below are free from error or suitable for particular applications,	!
!	and he disclaims all liability from any consequences arising from their		!
!	use.																		!
!==============================================================================*/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "global_defs.h"
#include "proto_defs.h"

local void forcecalc(bodyptr, int);  
//local void forcecalc(void);			// Opcion para solo correr barnes

void tree_ljforce(void)
{
    bodyptr p;

	DO_BODY(p,bodytab,bodytab+nbody)
        Update(p) = TRUE;
	if (forcemethod_int!=4 && forcemethod_int!=5 && forcemethod_int!=6)
		maketree(bodytab, nbody);		// Los metodos directos y cells no requiere del arbol
	forcecalc(bodytab, nbody);
//	forcecalc();						// Opcion para solo correr barnes
}

void stepsystem(void)
{
    bodyptr p;
    real velsq, Ekin, fEkin;
	
	Ekin=0.0;

	DO_BODY(p,bodytab,bodytab+nbody) {
        ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);     
        ADDMULVS(Pos(p), Vel(p), dtime);
    }

	DO_BODY(p,bodytab,bodytab+nbody) {
		VWrapAll (Pos(p));					// Periodic boundary condition
	}

    tree_ljforce(); 

	DO_BODY(p,bodytab,bodytab+nbody) {
        ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);     
        DOTVP(velsq, Vel(p), Vel(p));           
        Ekin += 0.5 * Mass(p) * velsq;
    }

#if (NDIM==3)
		fEkin=2.0*Ekin/((real)(3*nbody-3));
#else
#if (NDIM==2)
		fEkin=2.0*Ekin/((real)(2*nbody-2));
#endif
#endif
	Scale = rsqrt(temperature/fEkin);

	DO_BODY(p,bodytab,bodytab+nbody) {
        MULVS(Vel(p), Vel(p), Scale);		// Scaling can be done in EvalProps
	}

    nstep++;
    tnow = tnow + dtime;
}

#define BARNES		0
#define NULLMETHOD	1
#define NORMAL		2
#define NBLIST		3
#define DIRECT		4
#define CELLSMETHOD	5
#define DIRECT2		6
#define NORMAL2		7
#define BARNES2		8


//local void forcecalc(void)						// Opcion para solo correr barnes
local void forcecalc(bodyptr btab, int nbody)
{
    switch(forcemethod_int) {
        case BARNES:
            ljforcecalc_barnes(btab, nbody); break;
//            ljforcecalc_barnes(); break;			// Opcion para solo correr barnes
        case BARNES2:
            ljforcecalc_barnes2(btab, nbody); break;
        case NULLMETHOD:
            printf("\n\trunning default method (Normal)...\n");
            ljforcecalc_normal(btab, nbody); break;
        case NORMAL:
            ljforcecalc_normal(btab, nbody); break;
        case NORMAL2:
            ljforcecalc_normal2(btab, nbody); break;
        case NBLIST:
            ljforcecalc_nblist(btab, nbody); break;
        case DIRECT:
            ljforcecalc_direct(btab, nbody); break;
        case DIRECT2:
            ljforcecalc_direct2(btab, nbody); break;
        case CELLSMETHOD:
            ljforcecalc_cellsmethod(btab, nbody); break;
        default:
            printf("\n\tforcecalc_method: Unknown method...");
            printf("\n\trunning default force calculation method (Normal)...\n"); 
            ljforcecalc_normal(btab, nbody); break;
    }

}

#undef BARNES
#undef NULLMETHOD
#undef NORMAL
#undef NBLIST
#undef DIRECT 
#undef CELLSMETHOD
#undef DIRECT2
#undef NORMAL2
#undef BARNES2

