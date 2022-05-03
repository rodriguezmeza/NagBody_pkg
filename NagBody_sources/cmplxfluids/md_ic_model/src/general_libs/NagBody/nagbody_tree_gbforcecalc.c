/* ==============================================================================
	MODULE: nagbody_tree_gbforcecalc.c		[NagBody]
	Written by: M.A. Rodriguez-Meza
	Starting date:
	Purpose: Lennard-Jones force computation (Several methods)
	Language: C
	Use: forcecalc();
	Routines and functions:	tree_ljforce, stepsystem
	External modules, routines and headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ,
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

/*
NOTA: El orden de los terminos en el calculo de la aceleracion de cada particula
	es importante. Diferentes ordenamientos (diferentes metodos de calculo)
	dan diferentes resultados despues de un numero grande de iteraciones.
	Por ejemplo, un calculo con 
		md_lj_tree forcecalc_method=direct nbody=512 tstop=3.6
	da un resultado diferente a
		md_lj_tree forcecalc_method=normal nbody=512 tstop=3.6
	debido a que el orden de aparicion de los vecinos cercanos a cada particula
	es diferente. Despues de unas 3600 iteraciones el momento angular comienza
	a mostrarse diferente en los resultados del momento angular de los dos casos.
*/

#include "../general/stdinc.h"
#include "../math/mathfns.h"
#include "../math/vectmath.h"
#include "../math/mathutil.h"
#include "nagbody.h"
#include "../general/constant.h"

// Bloque para normal
local void normal_walktree(bodyptr, nodeptr, real, vector, 
						   real *, vector, global_data_tree_tljforcecalc *);
local void sumnode(bodyptr, cellptr, cellptr, vector, real *, vector, 
					global_data_tree_tljforcecalc *);

// Bloque para normal2
local void normal_walktree2(bodyptr, nodeptr, real, 
							global_data_tree_tljforcecalc *);
local void sumnode2(bodyptr, cellptr, cellptr, global_data_tree_tljforcecalc *);

// Bloque para normal3
local void normal_walktree3(bodyptr, nodeptr, real, 
							global_data_tree_tljforcecalc *);
local void sumnode3(bodyptr, cellptr, cellptr, global_data_tree_tljforcecalc *);

// Bloque para barnes
local void walktree_barnes(nodeptr *, nodeptr *, cellptr, cellptr,
                    nodeptr, real, vector, global_data_tree_tljforcecalc *);
local bool accept_blj(nodeptr, nodeptr,	global_data_tree_tljforcecalc *);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector, global_data_tree_tljforcecalc *);
local void gravsum(bodyptr, cellptr, cellptr, global_data_tree_tljforcecalc *);
local void sumnode_barnes(bodyptr, cellptr, cellptr, vector, real *, vector, 
							global_data_tree_tljforcecalc *);

// Bloque para barnes2
local void walktree_barnes2(nodeptr *, nodeptr *, cellptr *, bodyptr *,
                    nodeptr, real, vector, global_data_tree_tljforcecalc *);
local void walksub2(nodeptr *, nodeptr *, cellptr *, bodyptr *,
                   nodeptr, real, vector, global_data_tree_tljforcecalc *);
local void sumnode_barnes2(bodyptr, bodyptr *, bodyptr *, 
							global_data_tree_tljforcecalc *);

// Bloque para nblist
local void normal_walktree_nblist(bodyptr, nodeptr, real, vector, 
						   real *, vector, global_data_tree_tljforcecalc *);
local void sumnode_nblist_01(bodyptr, cellptr, cellptr, vector, real *, vector, 
								global_data_tree_tljforcecalc *);
local void sumnode_nblist_02(bodyptr, vector, real *, vector, 
								global_data_tree_tljforcecalc *);

local void bljFactors(bodyptr, bodyptr, real *, real *, real *,
		real *, real *, real *, global_data_tree_tljforcecalc *);

// Bloque que define la fuerza y el potencial de interaccion ===================

local void ForcePotentialCalc(bodyptr, bodyptr, int, real *, real *, real, 
								real, real, real, real, real);
local void PotentialCalc(bodyptr, bodyptr, int, real *, real, real, real, real, real);
local real PotentialFunction(int, real, real , real, real);
local real DPotentialFunction(int, real, real , real, real);

// Bloque Lennard-Jones ....
local void ForcePotentialCalc_LJ(real *, real *, real, 
								real, real, real, real, real);
local void PotentialCalc_LJ(real *, real, real, real, real, real);
local real PotentialFunction_LJ(real, real , real, real);
local real DPotentialFunction_LJ(real, real , real, real);
local void PotentialParameters_LJ(real, real, real, real, real, real,
	real, real, real, real, real, real, 
	real, real, real, real, real, real, global_data_tree_tljforcecalc *);

// Bloque Logaritmo natural ....
local void ForcePotentialCalc_LN(real *, real *, real, 
								real, real, real, real, real);
local void PotentialCalc_LN(real *, real, real, real, real, real);
local real PotentialFunction_LN(real, real , real, real);
local real DPotentialFunction_LN(real, real , real, real);
local void PotentialParameters_LN(real, real, real, real, real, real,
	real, real, real, real, real, real, 
	real, real, real, real, real, real, global_data_tree_tljforcecalc *);

// Bloque de fuerza de Lorentz ....
local void ForcePotentialCalc_Lorentz(bodyptr, bodyptr, real *, real *, real, 
								real, real, real, real, real);
local void PotentialCalc_Lorentz(bodyptr, bodyptr, real *, real, real, real, real, real);
local real PotentialFunction_Lorentz(real, real , real, real);
local real DPotentialFunction_Lorentz(real, real , real, real);
local void PotentialParameters_Lorentz(real, real, real, real, real, real,
	real, real, real, real, real, real, 
	real, real, real, real, real, real, global_data_tree_tljforcecalc *);

// Fin del bloque de fuerzas y potenciales =====================================

local void TransportCalc(bodyptr, vector, real, real);
local void TransportCalc2(bodyptr, bodyptr, vector, real, real);

#define DIAGNOSTICS					// Util para diagnosticar o expulgar... 
#undef DIAGNOSTICS					// Esta definicion es local ...
									// Remover cuando ya no sea util...

// COMIENZA METODO NORMAL DE CALCULO DE LA FUERZA

void ind_tljforcecalc_normal(bodyptr btab, int nbody, 
	global_data_tree_tljforcecalc *gdforce, global_data_tree *gdtree, 
	bodyptr p, real *uSum, real *virSum)
{
    double cpustart;
	vector pos0, acc0;
	real phi0;

    cpustart = cputime();

    gdforce->nbbcalc = gdforce->nbccalc = 0;
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;

	SETV(pos0, Pos(p));
	phi0 = 0.0;
	CLRV(acc0);
	normal_walktree(p, ((nodeptr) gdtree->root), gdtree->rsize,
		pos0, &phi0, acc0, gdforce);
	DIVVS(Acc(p), acc0, Mass(p));
	Phi(p) += phi0; gdforce->uSum += Phi(p);

	*virSum = gdforce->virSum;
	*uSum = gdforce->uSum;

    gdforce->cpuindforce = cputime() - cpustart;
}

void tljforcecalc_normal(bodyptr btab, int nbody, 
	global_data_tree_tljforcecalc *gdforce, global_data_tree *gdtree)
{
    bodyptr p;
    double cpustart;
	vector pos0, acc0;
	real phi0;

    cpustart = cputime();                       

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    gdforce->nbbcalc = gdforce->nbccalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

    for (p = btab; p < btab+nbody; p++) {
		SETV(pos0, Pos(p));
		phi0 = 0.0;
		CLRV(acc0);
		normal_walktree(p, ((nodeptr) gdtree->root), gdtree->rsize, 
			pos0, &phi0, acc0, gdforce);
//		DIVVS(Acc(p), acc0, Mass(p));
		ADDMULVS(Acc(p), acc0, 1.0/Mass(p));
		Phi(p) += phi0; gdforce->uSum += Phi(p);
	}
	gdforce->virSum = 0.5*gdforce->virSum;
	gdforce->uSum = 0.5*gdforce->uSum;

    gdforce->cpuforce = cputime() - cpustart;
}

local void normal_walktree(bodyptr p, nodeptr q, real qsize, vector pos0, 
			real *phi0, vector acc0, global_data_tree_tljforcecalc *gdforce)
{
    nodeptr l;
	real drpq, drpq2, rcell;
    vector dr;
	int k;
	real Rcut;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
//				for (k = 0; k < NDIM; k++)
//					dr[k]=dr[k] - 
//						((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
				DOTVP(drpq2, dr, dr);
				drpq = rsqrt(drpq2);
				rcell= qsize * rsqrt((real)(NDIM))/2.0;

				if (Type(p) == BODY1 )
					Rcut = gdforce->Rcut11Max;
				else
					Rcut = gdforce->Rcut22Max;

                if ( drpq >= Rcut+rcell ) {
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree(p,l,qsize/2,pos0,phi0,acc0, gdforce);
					}
				}
			} else {
					sumnode(p, ((cellptr) q),( (cellptr) q+1),
						pos0,phi0,acc0, gdforce);
			}
		}
	}
}

local void sumnode(bodyptr p, cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0, 
				   global_data_tree_tljforcecalc *gdforce)
{
    cellptr q;
    real dr2, fPot, fAcc;
    vector dr;
	real RcutSq, ssq, fphi, fa, vc, dvc;

    for (q = start; q < finish; q++) {
		bljFactors(p, (bodyptr)q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			ForcePotentialCalc(p, (bodyptr)q, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, vc, dvc);
			*phi0 += fPot;
			ADDMULVS(acc0, dr, fAcc);
			gdforce->virSum += fAcc*dr2;
			gdforce->nbbcalc += 1;
printf("\n count=%d Pot, fAcc: %g %g %g %g %g \t %g",gdforce->nbbcalc, fPot, fAcc,acc0[0],acc0[1],acc0[2],*phi0);
			if (gdforce->computeTransport)
				TransportCalc(p, dr, fPot, fAcc);
		}
    }
}

// TERMINA METODO NORMAL DE CALCULO DE LA FUERZA


// COMIENZA METODO NORMAL2 DE CALCULO DE LA FUERZA

void tljforcecalc_normal2(bodyptr btab, int nbody, 
			global_data_tree_tljforcecalc *gdforce, global_data_tree *gdtree)
{
    bodyptr p;
    double cpustart;
//	vector pos0, acc0;
//	real phi0;
//	int id0;
 
    cpustart = cputime();                       

	DO_BODY(p, btab, btab+nbody) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    gdforce->nbbcalc = gdforce->nbccalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

	DO_BODY(p, btab, btab+nbody) {
		normal_walktree2(p, ((nodeptr) gdtree->root), gdtree->rsize, gdforce);
		DIVVS(Acc(p), Acc(p), Mass(p));
	}

    gdforce->cpuforce = cputime() - cpustart;
}

local void normal_walktree2(bodyptr p, nodeptr q, real qsize, 
							global_data_tree_tljforcecalc *gdforce)
{
    nodeptr l;
	real drpq, drpq2, rcell;
    vector dr;
	int k;
	real Rcut;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
//				VWrapAll_ptr(dr);
				DOTVP(drpq2, dr, dr);
				drpq = rsqrt(drpq2);
				rcell= qsize * rsqrt((real)(NDIM))/2.0;

				if (Type(p) == BODY1 )
					Rcut = gdforce->Rcut11Max;
				else
					Rcut = gdforce->Rcut22Max;

                if ( drpq < Rcut+rcell )
					DO_DESCENDENTS(l,q)
						normal_walktree2(p,l,qsize/2, gdforce);
			} else {
				if ( Id(p) < Id(q) ) {
					sumnode2(p,((cellptr) q),( (cellptr) q+1), gdforce);
					}
			}
		}
	}
}

local void sumnode2(bodyptr p, cellptr start, cellptr finish, 
					global_data_tree_tljforcecalc *gdforce)
{
    cellptr q;
    real dr2, fPot, fAcc;
    vector dr;
//	real rri, rri3;
//	int k;
	bodyptr qb;
	real RcutSq, ssq, fphi, fa, vc, dvc;

	DO_BODY(q, start, finish) {
		qb = ((bodyptr) q);
		bljFactors(p, qb, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
//		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			ForcePotentialCalc(p, qb, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, vc, dvc);
/*
			rri=ssq/dr2; rri3=rri*rri*rri;
			uVal = fphi*(rri3-1.0)*rri3;
			fVal = fa*rri3*(rri3-0.5)*rri;
*/
			ADDMULVS(Acc(p), dr, fAcc);
			ADDMULVS(Acc(qb), dr, -fAcc);
			Phi(p) += fPot;
			Phi(qb) += fPot;
			gdforce->virSum += fAcc*dr2;
			gdforce->uSum += fPot;
			gdforce->nbbcalc += 2;
			if (gdforce->computeTransport)
				TransportCalc2(p, qb, dr, fPot, fAcc);
		}
    }
}

// TERMINA METODO NORMAL2 DE CALCULO DE LA FUERZA

// COMIENZA METODO NORMAL3 DE CALCULO DE LA FUERZA

void tljforcecalc_normal3(bodyptr btab, int nbody, 
			global_data_tree_tljforcecalc *gdforce, global_data_tree *gdtree)
{
    bodyptr p;
    double cpustart;
//	vector pos0, acc0;
//	real phi0;
//	int id0;
 
    cpustart = cputime();                       

	DO_BODY(p, btab, btab+nbody) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
        Phi11(p) = 0.0;
        CLRV(Acc11(p));
        Phi12(p) = 0.0;
        CLRV(Acc12(p));
        Phi22(p) = 0.0;
        CLRV(Acc22(p));
    }

    gdforce->nbbcalc = gdforce->nbccalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
//	if (gdforce->computeTransport) {
//		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
//	}

	DO_BODY(p, btab, btab+nbody) {
		gdforce->virSum11=0.0; gdforce->virSum12=0.0; gdforce->virSum22=0.0;
		normal_walktree3(p, ((nodeptr) gdtree->root), gdtree->rsize, gdforce);
		ADDMULVS(Acc(p), Acc11(p), gdforce->fa11);
		ADDMULVS(Acc(p), Acc12(p), gdforce->fa12);
		ADDMULVS(Acc(p), Acc22(p), gdforce->fa22);
		Phi(p) += gdforce->fphi11*Phi11(p)+gdforce->fphi12*Phi12(p)
					+gdforce->fphi22*Phi22(p);
		gdforce->virSum += gdforce->fa11*gdforce->virSum11
							+gdforce->fa12*gdforce->virSum12
							+gdforce->fa22*gdforce->virSum22;
		gdforce->uSum += Phi(p);
		DIVVS(Acc(p), Acc(p), Mass(p));
	}
	gdforce->uSum *= 0.5;

    gdforce->cpuforce = cputime() - cpustart;
}

local void normal_walktree3(bodyptr p, nodeptr q, real qsize, 
							global_data_tree_tljforcecalc *gdforce)
{
    nodeptr l;
	real drpq, drpq2, rcell;
    vector dr;
	int k;
	real Rcut;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
//				VWrapAll_ptr(dr);
				DOTVP(drpq2, dr, dr);
				drpq = rsqrt(drpq2);
				rcell= qsize * rsqrt((real)(NDIM))/2.0;

				if (Type(p) == BODY1 )
					Rcut = gdforce->Rcut11Max;
				else
					Rcut = gdforce->Rcut22Max;

                if ( drpq < Rcut+rcell )
					DO_DESCENDENTS(l,q)
						normal_walktree3(p,l,qsize/2, gdforce);
			} else {
				if ( Id(p) < Id(q) ) {
					sumnode3(p,((cellptr) q),( (cellptr) q+1), gdforce);
					}
			}
		}
	}
}

local void sumnode3(bodyptr p, cellptr start, cellptr finish, 
					global_data_tree_tljforcecalc *gdforce)
{
    cellptr q;
    real dr2, fPot, fAcc;
    vector dr;
	real rri, rri3;
//	int k;
	bodyptr qb;
	real RcutSq, ssq, fphi, fa, vc, dvc;
//	matrix w;

	DO_BODY(q, start, finish) {
		qb = ((bodyptr) q);
		if ( Type(p) == Type(qb) ) {
			if (Type(p)==BODY1) {
				RcutSq=gdforce->RcutSq11; 
				ssq=gdforce->ssq11; fphi=gdforce->fphi11; fa=gdforce->fa11;
				vc = gdforce->vc11; dvc = gdforce->dvc11;
				DOTPSUBV(dr2, dr, Pos(p), Pos(qb));
//				VWrapAll_ptr(dr);
				DOTVP(dr2, dr, dr);
				if (dr2<RcutSq) {
//					ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
//						fa, ssq, vc, dvc);

					rri=ssq/dr2; rri3=rri*rri*rri;
					fPot = (rri3-1.0)*rri3;
					fAcc = rri3*(rri3-0.5)*rri;

					ADDMULVS(Acc11(p), dr, fAcc);
					ADDMULVS(Acc11(qb), dr, -fAcc);
					Phi11(p) += fPot;
					Phi11(qb) += fPot;
					gdforce->virSum11 += fAcc*dr2;
					gdforce->nbbcalc += 2;
//					if (gdforce->computeTransport)
//						TransportCalc2(p, qb, dr, fPot, fAcc);
				}
			} else {
				RcutSq=gdforce->RcutSq22; 
				ssq=gdforce->ssq22; fphi=gdforce->fphi22; fa=gdforce->fa22;
				vc = gdforce->vc22; dvc = gdforce->dvc22;
				DOTPSUBV(dr2, dr, Pos(p), Pos(qb));
//				VWrapAll_ptr(dr);
				DOTVP(dr2, dr, dr);
				if (dr2<RcutSq) {
//					ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
//						fa, ssq, vc, dvc);

					rri=ssq/dr2; rri3=rri*rri*rri;
					fPot = (rri3-1.0)*rri3;
					fAcc = rri3*(rri3-0.5)*rri;

					ADDMULVS(Acc22(p), dr, fAcc);
					ADDMULVS(Acc22(qb), dr, -fAcc);
					Phi22(p) += fPot;
					Phi22(qb) += fPot;
					gdforce->virSum22 += fAcc*dr2;
					gdforce->nbbcalc += 2;
//					if (gdforce->computeTransport)
//						TransportCalc2(p, qb, dr, fPot, fAcc);
				}
			}
		} else {
			RcutSq=gdforce->RcutSq12; 
			ssq=gdforce->ssq12; fphi=gdforce->fphi12; fa=gdforce->fa12;
			vc = gdforce->vc12; dvc = gdforce->dvc12;
			DOTPSUBV(dr2, dr, Pos(p), Pos(qb));
//			VWrapAll_ptr(dr);
			DOTVP(dr2, dr, dr);
			if (dr2<RcutSq) {
//				ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
//					fa, ssq, vc, dvc);

				rri=ssq/dr2; rri3=rri*rri*rri;
				fPot = (rri3-1.0)*rri3;
				fAcc = rri3*(rri3-0.5)*rri;

				ADDMULVS(Acc12(p), dr, fAcc);
				ADDMULVS(Acc12(qb), dr, -fAcc);
				Phi12(p) += fPot;
				Phi12(qb) += fPot;
				gdforce->virSum12 += fAcc*dr2;
				gdforce->nbbcalc += 2;
//				if (gdforce->computeTransport)
//					TransportCalc2(p, qb, dr, fPot, fAcc);
			}
		}
	}
}

// TERMINA METODO NORMAL3 DE CALCULO DE LA FUERZA


#  define FACTIVE  0.75                         
local int actlen;
local int *activenb;
local int nblist;


/* COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE LISTA DE VECINOS */

void tljforcecalc_nblist(bodyptr btab, int nbody, 
			global_data_tree_tljforcecalc *gdforce, global_data_tree *gdtree)
{
    bodyptr p;
    double cpustart;
	vector pos0, acc0;
	real phi0;
 
    cpustart = cputime();                       

    actlen = FACTIVE * 216 * gdtree->tdepth;
    activenb = (int *) allocate(actlen * sizeof(int));

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    gdforce->nbbcalc = gdforce->nbccalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

    for (p = btab; p < btab+nbody; p++) {
		SETV(pos0, Pos(p));
		phi0 = 0.0;
		CLRV(acc0);
		nblist=0;
		normal_walktree_nblist(p, ((nodeptr) gdtree->root), gdtree->rsize, 
								pos0,&phi0,acc0, gdforce);
		int_piksrt(nblist, activenb);
		sumnode_nblist_02(p, pos0, &phi0, acc0, gdforce);
		ADDMULVS(Acc(p), acc0, 1.0/Mass(p));
		Phi(p) += phi0; gdforce->uSum += Mass(p)*Phi(p);
	}
	gdforce->virSum = 0.5*gdforce->virSum;
	gdforce->uSum = 0.5*gdforce->uSum;

	free(activenb);
    gdforce->cpuforce = cputime() - cpustart;
}

#undef FACTIVE

local void normal_walktree_nblist(bodyptr p, nodeptr q, real qsize, 
							vector pos0, real *phi0, vector acc0, 
							global_data_tree_tljforcecalc *gdforce)
{
    nodeptr l;
	real drpq, drpq2, rcell;
    vector dr;
	int k;
	real Rcut;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
//				for (k = 0; k < NDIM; k++)
//					dr[k]=dr[k] -
//						((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
				DOTVP(drpq2, dr, dr);
				drpq = rsqrt(drpq2);
				rcell= qsize * rsqrt((real)(NDIM))/2.0;

				if (Type(p) == BODY1 )
					Rcut = gdforce->Rcut11Max;
				else
					Rcut = gdforce->Rcut22Max;

                if ( drpq >= Rcut+rcell ) { 
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree_nblist(p,l,qsize/2,pos0,
												phi0,acc0,gdforce);
					}
				}
			} else {
					sumnode_nblist_01(p, ((cellptr) q),( (cellptr) q+1),
									pos0,phi0,acc0, gdforce);
			}
		}
	}
}

local void sumnode_nblist_01(bodyptr p0, cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0, 
					global_data_tree_tljforcecalc *gdforce)
{
    cellptr p;
    real dr2, phi_p, mr3i;
    vector dr;
	real rri, rri3;
	int k, ip;
	real RcutSq;
 
    for (p = start; p < finish; p++) {          
		DOTPSUBV(dr2, dr, pos0, Pos(p));
//		for (k = 0; k < NDIM; k++)
//			dr[k]=dr[k] - ((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
		DOTVP(dr2, dr, dr);

		if ( Type(p0) == Type(p) ) {
			if (Type(p0)==BODY1) {
				RcutSq=gdforce->RcutSq11; 
			} else {
				RcutSq=gdforce->RcutSq22;
			}
		} else {
			RcutSq=gdforce->RcutSq12;
		}

		if (dr2<RcutSq) {
			gdforce->nbbcalc += 1;
			ip = (bodyptr)p-bodytab;
			activenb[nblist]=ip;
			nblist +=1;
		}
    }
}

local void sumnode_nblist_02(bodyptr p0, vector pos0, real *phi0, vector acc0, 
							global_data_tree_tljforcecalc *gdforce)
{
    bodyptr p;
    real dr2, fPot, fAcc;
    vector dr;
//	real rri, rri3;
//	int k, ip;
	int ip;
	real RcutSq, ssq, fphi, fa, vc, dvc;
//	real RcutSq;
//	matrix w;
 
    for (ip = 0; ip < nblist; ip++) {
		p = bodytab + activenb[ip];

		bljFactors(p0, p, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

		DOTPSUBV(dr2, dr, pos0, Pos(p));
//		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

/*
		if ( Type(p0) == Type(p) ) {
			if (Type(p)==BODY1) {
				ssq=gdforce->ssq11; fphi=gdforce->fphi11; fa=gdforce->fa11;
			} else {
				ssq=gdforce->ssq22; fphi=gdforce->fphi22; fa=gdforce->fa22;
			}
		} else {
			ssq=gdforce->ssq12; fphi=gdforce->fphi12; fa=gdforce->fa12;
		}

		DOTPSUBV(dr2, dr, pos0, Pos(p));
		for (k = 0; k < NDIM; k++)
			dr[k]=dr[k] - ((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
		DOTVP(dr2, dr, dr);
*/
		ForcePotentialCalc(p0, p, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
			fa, ssq, vc, dvc);
/*
		rri=ssq/dr2; rri3=rri*rri*rri;
		phi_p = fphi*(rri3-1.0)*rri3;
		mr3i = fa*rri3*(rri3-0.5)*rri;
*/
		*phi0 += fPot;
		ADDMULVS(acc0, dr, fAcc);
		gdforce->virSum += fAcc*dr2;
		if (gdforce->computeTransport)
			TransportCalc(p0, dr, fPot, fAcc);
    }
}

/* TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE LISTA DE VECINOS */


/* COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CALCULO DIRECTO DE LA FUERZA */

void ind_tljforcecalc_direct(bodyptr btab, int nbody, 
						global_data_tree_tljforcecalc *gdforce, 
						bodyptr p, real *uSum, real *virSum, short *flag)
{
    bodyptr q;
    double cpustart;
    vector dr;
	real rri, rri3, fPot, fAcc, dr2;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	real dr2min=0.5;
 
    cpustart = cputime();

	Phi(p) = 0.0;
	CLRV(Acc(p));

    gdforce->nbbcalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	*flag = 1;

	DO_BODY(q, btab, btab+nbody) {
		if (p==q) continue;

		bljFactors(p, q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
//		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		if (dr2 < dr2min) {
			*flag = 0;
			break;
		}

		if (dr2 < RcutSq) {
			ForcePotentialCalc(p, q, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, vc, dvc);
/*
			rri=ssq/dr2; rri3=rri*rri*rri;
			fPot = fphi*(rri3-1.0)*rri3;
			fAcc = fa*rri3*(rri3-0.5)*rri;
*/
			Phi(p) += fPot;
			ADDMULVS(Acc(p), dr, fAcc);
			gdforce->virSum += fAcc*dr2;
			gdforce->nbbcalc += 1;
		}
	}
	DIVVS(Acc(p), Acc(p), Mass(p));
	gdforce->uSum += Phi(p);
	*virSum = gdforce->virSum;
	*uSum = gdforce->uSum;

    gdforce->cpuindforce = cputime() - cpustart;            
}

void tljforcecalc_direct(bodyptr btab, int nbody, 
						global_data_tree_tljforcecalc *gdforce)
{
    bodyptr p, q;
    double cpustart;
    vector dr;
//	real rri, rri3;
	real fPot, fAcc, dr2;
	real RcutSq, ssq, fphi, fa, vc, dvc;
//	matrix w;
 
    cpustart = cputime();

	DO_BODY(p, btab, btab+nbody) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    gdforce->nbbcalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

	DO_BODY(p, btab, btab+nbody) {
		DO_BODY(q, btab, btab+nbody) {
			if (p==q) continue;

			bljFactors(p, q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

			DOTPSUBV(dr2, dr, Pos(p), Pos(q));
//			VWrapAll_ptr(dr);
			DOTVP(dr2, dr, dr);

			if (dr2 < RcutSq) {
				ForcePotentialCalc(p, q, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
					fa, ssq, vc, dvc);
/*
				rri=ssq/dr2; rri3=rri*rri*rri;
				fPot = fphi*(rri3-1.0)*rri3;
				fAcc = fa*rri3*(rri3-0.5)*rri;
*/
				Phi(p) += fPot;
				ADDMULVS(Acc(p), dr, fAcc);
				gdforce->virSum += fAcc*dr2;
				gdforce->nbbcalc += 1;
				if (gdforce->computeTransport)
					TransportCalc(p, dr, fPot, fAcc);
			}
		}
		DIVVS(Acc(p), Acc(p), Mass(p));
		gdforce->uSum += Phi(p);
	}
	gdforce->virSum = 0.5*gdforce->virSum;
	gdforce->uSum = 0.5*gdforce->uSum;

    gdforce->cpuforce = cputime() - cpustart;            
}

void tljforcecalc_direct2(bodyptr btab, int nbody, 
						global_data_tree_tljforcecalc *gdforce)
{
    bodyptr p, q;
    double cpustart;
    vector dr;
//	int k;
//	real rri, rri3;
	real fPot, fAcc, dr2;
	real RcutSq, ssq, fphi, fa, vc, dvc;
//	matrix w;
 
    cpustart = cputime(); 

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    gdforce->nbbcalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

	DO_BODY(p, btab, btab+nbody-1) {
		DO_BODY(q, p+1, btab+nbody) {
			bljFactors(p, q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

			DOTPSUBV(dr2, dr, Pos(p), Pos(q));
//			VWrapAll_ptr(dr);
			DOTVP(dr2, dr, dr);

			if (dr2<RcutSq) {
				ForcePotentialCalc(p, q, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
					fa, ssq, vc, dvc);
/*
				rri=ssq/dr2; rri3=rri*rri*rri;
				fPot = fphi*(rri3-1.0)*rri3;
				fAcc = fa*rri3*(rri3-0.5)*rri;
*/
				Phi(p) += fPot; Phi(q) += fPot;
				ADDMULVS(Acc(p), dr, fAcc);
				ADDMULVS(Acc(q), dr, -fAcc);
				gdforce->uSum += fPot;
				gdforce->virSum += fAcc*dr2;
				gdforce->nbbcalc += 2;
				if (gdforce->computeTransport)
					TransportCalc2(p, q, dr, fPot, fAcc);
			}
		}
	}

	DO_BODY(p,btab,btab+nbody)
		DIVVS(Acc(p), Acc(p), Mass(p));

    gdforce->cpuforce = cputime() - cpustart;            
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE CALCULO DIRECTO DE LA FUERZA


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES PARA EL CALCULO DE LA FUERZA

#  define FACTIVE  2
local nodeptr *active;
local cellptr interact;

#if defined (DIAGNOSTICS)
#define INTERFNAME		"interactions"
local stream interout;
#endif

void tljforcecalc_barnes(bodyptr btab, int nbody, 
			global_data_tree_tljforcecalc *gdforce, global_data_tree *gdtree)
{
    double cpustart;
    vector rmid;
	bodyptr p;
 
#if defined (DIAGNOSTICS)
interout=fopen(INTERFNAME,"w!");
#endif

    actlen = FACTIVE * 216 * gdtree->tdepth;
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    interact = (cellptr) allocate(actlen * sizeof(cell));

    cpustart = cputime();                       
    gdforce->actmax = gdforce->nbbcalc = gdforce->nbccalc = 0;             
    active[0] = (nodeptr) gdtree->root;
    CLRV(rmid);

	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

    walktree_barnes(active, active + 1, interact, interact + actlen,
             (nodeptr) gdtree->root, gdtree->rsize, rmid, gdforce);      

	gdforce->virSum = 0.5*gdforce->virSum;
	gdforce->uSum = 0.5*gdforce->uSum;

    gdforce->cpuforce = cputime() - cpustart;            
    free(active);
    free(interact);
	
#if defined (DIAGNOSTICS)
fclose(interout);
#endif
}

#undef FACTIVE

local void walktree_barnes(nodeptr *aptr, nodeptr *nptr, cellptr cptr, 
					cellptr bptr, nodeptr p, real psize, vector pmid, 
					global_data_tree_tljforcecalc *gdforce)
{
    nodeptr *np, *ap, q;
    int actsafe;

    if (Update(p)) {                            
        np = nptr;                              
        actsafe = actlen - NSUB;                
        for (ap = aptr; ap < nptr; ap++)        
            if (Type(*ap) == CELL)
                if (accept_blj(p, *ap, gdforce)) {
				} else {
                    if (np - active >= actsafe)
                        error("walktree_barnes: active list overflow\n");
                    for (q = More(*ap); q != Next(*ap); q = Next(q))
						*np++= q;
				}
            else
                if (*ap != p) {
//					if (! accept_blj(p, *ap, gdforce) ) { 
						--bptr;
						Mass(bptr) = Mass(*ap);
						SETV(Pos(bptr), Pos(*ap));
						Id(bptr) = Id(*ap);
//					}
                }
        gdforce->actmax = MAX(gdforce->actmax, np - active);
        if (np != nptr)
            walksub(nptr, np, cptr, bptr, p, psize, pmid, gdforce);
        else {
            if (Type(p) == CELL)
                error("walktree_barnes: recursion terminated with cell\n");
            gravsum((bodyptr) p, cptr, bptr, gdforce);
        }
    }
}

local bool accept_blj(nodeptr p, nodeptr q, 
						global_data_tree_tljforcecalc *gdforce)
{
	real drpq, drpq2;
    vector dr;
	int k;
	vector posp, posq;

	if (Type(p)==CELL) {
		SETV(posp, GPos(p));
	} else
		SETV(posp, Pos(p));

	if (Type(q)==CELL) {
		SETV(posq, GPos(q));
	} else
		SETV(posq, Pos(q));

	DOTPSUBV(drpq2, dr, posp, posq);
//	DO_COORD(k)
//		dr[k]=dr[k] - ( (real) (nint(dr[k]/gdforce->Box[k])) )*gdforce->Box[k];
	DOTVP(drpq2, dr, dr);
	drpq = rsqrt(drpq2);

	if ( drpq >= Rcut(p)+Rcut(q) )
		return (TRUE);
	else
		return (FALSE);
}

local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid, 
				   global_data_tree_tljforcecalc *gdforce)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;
 
    poff = psize / 4;                           
    if (Type(p) == CELL) {                       
        for (q = More(p); q != Next(p); q = Next(q)) {
            for (k = 0; k < NDIM; k++)          
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walktree_barnes(nptr, np, cptr, bptr, q, psize / 2, nmid, gdforce);
        }
    } else {
        for (k = 0; k < NDIM; k++)              
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree_barnes(nptr, np, cptr, bptr, p, psize / 2, nmid, gdforce);
                                                
    }
}
 
local void gravsum(bodyptr p0, cellptr cptr, cellptr bptr, 
					global_data_tree_tljforcecalc *gdforce)
{
    vector pos0, acc0;
    real phi0;
 
    SETV(pos0, Pos(p0));
    phi0 = 0.0;
    CLRV(acc0);

    sumnode_barnes(p0, bptr, interact + actlen, pos0, &phi0, acc0, gdforce);

    Phi(p0) = phi0; gdforce->uSum += Mass(p0)*Phi(p0);
	DIVVS(Acc(p0), acc0, Mass(p0));
}

local void sumnode_barnes(bodyptr p0, cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0, 
				   global_data_tree_tljforcecalc *gdforce)
{
    cellptr p;
    real dr2, fPot, fAcc;
    vector dr;
//	real rri, rri3;
//	int k;
	real RcutSq, ssq, fphi,fa, vc, dvc;
//	matrix w;

#if defined (DIAGNOSTICS)
	vector acc;
#endif

#if defined (DIAGNOSTICS)
fprintf(interout,"\n[%d] :",Id(p0));
#endif

    for (p = start; p < finish; p++) {
		bljFactors(p0, (bodyptr)p, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

		DOTPSUBV(dr2, dr, pos0, Pos(p));
//		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			ForcePotentialCalc(p0, (bodyptr)p, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, vc, dvc);
/*
			rri=ssq/dr2; rri3=rri*rri*rri;
			phi_p = fphi*(rri3-1.0)*rri3;
			mr3i = fa*rri3*(rri3-0.5)*rri;
*/
			*phi0 += fPot;
			ADDMULVS(acc0, dr, fAcc);

#if defined (DIAGNOSTICS)
MULVS(acc, dr, mr3i);
fprintf(interout,"\n(%d, %d) : %g %g %g", Id(p0), Id(p), acc[0], acc[1], acc[2]);
#endif

			gdforce->virSum += fAcc*dr2;
			gdforce->nbbcalc += 1;
			if (gdforce->computeTransport)
				TransportCalc(p0, dr, fPot, fAcc);
		}
    }

#if defined (DIAGNOSTICS)
fprintf(interout,"\n");
#endif
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES PARA EL CALCULO DE LA FUERZA


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS PARA EL CALCULO DE LA FUERZA

void tljforcecalc_cellsmethod(bodyptr btab, int nbody, 
							global_data_tree_tljforcecalc *gdforce)
{
	bodyptr p, p1, p2;
//	matrix w;
//	int k;
	vector dr;
	vector invWid, rs, shift;
	vectorI cc, m1v, m2v, vOff[] = OFFSET_VALS;
//	real fcVal, rr, rri, rri3, uVal;
	real fPot, fAcc, rr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
	double cpustart;
	real RcutSq, ssq, fphi, fa, vc, dvc;

	cpustart = cputime();                       

	VDiv (invWid, gdforce->cells, gdforce->Box);
	for (n = nbody; n < nbody + VProd (gdforce->cells); n ++) 
		gdforce->cellList[n] = -1;

	DO_BODY(p,btab,btab+nbody) {
		n = p-btab;										// n indice de cuerpo...
		ADD2VS(rs, Pos(p), gdforce->Box, 0.5);
		VMul (cc, rs, invWid);
		c = VLinear (cc, gdforce->cells) + nbody;		// c indice de celda...
		gdforce->cellList[n] = gdforce->cellList[c];
		gdforce->cellList[c] = n;
	}

	DO_BODY(p,btab,btab+nbody) {CLRV(Acc(p)); Phi(p) = 0.0;}
	gdforce->uSum = 0.; gdforce->virSum = 0.; gdforce->nbbcalc = 0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

#if (NDIM==3)
	for (m1z = 0; m1z < gdforce->cells[2]; m1z ++) {
#endif
		for (m1y = 0; m1y < gdforce->cells[1]; m1y ++) {
			for (m1x = 0; m1x < gdforce->cells[0]; m1x ++) {
#if (NDIM==3)
				VSet (m1v, m1x, m1y, m1z);
#else
#if (NDIM==2)
				VSet (m1v, m1x, m1y);
#endif
#endif
				m1 = VLinear (m1v, gdforce->cells) + nbody;
				for (offset = 0; offset < N_OFFSET; offset ++) {
					ADDV(m2v, m1v, vOff[offset]);
					CLRV(shift);
//					VCellWrapAll_ptr ();
					m2 = VLinear (m2v, gdforce->cells) + nbody;
					DO_CELL_ptr (j1, m1) {
						DO_CELL_ptr (j2, m2) {
							if (m1 != m2 || j2 < j1) {
								p1 = btab+j1; p2 = btab+j2;
								bljFactors(p1,p2,&RcutSq,&ssq,&fphi,&fa,&vc,&dvc,gdforce);
								SUBV(dr, Pos(p1), Pos(p2));
								VVSub (dr, shift);
								rr = VLenSq (dr);
								if (rr < RcutSq) {

									ForcePotentialCalc(p1, p2, gdforce->potType, &fPot, 
										&fAcc, rr, fphi, fa, ssq, vc, dvc);
/*
									rri = ssq / rr;
									rri3 = Cube (rri);
									fcVal = fa * rri3 * (rri3 - 0.5) * rri;
									uVal = fphi * rri3 * (rri3 - 1.);
*/
									VVSAdd(Acc(p1), fAcc/Mass(p1), dr);
									VVSAdd (Acc(p2), - fAcc/Mass(p2), dr);
									Phi(p1) += fPot; Phi(p2) +=fPot;
									gdforce->uSum += fPot;
									gdforce->virSum += fAcc * rr;
									gdforce->nbbcalc+=2;
									if (gdforce->computeTransport) {
										TransportCalc2(p1, p2, dr, fPot, fAcc);
/*										en(p1) += uVal;
										en(p2) += uVal;
										OUTVP(w, dr, dr);
										MULMS(w, w, fcVal);
										DO_COORD(k) {
											VVAdd(rf(p1)[k], w[k]);
											VVAdd(rf(p2)[k], w[k]);
										} */
									}
								}
							}
						}
					}
				}
			}
		}
#if (NDIM==3)
	}
#endif

	gdforce->cpuforce = cputime() - cpustart;
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS PARA EL CALCULO DE LA FUERZA

// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS3 PARA EL CALCULO DE LA FUERZA

void tljforcecalc_cellsmethod3(bodyptr btab, int nbody, 
								global_data_tree_tljforcecalc *gdforce)
{
	bodyptr p, p1, p2;
	vector dr, invWid, rs, shift;
	vectorI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	real fcVal, rr, rri, rri3, uVal;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
	double cpustart;
	real RcutSq, ssq,fphi,fa;
//	matrix w;

	cpustart = cputime();                       

	VDiv (invWid, gdforce->cells, gdforce->Box);
	for (n = nbody; n < nbody + VProd (gdforce->cells); n ++) 
		gdforce->cellList[n] = -1;

	DO_BODY(p,btab,btab+nbody) {
		n = p-btab;										// n indice de cuerpo...
		ADD2VS(rs, Pos(p), gdforce->Box, 0.5);
		VMul (cc, rs, invWid);
		c = VLinear (cc, gdforce->cells) + nbody;		// c indice de celda...
		gdforce->cellList[n] = gdforce->cellList[c];
		gdforce->cellList[c] = n;
	}

	DO_BODY(p,btab,btab+nbody) {
		CLRV(Acc(p)); Phi(p) = 0.0;
		CLRV(Acc11(p)); Phi11(p) = 0.0;
		CLRV(Acc12(p)); Phi12(p) = 0.0;
		CLRV(Acc22(p)); Phi22(p) = 0.0;
	}
	gdforce->uSum = 0.;
	gdforce->virSum = 0.;
	gdforce->virSum11 = 0.;
	gdforce->virSum12 = 0.;
	gdforce->virSum22 = 0.;
	gdforce->nbbcalc = 0;
//	if (gdforce->computeTransport) {
//		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
//	}

#if (NDIM==3)
	for (m1z = 0; m1z < gdforce->cells[2]; m1z ++) {
#endif
		for (m1y = 0; m1y < gdforce->cells[1]; m1y ++) {
			for (m1x = 0; m1x < gdforce->cells[0]; m1x ++) {
#if (NDIM==3)
				VSet (m1v, m1x, m1y, m1z);
#else
#if (NDIM==2)
				VSet (m1v, m1x, m1y);
#endif
#endif
				m1 = VLinear (m1v, gdforce->cells) + nbody;
				for (offset = 0; offset < N_OFFSET; offset ++) {
					ADDV(m2v, m1v, vOff[offset]);
					CLRV(shift);
//					VCellWrapAll_ptr ();
					m2 = VLinear (m2v, gdforce->cells) + nbody;
					DO_CELL_ptr (j1, m1) {
						DO_CELL_ptr (j2, m2) {
							if (m1 != m2 || j2 < j1) {
								p1 = btab+j1; p2 = btab+j2;

								if ( Type(p1) == Type(p2) ) {
									if (Type(p1)==BODY1) {
										RcutSq=gdforce->RcutSq11; 
										ssq=gdforce->ssq11;
/*
										fphi=gdforce->fphi11;
										fa=gdforce->fa11;
										vc=gdforce->vc11;
										dvc=gdforce->dvc11;
*/
										SUBV(dr, Pos(p1), Pos(p2));
										VVSub (dr, shift);
										rr = VLenSq (dr);
										if (rr < RcutSq) {
//											ForcePotentialCalc(gdforce->potType, &uVal, 
//												&fcVal, rr, fphi, fa, ssq, vc, dvc);

											rri = ssq / rr;
											rri3 = Cube (rri);
											fcVal = rri3 * (rri3 - 0.5) * rri;
											uVal = rri3 * (rri3 - 1.);

											VVSAdd(Acc11(p1), fcVal, dr);
											VVSAdd (Acc11(p2), - fcVal, dr);
											Phi11(p1) += uVal; Phi11(p2) +=uVal;
											gdforce->virSum11 += fcVal * rr;
											gdforce->nbbcalc+=2;
//											if (gdforce->computeTransport)
//												TransportCalc2(p1, p2, dr, uVal, fcVal);
										}
									} else {
										RcutSq=gdforce->RcutSq22; 
										ssq=gdforce->ssq22; 
/*
										fphi=gdforce->fphi22;
										fa=gdforce->fa22;
										vc=gdforce->vc22;
										dvc=gdforce->dvc22;
*/
										SUBV(dr, Pos(p1), Pos(p2));
										VVSub (dr, shift);
										rr = VLenSq (dr);
										if (rr < RcutSq) {
//											ForcePotentialCalc(gdforce->potType, &uVal, 
//												&fcVal, rr, fphi, fa, ssq, vc, dvc);

											rri = ssq / rr;
											rri3 = Cube (rri);
											fcVal = rri3 * (rri3 - 0.5) * rri;
											uVal = rri3 * (rri3 - 1.);

											VVSAdd(Acc22(p1), fcVal, dr);
											VVSAdd (Acc22(p2), - fcVal, dr);
											Phi22(p1) += uVal; Phi22(p2) +=uVal;
											gdforce->virSum22 += fcVal * rr;
											gdforce->nbbcalc+=2;
//											if (gdforce->computeTransport)
//												TransportCalc2(p1, p2, dr, uVal, fcVal);
										}
									}
								} else {
									RcutSq=gdforce->RcutSq12; 
									ssq=gdforce->ssq12; 
/*
									fphi=gdforce->fphi12;
									fa=gdforce->fa12;
									vc=gdforce->vc12;
									dvc=gdforce->dvc12;
*/
									SUBV(dr, Pos(p1), Pos(p2));
									VVSub (dr, shift);
									rr = VLenSq (dr);
									if (rr < RcutSq) {
//										ForcePotentialCalc(gdforce->potType, &uVal, 
//											&fcVal, rr, fphi, fa, ssq, vc, dvc);

										rri = ssq / rr;
										rri3 = Cube (rri);
										fcVal = rri3 * (rri3 - 0.5) * rri;
										uVal = rri3 * (rri3 - 1.);

										VVSAdd(Acc12(p1), fcVal, dr);
										VVSAdd (Acc12(p2), - fcVal, dr);
										Phi12(p1) += uVal; Phi12(p2) +=uVal;
										gdforce->virSum12 += fcVal * rr;
										gdforce->nbbcalc+=2;
//										if (gdforce->computeTransport)
//											TransportCalc2(p1, p2, dr, uVal, fcVal);
									}
								}
							}
						}
					}
				}
			}
		}
#if (NDIM==3)
	}
#endif

	DO_BODY(p, btab, btab+nbody) {
		ADDMULVS(Acc(p), Acc11(p), gdforce->fa11);
		ADDMULVS(Acc(p), Acc12(p), gdforce->fa12);
		ADDMULVS(Acc(p), Acc22(p), gdforce->fa22);
		Phi(p) += gdforce->fphi11*Phi11(p)+gdforce->fphi12*Phi12(p)
					+gdforce->fphi22*Phi22(p);
		gdforce->uSum += Phi(p);
		DIVVS(Acc(p), Acc(p), Mass(p));
	}
	gdforce->uSum *= 0.5;
	gdforce->virSum = gdforce->fa11*gdforce->virSum11
							+gdforce->fa12*gdforce->virSum12
							+gdforce->fa22*gdforce->virSum22;

  gdforce->cpuforce = cputime() - cpustart;
}

/* TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS3 PARA EL CALCULO DE LA FUERZA */


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS PARA EL CALCULO DEL POTENCIAL

void tljpotcalc_cellsmethod(bodyptr btab, int nbody, 
							global_data_tree_tljforcecalc *gdforce)
{
	bodyptr p, p1, p2;
	vector dr, invWid, rs, shift;
	vectorI cc, m1v, m2v, vOff[] = OFFSET_VALS;
//	real fcVal; 
	real rr, fPot;
//	real rri, rri3;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
	double cpustart;
	real RcutSq, ssq, fphi,	fa, vc, dvc;

	cpustart = cputime();                       

	VDiv (invWid, gdforce->cells, gdforce->Box);
	for (n = nbody; n < nbody + VProd (gdforce->cells); n ++) 
		gdforce->cellList[n] = -1;

	DO_BODY(p,btab,btab+nbody) {
		n = p-btab;										// n indice de cuerpo...
		ADD2VS(rs, Pos(p), gdforce->Box, 0.5);
		VMul (cc, rs, invWid);
		c = VLinear (cc, gdforce->cells) + nbody;		// c indice de celda...
		gdforce->cellList[n] = gdforce->cellList[c];
		gdforce->cellList[c] = n;
	}

	DO_BODY(p,btab,btab+nbody) {Phi(p) = 0.0;}
	gdforce->uSum = 0.; gdforce->nbbcalc = 0;

#if (NDIM==3)
	for (m1z = 0; m1z < gdforce->cells[2]; m1z ++) {
#endif
		for (m1y = 0; m1y < gdforce->cells[1]; m1y ++) {
			for (m1x = 0; m1x < gdforce->cells[0]; m1x ++) {
#if (NDIM==3)
				VSet (m1v, m1x, m1y, m1z);
#else
#if (NDIM==2)
				VSet (m1v, m1x, m1y);
#endif
#endif
				m1 = VLinear (m1v, gdforce->cells) + nbody;
				for (offset = 0; offset < N_OFFSET; offset ++) {
					ADDV(m2v, m1v, vOff[offset]);
					CLRV(shift);
//					VCellWrapAll_ptr ();
					m2 = VLinear (m2v, gdforce->cells) + nbody;
					DO_CELL_ptr (j1, m1) {
						DO_CELL_ptr (j2, m2) {
							if (m1 != m2 || j2 < j1) {
								p1 = btab+j1; p2 = btab+j2;
								bljFactors(p1,p2,&RcutSq,&ssq,&fphi,&fa,&vc,&dvc,gdforce);
								SUBV(dr, Pos(p1), Pos(p2));
								VVSub (dr, shift);
								rr = VLenSq (dr);
								if (rr < RcutSq) {
									PotentialCalc(p1, p2, gdforce->potType, &fPot, 
										rr, fphi, ssq, vc, dvc);
/*
									rri = ssq / rr;
									rri3 = Cube (rri);
									uVal = fphi * rri3 * (rri3 - 1.);
*/
//									fcVal = fa * rri3 * (rri3 - 0.5) * rri;
//									VVSAdd(Acc(p1), fcVal/Mass(p1), dr);
//									VVSAdd (Acc(p2), - fcVal/Mass(p2), dr);
									Phi(p1) += fPot; Phi(p2) +=fPot;
									gdforce->uSum += fPot;
//									gdforce->virSum += fcVal * rr;
									gdforce->nbbcalc+=2;
								}
							}
						}
					}
				}
			}
		}
#if (NDIM==3)
	}
#endif

	gdforce->cpupot = cputime() - cpustart;
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS PARA EL CALCULO DEL POTENCIAL

// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES2 PARA EL CALCULO DE LA FUERZA

#define FACTIVE	2
local bodyptr *interactbody;
local cellptr *interactcell;

void tljforcecalc_barnes2(bodyptr btab, int nbody, 
			global_data_tree_tljforcecalc *gdforce, global_data_tree *gdtree)
{
    double cpustart;
    vector rmid;
	bodyptr p;

#if defined (DIAGNOSTICS)
interout=fopen(INTERFNAME,"w!");
#endif

    actlen = FACTIVE * 216 * gdtree->tdepth;
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    interactcell = (cellptr *) allocate(actlen * sizeof(cellptr));

    interactbody = (bodyptr *) allocate(actlen * sizeof(bodyptr));

    cpustart = cputime();                       
    gdforce->actmax = gdforce->nbbcalc = gdforce->nbccalc = 0;             
    active[0] = (nodeptr) gdtree->root;
    CLRV(rmid);

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	if (gdforce->computeTransport) {
		DO_BODY(p,btab,btab+nbody) {CLRM(rf(p)); en(p) = 0.0;}
	}

    walktree_barnes2(active, active + 1, interactcell, interactbody + actlen,
             (nodeptr) gdtree->root, gdtree->rsize, rmid, gdforce);      

	DO_BODY(p,btab,btab+nbody) {
		DIVVS(Acc(p), Acc(p), Mass(p));
	}

    gdforce->cpuforce = cputime() - cpustart;            
    free(active);
    free(interact);
	
#if defined (DIAGNOSTICS)
fclose(interout);
#endif
}

#undef FACTIVE

local void walktree_barnes2(nodeptr *aptr, nodeptr *nptr, cellptr *cptr, 
					bodyptr *bptr, nodeptr p, real psize, vector pmid, 
					global_data_tree_tljforcecalc *gdforce)
{
    nodeptr *np, *ap, q;
    int actsafe;

    if (Update(p)) {
        np = nptr;                              
        actsafe = actlen - NSUB;                
        for (ap = aptr; ap < nptr; ap++)        
            if (Type(*ap) == CELL) {
                if (accept_blj(p, *ap, gdforce)) { 
                } else {
                    if (np - active >= actsafe)
                        error("walktree_barnes2: active list overflow\n");
                    for (q = More(*ap); q != Next(*ap); q = Next(q))
						*np++= q;
                }
            } else
                if (*ap != p) {
					*--bptr = (bodyptr)*ap;
                }
        gdforce->actmax = MAX(gdforce->actmax, np - active);
        if (np != nptr)
            walksub2(nptr, np, cptr, bptr, p, psize, pmid, gdforce);

        else {
            if (Type(p) == CELL)
                error("walktree_barnes2: recursion terminated with cell\n");
			sumnode_barnes2((bodyptr)p, bptr, interactbody + actlen, gdforce);
        }
    }
}

local void walksub2(nodeptr *nptr, nodeptr *np, cellptr *cptr, bodyptr *bptr,
                   nodeptr p, real psize, vector pmid, 
				   global_data_tree_tljforcecalc *gdforce)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;
 
    poff = psize / 4;                           
    if (Type(p) == CELL) {                       
        for (q = More(p); q != Next(p); q = Next(q)) {
            for (k = 0; k < NDIM; k++)          
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walktree_barnes2(nptr, np, cptr, bptr, q, psize / 2, nmid, gdforce);
        }
    } else {
        for (k = 0; k < NDIM; k++)              
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree_barnes2(nptr, np, cptr, bptr, p, psize / 2, nmid, gdforce);

    }
}

local void sumnode_barnes2(bodyptr p, bodyptr *start, bodyptr *finish,
				   global_data_tree_tljforcecalc *gdforce)
{
    cellptr c;
    real dr2, fPot, fAcc;
    vector dr;
//	real rri, rri3;
	int k;
	real RcutSq, ssq, fphi,fa, vc, dvc;
	bodyptr *q;
//	matrix w;

#if defined (DIAGNOSTICS)
	vector acc;
#endif

#if defined (DIAGNOSTICS)
fprintf(interout,"\n[%d] :",Id(p));
fprintf(interout,"\nBEFORE: Acc(%d) : %g %g %g", Id(p), Acc(p)[0], Acc(p)[1], Acc(p)[2]);
#endif

	DO_BODY(q, start, finish) {
		if ( Id(p) < Id(*q) ) {
			bljFactors(p, *q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, gdforce);

			DOTPSUBV(dr2, dr, Pos(p), Pos(*q));
//			VWrapAll_ptr(dr);
			DOTVP(dr2, dr, dr);

			if (dr2<RcutSq) {
#if defined (DIAGNOSTICS)
fprintf(interout,"\nBEFORE: Acc(%d) : %g %g %g", Id(*q), Acc(*q)[0], Acc(*q)[1], Acc(*q)[2]);
#endif
				ForcePotentialCalc(p, *q, gdforce->potType, &fPot, &fAcc, dr2, fphi, 
					fa, ssq, vc, dvc);
/*
				rri=ssq/dr2; rri3=rri*rri*rri;
				fPot = fphi*(rri3-1.0)*rri3;
				fAcc = fa*rri3*(rri3-0.5)*rri;
*/
				ADDMULVS(Acc(p), dr, fAcc);

#if defined (DIAGNOSTICS)
MULVS(acc, dr, fAcc);
fprintf(interout,"\n(%d, %d) : %g %g %g", Id(p), Id(*q), acc[0], acc[1], acc[2]);
#endif

				ADDMULVS(Acc(*q), dr, -fAcc);

#if defined (DIAGNOSTICS)
MULVS(acc, dr, -fAcc);
fprintf(interout,"\n(%d, %d) : %g %g %g", Id(*q), Id(p), acc[0], acc[1], acc[2]);
#endif
#if defined (DIAGNOSTICS)
fprintf(interout,"\nAFTER: Acc(%d) : %g %g %g\n", Id(*q), Acc(*q)[0], Acc(*q)[1], Acc(*q)[2]);
#endif
				Phi(p) += fPot;
				Phi(*q) += fPot;
				gdforce->uSum += fPot;
				gdforce->virSum += fAcc*dr2;
				gdforce->nbbcalc += 2;
				if (gdforce->computeTransport)
					TransportCalc2(p, *q, dr, fPot, fAcc);
			}
		}
	}

#if defined (DIAGNOSTICS)
fprintf(interout,"\nAFTER: Acc(%d) : %g %g %g", Id(p), Acc(p)[0], Acc(p)[1], Acc(p)[2]);
fprintf(interout,"\n");
#endif
}

#undef INTERFNAME

/* TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES2 PARA EL CALCULO DE LA FUERZA */


// COMIENZA RUTINAS AUXILIARES PARA EL CALCULO DE LA FUERZA ...

local void bljFactors(bodyptr p, bodyptr q, real *RcutSq, real *ssq, real *fphi,
		real *fa, real *vc, real *dvc, global_data_tree_tljforcecalc *gdforce)
{
	if ( Type(p) == Type(q) ) {
		if (Type(p)==BODY1) {
			*RcutSq=gdforce->RcutSq11; 
			*ssq=gdforce->ssq11; 
			*fphi=gdforce->fphi11; 
			*fa=gdforce->fa11;
			*vc=gdforce->vc11;
			*dvc=gdforce->dvc11;
		} else {
			if (Type(p)==BODY2) {
				*RcutSq=gdforce->RcutSq22; 
				*ssq=gdforce->ssq22; 
				*fphi=gdforce->fphi22; 
				*fa=gdforce->fa22;
				*vc=gdforce->vc22;
				*dvc=gdforce->dvc22;
			} else {
				*RcutSq=gdforce->RcutSq33; 
				*ssq=gdforce->ssq33;
				*fphi=gdforce->fphi33; 
				*fa=gdforce->fa33;
				*vc=gdforce->vc33;
				*dvc=gdforce->dvc33;
			}
		}
	} else {
		if ( (Type(p)==BODY1 && Type(q)==BODY2) || (Type(p)==BODY2 && Type(q)==BODY1) ) {
			*RcutSq=gdforce->RcutSq12; 
			*ssq=gdforce->ssq12; 
			*fphi=gdforce->fphi12; 
			*fa=gdforce->fa12;
			*vc=gdforce->vc12;
			*dvc=gdforce->dvc12;
		} else
			if (Type(p)==BODY1 && Type(q)==STATICBODY) {
				*RcutSq=gdforce->RcutSq13; 
				*ssq=gdforce->ssq13; 
				*fphi=gdforce->fphi13; 
				*fa=gdforce->fa13;
				*vc=gdforce->vc13;
				*dvc=gdforce->dvc13;
			} else {
				*RcutSq=gdforce->RcutSq23; 
				*ssq=gdforce->ssq23; 
				*fphi=gdforce->fphi23; 
				*fa=gdforce->fa23;
				*vc=gdforce->vc23;
				*dvc=gdforce->dvc23;
			}
	}
}

// COMIENZA RUTINAS PARA EL CALCULO DE LA FUERZA Y EL POTENCIAL ================

local void ForcePotentialCalc(bodyptr p, bodyptr q, int type, real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real vc, real dvc)
{
    switch(type) {
        case 0:
			ForcePotentialCalc_LJ(fPot, fAcc, dr2, fphi, fa, ssq, vc, dvc);
			break;
        default:
			ForcePotentialCalc_LJ(fPot, fAcc, dr2, fphi, fa, ssq, vc, dvc);
			break;
        case 1:
			ForcePotentialCalc_LN(fPot, fAcc, dr2, fphi, fa, ssq, vc, dvc);
			break;
        case 2:
			ForcePotentialCalc_Lorentz(p, q, fPot, fAcc, dr2, fphi, fa, ssq, vc, dvc);
			break;
    }
}

local void PotentialCalc(bodyptr p, bodyptr q, int type, real *fPot, 
	real dr2, real fphi, real ssq, real vc, real dvc)
{
    switch(type) {
        case 0:
			PotentialCalc_LJ(fPot, dr2, fphi, ssq, vc, dvc);
			break;
        default:
			PotentialCalc_LJ(fPot, dr2, fphi, ssq, vc, dvc);
			break;
        case 1:
			PotentialCalc_LN(fPot, dr2, fphi, ssq, vc, dvc);
        case 2:
			PotentialCalc_Lorentz(p, q, fPot, dr2, fphi, ssq, vc, dvc);
			break;
    }
}

local real PotentialFunction(int type, real r, real fphi, real sigma, real eps)
{
	real fPot;

    switch(type) {
        case 0:
			fPot=PotentialFunction_LJ(r, fphi, sigma, eps);
			break;
        default:
			fPot=PotentialFunction_LJ(r, fphi, sigma, eps);
			break;
        case 1:
			fPot=PotentialFunction_LN(r, fphi, sigma, eps);
        case 2:
			fPot=PotentialFunction_Lorentz(r, fphi, sigma, eps);
			break;
    }
	return (fPot);
}

local real DPotentialFunction(int type, real r, real fphi, real sigma, real eps)
{
	real dfPot;

    switch(type) {
        case 0:
			dfPot=DPotentialFunction_LJ(r, fphi, sigma, eps);
			break;
        default:
			dfPot=DPotentialFunction_LJ(r, fphi, sigma, eps);
			break;
        case 1:
			dfPot=DPotentialFunction_LN(r, fphi, sigma, eps);
			break;
        case 2:
			dfPot=DPotentialFunction_Lorentz(r, fphi, sigma, eps);
			break;
    }
	return (dfPot);
}

void PotentialParameters_t(int type, 
	real sigma11, real sigma12, real sigma13, real sigma22, real sigma23, real sigma33,
	real eps11, real eps12, real eps13, real eps22, real eps23, real eps33, 
	real Rcut11, real Rcut12, real Rcut13, real Rcut22, real Rcut23, real Rcut33, 
	global_data_tree_tljforcecalc *gdforce)
{
    switch(type) {
        case 0:
			PotentialParameters_LJ(sigma11, sigma12, sigma13, sigma22, sigma23, sigma33,
				eps11, eps12, eps13, eps22, eps23, eps33, 
				Rcut11, Rcut12, Rcut13, Rcut22, Rcut23, Rcut33, gdforce);
			break;
        default:
			PotentialParameters_LJ(sigma11, sigma12, sigma13, sigma22, sigma23, sigma33,
				eps11, eps12, eps13, eps22, eps23, eps33, 
				Rcut11, Rcut12, Rcut13, Rcut22, Rcut23, Rcut33, gdforce);
			break;
        case 1:
			PotentialParameters_LN(sigma11, sigma12, sigma13, sigma22, sigma23, sigma33,
				eps11, eps12, eps13, eps22, eps23, eps33, 
				Rcut11, Rcut12, Rcut13, Rcut22, Rcut23, Rcut33, gdforce);
			break;
        case 2:
			PotentialParameters_Lorentz(sigma11, sigma12, sigma13,
				sigma22, sigma23, sigma33,
				eps11, eps12, eps13, eps22, eps23, eps33, 
				Rcut11, Rcut12, Rcut13, Rcut22, Rcut23, Rcut33, gdforce);
			break;
    }
}

// COMIENZAN RUTINAS PARA EL POTENCIAL LENNARD-JONES ---------------------------
local void ForcePotentialCalc_LJ(real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real vc, real dvc)
{
	real rri, rri3;

	rri=ssq/dr2; rri3=rri*rri*rri;
	*fPot = fphi*(rri3-1.0)*rri3;
	*fAcc = fa*rri3*(rri3-0.5)*rri;
}

local void PotentialCalc_LJ(real *fPot, 
	real dr2, real fphi, real ssq, real vc, real dvc)
{
	real rri, rri3;

	rri=ssq/dr2; rri3=rri*rri*rri;
	*fPot = fphi*(rri3-1.0)*rri3;
}

local real PotentialFunction_LJ(real r, real fphi, real sigma, real eps)
{
	real s2, r2, rri, rri3, fPot;

	s2 = sigma*sigma;
	r2 = r*r;
	rri=s2/r2; rri3=rri*rri*rri;
	fPot = fphi*(rri3-1.0)*rri3;
	return (fPot);
}

local real DPotentialFunction_LJ(real r, real fa, real sigma, real eps)
{
	real s2, r2, rri, rri3, dfPot;

	s2 = sigma*sigma;
	r2 = r*r;
	rri=s2/r2; rri3=rri*rri*rri;
	dfPot = -fa*(rri3-0.5)*rri3*rri*r;
	return (dfPot);
}

local void PotentialParameters_LJ(real sigma11, real sigma12, real sigma13, 
	real sigma22, real sigma23, real sigma33,
	real eps11, real eps12, real eps13, real eps22, real eps23, real eps33, 
	real Rcut11, real Rcut12, real Rcut13, real Rcut22, real Rcut23, real Rcut33, 
	global_data_tree_tljforcecalc *gdforce)
{

	gdforce->fphi11 = 4.0*eps11;
	gdforce->ssq11 = sigma11*sigma11;
	gdforce->fa11 = 48.0*eps11/gdforce->ssq11;
	gdforce->fphi12 = 4.0*eps12;
	gdforce->ssq12 = sigma12*sigma12;
	gdforce->fa12 = 48.0*eps12/gdforce->ssq12;

	gdforce->fphi13 = 4.0*eps13;
	gdforce->ssq13 = sigma13*sigma13;
	gdforce->fa13 = 48.0*eps13/gdforce->ssq13;

	gdforce->fphi22 = 4.0*eps22;
	gdforce->ssq22 = sigma22*sigma22;
	gdforce->fa22 = 48.0*eps22/gdforce->ssq22;

	gdforce->fphi23 = 4.0*eps23;
	gdforce->ssq23 = sigma23*sigma23;
	gdforce->fa23 = 48.0*eps23/gdforce->ssq23;

	gdforce->fphi33 = 4.0*eps33;
	gdforce->ssq33 = sigma33*sigma33;
	gdforce->fa33 = 48.0*eps33/gdforce->ssq33;

	gdforce->vc11 = PotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fphi11, sigma11, eps11);
	gdforce->dvc11 = DPotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fa11, sigma11, eps11);
	gdforce->vc12 = PotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fphi12, sigma12, eps12);
	gdforce->dvc12 = DPotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fa12, sigma12, eps12);

	gdforce->vc13 = PotentialFunction(gdforce->potType, Rcut13, 
					gdforce->fphi13, sigma13, eps13);
	gdforce->dvc13 = DPotentialFunction(gdforce->potType, Rcut13, 
					gdforce->fa13, sigma13, eps13);

	gdforce->vc22 = PotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fphi22, sigma22, eps22);
	gdforce->dvc22 = DPotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fa22, sigma22, eps22);

	gdforce->vc23 = PotentialFunction(gdforce->potType, Rcut23, 
					gdforce->fphi23, sigma23, eps23);
	gdforce->dvc23 = DPotentialFunction(gdforce->potType, Rcut23, 
					gdforce->fa23, sigma23, eps23);

	gdforce->vc33 = PotentialFunction(gdforce->potType, Rcut33, 
					gdforce->fphi33, sigma33, eps33);
	gdforce->dvc33 = DPotentialFunction(gdforce->potType, Rcut33, 
					gdforce->fa33, sigma33, eps33);
}
// TERMINAN RUTINAS PARA EL POTENCIAL LENNARD-JONES ----------------------------

// COMIENZAN RUTINAS PARA EL POTENCIAL LOGARITMO ---------------------------
local void ForcePotentialCalc_LN(real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real vc, real dvc)
{
	real rri;

	rri=dr2/ssq;
	*fPot = fphi*rlog(rri+1.0);
	*fAcc = -fa/(rri+1.0);
}

local void PotentialCalc_LN(real *fPot, 
	real dr2, real fphi, real ssq, real vc, real dvc)
{
	real rri;

	rri=dr2/ssq;
	*fPot = fphi*rlog(rri+1.0);
}

local real PotentialFunction_LN(real r, real fphi, real sigma, real eps)
{
	real s2, r2, rri, fPot;

	s2 = sigma*sigma;
	r2 = r*r;
	rri=r2/s2;;
	fPot = fphi*rlog(rri+1.0);
	return (fPot);
}

local real DPotentialFunction_LN(real r, real fa, real sigma, real eps)
{
	real s2, r2, rri, rri3, dfPot;

	s2 = sigma*sigma;
	r2 = r*r;
	rri=r2/s2;
	dfPot = fa*r/(rri+1.0);
	return (dfPot);
}

local void PotentialParameters_LN(real sigma11, real sigma12, real sigma13, 
	real sigma22, real sigma23, real sigma33,
	real eps11, real eps12, real eps13, real eps22, real eps23, real eps33, 
	real Rcut11, real Rcut12, real Rcut13, real Rcut22, real Rcut23, real Rcut33, 
	global_data_tree_tljforcecalc *gdforce)
{

	gdforce->fphi11 = eps11/LNE;
	gdforce->ssq11 = sigma11*sigma11;
	gdforce->fa11 = 2.0*eps11/gdforce->ssq11;
	gdforce->fphi12 = eps12/LNE;
	gdforce->ssq12 = sigma12*sigma12;
	gdforce->fa12 = 2.0*eps12/gdforce->ssq12;

	gdforce->fphi13 = eps13/LNE;
	gdforce->ssq13 = sigma13*sigma13;
	gdforce->fa13 = 2.0*eps13/gdforce->ssq13;

	gdforce->fphi22 = eps22/LNE;
	gdforce->ssq22 = sigma22*sigma22;
	gdforce->fa22 = 2.0*eps22/gdforce->ssq22;

	gdforce->fphi23 = eps23/LNE;
	gdforce->ssq23 = sigma23*sigma23;
	gdforce->fa23 = 2.0*eps23/gdforce->ssq23;

	gdforce->fphi33 = eps33/LNE;
	gdforce->ssq33 = sigma33*sigma33;
	gdforce->fa33 = 2.0*eps33/gdforce->ssq33;

	gdforce->vc11 = PotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fphi11, sigma11, eps11);
	gdforce->dvc11 = DPotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fa11, sigma11, eps11);
	gdforce->vc12 = PotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fphi12, sigma12, eps12);
	gdforce->dvc12 = DPotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fa12, sigma12, eps12);

	gdforce->vc13 = PotentialFunction(gdforce->potType, Rcut13, 
					gdforce->fphi13, sigma13, eps13);
	gdforce->dvc13 = DPotentialFunction(gdforce->potType, Rcut13, 
					gdforce->fa13, sigma13, eps13);

	gdforce->vc22 = PotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fphi22, sigma22, eps22);
	gdforce->dvc22 = DPotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fa22, sigma22, eps22);

	gdforce->vc23 = PotentialFunction(gdforce->potType, Rcut23, 
					gdforce->fphi23, sigma23, eps23);
	gdforce->dvc23 = DPotentialFunction(gdforce->potType, Rcut23, 
					gdforce->fa23, sigma23, eps23);

	gdforce->vc33 = PotentialFunction(gdforce->potType, Rcut33, 
					gdforce->fphi33, sigma33, eps33);
	gdforce->dvc33 = DPotentialFunction(gdforce->potType, Rcut33, 
					gdforce->fa33, sigma33, eps33);
}
// TERMINAN RUTINAS PARA EL POTENCIAL LOGARITMO ----------------------------

// COMIENZAN RUTINAS PARA LA FUERZA DE LORENTZ ---------------------------
local void ForcePotentialCalc_Lorentz(bodyptr p, bodyptr q,
	real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real vc, real dvc)
{
	real rri, dr;

	dr2 += ssq;
	dr=rsqrt(dr2);
	rri=Mass(p)*Mass(q)/dr;
	*fPot = -fphi*rri;
	*fAcc = -fa*rri/dr2;
	Acc(p)[0] += fa*Vel(p)[1]/Mass(p);
	Acc(p)[1] += -fa*Vel(p)[0]/Mass(p);
}

local void PotentialCalc_Lorentz(bodyptr p, bodyptr q,
	real *fPot, 
	real dr2, real fphi, real ssq, real vc, real dvc)
{
	real rri, dr;

	dr2 += ssq;
	dr=rsqrt(dr2);
	rri=Mass(q)/dr2;
	*fPot = -fphi*rri;
}

local real PotentialFunction_Lorentz(real r, real fphi, real sigma, real eps)
{
	real s2, r2, dr, rri, fPot;

	s2 = sigma*sigma;
	r2 = r*r+s2;
	dr=rsqrt(r2);
	rri=1.0/dr;
	fPot = -fphi*rri;
	return (fPot);
}

local real DPotentialFunction_Lorentz(real r, real fa, real sigma, real eps)
{
	real s2, r2, rri, dr, dfPot;

	s2 = sigma*sigma;
	r2 = r*r+s2;
	dr=rsqrt(r2);
	rri=1.0/r2;
	dfPot = fa*rri;
	return (dfPot);
}

local void PotentialParameters_Lorentz(real sigma11, real sigma12, real sigma13, 
	real sigma22, real sigma23, real sigma33,
	real eps11, real eps12, real eps13, real eps22, real eps23, real eps33, 
	real Rcut11, real Rcut12, real Rcut13, real Rcut22, real Rcut23, real Rcut33, 
	global_data_tree_tljforcecalc *gdforce)
{

	gdforce->fphi11 = eps11;
	gdforce->ssq11 = sigma11*sigma11;
	gdforce->fa11 = eps11;

	gdforce->fphi12 = eps12;
	gdforce->ssq12 = sigma12*sigma12;
	gdforce->fa12 = eps12;

	gdforce->fphi13 = eps13;
	gdforce->ssq13 = sigma13*sigma13;
	gdforce->fa13 = eps13;

	gdforce->fphi22 = eps22;
	gdforce->ssq22 = sigma22*sigma22;
	gdforce->fa22 = eps22;

	gdforce->fphi23 = eps23;
	gdforce->ssq23 = sigma23*sigma23;
	gdforce->fa23 = eps23;

	gdforce->fphi33 = eps33;
	gdforce->ssq33 = sigma33*sigma33;
	gdforce->fa33 = eps33;

	gdforce->vc11 = PotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fphi11, sigma11, eps11);
	gdforce->dvc11 = DPotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fa11, sigma11, eps11);
	gdforce->vc12 = PotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fphi12, sigma12, eps12);
	gdforce->dvc12 = DPotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fa12, sigma12, eps12);

	gdforce->vc13 = PotentialFunction(gdforce->potType, Rcut13, 
					gdforce->fphi13, sigma13, eps13);
	gdforce->dvc13 = DPotentialFunction(gdforce->potType, Rcut13, 
					gdforce->fa13, sigma13, eps13);

	gdforce->vc22 = PotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fphi22, sigma22, eps22);
	gdforce->dvc22 = DPotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fa22, sigma22, eps22);

	gdforce->vc23 = PotentialFunction(gdforce->potType, Rcut23, 
					gdforce->fphi23, sigma23, eps23);
	gdforce->dvc23 = DPotentialFunction(gdforce->potType, Rcut23, 
					gdforce->fa23, sigma23, eps23);

	gdforce->vc33 = PotentialFunction(gdforce->potType, Rcut33, 
					gdforce->fphi33, sigma33, eps33);
	gdforce->dvc33 = DPotentialFunction(gdforce->potType, Rcut33, 
					gdforce->fa33, sigma33, eps33);
}
// TERMINAN RUTINAS PARA LA FUERZA DE LORENTZ ----------------------------

// TERMINAN RUTINAS PARA EL CALCULO DE LA FUERZA Y EL POTENCIAL ================

local void TransportCalc(bodyptr p, vector dr, real fPot, real fAcc)
{
	int k;
	matrix w;

	en(p) += fPot;
	OUTVP(w, dr, dr);
	MULMS(w, w, fAcc);
	DO_COORD(k) {
		VVAdd(rf(p)[k], w[k]);
	}
}

local void TransportCalc2(bodyptr p, bodyptr q, vector dr, real fPot, real fAcc)
{
	int k;
	matrix w;

	en(p) += fPot;
	en(q) += fPot;
	OUTVP(w, dr, dr);
	MULMS(w, w, fAcc);
	DO_COORD(k) {
		VVAdd(rf(p)[k], w[k]);
		VVAdd(rf(q)[k], w[k]);
	}
}

