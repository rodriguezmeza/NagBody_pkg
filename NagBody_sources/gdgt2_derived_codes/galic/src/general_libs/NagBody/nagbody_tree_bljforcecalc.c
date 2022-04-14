/* ==============================================================================
	MODULE: nagbody_tree_bljforcecalc.c		[NagBody]
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
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
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

// Bloque para normal
local void normal_walktree(bodyptr, nodeptr, real, vector, 
						   real *, vector, global_data_tree_bljforcecalc *);
local void sumnode(bodyptr, cellptr, cellptr, vector, real *, vector, 
					global_data_tree_bljforcecalc *);

// Bloque para normal2
local void normal_walktree2(bodyptr, nodeptr, real, 
							global_data_tree_bljforcecalc *);
local void sumnode2(bodyptr, cellptr, cellptr, global_data_tree_bljforcecalc *);

// Bloque para normal3
local void normal_walktree3(bodyptr, nodeptr, real, 
							global_data_tree_bljforcecalc *);
local void sumnode3(bodyptr, cellptr, cellptr, global_data_tree_bljforcecalc *);

// Bloque para barnes
local void walktree_barnes(nodeptr *, nodeptr *, cellptr, cellptr,
                    nodeptr, real, vector, global_data_tree_bljforcecalc *);
local bool accept_blj(nodeptr, nodeptr,	global_data_tree_bljforcecalc *);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector, global_data_tree_bljforcecalc *);
local void gravsum(bodyptr, cellptr, cellptr, global_data_tree_bljforcecalc *);
local void sumnode_barnes(bodyptr, cellptr, cellptr, vector, real *, vector, 
							global_data_tree_bljforcecalc *);

// Bloque para barnes2
local void walktree_barnes2(nodeptr *, nodeptr *, cellptr *, bodyptr *,
                    nodeptr, real, vector, global_data_tree_bljforcecalc *);
local void walksub2(nodeptr *, nodeptr *, cellptr *, bodyptr *,
                   nodeptr, real, vector, global_data_tree_bljforcecalc *);
local void sumnode_barnes2(bodyptr, bodyptr *, bodyptr *, 
							global_data_tree_bljforcecalc *);

// Bloque para nblist
local void normal_walktree_nblist(bodyptr, nodeptr, real, vector, 
						   real *, vector, global_data_tree_bljforcecalc *);
local void sumnode_nblist_01(bodyptr, cellptr, cellptr, vector, real *, vector, 
								global_data_tree_bljforcecalc *);
local void sumnode_nblist_02(bodyptr, vector, real *, vector, 
								global_data_tree_bljforcecalc *);

local void bljFactors(bodyptr, bodyptr, real *, real *, real *,
		real *, real *, real *, int *, global_data_tree_bljforcecalc *);

// Bloque que define la fuerza y el potencial de interaccion ===================

local void ForcePotentialCalc(int, real *, real *, real, 
								real, real, real, real, real, real, int);
local void PotentialCalc(int, real *, real, real, real, real, real, real, int);
local real PotentialFunction(int, real, real , real, real, int);
local real DPotentialFunction(int, real, real , real, real, int);

// Bloque Lennard-Jones ....
local void ForcePotentialCalc_LJ(real *, real *, real, 
								real, real, real, real, real, real, int);
local void PotentialCalc_LJ(real *, real, real, real, real, real, real, int);
local real PotentialFunction_LJ(real, real , real, real, int);
local real DPotentialFunction_LJ(real, real , real, real, int);
local void PotentialParameters_LJ(real, real, real,
	real, real, real, real, real, real, global_data_tree_bljforcecalc *);

// Bloque Shifted-Lennard-Jones ....
local void ForcePotentialCalc_SLJ(real *, real *, real, 
								real, real, real, real, real, real, int);
local void PotentialCalc_SLJ(real *, real, real, real, real, real, real, int);
local void PotentialParameters_SLJ(real, real, real,
	real, real, real, real, real, real, global_data_tree_bljforcecalc *);

// Bloque Pot table in file ....
local void ForcePotentialCalc_File(real *, real *, real, 
								real, real, real, real, real, real, int);
local void PotentialCalc_File(real *, real, real, real, real, real, real, int);
local void PotentialParameters_File(real, real, real,
	real, real, real, real, real, real, global_data_tree_bljforcecalc *);

// Fin del bloque de fuerzas y potenciales =====================================

local void TransportCalc(bodyptr, vector, real, real);
local void TransportCalc2(bodyptr, bodyptr, vector, real, real);


// COMIENZA METODO NORMAL DE CALCULO DE LA FUERZA

void ind_ljforcecalc_normal(bodyptr btab, int nbody, 
	global_data_tree_bljforcecalc *gdforce, global_data_tree *gdtree, 
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

void ljforcecalc_normal(bodyptr btab, int nbody, 
	global_data_tree_bljforcecalc *gdforce, global_data_tree *gdtree)
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
		DIVVS(Acc(p), acc0, Mass(p));
		Phi(p) += phi0; gdforce->uSum += Phi(p);
	}
	gdforce->virSum = 0.5*gdforce->virSum;
	gdforce->uSum = 0.5*gdforce->uSum;

    gdforce->cpuforce = cputime() - cpustart;
}

local void normal_walktree(bodyptr p, nodeptr q, real qsize, vector pos0, 
			real *phi0, vector acc0, global_data_tree_bljforcecalc *gdforce)
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
				for (k = 0; k < NDIM; k++)
					dr[k]=dr[k] - 
						((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
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
				   global_data_tree_bljforcecalc *gdforce)
{
    cellptr q;
    real dr2, fPot, fAcc;
    vector dr;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	int ij;

    for (q = start; q < finish; q++) {
		bljFactors(p,(bodyptr)q,&RcutSq,&ssq,&fphi,&fa,&vc,&dvc,&ij,gdforce);

		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, RcutSq, vc, dvc, ij);
			*phi0 += fPot;
			ADDMULVS(acc0, dr, fAcc);
			gdforce->virSum += fAcc*dr2;
			gdforce->nbbcalc += 1;

			if (gdforce->computeTransport)
				TransportCalc(p, dr, fPot, fAcc);
		}
    }
}

// TERMINA METODO NORMAL DE CALCULO DE LA FUERZA


// COMIENZA METODO NORMAL2 DE CALCULO DE LA FUERZA

void ljforcecalc_normal2(bodyptr btab, int nbody, 
			global_data_tree_bljforcecalc *gdforce, global_data_tree *gdtree)
{
    bodyptr p;
    double cpustart;

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
							global_data_tree_bljforcecalc *gdforce)
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
				VWrapAll_ptr(dr);
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
					global_data_tree_bljforcecalc *gdforce)
{
    cellptr q;
    real dr2, fPot, fAcc;
    vector dr;
	bodyptr qb;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	int ij;

	DO_BODY(q, start, finish) {
		qb = ((bodyptr) q);
		bljFactors(p, qb, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, &ij, gdforce);

		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, RcutSq, vc, dvc, ij);
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

void ljforcecalc_normal3(bodyptr btab, int nbody, 
			global_data_tree_bljforcecalc *gdforce, global_data_tree *gdtree)
{
    bodyptr p;
    double cpustart;
 
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
							global_data_tree_bljforcecalc *gdforce)
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
				VWrapAll_ptr(dr);
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
					global_data_tree_bljforcecalc *gdforce)
{
    cellptr q;
    real dr2, fPot, fAcc;
    vector dr;
	real rri, rri3;
	bodyptr qb;
	real RcutSq, ssq, fphi, fa, vc, dvc;

	DO_BODY(q, start, finish) {
		qb = ((bodyptr) q);
		if ( Type(p) == Type(qb) ) {
			if (Type(p)==BODY1) {
				RcutSq=gdforce->RcutSq11; 
				ssq=gdforce->ssq11; fphi=gdforce->fphi11; fa=gdforce->fa11;
				vc = gdforce->vc11; dvc = gdforce->dvc11;
				DOTPSUBV(dr2, dr, Pos(p), Pos(qb));
				VWrapAll_ptr(dr);
				DOTVP(dr2, dr, dr);
				if (dr2<RcutSq) {
//					ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, 
//						fphi, fa, ssq, vc, dvc);

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
				VWrapAll_ptr(dr);
				DOTVP(dr2, dr, dr);
				if (dr2<RcutSq) {
//					ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, 
//						fphi, fa, ssq, vc, dvc);

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
			VWrapAll_ptr(dr);
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


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE LISTA DE VECINOS

void ljforcecalc_nblist(bodyptr btab, int nbody, 
			global_data_tree_bljforcecalc *gdforce, global_data_tree *gdtree)
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
							global_data_tree_bljforcecalc *gdforce)
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
				for (k = 0; k < NDIM; k++)
					dr[k]=dr[k] -
						((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
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
					global_data_tree_bljforcecalc *gdforce)
{
    cellptr p;
    real dr2, phi_p, mr3i;
    vector dr;
	real rri, rri3;
	int k, ip;
	real RcutSq;
 
    for (p = start; p < finish; p++) {          
		DOTPSUBV(dr2, dr, pos0, Pos(p));
		for (k = 0; k < NDIM; k++)
			dr[k]=dr[k] - ((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
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
							global_data_tree_bljforcecalc *gdforce)
{
    bodyptr p;
    real dr2, fPot, fAcc;
    vector dr;
	int ip;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	int ij;

    for (ip = 0; ip < nblist; ip++) {
		p = bodytab + activenb[ip];

		bljFactors(p0, p, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, &ij, gdforce);

		DOTPSUBV(dr2, dr, pos0, Pos(p));
		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
			fa, ssq, RcutSq, vc, dvc, ij);
		*phi0 += fPot;
		ADDMULVS(acc0, dr, fAcc);
		gdforce->virSum += fAcc*dr2;
		if (gdforce->computeTransport)
			TransportCalc(p0, dr, fPot, fAcc);
    }
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE LISTA DE VECINOS


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CALCULO DIRECTO DE LA FUERZA

void ind_ljforcecalc_direct(bodyptr btab, int nbody, 
						global_data_tree_bljforcecalc *gdforce, 
						bodyptr p, real *uSum, real *virSum, short *flag)
{
    bodyptr q;
    double cpustart;
    vector dr;
	real rri, rri3, fPot, fAcc, dr2;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	real dr2min=0.5;
	int ij;
 
    cpustart = cputime();

	Phi(p) = 0.0;
	CLRV(Acc(p));

    gdforce->nbbcalc = 0;   
	gdforce->uSum = 0.0; gdforce->virSum = 0.0;
	*flag = 1;

	DO_BODY(q, btab, btab+nbody) {
		if (p==q) continue;

		bljFactors(p, q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, &ij, gdforce);

		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		if (dr2 < dr2min) {
			*flag = 0;
			break;
		}

		if (dr2 < RcutSq) {
			ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, RcutSq, vc, dvc, ij);
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

void ljforcecalc_direct(bodyptr btab, int nbody, 
						global_data_tree_bljforcecalc *gdforce)
{
    bodyptr p, q;
    double cpustart;
    vector dr;
	real fPot, fAcc, dr2;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	int ij;
 
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

			bljFactors(p, q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, &ij, gdforce);

			DOTPSUBV(dr2, dr, Pos(p), Pos(q));
			VWrapAll_ptr(dr);
			DOTVP(dr2, dr, dr);

			if (dr2 < RcutSq) {
				ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
					fa, ssq, RcutSq, vc, dvc, ij);
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

void ljforcecalc_direct2(bodyptr btab, int nbody, 
						global_data_tree_bljforcecalc *gdforce)
{
    bodyptr p, q;
    double cpustart;
    vector dr;
	real fPot, fAcc, dr2;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	int ij;

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
			bljFactors(p, q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, &ij, gdforce);

			DOTPSUBV(dr2, dr, Pos(p), Pos(q));
			VWrapAll_ptr(dr);
			DOTVP(dr2, dr, dr);

			if (dr2<RcutSq) {
				ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
					fa, ssq, RcutSq, vc, dvc, ij);
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

void ljforcecalc_barnes(bodyptr btab, int nbody, 
			global_data_tree_bljforcecalc *gdforce, global_data_tree *gdtree)
{
    double cpustart;
    vector rmid;
	bodyptr p;
 
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
					global_data_tree_bljforcecalc *gdforce)
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
						global_data_tree_bljforcecalc *gdforce)
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
	DO_COORD(k)
		dr[k]=dr[k] - ( (real) (nint(dr[k]/gdforce->Box[k])) )*gdforce->Box[k];
	DOTVP(drpq2, dr, dr);
	drpq = rsqrt(drpq2);

	if ( drpq >= Rcut(p)+Rcut(q) )
		return (TRUE);
	else
		return (FALSE);
}

local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid, 
				   global_data_tree_bljforcecalc *gdforce)
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
					global_data_tree_bljforcecalc *gdforce)
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
				   global_data_tree_bljforcecalc *gdforce)
{
    cellptr p;
    real dr2, fPot, fAcc;
    vector dr;
	real RcutSq, ssq, fphi,fa, vc, dvc;
	int ij;

    for (p = start; p < finish; p++) {
		bljFactors(p0, (bodyptr)p, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, &ij, gdforce);

		DOTPSUBV(dr2, dr, pos0, Pos(p));
		VWrapAll_ptr(dr);
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
				fa, ssq, RcutSq, vc, dvc, ij);
			*phi0 += fPot;
			ADDMULVS(acc0, dr, fAcc);

			gdforce->virSum += fAcc*dr2;
			gdforce->nbbcalc += 1;
			if (gdforce->computeTransport)
				TransportCalc(p0, dr, fPot, fAcc);
		}
    }
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES PARA EL CALCULO DE LA FUERZA


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS PARA EL CALCULO DE LA FUERZA

void ljforcecalc_cellsmethod(bodyptr btab, int nbody, 
							global_data_tree_bljforcecalc *gdforce)
{
	bodyptr p, p1, p2;
	vector dr;
	vector invWid, rs, shift;
	vectorI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	real fPot, fAcc, rr;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
	double cpustart;
	real RcutSq, ssq, fphi, fa, vc, dvc;
	int ij;

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
					VCellWrapAll_ptr ();
					m2 = VLinear (m2v, gdforce->cells) + nbody;
					DO_CELL_ptr (j1, m1) {
						DO_CELL_ptr (j2, m2) {
							if (m1 != m2 || j2 < j1) {
								p1 = btab+j1; p2 = btab+j2;
								bljFactors(p1,p2,&RcutSq,&ssq,&fphi,&fa,&vc,
										&dvc,&ij,gdforce);
								SUBV(dr, Pos(p1), Pos(p2));
								VVSub (dr, shift);
								rr = VLenSq (dr);
								if (rr < RcutSq) {
									ForcePotentialCalc(gdforce->potType, &fPot, 
										&fAcc, rr, fphi, fa, ssq, RcutSq, vc, dvc, ij);
									VVSAdd(Acc(p1), fAcc/Mass(p1), dr);
									VVSAdd (Acc(p2), - fAcc/Mass(p2), dr);
									Phi(p1) += fPot; Phi(p2) +=fPot;
									gdforce->uSum += fPot;
									gdforce->virSum += fAcc * rr;
									gdforce->nbbcalc+=2;
									if (gdforce->computeTransport) {
										TransportCalc2(p1, p2, dr, fPot, fAcc);
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

void ljforcecalc_cellsmethod3(bodyptr btab, int nbody, 
								global_data_tree_bljforcecalc *gdforce)
{
	bodyptr p, p1, p2;
	vector dr, invWid, rs, shift;
	vectorI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	real fcVal, rr, rri, rri3, uVal;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
	double cpustart;
	real RcutSq, ssq,fphi,fa;

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
					VCellWrapAll_ptr ();
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
//											ForcePotentialCalc(gdforce->potType,
//												&uVal,&fcVal,rr,fphi,
//												fa,ssq,vc,dvc);
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
//												TransportCalc2(p1,p2,dr,
//															uVal,fcVal);
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
//											ForcePotentialCalc(gdforce->potType,
//												&uVal,&fcVal,rr,fphi,
//												fa,ssq,vc,dvc);
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
//												TransportCalc2(p1,p2,dr,
//																uVal,fcVal);
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
//										ForcePotentialCalc(gdforce->potType,
//											&uVal,&fcVal,rr,fphi,fa,ssq,vc,dvc);
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
//											TransportCalc2(p1,p2,dr,uVal,fcVal);
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

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS3 PARA EL CALCULO DE LA FUERZA


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS PARA EL CALCULO DEL POTENCIAL

void ljpotcalc_cellsmethod(bodyptr btab, int nbody, 
							global_data_tree_bljforcecalc *gdforce)
{
	bodyptr p, p1, p2;
	vector dr, invWid, rs, shift;
	vectorI cc, m1v, m2v, vOff[] = OFFSET_VALS;
	real rr, fPot;
	int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
	double cpustart;
	real RcutSq, ssq, fphi,	fa, vc, dvc;
	int ij;

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
					VCellWrapAll_ptr ();
					m2 = VLinear (m2v, gdforce->cells) + nbody;
					DO_CELL_ptr (j1, m1) {
						DO_CELL_ptr (j2, m2) {
							if (m1 != m2 || j2 < j1) {
								p1 = btab+j1; p2 = btab+j2;
								bljFactors(p1,p2,&RcutSq,&ssq,&fphi,
											&fa,&vc,&dvc,&ij,gdforce);
								SUBV(dr, Pos(p1), Pos(p2));
								VVSub (dr, shift);
								rr = VLenSq (dr);
								if (rr < RcutSq) {
									PotentialCalc(gdforce->potType, &fPot, 
										rr, fphi, ssq, RcutSq, vc, dvc, ij);
									Phi(p1) += fPot; Phi(p2) +=fPot;
									gdforce->uSum += fPot;
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

void ljforcecalc_barnes2(bodyptr btab, int nbody, 
			global_data_tree_bljforcecalc *gdforce, global_data_tree *gdtree)
{
    double cpustart;
    vector rmid;
	bodyptr p;

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
}

#undef FACTIVE

local void walktree_barnes2(nodeptr *aptr, nodeptr *nptr, cellptr *cptr, 
					bodyptr *bptr, nodeptr p, real psize, vector pmid, 
					global_data_tree_bljforcecalc *gdforce)
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
				   global_data_tree_bljforcecalc *gdforce)
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
				   global_data_tree_bljforcecalc *gdforce)
{
    cellptr c;
    real dr2, fPot, fAcc;
    vector dr;
	int k;
	real RcutSq, ssq, fphi,fa, vc, dvc;
	int ij;
	bodyptr *q;

	DO_BODY(q, start, finish) {
		if ( Id(p) < Id(*q) ) {
			bljFactors(p, *q, &RcutSq, &ssq, &fphi, &fa, &vc, &dvc, &ij, gdforce);

			DOTPSUBV(dr2, dr, Pos(p), Pos(*q));
			VWrapAll_ptr(dr);
			DOTVP(dr2, dr, dr);

			if (dr2<RcutSq) {
				ForcePotentialCalc(gdforce->potType, &fPot, &fAcc, dr2, fphi, 
					fa, ssq, RcutSq, vc, dvc, ij);
				ADDMULVS(Acc(p), dr, fAcc);
				ADDMULVS(Acc(*q), dr, -fAcc);
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
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES2 PARA EL CALCULO DE LA FUERZA


// COMIENZA RUTINAS AUXILIARES PARA EL CALCULO DE LA FUERZA ...

local void bljFactors(bodyptr p, bodyptr q, real *RcutSq, real *ssq, real *fphi,
		real *fa, real *vc, real *dvc, int *ij, global_data_tree_bljforcecalc *gdforce)
{
	if ( Type(p) == Type(q) ) {
		if (Type(p)==BODY1) {
			*RcutSq=gdforce->RcutSq11; 
			*ssq=gdforce->ssq11; 
			*fphi=gdforce->fphi11; 
			*fa=gdforce->fa11;
			*vc=gdforce->vc11;
			*dvc=gdforce->dvc11;
			*ij=11;
		} else {
			*RcutSq=gdforce->RcutSq22; 
			*ssq=gdforce->ssq22; 
			*fphi=gdforce->fphi22; 
			*fa=gdforce->fa22;
			*vc=gdforce->vc22;
			*dvc=gdforce->dvc22;
			*ij=22;
		}
	} else {
		*RcutSq=gdforce->RcutSq12; 
		*ssq=gdforce->ssq12; 
		*fphi=gdforce->fphi12; 
		*fa=gdforce->fa12;
		*vc=gdforce->vc12;
		*dvc=gdforce->dvc12;
		*ij=12;
	}
}

// COMIENZA RUTINAS PARA EL CALCULO DE LA FUERZA Y EL POTENCIAL ================

local void ForcePotentialCalc(int type, real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real rc, real vc, real dvc, int ij)
{
    switch(type) {
        case POTTYPE_LJ:
			ForcePotentialCalc_LJ(fPot, fAcc, dr2, fphi, fa, ssq, rc, vc, dvc, ij);
			break;
        case POTTYPE_SLJ:
			ForcePotentialCalc_SLJ(fPot, fAcc, dr2, fphi, fa, ssq, rc, vc, dvc, ij);
			break;
        case POTTYPE_FILE:
			ForcePotentialCalc_File(fPot, fAcc, dr2, fphi, fa, ssq, rc, vc, dvc, ij);
			break;
        default:
			ForcePotentialCalc_LJ(fPot, fAcc, dr2, fphi, fa, ssq, rc, vc, dvc, ij);
			break;
    }
}

local void PotentialCalc(int type, real *fPot, 
	real dr2, real fphi, real ssq, real rc, real vc, real dvc, int ij)
{
    switch(type) {
        case POTTYPE_LJ:
			PotentialCalc_LJ(fPot, dr2, fphi, ssq, rc, vc, dvc, ij);
			break;
        case POTTYPE_SLJ:
			PotentialCalc_SLJ(fPot, dr2, fphi, ssq, rc, vc, dvc, ij);
			break;
        case POTTYPE_FILE:
			PotentialCalc_File(fPot, dr2, fphi, ssq, rc, vc, dvc, ij);
			break;
        default:
			PotentialCalc_LJ(fPot, dr2, fphi, ssq, rc, vc, dvc, ij);
			break;
    }
}

local real PotentialFunction(int type, real r, real fphi, real sigma, real eps, int ij)
{
	real fPot;

    switch(type) {
        case POTTYPE_LJ:
			fPot=PotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
        case POTTYPE_SLJ:
			fPot=PotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
        case POTTYPE_FILE:
			fPot=PotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
        default:
			fPot=PotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
    }
	return (fPot);
}

local real DPotentialFunction(int type, real r, real fphi, real sigma, real eps, int ij)
{
	real dfPot;

    switch(type) {
        case POTTYPE_LJ:
			dfPot=DPotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
        case POTTYPE_SLJ:
			dfPot=DPotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
        case POTTYPE_FILE:
			dfPot=DPotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
        default:
			dfPot=DPotentialFunction_LJ(r, fphi, sigma, eps, ij);
			break;
    }
	return (dfPot);
}

void PotentialParameters(int type, real sigma11, real sigma12, real sigma22,
	real eps11, real eps12, real eps22, 
	real Rcut11, real Rcut12, real Rcut22,
	global_data_tree_bljforcecalc *gdforce)
{
    switch(type) {
        case POTTYPE_LJ:
			PotentialParameters_LJ(sigma11, sigma12, sigma22,
				eps11, eps12, eps22, Rcut11, Rcut12, Rcut22, gdforce);
			break;
        case POTTYPE_SLJ:
			PotentialParameters_SLJ(sigma11, sigma12, sigma22,
				eps11, eps12, eps22, Rcut11, Rcut12, Rcut22, gdforce);
			break;
        case POTTYPE_FILE:
			PotentialParameters_File(sigma11, sigma12, sigma22,
				eps11, eps12, eps22, Rcut11, Rcut12, Rcut22, gdforce);
			break;
        default:
			PotentialParameters_LJ(sigma11, sigma12, sigma22,
				eps11, eps12, eps22, Rcut11, Rcut12, Rcut22, gdforce);
			break;
    }
}

// COMIENZAN RUTINAS PARA EL POTENCIAL LENNARD-JONES ---------------------------
local void ForcePotentialCalc_LJ(real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real rc2, real vc, real dvc, int ij)
{
	real rri, rri3;

	rri=ssq/dr2; rri3=rri*rri*rri;
	*fPot = fphi*(rri3-1.0)*rri3;
	*fAcc = fa*rri3*(rri3-0.5)*rri;
}

local void PotentialCalc_LJ(real *fPot, 
	real dr2, real fphi, real ssq, real rc2, real vc, real dvc, int ij)
{
	real rri, rri3;

	rri=ssq/dr2; rri3=rri*rri*rri;
	*fPot = fphi*(rri3-1.0)*rri3;
}

local real PotentialFunction_LJ(real r, real fphi, real sigma, real eps, int ij)
{
	real s2, r2, rri, rri3, fPot;

	s2 = sigma*sigma;
	r2 = r*r;
	rri=s2/r2; rri3=rri*rri*rri;
	fPot = fphi*(rri3-1.0)*rri3;
	return (fPot);
}

local real DPotentialFunction_LJ(real r, real fa, real sigma, real eps, int ij)
{
	real s2, r2, rri, rri3, dfPot;

	s2 = sigma*sigma;
	r2 = r*r;
	rri=s2/r2; rri3=rri*rri*rri;
	dfPot = -fa*(rri3-0.5)*rri3*rri*r;
	return (dfPot);
}

local void PotentialParameters_LJ(real sigma11, real sigma12, real sigma22,
	real eps11, real eps12, real eps22, real Rcut11, real Rcut12, real Rcut22,
	global_data_tree_bljforcecalc *gdforce)
{

	gdforce->fphi11 = 4.0*eps11;
	gdforce->ssq11 = sigma11*sigma11;
	gdforce->fa11 = 48.0*eps11/gdforce->ssq11;
	gdforce->fphi12 = 4.0*eps12;
	gdforce->ssq12 = sigma12*sigma12;
	gdforce->fa12 = 48.0*eps12/gdforce->ssq12;
	gdforce->fphi22 = 4.0*eps22;
	gdforce->ssq22 = sigma22*sigma22;
	gdforce->fa22 = 48.0*eps22/gdforce->ssq22;

	gdforce->vc11 = PotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fphi11, sigma11, eps11, 11);
	gdforce->dvc11 = DPotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fa11, sigma11, eps11, 11);
	gdforce->vc12 = PotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fphi12, sigma12, eps12, 12);
	gdforce->dvc12 = DPotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fa12, sigma12, eps12, 12);
	gdforce->vc22 = PotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fphi22, sigma22, eps22, 22);
	gdforce->dvc22 = DPotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fa22, sigma22, eps22, 22);
}
// TERMINAN RUTINAS PARA EL POTENCIAL LENNARD-JONES ----------------------------

// COMIENZAN RUTINAS PARA EL POTENCIAL SHIFTED-LENNARD-JONES -------------------
local void ForcePotentialCalc_SLJ(real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real rc2, real vc, real dvc, int ij)
{
	real rri, rri3, r, rc;

	rc = rsqrt(rc2);
	r = rsqrt(dr2);
	rri=ssq/dr2; rri3=rri*rri*rri;
	*fPot = fphi*(rri3-1.0)*rri3-vc-dvc*(r-rc);
	*fAcc = fa*rri3*(rri3-0.5)*rri-(rc/r)*dvc;
}

local void PotentialCalc_SLJ(real *fPot, 
	real dr2, real fphi, real ssq, real rc2, real vc, real dvc, int ij)
{
	real rri, rri3, r, rc;

	rc = rsqrt(rc2);
	r = rsqrt(dr2);
	rri=ssq/dr2; rri3=rri*rri*rri;
	*fPot = fphi*(rri3-1.0)*rri3-vc-dvc*(r-rc);
}

local void PotentialParameters_SLJ(real sigma11, real sigma12, real sigma22,
	real eps11, real eps12, real eps22, real Rcut11, real Rcut12, real Rcut22,
	global_data_tree_bljforcecalc *gdforce)
{

	gdforce->fphi11 = 4.0*eps11;
	gdforce->ssq11 = sigma11*sigma11;
	gdforce->fa11 = 48.0*eps11/gdforce->ssq11;
	gdforce->fphi12 = 4.0*eps12;
	gdforce->ssq12 = sigma12*sigma12;
	gdforce->fa12 = 48.0*eps12/gdforce->ssq12;
	gdforce->fphi22 = 4.0*eps22;
	gdforce->ssq22 = sigma22*sigma22;
	gdforce->fa22 = 48.0*eps22/gdforce->ssq22;

	gdforce->vc11 = PotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fphi11, sigma11, eps11, 11);
	gdforce->dvc11 = DPotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fa11, sigma11, eps11, 11);
	gdforce->vc12 = PotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fphi12, sigma12, eps12, 12);
	gdforce->dvc12 = DPotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fa12, sigma12, eps12, 12);
	gdforce->vc22 = PotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fphi22, sigma22, eps22, 22);
	gdforce->dvc22 = DPotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fa22, sigma22, eps22, 22);
}
// TERMINAN RUTINAS PARA EL POTENCIAL SHIFTED-LENNARD-JONES --------------------

// COMIENZAN RUTINAS PARA EL POTENCIAL IN FILE ---------------------------------
local void ForcePotentialCalc_File(real *fPot, real *fAcc, 
	real dr2, real fphi, real fa, real ssq, real rc2, real vc, real dvc, int ij)
{
	pointForcePotptr p, pf, pi;
	int j;
	real dr, r;

	r = rsqrt(dr2);

	pi = forcepottab;
	if ( r <= rPos(pi) ) {
		if (ij==11) {
			*fPot = Pot11(pi);
			*fAcc = Force11(pi);
		} else
			if (ij==12) {
				*fPot = Pot12(pi);
				*fAcc = Force12(pi);
			} else {
				*fPot = Pot22(pi);
				*fAcc = Force22(pi);
			}
		return;
	}

	pf = forcepottab+nforcepot-1;
	if ( r >= rPos(pf) ) {
		if (ij==11) {
			*fPot = Pot11(pf);
			*fAcc = Force11(pf);
		} else
			if (ij==12) {
				*fPot = Pot12(pf);
				*fAcc = Force12(pf);
			} else {
				*fPot = Pot22(pf);
				*fAcc = Force22(pf);
			}
		return;
	}

	dr = ( rPos(pf) - rPos(pi) ) / ((real)(nforcepot-1));
	j = (int) ( ( r - rPos(pi) ) / dr );
	p = forcepottab + j;
	if (ij==11) {
		*fPot = Pot11(p) + ( Pot11(p+1) - Pot11(p) ) * ( r - rPos(p) ) / dr;
		*fAcc = Force11(p) + ( Force11(p+1) - Force11(p) ) * ( r - rPos(p) ) / dr;
	} else
		if (ij==12) {
			*fPot = (Pot12(p)+(Pot12(p+1)-Pot12(p))*(r-rPos(p))/dr);
			*fAcc = (Force12(p)+(Force12(p+1)-Force12(p))*(r-rPos(p))/dr);
		} else {
			*fPot = (Pot22(p)+(Pot22(p+1)-Pot22(p))*(r-rPos(p))/dr);
			*fAcc = (Force22(p)+(Force22(p+1)-Force22(p))*(r-rPos(p))/dr);
		}
}

local void PotentialCalc_File(real *fPot, 
	real dr2, real fphi, real ssq, real rc2, real vc, real dvc, int ij)
{
	pointForcePotptr p, pf, pi;
	int j;
	real dr;

	pi = forcepottab;
	if ( dr2 <= rPos(pi) ) {
		if (ij==11)
			*fPot = Pot11(pi);
		else
			if (ij==12)
				*fPot = Pot12(pi);
			else
				*fPot = Pot22(pi);
		return;
	}

	pf = forcepottab+nforcepot-1;
	if ( dr2 >= rPos(pf) ) {
		if (ij==11)
			*fPot = Pot11(pf);
		else
			if (ij==12)
				*fPot = Pot12(pf);
			else
				*fPot = Pot22(pf);
		return;
	}

	dr = ( rPos(pf) - rPos(pi) ) / ((real)(nforcepot-1));
	j = (int) ( ( dr2 - rPos(pi) ) / dr );
	p = forcepottab + j;
	if (ij==11)
		*fPot = (Pot11(p)+(Pot11(p+1)-Pot11(p))*(dr2-rPos(p))/dr);
	else
		if (ij==12)
			*fPot = (Pot12(p)+(Pot12(p+1)-Pot12(p))*(dr2-rPos(p))/dr);
		else
			*fPot = (Pot22(p)+(Pot22(p+1)-Pot22(p))*(dr2-rPos(p))/dr);
}

local void PotentialParameters_File(real sigma11, real sigma12, real sigma22,
	real eps11, real eps12, real eps22, real Rcut11, real Rcut12, real Rcut22,
	global_data_tree_bljforcecalc *gdforce)
{
	gdforce->vc11 = PotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fphi11, sigma11, eps11, 11);
	gdforce->dvc11 = DPotentialFunction(gdforce->potType, Rcut11, 
					gdforce->fa11, sigma11, eps11, 11);
	gdforce->vc12 = PotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fphi12, sigma12, eps12, 12);
	gdforce->dvc12 = DPotentialFunction(gdforce->potType, Rcut12, 
					gdforce->fa12, sigma12, eps12, 12);
	gdforce->vc22 = PotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fphi22, sigma22, eps22, 22);
	gdforce->dvc22 = DPotentialFunction(gdforce->potType, Rcut22, 
					gdforce->fa22, sigma22, eps22, 22);
}
// TERMINAN RUTINAS PARA EL POTENCIAL IN FILE ----------------------------------


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

