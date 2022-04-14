/* ==============================================================================
	MODULE: nagbody_tree_bljforcecalc.c		[NagBody]
	Written by: Mario A. Rodriguez-Meza
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

#include "globaldefs.h"



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

// COMIENZA BLOQUE METODO DIRECTO E INEFICIENTE PARA CALCULAR LA FUERZA

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

//

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

