/* ==============================================================================
!	MODULE: tree_forcecalc.c		[md_lj_tree]
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date:																!
!	Purpose: Lennard-Jones force computation (Several methods)					!
!	Language: C																	!
!	Use: forcecalc();															!
!	Routines and functions:	tree_ljforce, stepsystem							!
!	External modules, routines and headers:										!
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

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "mathutil.h"
#include "global_defs.h"
#include "proto_defs.h"


// Bloque para normal
local void normal_walktree(bodyptr, nodeptr, real, vector, 
						   real *, vector);
local void sumnode(cellptr, cellptr, vector, real *, vector);

// Bloque para normal2
local void normal_walktree2(bodyptr, nodeptr, real);
local void sumnode2(bodyptr, cellptr, cellptr);

// Bloque para barnes
local void walktree_barnes(nodeptr *, nodeptr *, cellptr, cellptr,
                    nodeptr, real, vector);
local bool accept_lj(nodeptr, nodeptr, real);
local bool accept_body(nodeptr, nodeptr);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local void gravsum(bodyptr, cellptr, cellptr);
local void sumnode_barnes(cellptr, cellptr, vector, real *, vector);

// Bloque para barnes2
local void walktree_barnes2(nodeptr *, nodeptr *, cellptr, cellptr,
                    nodeptr, real, vector);
local void walksub2(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local bool accept_lj2(nodeptr, nodeptr, real);
local void sumnode_barnes2(bodyptr, cellptr, cellptr);

// Bloque para nblist
local void normal_walktree_nblist(bodyptr, nodeptr, real, vector, 
						   real *, vector);
local void sumnode_nblist_01(cellptr, cellptr, vector, real *, vector);
local void sumnode_nblist_02(vector, real *, vector);

#if !defined(FACTIVE)
#  define FACTIVE  0.75                         
#endif
 
local int actlen; 
local nodeptr *active;                          
local cellptr interact;                         
local int *activenb;                          
local int nblist;


// COMIENZA METODO NORMAL DE CALCULO DE LA FUERZA

void ljforcecalc_normal(bodyptr btab, int nbody)
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

    nbbcalc = nbccalc = 0;   
	uSum = 0.0; virSum = 0.0;

    for (p = btab; p < btab+nbody; p++) {
		SETV(pos0, Pos(p));
		phi0 = 0.0;
		CLRV(acc0);
#if defined(DIAGNOSTICS)
		printf("\nParticle %d\n",p-btab+1);
#endif
		normal_walktree(p, ((nodeptr) root), rsize, pos0,&phi0,acc0);
		ADDMULVS(Acc(p), acc0, fa);
		Phi(p) += phi0; uSum += Mass(p)*Phi(p);
	}
	virSum = 0.5*virSum*fa;
	uSum = 0.5*uSum;

    cpuforce = cputime() - cpustart;            
}

local void normal_walktree(bodyptr p, nodeptr q, real qsize, 
							vector pos0, real *phi0, vector acc0)
{
    nodeptr l;
	real drpq, drpq2, rcell;
    vector dr;
	int k;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
				for (k = 0; k < NDIM; k++)
					dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
				DOTVP(drpq2, dr, dr);
				drpq = rsqrt(drpq2);
				rcell= qsize * rsqrt((real)(NDIM))/2.0;

                if ( drpq >= Rcut+rcell ) { 
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree(p,l,qsize/2,pos0,phi0,acc0);
					}
				}
			} else {
					sumnode(((cellptr) q),( (cellptr) q+1),
						pos0,phi0,acc0);
			}
		}
	}
}

local void sumnode(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real dr2, phi_p, mr3i;
    vector dr;
	real rri, rri3;
	int k, ip;
 
    for (p = start; p < finish; p++) {          
		DOTPSUBV(dr2, dr, pos0, Pos(p));
		for (k = 0; k < NDIM; k++)
			dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			rri=ssq/dr2; rri3=rri*rri*rri;
			phi_p = fphi*(rri3-1.0)*rri3+eps;
			*phi0 += phi_p;
			mr3i = rri3*(rri3-0.5)*rri;
			ADDMULVS(acc0, dr, mr3i);
			virSum += mr3i*dr2;
			nbbcalc += 1;
			ip = (bodyptr)p-bodytab+1;
#if defined(DIAGNOSTICS)
	printf("\nInteracting particle %d\n",ip);
#endif
		}
    }
}

// TERMINA METODO NORMAL DE CALCULO DE LA FUERZA


// COMIENZA METODO NORMAL2 DE CALCULO DE LA FUERZA

void ljforcecalc_normal2(bodyptr btab, int nbody)
{
    bodyptr p;
    double cpustart;
	vector pos0, acc0;
	real phi0;
	int id0;
 
    cpustart = cputime();                       

	DO_BODY(p, btab, btab+nbody) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    nbbcalc = nbccalc = 0;   
	uSum = 0.0; virSum = 0.0;

	DO_BODY(p, btab, btab+nbody)
		normal_walktree2(p, ((nodeptr) root), rsize);

    cpuforce = cputime() - cpustart;            
}

local void normal_walktree2(bodyptr p, nodeptr q, real qsize)
{
    nodeptr l;
	real drpq, drpq2, rcell;
    vector dr;
	int k;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
				VWrapAll(dr);
				DOTVP(drpq2, dr, dr);
				drpq = rsqrt(drpq2);

                if ( drpq < Rcut+Rcrit2(q) )
					DO_DESCENDENTS(l,q)
						normal_walktree2(p,l,qsize/2);
			} else
				if ( Id(p) < Id(q) )
					sumnode2(p,((cellptr) q),( (cellptr) q+1));
		}
	}
}

local void sumnode2(bodyptr p, cellptr start, cellptr finish)
{
    cellptr q;
    real dr2, uVal, fVal;
    vector dr;
	real rri, rri3;
	int k;
	bodyptr qb;
 
	DO_BODY(q, start, finish) {
		qb = ((bodyptr) q);
		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
		VWrapAll(dr);
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			rri=ssq/dr2; rri3=rri*rri*rri;
			uVal = 4.0*(rri3-1.0)*rri3+1.0;
			fVal = 48.0*rri3*(rri3-0.5)*rri;
			ADDMULVS(Acc(p), dr, fVal);
			ADDMULVS(Acc(qb), dr, -fVal);
			Phi(p) += uVal;
			Phi(qb) += uVal;
			virSum += fVal*dr2;
			uSum += uVal;
			nbbcalc += 2;
		}
    }
}

// TERMINA METODO NORMAL2 DE CALCULO DE LA FUERZA


/* COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE LISTA DE VECINOS */

void ljforcecalc_nblist(bodyptr btab, int nbody)
{
    bodyptr p;
    double cpustart;
	vector pos0, acc0;
	real phi0;
 
    cpustart = cputime();                       

    actlen = FACTIVE * 216 * tdepth;            
    activenb = (int *) allocate(actlen * sizeof(int));

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    nbbcalc = nbccalc = 0;   
	uSum = 0.0; virSum = 0.0;

    for (p = btab; p < btab+nbody; p++) {
		SETV(pos0, Pos(p));
		phi0 = 0.0;
		CLRV(acc0);
		nblist=0;
#if defined(DIAGNOSTICS)
		printf("\nParticle %d\n",p-btab+1);
#endif
		normal_walktree_nblist(p, ((nodeptr) root), rsize, pos0,&phi0,acc0);
		int_piksrt(nblist, activenb);
		sumnode_nblist_02(pos0,&phi0,acc0);
       ADDMULVS(Acc(p), acc0, fa);
       Phi(p) += phi0; uSum += Mass(p)*Phi(p);
	}
	virSum = 0.5*virSum*fa;
	uSum = 0.5*uSum;

    cpuforce = cputime() - cpustart;            
}

local void normal_walktree_nblist(bodyptr p, nodeptr q, real qsize, 
							vector pos0, real *phi0, vector acc0)
{
    nodeptr l;
	real drpq, drpq2, rcell;
    vector dr;
	int k;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
				for (k = 0; k < NDIM; k++)
					dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
				DOTVP(drpq2, dr, dr);
				drpq = rsqrt(drpq2);
				rcell= qsize * rsqrt((real)(NDIM))/2.0;

                if ( drpq >= Rcut+rcell ) { 
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree_nblist(p,l,qsize/2,pos0,phi0,acc0);
					}
				}
			} else {
					sumnode_nblist_01(((cellptr) q),( (cellptr) q+1),
						pos0,phi0,acc0);
			}
		}
	}
}

local void sumnode_nblist_01(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real dr2, phi_p, mr3i;
    vector dr;
	real rri, rri3;
	int k, ip;
 
    for (p = start; p < finish; p++) {          
		DOTPSUBV(dr2, dr, pos0, Pos(p));
		for (k = 0; k < NDIM; k++)
			dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			nbbcalc += 1;
			ip = (bodyptr)p-bodytab;
			activenb[nblist]=ip;
			nblist +=1;
		}
    }
}

local void sumnode_nblist_02(vector pos0, real *phi0, vector acc0)
{
    bodyptr p;
    real dr2, phi_p, mr3i;
    vector dr;
	real rri, rri3;
	int k, ip;
 
    for (ip = 0; ip < nblist; ip++) {
		p = bodytab + activenb[ip];
#if defined(DIAGNOSTICS)
		printf("\nInteracting particle %d\n",p-bodytab+1);
#endif
		DOTPSUBV(dr2, dr, pos0, Pos(p));
		for (k = 0; k < NDIM; k++)
			dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
		DOTVP(dr2, dr, dr);

			rri=ssq/dr2; rri3=rri*rri*rri;
			phi_p = fphi*(rri3-1.0)*rri3+eps;
			*phi0 += phi_p;
			mr3i = rri3*(rri3-0.5)*rri;
			ADDMULVS(acc0, dr, mr3i);
			virSum += mr3i*dr2;
    }
}

/* TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE LISTA DE VECINOS */


/* COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE CALCULO DIRECTO DE LA FUERZA */

void ljforcecalc_direct(bodyptr btab, int nbody)
{
    bodyptr p, q;
    double cpustart;
    vector pos0, acc0;
    real phi0;
    vector dr;
	int k;
	real rri, rri3, phi_q, mr3i, dr2;
 
    cpustart = cputime(); 

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    nbbcalc = 0;   
	uSum = 0.0; virSum = 0.0;
    
    for (p = btab; p < btab+nbody; p++) {
       SETV(pos0, Pos(p));
       phi0 = 0.0;
       CLRV(acc0);
#if defined(DIAGNOSTICS)
		printf("\nParticle %d\n",p-btab+1);
#endif
       for (q = btab; q < btab+nbody; q++) {   
          if (p==q) continue;

          DOTPSUBV(dr2, dr, pos0, Pos(q));
		  VWrapAll(dr);
		  DOTVP(dr2, dr, dr);

		  if (dr2<RcutSq) {
#if defined(DIAGNOSTICS)
		printf("\nInteracting article %d\n",q-btab+1);
#endif
			rri=ssq/dr2; rri3=rri*rri*rri;
			phi_q = fphi*(rri3-1.0)*rri3+eps;
			phi0 += phi_q;
			mr3i = rri3*(rri3-0.5)*rri;
			ADDMULVS(acc0, dr, mr3i);
			virSum += mr3i*dr2;
			nbbcalc += 1;
		  }
       }
       ADDMULVS(Acc(p), acc0, fa);
       Phi(p) += phi0; uSum += Mass(p)*Phi(p);
    }
	virSum = 0.5*virSum*fa;
	uSum = 0.5*uSum;

    cpuforce = cputime() - cpustart;            
}

void ljforcecalc_direct2(bodyptr btab, int nbody)
{
    bodyptr p, q;
    double cpustart;
    vector posp, accp;
    real phip;
    vector dr;
	int k;
	real rri, rri3, uVal, fcVal, dr2;
 
    cpustart = cputime(); 

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    nbbcalc = 0;   
	uSum = 0.0; virSum = 0.0;
    
    for (p = btab; p < btab+nbody-1; p++) {
       for (q = p+1; q < btab+nbody; q++) {   
          DOTPSUBV(dr2, dr, Pos(p), Pos(q));
		  VWrapAll(dr);
		  DOTVP(dr2, dr, dr);

		  if (dr2<RcutSq) {
			rri=1./dr2; rri3=rri*rri*rri;
			uVal = 4.*(rri3-1.0)*rri3+1.;
			Phi(p) += uVal; Phi(q) += uVal;
			fcVal = 48.*rri3*(rri3-0.5)*rri;
			ADDMULVS(Acc(p), dr, fcVal);
			ADDMULVS(Acc(q), dr, -fcVal);
			uSum += uVal;
			virSum += fcVal*dr2;
			nbbcalc += 1;
		  }
       }
    }

    cpuforce = cputime() - cpustart;            
}

// TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE CALCULO DIRECTO DE LA FUERZA


// COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES PARA EL CALCULO DE LA FUERZA

//void ljforcecalc_barnes(void)					// Opcion para solo correr barnes
void ljforcecalc_barnes(bodyptr btab, int nbody)
{
    bodyptr p;
    double cpustart;
    vector rmid;
 
    actlen = FACTIVE * 216 * tdepth;            
    actlen = actlen * rpow(theta, -2.5);        
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    interact = (cellptr) allocate(actlen * sizeof(cell));

    cpustart = cputime();                       
    actmax = nbbcalc = nbccalc = 0;             
    active[0] = (nodeptr) root;                 
    CLRV(rmid);                                 

//    for (p = btab; p < btab+nbody; p++) {
//        Phi(p) = 0.0;
//        CLRV(Acc(p));
//    }
	uSum = 0.0; virSum = 0.0;

    walktree_barnes(active, active + 1, interact, interact + actlen,
             (nodeptr) root, rsize, rmid);      

	virSum = 0.5*virSum*fa;
	uSum = 0.5*uSum;

    cpuforce = cputime() - cpustart;            
    free(active);
    free(interact);
}

local void walktree_barnes(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
    nodeptr *np, *ap, q;
    int actsafe;

    if (Update(p)) {                            
        np = nptr;                              
        actsafe = actlen - NSUB;                
        for (ap = aptr; ap < nptr; ap++)        
            if (Type(*ap) == CELL) {
                if (accept_lj(p, *ap, psize)) { 
                } else {                        
                    if (np - active >= actsafe) 
                        error("walktree: active list overflow\n");
                    for (q = More(*ap); q != Next(*ap); q = Next(q))
						*np++= q;               
                }
            } else                              
                if (*ap != p) {                 
                    --bptr;                     
                    Mass(bptr) = Mass(*ap);     
                    SETV(Pos(bptr), Pos(*ap));
					Id(bptr) = Id(*ap);
#if defined(DIAGNOSTICSQUES)
					if (accept_body(p,*ap)) {
						printf("\n (p,q) = %d %d %d",Id(p),Id(bptr));
					}
#endif
                }
        actmax = MAX(actmax, np - active);      
        if (np != nptr)                         
            walksub(nptr, np, cptr, bptr, p, psize, pmid);
                                                
        else {                                  
            if (Type(p) != BODY)                
                error("walktree: recursion terminated with cell\n");
#if defined(DIAGNOSTICSQUES)
						printf("\n p= %d", (bodyptr)p-bodytab+1);
#endif
            gravsum((bodyptr) p, cptr, bptr);
        }
    }
}

local bool accept_lj(nodeptr p, nodeptr q, real psize)
{
	real drpq, drpq2, rcell;
    vector dr;
	int k;

	DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
	for (k = 0; k < NDIM; k++)
		dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
	DOTVP(drpq2, dr, dr);
	drpq = rsqrt(drpq2);
	rcell= psize * rsqrt((real)(NDIM))/2.0;
	if ( drpq >= Rcut+rcell+Rcrit2(q)/2. )
		return (TRUE);
	else
		return (FALSE);
}

local bool accept_body(nodeptr p, nodeptr q)
{
	real drpq2;
    vector dr;
	int k;

	DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
	for (k = 0; k < NDIM; k++)
		dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
	DOTVP(drpq2, dr, dr);
	if ( drpq2 < RcutSq )
		return (TRUE);
	else
		return (FALSE);
}

local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
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
            walktree_barnes(nptr, np, cptr, bptr, q, psize / 2, nmid);
        }
    } else {
        for (k = 0; k < NDIM; k++)              
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree_barnes(nptr, np, cptr, bptr, p, psize / 2, nmid);
                                                
    }
}
 
local void gravsum(bodyptr p0, cellptr cptr, cellptr bptr)
{
    vector pos0, acc0;
    real phi0;
 
    SETV(pos0, Pos(p0));                        
    phi0 = 0.0;                                 
    CLRV(acc0);

    sumnode_barnes(bptr, interact + actlen, pos0, &phi0, acc0);
                                                
    Phi(p0) = phi0; uSum += Mass(p0)*Phi(p0);                             
	SETV(Acc(p0), acc0);
}

local void sumnode_barnes(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real dr2, phi_p, mr3i;
    vector dr;
	real rri, rri3;
	int k;

    for (p = start; p < finish; p++) {          
		DOTPSUBV(dr2, dr, pos0, Pos(p));
		for (k = 0; k < NDIM; k++)
			dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			rri=1.0/dr2; rri3=rri*rri*rri;
			phi_p = 4.0*(rri3-1.0)*rri3+1.0;
			*phi0 += phi_p;
			mr3i = 48.0*rri3*(rri3-0.5)*rri;
			ADDMULVS(acc0, dr, mr3i);
			virSum += mr3i*dr2;
			nbbcalc += 1;
		}
    }
}

/* TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES PARA EL CALCULO DE LA FUERZA */


/* COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE RAPAPORT PARA EL CALCULO DE LA FUERZA */

void ljforcecalc_cellsmethod(bodyptr btab, int nbody)
{
  bodyptr p, p1, p2;
  vector dr, invWid, rs, shift;
  vectorI cc, m1v, m2v, vOff[] = OFFSET_VALS;
  real fcVal, rr, rri, rri3, uVal;
  int c, j1, j2, m1, m1x, m1y, m1z, m2, n, offset;
  double cpustart;

  cpustart = cputime();                       

  VDiv (invWid, cells, Box);
  for (n = nbody; n < nbody + VProd (cells); n ++) cellList[n] = -1;

  DO_BODY(p,btab,btab+nbody) {
	n = p-btab;
	ADD2VS(rs, Pos(p), Box, 0.5);
    VMul (cc, rs, invWid);
    c = VLinear (cc, cells) + nbody;
    cellList[n] = cellList[c];
    cellList[c] = n;
  }

  DO_BODY(p,btab,btab+nbody) {CLRV(Acc(p)); Phi(p) = 0.0;}
  uSum = 0.;
  virSum = 0.;
  nbbcalc = 0;

#if (NDIM==3)
  for (m1z = 0; m1z < cells[2]; m1z ++) {		// Este loop es 3D
#endif
    for (m1y = 0; m1y < cells[1]; m1y ++) {
      for (m1x = 0; m1x < cells[0]; m1x ++) {
#if (NDIM==3)
        VSet (m1v, m1x, m1y, m1z);				// Vector en 3D
#else
#if (NDIM==2)
        VSet (m1v, m1x, m1y);					// Vector en 2D
#endif
#endif
        m1 = VLinear (m1v, cells) + nbody;
        for (offset = 0; offset < N_OFFSET; offset ++) {
          ADDV(m2v, m1v, vOff[offset]);
          CLRV(shift);
          VCellWrapAll ();
          m2 = VLinear (m2v, cells) + nbody;
          DO_CELL (j1, m1) {
            DO_CELL (j2, m2) {
              if (m1 != m2 || j2 < j1) {
				p1 = btab+j1; p2 = btab+j2;
                SUBV(dr, Pos(p1), Pos(p2));
                VVSub (dr, shift);
                rr = VLenSq (dr);
                if (rr < RcutSq) {
                  rri = 1. / rr;
                  rri3 = Cube (rri);
                  fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                  uVal = 4. * rri3 * (rri3 - 1.) + 1.;
                  VVSAdd(Acc(p1), fcVal, dr);
                  VVSAdd (Acc(p2), - fcVal, dr);
				  Phi(p1) += uVal; Phi(p2) +=uVal;
                  uSum += uVal;
                  virSum += fcVal * rr;
				  nbbcalc+=2;
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
  cpuforce = cputime() - cpustart;            
}

/* TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE CELDAS DE RAPAPORT PARA EL CALCULO DE LA FUERZA */


/* COMIENZA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES2 PARA EL CALCULO DE LA FUERZA */

void ljforcecalc_barnes2(bodyptr btab, int nbody)
{
    bodyptr p;
    double cpustart;
    vector rmid;
 
    actlen = FACTIVE * 216 * tdepth;            
    actlen = actlen * rpow(theta, -2.5);        
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    interact = (cellptr) allocate(actlen * sizeof(cell));

    cpustart = cputime();                       
    actmax = nbbcalc = nbccalc = 0;             
    active[0] = (nodeptr) root;                 
    CLRV(rmid);                                 

    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }
	uSum = 0.0; virSum = 0.0;

    walktree_barnes2(active, active + 1, interact, interact + actlen,
             (nodeptr) root, rsize, rmid);      

    cpuforce = cputime() - cpustart;            
    free(active);
    free(interact);
}

local void walktree_barnes2(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
    nodeptr *np, *ap, q;
    int actsafe;

    if (Update(p)) {                            
        np = nptr;                              
        actsafe = actlen - NSUB;                
        for (ap = aptr; ap < nptr; ap++)        
            if (Type(*ap) == CELL) {
                if (!accept_lj2(p, *ap, psize)) { 
                    if (np - active >= actsafe) 
                        error("walktree: active list overflow\n");
                    for (q = More(*ap); q != Next(*ap); q = Next(q))
						*np++= q;               
                }
            } else                              
                if (*ap != p) {                 
                    --bptr;                     
                    Mass(bptr) = Mass(*ap);     
                    SETV(Pos(bptr), Pos(*ap));
					Id(bptr) = Id(*ap);
                }
        actmax = MAX(actmax, np - active);      
        if (np != nptr)                         
            walksub2(nptr, np, cptr, bptr, p, psize, pmid);
		else {                                  
            if (Type(p) != BODY)                
                error("walktree: recursion terminated with cell\n");
            sumnode_barnes2((bodyptr) p, cptr, bptr);   // sum force on the body
        }
    }
}

local bool accept_lj2(nodeptr p, nodeptr q, real psize) // Se puede unificar con
{														// accept_lj
	real drpq, drpq2, rcell;
    vector dr;
	int k;

	DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
	for (k = 0; k < NDIM; k++)
		dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
	DOTVP(drpq2, dr, dr);
	drpq = rsqrt(drpq2);
	rcell= psize * rsqrt((real)(NDIM))/2.0;
	if ( drpq >= Rcut+rcell+Rcrit2(q)/2. )	// denominator 2 is experimental
		return (TRUE);
	else
		return (FALSE);
}

local void walksub2(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;
 
    poff = psize / 4;                           
    if (Type(p) == CELL)
        for (q = More(p); q != Next(p); q = Next(q)) {
			for (k = 0; k < NDIM; k++)          
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walktree_barnes2(nptr, np, cptr, bptr, q, psize / 2, nmid);
        }
    else {                                    
        for (k = 0; k < NDIM; k++)              
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree_barnes2(nptr, np, cptr, bptr, p, psize / 2, nmid);
                                                
    }
}

local void sumnode_barnes2(bodyptr p, cellptr start, cellptr finish)
{
    cellptr q;
    real dr2, uVal, fVal;
    vector dr;
	real rri, rri3;
	int k;
	bodyptr qb;
 
    for (q = start; q < finish; q++) {
	  qb = ((bodyptr) q);
//	  if ( Id(p) < Id(qb) ) {
		DOTPSUBV(dr2, dr, Pos(p), Pos(q));
		for (k = 0; k < NDIM; k++)
			dr[k]=dr[k] - ( (real) (nint(dr[k]/Box[k])) )*Box[k];
		DOTVP(dr2, dr, dr);

		if (dr2<RcutSq) {
			rri=ssq/dr2; rri3=rri*rri*rri;
			uVal = 4.0*(rri3-1.0)*rri3+1.0;
			fVal = 48.0*rri3*(rri3-0.5)*rri;
			ADDMULVS(Acc(p), dr, fVal);
//			ADDMULVS(Acc(qb), dr, -fVal);
			Phi(p) += uVal;
//			Phi(qb) += uVal;
			virSum += 0.5*fVal*dr2;
			uSum += 0.5*uVal;
			nbbcalc += 1;
		}
//	  }
    }
}

/* TERMINA BLOQUE SIGUIENDO EL ESQUEMA DE BARNES2 PARA EL CALCULO DE LA FUERZA */

