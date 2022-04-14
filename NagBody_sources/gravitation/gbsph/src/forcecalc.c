/* =============================================================================
	MODULE: forcecalc.c					[gbsph]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: Compute forces in the N-body system
	Language: C
	Use:
	Routines and functions:
	Modules, routines and external headers: 
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 24, 2007; March 5, 2008;
	Copyright: (c) 1999-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "switches.h"

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/mathfns.h"
//#include "../../../General_libs/math/vectmath.h"

#include "globaldefs.h"
#include "protodefs.h"

local void normal_gravcalc(bodyptr, int);  
local void direct_gravcalc(bodyptr, int);  
 
local void walktree(nodeptr *, nodeptr *, cellptr, cellptr,
                    nodeptr, real, vector);
local bool accept(nodeptr, real, vector);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local void gravsum(bodyptr, cellptr, cellptr);

local void sumnode_gravity(cellptr, cellptr, vector, real *, vector);
local void sumnode_sf(cellptr, cellptr, vector, real *, vector);
local void sumcell_gravity(cellptr, cellptr, vector, real *, vector);
local void sumcell_gravity_spline(cellptr, cellptr, vector, real *, vector);
local void sumcell_sf(cellptr, cellptr, vector, real *, vector);

local void normal_walktree(bodyptr, nodeptr, real, vector, real *, vector);
local void normal_walktree_spline(bodyptr, nodeptr, real, vector, real *, vector);
local void normal_walktree_sf(bodyptr, nodeptr, real, vector, real *, vector);


// Definiciones para fuerza SPLINE =============================================
#define KERN_LEN   10000  
static real knlrad  [KERN_LEN+1],   
            knlforce[KERN_LEN+1],
            knlpot  [KERN_LEN+1],
            knlW2   [KERN_LEN+1],
            knlW3   [KERN_LEN+1],
            knlW4   [KERN_LEN+1];

local void sumnode_gravity_spline(cellptr, cellptr, vector, real *, vector);


// Definiciones para Movimiento Colectivo Sin Un Lider =========================
local void walktree_collective_motion_without_leader_01(bodyptr, nodeptr,
		real, vector, real *, vector);
local void walktree_collective_motion_without_leader_02(bodyptr, real, nodeptr, 
		real, vector, vector, real *, vector, int *);
local void sumnode_collective_motion_without_leader(cellptr, cellptr, 
		vector, real *, vector);
local int p_num_nearest_neighbours;
// -----------------------------------------------------------------------------

local void darkmatter_walktree(bodyptr, real, nodeptr, 
				real, vector, real *, vector, int *);
local void sum_dm_node(bodyptr, bodyptr, real, real *, vector);

// Definiciones para Potencial Externo =========================================
local real external_pot(real);
local real external_acc_factor(real, real);
local real external_pot_LJ(real);
local real external_acc_factor_LJ(real, real);
local real external_pot_n1_n2(real);
local real external_acc_factor_n1_n2(real, real);
local real external_pot_ALJ(real);
local real external_acc_factor_ALJ(real, real);
local real external_pot_point_halo(real);
local real external_acc_factor_point_halo(real, real);
local real external_pot_ln_halo(real);
local real external_acc_factor_ln_halo(real, real);
local real external_pot_bmn_halo(real);
local real external_acc_factor_bmn_halo(real, real);
// -----------------------------------------------------------------------------

// Definiciones para calculo de Transporte
local void TransportCalc(bodyptr, vector, real, real);
local void TransportCalc2(bodyptr, bodyptr, vector, real, real);
// -----------------------------------------------------------------------------

local bool traslape(bodyptr, real, vector, real);

// Lists of active nodes and interactions...

#if !defined(FACTIVE)
#  define FACTIVE  0.75                         
#endif
 
local int actlen;                               

local nodeptr *active;                          

local cellptr interact;                         

local void forcecalc_method_string_to_int(string,int *);
local int method_int;

#define BARNES		0
#define NULLMETHOD	1
#define NORMAL		2
#define DIRECT		3

void forcecalc(bodyptr btab, int nbody)
{
    forcecalc_method_string_to_int(cmd.forcecalc_method, &method_int);
    switch(method_int) {
        case BARNES:
            gravcalc(); break;
        case NULLMETHOD:
            error("\nNo method was given, we stop ...\n\n"); break;
        case NORMAL:
            normal_gravcalc(btab, nbody); break;
        case DIRECT:
            direct_gravcalc(btab, nbody); break;
        default:
            error("\nUnknown method, we stop ...\n\n"); break;
    }

}

local void forcecalc_method_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"barnes") == 0) *method_int = BARNES;
    if (strnull(method_str)) *method_int = NULLMETHOD;
    if (strcmp(method_str,"normal") == 0) *method_int = NORMAL;
    if (strcmp(method_str,"direct") == 0) *method_int = DIRECT;
}

#undef BARNES
#undef NULLMETHOD
#undef NORMAL
#undef DIRECT


void gravcalc(void)
{
    double cpustart;
    vector rmid;

    actlen = FACTIVE * 216 * gd.tdepth;            
#if !defined(QUICKSCAN)
    actlen = actlen * rpow(cmd.theta, -2.5);        
#endif
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));

    interact = (cellptr) allocate(actlen * sizeof(cell));
    cpustart = cputime();                       
    gd.actmax = gd.nbbcalc = gd.nbccalc = 0;             
    active[0] = (nodeptr) gd.root;                 
    CLRV(rmid);                                 

    walktree(active, active + 1, interact, interact + actlen,
             (nodeptr) gd.root, gd.rsize, rmid);      
    gd.cpuforce = cputime() - cpustart;            
    free(active);
    free(interact);
}

 
local void walktree(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
    nodeptr *np, *ap, q;
    int actsafe;
 
    if (Update(p)) {                            
        np = nptr;                              
        actsafe = actlen - NSUB;                
        for (ap = aptr; ap < nptr; ap++)        
            if (Type(*ap) == CELL) {            
                if (accept(*ap, psize, pmid)) { 
                    Mass(cptr) = Mass(*ap);     
                    dm_Mass(cptr) = dm_Mass(*ap);			// DARKMATTERACC
                    SETV(Pos(cptr), Pos(*ap));
                    SETM(Quad(cptr), Quad(*ap));
                    cptr++;                     
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
                    dm_Mass(bptr) = dm_Mass(*ap);			// DARKMATTERACC
                    SETV(Pos(bptr), Pos(*ap));
                }
        gd.actmax = MAX(gd.actmax, np - active);      
        if (np != nptr)                         
            walksub(nptr, np, cptr, bptr, p, psize, pmid);
                                                
        else {                                  
            if (Type(p) != BODY) {
				printf("Type=%d\n",Type(p));
                error("walktree: recursion terminated with cell\n");
			}
            gravsum((bodyptr) p, cptr, bptr);   
        }
    }
}


// Calculo de la interaccion entre particulas caminando el arbol de manera
// normal. Involucra a la rutina siguiente, la rutina recursiva normal_walktree
// y las dos rutinas comunes, sumnode y sumcell.

local void normal_gravcalc(bodyptr btab, int nbody)
{
    bodyptr p;
    double cpustart;
	vector pos0, acc0;
	real phi0;
 
    cpustart = cputime();                       
    gd.nbbcalc = gd.nbccalc = 0;   
    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

    if (gd.IncludeGrav) {
		if (gd.IncludeGrav_Plummer) {
			for (p = btab; p < btab+nbody; p++) {
				SETV(pos0, Pos(p));
				phi0 = 0.0;
				CLRV(acc0);
				normal_walktree(p, ((nodeptr) gd.root), gd.rsize, pos0,&phi0,acc0);
				Phi(p) += phi0;                       
				ADDMULVS(Acc(p), acc0, 1.0);
			}
		} else {
			for (p = btab; p < btab+nbody; p++) {
				SETV(pos0, Pos(p));
				phi0 = 0.0;
				CLRV(acc0);
				normal_walktree_spline(p, ((nodeptr) gd.root), gd.rsize, pos0,&phi0,acc0);
				Phi(p) += phi0;                       
				ADDMULVS(Acc(p), acc0, 1.0);
			}
		}
	}

    if (gd.IncludeSF_eps) {
		for (p = btab; p < btab+nbody; p++) {
			SETV(pos0, Pos(p));
			phi0 = 0.0;
			CLRV(acc0);
			normal_walktree_sf(p, ((nodeptr) gd.root), gd.rsize, pos0,&phi0,acc0);
			Phi(p) += phi0;                       
			ADDMULVS(Acc(p), acc0, 1.0);
		}
	}

    gd.cpuforce = cputime() - cpustart;            
}

local void normal_walktree(bodyptr p, nodeptr q, real qsize, 
							vector pos0, real *phi0, vector acc0)
{
    nodeptr l;
	real drpq,drpq2;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DISTV(drpq,Pos(p),Pos(q));
				drpq2 = drpq*drpq;
                if ( drpq2 >= Rcrit2(q) ) { 
					if (cmd.usequad)
						sumcell_gravity(((cellptr) q),( (cellptr) q+1),
							pos0,phi0,acc0);
					else 
						sumnode_gravity(((cellptr) q),( (cellptr) q+1),
							pos0,phi0,acc0);
					gd.nbccalc++;
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree(p,l,qsize/2,pos0,phi0,acc0);
					}
				}
			} else {
				sumnode_gravity(((cellptr) q),( (cellptr) q+1),
					pos0,phi0,acc0);
				gd.nbbcalc++;
			}
		}
	}
}

local void normal_walktree_spline(bodyptr p, nodeptr q, real qsize, 
							vector pos0, real *phi0, vector acc0)
{
    nodeptr l;
	real drpq,drpq2;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DISTV(drpq,Pos(p),Pos(q));
				drpq2 = drpq*drpq;
                if ( drpq2 >= Rcrit2(q) ) { 
					if (cmd.usequad)
						sumcell_gravity_spline(((cellptr) q),( (cellptr) q+1),
							pos0,phi0,acc0);
					else
						sumnode_gravity_spline(((cellptr) q),( (cellptr) q+1),
							pos0,phi0,acc0);
					gd.nbccalc++;
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree_spline(p,l,qsize/2,pos0,phi0,acc0);
					}
				}
			} else {
				sumnode_gravity_spline(((cellptr) q),( (cellptr) q+1),
					pos0,phi0,acc0);
				gd.nbbcalc++;
			}
		}
	}
}

local void normal_walktree_sf(bodyptr p, nodeptr q, real qsize, 
								vector pos0, real *phi0, vector acc0)
{
    nodeptr l;
	real drpq,drpq2;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DISTV(drpq,Pos(p),Pos(q));
				drpq2 = drpq*drpq;
                if ( drpq2 >= Rcrit2(q) ) { 
					if (cmd.usequad)
						sumcell_sf(((cellptr) q),( (cellptr) q+1),
							pos0,phi0,acc0);
					else
						sumnode_sf(((cellptr) q),( (cellptr) q+1),pos0,phi0,acc0);
					gd.nbccalc++;
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree_sf(p,l,qsize/2,pos0,phi0,acc0);
					}
				}
			} else {
				sumnode_sf(((cellptr) q),( (cellptr) q+1),pos0,phi0,acc0);
				gd.nbbcalc++;
			}
		}
	}
}

// Terminan rutinas para el c'alculo normal de la interacci'on -----------------


// COMIENZA rutina para el calculo directo de la interaccion entre particulas --

local void direct_gravcalc(bodyptr btab, int nbody)
{
    bodyptr p, q;
    double cpustart;
    vector pos0, acc0;
    real phi0;
    real eps2, dr2, drab, phi_q, mr3i;
    vector dr;
 
    cpustart = cputime();                       

    eps2 = cmd.eps * cmd.eps;                           
	gd.nbbcalc = gd.nbccalc = 0;
    for (p = btab; p < btab+nbody; p++) {
        Phi(p) = 0.0;
        CLRV(Acc(p));
    }

	gd.virSum = 0.0;
	if (cmd.computeTransport) {
		for (p = btab; p < btab+nbody; p++) {
			CLRM(rf(p)); en(p) = 0.0;
		}
	}

    if (gd.IncludeGrav) {
        for (p = btab; p < btab+nbody; p++) {
            SETV(pos0, Pos(p));
            phi0 = 0.0;
            CLRV(acc0);
            for (q = btab; q < btab+nbody; q++) {
                if (p==q) continue;
                DOTPSUBV(dr2, dr, Pos(q), pos0);
                dr2 += eps2;
                drab = rsqrt(dr2);
                phi_q = ( Mass(q) / drab )*cmd.G;
                phi0 -= phi_q;
                mr3i = phi_q / dr2;
                ADDMULVS(acc0, dr, mr3i);
                gd.nbbcalc += 1;

				gd.virSum += mr3i*dr2;
				if (cmd.computeTransport)				// Transport computation ...
					TransportCalc(p, dr, -phi_q, mr3i);
            }
            Phi(p) += phi0;                       
            ADDMULVS(Acc(p), acc0, 1.0);
        }
    }
    if (gd.IncludeSF_eps) {
        for (p = btab; p < btab+nbody; p++) {
            SETV(pos0, Pos(p));
            phi0 = 0.0;
            CLRV(acc0);
            for (q = btab; q < btab+nbody; q++) {
                if (p==q) continue;
                DOTPSUBV(dr2, dr, Pos(q), pos0);
                dr2 += eps2;
                drab = rsqrt(dr2);
                phi_q = ( cmd.dm_alpha*dm_Mass(q)
						* rexp(-drab/cmd.dm_lambda) / drab )*cmd.G;
                phi0 -= phi_q;
                phi_q = ( Mass(q) / drab )*cmd.G;
                mr3i = phi_q * ( cmd.dm_alpha*(dm_Mass(q)/Mass(q))
						*rexp(-drab/cmd.dm_lambda)*(1.0+drab/cmd.dm_lambda) ) / dr2;
                ADDMULVS(acc0, dr, mr3i);
                gd.nbbcalc += 1;
            }
            Phi(p) += phi0;                       
            ADDMULVS(Acc(p), acc0, 1.0);
        }
    }

	gd.virSum = 0.5*gd.virSum;	// Transport ...

    gd.cpuforce = cputime() - cpustart;            
}

// TERMINA rutina para el calculo directo de la interaccion entre particulas ---


#if defined(QUICKSCAN)

local bool accept(nodeptr c, real psize, vector pmid)
{
    real p15, dk;
 
    p15 = ((real) 1.5) * psize;                 
    dk = Pos(c)[0] - pmid[0];                   
    if (ABS(dk) > p15)                          
        return (TRUE);                          
    dk = Pos(c)[1] - pmid[1];                   
    if (ABS(dk) > p15)                          
        return (TRUE);                          
    dk = Pos(c)[2] - pmid[2];                   
    if (ABS(dk) > p15)                          
        return (TRUE);                          
    return (FALSE);                             
}
 
#else
  
local bool accept(nodeptr c, real psize, vector pmid)
{
    real dmax, dsq, dk;
    int k;
 
    dmax = psize;                               
    dsq = 0.0;                                  
    for (k = 0; k < NDIM; k++) {                
        dk = Pos(c)[k] - pmid[k];               
        if (dk < 0)                             
            dk = - dk;
        if (dk > dmax)                          
            dmax = dk;
        dk -= ((real) 0.5) * psize;             
        if (dk > 0)
            dsq += dk * dk;                     
    }
    return (dsq > Rcrit2(c) &&                  
              dmax > ((real) 1.5) * psize);     
}

#endif

local bool traslape(bodyptr p, real qsize, vector qmid, real Hpp)
{
    real p15, dk;
 
    p15 = ((real) 0.5)*(Hpp + qsize); 
    dk = Pos(p)[0] - qmid[0];           
    if (ABS(dk) > p15)             
        return (FALSE);             
#if defined(ONEDIM)
		return(TRUE);
#endif
    dk = Pos(p)[1] - qmid[1];      
    if (ABS(dk) > p15)             
        return (FALSE);             
#if defined(TWODIM)
		return(TRUE);
#endif
    dk = Pos(p)[2] - qmid[2];      
    if (ABS(dk) > p15)             
        return (FALSE);             
    return (TRUE);                
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
            walktree(nptr, np, cptr, bptr, q, psize / 2, nmid);
                                                
        }
    } else {                                    
        for (k = 0; k < NDIM; k++)              
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree(nptr, np, cptr, bptr, p, psize / 2, nmid);
                                                
    }
}

 
local void gravsum(bodyptr p0, cellptr cptr, cellptr bptr)
{
    vector pos0, acc0;
    real phi0;
 
    SETV(pos0, Pos(p0));                        
    phi0 = 0.0;                                 
    CLRV(acc0);                                 

    if (cmd.usequad) {
        if (gd.IncludeGrav) {
			if (gd.IncludeGrav_Plummer)
				sumcell_gravity(interact, cptr, pos0, &phi0, acc0);
			else
				sumcell_gravity_spline(interact, cptr, pos0, &phi0, acc0);
		}
        if (gd.IncludeSF_eps)
            sumcell_sf(interact, cptr, pos0, &phi0, acc0);
    } else {
        if (gd.IncludeGrav)
            if (gd.IncludeGrav_Plummer)
                sumnode_gravity(interact, cptr, pos0, &phi0, acc0);
            else
                sumnode_gravity_spline(interact, cptr, pos0, &phi0, acc0);
        if (gd.IncludeSF_eps)
            sumnode_sf(interact, cptr, pos0, &phi0, acc0);
    }

    if (gd.IncludeGrav)
        if (gd.IncludeGrav_Plummer)
            sumnode_gravity(bptr, interact + actlen, pos0, &phi0, acc0);
        else
            sumnode_gravity_spline(bptr, interact + actlen, pos0, &phi0, acc0);        

    if (gd.IncludeSF_eps)
        sumnode_sf(bptr, interact + actlen, pos0, &phi0, acc0);

    Phi(p0) = phi0;                             
    SETV(Acc(p0), acc0);                        
    gd.nbbcalc += interact + actlen - bptr;        
    gd.nbccalc += cptr - interact;                 

}


// * sumnode - Gravity - PLUMMER 

local void sumnode_gravity(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i;
    vector dr;
 
    eps2 = cmd.eps * cmd.eps;                           
    for (p = start; p < finish; p++) {          
        DOTPSUBV(dr2, dr, Pos(p), pos0);        
        dr2 += eps2;                            
        drab = rsqrt(dr2);
        phi_p = ( Mass(p) / drab )*cmd.G;
        *phi0 -= phi_p;
        mr3i = phi_p / dr2;
        ADDMULVS(acc0, dr, mr3i); 

//		if (cmd.computeTransport)				// Transport computation ...
//			TransportCalc(p, dr, -phi_p, mr3i);
    }
}


// * sumnode - Gravity - SPLINE

local void sumnode_gravity_spline(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real dr2,r,r_inv,fac,u,h,h_inv,h3_inv,ff,wf,wp;
    vector dr;
    int ii;

    h = 2.8*cmd.eps;
    h_inv=1/h;
    h3_inv=h_inv*h_inv*h_inv;

    for (p = start; p < finish; p++) {          

        DOTPSUBV(dr2, dr, Pos(p), pos0);        

        r = rsqrt(dr2);
        u=r*h_inv;

        if(u>=1) {
            r_inv=1/r;
            fac= Mass(p)*r_inv*r_inv*r_inv;
            ADDMULVS(acc0, dr, fac); 
            *phi0 -= Mass(p)*r_inv;
        } else {
            ii = (int)(u*KERN_LEN); 
            ff=(u-knlrad[ii])*KERN_LEN;
            wf=knlforce[ii]+(knlforce[ii+1]-knlforce[ii])*ff;
            wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;

            if(u>1.0e-4) {
                fac=Mass(p)*h3_inv*wf;
                ADDMULVS(acc0, dr, fac); 
            }
            *phi0 += Mass(p)*h_inv*wp;					// Correcto el signo?
        }
    }
}



// * sumnode. Only scalar field potential (smoothed)

local void sumnode_sf(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i;
    vector dr;
 
    eps2 = cmd.eps * cmd.eps;                           

    for (p = start; p < finish; p++) {
        DOTPSUBV(dr2, dr, Pos(p), pos0);
        dr2 += eps2;
        drab = rsqrt(dr2);
        phi_p = ( cmd.dm_alpha*dm_Mass(p) 
				* rexp(-drab/cmd.dm_lambda) / drab )*cmd.G;
        *phi0 -= phi_p;
        phi_p = ( Mass(p) / drab )*cmd.G;
        mr3i = phi_p * ( cmd.dm_alpha*(dm_Mass(p)/Mass(p))
			*rexp(-drab/cmd.dm_lambda)*(1.0+drab/cmd.dm_lambda) ) / dr2;
        ADDMULVS(acc0, dr, mr3i); 
    }
}


// * sumcell - Gravity - PLUMMER

local void sumcell_gravity(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i, drqdr, dr5i, phi_q;
    vector dr, qdr;

    real phi_p_sf;
 
    eps2 = cmd.eps * cmd.eps;

    for (p = start; p < finish; p++) {          
        DOTPSUBV(dr2, dr, Pos(p), pos0);        
        dr2 += eps2;
        drab = rsqrt(dr2);
        phi_p = Mass(p) / drab;
        mr3i = phi_p / dr2;
        DOTPMULMV(drqdr, qdr, Quad(p), dr);     
        dr5i = ((real) 1.0) / (dr2 * dr2 * drab);
        phi_q = ((real) 0.5) * dr5i * drqdr;
        mr3i += ((real) 5.0) * phi_q / dr2;
        *phi0 -= ( phi_p + phi_q )*cmd.G;
        mr3i *= cmd.G;
        ADDMULVS2(acc0, dr, mr3i, qdr, -dr5i);  
    }
}


// * sumcell - Gravity - SPLINE
 
// CHECAR ESTA RUTINA - E INTEGRARLA AL CODIGO 
// SU NOMBRE PREVIO ERA sumcell_gravity ....

local void sumcell_gravity_spline(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real dr2;
    vector dr;

	double r2,dx,dy,dz,r,fac,u,h,h2_inv,h3_inv,h5_inv,h4_inv,h6_inv,ff;
	double q11dx,q12dy,q13dz,q12dx,q22dy,q23dz,q13dx,q23dy,q33dz;
	double r_inv,r2_inv,r3_inv,r5_inv;
	double h_inv;
	double wf,wp,w2,w3,w4,potq;
	int ii;

	h = 2.8*cmd.eps;
	h_inv=1/h;
	h2_inv=h_inv*h_inv;
	h3_inv=h2_inv*h_inv;
	h4_inv=h2_inv*h2_inv;
	h5_inv=h2_inv*h3_inv;
	h6_inv=h3_inv*h3_inv;

    for (p = start; p < finish; p++) {          
		DOTPSUBV(dr2, dr, Pos(p), pos0);        
        r = rsqrt(dr2);
        u=r*h_inv;

        dx = dr[0];
        dy = dr[1];
        dz = dr[2];

		if(u>=1) {
			r_inv=1/r;
			r2_inv=r_inv*r_inv;
			r3_inv=r2_inv*r_inv;
			r5_inv=r2_inv*r3_inv;

			q11dx=QuadQ(p)[0][0]*dx;
			q12dy=QuadQ(p)[0][1]*dy;
			q13dz=QuadQ(p)[0][2]*dz;
			q12dx=QuadQ(p)[0][1]*dx;
			q22dy=QuadQ(p)[1][1]*dy;
			q23dz=QuadQ(p)[1][2]*dz;
			q13dx=QuadQ(p)[0][2]*dx;
			q23dy=QuadQ(p)[1][2]*dy;
			q33dz=QuadQ(p)[2][2]*dz;

			potq = 0.5*(q11dx*dx + q22dy*dy + q33dz*dz) 
					+ q12dx*dy + q13dx*dz + q23dy*dz;

			*phi0 += -Mass(p)*r_inv + r3_inv*( -3*potq*r2_inv + 0.5*QuadP(p));  

			fac= Mass(p)*r3_inv + (15*potq*r2_inv -1.5*QuadP(p) )*r5_inv;

			acc0[0]+=dx*fac;
			acc0[1]+=dy*fac;
			acc0[2]+=dz*fac;

			ff=-3*r5_inv;
			acc0[0] += ff*(q11dx + q12dy + q13dz);
			acc0[1] += ff*(q12dx + q22dy + q23dz);
			acc0[2] += ff*(q13dx + q23dy + q33dz);
		} else {
			ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
			wf=knlforce[ii] + (knlforce[ii+1]-knlforce[ii])*ff;
			wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;
			w2=knlW2[ii]    + (knlW2[ii+1]-knlW2[ii])*ff;
			w3=knlW3[ii]    + (knlW3[ii+1]-knlW3[ii])*ff;
			w4=knlW4[ii]    + (knlW4[ii+1]-knlW4[ii])*ff;

			q11dx=QuadQ(p)[0][0]*dx;
			q12dy=QuadQ(p)[0][1]*dy;
			q13dz=QuadQ(p)[0][2]*dz;
			q12dx=QuadQ(p)[0][1]*dx;
			q22dy=QuadQ(p)[1][1]*dy;
			q23dz=QuadQ(p)[1][2]*dz;
			q13dx=QuadQ(p)[0][2]*dx;
			q23dy=QuadQ(p)[1][2]*dy;
			q33dz=QuadQ(p)[2][2]*dz;

			potq = 0.5*(q11dx*dx + q22dy*dy + q33dz*dz)
					+ q12dx*dy + q13dx*dz + q23dy*dz;

			*phi0 += Mass(p)*h_inv*wp + potq*w2*h5_inv + 0.5*QuadP(p)*wf*h2_inv*h_inv;

			if(u>1.0e-4) {
				r_inv=1/r;
				fac= Mass(p)*h2_inv*h_inv*wf + 
					+ potq*h6_inv * w3*r_inv + 0.5*QuadP(p) * w4 *h4_inv*r_inv; 

				acc0[0]+=dx*fac;
				acc0[1]+=dy*fac;
				acc0[2]+=dz*fac;

				ff=w2*h5_inv;
				acc0[0] += ff*(q11dx + q12dy + q13dz);
				acc0[1] += ff*(q12dx + q22dy + q23dz);
				acc0[2] += ff*(q13dx + q23dy + q33dz);
		    }
		}
    }
}


// * sumcell with scalar field // NO ES LA VERSION SUAVIZADA ... NO USAR!!...

local void sumcell_sf(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i, drqdr, dr5i, phi_q;
    vector dr, qdr;

    real phi_p_sf;
 
    eps2 = cmd.eps * cmd.eps;
    for (p = start; p < finish; p++) {          
        DOTPSUBV(dr2, dr, Pos(p), pos0);        
        dr2 += eps2;
        drab = rsqrt(dr2);
        phi_p = Mass(p) / drab;
        phi_p_sf = cmd.dm_alpha*dm_Mass(p) * rexp(-drab/cmd.dm_lambda) / drab;
        mr3i = phi_p * ( cmd.dm_alpha*(dm_Mass(p)/Mass(p))
				*rexp(-drab/cmd.dm_lambda)*(1.0+drab/cmd.dm_lambda) ) / dr2;
        *phi0 -= ( phi_p_sf )*cmd.G;
        mr3i *= cmd.G;
        ADDMULVS(acc0, dr, mr3i);  
    }
}


// COMIENZAN DEFINICIIONES POTENCIALEXTERNO ------------------------------------

/*
Estas rutinas tiene el proposito de calcular el efecto de un potencial
central externo en el sistema de N-cuerpos. Se supone que el origen del
potencial esta en el origen del sistema.

En un futuro se debera considerar incluir a este potencial como una 
"particula" miembro del sistema. Esto nos dara la posibilidad de
generalizar al caso de varios centros de dispersion o pozos de potencial,
para aplicaciones en la dinamica de sistemas astrofisicos inmersos en 
una distribuciones de pozos de potencial, producto del colapso de
campos escalares.

Agregamos dos rutinas: external_pot_point_halo, external_acc_factor_point_halo.
para tomar en cuenta el potencial de un halo oscuro puntual descrito para
comenzar por un modelo de Dehnen con gamma=0 y centrado en el origen. [Nov. 19, 2004].
*/

void external_forcecalc(bodyptr btab, int nbody)
{
    bodyptr p;
    real dr2, drab, mr3i;
    vector dr, pos0;
    double cpustart, cpu_ext_pot;

    cpustart = cputime();                       
 
    for (p = btab; p < btab+nbody; p++) {
        DOTPSUBV(dr2, dr, gd.pos_pot, Pos(p));	// dr esta va en la direcci'on de pos_pot
        drab = rsqrt(dr2);
		Phi(p) += external_pot(drab);			// es atractivo si pot es negativo
        mr3i = external_acc_factor(drab,dr2)/Mass(p);
        ADDMULVS(Acc(p), dr, mr3i);			// la particula se mueve hacia pos_pot
											// cuando pot es negativo
    }

    cpu_ext_pot = cputime() - cpustart;            
    printf("\ncpu_ext_pot= %8.3f\n",cpu_ext_pot);
}


/*
Se considera para comenzar, como potencial externo, uno de Lennard-Jones

v(r) = 4e [ (s/r)^12 - (s/r)^6  ]

s -- es una longitud efectiva
e -- es una energia tipica
*/

local real external_pot(real rp)
{    
   real pot;
//   pot = external_pot_LJ(rp);
//   pot = external_pot_ALJ(rp);
//   pot = external_pot_n1_n2(rp);
//   pot = external_pot_point_halo(rp);
   pot = external_pot_ln_halo(rp);
   return pot;
}

local real external_acc_factor(real drab, real dr2)
{
   real acc;
//   acc=external_acc_factor_LJ(drab,dr2);
//   acc=external_acc_factor_ALJ(drab,dr2);
//   acc=external_acc_factor_n1_n2(drab,dr2);
//   acc=external_acc_factor_point_halo(drab,dr2); 
//   acc=external_acc_factor_ln_halo(drab,dr2); 
   acc=external_acc_factor_bmn_halo(drab,dr2); 
   return acc;
}


local real external_pot_LJ(real rp)
{    
	real pot;
	pot = 4.0*cmd.eps_pot
			*( rpow((cmd.sigma_pot/rp), 12.0) - rpow((cmd.sigma_pot/rp), 6.0) );
	return pot;
}

local real external_acc_factor_LJ(real drab, real dr2)
{
   real acc;
   acc = -4.0*cmd.eps_pot
		*( 12.0*rpow((cmd.sigma_pot/drab), 12.0) - 6.0*rpow((cmd.sigma_pot/drab), 6.0) )
		/ dr2;
   return acc;
}

#define n1	12.0
#define n2	6.0
#define sigma_pot_1	1.0
#define sigma_pot_2	1.0

local real external_pot_n1_n2(real rp)
{    
	real pot;
	pot = 4.0*cmd.eps_pot
			*( rpow((sigma_pot_1/rp), n1) - rpow((sigma_pot_2/rp), n2) );
	return pot;
}

local real external_acc_factor_n1_n2(real drab, real dr2)
{
   real acc;
   acc = -4.0*cmd.eps_pot
		*( n1*rpow((sigma_pot_1/drab), n1) - n2*rpow((sigma_pot_2/drab), n2) )
		/ dr2;
   return acc;
}

#undef n1
#undef n2
#undef sigma_pot_1
#undef sigma_pot_2

local real external_pot_ALJ(real rp)
{    
	real pot;
	real rmin;

	rmin=rpow(2.0, 1.0/6.0) * cmd.sigma_pot;
	if (rp < rmin)
		pot=-cmd.eps_pot;
	else
	pot = external_pot_LJ(rp);
	return pot;
}

local real external_acc_factor_ALJ(real drab, real dr2)
{
	real acc;
	real rmin;

	rmin=rpow(2.0, 1.0/6.0) * cmd.sigma_pot;
	if (drab<rmin)
		acc=0.0;
	else
		acc=external_acc_factor_LJ(drab,dr2);
	return acc;
}

/* Comienzan definiciones del potencial punto tipo Halo oscuro. 
	Consideramos el potencial puntual del modelo de Dehnen
	con gamma=0. Tambien, estamos considerando que G=1. 

	Phi(r) = -(GM/a) (1/(2-gamma)) [1 - (r/(r+a))^(2-gamma)] 	

	gamma neq 2.
*/

#define GAMMA_POT 0

// Para usar los parametros en linea, eps_pot y sigma_pot hacemos la 
//	correspondencia MH_POT -> eps_pot y AH_POT -> sigma_pot

local real external_pot_point_halo(real rp)
{    
	real pot;
	pot = -(cmd.eps_pot/cmd.sigma_pot)*(1./(2.-GAMMA_POT))
			*( 1.0 - rpow(rp/(rp+cmd.sigma_pot),2.0-GAMMA_POT) );
	return pot;
}

local real external_acc_factor_point_halo(real drab, real dr2)
{
	real acc;
	acc = cmd.eps_pot*rpow(drab,2.0-GAMMA_POT)
			/ rpow((drab+cmd.sigma_pot), 3-GAMMA_POT)  / dr2;
	return acc;
}

#undef GAMMA_POT

// Terminan definiciones del potencial punto tipo Halo oscuro ------------------

// Comienzan definiciones del potencial tipo Halo de campo escalar colapsado Ln --
local real external_pot_ln_halo(real drab)
{
	real pot;
	pot = cmd.eps_pot*rlog( rsqr(drab) + rsqr(cmd.sigma_pot) ) ;
	return pot;
}

local real external_acc_factor_ln_halo(real drab, real dr2)
{
	real acc;
	acc = 2.0*cmd.eps_pot/(dr2+rsqr(cmd.sigma_pot));
	return acc;
}
// Terminan definiciones del potencial tipo Halo de campo escalar colapsado Ln --

// Comienzan definiciones del potencial tipo Halo de campo escalar colapsado bmn --
#define PHI02		5.0e-5
local real external_pot_bmn_halo(real drab)
{
	real pot;
	real Phi0, k, omega, t, r;
	
	Phi0 = rsqrt(PHI02);
	k=2.0; omega=20.001; t=1.0;
	r = drab;

	pot = -(
			rpow(Phi0,2)
			*(
				rpow(omega,2)*rpow(rsin(k*r),2)
				+ rpow(rcos(t*omega),2)
			*(
				-(rpow(k,2)*rpow(rsin(k*r),2)) 
				+ rpow(-rcos(k*r) + rsin(k*r)/(k*r),2)
			)
			)
			)
			/(2.*rpow(r,2));
	
//	-( 
//			rpow(Phi0,2) * (
//			 rpow(omega,2)*rpow(rsin(k*r),2) + rpow(rcos(omega*t), 2) * (
//						-(rpow(k,2)*rpow(rsin(k*r),2)) 
//						+ rpow(-rcos(k*r)+rsin(k*r)/(k*r),2)
//					)
//				)
//			)/(2.0*rpow(r,2));
	return pot;
}

local real external_acc_factor_bmn_halo(real drab, real dr2)
{
	real acc;
	real Phi0, k, omega, t, r;

	Phi0 = rsqrt(PHI02);
	k=2.0; omega=20.001; t=1.0;
	r = drab;

	acc = -(
			rpow(Phi0,2)
			*(
				2.0*k*rpow(omega,2)*rcos(k*r)*rsin(k*r) 
				+ rpow(rcos(t*omega),2)
			*(
				-2*rpow(k,3)*rcos(k*r)*rsin(k*r) 
				+ 
				2*(
				rcos(k*r)/r + k*rsin(k*r) - rsin(k*r)/(k*rpow(r,2))
				)
				*(-rcos(k*r) + rsin(k*r)/(k*r))
				)
				)
			)/(2.*rpow(r,2)) 
			+ 
			(
				rpow(Phi0,2)
				*(
					rpow(omega,2)*rpow(rsin(k*r),2) 
					+ rpow(rcos(t*omega),2)
					*(
						-(rpow(k,2)*rpow(rsin(k*r),2)) 
						+ rpow(-rcos(k*r) + rsin(k*r)/(k*r),2)
					)
				)
			)
			/rpow(r,3);
	return acc;
}
// Terminan definiciones del potencial tipo Halo de campo escalar colapsado bmn --

// TERMINAN DEFINICIIONES POTENCIALEXTERNO -------------------------------------


// * Setting kernel table .... sumnode_gravity_spline...

void force_setkernel(void) 
{
  int i;
  double u;

  for(i=0;i<=KERN_LEN;i++)
    {
      u=((double)i)/KERN_LEN;

      knlrad[i] = u;

      if(u<=0.5)
	{
	  knlforce[i]=32 * (1.0/3 -6.0/5*pow(u,2) + pow(u,3));
	  knlpot[i]=16.0/3*pow(u,2)-48.0/5*pow(u,4)+32.0/5*pow(u,5)-14.0/5;

	  knlW2[i]= -384.0/5 +96.0*u;
	  knlW3[i]= 96.0;
	  knlW4[i]= 96.0/5*u*(5*u-4);
	}
      else
	{
	  knlforce[i]=64*(1.0/3 -3.0/4*u + 
		  3.0/5*pow(u,2)-pow(u,3)/6) - 1.0/15/pow(u,3);
	  knlpot[i]=1.0/15/u +32.0/3*pow(u,2)-16.0*pow(u,3)+
		  48.0/5*pow(u,4)-32.0/15*pow(u,5)-16.0/5;

	  knlW2[i]= 384.0/5 + 1/(5.0*pow(u,5)) -48.0/u -32*u;
	  knlW3[i]= -32-1/pow(u,6)+48/pow(u,2);
	  knlW4[i]= -48 +1/(5*pow(u,4)) +384.0/5*u -32*pow(u,2);
	}
    }
}


// COMIENZAN rutinas para el calculo del movimiento colectivo sin un lider =====

/*
AUN FALTA TERMINAR ESTAS RUTINAS.... ESTAN EN DESARROLLO USENSE CON CUIDADO!!!!!!!

1. Las rutinas de recorrido del arbol ya han sido mejoradas, puede usarse una o la otra:

a. walktree_collective_motion_without_leader_01
b. walktree_collective_motion_without_leader_02

Pero hay que unificarlas en terminos de las variables y parametros que usan.

2. Hay que checar en si la acelaracion es un promedio de las
velocidades de las particulas que rodean a una dada. Sucede que se
la velocidad de las particulas va disminuyendo con el tiempo, hasta
que su velocidad es cero. El promedio que hemos implementado es:

vx = (1/Num. Vecinos) Sum_i [ v_x_i]

vy = (1/Num. Vecinos) Sum_i [ v_y_i]

vz = (1/Num. Vecinos) Sum_i [ v_z_i]

Se deben checar estas tres rutinas de abajo y stepsystem.
*/

#define VIEWLENGTH	0.1

void forcecalc_collective_motion_without_leader(bodyptr btab, int nbody)
{
    bodyptr p;
    double cpustart;
	vector pos0, acc0;
	real phi0, mr3i;

// Por considerar en vez de las anteriores lineas
/*	int nndm;
	vector rmid;
	real lambda_p; */
// Termina la consideracion ...

    cpustart = cputime();                       
    gd.nbbcalc = 0;

    for (p = btab; p < btab+nbody; p++) {
        p_num_nearest_neighbours = 0;
        SETV(pos0, Pos(p));
        phi0 = 0.0;
        CLRV(Acc(p));
        CLRV(acc0);
        walktree_collective_motion_without_leader_01(p, ((nodeptr) gd.root),
				gd.rsize, pos0,&phi0,acc0);

// Por considerar en vez de las anteriores lineas ...
/*		lambda_p = dm_lambda*rsize;
		CLRV(rmid);
		phi0 = 0.0;
		CLRV(acc0);
		nndm=0;
        walktree_collective_motion_without_leader_02(p, lambda_p, ((nodeptr) root), rsize, pos0, rmid,
			&phi0, acc0, &nndm); */
// Termina la consideracion...

        Phi(p) = phi0;                       
//        SETV(Acc(p), acc0);
        mr3i = 1.0/((real) p_num_nearest_neighbours);
        ADDMULVS(Acc(p), acc0, mr3i);
//        DIVVS(Acc(p), acc0, ((real) p_num_nearest_neighbours));
//        printf("%d -- nearest neighbours = %d\n", p-btab+1, p_num_nearest_neighbours);
    }
    gd.cpuforce = cputime() - cpustart;            
}


local void walktree_collective_motion_without_leader_01(bodyptr p, nodeptr q,
		real qsize, vector pos0, real *phi0, vector acc0)
{
    nodeptr l;
    real drpq;
    real length;

    if (Update(p)) {
        if (Type(q) == CELL) {
            DISTV(drpq,Pos(p),Pos(q));
            length=VIEWLENGTH+qsize*rsqrt(3.0)/2.0;
            if ( drpq <= length )
                for (l = More(q); l != Next(q); l = Next(l))
                    walktree_collective_motion_without_leader_01(p,l,qsize/2,pos0,phi0,acc0);
        } else {
            DISTV(drpq,Pos(p),Pos(q));
            if ( drpq <= VIEWLENGTH ) {
                sumnode_collective_motion_without_leader(((cellptr) q),( (cellptr) q+1),
                                    pos0,phi0,acc0);
                gd.nbbcalc++;
            }
        }
    }
}

local void walktree_collective_motion_without_leader_02(bodyptr p, real hpp, nodeptr q, 
				real qsize, vector pos0, vector qmid, real *phi0, vector acc0, int *nndm)
{
    nodeptr l;
    real drpq;
    real poff;
    int k;
    vector nmid;
  
    if (Type(q) == CELL) {            
        if (traslape(p, qsize, qmid, hpp) ) { 
            poff = qsize / 4; 
            for (l = More(q); l != Next(q); l = Next(l)) {
                for (k = 0; k < NDIM; k++) 
                    nmid[k] = qmid[k] + (Pos(l)[k] < qmid[k] ? - poff : poff);
                walktree_collective_motion_without_leader_02(p,hpp,l,qsize/2,pos0, nmid,phi0,acc0,nndm);
            }
        }
    } else {
        DISTV(drpq,Pos(p),Pos(q));
        if ( (drpq <= hpp) ) {
            sumnode_collective_motion_without_leader(((cellptr) q),( (cellptr) q+1),pos0,phi0,acc0);
            (*nndm)++;
        }
        return;
    }
}

local void sumnode_collective_motion_without_leader(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real mr3i;
 
    for (p = start; p < finish; p++) {          
        *phi0 -= 0.0;
        mr3i = 1.0;
        ADDMULVS(acc0, Vel(p), mr3i);               
        p_num_nearest_neighbours++;
    }
}

#undef VIEWLENGTH

//TERMINAN rutinas para el calculo del movimiento colectivo sin un lider -------

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

