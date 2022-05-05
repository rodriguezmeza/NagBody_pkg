/*==============================================================================
	MODULE: nagbody_tree.c			[NagBody]
	Written by: M.A. Rodriguez-Meza
	Starting date:	January, 2005
	Purpose: Routines to drive input and output data
	Language: C
	Use:
	Routines and functions:
	External modules, routines and headers:	stdinc.h, mathfns.h, vectmath
						vectmath.h, getparam.h
						types.h, stat.h, inout.h
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
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

#include "../general/stdinc.h"
#include "../math/mathfns.h"
#include "../io/inout.h"
#include "../math/vectdefs.h"
#include "../math/vectmath.h"
#include "nagbody.h"
#include "../physics/physconstants.h"

#include <string.h>
// #include <strings.h>							// For unix
//#include "../../../../General_libs/strings.h"	// For Visual C

#include <sys/stat.h>

//--------------------------- COMIENZA BLOQUE ARBOL ----------------------------

local void newtree(global_data_tree *);
local cellptr makecell(global_data_tree *);
local void expandbox(bodyptr, int, global_data_tree *);
local void loadbody(bodyptr, global_data_tree *);
local int subindex(bodyptr, cellptr);
local void hackcofm(cellptr, real, int, global_data_tree *, 
					global_data_tree_bljforcecalc *);
local void setrcrit(cellptr, vector, real, global_data_tree *);
local void threadtree(nodeptr, nodeptr);
local void hackquad(cellptr);					// Not in use by now
 
local bool bh86, sw94;
local nodeptr freecell = NULL;

#define MAXLEVEL  32  

local int cellhist[MAXLEVEL];
local int subnhist[MAXLEVEL];
														// CHECK 2D --- OK!!
void maketree(bodyptr btab, int nbody, char *options, 
	global_data_tree *gdtree, global_data_tree_bljforcecalc *gdforce)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = cputime();

    newtree(gdtree);
    gdtree->root = makecell(gdtree);
    CLRV(Pos(gdtree->root));
	SETV( GPos(gdtree->root), Pos(gdtree->root) );	// Set Gemetric Center...

    expandbox(btab, nbody, gdtree);

	DO_BODY(p, btab, btab+nbody)
        loadbody(p, gdtree);

    bh86 = scanopt(options, "bh86");
    sw94 = scanopt(options, "sw94");
    if (bh86 && sw94)
        error("maketree: incompatible options bh86 and sw94\n");
    gdtree->tdepth = 0;
    for (i = 0; i < MAXLEVEL; i++)             
        cellhist[i] = subnhist[i] = 0;
    hackcofm(gdtree->root, gdtree->rsize, 0, gdtree, gdforce);                  
    threadtree((nodeptr) gdtree->root, NULL);          
    if (gdtree->usequad)	// Si esta corriendo md_blj_tree inicializar usequad
        hackquad(gdtree->root);							// Not in use by now
    gdtree->cputree = cputime() - cpustart;
}

local void newtree(global_data_tree *gdtree)			// CHECK 2D --- OK!!
{
    static bool firstcall = TRUE;
    nodeptr p;
 
    if (! firstcall) {                          
        p = (nodeptr) gdtree->root;                     
        while (p != NULL)                       
            if (Type(p) == CELL) {              
                Next(p) = freecell;             
                freecell = p;                  
                p = More(p);                    
            } else                             
                p = Next(p);                    
    } else                                      
        firstcall = FALSE;                     
    gdtree->root = NULL;                               
    gdtree->ncell = 0;                                  
}

local cellptr makecell(global_data_tree *gdtree)		// CHECK 2D --- OK!!
{
    cellptr c;
    int i;
 
    if (freecell == NULL)                       
        c = (cellptr) allocate(sizeof(cell));   
    else {                                     
        c = (cellptr) freecell;                 
        freecell = Next(c);                     
    }
    Type(c) = CELL;                             
    Update(c) = FALSE;                          
    for (i = 0; i < NSUB; i++)                  
        Subp(c)[i] = NULL;                     
    gdtree->ncell++;                                    
    return (c);                                 
}
														// CHECK 2D --- OK!!
local void expandbox(bodyptr btab, int nbody, global_data_tree *gdtree)
{
    real dmax, d;
    bodyptr p;
    int k;
 
    dmax = 0.0;                                 
	DO_BODY(p, btab, btab+nbody)
		DO_COORD(k) {
            d = rabs(Pos(p)[k] - Pos(gdtree->root)[k]); 
            if (d > dmax)                      
                dmax = d;                       
        }
    while (gdtree->rsize < 2 * dmax)                   
        gdtree->rsize = 2 * gdtree->rsize;                      
}

local void loadbody(bodyptr p, global_data_tree *gdtree)	// CHECK 2D --- OK!!
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;
 
    q = gdtree->root;                                   
    qind = subindex(p, q);                      
    qsize = gdtree->rsize;
    while (Subp(q)[qind] != NULL) {
        if (Type(Subp(q)[qind]) != CELL) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0)
                error("loadbody: two bodies have same position\n");
            c = makecell(gdtree);
			DO_COORD(k)
                Pos(c)[k] = Pos(q)[k] +         
                    (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;

			SETV( GPos(c), Pos(c) );		// Set GC ...

			Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
                                               
            Subp(q)[qind] = (nodeptr) c;        
        }
        q = (cellptr) Subp(q)[qind];            
        qind = subindex(p, q);                  
        qsize = qsize / 2;                      
    }
    Subp(q)[qind] = (nodeptr) p;                
}

local int subindex(bodyptr p, cellptr q)				// CHECK 2D --- OK!!
{
    int ind, k;
 
    ind = 0;                                    
	DO_COORD(k)
        if (Pos(q)[k] <= Pos(p)[k])             
            ind += NSUB >> (k + 1);             
    return (ind);
}
														// CHECK 2D --- OK!!
local void hackcofm(cellptr p, real psize, int lev, global_data_tree *gdtree, 
					global_data_tree_bljforcecalc *gdforce)
{
    vector cmpos, tmpv;
    int i, k;
    nodeptr q;
	vector cmvel;					// Agregado para contener el num. de cuerpos
	vector dr;
	real drpq2, drpq;

    gdtree->tdepth = MAX(gdtree->tdepth, lev);
    cellhist[lev]++;
    Mass(p) = 0.0;
	NBodies(p) = 0;					// Agregado para contener el num. de cuerpos
    CLRV(cmpos);
    CLRV(cmvel);					// Agregado para contener el num. de cuerpos
	Rcut(p) = psize * rsqrt((real)(NDIM))/2.0;
    for (i = 0; i < NSUB; i++)
        if ((q = Subp(p)[i]) != NULL) {         
            subnhist[lev]++;                    
            if (Type(q) == CELL) {
                hackcofm((cellptr) q, psize/2, lev+1, gdtree, gdforce);
				DOTPSUBV(drpq2, dr, GPos(p), GPos(q));
			} else
				DOTPSUBV(drpq2, dr, GPos(p), Pos(q));

			DO_COORD(k)
				dr[k]=dr[k]
						-((real)(nint(dr[k]/gdforce->Box[k])))*gdforce->Box[k];
			DOTVP(drpq2, dr, dr);
			drpq = rsqrt(drpq2);
			if (drpq+Rcut(q)>Rcut(p))
				Rcut(p)=drpq+Rcut(q);

            Update(p) |= Update(q);
            Mass(p) += Mass(q);
			NBodies(p) += NBodies(q);	// Agregado para cont el num. de cuerpos
            MULVS(tmpv, Pos(q), Mass(q));
            ADDV(cmpos, cmpos, tmpv);
            MULVS(tmpv, Vel(q), Mass(q)); // Agregado para cont num. de cuerpos
            ADDV(cmvel, cmvel, tmpv);	 // Agregado para cont num. de cuerpos
        }
    if (Mass(p) > 0.0) {                        
        DIVVS(cmpos, cmpos, Mass(p));
        DIVVS(cmvel, cmvel, Mass(p));	// Agregado para cont num. de cuerpos
    } else {
        SETV(cmpos, Pos(p));
        SETV(cmvel, Vel(p));
    }
	DO_COORD(k)
        if (cmpos[k] < Pos(p)[k] - psize/2 ||   
              Pos(p)[k] + psize/2 <= cmpos[k])
            error("hackcofm: tree structure error\n");

    setrcrit(p, cmpos, psize, gdtree);                  
    SETV(Pos(p), cmpos);
    SETV(Vel(p), cmvel);				// Agregado para cont num. de cuerpos
}
														// CHECK 2D --- OK!!
local void setrcrit(cellptr p, vector cmpos, real psize, global_data_tree *gdtree)
{
    real bmax2, d;
    int k;

    if (gdtree->theta == 0.0)
        Rcrit2(p) = rsqr(2 * gdtree->rsize);  // ES CORRECTO EN 2D? ...
    else if (sw94) {
        bmax2 = 0.0;
		DO_COORD(k) {
            d = cmpos[k] - Pos(p)[k] + psize/2;
            bmax2 += rsqr(MAX(d, psize - d));
        }
        Rcrit2(p) = bmax2 / rsqr(gdtree->theta);
    } else if (bh86)
        Rcrit2(p) = rsqr(psize / gdtree->theta);
    else {
		Rcrit2(p) = psize * rsqrt((real)(NDIM))/2.0 ;	// Usado en met. normal2
    }
}

local void threadtree(nodeptr p, nodeptr n)				// CHECK 2D --- OK!!
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    Next(p) = n;
    if (Type(p) == CELL) {                      
        ndesc = 0;                              
        for (i = 0; i < NSUB; i++)              
            if (Subp(p)[i] != NULL)             
                desc[ndesc++] = Subp(p)[i];     
        More(p) = desc[0];                      
        desc[ndesc] = n;                       
        for (i = 0; i < ndesc; i++)             
            threadtree(desc[i], desc[i+1]);     
    }
}
														// CHECK 2D --- OK!!
local void hackquad(cellptr p)				// Not in use by now
{
    int ndesc, i;
    nodeptr desc[NSUB], q;
    vector dr;
    real drsq;
    matrix drdr, Idrsq, tmpm;
 
    ndesc = 0;                                  
    for (i = 0; i < NSUB; i++)                  
        if (Subp(p)[i] != NULL)                 
            desc[ndesc++] = Subp(p)[i];         
    CLRM(Quad(p));                              
    for (i = 0; i < ndesc; i++) {               
        q = desc[i];                            
        if (Type(q) == CELL)                    
            hackquad((cellptr) q);              
        SUBV(dr, Pos(q), Pos(p));               
        OUTVP(drdr, dr, dr);                    
        DOTVP(drsq, dr, dr);                   
        SETMI(Idrsq);                           
        MULMS(Idrsq, Idrsq, drsq);              
        MULMS(tmpm, drdr, 3.0);                
        SUBM(tmpm, tmpm, Idrsq);                
        MULMS(tmpm, tmpm, Mass(q));             
        if (Type(q) == CELL)                    
            ADDM(tmpm, tmpm, Quad(q));          
        ADDM(Quad(p), Quad(p), tmpm);           
    }
}

//---------------------------- TERMINA BLOQUE ARBOL ----------------------------


