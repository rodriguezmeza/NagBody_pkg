/* ==============================================================================
!	MODULE: treeload.c															!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date:																!
!	Purpose: Lennard-Jones force computation (Several methods)					!
!	Language: C																	!
!	Use: maketree();															!
!	Routines and functions:	tree_ljforce										!
!	External modules, routines and headers:										!
!	Comments and notes:															!
!	Info: M.A. Rodriguez-Meza,													!
!		Depto. de Fisica, ININ,													!
!		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico.							!
!		e-mail: marioalberto@inin.gob.mx
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
#include "data_struc_defs.h"
#include "global_defs.h"
#include "proto_defs.h"

local void newtree(void);
local cellptr makecell(void);
local void expandbox(bodyptr, int);
local void loadbody(bodyptr);
local int subindex(bodyptr, cellptr);
local void hackcofm(cellptr, real, int);
local void setrcrit(cellptr, vector, real);
local void threadtree(nodeptr, nodeptr);
local void hackquad(cellptr);					// Not in use by now
 
local bool bh86, sw94;
local nodeptr freecell = NULL;

#define MAXLEVEL  32  

local int cellhist[MAXLEVEL];
local int subnhist[MAXLEVEL];

void maketree(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = cputime();
    newtree();
    root = makecell();
    CLRV(Pos(root));
    expandbox(btab, nbody);
	DO_BODY(p, btab, btab+nbody)
        loadbody(p);                           
    bh86 = scanopt(options, "bh86");           
    sw94 = scanopt(options, "sw94");
    if (bh86 && sw94)                          
        error("maketree: incompatible options bh86 and sw94\n");
    tdepth = 0;                                
    for (i = 0; i < MAXLEVEL; i++)             
        cellhist[i] = subnhist[i] = 0;
    hackcofm(root, rsize, 0);                  
    threadtree((nodeptr) root, NULL);          
    if (usequad)                               
        hackquad(root);                        // Not in use by now
    cputree = cputime() - cpustart;            
}

local void newtree(void)
{
    static bool firstcall = TRUE;
    nodeptr p;
 
    if (! firstcall) {                          
        p = (nodeptr) root;                     
        while (p != NULL)                       
            if (Type(p) == CELL) {              
                Next(p) = freecell;             
                freecell = p;                  
                p = More(p);                    
            } else                             
                p = Next(p);                    
    } else                                      
        firstcall = FALSE;                     
    root = NULL;                               
    ncell = 0;                                  
}

local cellptr makecell(void)
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
    ncell++;                                    
    return (c);                                 
}

local void expandbox(bodyptr btab, int nbody)
{
    real dmax, d;
    bodyptr p;
    int k;
 
    dmax = 0.0;                                 
	DO_BODY(p, btab, btab+nbody)
		DO_COORD(k) {
            d = rabs(Pos(p)[k] - Pos(root)[k]); 
            if (d > dmax)                      
                dmax = d;                       
        }
    while (rsize < 2 * dmax)                   
        rsize = 2 * rsize;                      
}

local void loadbody(bodyptr p)
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;
 
    q = root;                                   
    qind = subindex(p, q);                      
    qsize = rsize;                              
    while (Subp(q)[qind] != NULL) {            
        if (Type(Subp(q)[qind]) == BODY) {      
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0)			
                error("loadbody: two bodies have same position\n");
            c = makecell();                     
			DO_COORD(k)
                Pos(c)[k] = Pos(q)[k] +         
                    (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;
            Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
                                               
            Subp(q)[qind] = (nodeptr) c;        
        }
        q = (cellptr) Subp(q)[qind];            
        qind = subindex(p, q);                  
        qsize = qsize / 2;                      
    }
    Subp(q)[qind] = (nodeptr) p;                
}

local int subindex(bodyptr p, cellptr q)
{
    int ind, k;
 
    ind = 0;                                    
	DO_COORD(k)
        if (Pos(q)[k] <= Pos(p)[k])             
            ind += NSUB >> (k + 1);             
    return (ind);
}

local void hackcofm(cellptr p, real psize, int lev)
{
    vector cmpos, tmpv;
    int i, k;
    nodeptr q;
 
    tdepth = MAX(tdepth, lev);                  
    cellhist[lev]++;                            
    Mass(p) = 0.0;                              
    CLRV(cmpos);                                
    for (i = 0; i < NSUB; i++)                  
        if ((q = Subp(p)[i]) != NULL) {         
            subnhist[lev]++;                    
            if (Type(q) == CELL)                
                hackcofm((cellptr) q, psize/2, lev+1);
                                                
            Update(p) |= Update(q);             
            Mass(p) += Mass(q);                 
            MULVS(tmpv, Pos(q), Mass(q));       
            ADDV(cmpos, cmpos, tmpv);           
        }
    if (Mass(p) > 0.0) {                        
        DIVVS(cmpos, cmpos, Mass(p));           
    } else {                                    
        SETV(cmpos, Pos(p));                    
    }
	DO_COORD(k)
        if (cmpos[k] < Pos(p)[k] - psize/2 ||   
              Pos(p)[k] + psize/2 <= cmpos[k])
            error("hackcofm: tree structure error\n");

    setrcrit(p, cmpos, psize);                  
    SETV(Pos(p), cmpos);                        
}


local void setrcrit(cellptr p, vector cmpos, real psize)
{
    real bmax2, d;
    int k;

    if (theta == 0.0)                           
        Rcrit2(p) = rsqr(2 * rsize);            
    else if (sw94) {                            
        bmax2 = 0.0;                            
		DO_COORD(k) {
            d = cmpos[k] - Pos(p)[k] + psize/2; 
            bmax2 += rsqr(MAX(d, psize - d));   
        }
        Rcrit2(p) = bmax2 / rsqr(theta);        
    } else if (bh86)                            
        Rcrit2(p) = rsqr(psize / theta);        
    else {
        DISTV(d, cmpos, Pos(p));                
//      Rcrit2(p) = rsqr(psize / theta + d);				// original line
//      Rcrit2(p) = rsqr(psize / theta);					// Activa cuando metodo barnes
//		Rcrit2(p) = rsqr( psize ) * ((real) NDIM)/4.0;		// Square of the Radius of 
															// the sphere containing the cell
//		Rcrit2(p) = rsqr(psize * rsqrt((real)(NDIM))/2.0);	// Usado en met. barnes new
		Rcrit2(p) = psize * rsqrt((real)(NDIM))/2.0 ;		// Usado en met. normal2
    }
}


local void threadtree(nodeptr p, nodeptr n)
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

