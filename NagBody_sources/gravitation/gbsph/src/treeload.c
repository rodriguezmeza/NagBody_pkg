/* =============================================================================
	MODULE: treeload.c					[gbsph]
	Written by: M.A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Make the tree of the N-body system
	Language: C
	Use: 'maketree();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 24, 2007;
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

//#include "treedefs.h"

//#include <sys/stat.h>

// Incluido para depurar dm_hackcofm
#include "globaldefs.h"
#include "protodefs.h"

//#include "forcedefs.h"

 
local void newtree(void);                       
local cellptr makecell(void);                   
local void expandbox(bodyptr, int);             
local void loadbody(bodyptr,bodyptr);           
local int subindex(bodyptr, cellptr);           
local void hackcofm(cellptr, real, int);        
//#if defined(DARKMATTERACC)
local void dm_hackcofm(cellptr, real, int);      
//#endif  
local void setrcrit(cellptr, vector, real);     
local void threadtree(nodeptr, nodeptr);        
local void hackquad(cellptr);                   
 
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
    gd.root = makecell();                          

    CLRV(Pos(gd.root));                            
    expandbox(btab, nbody);                     

    for (p = btab; p < btab+nbody; p++)         
        loadbody(p,btab);                       

    bh86 = scanopt(cmd.options, "bh86");            
    sw94 = scanopt(cmd.options, "sw94");            
    if (bh86 && sw94)                           
        error("maketree: incompatible options bh86 and sw94\n");
    gd.tdepth = 0;                                 
    for (i = 0; i < MAXLEVEL; i++)              
        cellhist[i] = subnhist[i] = 0;
    hackcofm(gd.root, gd.rsize, 0);

    if (scanopt(cmd.force_models, "scalar_field_potential_eps") 
		|| scanopt(cmd.force_models, "scalar_field_potential_no_eps") ) {
printf("Entering dm_hackcofm: (rsize) %8.3f\n",gd.rsize);	// Se puede borrar esta linea de impresion
        dm_hackcofm(gd.root, gd.rsize, 0);  
printf("Saliendo dm_hackcofm: (dm_Mass y ncell) %8.3f %d\n",dm_Mass(gd.root),gd.ncell);	// Se puede borrar esta linea de impresion
    }


    threadtree((nodeptr) gd.root, NULL);           
    if (cmd.usequad)                                
        hackquad(gd.root);                         
    gd.cputree = cputime() - cpustart;             
}

 
local void newtree(void)
{
    static bool firstcall = TRUE;
    nodeptr p;
 
    if (! firstcall) {                          
        p = (nodeptr) gd.root;                     
        while (p != NULL)                       
            if (Type(p) == CELL) {              
                Next(p) = freecell;             
                freecell = p;                   
                p = More(p);                    
            } else                              
                p = Next(p);                    
    } else                                      
        firstcall = FALSE;                      
    gd.root = NULL;                                
    gd.ncell = 0;                                  
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
    gd.ncell++;                                    
    return (c);                                 
}

 
local void expandbox(bodyptr btab, int nbody)
{
    real dmax, d;
    bodyptr p;
    int k;
 
    dmax = 0.0;                                 
    for (p = btab; p < btab+nbody; p++)         
        for (k = 0; k < NDIM; k++) {            
            d = rabs(Pos(p)[k] - Pos(gd.root)[k]); 
            if (d > dmax)                       
                dmax = d;                       
        }
    while (gd.rsize < 2 * dmax)                    
        gd.rsize = 2 * gd.rsize;                      
}

 
local void loadbody(bodyptr p, bodyptr btab)
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;
 
    q = gd.root;                                   
    qind = subindex(p, q);                      
    qsize = gd.rsize;                              
    while (Subp(q)[qind] != NULL) {             
        if (Type(Subp(q)[qind]) == BODY) {      
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0)			
                error("loadbody: two bodies have same position ( %d, %d ) %lf %lf %lf %lf %lf %lf\n",
				p-btab+1, ((bodyptr) Subp(q)[qind])-btab+1, Pos(p)[0], Pos(p)[1], Pos(p)[1],
				Pos(Subp(q)[qind])[0], Pos(Subp(q)[qind])[1], Pos(Subp(q)[qind])[2]);
            c = makecell();             
            for (k = 0; k < NDIM; k++)          
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
    for (k = 0; k < NDIM; k++)                  
        if (Pos(q)[k] <= Pos(p)[k])             
            ind += NSUB >> (k + 1);             
    return (ind);
}

 
local void hackcofm(cellptr p, real psize, int lev)
{
    vector cmpos, tmpv;
    int i, k;
    nodeptr q;
	
    gd.tdepth = MAX(gd.tdepth, lev);                  
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
    for (k = 0; k < NDIM; k++)                  
        if (cmpos[k] < Pos(p)[k] - psize/2 ||   
              Pos(p)[k] + psize/2 <= cmpos[k])
            error("hackcofm: tree structure error\n");
#if !defined(QUICKSCAN)
    setrcrit(p, cmpos, psize);                  
#endif
    SETV(Pos(p), cmpos);                        
}


local void dm_hackcofm(cellptr p, real psize, int lev)
{
    int i, k;
    nodeptr q;
    real dr2, drab;

    vector dr;
	
    gd.tdepth = MAX(gd.tdepth, lev);                  
    cellhist[lev]++;                            
    dm_Mass(p) = 0.0;                              
    for (i = 0; i < NSUB; i++)                  
        if ((q = Subp(p)[i]) != NULL) {         
            subnhist[lev]++;                    
            if (Type(q) == CELL)                
                dm_hackcofm((cellptr) q, psize/2, lev+1);

            Update(p) |= Update(q);             
            DOTPSUBV(dr2, dr, Pos(p), Pos(q));                                                 
//            dr2 += eps*eps;
            drab = rsqrt(dr2);
	    if (drab==0.) {
	    	dm_Mass(p) += Mass(q);
	    } else {
            	dm_Mass(p) += Mass(q) 
								* sinh(drab/cmd.dm_lambda)/(drab/cmd.dm_lambda);
	    }                
        }
}


#if !defined(QUICKSCAN)


local void setrcrit(cellptr p, vector cmpos, real psize)
{
    real bmax2, d;
    int k;

    if (cmd.theta == 0.0)                           
        Rcrit2(p) = rsqr(2 * gd.rsize);            
    else if (sw94) {                            
        bmax2 = 0.0;                            
        for (k = 0; k < NDIM; k++) {            
            d = cmpos[k] - Pos(p)[k] + psize/2; 
            bmax2 += rsqr(MAX(d, psize - d));   
        }
        Rcrit2(p) = bmax2 / rsqr(cmd.theta);        
    } else if (bh86)                            
        Rcrit2(p) = rsqr(psize / cmd.theta);        
    else {                                      
        DISTV(d, cmpos, Pos(p));                
        Rcrit2(p) = rsqr(psize / cmd.theta + d);    
    }
}

#endif

 
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

 
local void hackquad(cellptr p)
{
    int ndesc, i;
    nodeptr desc[NSUB], q;
    vector dr;
    real drsq;
    matrix drdr, Idrsq, tmpm;
	matrix tmpmq;						// Added for Spline
	real tmpp;
 
    ndesc = 0;                                  
    for (i = 0; i < NSUB; i++)                  
        if (Subp(p)[i] != NULL)                 
            desc[ndesc++] = Subp(p)[i];         
    CLRM(Quad(p));
    CLRM(QuadQ(p));								// Added for Spline
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

        MULMS(tmpmq, tmpm, Mass(q));             // Added for Spline
		tmpp = drsq*Mass(q);					// Added for Spline

        SUBM(tmpm, tmpm, Idrsq);                
        MULMS(tmpm, tmpm, Mass(q));
        if (Type(q) == CELL) {
            ADDM(tmpm, tmpm, Quad(q));
            ADDM(tmpmq, tmpmq, QuadQ(q));		// Added for Spline
			tmpp += QuadP(q);					// Added for Spline
		}
        ADDM(Quad(p), Quad(p), tmpm);
        ADDM(QuadQ(p), QuadQ(p), tmpmq);		// Added for Spline
		QuadP(p) += tmpp;						// Added for Spline
    }
}

