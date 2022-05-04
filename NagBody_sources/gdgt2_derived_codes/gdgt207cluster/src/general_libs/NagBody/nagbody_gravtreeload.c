
#include "machines.h"						// Verificar su uso...

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"

#include "nagbody.h"

//#include "globaldefs.h"
//#include "protodefs.h"

#include <sys/stat.h>
 
local void newtree(global_data_treegrav *);
local cellptr makecell(global_data_treegrav *);
local void expandbox(bodyptr, int, global_data_treegrav *);
local void loadbody(bodyptr,bodyptr,global_data_treegrav *);
local int subindex(bodyptr, cellptr);
local void hackcofm(cellptr, real, int, global_data_treegrav *);
local void setrcrit(cellptr, vector, real, global_data_treegrav *);
local void threadtree(nodeptr, nodeptr);
local void hackquad(cellptr);
 
local bool bh86, sw94;                          
local nodeptr freecell = NULL;                  

#define MAXLEVEL  32                            

local int cellhist[MAXLEVEL];                   
local int subnhist[MAXLEVEL];                   

 
void maketree_grav(bodyptr btab, int nbody, global_data_treegrav *gdtreegrav)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = cputime();
    newtree(gdtreegrav);
    gdtreegrav->root = makecell(gdtreegrav);

    CLRV(Pos(gdtreegrav->root));
    expandbox(btab, nbody, gdtreegrav);

    for (p = btab; p < btab+nbody; p++)
        loadbody(p,btab,gdtreegrav);

    bh86 = scanopt(gdtreegrav->options, "bh86");
    sw94 = scanopt(gdtreegrav->options, "sw94");
    if (bh86 && sw94)
        error("maketree: incompatible options bh86 and sw94\n");
    gdtreegrav->tdepth = 0;
    for (i = 0; i < MAXLEVEL; i++)
        cellhist[i] = subnhist[i] = 0;
    hackcofm(gdtreegrav->root, gdtreegrav->rsize, 0, gdtreegrav);
    threadtree((nodeptr) gdtreegrav->root, NULL);
    if (gdtreegrav->usequad)
        hackquad(gdtreegrav->root);
    gdtreegrav->cputree = cputime() - cpustart;
}
 
local void newtree(global_data_treegrav *gdtreegrav)
{
    static bool firstcall = TRUE;
    nodeptr p;
 
    if (! firstcall) {
        p = (nodeptr) gdtreegrav->root;
        while (p != NULL)
            if (Type(p) == CELL) {              
                Next(p) = freecell;             
                freecell = p;                   
                p = More(p);                    
            } else                              
                p = Next(p);                    
    } else                                      
        firstcall = FALSE;                      
		
//firstcall = TRUE;					// AGREGADO PARA RECUPERAR MEMORIA ...
									// CHECAR CUIDADOSAMENTE Y BORRAR !!!!!!!!!!

    gdtreegrav->root = NULL;
    gdtreegrav->ncell = 0;
}

local cellptr makecell(global_data_treegrav *gdtreegrav)
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
    gdtreegrav->ncell++;
    return (c);
}
 
local void expandbox(bodyptr btab, int nbody, global_data_treegrav *gdtreegrav)
{
    real dmax, d;
    bodyptr p;
    int k;
 
    dmax = 0.0;
    for (p = btab; p < btab+nbody; p++)
        for (k = 0; k < NDIM; k++) {
            d = rabs(Pos(p)[k] - Pos(gdtreegrav->root)[k]);
            if (d > dmax)
                dmax = d;
        }
    while (gdtreegrav->rsize < 2 * dmax)
        gdtreegrav->rsize = 2 * gdtreegrav->rsize;
}

local void loadbody(bodyptr p, bodyptr btab, global_data_treegrav *gdtreegrav)
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;

    q = gdtreegrav->root;
    qind = subindex(p, q);                      
    qsize = gdtreegrav->rsize;
    while (Subp(q)[qind] != NULL) {
        if (Type(Subp(q)[qind]) == BODY) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0)
                error("loadbody: two bodies have same position ( %d, %d ) %lf %lf %lf %lf %lf %lf\n",
				p-btab+1, ((bodyptr) Subp(q)[qind])-btab+1, Pos(p)[0], Pos(p)[1], Pos(p)[1],
				Pos(Subp(q)[qind])[0], Pos(Subp(q)[qind])[1], Pos(Subp(q)[qind])[2]);
            c = makecell(gdtreegrav);
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

local void hackcofm(cellptr p, real psize, int lev, 
	global_data_treegrav *gdtreegrav)
{
    vector cmpos, tmpv;
    int i, k;
    nodeptr q;

    gdtreegrav->tdepth = MAX(gdtreegrav->tdepth, lev);                  
    cellhist[lev]++;                            
    Mass(p) = 0.0;                              
    CLRV(cmpos);                                
    for (i = 0; i < NSUB; i++)
        if ((q = Subp(p)[i]) != NULL) {         
            subnhist[lev]++;                    
            if (Type(q) == CELL)
                hackcofm((cellptr) q, psize/2, lev+1, gdtreegrav);
                                                
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
    setrcrit(p, cmpos, psize, gdtreegrav);
    SETV(Pos(p), cmpos);
}

local void setrcrit(cellptr p, vector cmpos, real psize, global_data_treegrav *gdtreegrav)
{
    real bmax2, d;
    int k;

    if (gdtreegrav->theta == 0.0)                           
        Rcrit2(p) = rsqr(2 * gdtreegrav->rsize);
    else if (sw94) {
        bmax2 = 0.0;
        for (k = 0; k < NDIM; k++) {
            d = cmpos[k] - Pos(p)[k] + psize/2;
            bmax2 += rsqr(MAX(d, psize - d));
        }
        Rcrit2(p) = bmax2 / rsqr(gdtreegrav->theta);
    } else if (bh86)                            
        Rcrit2(p) = rsqr(psize / gdtreegrav->theta);        
    else {                                      
        DISTV(d, cmpos, Pos(p));                
        Rcrit2(p) = rsqr(psize / gdtreegrav->theta + d);    
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

local void hackquad(cellptr p)
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

