//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)

 
#include "../../../General_libs/zeno/clib/stdinc.h"
#include "../../../General_libs/zeno/clib/mathfns.h"
#include "../../../General_libs/zeno/clib/vectmath.h"
#include "sphdefs.h"
 
 
local void walkgrav(nodeptr *, nodeptr *, cellptr, cellptr,
		    nodeptr, real, vector);
local bool accept(nodeptr, real, vector);
local void walksub(nodeptr *, nodeptr *, cellptr, cellptr,
                   nodeptr, real, vector);
local void gravsum(bodyptr, cellptr, cellptr);
local void sumnode(cellptr, cellptr, vector, real *, vector);
local void sumcell(cellptr, cellptr, vector, real *, vector);
 
 
#if !defined(FACTIVE)
#define FACTIVE  1
#endif
 
local int actlen;

local nodeptr *active;
local cellptr interact;

 
void gravforce(void)
{
    double cpustart;
    vector rmid;
 
    actlen = FACTIVE * 216 * tdepth;
#if !defined(QUICKSCAN)
    actlen = actlen * rpow(theta, -2.5);
#endif
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    interact = (cellptr) allocate(actlen * sizeof(cell));
    cpustart = cputime();
    actmax = nfcalc = nbbcalc = nbccalc = 0;
    active[0] = (nodeptr) root;
    CLRV(rmid);
    walkgrav(active, active + 1, interact, interact + actlen,
             (nodeptr) root, rsize, rmid);
    cpuforce = cputime() - cpustart;
    free(active);
    free(interact);
}


void report_force(stream ostr, int nbody)
{
    fprintf(ostr, "\n %9s %9s %9s %9s %13s %13s %9s\n",
	    "rsize", "tdepth", "ftree",
	    "actmax", "nbbint", "nbcint", "CPUhfc");
    fprintf(ostr, " %9.1f %9d %9.3f %9d %13d %13d %9.3f\n",
	    rsize, tdepth, (nbody + ncell - 1) / ((real) ncell),
	    actmax, nbbcalc, nbccalc, cpuforce);
    fflush(ostr);
}

 
local void walkgrav(nodeptr *aptr, nodeptr *nptr, cellptr cptr, cellptr bptr,
                    nodeptr p, real psize, vector pmid)
{
    nodeptr *np, *ap, q;
    int actsafe;
 
    if (Update(p)) {
	np = nptr;
	actsafe = actlen - NSUB;
	for (ap = aptr; ap < nptr; ap++)
	    if (Cell(*ap)) {
		if (accept(*ap, psize, pmid)) {
		    Mass(cptr) = Mass(*ap);
		    SETV(Pos(cptr), Pos(*ap));
		    SETM(Quad(cptr), Quad(*ap));
		    cptr++;
		} else {
		    if (np - active >= actsafe)
			error("walkgrav: active list overflow\n");
		    for (q = More(*ap); q != Next(*ap); q = Next(q))

			*np++= q;
		}
	    } else
		if (*ap != p) {
		    --bptr;
		    Mass(bptr) = Mass(*ap);
		    SETV(Pos(bptr), Pos(*ap));
		}
	actmax = MAX(actmax, np - active);
	if (np != nptr)
	    walksub(nptr, np, cptr, bptr, p, psize, pmid);

	else {
	    if (! Body(p))
		error("walkgrav: recursion terminated with cell\n");
	    gravsum((bodyptr) p, cptr, bptr);
	}
    }
}

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

 
local void walksub(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
                   nodeptr p, real psize, vector pmid)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;
 
    poff = psize / 4;
    if (Cell(p)) { 
        for (q = More(p); q != Next(p); q = Next(q)) {

            for (k = 0; k < NDIM; k++)
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walkgrav(nptr, np, cptr, bptr, q, psize / 2, nmid);

        }
    } else {
        for (k = 0; k < NDIM; k++)
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walkgrav(nptr, np, cptr, bptr, p, psize / 2, nmid);

    }
}

 
local void gravsum(bodyptr p0, cellptr cptr, cellptr bptr)
{
    vector pos0, acc0;
    real phi0;
 
    SETV(pos0, Pos(p0));
    phi0 = 0.0;
    CLRV(acc0);
    if (usequad)
        sumcell(interact, cptr, pos0, &phi0, acc0);

    else
        sumnode(interact, cptr, pos0, &phi0, acc0);

    sumnode(bptr, interact + actlen, pos0, &phi0, acc0);

    Phi(p0) = phi0;
    SETV(Acc(p0), acc0);
    nfcalc++;
    nbbcalc += interact + actlen - bptr;
    nbccalc += cptr - interact;
}

 
local void sumnode(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i;
    vector dr;
 
    eps2 = eps * eps;
    for (p = start; p < finish; p++) {
        DOTPSUBV(dr2, dr, Pos(p), pos0);

        dr2 += eps2;
        drab = rsqrt(dr2);
        phi_p = Mass(p) / drab;
        *phi0 -= phi_p;
        mr3i = phi_p / dr2;
        ADDMULVS(acc0, dr, mr3i);
    }
}

 
local void sumcell(cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i, drqdr, dr5i, phi_q;
    vector dr, qdr;
 
    eps2 = eps * eps;
    for (p = start; p < finish; p++) {
        DOTPSUBV(dr2, dr, Pos(p), pos0);
        dr2 += eps2;
        drab = rsqrt(dr2);
        phi_p = Mass(p) / drab;
        mr3i = phi_p / dr2;
        DOTPMULMV(drqdr, qdr, Quad(p), dr);
        dr5i = ((real) 1.0) / (dr2 * dr2 * drab);
        phi_q = ((real) 0.5) * dr5i * drqdr;
        *phi0 -= phi_p + phi_q;
        mr3i += ((real) 5.0) * phi_q / dr2;
        ADDMULVS2(acc0, dr, mr3i, qdr, -dr5i);
    }
}
