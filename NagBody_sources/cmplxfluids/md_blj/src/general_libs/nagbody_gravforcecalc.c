#include "machines.h"						// Verificar su uso...

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"

#include "nagbody.h"

//#include "globaldefs.h"
//#include "protodefs.h"

local void sumnode(bodyptr, cellptr, cellptr, vector, real *, vector,
	global_data_treegrav *gdtreegrav);
local void sumcell(bodyptr, cellptr, cellptr, vector, real *, vector,
	global_data_treegrav *gdtreegrav);
local void normal_walktree(bodyptr, nodeptr, real, vector, real *, vector,
			global_data_treegrav *gdtreegrav);

void ind_normal_gravcalc(bodyptr btab, int nbody, bodyptr p, 
	global_data_treegrav *gdtreegrav)
{
    double cpustart;
	vector pos0, acc0;
	real phi0;
 
    cpustart = cputime();                       
    gdtreegrav->nbbcalc = gdtreegrav->nbccalc = 0;   

	SETV(pos0, Pos(p));
	phi0 = 0.0;
	CLRV(acc0);
//	if (scanopt(cmd.data_analysis_type, "acf-anim")) {
	if (scanopt(gdtreegrav->options, "acf-anim")) {
		CLRM(rf(p));
		en(p) = 0.0;
	}
	normal_walktree(p, ((nodeptr) gdtreegrav->root), gdtreegrav->rsize, 
		pos0,&phi0,acc0, gdtreegrav);
	Phi(p) = phi0;                       
	SETV(Acc(p), acc0);                  
//	if (scanopt(cmd.data_analysis_type, "acf-anim")) {
	if (scanopt(gdtreegrav->options, "acf-anim")) {
		en(p) *= Mass(p);
	}
    gdtreegrav->cpuforce += cputime() - cpustart;            
}

void normal_gravcalc(bodyptr btab, int nbody, global_data_treegrav *gdtreegrav)
{
    bodyptr p;
    double cpustart;
	vector pos0, acc0;
	real phi0;
 
    cpustart = cputime();                       
    gdtreegrav->nbbcalc = gdtreegrav->nbccalc = 0;   

    for (p = btab; p < btab+nbody; p++) {
		SETV(pos0, Pos(p));
		phi0 = 0.0;
		CLRV(acc0);
		normal_walktree(p, ((nodeptr) gdtreegrav->root), gdtreegrav->rsize, 
			pos0,&phi0,acc0, gdtreegrav);
		Phi(p) = phi0;                       
		SETV(Acc(p), acc0);                  
	}
    gdtreegrav->cpuforce = cputime() - cpustart;            
}

local void normal_walktree(bodyptr p, nodeptr q, real qsize, 
							vector pos0, real *phi0, vector acc0,
							global_data_treegrav *gdtreegrav)
{
    nodeptr l;
	real drpq,drpq2;

	if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
			if (Type(q) == CELL) {
				DISTV(drpq,Pos(p),Pos(q));
				drpq2 = drpq*drpq;
                if ( drpq2 >= Rcrit2(q) ) { 
					if (gdtreegrav->usequad)
						sumcell(p, ((cellptr) q),( (cellptr) q+1),
							pos0,phi0,acc0, gdtreegrav);
					else 
						sumnode(p, ((cellptr) q),( (cellptr) q+1),
							pos0,phi0,acc0, gdtreegrav);
					gdtreegrav->nbccalc++;
				} else {
					for (l = More(q); l != Next(q); l = Next(l)) {
						normal_walktree(p,l,qsize/2,pos0,phi0,acc0, gdtreegrav);
					}
				}
			} else {
					sumnode(p, ((cellptr) q),( (cellptr) q+1),
						pos0,phi0,acc0, gdtreegrav);
					gdtreegrav->nbbcalc++;
			}
		}
	}
}

local void sumnode(bodyptr p, cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0,
				   global_data_treegrav *gdtreegrav)
{
    cellptr q;
    real eps2, dr2, drab, phi_q, mr3i;
    vector dr;
	matrix w;
	int k;
 
//    eps2 = gdtreegrav->eps * gdtreegrav->eps;                           

    for (q = start; q < finish; q++) {          
        DOTPSUBV(dr2, dr, Pos(q), pos0);        
        dr2 += gdtreegrav->eps2;                            
        drab = rsqrt(dr2);   
        phi_q = ( Mass(q) / drab );
        *phi0 -= phi_q;
        mr3i = phi_q / dr2;
        ADDMULVS(acc0, dr, mr3i); 
//		if (scanopt(cmd.data_analysis_type, "acf-anim")) {
		if (scanopt(gdtreegrav->options, "acf-anim")) {
			en(p) -= phi_q;
			OUTVP(w, dr, dr);
			MULMS(w, w, mr3i);
			DO_COORD(k) {
				VVAdd(rf(p)[k], w[k]);
			}
		}
    }
}
 
local void sumcell(bodyptr p0, cellptr start, cellptr finish,
                   vector pos0, real *phi0, vector acc0,
				   global_data_treegrav *gdtreegrav)
{
    cellptr p;
    real eps2, dr2, drab, phi_p, mr3i, drqdr, dr5i, phi_q;
    vector dr, qdr;
 
//    eps2 = gdtreegrav->eps * gdtreegrav->eps;

    for (p = start; p < finish; p++) {          
        DOTPSUBV(dr2, dr, Pos(p), pos0);        
        dr2 += gdtreegrav->eps2;
        drab = rsqrt(dr2);
        phi_p = Mass(p) / drab;
        mr3i = phi_p / dr2;
        DOTPMULMV(drqdr, qdr, Quad(p), dr);     
        dr5i = ((real) 1.0) / (dr2 * dr2 * drab);
        phi_q = ((real) 0.5) * dr5i * drqdr;
        mr3i += ((real) 5.0) * phi_q / dr2;
        *phi0 -= ( phi_p + phi_q );
        ADDMULVS2(acc0, dr, mr3i, qdr, -dr5i);
    }
}

