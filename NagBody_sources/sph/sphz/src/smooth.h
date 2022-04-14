//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)


#ifndef _smooth_h
#define _smooth_h


#define MAXLEV   20
#define NCOEFS    5
#define NRPROF   11
#define EXTLIST  10


typedef struct _pqnode {
    real pqkey;
    int pind;
    struct _pqnode *pqLoser;
    struct _pqnode *pqFromInt;
    struct _pqnode *pqFromExt;
    struct _pqnode *pqWinner;
} pqnode;


typedef struct {
    kdxptr kd;
    int  nsmooth;
    real coefs[NCOEFS];
    int  rprof[NRPROF];
    real freqmax;
    real freqsum;
    pqnode *pque;
    pqnode *pqhead;
    real *r2ball;
    real *r2list;
    int *inlist;
    real cpustart;
} smcontext, *smxptr;


#define WSmooth(wsm, wsc, x2, cf)					\
{									\
    real _x = rsqrt(x2);                                                \
    if (x2 < 1)                                                         \
        wsm = wsc * (cf[0] + cf[1]*_x + cf[2]*x2 + cf[3]*_x*x2);        \
    else                                                                \
        wsm = wsc * cf[4] * (2-_x)*(2-_x)*(2-_x);                       \
}

#define dWSmooth(dsm, dsc, x2, cf)					\
{					                                \
    real _x = rsqrt(x2);                                                \
    if (x2 < 1)                                                         \
        dsm = dsc * (cf[1] + 2*cf[2]*_x + 3*cf[3]*x2);                  \
    else                                                                \
        dsm = dsc * cf[4] * -3 * (2-_x)*(2-_x);                         \
}


smxptr init_smooth(kdxptr, int, real);
void smooth(smxptr, void (*)(smxptr, int, int));
void resmooth(smxptr, void (*)(smxptr, int, int));
void report_smooth(smxptr, stream, string);
void finish_smooth(smxptr);

#endif  // ! _smooth_h
