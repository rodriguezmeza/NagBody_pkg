//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)

 
#ifndef _sphdefs_h
#define _sphdefs_h
 

typedef struct _node {
    byte type;
    byte flags;
    byte curlevel;
    byte newlevel;
    struct _node *next;   
    real mass;
    vector pos;
} node, *nodeptr;
 

#define Type(x)      (((nodeptr) (x))->type)
#define Flags(x)     (((nodeptr) (x))->flags)
#define CurLevel(x)  (((nodeptr) (x))->curlevel)
#define NewLevel(x)  (((nodeptr) (x))->newlevel)
#define Next(x)      (((nodeptr) (x))->next)
#define Mass(x)      (((nodeptr) (x))->mass)
#define Pos(x)       (((nodeptr) (x))->pos)


#define CELL  001
#define BODY  002
#define GAS   004
#define STAR  010

#define Cell(x)  ((Type(x) & CELL) != 0)
#define Body(x)  ((Type(x) & BODY) != 0)
#define Gas(x)   ((Type(x) & GAS) != 0)
#define Star(x)  ((Type(x) & STAR) != 0)


#define INCLUDE  001
#define UPDATE   002
#define INQUE    004
#define DONE     010
#define SURFACE  020

#define Include(x)  ((Flags(x) & INCLUDE) != 0)
#define Update(x)   ((Flags(x) & UPDATE) != 0)
#define InQue(x)    ((Flags(x) & INQUE) != 0)
#define Done(x)     ((Flags(x) & DONE) != 0)
#define Surface(x)  ((Flags(x) & SURFACE) != 0)

#define SetFlag(x,f)  (Flags(x) |= (f))
#define ClrFlag(x,f)  (Flags(x) &= ~(f))

 
typedef struct {
    node bodynode;
    vector vel;
    vector vmid;
    vector acc;
    real smooth;
    real phi;
    real rho;
#if defined(ENTROPY)
    real entf;
#else
    real uint;
#endif
    real udotint;
#if defined(RADIATING)
    real udotrad;
#endif
#if defined(COMPVISC)
    real udotvis;
#endif
    real press;
    real frequency;
#if defined(DIFFUSING) || defined(OPAQUE)
    real tau;
#endif
#if defined(STARFORM)
    real birth;
#endif
} body, *bodyptr;


#define Vel(x)       (((bodyptr) (x))->vel)
#define Vmid(x)      (((bodyptr) (x))->vmid)
#define Acc(x)       (((bodyptr) (x))->acc)
#define Phi(x)       (((bodyptr) (x))->phi)
#define Smooth(x)    (((bodyptr) (x))->smooth)
#define Rho(x)       (((bodyptr) (x))->rho)
#define EntFunc(x)   (((bodyptr) (x))->entf)
#define Uintern(x)   (((bodyptr) (x))->uint)
#define UdotInt(x)   (((bodyptr) (x))->udotint)
#define UdotRad(x)   (((bodyptr) (x))->udotrad)
#define UdotVis(x)   (((bodyptr) (x))->udotvis)
#define Press(x)     (((bodyptr) (x))->press)
#define Frequency(x) (((bodyptr) (x))->frequency)
#define Tau(x)       (((bodyptr) (x))->tau)
#define Birth(x)     (((bodyptr) (x))->birth)

#define NthBody(bp,n)  ((bp) + (n))

  
#define NSUB (1 << NDIM)
 
typedef struct {
    node cellnode;
#if !defined(QUICKSCAN)
    real rcrit2;
#endif
    nodeptr more;  
    union {
	nodeptr subp[NSUB];
	matrix quad;
    } sorq;
} cell, *cellptr;
 

#define Rcrit2(x)    (((cellptr) (x))->rcrit2)
#define More(x)      (((cellptr) (x))->more)
#define Subp(x)      (((cellptr) (x))->sorq.subp)
#define Quad(x)      (((cellptr) (x))->sorq.quad)


#if !defined(global)
#define global extern
#endif
 
 
#if !defined(QUICKSCAN)
global real theta;
#endif

global string options;
global bool usequad;
global real eps;

 
void maketree(bodyptr, int);

global cellptr root;
global real rsize;
global int ncell;
global int tdepth;
global real cputree;
 

void gravforce(void);
void report_force(stream, int);

global int actmax;
global int nfcalc;
global int nbbcalc;
global int nbccalc;
global real cpuforce;

#endif // ! _sphdefs_h
