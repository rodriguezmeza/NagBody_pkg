/* ==============================================================================
	HEADER: nagbody_struct.h			[NagBody]
	Written by: Mario A. Rodriguez-Meza
	Starting date: February 2005
	Purpose: Definitions of global variables and parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main), start_run, time_step,
					forcecalc, md_lj_tree_io
	External headers:
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _nagbody_struct_h
#define _nagbody_struct_h

#include <stdio.h>
//#include <gsl/gsl_rng.h>
//#include "tags.h"

#include "globaldefs.h"

// ------------------------START STRUCTURE DEFINITIONS--------------------------

typedef struct _node {
    short type;                 
    bool update;                
    real mass;                  
    vector pos;                 
    struct _node *next;
    vector vel;					// Agregado para manejar num. de cuerpos en nodo   
	int nbodies;				// Agregado para manejar num. de cuerpos en nodo
	real rcut;					// Agregado para metodo de Barnes...
} node, *nodeptr;
 
#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Mass(x)   (((nodeptr) (x))->mass)
#define Pos(x)    (((nodeptr) (x))->pos)
#define Next(x)   (((nodeptr) (x))->next)
#define Vel(x)    (((nodeptr) (x))->vel)		// Agregado para manejar num. 
												// de cuerpos en celda
#define NBodies(x)   (((nodeptr) (x))->nbodies)	// Agregado para manejar num. 
												// de cuerpos en celda
#define Rcut(x)   (((nodeptr) (x))->rcut)

#define BODY 0
#define BODY1 01
#define BODY2 02							// Binary correction
#define CELL 03								// Binary correction
#define SPHBODY 03							// Revisar su uso...
#define DMBODY 04				
#define LNNMAX 1000

#define TESTBODYMU		100
#define TESTBODYSPH		101
#define STATICBODY		102

typedef struct {
    node bodynode;
//    vector vel;			// Comentado para manejar num. de cuerpos en celda
    vector acc;
    real phi;
	real rho;
	real up;
	int Id;
	int IdG;
    vector acc11;
    vector acc12;
    vector acc22;
    real phi11;
    real phi12;
    real phi22;
	
	matrix rf;							// Transport coefficients
	real eng;
} body, *bodyptr;

//#define Vel(x)    (((bodyptr) (x))->vel)	// Comentado para manejar num. de 
											// cuerpos en celda
#define Acc(x)    (((bodyptr) (x))->acc)
#define Phi(x)    (((bodyptr) (x))->phi)
#define Acc11(x)    (((bodyptr) (x))->acc11)
#define Phi11(x)    (((bodyptr) (x))->phi11)
#define Acc12(x)    (((bodyptr) (x))->acc12)
#define Phi12(x)    (((bodyptr) (x))->phi12)
#define Acc22(x)    (((bodyptr) (x))->acc22)
#define Phi22(x)    (((bodyptr) (x))->phi22)
#define Rho(x)    (((bodyptr) (x))->rho)
#define Up(x)    (((bodyptr) (x))->up)
#define Id(x)    (((bodyptr) (x))->Id)
#define IdG(x)    (((bodyptr) (x))->IdG)
#define rf(x)	  (((bodyptr) (x))->rf)			// Transport coefficients
#define en(x)    (((bodyptr) (x))->eng)

// Particle data structure to manipulate I/O
// N > 10^6 purpose...
typedef struct _node_long {
    short type;                 
    bool update;                
    real mass;                  
    vector pos;                 
    struct _node *next;
//    vector vel;					// Agregado para manejar num. de cuerpos en nodo   
//	int nbodies;				// Agregado para manejar num. de cuerpos en nodo
//	real rcut;					// Agregado para metodo de Barnes...
} node_long, *nodeptr_long;

#define Type_long(x)   (((nodeptr_long) (x))->type)
#define Update_long(x) (((nodeptr_long) (x))->update)
#define Mass_long(x)   (((nodeptr_long) (x))->mass)
#define Pos_long(x)    (((nodeptr_long) (x))->pos)
#define Next_long(x)   (((nodeptr_long) (x))->next)
//#define Vel(x)    (((nodeptr) (x))->vel)		// Agregado para manejar num. 
												// de cuerpos en celda
//#define NBodies(x)   (((nodeptr) (x))->nbodies)	// Agregado para manejar num. 
												// de cuerpos en celda
//#define Rcut(x)   (((nodeptr) (x))->rcut)

//#define BODY 0
//#define BODY1 01
//#define BODY2 02							// Binary correction
//#define CELL 03								// Binary correction
//#define SPHBODY 03							// Revisar su uso...
//#define DMBODY 04				
//#define LNNMAX 1000

//#define TESTBODYMU		100
//#define TESTBODYSPH		101
//#define STATICBODY		102

typedef struct {
    node_long bodynode;
    vector vel;
//    vector acc;
//    real phi;
//	real rho;
//	real up;
	int Id;
//	int IdG;
//    vector acc11;
//    vector acc12;
//    vector acc22;
//    real phi11;
//    real phi12;
//    real phi22;
	
//	matrix rf;							// Transport coefficients
//	real eng;
} body_long, *bodyptr_long;

#define Vel_long(x)    (((bodyptr_long) (x))->vel)	// Comentado para manejar num. de 
											// cuerpos en celda
//#define Acc(x)    (((bodyptr) (x))->acc)
//#define Phi(x)    (((bodyptr) (x))->phi)
//#define Acc11(x)    (((bodyptr) (x))->acc11)
//#define Phi11(x)    (((bodyptr) (x))->phi11)
//#define Acc12(x)    (((bodyptr) (x))->acc12)
//#define Phi12(x)    (((bodyptr) (x))->phi12)
//#define Acc22(x)    (((bodyptr) (x))->acc22)
//#define Phi22(x)    (((bodyptr) (x))->phi22)
//#define Rho(x)    (((bodyptr) (x))->rho)
//#define Up(x)    (((bodyptr) (x))->up)
#define Id_long(x)    (((bodyptr_long) (x))->Id)
//#define IdG(x)    (((bodyptr) (x))->IdG)
//#define rf(x)	  (((bodyptr) (x))->rf)			// Transport coefficients
//#define en(x)    (((bodyptr) (x))->eng)

//


#define NSUB (1 << NDIM) 
 
typedef struct {
    node cellnode;       
    real rcrit2;         
    nodeptr more;        
    union {
        nodeptr subp[NSUB]; 
        matrix quad;        
    } sorq;

    vector gpos;						// Geometric Center ...
} cell, *cellptr;
 
#define Rcrit2(x) (((cellptr) (x))->rcrit2)

#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->sorq.subp)
#define Quad(x)   (((cellptr) (x))->sorq.quad)
#define GPos(x)    (((cellptr) (x))->gpos)

//------------------------------------------------------------------------------


#if !defined(global)			// Posicion original de la definicion "global"
#define global extern
#endif

// ------------------------START MACROS DEFINITIONS--------------------------


#define DO_BODY(p,start,finish)  for (p = start; p < finish; p++)
#define DO_DESCENDENTS(q,p)		for (q = More(p); q != Next(p); q = Next(q))
#define DO_COORD(k)				for (k=0; k<NDIM; k++)

// ------------------------START GLOBAL DEFINITIONS--------------------------
//global int nbody;								// Command line... normaly
global bodyptr bodytab;
// Particle data structure to manipulate I/O
// N > 10^6 purpose...
global bodyptr_long bodytab_long;
//


// i/o definitions...
//global string in;			// Command line... normaly DEBE SER LOCAL EN CODIGO
//global string infmt;							// Command line... normaly
//global string out;							// Command line... normaly
//global string outfmt;							// Command line... normaly
//global string basedir;						// Command line... normaly

//global string options;						// Command line... normaly

//global real cpuinit;
//global real tnow;						// DEBE SER LOCAL EN EL CODIGO...

//global FILE *outlog;

// Specific to tree construction

/*
global cellptr root;
global real rsize;
global int ncell;
global int tdepth;
global real cputree;

global real theta;
global bool usequad;
*/

typedef struct {
	cellptr root;
	real rsize;
	int ncell;
	int tdepth;
	real cputree;

	real theta;
	bool usequad;
} global_data_tree;

//global_data_tree gdtree;

typedef struct {
	int actmax;
	int nbccalc;
	int nbbcalc;

	real cpuforce;
	real cpuindforce;
	real cpupot;

	real virSum;
	real virSum11;
	real virSum12;
	real virSum22;
	real uSum;

	real Rcut11Max;
	real Rcut22Max;

	vector Box;

	real RcutSq11;
	real RcutSq12;
	real RcutSq22;

	real ssq11;
	real ssq12;
	real ssq22;

	real fphi11;
	real fphi12;
	real fphi22;

	real fa11;
	real fa12;
	real fa22;

	real vc11;
	real vc12;
	real vc22;
	real dvc11;
	real dvc12;
	real dvc22;

	real Rcut;
	real RcutSq;
	real fphi;
	real ssq;
	real fa;

	vectorI cells;
	int *cellList;

	bool computeTransport;

	int potType;
} global_data_tree_bljforcecalc;

//global_data_tree_bljforcecalc gdbljforce;

typedef struct {
	int actmax;
	int nbccalc;
	int nbbcalc;

	real cpuforce;
	real cpuindforce;
	real cpupot;

	real virSum;
	real virSum11;
	real virSum12;
	real virSum22;
	real uSum;

	real Rcut11Max;
	real Rcut22Max;
	real Rcut33Max;

	vector Box;

	real RcutSq11;
	real RcutSq12;
	real RcutSq13;
	real RcutSq22;
	real RcutSq23;
	real RcutSq33;

	real ssq11;
	real ssq12;
	real ssq13;
	real ssq22;
	real ssq23;
	real ssq33;

	real fphi11;
	real fphi12;
	real fphi13;
	real fphi22;
	real fphi23;
	real fphi33;

	real fa11;
	real fa12;
	real fa13;
	real fa22;
	real fa23;
	real fa33;

	real vc11;
	real vc12;
	real vc13;
	real vc22;
	real vc23;
	real vc33;
	real dvc11;
	real dvc12;
	real dvc13;
	real dvc22;
	real dvc23;
	real dvc33;

	real Rcut;
	real RcutSq;
	real fphi;
	real ssq;
	real fa;

	vectorI cells;
	int *cellList;

	bool computeTransport;

	int potType;
} global_data_tree_tljforcecalc;


typedef struct _pointForcePot {
	int id;
	real r;
	real pot11;
	real pot12;
	real pot22;
	real force11;
	real force12;
	real force22;
} pointForcePot, *pointForcePotptr;

global int nforcepot;
global pointForcePotptr forcepottab;

#define POTTYPE_LJ			0
#define POTTYPE_SLJ			1
#define POTTYPE_FILE		2

#define idPos(x)    (((pointForcePotptr) (x))->id)
#define rPos(x)    (((pointForcePotptr) (x))->r)
#define Pot11(x)    (((pointForcePotptr) (x))->pot11)
#define Pot12(x)    (((pointForcePotptr) (x))->pot12)
#define Pot22(x)    (((pointForcePotptr) (x))->pot22)
#define Force11(x)    (((pointForcePotptr) (x))->force11)
#define Force12(x)    (((pointForcePotptr) (x))->force12)
#define Force22(x)    (((pointForcePotptr) (x))->force22)



// Headers for the snap file formats...

typedef struct {
	int nbody;
	int nbody1;
	int nbody2;
	int ndim;
	real tnow;
	real temperature;
	real density;
	real mass1;
	real mass2;
	real Lx;
	real Ly;
#ifdef THREEDIM
	real Lz;
#endif
	real eps11;
	real eps12;
	real eps22;
	real sigma11;
	real sigma12;
	real sigma22;
	real Rcut11;
	real Rcut12;
	real Rcut22;
} io_header_blj, *io_header_blj_ptr;

typedef struct {
	int nbody;
	int nbody1;
	int nbody2;
	int nbody3;
	int ndim;
	real tnow;
	real temperature;
	real density;
	real mass1;
	real mass2;
	real mass3;
	real Lx;
	real Ly;
#ifdef THREEDIM
	real Lz;
#endif
	real eps11;
	real eps12;
	real eps13;
	real eps22;
	real eps23;
	real eps33;
	real sigma11;
	real sigma12;
	real sigma13;
	real sigma22;
	real sigma23;
	real sigma33;
	real Rcut11;
	real Rcut12;
	real Rcut13;
	real Rcut22;
	real Rcut23;
	real Rcut33;
} io_header_tlj, *io_header_tlj_ptr;


//------------------------------------------------------------------------------


#endif	// ! _nagbody_struct_h
