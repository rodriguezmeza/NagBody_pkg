/* =============================================================================
	HEADER: globaldefs.h				[gbsph]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: Global definitions and structures for command line parameters
	Language: C
	Use: '#include "...."
	Use in routines and functions:
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 24, 2007; October 04, 2007;
	Copyright: (c) 1999-2008 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "switches.h"

#ifndef _globaldefs_h
#define _globaldefs_h

#if !defined(global)
#  define global extern
#endif


#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#ifndef NOGNU
#include "./general_libs/general/stdinc.h"
#include "./general_libs/math/mathfns.h"
#include "./general_libs/math/vectmath.h"
#include "./general_libs/general/getparam.h"
#include "./general_libs/general/strings.h"
#include "./general_libs/io/inout.h"
#include "./general_libs/math/numrec.h"
#include "./general_libs/general/lic.h"
#include "./general_libs/general/machines.h"
#else
#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "strings.h"
#include "inout.h"
#include "numrec.h"
#include "lic.h"
#include "machines.h"
#endif



// COMIENZO DE DEFINICION DE ESTRUCTURAS DE NODOS, CUERPOS Y CELDAS ------------

typedef struct _node {
    short type;                 
    bool update;                
    real mass;                  
    real dm_mass;            
    vector pos;                 
    struct _node *next;         
} node, *nodeptr;
 
#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Mass(x)   (((nodeptr) (x))->mass)
#define dm_Mass(x)   (((nodeptr) (x))->dm_mass)
#define Pos(x)    (((nodeptr) (x))->pos)
#define Next(x)   (((nodeptr) (x))->next)

#define BODY 01                 
#define CELL 02                 
#define SPHBODY 03              
#define DMBODY 04				
#define LNNMAX 1000

typedef struct {
    node bodynode; 
    vector vel;    
    vector acc;    
    real phi;      
	real rho;			
	real up;

	matrix rf;							// Transport coefficients
	real eng;
} body, *bodyptr;

#define Vel(x)    (((bodyptr) (x))->vel)
#define Acc(x)    (((bodyptr) (x))->acc)
#define Phi(x)    (((bodyptr) (x))->phi)
#define Rho(x)    (((bodyptr) (x))->rho)
#define Up(x)    (((bodyptr) (x))->up)

#define rf(x)    (((bodyptr) (x))->rf)		// Transport coefficients
#define en(x)    (((bodyptr) (x))->eng)
  
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
	
	matrix quadq;						// Added for Spline
	real quadp;							// Added for Spline
} cell, *cellptr;
 
#if !defined(QUICKSCAN)
#define Rcrit2(x) (((cellptr) (x))->rcrit2)
#endif

#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->sorq.subp)
#define Quad(x)   (((cellptr) (x))->sorq.quad)
#define QuadQ(x)   (((cellptr) (x))->quadq)
#define QuadP(x)   (((cellptr) (x))->quadp)

// FIN DE DEFINICION DE ESTRUCTURAS DE NODOS, CUERPOS Y CELDAS -----------------


typedef struct {
	string paramfile;
	string forcecalc_method;
	string force_models;

	string icfile;
	string icfilefmt;
	string snapoutfile;
	string snapoutfilefmt;
	string statefile;
	string restorefile;

	string dtimestr;
	string dtoutstr;
	string dtoutinfostr;
	int stepState;
	real tstop;                      

	bool computeTransport;					// Transport properties ...

	string options;                  
 
	bool usequad;                    
#if !defined(QUICKSCAN)
	real theta; 
#endif
	real eps;                        

	int seed;

	real dm_lambda;
	real dm_alpha;
//	real dm_inv_avgphi;
	real G;
	real dm_a;
	real dm_time;

	real eps_pot;
	real sigma_pot;
	real x_pot;
	real y_pot;
	real z_pot;

#if defined(PPNACC)
	real ppn_a1;
#endif

	int nbody;                       

} cmdline_data, *cmdline_data_ptr;


typedef struct {

	int stepSnapInit;
//	int stepSnap;

	real cpuinit;
	real cputotout;
	real cputotal;

	real dtime;
	real dtout;
	real dtoutinfo;

	string headline;
	string headline0;
	string headline1;
	string headline2;
	string headline3;

	string model_comment;			

	real tnow;
	real tout;
	real toutinfo;

	int nstep;                       
	int nstep_grav;					

	vector pos_pot;

	char mode[2];

	cellptr root;                    
	real rsize;                      
	int ncell;                       
	int tdepth;                      
	real cputree;

	int actmax;
	int nbbcalc;
	int nbccalc;
	real cpuforce;

//	int OnlyGrav;
	int IncludeGrav;
	int IncludeGrav_Plummer;
	int IncludeGrav_Spline;
//	int Only2G;
//	int OnlySF_eps;
	int IncludeSF_eps;
//	int OnlySF_noeps;

	bool Plummer;
	
	FILE *outlog;

	int stopflag;

	real dm_alphaTime;

	real virSum;			// Transport ...

} global_data, *global_data_ptr;

global global_data gd;
global cmdline_data cmd;

global bodyptr bodytab;                 


// MACROS - Definitions to add a parameter in the scheme of parameterfile ------

#define IPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&(param);												\
  id[nt++]=INT;}

#define RPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=DOUBLE;}

#define BPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)									\
  {strcpy(tag[nt],paramtext);										\
  param=(string) malloc(n);											\
  addr[nt]=param;													\
  id[nt++]=STRING;}

// -----------------------------------------------------------------------------

#define DO_BODY(p,start,finish)  for (p = start; p < finish; p++)
#define DO_DESCENDENTS(q,p)		for (q = More(p); q != Next(p); q = Next(q))
#define DO_COORD(k)				for (k=0; k<NDIM; k++)

// STATIC problem: gcc version 11
// From inout.h
global real *inout_xval;
global real *inout_yval;
global real *inout_zval;
global real *inout_wval;

// STATIC problem: gcc version 11
// From diffeqs.h
global double dxsav,*xp,**yp;
global int kmax,kount;
global int nrhs;

// STATIC problem: gcc version 11
// From stdinc.h
global long idum;                // seed for random generators

#endif // !_globaldefs_h

