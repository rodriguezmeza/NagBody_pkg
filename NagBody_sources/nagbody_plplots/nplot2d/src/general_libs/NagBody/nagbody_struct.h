/* ==============================================================================
	HEADER: nagbody_struct.h			[NagBody]
	Written by: M.A. Rodriguez-Meza
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
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
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
// gcc11 :: To avoid Error :: duplicate symbol '_bodytab' in:
//global bodyptr bodytab;
local bodyptr bodytab;
// Particle data structure to manipulate I/O
// N > 10^6 purpose...
//global bodyptr_long bodytab_long;
// gcc11 :: To avoid Error :: duplicate symbol '_bodytab_long' in:
local bodyptr_long bodytab_long;
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


typedef struct {
	bool usequad;
	real eps;
	real theta;

	int nbbcalc;
	int nbccalc;
	cellptr root;
	real rsize;
	double cpuforce;
	double cputree;
	int tdepth;
	int ncell;

	char options[200];

	real eps2;
} global_data_treegrav;

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

// gcc11 :: To avoid Error :: duplicate symbol '_nforcepot' in:
//global int nforcepot;
//global pointForcePotptr forcepottab;
local int nforcepot;
local pointForcePotptr forcepottab;

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


// Structures coming from nplot2d ... -------------------------------------------

// Desactivadas por ahora....
/*
real *row;

typedef struct _point {
	real x;
	real y;
} point, *pointptr;

#define Xval(i)   (((pointptr) (i))->x)
#define Yval(i)   (((pointptr) (i))->y)

#define DO_POINT(p,start,finish)  for (p = start; p < finish; p++)

// Block of interface variables with plplot routines
global real *xval;
global real *yval;
global real *zval;
global real *wval;
global int *Typeval;

*/

typedef struct {
	real *xval;
	real *yval;
	real *zval;
	real *wval;
	int *Typeval;
} nplt, *npltptr;

// gcc11 :: To avoid Error :: duplicate symbol '_npltd' in:
//global nplt npltd;
local nplt npltd;

// -----------------------------------------------------------------------------


// gadget global definitions...
#ifdef T3E
  typedef short int int4byte;   // Note: int has 8 Bytes on the T3E !
#else
  typedef int int4byte;
#endif

#define  KERNEL_TABLE 1000

#define  MAX_NGB  20000			// defines maximum length of neighbour list

#define  MAXLEN_OUTPUTLIST 350 // maxmimum number of entries in output list

#define  TIMESTEP_INCREASE_FACTOR 1.3


extern int    NumForceUpdate, NumSphUpdate, IndFirstUpdate;
extern int    TimeTreeRoot;

extern int    Num_nodeupdates, Num_nodeupdate_particles;

// gcc11 :: To avoid Error :: duplicate symbol '_RestartFlag' in:
//global int    RestartFlag;        // Se inicializa en 0, el valor por defecto...
local int    RestartFlag;        // Se inicializa en 0, el valor por defecto...

#define MAX_PARTICLE_TYPES	6		// MAX_PARTICLE_TYPES>=MaxParticlesTypes
// gcc11 :: To avoid Error :: duplicate symbol '_NumPart' in:
/*
global int     NumParticleTypes;        // Number of particles types

global int    NumPart,        // Note: these are the LOCAL process values
              N_gas;
*/

local int     NumParticleTypes;        // Number of particles types
local int    NumPart,        // Note: these are the LOCAL process values
              N_gas;


// variables for input/output ,  usually only used on process 0

extern  char   ParameterFile[100];
extern  FILE  *FdInfo,
              *FdEnergy,
              *FdTimings,
              *FdCPU;


// tabulated smoothing kernel

extern double  Kernel[KERNEL_TABLE+2],
               KernelDer[KERNEL_TABLE+2],
               KernelRad[KERNEL_TABLE+2];


extern double  CPUThisRun;

//#define MAX_PARTICLE_TYPES	6		// MAX_PARTICLE_TYPES>=MaxParticlesTypes
//extern int 	NumParticleTypes;		// Number of particles types

//extern int    NumPart,        // Note: these are the LOCAL process values
//              N_gas;

typedef struct {	// this struct contains data which is the same for all tasks 
					//  (mostly code parameters read from the parameter file)

  int   TotNumPart,         //  particle numbers
        TotN_gas,
        TotN_p0,	// Agregado para tomar en cuenta formato reducido de snaps
        TotN_halo,
        TotN_disk,
        TotN_bulge,
        TotN_stars,
		TotN_dm;			// Dark matter Particle type 6


  int   MaxPart,				// These numbers give the max number of 
//  long int   MaxPart,				// These numbers give the max number of 
								// particles that can be hold
        MaxPartSph;				// on the current processor. (Set it to ~2-3 
								// times the average load)

  int   ICFormat;


  double PartAllocFactor;  // in order to maintain work-load balance,  
//			      the particle load is usually NOT balanced 
//			      each processor allocates memory for PartAllocFactor times  
//			      the average number of particles

  double TreeAllocFactor;  // similarly for the tree:
//			     each processor allocates a number of nodes which is 
//			     TreeAllocFactor times the maximum(!) number of particles 
//			     Note, that a typical local tree for N particles needs typically
//			     1.5-2 N nodes
 
 // some SPH parameters

  int    DesNumNgb;
  double ArtBulkViscConst;
  double InitGasTemp;     // may be used to set the temperature in the IC's  
  double MinGasTemp;      // may be used to set a floor for the gas temperature
  double MinEgySpec;


  // diagnostics

  int    TotNumOfForces;   // counts total number of force computations

  // system of units

  double UnitTime_in_s,
         UnitMass_in_g,
         UnitVelocity_in_cm_per_s,
         UnitLength_in_cm,
         UnitPressure_in_cgs,
         UnitDensity_in_cgs,
         UnitCoolingRate_in_cgs,
         UnitEnergy_in_cgs,
         UnitTime_in_Megayears,
         GravityConstantInternal,
         G;

// Cosmology
  double Hubble;
  double BoxSize, BoxHalf;
  double Omega0,        
         OmegaLambda,
         OmegaBaryon,
         HubbleParam; // little `h', i.e. Hubble const in units of 100 km/s/Mpc. 
					  // Only needed to get absolute physical values 
					  // for cooling physics 

// Code options

  int    ComovingIntegrationOn;   // enables comoving integration
  int    PeriodicBoundariesOn;
  int    ResubmitOn;
  int    TypeOfOpeningCriterion;
  int    TypeOfTimestepCriterion;
  int    OutputListOn;
  int    CoolingOn;


// parameters determining output frequency

  int    SnapshotFileCount;
  double TimeBetSnapshot,
         TimeOfFirstSnapshot,
         CpuTimeBetRestartFile,
         TimeLastRestartFile,
         TimeBetStatistics,
         TimeLastStatistics;

// Current time of the simulation

  int     NumCurrentTiStep;

  double  Time, 
          TimeStep,
          TimeBegin,
          TimeMax;   // marks end of the simulation

  // variables that keep track of cumulative CPU consumption
  double  TimeLimitCPU;
  double  CPU_TreeConstruction;
  double  CPU_TreeWalk;
  double  CPU_Gravity;
  double  CPU_Potential;
  double  CPU_Snapshot;
  double  CPU_Total;
  double  CPU_Hydro;
  double  CPU_Predict;
  double  CPU_TimeLine;


  // tree code opening criterion
  double  ErrTolTheta;
  double  ErrTolForceAcc;


 // adjusts accuracy of time-integration

  double  ErrTolIntAccuracy;  // for 1/a^{1/2} collisionless timestep criterion
  double  ErrTolVelScale;     // for 1/a  collisionless timestep criterion
  double  MinSizeTimestep,
          MaxSizeTimestep;

  double  CourantFac;      // SPH-Courant factor
 

  // frequency of tree reconstruction

  double  MaxNodeMove;
  double  TreeUpdateFrequency;
  int     NumForcesSinceLastTreeConstruction;


// gravitational and hydrodynamical softening lengths 
//   * (given in terms of an `equivalent' Plummer softening length) 
//   *
//   * five groups of particles are supported 
//   * 0=gas,1=halo,2=disk,3=bulge,4=stars 

  double  MinGasHsmlFractional, MinGasHsml;


  double  SofteningGas,
          SofteningHalo,
          SofteningDisk,
          SofteningBulge,
          SofteningStars;

  double  SofteningGasMaxPhys,
          SofteningHaloMaxPhys,
          SofteningDiskMaxPhys,
          SofteningBulgeMaxPhys,
          SofteningStarsMaxPhys;



  double  SofteningTable[6];    
  double  SofteningTableMaxPhys[6];


// If particle masses are all equal for one type, 
//   *  the corresponding entry in MassTable is set to this value,
//   * allowing the size of the snapshot files to be reduced

  double  MassTable[6];   
  
  char    InitCondFile[100],
          OutputDir[100],
          SnapshotFileBase[100],
          EnergyFile[100],
          InfoFile[100],
          TimingsFile[100],
          CpuFile[100],
          RestartFile[100],
          ResubmitCommand[100],
          OutputListFilename[100];

  double  OutputListTimes[MAXLEN_OUTPUTLIST]; // was 200 in earlier version
  int     OutputListLength;



} global_data_all_processes, *global_data_all_processes_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_All' in:
//global global_data_all_processes All;
local global_data_all_processes All;


// The following structure holds all the information that is
// * stored for each particle of the simulation.

typedef struct {
  float     Pos[3];			// particle position at its current time
  float     Vel[3];			// particle velocity at its current time
  float     Mass;			// particle mass
  int4byte  ID;				// unique particle identifier
  int4byte  Type;			// flags particle type. 0=gas, 1=halo, 2=disk, 
							// 3=bulge, 4=stars

  float     CurrentTime;	// current time of the particle
  float     MaxPredTime;	// current time plus half the particles allowed
							// timestep
  float     PosPred[3];		// particle position at the global prediction time   
  float     VelPred[3];		// particle velocity at the global prediction time

  float     Accel[3];		// particle acceleration
  float     Potential;		// particle potential

  float     OldAcc;			// magnitude of old force. Used in new relative 
							// opening criterion

  int4byte  ForceFlag;		// points to next active particle

#ifdef VELDISP  
  float     VelDisp;
  float     HsmlVelDisp;
  float     DensVelDisp;
#endif
} particle_data, *particle_data_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_P_data' in:
//global particle_data_ptr P, P_data;
local particle_data_ptr P, P_data;






typedef struct {
  double     Pos[3];		// particle position at its current time  
  double     Vel[3];		// particle velocity at its current time  
  double     Mass;			// particle mass
  int4byte  ID;				// unique particle identifier
  int4byte  Type;			// flags particle type. 0=gas, 1=halo, 2=disk, 
							// 3=bulge, 4=stars

  double     CurrentTime;	// current time of the particle
  double     MaxPredTime;	// current time plus half the particles 
							// allowed timestep
  double     PosPred[3];	// particle position at the global prediction time
  double     VelPred[3];	// particle velocity at the global prediction time

  double     Accel[3];		// particle acceleration
  double     Potential;		// particle potential

  double     OldAcc;		// magnitude of old force. Used in new relative 
							// opening criterion

  int4byte  ForceFlag;		// points to next active particle

#ifdef VELDISP  
  double     VelDisp;
  double     HsmlVelDisp;
  double     DensVelDisp;
#endif
} particle_data_double, *particle_data_double_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_P_data_double' in:
//global particle_data_ptr P_double, P_data_double;
local particle_data_ptr P_double, P_data_double;


// Particle data structure to manipulate I/O
// N > 10^6 purpose...
typedef struct {
  float     Pos[3];			// particle position at its current time
  float     Vel[3];			// particle velocity at its current time
  float     Mass;			// particle mass
  int4byte  ID;				// unique particle identifier
  int4byte  Type;			// flags particle type. 0=gas, 1=halo, 2=disk, 
							// 3=bulge, 4=stars

//  float     CurrentTime;	// current time of the particle
//  float     MaxPredTime;	// current time plus half the particles allowed
							// timestep
//  float     PosPred[3];		// particle position at the global prediction time   
//  float     VelPred[3];		// particle velocity at the global prediction time

//  float     Accel[3];		// particle acceleration
//  float     Potential;		// particle potential

//  float     OldAcc;			// magnitude of old force. Used in new relative 
							// opening criterion

//  int4byte  ForceFlag;		// points to next active particle

//#ifdef VELDISP  
//  float     VelDisp;
//  float     HsmlVelDisp;
//  float     DensVelDisp;
//#endif
} particle_data_long, *particle_data_long_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_P_data_long' in:
//global particle_data_long_ptr P_long, P_data_long;
local particle_data_long_ptr P_long, P_data_long;
//

typedef struct {
  float  Density;         // particle density at its current time
  float  DtDensity;       // rate of change of density
  float  DensityPred;     // predicted particle density

  float  EgySpec;         // internal energy per unit mass
  float  DtEgySpec;       // rate of change of the internal energy
  float  EgySpecPred;     // predicted internal energy per unit mass

  float  Pressure;        // pressure

  float  Hsml;            // smoothing length
  float  DtHsml;          // rate of change of smoothing length

  int    NumNgb;          // number of SPH neighbours

  float  DivVel;          // local velocity divergence
  float  CurlVel;         // local velocity curl

#ifdef COOLING
  float  Ne;              // electron fraction. Gives indirectly ionization 
                          // state and mean molecular weight.
#endif
} sph_particle_data, *sph_particle_data_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_SphP_data' in:
//global sph_particle_data_ptr SphP, SphP_data;
local sph_particle_data_ptr SphP, SphP_data;

typedef struct {
  double  Density;         // particle density at its current time
  double  DtDensity;       // rate of change of density
  double  DensityPred;     // predicted particle density

  double  EgySpec;         // internal energy per unit mass
  double  DtEgySpec;       // rate of change of the internal energy
  double  EgySpecPred;     // predicted internal energy per unit mass

  double  Pressure;        // pressure

  double  Hsml;            // smoothing length
  double  DtHsml;          // rate of change of smoothing length

  int    NumNgb;          // number of SPH neighbours

  double  DivVel;          // local velocity divergence
  double  CurlVel;         // local velocity curl

#ifdef COOLING
  double  Ne;              // electron fraction. Gives indirectly ionization 
                           // state and mean molecular weight.
#endif
} sph_particle_data_double, *sph_particle_data_double_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_SphP_data_double' in:
//global sph_particle_data_double_ptr SphP_double, SphP_data_double;
local sph_particle_data_double_ptr SphP_double, SphP_data_double;

// this structure holds nodes for the ordered binary tree of the timeline.
typedef struct {
  int4byte left,right;
} timetree_data, *timetree_data_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_PTimeTree' in:
//global timetree_data_ptr PTimeTree;        // for ordered binary tree of max pred.
local timetree_data_ptr PTimeTree;        // for ordered binary tree of max pred.
										// times

// state of total system
typedef struct {
  double  Mass,
          EnergyKin,
          EnergyPot,
          EnergyInt,
          EnergyTot,
          Momentum[4],
          AngMomentum[4],
          CenterOfMass[4],

          MassComp[5],
          EnergyKinComp[5],
          EnergyPotComp[5],
          EnergyIntComp[5],
          EnergyTotComp[5],
          MomentumComp[5][4],
          AngMomentumComp[5][4],
          CenterOfMassComp[5][4];
} state_of_system, *state_of_system_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_SysStateAtStart' in:
//global  state_of_system SysState,SysStateAtStart,SysStateAtEnd;
local  state_of_system SysState,SysStateAtStart,SysStateAtEnd;


// Headers for the snap file formats...

typedef struct {
  int4byte npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int4byte flag_sfr;
  int4byte flag_feedback;
  int4byte npartTotal[6];
  int4byte flag_cooling;
  int4byte num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  // fills to 256 Bytes
} io_header_1, *io_header_1_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_header1' in:
//global io_header_1 header1;
local io_header_1 header1;

typedef struct {
  int4byte npart[6];
  double   mass[6];
  double   time;
  int4byte npartTotal[6];
  int4byte num_files;
  double   BoxSize;
  char     fill[256- 6*4- 6*8- 1*8- 0*4- 6*4- 1*4 - 1*8];
                                    // Completa 256 Bytes
} io_header_reducido, *io_header_reducido_ptr;


// gcc11 :: To avoid Error :: duplicate symbol '_header_reducido' in:
//global io_header_reducido header_reducido;
local io_header_reducido header_reducido;

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


typedef struct {
	int nbody;
	int ndim;
	real tnow;
} io_header_standard;


// Comienzo de TIPSY STRUCTURE -------------------------------------------------

#define MAXDIM 3
#define forever for(;;)

typedef float Real;

typedef struct  {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} gas_particle, *gas_particle_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_gas_particles' in:
//global gas_particle_ptr gas_particles;
local gas_particle_ptr gas_particles;

typedef struct  {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} dark_particle, *dark_particle_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_dark_particles' in:
//global dark_particle_ptr dark_particles;
local dark_particle_ptr dark_particles;

typedef struct {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} star_particle, *star_particle_ptr;

// gcc11 :: To avoid Error :: duplicate symbol '_star_particles' in:
//global star_particle_ptr star_particles;
local star_particle_ptr star_particles;

typedef struct {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} dump, *dumptr;

// gcc11 :: To avoid Error :: duplicate symbol '_header' in:
//global dump header ;
local dump header ;

// FIN de TIPSY STRUCTURE ------------------------------------------------------

//------------------------------------------------------------------------------


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////// MODIFICACIONES PARA LEER Y ESCRIBIR EN FORMATO GADGET207
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// GADGET207 :: COMIENZO //////////////////////////////////////////////////////////
// ESTRUCTURAS DE DATOS PARA INTERACTUAR CON EL FORMATO GADGET 207
#ifdef DOUBLEPRECISION   /*!< If defined, the variable type FLOAT is set to "double", otherwise to FLOAT */
#define FLOAT double
#else
#define FLOAT float
#endif

#define  MAXLEN_FILENAME  100    /*!< Maximum number of characters for filenames (including the full path) */


typedef struct
{
    long long TotNumPart;		/*!< total particle numbers (global value) */
    long long TotN_gas;		/*!< total gas particle number (global value) */
    
    int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one processor. */
    int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one processor. */
    
    double BoxSize;               /*!< Boxsize in case periodic boundary conditions are used */
    
    int ICFormat;			/*!< selects different versions of IC file-format */
    
    int SnapFormat;		/*!< selects different versions of snapshot file-formats */
    
    int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
    int NumFilesWrittenInParallel;/*!< maximum number of files that may be written simultaneously when
                                   writing/reading restart-files, or when writing snapshot files */ 
    
    int BufferSize;		/*!< size of communication buffer in MB */
    int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force algorithm  */
    int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density computation */
    int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH hydrodynamical force computation */
    int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain decomposition */
    
    double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
                             NOT be balanced.  Each processor allocates memory for PartAllocFactor times
                             the average number of particles to allow for that */
    
    double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
                             the maximum(!) number of particles.  Note: A typical local tree for N
                             particles needs usually about ~0.65*N nodes. */
    
    /* some SPH parameters */
    
    double DesNumNgb;             /*!< Desired number of SPH neighbours */
    double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */
    
    double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
    double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
    double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
    double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */
    
    
    /* some force counters  */
    
    long long TotNumOfForces;	             /*!< counts total number of force computations  */
    long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */
    
    
    /* system of units  */
    
    double G;                        /*!< Gravity-constant in internal units */
    double UnitTime_in_s;   	   /*!< factor to convert internal time unit to seconds/h */
    double UnitMass_in_g;            /*!< factor to convert internal mass unit to grams/h */
    double UnitVelocity_in_cm_per_s; /*!< factor to convert intqernal velocity unit to cm/sec */
    double UnitLength_in_cm;         /*!< factor to convert internal length unit to cm/h */
    double UnitPressure_in_cgs;      /*!< factor to convert internal pressure unit to cgs units (little 'h' still around!) */
    double UnitDensity_in_cgs;       /*!< factor to convert internal length unit to g/cm^3*h^2 */
    double UnitCoolingRate_in_cgs;   /*!< factor to convert internal cooling rate to cgs units */
    double UnitEnergy_in_cgs;        /*!< factor to convert internal energy to cgs units */
    double UnitTime_in_Megayears;    /*!< factor to convert internal time to megayears/h */
    double GravityConstantInternal;  /*!< If set to zero in the parameterfile, the internal value of the
                                      gravitational constant is set to the Newtonian value based on the system of
                                      units specified. Otherwise the value provided is taken as internal gravity constant G. */
    
    
    /* Cosmological parameters */
    
    double Hubble;       /*!< Hubble-constant in internal units */
    double Omega0;       /*!< matter density in units of the critical density (at z=0)*/
    double OmegaLambda;  /*!< vaccum energy density relative to crictical density (at z=0) */
    double OmegaBaryon;  /*!< baryon density in units of the critical density (at z=0)*/
    double HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute physical values for cooling physics */
    
    
    /* Code options */
    
    int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
    int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
    int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
    int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative criterion */
    int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike gadget-1.1) */
    int OutputListOn;             /*!< flags that output times are listed in a specified file */
    
    
    /* Parameters determining output frequency */
    
    int SnapshotFileCount;        /*!< number of snapshot that is written next */
    double TimeBetSnapshot;       /*!< simulation time interval between snapshot files */
    double TimeOfFirstSnapshot;   /*!< simulation time of first snapshot files */
    double CpuTimeBetRestartFile; /*!< cpu-time between regularly generated restart files */
    double TimeLastRestartFile;   /*!< cpu-time when last restart-file was written */
    double TimeBetStatistics;     /*!< simulation time interval between computations of energy statistics */
    double TimeLastStatistics;    /*!< simulation time when the energy statistics was computed the last time */
    int NumCurrentTiStep;         /*!< counts the number of system steps taken up to this point */
    
    
    /* Current time of the simulation, global step, and end of simulation */
    
    double Time;                  /*!< current time of the simulation */
    double TimeBegin;             /*!< time of initial conditions of the simulation */
    double TimeStep;              /*!< difference between current times of previous and current timestep */
    double TimeMax;	        /*!< marks the point of time until the simulation is to be evolved */
    
    
    /* variables for organizing discrete timeline */
    
    double Timebase_interval;     /*!< factor to convert from floating point time interval to integer timeline */
    int Ti_Current;               /*!< current time on integer timeline */ 
    int Ti_nextoutput;            /*!< next output time on integer timeline */
#ifdef FLEXSTEPS
    int PresentMinStep;           /*!< If FLEXSTEPS is used, particle timesteps are chosen as multiples of the present minimum timestep. */
    int PresentMaxStep;		/*!< If FLEXSTEPS is used, this is the maximum timestep in timeline units, rounded down to the next power 2 division */
#endif
#ifdef PMGRID
    int PM_Ti_endstep;            /*!< begin of present long-range timestep */
    int PM_Ti_begstep;            /*!< end of present long-range timestep */
#endif
    
    
    /* Placement of PM grids */
    
#ifdef PMGRID
    double Asmth[2];              /*!< Gives the scale of the long-range/short-range split (in mesh-cells), both for the coarse and the high-res mesh */
    double Rcut[2];               /*!< Gives the maximum radius for which the short-range force is evaluated with the tree (in mesh-cells), both for the coarse and the high-res mesh */
    double Corner[2][3];          /*!< lower left corner of coarse and high-res PM-mesh */
    double UpperCorner[2][3];     /*!< upper right corner of coarse and high-res PM-mesh */
    double Xmintot[2][3];         /*!< minimum particle coordinates both for coarse and high-res PM-mesh */
    double Xmaxtot[2][3];         /*!< maximum particle coordinates both for coarse and high-res PM-mesh */
    double TotalMeshSize[2];      /*!< total extension of coarse and high-res PM-mesh */
#endif
    
    
    /* Variables that keep track of cumulative CPU consumption */
    
    double TimeLimitCPU;          /*!< CPU time limit as defined in parameterfile */
    double CPU_TreeConstruction;  /*!< time spent for constructing the gravitational tree */
    double CPU_TreeWalk;          /*!< actual time spent for pure tree-walks */
    double CPU_Gravity;           /*!< cumulative time used for gravity computation (tree-algorithm only) */
    double CPU_Potential;         /*!< time used for computing gravitational potentials */
    double CPU_Domain;            /*!< cumulative time spent for domain decomposition */
    double CPU_Snapshot;          /*!< time used for writing snapshot files */
    double CPU_Total;             /*!< cumulative time spent for domain decomposition */
    double CPU_CommSum;           /*!< accumulated time used for communication, and for collecting partial results, in tree-gravity */
    double CPU_Imbalance;         /*!< cumulative time lost accross all processors as work-load imbalance in gravitational tree */
    double CPU_HydCompWalk;       /*!< time used for actual SPH computations, including neighbour search */
    double CPU_HydCommSumm;       /*!< cumulative time used for communication in SPH, and for collecting partial results */
    double CPU_HydImbalance;      /*!< cumulative time lost due to work-load imbalance in SPH */
    double CPU_Hydro;             /*!< cumulative time spent for SPH related computations */
    double CPU_EnsureNgb;         /*!< time needed to iterate on correct neighbour numbers */
    double CPU_Predict;           /*!< cumulative time to drift the system forward in time, including dynamic tree updates */
    double CPU_TimeLine;          /*!< time used for determining new timesteps, and for organizing the timestepping, including kicks of active particles */
    double CPU_PM;                /*!< time used for long-range gravitational force */
    double CPU_Peano;             /*!< time required to establish Peano-Hilbert order */
    
    /* tree code opening criterion */
    
    double ErrTolTheta;		/*!< BH tree opening angle */
    double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */
    
    
    /* adjusts accuracy of time-integration */
    
    double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                 timestep is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */
    
    double MinSizeTimestep;       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                   timestep determined by the timestep criteria falls below this limit. */ 
    double MaxSizeTimestep;       /*!< maximum allowed timestep */
    
    double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
                                   in comoving coordinates.  To this end, the code computes the rms velocity
                                   of all particles, and limits the timestep such that the rms displacement
                                   is a fraction of the mean particle separation (determined from the
                                   particle mass and the cosmological parameters). This parameter specifies
                                   this fraction. */
    
    double CourantFac;		/*!< SPH-Courant factor */
    
    
    /* frequency of tree reconstruction/domain decomposition */
    
    double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */
    
    
    /* Gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening length).
     * Five groups of particles are supported 0="gas", 1="halo", 2="disk", 3="bulge", 4="stars", 5="bndry"
     */
    
    double MinGasHsmlFractional;  /*!< minimum allowed SPH smoothing length in units of SPH gravitational softening length */
    double MinGasHsml;            /*!< minimum allowed SPH smoothing length */
    
    
    double SofteningGas;          /*!< comoving gravitational softening lengths for type 0 */ 
    double SofteningHalo;         /*!< comoving gravitational softening lengths for type 1 */ 
    double SofteningDisk;         /*!< comoving gravitational softening lengths for type 2 */ 
    double SofteningBulge;        /*!< comoving gravitational softening lengths for type 3 */ 
    double SofteningStars;        /*!< comoving gravitational softening lengths for type 4 */ 
    double SofteningBndry;        /*!< comoving gravitational softening lengths for type 5 */ 
    
    double SofteningGasMaxPhys;   /*!< maximum physical softening length for type 0 */ 
    double SofteningHaloMaxPhys;  /*!< maximum physical softening length for type 1 */ 
    double SofteningDiskMaxPhys;  /*!< maximum physical softening length for type 2 */ 
    double SofteningBulgeMaxPhys; /*!< maximum physical softening length for type 3 */ 
    double SofteningStarsMaxPhys; /*!< maximum physical softening length for type 4 */ 
    double SofteningBndryMaxPhys; /*!< maximum physical softening length for type 5 */ 
    
    double SofteningTable[6];     /*!< current (comoving) gravitational softening lengths for each particle type */
    double ForceSoftening[6];     /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */
    
    
    double MassTable[6];          /*!< Table with particle masses for particle types with equal mass.
                                   If particle masses are all equal for one type, the corresponding entry in MassTable 
                                   is set to this value, allowing the size of the snapshot files to be reduced. */
    
    
    
    /* some filenames */
    
    char InitCondFile[MAXLEN_FILENAME];          /*!< filename of initial conditions */
    char OutputDir[MAXLEN_FILENAME];             /*!< output directory of the code */
    char SnapshotFileBase[MAXLEN_FILENAME];      /*!< basename to construct the names of snapshotf files */
    char EnergyFile[MAXLEN_FILENAME];            /*!< name of file with energy statistics */
    char CpuFile[MAXLEN_FILENAME];               /*!< name of file with cpu-time statistics */
    char InfoFile[MAXLEN_FILENAME];              /*!< name of log-file with a list of the timesteps taken */
    char TimingsFile[MAXLEN_FILENAME];           /*!< name of file with performance metrics of gravitational tree algorithm */
    char RestartFile[MAXLEN_FILENAME];           /*!< basename of restart-files */
    char ResubmitCommand[MAXLEN_FILENAME];       /*!< name of script-file that will be executed for automatic restart */
    char OutputListFilename[MAXLEN_FILENAME];    /*!< name of file with list of desired output times */
    
    double OutputListTimes[MAXLEN_OUTPUTLIST];   /*!< table with desired output times */
    int OutputListLength;                        /*!< number of output times stored in the table of desired output times */
    
} global_data_all_processes_GADGET207;

// gcc11 :: To avoid Error :: duplicate symbol '_All_GADGET' in:
//global global_data_all_processes_GADGET207 All_GADGET;
local global_data_all_processes_GADGET207 All_GADGET;

// BORRAR DESPUES DE REVISAR USO DE LA ESTRUCTURA ANTERIOR ...
//All;                                          /*!< a container variable for global variables that are equal on all processors */


typedef struct
{
    FLOAT Pos[3];			/*!< particle position at its current time */
    FLOAT Mass;			/*!< particle mass */
    FLOAT Vel[3];			/*!< particle velocity at its current time */
    FLOAT GravAccel[3];		/*!< particle acceleration due to gravity */
#ifdef PMGRID
    FLOAT GravPM[3];		/*!< particle acceleration due to long-range PM gravity force*/
#endif
#ifdef FORCETEST
    FLOAT GravAccelDirect[3];	/*!< particle acceleration when computed with direct summation */
#endif
    FLOAT Potential;		/*!< gravitational potential */
    FLOAT OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening criterion */
#ifndef LONGIDS
    unsigned int ID;		/*!< particle identifier */
#else
    unsigned long long ID;        /*!< particle identifier */
#endif
    
    int Type;		        /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
    int Ti_endstep;               /*!< marks start of current timestep of particle on integer timeline */ 
    int Ti_begstep;               /*!< marks end of current timestep of particle on integer timeline */
#ifdef FLEXSTEPS
    int FlexStepGrp;		/*!< a random 'offset' on the timeline to create a smooth groouping of particles */
#endif
    float GravCost;		/*!< weight factor used for balancing the work-load */
#ifdef PSEUDOSYMMETRIC
    float AphysOld;               /*!< magnitude of acceleration in last timestep. Used to make a first order
                                   prediction of the change of acceleration expected in the future, thereby
                                   allowing to guess whether a decrease/increase of the timestep should occur
                                   in the timestep that is started. */
#endif
} particle_data_GADGET207, *particle_data_ptr_GADGET207;

// gcc11 :: To avoid Error :: duplicate symbol '_P_GADGET' in:
//global particle_data_ptr_GADGET207 P_GADGET;
local particle_data_ptr_GADGET207 P_GADGET;

// BORRAR DESPUES DE REVISAR USO DE LA ESTRUCTURA ANTERIOR ...
//*P,              /*!< holds particle data on local processor */
//*DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */

typedef struct 
{
    int npart[6];                        /*!< number of particles of each type in this file */
    double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                          stored in the mass-block of the snapshot file, otherwise they are omitted */
    double time;                         /*!< time of snapshot file */
    double redshift;                     /*!< redshift of snapshot file */
    int flag_sfr;                        /*!< flags whether the simulation was including star formation */
    int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
    unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
                                          different from npart if one is dealing with a multi-file snapshot. */
    int flag_cooling;                    /*!< flags whether cooling was included  */
    int num_files;                       /*!< number of files in multi-file snapshot */
    double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
    double Omega0;                       /*!< matter density in units of critical density */
    double OmegaLambda;                  /*!< cosmological constant parameter */
    double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
    int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
    int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
    unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
    int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
    char fill[60];	               /*!< fills to 256 Bytes */
} io_header_GADGET207, *io_header_ptr_GADGET207;

// gcc11 :: To avoid Error :: duplicate symbol '_header_GADGET' in:
//global io_header_ptr_GADGET207 header_GADGET;
local io_header_ptr_GADGET207 header_GADGET;

// BORRAR DESPUES DE REVISAR USO DE LA ESTRUCTURA ANTERIOR ...
//header;                               /*!< holds header for snapshot files */

// GADGET207 FIN ///////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/////// FIN FIN FIN MODIFICACIONES PARA LEER Y ESCRIBIR EN FORMATO GADGET207
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/////// MODIFICACIONES PARA LEER Y ESCRIBIR EN FORMATO GADGET207 version IBERO
////////////////////////////////////////////////////////////////////////////////
// IBERO :: COMIENZO //////////////////////////////////////////////////////////
typedef struct {
    float     Pos[3];			// particle position at its current time
    float     Vel[3];			// particle velocity at its current time
    float     Mass;			// particle mass
    int  ID;				// unique particle identifier
    int  Type;			// flags particle type. 0=gas, 1=halo, 2=disk, 
    // 3=bulge, 4=stars
    
    float     CurrentTime;	// current time of the particle
    float     MaxPredTime;	// current time plus half the particles allowed
    // timestep
    float     PosPred[3];		// particle position at the global prediction time   
    float     VelPred[3];		// particle velocity at the global prediction time
    
    float     Accel[3];		// particle acceleration
    float     Potential;		// particle potential
    
    float     OldAcc;			// magnitude of old force. Used in new relative 
    // opening criterion
    
    int  ForceFlag;		// points to next active particle
    
#ifdef VELDISP  
    float     VelDisp;
    float     HsmlVelDisp;
    float     DensVelDisp;
#endif
} particle_data_IBERO, *particle_data_ptr_IBERO;

// gcc11 :: To avoid Error :: duplicate symbol '_P_IBERO_data' in:
//global particle_data_ptr_IBERO P_IBERO, P_IBERO_data;
local particle_data_ptr_IBERO P_IBERO, P_IBERO_data;


typedef struct {
    int4byte npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int4byte flag_sfr;
    int4byte flag_feedback;
    int4byte npartTotal[6];
    int4byte flag_cooling;
    int4byte num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  // fills to 256 Bytes
} io_header_1_IBERO, *io_header_1_ptr_IBERO;

// gcc11 :: To avoid Error :: duplicate symbol '_header1_IBERO' in:
//global io_header_1_IBERO header1_IBERO;
local io_header_1_IBERO header1_IBERO;

// IBERO FIN ///////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/////// FIN FIN FIN MODIFICACIONES PARA LEER Y ESCRIBIR EN FORMATO GADGET207 version IBERO
////////////////////////////////////////////////////////////////////////////////

#endif	// ! _nagbody_struct_h
