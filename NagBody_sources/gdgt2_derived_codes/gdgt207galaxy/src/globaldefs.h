
//
// QUITAR LAS VARIABLES Y CONSTANTES ASOCIADAS A ComovingIntegrationOn AND PeriodicBoundariesOn
// BoxSize, Hubble, Omegas ...
//


#ifndef ALLVARS_H
#define ALLVARS_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>

//#include "globaldefs.h"
//#include "protodefs.h"



//                      Poner hasta que se decida la estrategia para unificar
//                      gassphere, galaxy, cluster, lcdm_gas
//#include "switches.h"
//

#ifndef NOGNU
#include "./general_libs/general/constant.h"
#include "./general_libs/general/stdinc.h"
#include "./general_libs/general/getparam.h"
#include "./general_libs/math/mathfns.h"
#include "./general_libs/general/lic.h"
#include "./general_libs/io/inout.h"
#include "./general_libs/physics/physcons.h"
#include "./general_libs/physics/eosparam.h"
#include "./general_libs/physics/units.h"
#include "./general_libs/general/machines.h"
#include "./general_libs/mpi/mpi_proto.h"
#else
#include "constant.h"
#include "stdinc.h"
#include "getparam.h"
#include "mathfns.h"
#include "lic.h"
#include "inout.h"
#include "physcons.h"
#include "eosparam.h"
#include "units.h"
#include "machines.h"
#include "mpi_proto.h"
#endif

#include "tags.h"


//
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <time.h>
//#include <mpi.h>

//#include "../../../General_libs/general/constant.h"
//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/mathfns.h"

//#include "../../../General_libs/io/inout.h"

//#include "../../../General_libs/physics/physcons.h"
//#include "../../../General_libs/physics/eosparam.h"
//#include "../../../General_libs/physics/units.h"

//#include "../../../General_libs/general/machines.h"
//



#define  GADGETVERSION   "2.0"

#define  TIMEBASE        (1<<28)

#define  MAXTOPNODES     200000
//#define  MAXTOPNODES     2000000


typedef  long long  peanokey;

#define  BITS_PER_DIMENSION 18

#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))


#define  RNDTABLE         3000
#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  MAXLEN_FILENAME  100

#ifdef   ISOTHERM_EQS
#define  GAMMA         (1.0)
#else
#define  GAMMA         (5.0/3)
#endif

#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.76

#define  GRAVITY           6.672e-8
#define  SOLAR_MASS        1.989e33
#define  SOLAR_LUM         3.826e33
#define  RAD_CONST         7.565e-15
#define  AVOGADRO          6.0222e23
#define  BOLTZMANN         1.3806e-16
#define  GAS_CONST         8.31425e7
#define  C                 2.9979e10
#define  PLANCK            6.6262e-27
#define  CM_PER_MPC        3.085678e24
#define  PROTONMASS        1.6726e-24
#define  ELECTRONMASS      9.10953e-28
#define  THOMPSON          6.65245e-25
#define  ELECTRONCHARGE    4.8032e-10

// SOLO SE NECESITA EN COSMOLOGIA
#define  HUBBLE            3.2407789e-18

#define  SEC_PER_MEGAYEAR  3.155e13
#define  SEC_PER_YEAR      3.155e7

#ifndef ASMTH
#define ASMTH 1.25
#endif

#ifndef RCUT
#define RCUT  4.5
#endif

#define MAX_NGB             20000

#define MAXLEN_OUTPUTLIST   500

#define DRIFT_TABLE_LENGTH  1000

#define MAXITER             150


#ifdef DOUBLEPRECISION
#define FLOAT double
#else
#define FLOAT float
#endif


#ifndef  TWODIMS
#define  NUMDIMS 3
#define  KERNEL_COEFF_1  2.546479089470
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786
#else
#define  NUMDIMS 2
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI
#endif


// ==================== Data structures ========================================

typedef struct {
    
//	char	paramfile[100];
	int		RestartFlag;
//	int		nbody;

//	int		snapoutfmt;
//	char    icfile[100];
    
//	char	outputdir[100],
//    snapout[100],
//    energyfile[100],
//    cpufile[100],
//    infofile[100],
//    timingfile[100],
//    restorefile[100],
//    resubmitcmd[100],
//    outputlistfile[100];
    
//	int		outputlist;
    
//	real	boxsize;
//	int		periodicboundaries;
//	real	timefirstsnap;
//	real	cputimerestore;
//	real	timestatistics;
//	real	timebegin;
//	real	timestop;
//	real	timesnap;
//	real	unitmass,
//    unitvelocity,
//    unitlength;
//	real	gravityconstant;
//	real	maxnodemove;
//	real	treeupdatefreq;
//	real	errtolintegration;
//	real	errtolvelocity;
//	real	errtoltheta;
//	real	errtolforce;
//	real	mintimestep,
//    maxtimestep;
//	int		icfilefmt;
//	int		nfilessnap;
//	int		nfileswparallel;
//	int		resubmit;
//	int		typeopening;
//	int		typetimestep;
//	real	timemaxcpu;
//	real	domainupdatefreq;

//	real	epsbody0;
//	real	epsbody1,
//    epsbody2,
//    epsbody3,
//    epsbody4,
//    epsbody5;
//	int		buffersize;
//	real	pallocfactor;
//	real	treeallocfactor;
    
//	real	dm_alpha,
//    dm_lambda;
//	int		typesoftening;
//	char	options[100];
    
//	int 	bodytypes;
} cmdline_data;

typedef struct
{
  FLOAT s[3];
  FLOAT vs[3];
  FLOAT mass;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int   bitflags;
#else
  FLOAT maxsoft;
#endif
#endif
} DomainNODE;


typedef struct {
  int Daughter;  
  int Pstart;
  int Blocks;
  int Leaf;
  peanokey Size;
  peanokey StartKey;
  long long Count;
} topnode_data;
    

typedef struct {
  long long TotNumPart;
  long long TotN_gas;

  int MaxPart;
  int MaxPartSph;

// SOLO SE NECESITA EN COSMOLOGIA
  double BoxSize;

  int ICFormat;

  int SnapFormat;

  int NumFilesPerSnapshot;
  int NumFilesWrittenInParallel;

  int BufferSize;
  int BunchSizeForce;
  int BunchSizeDensity;
  int BunchSizeHydro;
  int BunchSizeDomain;

  double PartAllocFactor;

  double TreeAllocFactor;

  double DesNumNgb;
  double MaxNumNgbDeviation;

  double ArtBulkViscConst;
  double InitGasTemp;
  double MinGasTemp;
  double MinEgySpec;

  long long TotNumOfForces;
  long long NumForcesSinceLastDomainDecomp;


  double G;
  double UnitTime_in_s;
  double UnitMass_in_g;
  double UnitVelocity_in_cm_per_s;
  double UnitLength_in_cm;
  double UnitPressure_in_cgs;
  double UnitDensity_in_cgs;
  double UnitCoolingRate_in_cgs;
  double UnitEnergy_in_cgs;
  double UnitTime_in_Megayears;
  double GravityConstantInternal;



// SOLO SE NECESITA EN COSMOLOGIA
  double Hubble;
  double Omega0;
  double OmegaLambda;
  double OmegaBaryon;
  double HubbleParam;
  


// SOLO SE NECESITA EN COSMOLOGIA
  int ComovingIntegrationOn;
  int PeriodicBoundariesOn;

  int ResubmitOn;
  int TypeOfOpeningCriterion;
  int TypeOfTimestepCriterion;
  int OutputListOn;



  int SnapshotFileCount;
  double TimeBetSnapshot;
  double TimeOfFirstSnapshot;
  double CpuTimeBetRestartFile;
  double TimeLastRestartFile;
  double TimeBetStatistics;
  double TimeLastStatistics;
  int NumCurrentTiStep;



  double Time;
  double TimeBegin;
  double TimeStep;
  double TimeMax;



  double Timebase_interval;
  int Ti_Current;
  int Ti_nextoutput;
#ifdef FLEXSTEPS
  int PresentMinStep;
  int PresentMaxStep;
#endif
#ifdef PMGRID
  int PM_Ti_endstep;
  int PM_Ti_begstep;
#endif



#ifdef PMGRID
  double Asmth[2];
  double Rcut[2];
  double Corner[2][3];
  double UpperCorner[2][3];
  double Xmintot[2][3];
  double Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif



  double TimeLimitCPU;
  double CPU_TreeConstruction;
  double CPU_TreeWalk;
  double CPU_Gravity;
  double CPU_Potential;
  double CPU_Domain;
  double CPU_Snapshot;
  double CPU_Total;
  double CPU_CommSum;
  double CPU_Imbalance;
  double CPU_HydCompWalk;
  double CPU_HydCommSumm;
  double CPU_HydImbalance;
  double CPU_Hydro;
  double CPU_EnsureNgb;
  double CPU_Predict;
  double CPU_TimeLine;
  double CPU_PM;
  double CPU_Peano;


  double ErrTolTheta;
  double ErrTolForceAcc;


  double ErrTolIntAccuracy;

  double MinSizeTimestep;
  double MaxSizeTimestep;

  double MaxRMSDisplacementFac;

  double CourantFac;


  double TreeDomainUpdateFrequency;


  double MinGasHsmlFractional;
  double MinGasHsml;


  double SofteningGas;
  double SofteningHalo;
  double SofteningDisk;
  double SofteningBulge;
  double SofteningStars;
  double SofteningBndry;

  double SofteningGasMaxPhys;
  double SofteningHaloMaxPhys;
  double SofteningDiskMaxPhys;
  double SofteningBulgeMaxPhys;
  double SofteningStarsMaxPhys;
  double SofteningBndryMaxPhys;

  double SofteningTable[6];
  double ForceSoftening[6];


  double MassTable[6];
  



  char InitCondFile[MAXLEN_FILENAME];
  char OutputDir[MAXLEN_FILENAME];
  char SnapshotFileBase[MAXLEN_FILENAME];
  char EnergyFile[MAXLEN_FILENAME];
  char CpuFile[MAXLEN_FILENAME];
  char InfoFile[MAXLEN_FILENAME];
  char TimingsFile[MAXLEN_FILENAME];
  char RestartFile[MAXLEN_FILENAME];
  char ResubmitCommand[MAXLEN_FILENAME];
  char OutputListFilename[MAXLEN_FILENAME];

  double OutputListTimes[MAXLEN_OUTPUTLIST];
  int OutputListLength;

  string	headline0;
  string	headline1;
  string	headline2;
  string	headline3;

} global_data_all_processes;



typedef struct {
  FLOAT Pos[3];
  FLOAT Mass;
  FLOAT Vel[3];
  FLOAT GravAccel[3];
#ifdef PMGRID
  FLOAT GravPM[3];
#endif
#ifdef FORCETEST
  FLOAT GravAccelDirect[3];
#endif
  FLOAT Potential;
  FLOAT OldAcc;
#ifndef LONGIDS
  unsigned int ID;
#else
  unsigned long long ID;
#endif

  int Type;
  int Ti_endstep;
  int Ti_begstep;
#ifdef FLEXSTEPS
  int FlexStepGrp;
#endif
  float GravCost;
#ifdef PSEUDOSYMMETRIC
  float AphysOld;
#endif
} particle_data;


typedef struct
{
  FLOAT Entropy;
  FLOAT Density;
  FLOAT Hsml;
  FLOAT Left;
  FLOAT Right;
  FLOAT NumNgb;
  FLOAT Pressure;
  FLOAT DtEntropy;
  FLOAT HydroAccel[3];
  FLOAT VelPred[3];
  FLOAT DivVel;
  FLOAT CurlVel;
  FLOAT Rot[3];
  FLOAT DhsmlDensityFactor;
  FLOAT MaxSignalVel;
} sph_particle_data;



typedef struct
{
  FLOAT len;
  FLOAT center[3];
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  FLOAT maxsoft;
#endif
  union
  {
    int suns[8];
    struct
    {
      FLOAT s[3];
      FLOAT mass;
      int bitflags;
      int sibling;
      int nextnode;
      int father;
    }
    d;
  }
  u;
} NODE;



typedef struct
{
  FLOAT hmax;
  FLOAT vs[3];
} extNODE;


typedef struct
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int  flag_entropy_instead_u;
  char fill[60];
} io_header;


typedef struct
{
  double Mass;
  double EnergyKin;
  double EnergyPot;
  double EnergyInt;
  double EnergyTot;
  double Momentum[4];
  double AngMomentum[4];
  double CenterOfMass[4];
  double MassComp[6];
  double EnergyKinComp[6];
  double EnergyPotComp[6];
  double EnergyIntComp[6];
  double EnergyTotComp[6];
  double MomentumComp[6][4]; 
  double AngMomentumComp[6][4]; 
  double CenterOfMassComp[6][4];
} state_of_system;
 


typedef struct
{
  union
  {
    FLOAT Pos[3];
    FLOAT Acc[3];
    FLOAT Potential;
  }
  u;
#ifdef UNEQUALSOFTENINGS
  int Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  FLOAT Soft;
#endif
#endif
  union
  {
    FLOAT OldAcc;
    int Ninteractions;
  }
  w;
} gravdata_in;


typedef struct
{
  int Task;
  int Index;
  int SortIndex;
} gravdata_index;



typedef struct
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  int Index;
  int Task;
} densdata_in;

typedef struct
{
  FLOAT Rho;
  FLOAT Div, Rot[3];
  FLOAT DhsmlDensity;
  FLOAT Ngb;
} densdata_out;


typedef struct
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT Density;
  FLOAT Pressure;
  FLOAT F1;
  FLOAT DhsmlDensityFactor;
  int   Timestep;
  int   Task;
  int   Index;
} hydrodata_in;


typedef struct
{
  FLOAT Acc[3];
  FLOAT DtEntropy;
  FLOAT MaxSignalVel;
} hydrodata_out;


#if !defined(global)
#define global extern
#endif

#ifndef NDIM
#define NDIM	3
#endif
    
#define DO_BODY(p,start,finish)  for (p = start; p < finish; p++)
#define DO_COORD(k)				for (k=0; k<NDIM; k++)
    
#define IPName(param,paramtext)										\
{strcpy(tag[nt],paramtext);                                         \
addr[nt]=&(param);                                                  \
id[nt++]=INT;}
    
#define RPName(param,paramtext)										\
{strcpy(tag[nt],paramtext);                                         \
addr[nt]=&(param);                                                  \
id[nt++]=DOUBLE;}
    
#define BPName(param,paramtext)										\
{strcpy(tag[nt],paramtext);                                         \
addr[nt]=&(param);                                                  \
id[nt++]=BOOLEAN;}
    
#define SPName(param,paramtext)										\
{strcpy(tag[nt],paramtext);                                         \
addr[nt]=param;                                                     \
id[nt++]=STRING;}
    
#define GetCharParam(param,paramtext)								\
{strtmp=GetParam(paramtext);										\
strcpy(param,strtmp);}


global int ThisTask;
global int NTask;
global int PTask;
    
global int NumPart;
global int N_gas;
global long long Ntype[6];
global int NtypeLocal[6];
    
global int NumForceUpdate;
global int NumSphUpdate;
    
global double CPUThisRun;
    

//global int RestartFlag;
    
global char *Exportflag;
    
global int  *Ngblist;
    
global int TreeReconstructFlag;
    
global int Flag_FullStep;
    
    
global gsl_rng *random_generator;
    
global double RndTable[RNDTABLE];
    
    
global double DomainCorner[3];
global double DomainCenter[3];
global double DomainLen;
global double DomainFac;
global int    DomainMyStart;
global int    DomainMyLast;
global int    *DomainStartList;
global int    *DomainEndList;
global double *DomainWork;
global int    *DomainCount;
global int    *DomainCountSph;
    
global int    *DomainTask;
global int    *DomainNodeIndex;
global FLOAT  *DomainTreeNodeLen;
global FLOAT  *DomainHmax;


global peanokey *DomainKeyBuf;
    
global peanokey *Key;
global peanokey *KeySorted;
    
    
global int NTopnodes;
global int NTopleaves;


global DomainNODE *DomainMoment;

global topnode_data *TopNodes;

global double TimeOfLastTreeConstruction;
    
    
    
global char ParameterFile[MAXLEN_FILENAME];
    
global FILE *FdInfo;
global FILE *FdEnergy;
global FILE *FdEnergy2;
global FILE *FdTimings;
global FILE *FdCPU;
    
#ifdef FORCETEST
global FILE *FdForceTest;
#endif
    
    
global double DriftTable[DRIFT_TABLE_LENGTH];
global double GravKickTable[DRIFT_TABLE_LENGTH];
global double HydroKickTable[DRIFT_TABLE_LENGTH];
    
global void *CommBuffer;

//global global_data_all_processes All;
global global_data_all_processes gd;
global cmdline_data cmd;


global particle_data *P,
    *DomainPartBuf;

global sph_particle_data *SphP,
    *DomainSphBuf;

global int MaxNodes;
global int Numnodestree;

global NODE *Nodes_base,
    *Nodes;

global int *Nextnode;
global int *Father;
    
global extNODE *Extnodes_base,
    *Extnodes;

global io_header header;

#define IO_NBLOCKS 11
    
enum iofields
{
        IO_POS,
        IO_VEL,
        IO_ID,
        IO_MASS,
        IO_U,
        IO_RHO,
        IO_HSML,
        IO_POT,
        IO_ACCEL,
        IO_DTENTR,
        IO_TSTP,
};
    
    
global char Tab_IO_Labels[IO_NBLOCKS][4];

global state_of_system SysState;

global gravdata_in *GravDataIn,
    *GravDataGet,
    *GravDataResult,
    *GravDataOut;

global gravdata_index *GravDataIndexTable;

global densdata_in *DensDataIn,
    *DensDataGet;

global densdata_out *DensDataResult,
    *DensDataPartialResult;

global hydrodata_in *HydroDataIn,
    *HydroDataGet;

global hydrodata_out *HydroDataResult,
	*HydroDataPartialResult;


// ===================== globals ===============================================

global  FILE  *stdoutstrm;


// ==================== switches ===============================================

// ... Added to manipulate stdout ...
#define STDOUT	stdoutstrm
#define STDOUTFILE	"stdoutput.log"


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

#endif      // ! ALLVARS_H
