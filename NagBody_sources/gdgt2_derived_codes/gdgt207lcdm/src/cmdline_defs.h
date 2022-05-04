
// DEJAR COMOVING INTEGRATION ON; DEJAR PERIODIC; DEJAR MAKEGLASS
//
// DEJAR LOS IF´S DE ComovingIntegrationOn,
//
//
// DEJAR LAS VARIABLES Y CONSTANTES ASOCIADAS A ComovingIntegrationOn AND PeriodicBoundariesOn
// BoxSize, Hubble, Omegas ...
//

/*==============================================================================
	HEADER: cmdline_defs.h			[gdgt207lcdm]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 2007
	Purpose: Definitions for importing arguments from the command line
	Language: C
	Use: '#include "...."
	Use in routines and functions: main
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Major revisions: October 2008
	Copyright: (c) 2005-2014 Mar.  All Rights Reserved
==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

//#define HEAD1	"NagBody"
//#define HEAD2	"Parallel code for the evolution of an N-Body selfgravitating system"
//#define HEAD3	"hierarchical force calculation"

//
// COMMAND LINE DEFINITONS APROPRIATE FOR GASSPHERE CASE ...
//

//string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,

#define HEAD1	"NagBody"
#define HEAD2	"Parallel code for the evolution of an N-Body selfgravitating system"
#define HEAD3	"hierarchical force calculation"

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t ", //HEAD3,
    "ParameterFile=",					";Parameter input file. Overwrite what follows", ":pfile",
	"RestartFlag=0",					";Flag to restart (0, 1, 2)", ":restart",
//
    "SnapFormat=1",						";Output format of snaps (1, 6, 7)", ":ofmt",
    "InitCondFile=./ICs/lcdm_gas_littleendian.dat",		";File with IC", ":ic",
    "OutputDir=lcdm_gas/",						";Simulation working directory",
    "SnapshotFileBase=snapshot",			";Snapshot name pattern",
    "EnergyFile=energy.dat",			";Name file to save simulation statistics",
    "CpuFile=cpu.dat",					";Name file to save cpu statistics",
    "InfoFile=info.dat",				";Name file to save general info",
    "TimingsFile=timings.dat",			";Name file to save timing info",
    "RestartFile=restart",				";Name file to save state run",
//
    "ResubmitCommand=my-scriptfile",				";Resubmit command",
    "OutputListFilename=./parameterfiles/outputs_lcdm_gas.txt",	";List with output times",
    "OutputListOn=1",					";Output list flag",
    "BoxSize=50000.0",						";Simulation box size",
    "PeriodicBoundariesOn=1",			";Flag to use periodic boundary condition",
    "TimeOfFirstSnapshot=0",			";Time to write the first snapshot",
    "CpuTimeBetRestartFile=36000",		";Cpu time to save a state run (here in seconds)",
    "TimeBetStatistics=0.05",			";Time to save statistics info",
    "TimeBegin=0.090909091",					";a to start simulation (z=10)",
    "TimeMax=1.0",						";Time to stop simulation",
    "TimeBetSnapshot=0.5",			";Time to save a snapshot",
    "UnitVelocity_in_cm_per_s=1e5",		";Velocity unit in cm/s (1 km/sec)",
    "UnitLength_in_cm=3.085678e21",				";Length unit in cm (1 kpc)",
    "UnitMass_in_g=1.989e43",				";Mass unit in 1 x 10^10 solar masses",
    "GravityConstantInternal=0",		";Constant G internal value",
    "ErrTolIntAccuracy=0.025",			";Internal tolerance accuracy",
    "ErrTolTheta=0.5",					";Opening parameter value",
    "ErrTolForceAcc=0.005",				";Tolerance for force-acc",
    "MaxSizeTimestep=0.03",			";Maximum time step",
    "MinSizeTimestep=0",				";Minimum time step",
    "MaxRMSDisplacementFac=0.2",				";Maximum RMS displacement factor",
    "DesNumNgb=33",						";Desired number of neighbours",
    "MaxNumNgbDeviation=2",				";Maximum number of neighbours deviation",
    "ICFormat=1",						";Initial condition file format",
    "Omega0=0.3",						";Dark matter content",
    "OmegaBaryon=0.04",					";Baryons content",
    "OmegaLambda=0.7",					";Dark energy content",
    "HubbleParam=0.7",					";Hubble constant present epoch in 100 km/s/Mpc",
//
    "MinGasHsmlFractional=0.25",			";Minimum fraction in gas smoothing",
    "ArtBulkViscConst=0.8",				";Artificial bulk viscosity constant",
    "CourantFac=0.15",						";Courant factor",
    "ComovingIntegrationOn=1",			";Flag for comoving integration",
//
    "NumFilesPerSnapshot=1",			";Number of files per snapshot",
    "NumFilesWrittenInParallel=1",		";Number of files written in parallel",
    "ResubmitOn=0",						";Resubmit flag",
    "TypeOfTimestepCriterion=0",		";Type of time step criterion",
    "TypeOfOpeningCriterion=1",			";Type of cell opening criterion",
    "TimeLimitCPU=36000",				";Cpu time limit (10 hours)",
    "TreeDomainUpdateFrequency=0.1",	";Frequency to update domain",
//
    "SofteningHalo=600.0",				";Softening force length for halo particles",
    "SofteningDisk=0",					";Softening for disk particles",
    "SofteningBulge=0",					";Softening for bulge particles",
    "SofteningGas=600.0",					";Softening for gas particles",
    "SofteningStars=0",					";Softening for star particles",
    "SofteningBndry=0.0",					";Softening for dm particles",
//
    "SofteningHaloMaxPhys=600.0",		";Maximum physical softening for halo particles",
    "SofteningDiskMaxPhys=0",			";Maximum physical softening for disk particles",
    "SofteningBulgeMaxPhys=0",			";Maximum physical softening for bulge particles",
    "SofteningGasMaxPhys=600.0",			";Maximum physical softening for gas particles",
    "SofteningStarsMaxPhys=0",			";Maximum physical softening for star particles",
    "SofteningBndryMaxPhys=0",				";Maximum physical softening for dm particles",
//
    "BufferSize=30",					";Buffer size",
    "PartAllocFactor=1.6",				";Memory allocation factor for particles",
    "TreeAllocFactor=0.8",				";Memory allocation factor for the tree",
//
    "InitGasTemp=1000",					";Initial gas temperature",
    "MinGasTemp=50",						";Minimum gas temperature",
//
    "Version=0.1",						";Mar 2005-2014",
    NULL,
};

#endif /* ! _cmdline_defs_h */
