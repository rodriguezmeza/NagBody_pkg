/*==============================================================================
	HEADER: cmdline_defs.h			[gdgt207]
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
	Copyright: (c) 2005-2008 Mar.  All Rights Reserved
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
    "InitCondFile=./ICs/galaxy_littleendian.dat",		";File with IC", ":ic",
    "OutputDir=galaxy/",				";Simulation working directory",
    "SnapshotFileBase=snapshot",			";Snapshot name pattern",
    "EnergyFile=energy.dat",			";Name file to save simulation statistics",
    "CpuFile=cpu.dat",					";Name file to save cpu statistics",
    "InfoFile=info.dat",				";Name file to save general info",
    "TimingsFile=timings.dat",			";Name file to save timing info",
    "RestartFile=restart",				";Name file to save state run",
//
    "ResubmitCommand=my-scriptfile",				";Resubmit command",
    "OutputListFilename=parameterfiles/output_list.txt",	";List with output times",
    "OutputListOn=0",					";Output list flag",
    "BoxSize=0.0",						";Simulation box size",
    "PeriodicBoundariesOn=0",			";Flag to use periodic boundary condition",
    "TimeOfFirstSnapshot=0.",			";Time to write the first snapshot",
    "CpuTimeBetRestartFile=36000",		";Cpu time to save a state run",
    "TimeBetStatistics=0.05",			";Time to save statistics info",
    "TimeBegin=0.0",					";Time/Redshift to start simulation",
    "TimeMax=3.0",						";Time to stop simulation",
    "TimeBetSnapshot=0.5",			";Time to save a snapshot",
    "UnitVelocity_in_cm_per_s=1.0e5",		";Velocity unit in cm/s [1 km/sec default]",
    "UnitLength_in_cm=3.085678e21",				";Length unit in cm [1 kpc default]",
    "UnitMass_in_g=1.989e43",				";Mass unit in 1 x 10^10 solar masses",
    "GravityConstantInternal=0",		";Constant G internal value",
    "ErrTolIntAccuracy=0.025",			";Internal tolerance accuracy",
    "ErrTolTheta=0.5",					";Opening parameter value",
    "ErrTolForceAcc=0.005",				";Tolerance for force-acc",
    "MaxSizeTimestep=0.01",			";Maximum time step",
    "MinSizeTimestep=0",				";Minimum time step",
    "MaxRMSDisplacementFac=0.2",				";Maximum RMS displacement factor",
    "DesNumNgb=50",						";Desired number of neighbours",
    "MaxNumNgbDeviation=2",				";Maximum number of neighbours deviation",
    "ICFormat=1",						";Initial condition file format",
    "Omega0=0",						";Dark matter content",
    "OmegaBaryon=0",					";Baryons content",
    "OmegaLambda=0",					";Dark energy content",
    "HubbleParam=1.0",					";Hubble constant present epoch in 100 km/s/Mpc",
//
    "MinGasHsmlFractional=0.25",			";Minimum fraction in gas smoothing",
    "ArtBulkViscConst=0.8",				";Artificial bulk viscosity constant",
    "CourantFac=0.15",						";Courant factor",
    "ComovingIntegrationOn=0",			";Flag for comoving integration",
//
    "NumFilesPerSnapshot=1",			";Number of files per snapshot",
    "NumFilesWrittenInParallel=1",		";Number of files written in parallel",
    "ResubmitOn=0",						";Resubmit flag",
    "TypeOfTimestepCriterion=0",		";Type of time step criterion",
    "TypeOfOpeningCriterion=1",			";Type of cell opening criterion",
    "TimeLimitCPU=36000",				";Cpu time limit",
    "TreeDomainUpdateFrequency=0.1",	";Frequency to update domain",
//
    "SofteningHalo=1.0",				";Softening force length for halo particles",
    "SofteningDisk=0.4",					";Softening for disk particles",
    "SofteningBulge=0",					";Softening for bulge particles",
    "SofteningGas=0.0",					";Softening for gas particles",
    "SofteningStars=0",					";Softening for star particles",
    "SofteningBndry=0.0",					";Softening for dm particles",
//
    "SofteningHaloMaxPhys=1.0",		";Maximum physical softening for halo particles",
    "SofteningDiskMaxPhys=0.4",			";Maximum physical softening for disk particles",
    "SofteningBulgeMaxPhys=0",			";Maximum physical softening for bulge particles",
    "SofteningGasMaxPhys=0.0",			";Maximum physical softening for gas particles",
    "SofteningStarsMaxPhys=0",			";Maximum physical softening for star particles",
    "SofteningBndryMaxPhys=0",				";Maximum physical softening for dm particles",
//
    "BufferSize=25",					";Buffer size (in MByte)",
    "PartAllocFactor=1.5",				";Memory allocation factor for particles",
    "TreeAllocFactor=0.8",				";Memory allocation factor for the tree",
//
    "InitGasTemp=0",					";Initial gas temperature",
    "MinGasTemp=0",						";Minimum gas temperature",
//
    "Version=0.1",						";Mar 2005-2014",
    NULL,
};

#endif /* ! _cmdline_defs_h */
