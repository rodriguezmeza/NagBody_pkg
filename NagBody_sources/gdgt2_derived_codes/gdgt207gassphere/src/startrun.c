
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS




#include "globaldefs.h"
#include "protodefs.h"

local int startrun_parameterfile(char *);
local int startrun_cmdline(void);
local void PrintParameterFile(char *);

//void begrun(void)
void startrun(string head0, string head1, string head2, string head3)
{
	global_data_all_processes all;
	cmdline_data cmd_old;

  if(ThisTask == 0)
    {
      printf("\nThis is Gadget, version `%s'.\n", GADGETVERSION);
      printf("\nRunning on %d processors.\n", NTask);
    }

  read_parameter_file(ParameterFile);

// QUITAR ESTA LINEA CUANDO YA SE HALLA EXPULGADO EL CODIGO (RESTARTFLAG)
// HAY OTRAS LINEAS SIMILARES EN OTRAS PARTES (BUSCAR RESTARTFLAG)
//fprintf_mpi_flush(STDOUT, ThisTask, "\nbegrun:: RestartFlag %d\n", cmd.RestartFlag);

  allocate_commbuffers();
  set_units();

/* #if defined(PERIODIC) && (!defined(PMGRID) || defined(FORCETEST))
  ewald_init();
#endif */

  open_outputfiles();

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);

#ifdef PMGRID
  long_range_init();
#endif

  gd.TimeLastRestartFile = CPUThisRun;

//  if(RestartFlag == 0 || RestartFlag == 2)
  if(cmd.RestartFlag==0 || cmd.RestartFlag==2)
    {
      set_random_numbers();

      init();
    }
  else
    {
      all = gd;
      cmd_old	=	cmd;

//        restart(RestartFlag);
        restart(cmd.RestartFlag);

      gd.MinSizeTimestep = gd.MinSizeTimestep;
      gd.MaxSizeTimestep = gd.MaxSizeTimestep;
      gd.BufferSize = gd.BufferSize;
      gd.BunchSizeForce = gd.BunchSizeForce;
      gd.BunchSizeDensity = gd.BunchSizeDensity;
      gd.BunchSizeHydro = gd.BunchSizeHydro;
      gd.BunchSizeDomain = gd.BunchSizeDomain;

      gd.TimeLimitCPU = gd.TimeLimitCPU;
      gd.ResubmitOn = gd.ResubmitOn;
      gd.TimeBetSnapshot = gd.TimeBetSnapshot;
      gd.TimeBetStatistics = gd.TimeBetStatistics;
      gd.CpuTimeBetRestartFile = gd.CpuTimeBetRestartFile;
      gd.ErrTolIntAccuracy = gd.ErrTolIntAccuracy;
      gd.MaxRMSDisplacementFac = gd.MaxRMSDisplacementFac;

      gd.ErrTolForceAcc = gd.ErrTolForceAcc;

      gd.TypeOfTimestepCriterion = gd.TypeOfTimestepCriterion;
      gd.TypeOfOpeningCriterion = gd.TypeOfOpeningCriterion;
      gd.NumFilesWrittenInParallel = gd.NumFilesWrittenInParallel;
      gd.TreeDomainUpdateFrequency = gd.TreeDomainUpdateFrequency;

      gd.SnapFormat = gd.SnapFormat;
      gd.NumFilesPerSnapshot = gd.NumFilesPerSnapshot;
      gd.MaxNumNgbDeviation = gd.MaxNumNgbDeviation;
      gd.ArtBulkViscConst = gd.ArtBulkViscConst;


      gd.OutputListOn = gd.OutputListOn;
      gd.CourantFac = gd.CourantFac;

      gd.OutputListLength = gd.OutputListLength;
      memcpy(gd.OutputListTimes, gd.OutputListTimes, sizeof(double) * gd.OutputListLength);


      strcpy(gd.ResubmitCommand, gd.ResubmitCommand);
      strcpy(gd.OutputListFilename, gd.OutputListFilename);
      strcpy(gd.OutputDir, gd.OutputDir);
      strcpy(gd.RestartFile, gd.RestartFile);
      strcpy(gd.EnergyFile, gd.EnergyFile);
      strcpy(gd.InfoFile, gd.InfoFile);
      strcpy(gd.CpuFile, gd.CpuFile);
      strcpy(gd.TimingsFile, gd.TimingsFile);
      strcpy(gd.SnapshotFileBase, gd.SnapshotFileBase);

      if(gd.TimeMax != gd.TimeMax)
	readjust_timebase(gd.TimeMax, gd.TimeMax);
    }

#ifdef PMGRID
  long_range_init_regionsize();
#endif

  if(gd.ComovingIntegrationOn)
    init_drift_table();

// if(RestartFlag == 2)
  if(cmd.RestartFlag == 2)
    gd.Ti_nextoutput = find_next_outputtime(gd.Ti_Current + 1);
  else
    gd.Ti_nextoutput = find_next_outputtime(gd.Ti_Current);


  gd.TimeLastRestartFile = CPUThisRun;
}



void set_units(void)
{
  double meanweight;

  gd.UnitTime_in_s = gd.UnitLength_in_cm / gd.UnitVelocity_in_cm_per_s;
  gd.UnitTime_in_Megayears = gd.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(gd.GravityConstantInternal == 0)
    gd.G = GRAVITY / pow(gd.UnitLength_in_cm, 3) * gd.UnitMass_in_g * pow(gd.UnitTime_in_s, 2);
  else
    gd.G = gd.GravityConstantInternal;

  gd.UnitDensity_in_cgs = gd.UnitMass_in_g / pow(gd.UnitLength_in_cm, 3);
  gd.UnitPressure_in_cgs = gd.UnitMass_in_g / gd.UnitLength_in_cm / pow(gd.UnitTime_in_s, 2);
  gd.UnitCoolingRate_in_cgs = gd.UnitPressure_in_cgs / gd.UnitTime_in_s;
  gd.UnitEnergy_in_cgs = gd.UnitMass_in_g * pow(gd.UnitLength_in_cm, 2) / pow(gd.UnitTime_in_s, 2);


  gd.Hubble = HUBBLE * gd.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", gd.Hubble);
      printf("G (internal units) = %g\n", gd.G);
      printf("UnitMass_in_g = %g \n", gd.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", gd.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", gd.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", gd.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", gd.UnitEnergy_in_cgs);
      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);

#ifdef ISOTHERM_EQS
  gd.MinEgySpec = 0;
#else
  gd.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * gd.MinGasTemp;
  gd.MinEgySpec *= gd.UnitMass_in_g / gd.UnitEnergy_in_cgs;
#endif

}



void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(ThisTask != 0)
    return;

//  if(RestartFlag == 0)
  if(cmd.RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");


  sprintf(buf, "%s%s", gd.OutputDir, gd.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun_mpi(0,1);
    }

  sprintf(buf, "%s%s", gd.OutputDir, gd.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun_mpi(0,1);
    }

  sprintf(buf, "%s%s", gd.OutputDir, gd.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun_mpi(0,1);
    }

  sprintf(buf, "%s%s", gd.OutputDir, gd.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun_mpi(0,1);
    }

#ifdef FORCETEST
//  if(RestartFlag == 0)
  if(cmd.RestartFlag == 0)
    {
      sprintf(buf, "%s%s", gd.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun_mpi(0,1);
	}
      fclose(FdForceTest);
    }
#endif
}


void close_outputfiles(void)
{
  if(ThisTask != 0)
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
#ifdef FORCETEST
  fclose(FdForceTest);
#endif
}


//-------------------------------------------------------------------
void read_parameter_file(char *fname)
{
	int  errorFlag = 0;
	string strtmp;
    
	if(sizeof(long long) != 8) {
		if(ThisTask == 0)
			printf("\n%s. %s.\n\n",
                   "Type `long long' is not 64 bit on this platform", "Stopping");
		endrun_mpi(ThisTask,0);
    }
    
	if(sizeof(int) != 4) {
		if(ThisTask == 0)
			printf("\n%s. %s.\n\n",
                   "Type `int' is not 32 bit on this platform", "Stopping");
		endrun_mpi(ThisTask,0);
    }
    
	if(sizeof(float) != 4) {
		if(ThisTask == 0)
			printf("\n%s. %s.\n\n",
                   "Type `float' is not 32 bit on this platform", "Stopping");
		endrun_mpi(ThisTask,0);
    }
    
	if(sizeof(double) != 8) {
		if(ThisTask == 0)
			printf("\n%s. %s.\n\n",
                   "Type `double' is not 64 bit on this platform", "Stopping");
		endrun_mpi(ThisTask,0);
    }
    
	if(ThisTask == 0) {
		GetCharParam(ParameterFile,"ParameterFile");
		if (!strnull(ParameterFile)) {
			fprintf(stdout,"\nParameter file: %s\n",ParameterFile);
//			RestartFlag = GetiParam("RestartFlag");
//fprintf_mpi_flush(STDOUT, ThisTask, "\nread_parameter_file:: RestartFlag %d\n", cmd.RestartFlag);
			cmd.RestartFlag = GetiParam("RestartFlag");
//fprintf_mpi_flush(STDOUT, ThisTask, "\nread_parameter_file:: RestartFlag %d\n", cmd.RestartFlag);
			errorFlag = startrun_parameterfile(ParameterFile);
		} else {
			errorFlag = startrun_cmdline();
		}
	}
    
//fprintf_mpi_flush(STDOUT, ThisTask, "\nread_parameter_file:: RestartFlag %d\n", cmd.RestartFlag);
//	MPI_Bcast(&RestartFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cmd.RestartFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
//fprintf_mpi_flush(STDOUT, ThisTask, "\nread_parameter_file:: RestartFlag %d\n", cmd.RestartFlag);

	MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
	if(errorFlag) {
		MPI_Finalize();
		exit(0);
    }
    
	MPI_Bcast(&gd, sizeof(global_data_all_processes),
              MPI_BYTE, 0, MPI_COMM_WORLD);
    
	if (strnull(gd.InitCondFile)) {
		if (ThisTask==0)
			fprintf(stdout,"\n\nYou must supply an IC file\n\n");
		MPI_Finalize();
		exit(0);
    }
    
	if(gd.NumFilesWrittenInParallel < 1) {
		if(ThisTask == 0)
			printf("NumFilesWrittenInParallel MUST be at least 1\n");
		endrun_mpi(ThisTask,0);
    }
    
	if(gd.NumFilesWrittenInParallel > NTask) {
		if(ThisTask == 0)
			printf("NumFilesWrittenInParallel MUST %s\n",
                   "be smaller than number of processors");
		endrun_mpi(ThisTask,0);
    }
    
/* #ifdef PERIODIC
	if(gd.PeriodicBoundariesOn == 0) {
		if(ThisTask == 0) {
			printf("%s %s.\n",
                   "Code was compiled with",
                   "periodic boundary conditions switched on");
			printf("You must set `PeriodicBoundariesOn=1', %s.\n",
                   "or recompile the code");
		}
		endrun_mpi(ThisTask,0);
    }
#else */
	if(gd.PeriodicBoundariesOn == 1) {
		if(ThisTask == 0) {
			printf("%s %s.\n",
                   "Code was compiled with",
                   "periodic boundary conditions switched off");
			printf("You must set `PeriodicBoundariesOn=0', %s.\n",
                   "or recompile the code");
		}
		endrun_mpi(ThisTask,0);
    }
//#endif
    
	if(gd.TypeOfTimestepCriterion >= 1) {
		if(ThisTask == 0) {
			printf("The specified timestep criterion\n");
			printf("is not valid\n");
		}
		endrun_mpi(ThisTask,0);
    }
    
#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
	if(ThisTask == 0) {
		printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
		printf("%s.\n",
               "Stretched periodic boxes are not implemented for gravity yet");
    }
	endrun_mpi(ThisTask,0);
#endif
#endif
}
//-------------------------------------------------------------------

local int startrun_parameterfile(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag = 0;

      nt = 0;

      SPName(gd.InitCondFile, "InitCondFile");
        SPName(gd.OutputDir, "OutputDir");
        SPName(gd.SnapshotFileBase, "SnapshotFileBase");
        SPName(gd.EnergyFile, "EnergyFile");
        SPName(gd.CpuFile, "CpuFile");
        SPName(gd.InfoFile, "InfoFile");
        SPName(gd.TimingsFile, "TimingsFile");
        SPName(gd.RestartFile, "RestartFile");
        SPName(gd.ResubmitCommand, "ResubmitCommand");
        SPName(gd.OutputListFilename, "OutputListFilename");
        IPName(gd.OutputListOn, "OutputListOn");
        RPName(gd.Omega0, "Omega0");
        RPName(gd.OmegaBaryon, "OmegaBaryon");
        RPName(gd.OmegaLambda, "OmegaLambda");
        RPName(gd.HubbleParam, "HubbleParam");
        RPName(gd.BoxSize, "BoxSize");
        IPName(gd.PeriodicBoundariesOn, "PeriodicBoundariesOn");
        RPName(gd.TimeOfFirstSnapshot,"TimeOfFirstSnapshot");
        RPName(gd.CpuTimeBetRestartFile,"CpuTimeBetRestartFile");
        RPName(gd.TimeBetStatistics,"TimeBetStatistics");
        RPName(gd.TimeBegin,"TimeBegin");
        RPName(gd.TimeMax,"TimeMax");
        RPName(gd.TimeBetSnapshot,"TimeBetSnapshot");
        RPName(gd.UnitVelocity_in_cm_per_s,"UnitVelocity_in_cm_per_s");
        RPName(gd.UnitLength_in_cm,"UnitLength_in_cm");
        RPName(gd.UnitMass_in_g,"UnitMass_in_g");
        RPName(gd.TreeDomainUpdateFrequency,"TreeDomainUpdateFrequency");
        RPName(gd.ErrTolIntAccuracy,"ErrTolIntAccuracy");
        RPName(gd.ErrTolTheta,"ErrTolTheta");
        RPName(gd.ErrTolForceAcc,"ErrTolForceAcc");
        RPName(gd.MinGasHsmlFractional,"MinGasHsmlFractional");
        RPName(gd.MaxSizeTimestep,"MaxSizeTimestep");
        RPName(gd.MinSizeTimestep,"MinSizeTimestep");
        RPName(gd.MaxRMSDisplacementFac,"MaxRMSDisplacementFac");
        RPName(gd.ArtBulkViscConst,"ArtBulkViscConst");
        RPName(gd.CourantFac,"CourantFac");
        RPName(gd.DesNumNgb,"DesNumNgb");
        RPName(gd.MaxNumNgbDeviation,"MaxNumNgbDeviation");
        IPName(gd.ComovingIntegrationOn,"ComovingIntegrationOn");
        IPName(gd.ICFormat,"ICFormat");
        IPName(gd.SnapFormat,"SnapFormat");
        IPName(gd.NumFilesPerSnapshot,"NumFilesPerSnapshot");
        IPName(gd.NumFilesWrittenInParallel,"NumFilesWrittenInParallel");
        IPName(gd.ResubmitOn,"ResubmitOn");
        IPName(gd.TypeOfTimestepCriterion,"TypeOfTimestepCriterion");
        IPName(gd.TypeOfOpeningCriterion,"TypeOfOpeningCriterion");
        RPName(gd.TimeLimitCPU,"TimeLimitCPU");
        RPName(gd.SofteningHalo,"SofteningHalo");
        RPName(gd.SofteningDisk,"SofteningDisk");
        RPName(gd.SofteningBulge,"SofteningBulge");
        RPName(gd.SofteningGas,"SofteningGas");
        RPName(gd.SofteningStars,"SofteningStars");
        RPName(gd.SofteningBndry,"SofteningBndry");
        RPName(gd.SofteningHaloMaxPhys,"SofteningHaloMaxPhys");
        RPName(gd.SofteningDiskMaxPhys,"SofteningDiskMaxPhys");
        RPName(gd.SofteningBulgeMaxPhys,"SofteningBulgeMaxPhys");
        RPName(gd.SofteningGasMaxPhys,"SofteningGasMaxPhys");
        RPName(gd.SofteningStarsMaxPhys,"SofteningStarsMaxPhys");
        RPName(gd.SofteningBndryMaxPhys,"SofteningBndryMaxPhys");
        IPName(gd.BufferSize,"BufferSize");
        RPName(gd.PartAllocFactor,"PartAllocFactor");
        RPName(gd.TreeAllocFactor,"TreeAllocFactor");
        RPName(gd.GravityConstantInternal,"GravityConstantInternal");
        RPName(gd.InitGasTemp,"InitGasTemp");
        RPName(gd.MinGasTemp,"MinGasTemp");

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(gd.OutputDir);
	      if(i > 0)
		if(gd.OutputDir[i - 1] != '/')
		  strcat(gd.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", gd.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
	      system(buf3);
	    }
	}
      else
	{
	  printf("\nParameter file %s not found.\n\n", fname);
	  errorFlag = 2;
	}

      if(errorFlag != 2)
	for(i = 0; i < nt; i++)
	  {
	    if(*tag[i])
	      {
		printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
		errorFlag = 1;
	      }
	  }

      if(gd.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(gd.OutputListFilename);
      else
	gd.OutputListLength = 0;

return errorFlag;

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}


int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  gd.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &gd.OutputListTimes[gd.OutputListLength]) == 1)
	gd.OutputListLength++;
      else
	break;
    }
  while(gd.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", gd.OutputListLength);

  return 0;
}


void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(ThisTask == 0)
    {
      printf("\ngd.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce gd.TimeMax\n\n");
      endrun_mpi(ThisTask,556);
    }

  if(gd.ComovingIntegrationOn)
    ti_end = log(TimeMax_new / gd.TimeBegin) / gd.Timebase_interval;
  else
    ti_end = (TimeMax_new - gd.TimeBegin) / gd.Timebase_interval;

  while(ti_end > TIMEBASE)
    {
      gd.Timebase_interval *= 2.0;

      ti_end /= 2;
      gd.Ti_Current /= 2;

#ifdef PMGRID
      gd.PM_Ti_begstep /= 2;
      gd.PM_Ti_endstep /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_endstep /= 2;
	}
    }

  gd.TimeMax = TimeMax_new;
}

#define parameter_null	"parameters_null"

local int startrun_cmdline(void)
{
	int  errorFlag=0;
	string strtmp;
    
	if(ThisTask == 0)
        fprintf(stdout,"\nInput command line parameters ...\n");
    
//	RestartFlag = GetiParam("RestartFlag");
	cmd.RestartFlag = GetiParam("RestartFlag");
fprintf_mpi_flush(STDOUT, ThisTask, "\nstartrun_cmdline:: RestartFlag %d\n", cmd.RestartFlag);

	gd.SnapFormat = GetiParam("SnapFormat");
	GetCharParam(gd.InitCondFile,"InitCondFile");
	GetCharParam(gd.OutputDir,"OutputDir");
	GetCharParam(gd.SnapshotFileBase,"SnapshotFileBase");
	GetCharParam(gd.EnergyFile,"EnergyFile");
	GetCharParam(gd.CpuFile,"CpuFile");
	GetCharParam(gd.InfoFile,"InfoFile");
	GetCharParam(gd.TimingsFile,"TimingsFile");
	GetCharParam(gd.RestartFile,"RestartFile");
    
	GetCharParam(gd.ResubmitCommand,"ResubmitCommand");
	GetCharParam(gd.OutputListFilename,"OutputListFilename");
	gd.OutputListOn = GetiParam("OutputListOn");
	gd.BoxSize = GetdParam("BoxSize");
	gd.PeriodicBoundariesOn = GetiParam("PeriodicBoundariesOn");
	gd.TimeOfFirstSnapshot = GetdParam("TimeOfFirstSnapshot");
	gd.CpuTimeBetRestartFile = GetdParam("CpuTimeBetRestartFile");
	gd.TimeBetStatistics = GetdParam("TimeBetStatistics");
	gd.TimeBegin = GetdParam("TimeBegin");
	gd.TimeMax = GetdParam("TimeMax");
	gd.TimeBetSnapshot = GetdParam("TimeBetSnapshot");
	gd.UnitVelocity_in_cm_per_s = GetdParam("UnitVelocity_in_cm_per_s");
	gd.UnitLength_in_cm = GetdParam("UnitLength_in_cm");
	gd.UnitMass_in_g = GetdParam("UnitMass_in_g");
	gd.GravityConstantInternal = GetdParam("GravityConstantInternal");
    
	gd.ErrTolIntAccuracy = GetdParam("ErrTolIntAccuracy");
    
	gd.ErrTolTheta = GetdParam("ErrTolTheta");
	gd.ErrTolForceAcc = GetdParam("ErrTolForceAcc");
	gd.MaxSizeTimestep = GetdParam("MaxSizeTimestep");
	gd.MinSizeTimestep = GetdParam("MinSizeTimestep");
	gd.MaxRMSDisplacementFac = GetdParam("MaxRMSDisplacementFac");
	gd.DesNumNgb = GetiParam("DesNumNgb");
	gd.MaxNumNgbDeviation = GetiParam("MaxNumNgbDeviation");
	gd.ICFormat = GetiParam("ICFormat");
	gd.Omega0 = GetdParam("Omega0");
	gd.OmegaBaryon = GetdParam("OmegaBaryon");
	gd.OmegaLambda = GetdParam("OmegaLambda");
	gd.HubbleParam = GetdParam("HubbleParam");
    
	gd.MinGasHsmlFractional = GetdParam("MinGasHsmlFractional");
	gd.ArtBulkViscConst = GetdParam("ArtBulkViscConst");
	gd.CourantFac = GetdParam("CourantFac");
	gd.ComovingIntegrationOn = GetiParam("ComovingIntegrationOn");
    
	gd.NumFilesPerSnapshot = GetiParam("NumFilesPerSnapshot");
	gd.NumFilesWrittenInParallel = GetiParam("NumFilesWrittenInParallel");
	gd.ResubmitOn = GetiParam("ResubmitOn");
	gd.TypeOfTimestepCriterion = GetiParam("TypeOfTimestepCriterion");
	gd.TypeOfOpeningCriterion = GetiParam("TypeOfOpeningCriterion");
	gd.TimeLimitCPU = GetdParam("TimeLimitCPU");
	gd.TreeDomainUpdateFrequency = GetdParam("TreeDomainUpdateFrequency");
    
	gd.SofteningHalo = GetdParam("SofteningHalo");
	gd.SofteningDisk = GetdParam("SofteningDisk");
	gd.SofteningBulge = GetdParam("SofteningBulge");
	gd.SofteningGas = GetdParam("SofteningGas");
	gd.SofteningStars = GetdParam("SofteningStars");
	gd.SofteningBndry = GetdParam("SofteningBndry");
    
	gd.SofteningHaloMaxPhys = GetdParam("SofteningHaloMaxPhys");
	gd.SofteningDiskMaxPhys = GetdParam("SofteningDiskMaxPhys");
	gd.SofteningBulgeMaxPhys = GetdParam("SofteningBulgeMaxPhys");
	gd.SofteningGasMaxPhys = GetdParam("SofteningGasMaxPhys");
	gd.SofteningStarsMaxPhys = GetdParam("SofteningStarsMaxPhys");
	gd.SofteningBndryMaxPhys = GetdParam("SofteningBndryMaxPhys");
    
	gd.BufferSize = GetiParam("BufferSize");
	gd.PartAllocFactor = GetdParam("PartAllocFactor");
	gd.TreeAllocFactor = GetdParam("TreeAllocFactor");
    
	gd.InitGasTemp = GetdParam("InitGasTemp");
	gd.MinGasTemp = GetdParam("MinGasTemp");
    
	if(gd.OutputListOn && errorFlag==0)
		errorFlag+= read_outputlist(gd.OutputListFilename);
	else
		gd.OutputListLength=0;
    
	PrintParameterFile(parameter_null);
    
	return errorFlag;
}

#undef parameter_null

#define FMTT	"%-35s%s\n"
#define FMTI	"%-35s%d\n"
#define FMTR	"%-35s%g\n"

local void PrintParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];
    
    sprintf(buf,"%s%s",fname,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
                "%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Parameter input file for:",gd.headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
                gd.headline3);
        fprintf(fdout,"%s\n%s\n",
                "%-------------------------------------------------------------------",
                "%");
        
        fprintf(fdout,FMTI,"SnapFormat",gd.SnapFormat);
        fprintf(fdout,FMTT,"InitCondFile",gd.InitCondFile);
        fprintf(fdout,FMTT,"OutputDir",gd.OutputDir);
        fprintf(fdout,FMTT,"SnapshotFileBase",gd.SnapshotFileBase);
        fprintf(fdout,FMTT,"EnergyFile",gd.EnergyFile);
        fprintf(fdout,FMTT,"CpuFile",gd.CpuFile);
        fprintf(fdout,FMTT,"InfoFile",gd.InfoFile);
        fprintf(fdout,FMTT,"TimingsFile",gd.TimingsFile);
        fprintf(fdout,FMTT,"RestartFile",gd.RestartFile);
        
        fprintf(fdout,FMTT,"ResubmitCommand",gd.ResubmitCommand);
        fprintf(fdout,FMTT,"OutputListFilename",gd.OutputListFilename);
        fprintf(fdout,FMTI,"OutputListOn",gd.OutputListOn);
        fprintf(fdout,FMTR,"BoxSize",gd.BoxSize);
        fprintf(fdout,FMTI,"PeriodicBoundariesOn",gd.PeriodicBoundariesOn);
        fprintf(fdout,FMTR,"TimeOfFirstSnapshot",gd.TimeOfFirstSnapshot);
        fprintf(fdout,FMTR,"CpuTimeBetRestartFile",gd.CpuTimeBetRestartFile);
        fprintf(fdout,FMTR,"TimeBetStatistics",gd.TimeBetStatistics);
        fprintf(fdout,FMTR,"TimeBegin",gd.TimeBegin);
        fprintf(fdout,FMTR,"TimeMax",gd.TimeMax);
        fprintf(fdout,FMTR,"TimeBetSnapshot",gd.TimeBetSnapshot);
        fprintf(fdout,FMTR,"UnitVelocity_in_cm_per_s",gd.UnitVelocity_in_cm_per_s);
        fprintf(fdout,FMTR,"UnitLength_in_cm",gd.UnitLength_in_cm);
        fprintf(fdout,FMTR,"UnitMass_in_g",gd.UnitMass_in_g);
        fprintf(fdout,FMTR,"GravityConstantInternal",gd.GravityConstantInternal);
        fprintf(fdout,FMTR,"ErrTolIntAccuracy",gd.ErrTolIntAccuracy);
        fprintf(fdout,FMTR,"ErrTolTheta",gd.ErrTolTheta);
        fprintf(fdout,FMTR,"ErrTolForceAcc",gd.ErrTolForceAcc);
        fprintf(fdout,FMTR,"MaxSizeTimestep",gd.MaxSizeTimestep);
        fprintf(fdout,FMTR,"MinSizeTimestep",gd.MinSizeTimestep);
        fprintf(fdout,FMTR,"MaxRMSDisplacementFac",gd.MaxRMSDisplacementFac);
        fprintf(fdout,FMTR,"DesNumNgb",gd.DesNumNgb);
        fprintf(fdout,FMTR,"MaxNumNgbDeviation",gd.MaxNumNgbDeviation);
        fprintf(fdout,FMTI,"ICFormat",gd.ICFormat);
        fprintf(fdout,FMTR,"Omega0",gd.Omega0);
        fprintf(fdout,FMTR,"OmegaBaryon",gd.OmegaBaryon);
        fprintf(fdout,FMTR,"OmegaLambda",gd.OmegaLambda);
        fprintf(fdout,FMTR,"HubbleParam",gd.HubbleParam);
        
        fprintf(fdout,FMTR,"MinGasHsmlFractional",gd.MinGasHsmlFractional);
        fprintf(fdout,FMTR,"ArtBulkViscConst",gd.ArtBulkViscConst);
        fprintf(fdout,FMTR,"CourantFac",gd.CourantFac);
        fprintf(fdout,FMTI,"ComovingIntegrationOn",gd.ComovingIntegrationOn);
        
        fprintf(fdout,FMTI,"NumFilesPerSnapshot",gd.NumFilesPerSnapshot);
        fprintf(fdout,FMTI,"NumFilesWrittenInParallel",gd.NumFilesWrittenInParallel);
        fprintf(fdout,FMTI,"ResubmitOn",gd.ResubmitOn);
        fprintf(fdout,FMTI,"TypeOfTimestepCriterion",gd.TypeOfTimestepCriterion);
        fprintf(fdout,FMTI,"TypeOfOpeningCriterion",gd.TypeOfOpeningCriterion);
        fprintf(fdout,FMTR,"TimeLimitCPU",gd.TimeLimitCPU);
        fprintf(fdout,FMTR,"TreeDomainUpdateFrequency",gd.TreeDomainUpdateFrequency);
        
        fprintf(fdout,FMTR,"SofteningHalo",gd.SofteningHalo);
        fprintf(fdout,FMTR,"SofteningDisk",gd.SofteningDisk);
        fprintf(fdout,FMTR,"SofteningBulge",gd.SofteningBulge);
        fprintf(fdout,FMTR,"SofteningGas",gd.SofteningGas);
        fprintf(fdout,FMTR,"SofteningStars",gd.SofteningStars);
        fprintf(fdout,FMTR,"SofteningBndry",gd.SofteningBndry);
        
        fprintf(fdout,FMTR,"SofteningHaloMaxPhys",gd.SofteningHaloMaxPhys);
        fprintf(fdout,FMTR,"SofteningDiskMaxPhys",gd.SofteningDiskMaxPhys);
        fprintf(fdout,FMTR,"SofteningBulgeMaxPhys",gd.SofteningBulgeMaxPhys);
        fprintf(fdout,FMTR,"SofteningGasMaxPhys",gd.SofteningGasMaxPhys);
        fprintf(fdout,FMTR,"SofteningStarsMaxPhys",gd.SofteningStarsMaxPhys);
        fprintf(fdout,FMTR,"SofteningBndryMaxPhys",gd.SofteningBndryMaxPhys);
        
        fprintf(fdout,FMTI,"BufferSize",gd.BufferSize);
        fprintf(fdout,FMTR,"PartAllocFactor",gd.PartAllocFactor);
        fprintf(fdout,FMTR,"TreeAllocFactor",gd.TreeAllocFactor);
        
        fprintf(fdout,FMTR,"InitGasTemp",gd.InitGasTemp);
        fprintf(fdout,FMTR,"MinGasTemp",gd.MinGasTemp);
        
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR
