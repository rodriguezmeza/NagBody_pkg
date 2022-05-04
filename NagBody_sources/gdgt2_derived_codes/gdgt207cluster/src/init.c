
// DEJAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS

//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <mpi.h>

#include "globaldefs.h"
#include "protodefs.h"
//#include "../../../General_libs/mpi/mpi_proto.h"


void init(void)
{
  int i, j;
  double a3;

  gd.Time = gd.TimeBegin;

  switch (gd.ICFormat)
    {
    case 1:
#if (MAKEGLASS > 1)
      seed_glass();
#else
      read_ic(gd.InitCondFile);
#endif
      break;
    case 2:
    case 3:
      read_ic(gd.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", gd.ICFormat);
      endrun_mpi(0,0);
    }

  gd.Time = gd.TimeBegin;
  gd.Ti_Current = 0;

  if(gd.ComovingIntegrationOn)
    {
      gd.Timebase_interval = (log(gd.TimeMax) - log(gd.TimeBegin)) / TIMEBASE;
      a3 = gd.Time * gd.Time * gd.Time;
    }
  else
    {
      gd.Timebase_interval = (gd.TimeMax - gd.TimeBegin) / TIMEBASE;
      a3 = 1;
    }

  set_softenings();

  gd.NumCurrentTiStep = 0;
  gd.SnapshotFileCount = 0;
//    if(RestartFlag == 2)
    if(cmd.RestartFlag == 2)
    gd.SnapshotFileCount = atoi(gd.InitCondFile + strlen(gd.InitCondFile) - 3) + 1;

  gd.TotNumOfForces = 0;
  gd.NumForcesSinceLastDomainDecomp = 0;

  if(gd.ComovingIntegrationOn)
    if(gd.PeriodicBoundariesOn == 1)
      check_omega();

  gd.TimeLastStatistics = gd.TimeBegin - gd.TimeBetStatistics;

  if(gd.ComovingIntegrationOn)
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(gd.Time) * gd.Time;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].Potential = 0;
    }

#ifdef PMGRID
  gd.PM_Ti_endstep = gd.PM_Ti_begstep = 0;
#endif

#ifdef FLEXSTEPS
  gd.PresentMinStep = TIMEBASE;
  for(i = 0; i < NumPart; i++)
    {
      P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif


  for(i = 0; i < N_gas; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = 0;
	}

      SphP[i].DtEntropy = 0;

//        if(RestartFlag == 0)
        if(cmd.RestartFlag == 0)
	{
	  SphP[i].Hsml = 0;
	  SphP[i].Density = -1;
	}
    }

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(gd.TreeAllocFactor * gd.MaxPart, gd.MaxPart);

  gd.NumForcesSinceLastDomainDecomp = 1 + gd.TotNumPart * gd.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;

  domain_Decomposition();

  ngb_treebuild();

  setup_smoothinglengths();

  TreeReconstructFlag = 1;

#ifndef ISOTHERM_EQS
  if(header.flag_entropy_instead_u == 0)
    for(i = 0; i < N_gas; i++)
      SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif
}


void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (gd.BoxSize * gd.BoxSize * gd.BoxSize) / (3 * gd.Hubble * gd.Hubble / (8 * M_PI * gd.G));

  if(fabs(omega - gd.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, gd.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun_mpi(0,1);
    }
}


void setup_smoothinglengths(void)
{
  int i, no, p;

//    if(RestartFlag == 0)
    if(cmd.RestartFlag == 0)
    {

      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * gd.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }
#ifndef TWODIMS
	  SphP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * gd.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  SphP[i].Hsml =
	    pow(1.0 / (M_PI) * gd.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	}
    }

  density();
}


#if (MAKEGLASS > 1)
void seed_glass(void)
{
  int i, k, n_for_this_task;
  double Range[3], LowerBound[3];
  double drandom, partmass;
  long long IDstart;

  gd.TotNumPart = MAKEGLASS;
  partmass = gd.Omega0 * (3 * gd.Hubble * gd.Hubble / (8 * M_PI * gd.G))
    * (gd.BoxSize * gd.BoxSize * gd.BoxSize) / gd.TotNumPart;

  gd.MaxPart = gd.PartAllocFactor * (gd.TotNumPart / NTask);

  allocate_memory();

  header.npartTotal[1] = gd.TotNumPart;
  header.mass[1] = partmass;

  if(ThisTask == 0)
    {
      printf("\nGlass initialising\nPartMass= %g\n", partmass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (gd.TotNumPart / 1000000000), (int) (gd.TotNumPart % 1000000000));
    }

  n_for_this_task = gd.TotNumPart / NTask;

  if(ThisTask == NTask - 1)
    n_for_this_task = gd.TotNumPart - (NTask - 1) * n_for_this_task;

  NumPart = 0;
  IDstart = 1 + (gd.TotNumPart / NTask) * ThisTask;

  Range[0] = Range[1] = gd.BoxSize;
  Range[2] = gd.BoxSize / NTask;
  LowerBound[0] = LowerBound[1] = 0;
  LowerBound[2] = ThisTask * Range[2];

  srand48(ThisTask);

  for(i = 0; i < n_for_this_task; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  drandom = drand48();

	  P[i].Pos[k] = LowerBound[k] + Range[k] * drandom;
	  P[i].Vel[k] = 0;
	}

      P[i].Mass = partmass;
      P[i].Type = 1;
      P[i].ID = IDstart + i;

      NumPart++;
    }
}
#endif
