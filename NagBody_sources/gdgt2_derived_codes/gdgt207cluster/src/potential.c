
// DEJAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <float.h>
//#include <mpi.h>


#include "globaldefs.h"
#include "protodefs.h"

//#include "../../../General_libs/mpi/mpi_proto.h"


void compute_potential(void)
{
  int i;

#ifndef NOGRAVITY
  long long ntot, ntotleft;
  int j, k, level, sendTask, recvTask;
  int ndone;
  int maxfill, ngrp, place, nexport;
  int *nsend, *noffset, *nsend_local, *nbuffer, *ndonelist, *numlist;
  double fac;
  double t0, t1, tstart, tend;
  MPI_Status status;
  double r2;

  t0 = second();

  if(gd.ComovingIntegrationOn)
    set_softenings();

  if(ThisTask == 0)
    {
      printf("Start computation of potential for all particles...\n");
      fflush(stdout);
    }


  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

      force_treebuild(NumPart);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }
  tend = second();
  gd.CPU_TreeConstruction += timediff(tstart, tend);

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumPart, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);

  noffset = malloc(sizeof(int) * NTask);
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;
  ntotleft = ntot;

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      for(nexport = 0, ndone = 0; i < NumPart && nexport < gd.BunchSizeForce - NTask; i++)
	{
	  ndone++;

	  for(j = 0; j < NTask; j++)
	    Exportflag[j] = 0;

#ifndef PMGRID
	  force_treeevaluate_potential(i, 0);
#else
	  force_treeevaluate_potential_shortrange(i, 0);
#endif

	  for(j = 0; j < NTask; j++)
	    {
	      if(Exportflag[j])
		{
		  for(k = 0; k < 3; k++)
		    GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
#ifdef UNEQUALSOFTENINGS
		  GravDataGet[nexport].Type = P[i].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                  if(P[i].Type == 0)
                    GravDataGet[nexport].Soft = SphP[i].Hsml;
#endif
#endif
		  GravDataGet[nexport].w.OldAcc = P[i].OldAcc;

		  GravDataIndexTable[nexport].Task = j;
		  GravDataIndexTable[nexport].Index = i;
		  GravDataIndexTable[nexport].SortIndex = nexport;

		  nexport++;
		  nsend_local[j]++;
		}
	    }
	}

      qsort(GravDataIndexTable, nexport, sizeof(gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);


      for(level = 1; level < (1 << PTask); level++)
	{
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= gd.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(gravdata_in), MPI_BYTE,
				   recvTask, TAG_POTENTIAL_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(gravdata_in), MPI_BYTE,
				   recvTask, TAG_POTENTIAL_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
#ifndef PMGRID
	      force_treeevaluate_potential(j, 1);
#else
	      force_treeevaluate_potential_shortrange(j, 1);
#endif
	    }


	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= gd.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(gravdata_in),
				   MPI_BYTE, recvTask, TAG_POTENTIAL_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(gravdata_in),
				   MPI_BYTE, recvTask, TAG_POTENTIAL_B, MPI_COMM_WORLD, &status);

		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  P[place].Potential += GravDataOut[j + noffset[recvTask]].u.Potential;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);



  for(i = 0; i < NumPart; i++)
    {
      P[i].Potential += P[i].Mass / gd.SofteningTable[P[i].Type];

      if(gd.ComovingIntegrationOn)
	if(gd.PeriodicBoundariesOn)
	  P[i].Potential -= 2.8372975 * pow(P[i].Mass, 2.0 / 3) *
	    pow(gd.Omega0 * 3 * gd.Hubble * gd.Hubble / (8 * M_PI * gd.G), 1.0 / 3);
    }



  for(i = 0; i < NumPart; i++)
    P[i].Potential *= gd.G;


#ifdef PMGRID

/* #ifdef PERIODIC
  pmpotential_periodic();
#ifdef PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(1);
    }
  if(i == 1)
    endrun_mpi(0,88686);
#endif
#else */
  i = pmpotential_nonperiodic(0);
  if(i == 1)
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmpotential_nonperiodic(0);
    }
  if(i == 1)
    endrun_mpi(0,88687);
#ifdef PLACEHIGHRESREGION
  i = pmpotential_nonperiodic(1);
  if(i == 1)
    {
      pm_init_regionsize();

      i = pmpotential_nonperiodic(1);
    }
  if(i != 0)
    endrun_mpi(0,88688);
#endif
//#endif

#endif



  if(gd.ComovingIntegrationOn)
    {
// #ifndef PERIODIC
      fac = -0.5 * gd.Omega0 * gd.Hubble * gd.Hubble;

      for(i = 0; i < NumPart; i++)
	{
	  for(k = 0, r2 = 0; k < 3; k++)
	    r2 += P[i].Pos[k] * P[i].Pos[k];

	  P[i].Potential += fac * r2;
	}
//#endif
    }
  else
    {
      fac = -0.5 * gd.OmegaLambda * gd.Hubble * gd.Hubble;
      if(fac != 0)
	{
	  for(i = 0; i < NumPart; i++)
	    {
	      for(k = 0, r2 = 0; k < 3; k++)
		r2 += P[i].Pos[k] * P[i].Pos[k];

	      P[i].Potential += fac * r2;
	    }
	}
    }


  if(ThisTask == 0)
    {
      printf("potential done.\n");
      fflush(stdout);
    }

  t1 = second();

  gd.CPU_Potential += timediff(t0, t1);

#else
  for(i = 0; i < NumPart; i++)
    P[i].Potential = 0;
#endif
}
