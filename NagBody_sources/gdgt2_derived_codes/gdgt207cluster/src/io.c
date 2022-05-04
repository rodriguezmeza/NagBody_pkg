
// DEJAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <mpi.h>
//#include <errno.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "globaldefs.h"
#include "protodefs.h"
//#include "../../../General_libs/mpi/mpi_proto.h"


static int n_type[6];
static long long ntot_type_all[6];


void savepositions(int num)
{
  double t0, t1;
  char buf[500];
  int i, j, *temp, n, filenr, gr, ngroups, masterTask, lastTask;

  t0 = second();

  if(ThisTask == 0)
    printf("\nwriting snapshot file... \n");

#if defined(SFR) || defined(BLACK_HOLES)
  rearrange_particle_sequence();

  gd.NumForcesSinceLastDomainDecomp = 1 + gd.TreeDomainUpdateFrequency * gd.TotNumPart;
#endif

  if(NTask < gd.NumFilesPerSnapshot)
    {
      if(ThisTask == 0)
	printf("Fatal error.\nNumber of processors must be larger or equal than gd.NumFilesPerSnapshot.\n");
      endrun_mpi(0,0);
    }
  if(gd.SnapFormat < 1 || gd.SnapFormat > 3)
    {
      if(ThisTask == 0)
	printf("Unsupported File-Format\n");
      endrun_mpi(0,0);
    }
#ifndef  HAVE_HDF5
  if(gd.SnapFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun_mpi(0,0);
    }
#endif

  for(n = 0; n < 6; n++)
    n_type[n] = 0;

  for(n = 0; n < NumPart; n++)
    n_type[P[n].Type]++;

  temp = malloc(NTask * 6 * sizeof(int));
  MPI_Allgather(n_type, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
  for(i = 0; i < 6; i++)
    {
      ntot_type_all[i] = 0;
      for(j = 0; j < NTask; j++)
	ntot_type_all[i] += temp[j * 6 + i];
    }
  free(temp);


  distribute_file(gd.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

  fill_Tab_IO_Labels();

  if(gd.NumFilesPerSnapshot > 1)
    sprintf(buf, "%s%s_%03d.%d", gd.OutputDir, gd.SnapshotFileBase, num, filenr);
  else
    sprintf(buf, "%s%s_%03d", gd.OutputDir, gd.SnapshotFileBase, num);

  ngroups = gd.NumFilesPerSnapshot / gd.NumFilesWrittenInParallel;
  if((gd.NumFilesPerSnapshot % gd.NumFilesWrittenInParallel))
    ngroups++;

  for(gr = 0; gr < ngroups; gr++)
    {
      if((filenr / gd.NumFilesWrittenInParallel) == gr)
	write_file(buf, masterTask, lastTask);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(ThisTask == 0)
    printf("done with snapshot.\n");

  t1 = second();

  gd.CPU_Snapshot += timediff(t0, t1);

}


void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex;
  float *fp;

#ifdef LONGIDS
  long long *ip;
#else
  int *ip;
#endif

/* #ifdef PERIODIC
  FLOAT boxSize;
#endif */
#ifdef PMGRID
  double dt_gravkick_pm = 0;
#endif
  double dt_gravkick, dt_hydrokick, a3inv = 1, fac1, fac2;

  if(gd.ComovingIntegrationOn)
    {
      a3inv = 1 / (gd.Time * gd.Time * gd.Time);
      fac1 = 1 / (gd.Time * gd.Time);
      fac2 = 1 / pow(gd.Time, 3 * GAMMA - 2);
    }
  else
    a3inv = fac1 = fac2 = 1;

#ifdef PMGRID
  if(gd.ComovingIntegrationOn)
    dt_gravkick_pm =
      get_gravkick_factor(gd.PM_Ti_begstep,
			  gd.Ti_Current) -
      get_gravkick_factor(gd.PM_Ti_begstep, (gd.PM_Ti_begstep + gd.PM_Ti_endstep) / 2);
  else
    dt_gravkick_pm = (gd.Ti_Current - (gd.PM_Ti_begstep + gd.PM_Ti_endstep) / 2) * gd.Timebase_interval;
#endif

  fp = CommBuffer;
  ip = CommBuffer;

  pindex = *startindex;

  switch (blocknr)
    {
    case IO_POS:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Pos[k];
/*#ifdef PERIODIC
		boxSize = gd.BoxSize;
#ifdef LONG_X
		if(k == 0)
		  boxSize = gd.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
		if(k == 1)
		  boxSize = gd.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
		if(k == 2)
		  boxSize = gd.BoxSize * LONG_Z;
#endif
		while(fp[k] < 0)
		  fp[k] += boxSize;
		while(fp[k] >= boxSize)
		  fp[k] -= boxSize;
#endif */
	      }
	    n++;
	    fp += 3;
	  }
      break;

    case IO_VEL:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    if(gd.ComovingIntegrationOn)
	      {
		dt_gravkick =
		  get_gravkick_factor(P[pindex].Ti_begstep,
				      gd.Ti_Current) -
		  get_gravkick_factor(P[pindex].Ti_begstep,
				      (P[pindex].Ti_begstep + P[pindex].Ti_endstep) / 2);
		dt_hydrokick =
		  get_hydrokick_factor(P[pindex].Ti_begstep,
				       gd.Ti_Current) -
		  get_hydrokick_factor(P[pindex].Ti_begstep,
				       (P[pindex].Ti_begstep + P[pindex].Ti_endstep) / 2);
	      }
	    else
	      dt_gravkick = dt_hydrokick =
		(gd.Ti_Current - (P[pindex].Ti_begstep + P[pindex].Ti_endstep) / 2) * gd.Timebase_interval;

	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Vel[k] + P[pindex].GravAccel[k] * dt_gravkick;
		if(P[pindex].Type == 0)
		  fp[k] += SphP[pindex].HydroAccel[k] * dt_hydrokick;
	      }
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += P[pindex].GravPM[k] * dt_gravkick_pm;
#endif
	    for(k = 0; k < 3; k++)
	      fp[k] *= sqrt(a3inv);

	    n++;
	    fp += 3;
	  }
      break;

    case IO_ID:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *ip++ = P[pindex].ID;
	    n++;
	  }
      break;

    case IO_MASS:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Mass;
	    n++;
	  }
      break;

    case IO_U:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
#ifdef ISOTHERM_EQS
	    *fp++ = SphP[pindex].Entropy;
#else
	    *fp++ =
	      dmax(gd.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].Density * a3inv, GAMMA_MINUS1));
#endif
	    n++;
	  }
      break;

    case IO_RHO:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Density;
	    n++;
	  }
      break;

    case IO_HSML:
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].Hsml;
	    n++;
	  }
      break;

    case IO_POT:
#ifdef OUTPUTPOTENTIAL
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Potential;
	    n++;
	  }
#endif
      break;

    case IO_ACCEL:
#ifdef OUTPUTACCELERATION
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      fp[k] = fac1 * P[pindex].GravAccel[k];
#ifdef PMGRID
	    for(k = 0; k < 3; k++)
	      fp[k] += fac1 * P[pindex].GravPM[k];
#endif
	    if(P[pindex].Type == 0)
	      for(k = 0; k < 3; k++)
		fp[k] += fac2 * SphP[pindex].HydroAccel[k];
	    fp += 3;
	    n++;
	  }
#endif
      break;

    case IO_DTENTR:
#ifdef OUTPUTCHANGEOFENTROPY
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].DtEntropy;
	    n++;
	  }
#endif
      break;

    case IO_TSTP:
#ifdef OUTPUTTIMESTEP

      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = (P[pindex].Ti_endstep - P[pindex].Ti_begstep) * gd.Timebase_interval;
	    n++;
	  }
#endif
      break;

    }

  *startindex = pindex;
}


int get_bytes_per_blockelement(enum iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      bytes_per_blockelement = 3 * sizeof(float);
      break;

    case IO_ID:
#ifdef LONGIDS
      bytes_per_blockelement = sizeof(long long);
#else
      bytes_per_blockelement = sizeof(int);
#endif
      break;

    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
      bytes_per_blockelement = sizeof(float);
      break;
    }

  return bytes_per_blockelement;
}


int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifdef LONGIDS
      typekey = 2;
#else
      typekey = 0;
#endif
      break;

    default:
      typekey = 1;
      break;
    }

  return typekey;
}


int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
      values = 1;
      break;
    }

  return values;
}


int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(gd.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(gd.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_DTENTR:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngas;
      break;
    }

  endrun_mpi(0,212);
  return 0;
}


int blockpresent(enum iofields blocknr)
{

#ifndef OUTPUTPOTENTIAL
  if(blocknr == IO_POT)
    return 0;
#endif

#ifndef OUTPUTACCELERATION
  if(blocknr == IO_ACCEL)
    return 0;
#endif

#ifndef OUTPUTCHANGEOFENTROPY
  if(blocknr == IO_DTENTR)
    return 0;
#endif

#ifndef OUTPUTTIMESTEP
  if(blocknr == IO_TSTP)
    return 0;
#endif

  return 1;
}


void fill_Tab_IO_Labels(void)
{
  enum iofields i;

  for(i = 0; i < IO_NBLOCKS; i++)
    switch (i)
      {
      case IO_POS:
	strncpy(Tab_IO_Labels[IO_POS], "POS ", 4);
	break;
      case IO_VEL:
	strncpy(Tab_IO_Labels[IO_VEL], "VEL ", 4);
	break;
      case IO_ID:
	strncpy(Tab_IO_Labels[IO_ID], "ID  ", 4);
	break;
      case IO_MASS:
	strncpy(Tab_IO_Labels[IO_MASS], "MASS", 4);
	break;
      case IO_U:
	strncpy(Tab_IO_Labels[IO_U], "U   ", 4);
	break;
      case IO_RHO:
	strncpy(Tab_IO_Labels[IO_RHO], "RHO ", 4);
	break;
      case IO_HSML:
	strncpy(Tab_IO_Labels[IO_HSML], "HSML", 4);
	break;
      case IO_POT:
	strncpy(Tab_IO_Labels[IO_POT], "POT ", 4);
	break;
      case IO_ACCEL:
	strncpy(Tab_IO_Labels[IO_ACCEL], "ACCE", 4);
	break;
      case IO_DTENTR:
	strncpy(Tab_IO_Labels[IO_DTENTR], "ENDT", 4);
	break;
      case IO_TSTP:
	strncpy(Tab_IO_Labels[IO_TSTP], "TSTP", 4);
	break;
      }
}


void get_dataset_name(enum iofields blocknr, char *buf)
{

  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;
    case IO_HSML:
      strcpy(buf, "SmoothingLength");
      break;
    case IO_POT:
      strcpy(buf, "Potential");
      break;
    case IO_ACCEL:
      strcpy(buf, "Acceleration");
      break;
    case IO_DTENTR:
      strcpy(buf, "RateOfChangeOfEntropy");
      break;
    case IO_TSTP:
      strcpy(buf, "TimeStep");
      break;
    }
}


void write_file(char *fname, int writeTask, int lastTask)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[6];
  int n_for_this_task, ntask, n, p, pc, offset = 0, task;
  int blockmaxlen, ntot_type[6], nn[6];
  enum iofields blocknr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[6], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank, pcsum = 0;
  char buf[500];
#endif

#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 6; n++)
	ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
	  for(n = 0; n < 6; n++)
	    ntot_type[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }


  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      header.npartTotal[n] = (unsigned int) ntot_type_all[n];
      header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }

  for(n = 0; n < 6; n++)
    header.mass[n] = gd.MassTable[n];

  header.time = gd.Time;

  if(gd.ComovingIntegrationOn)
    header.redshift = 1.0 / gd.Time - 1;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

#ifdef COOLING
  header.flag_cooling = 1;
#endif
#ifdef SFR
  header.flag_sfr = 1;
  header.flag_feedback = 1;
#ifdef STELLARAGE
  header.flag_stellarage = 1;
#endif
#ifdef METALS
  header.flag_metals = 1;
#endif
#endif

  header.num_files = gd.NumFilesPerSnapshot;
  header.BoxSize = gd.BoxSize;
  header.Omega0 = gd.Omega0;
  header.OmegaLambda = gd.OmegaLambda;
  header.HubbleParam = gd.HubbleParam;

  if(ThisTask == writeTask)
    {
      if(gd.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  sprintf(buf, "%s.hdf5", fname);
	  hdf5_file = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	  hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gcreate(hdf5_file, buf, 0);
		}
	    }

	  write_header_attributes_in_hdf5(hdf5_headergrp);
#endif
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s' for writing snapshot.\n", fname);
	      endrun_mpi(0,123);
	    }

	  if(gd.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP;
	      my_fwrite("HEAD", sizeof(char), 4, fd);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, fd);
	      SKIP;
	    }

	  blksize = sizeof(header);
	  SKIP;
	  my_fwrite(&header, sizeof(header), 1, fd);
	  SKIP;
	}
    }

  ntask = lastTask - writeTask + 1;

  for(blocknr = 0; blocknr < IO_NBLOCKS; blocknr++)
    {
      if(blockpresent(blocknr))
	{
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr);

	  blockmaxlen = ((int) (gd.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == writeTask)
		{

		  if(gd.SnapFormat == 1 || gd.SnapFormat == 2)
		    {
		      if(gd.SnapFormat == 2)
			{
			  blksize = sizeof(int) + 4 * sizeof(char);
			  SKIP;
			  my_fwrite(Tab_IO_Labels[blocknr], sizeof(char), 4, fd);
			  nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
			  my_fwrite(&nextblock, sizeof(int), 1, fd);
			  SKIP;
			}

		      blksize = npart * bytes_per_blockelement;
		      SKIP;

		    }
		}

	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
		    {
#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && gd.SnapFormat == 3 && header.npart[type] > 0)
			{
			  switch (get_datatype_in_block(blocknr))
			    {
			    case 0:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
			      break;
			    case 1:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
			      break;
			    case 2:
			      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
			      break;
			    }

			  dims[0] = header.npart[type];
			  dims[1] = get_values_per_blockelement(blocknr);
			  if(dims[1] == 1)
			    rank = 1;
			  else
			    rank = 2;

			  get_dataset_name(blocknr, buf);

			  hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
			  hdf5_dataset =
			    H5Dcreate(hdf5_grp[type], buf, hdf5_datatype, hdf5_dataspace_in_file,
				      H5P_DEFAULT);
			  pcsum = 0;
			}
#endif

		      for(task = writeTask, offset = 0; task <= lastTask; task++)
			{
			  if(task == ThisTask)
			    {
			      n_for_this_task = n_type[type];

			      for(p = writeTask; p <= lastTask; p++)
				if(p != ThisTask)
				  MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
			    }
			  else
			    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD,
				     &status);

			  while(n_for_this_task > 0)
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == task)
				fill_write_buffer(blocknr, &offset, pc, type);

			      if(ThisTask == writeTask && task != writeTask)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					 TAG_PDATA, MPI_COMM_WORLD, &status);

			      if(ThisTask != writeTask && task == ThisTask)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
					  TAG_PDATA, MPI_COMM_WORLD);

			      if(ThisTask == writeTask)
				{
				  if(gd.SnapFormat == 3)
				    {
#ifdef HAVE_HDF5
				      start[0] = pcsum;
				      start[1] = 0;

				      count[0] = pc;
				      count[1] = get_values_per_blockelement(blocknr);
				      pcsum += pc;

				      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							  start, NULL, count, NULL);

				      dims[0] = pc;
				      dims[1] = get_values_per_blockelement(blocknr);
				      hdf5_dataspace_memory = H5Screate_simple(rank, dims, NULL);

				      hdf5_status =
					H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory,
						 hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

				      H5Sclose(hdf5_dataspace_memory);
#endif
				    }
				  else
				    my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
				}

			      n_for_this_task -= pc;
			    }
			}

#ifdef HAVE_HDF5
		      if(ThisTask == writeTask && gd.SnapFormat == 3 && header.npart[type] > 0)
			{
			  if(gd.SnapFormat == 3)
			    {
			      H5Dclose(hdf5_dataset);
			      H5Sclose(hdf5_dataspace_in_file);
			      H5Tclose(hdf5_datatype);
			    }
			}
#endif
		    }
		}

	      if(ThisTask == writeTask)
		{
		  if(gd.SnapFormat == 1 || gd.SnapFormat == 2)
		    SKIP;
		}
	    }
	}
    }

  if(ThisTask == writeTask)
    {
      if(gd.SnapFormat == 3)
	{
#ifdef HAVE_HDF5
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Gclose(hdf5_headergrp);
	  H5Fclose(hdf5_file);
#endif
	}
      else
	fclose(fd);
    }
}


#ifdef HAVE_HDF5
void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = { 6 };
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npart);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);


  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  header.flag_entropy_instead_u = 0;

  hdf5_dataspace = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL);
  hdf5_attribute = H5Acreate(handle, "Flag_Entropy_ICs", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &header.flag_entropy_instead_u);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
}
#endif


size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun_mpi(0,777);
    }
  return nwritten;
}

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun_mpi(0,778);
    }
  return nread;
}
