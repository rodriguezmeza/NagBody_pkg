
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS


#include "globaldefs.h"
#include "protodefs.h"


void read_ic(char *fname)
{
  int i, num_files, rest_files, ngroups, gr, filenr, masterTask, lastTask, groupMaster;
  double u_init;
  char buf[500];

#ifndef ISOTHERM_EQS
  double molecular_weight;
#endif
#ifdef SFR
  double original_gas_mass, mass, masstot;
#endif

  NumPart = 0;
  N_gas = 0;
  gd.TotNumPart = 0;

  num_files = find_files(fname);

  rest_files = num_files;

  fill_Tab_IO_Labels();

  while(rest_files > NTask)
    {
      sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
      if(gd.ICFormat == 3)
	sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));

      ngroups = NTask / gd.NumFilesWrittenInParallel;
      if((NTask % gd.NumFilesWrittenInParallel))
	ngroups++;
      groupMaster = (ThisTask / ngroups) * ngroups;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if(ThisTask == (groupMaster + gr))
	    read_file(buf, ThisTask, ThisTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      rest_files -= NTask;
    }


  if(rest_files > 0)
    {
      distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(num_files > 1)
	{
	  sprintf(buf, "%s.%d", fname, filenr);
	  if(gd.ICFormat == 3)
	    sprintf(buf, "%s.%d.hdf5", fname, filenr);
	}
      else
	{
	  sprintf(buf, "%s", fname);
	  if(gd.ICFormat == 3)
	    sprintf(buf, "%s.hdf5", fname);
	}

      ngroups = rest_files / gd.NumFilesWrittenInParallel;
      if((rest_files % gd.NumFilesWrittenInParallel))
	ngroups++;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if((filenr / gd.NumFilesWrittenInParallel) == gr)
	    read_file(buf, masterTask, lastTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }

  for(i = 0; i < NumPart; i++)
    {
      if(gd.MassTable[P[i].Type] != 0)
	P[i].Mass = gd.MassTable[P[i].Type];
    }

//    if(RestartFlag == 0)
    if(cmd.RestartFlag == 0)
    {
      if(gd.InitGasTemp > 0)
	{
	  u_init = (BOLTZMANN / PROTONMASS) * gd.InitGasTemp;
	  u_init *= gd.UnitMass_in_g / gd.UnitEnergy_in_cgs;

#ifdef ISOTHERM_EQS
	  u_init *= 1.0;
#else
	  u_init *= (1.0 / GAMMA_MINUS1);

	  if(gd.InitGasTemp > 1.0e4)
	    molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
	  else
	    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);

	  u_init /= molecular_weight;
#endif

	  for(i = 0; i < N_gas; i++)
	    {
	      if(SphP[i].Entropy == 0)
		SphP[i].Entropy = u_init;


	    }
	}
    }

  for(i = 0; i < N_gas; i++)
    SphP[i].Entropy = dmax(gd.MinEgySpec, SphP[i].Entropy);

  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("reading done.\n");
      fflush(stdout);
    }

  if(ThisTask == 0)
    {
      printf("Total number of particles :  %d%09d\n\n",
	     (int) (gd.TotNumPart / 1000000000), (int) (gd.TotNumPart % 1000000000));
      fflush(stdout);
    }
}


void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
  int n, k;
  float *fp;

#ifdef LONGIDS
  long long *ip;
#else
  int *ip;
#endif

  fp = CommBuffer;
  ip = CommBuffer;

  switch (blocknr)
    {
    case IO_POS:
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  P[offset + n].Pos[k] = *fp++;

      for(n = 0; n < pc; n++)
	P[offset + n].Type = type;
      break;

    case IO_VEL:
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  P[offset + n].Vel[k] = *fp++;
      break;

    case IO_ID:
      for(n = 0; n < pc; n++)
	P[offset + n].ID = *ip++;
      break;

    case IO_MASS:
      for(n = 0; n < pc; n++)
	P[offset + n].Mass = *fp++;
      break;

    case IO_U:
      for(n = 0; n < pc; n++)
	SphP[offset + n].Entropy = *fp++;
      break;

    case IO_RHO:
      for(n = 0; n < pc; n++)
	SphP[offset + n].Density = *fp++;
      break;


    case IO_HSML:
      for(n = 0; n < pc; n++)
	SphP[offset + n].Hsml = *fp++;
      break;




    case IO_POT:
    case IO_ACCEL:
    case IO_DTENTR:
    case IO_TSTP:
      break;
    }
}



void read_file(char *fname, int readTask, int lastTask)
{
  int blockmaxlen;
  int i, n_in_file, n_for_this_task, ntask, pc, offset = 0, task;
  int blksize1, blksize2;
  MPI_Status status;
  FILE *fd = 0;
  int nall;
  int type;
  char label[4];
  int nstart, bytes_per_blockelement, npart, nextblock, typelist[6];
  enum iofields blocknr;

#ifdef HAVE_HDF5
  char buf[500];
  int rank, pcsum;
  hid_t hdf5_file, hdf5_grp[6], hdf5_dataspace_in_file;
  hid_t hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataset;
  hsize_t dims[2], count[2], start[2];
#endif

#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}

  if(ThisTask == readTask)
    {
      if(gd.ICFormat == 1 || gd.ICFormat == 2)
	{
	  if(!(fd = fopen(fname, "r")))
	    {
	      printf("can't open file `%s' for reading initial conditions.\n", fname);
	      endrun_mpi(0,123);
	    }

	  if(gd.ICFormat == 2)
	    {
	      SKIP;
	      my_fread(&label, sizeof(char), 4, fd);
	      my_fread(&nextblock, sizeof(int), 1, fd);
	      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3],
		     nextblock);
	      SKIP2;
	    }

	  SKIP;
	  my_fread(&header, sizeof(header), 1, fd);
	  SKIP2;

	  if(blksize1 != 256 || blksize2 != 256)
	    {
	      printf("incorrect header format\n");
	      fflush(stdout);
	      endrun_mpi(0,890);
	    }
	}


#ifdef HAVE_HDF5
      if(gd.ICFormat == 3)
	{
	  read_header_attributes_in_hdf5(fname);

	  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

	  for(type = 0; type < 6; type++)
	    {
	      if(header.npart[type] > 0)
		{
		  sprintf(buf, "/PartType%d", type);
		  hdf5_grp[type] = H5Gopen(hdf5_file, buf);
		}
	    }
	}
#endif

      for(task = readTask + 1; task <= lastTask; task++)
	MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
    }
  else
    MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);


  if(gd.TotNumPart == 0)
    {
      if(header.num_files <= 1)
	for(i = 0; i < 6; i++)
	  header.npartTotal[i] = header.npart[i];

      gd.TotN_gas = header.npartTotal[0] + (((long long) header.npartTotalHighWord[0]) << 32);

      for(i = 0, gd.TotNumPart = 0; i < 6; i++)
	{
	  gd.TotNumPart += header.npartTotal[i];
	  gd.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
	}


      for(i = 0; i < 6; i++)
	gd.MassTable[i] = header.mass[i];

      gd.MaxPart = gd.PartAllocFactor * (gd.TotNumPart / NTask);
      gd.MaxPartSph = gd.PartAllocFactor * (gd.TotN_gas / NTask);
      allocate_memory();

//        if(RestartFlag == 2)
        if(cmd.RestartFlag == 2)
	gd.Time = gd.TimeBegin = header.time;
    }

  if(ThisTask == readTask)
    {
      for(i = 0, n_in_file = 0; i < 6; i++)
	n_in_file += header.npart[i];

      printf("\nreading file `%s' on task=%d (contains %d particles.)\n"
	     "distributing this file to tasks %d-%d\n"
	     "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 1 (halo):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 2 (disk):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 3 (bulge): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 4 (stars): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 5 (bndry): %8d  (tot=%6d%09d) masstab=%g\n\n", fname, ThisTask, n_in_file, readTask,
	     lastTask, header.npart[0], (int) (header.npartTotal[0] / 1000000000),
	     (int) (header.npartTotal[0] % 1000000000), gd.MassTable[0], header.npart[1],
	     (int) (header.npartTotal[1] / 1000000000), (int) (header.npartTotal[1] % 1000000000),
	     gd.MassTable[1], header.npart[2], (int) (header.npartTotal[2] / 1000000000),
	     (int) (header.npartTotal[2] % 1000000000), gd.MassTable[2], header.npart[3],
	     (int) (header.npartTotal[3] / 1000000000), (int) (header.npartTotal[3] % 1000000000),
	     gd.MassTable[3], header.npart[4], (int) (header.npartTotal[4] / 1000000000),
	     (int) (header.npartTotal[4] % 1000000000), gd.MassTable[4], header.npart[5],
	     (int) (header.npartTotal[5] / 1000000000), (int) (header.npartTotal[5] % 1000000000),
	     gd.MassTable[5]);
      fflush(stdout);
    }


  ntask = lastTask - readTask + 1;



  for(type = 0, nall = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;

      nall += n_for_this_task;
    }

  memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(particle_data));
  nstart = N_gas;



  for(blocknr = 0; blocknr < IO_NBLOCKS; blocknr++)
    {
      if(blockpresent(blocknr))
	{
//        if(RestartFlag == 0 && blocknr > IO_U)
        if(cmd.RestartFlag == 0 && blocknr > IO_U)
	    continue;

	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr);

	  blockmaxlen = ((int) (gd.BufferSize * 1024 * 1024)) / bytes_per_blockelement;

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == readTask)
		{
		  if(gd.ICFormat == 2)
		    {
		      SKIP;
		      my_fread(&label, sizeof(char), 4, fd);
		      my_fread(&nextblock, sizeof(int), 1, fd);
		      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2],
			     label[3], nextblock);
		      SKIP2;

		      if(strncmp(label, Tab_IO_Labels[blocknr], 4) != 0)
			{
			  printf("incorrect block-structure!\n");
			  printf("expected '%c%c%c%c' but found '%c%c%c%c'\n",
				 label[0], label[1], label[2], label[3],
				 Tab_IO_Labels[blocknr][0], Tab_IO_Labels[blocknr][1],
				 Tab_IO_Labels[blocknr][2], Tab_IO_Labels[blocknr][3]);
			  fflush(stdout);
			  endrun_mpi(0,1890);
			}
		    }

		  if(gd.ICFormat == 1 || gd.ICFormat == 2)
		    SKIP;
		}

	      for(type = 0, offset = 0; type < 6; type++)
		{
		  n_in_file = header.npart[type];
#ifdef HAVE_HDF5
		  pcsum = 0;
#endif
		  if(typelist[type] == 0)
		    {
		      n_for_this_task = n_in_file / ntask;
		      if((ThisTask - readTask) < (n_in_file % ntask))
			n_for_this_task++;

		      offset += n_for_this_task;
		    }
		  else
		    {
		      for(task = readTask; task <= lastTask; task++)
			{
			  n_for_this_task = n_in_file / ntask;
			  if((task - readTask) < (n_in_file % ntask))
			    n_for_this_task++;

			  if(task == ThisTask)
			    if(NumPart + n_for_this_task > gd.MaxPart)
			      {
				printf("too many particles\n");
				endrun_mpi(0,1313);
			      }


			  do
			    {
			      pc = n_for_this_task;

			      if(pc > blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == readTask)
				{
				  if(gd.ICFormat == 1 || gd.ICFormat == 2)
				    my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
#ifdef HAVE_HDF5
				  if(gd.ICFormat == 3)
				    {
				      get_dataset_name(blocknr, buf);
				      hdf5_dataset = H5Dopen(hdf5_grp[type], buf);

				      dims[0] = header.npart[type];
				      dims[1] = get_values_per_blockelement(blocknr);
				      if(dims[1] == 1)
					rank = 1;
				      else
					rank = 2;

				      hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);

				      dims[0] = pc;
				      hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

				      start[0] = pcsum;
				      start[1] = 0;

				      count[0] = pc;
				      count[1] = get_values_per_blockelement(blocknr);
				      pcsum += pc;

				      H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
							  start, NULL, count, NULL);

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

				      H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory,
					      hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

				      H5Tclose(hdf5_datatype);
				      H5Sclose(hdf5_dataspace_in_memory);
				      H5Sclose(hdf5_dataspace_in_file);
				      H5Dclose(hdf5_dataset);
				    }
#endif
				}

			      if(ThisTask == readTask && task != readTask)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA,
					  MPI_COMM_WORLD);

			      if(ThisTask != readTask && task == ThisTask)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask,
					 TAG_PDATA, MPI_COMM_WORLD, &status);

			      if(ThisTask == task)
				{
				  empty_read_buffer(blocknr, nstart + offset, pc, type);

				  offset += pc;
				}

			      n_for_this_task -= pc;
			    }
			  while(n_for_this_task > 0);
			}
		    }
		}
	      if(ThisTask == readTask)
		{
		  if(gd.ICFormat == 1 || gd.ICFormat == 2)
		    {
		      SKIP2;
		      if(blksize1 != blksize2)
			{
			  printf("incorrect block-sizes detected!\n");
			  printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, blocknr,
				 blksize1, blksize2);
			  fflush(stdout);
			  endrun_mpi(0,1889);
			}
		    }
		}
	    }
	}
    }


  for(type = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;

      NumPart += n_for_this_task;

      if(type == 0)
	N_gas += n_for_this_task;
    }

  if(ThisTask == readTask)
    {
      if(gd.ICFormat == 1 || gd.ICFormat == 2)
	fclose(fd);
#ifdef HAVE_HDF5
      if(gd.ICFormat == 3)
	{
	  for(type = 5; type >= 0; type--)
	    if(header.npart[type] > 0)
	      H5Gclose(hdf5_grp[type]);
	  H5Fclose(hdf5_file);
	}
#endif
    }
}




int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if(gd.ICFormat == 3)
    {
      sprintf(buf, "%s.%d.hdf5", fname, 0);
      sprintf(buf1, "%s.hdf5", fname);
    }

#ifndef  HAVE_HDF5
  if(gd.ICFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun_mpi(0,0);
    }
#endif

  header.num_files = 0;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "r")))
	{
	  if(gd.ICFormat == 1 || gd.ICFormat == 2)
	    {
	      if(gd.ICFormat == 2)
		{
		  fread(&dummy, sizeof(dummy), 1, fd);
		  fread(&dummy, sizeof(dummy), 1, fd);
		  fread(&dummy, sizeof(dummy), 1, fd);
		  fread(&dummy, sizeof(dummy), 1, fd);
		}

	      fread(&dummy, sizeof(dummy), 1, fd);
	      fread(&header, sizeof(header), 1, fd);
	      fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(gd.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf);
#endif
	}
    }

  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf1, "r")))
	{
	  if(gd.ICFormat == 1 || gd.ICFormat == 2)
	    {
	      if(gd.ICFormat == 2)
		{
		  fread(&dummy, sizeof(dummy), 1, fd);
		  fread(&dummy, sizeof(dummy), 1, fd);
		  fread(&dummy, sizeof(dummy), 1, fd);
		  fread(&dummy, sizeof(dummy), 1, fd);
		}

	      fread(&dummy, sizeof(dummy), 1, fd);
	      fread(&header, sizeof(header), 1, fd);
	      fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

#ifdef HAVE_HDF5
	  if(gd.ICFormat == 3)
	    read_header_attributes_in_hdf5(buf1);
#endif
	  header.num_files = 1;
	}
    }

  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      printf("\nCan't find initial conditions file.");
      printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
      fflush(stdout);
    }

  endrun_mpi(0,0);
  return 0;
}



void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last)
{
  int ntask, filesleft, filesright, tasksleft, tasksright;

  if(nfiles > 1)
    {
      ntask = lasttask - firsttask + 1;

      filesleft = (((double) (ntask / 2)) / ntask) * nfiles;
      if(filesleft <= 0)
	filesleft = 1;
      if(filesleft >= nfiles)
	filesleft = nfiles - 1;

      filesright = nfiles - filesleft;

      tasksleft = ntask / 2;
      tasksright = ntask - tasksleft;

      distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
      distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr, master,
		      last);
    }
  else
    {
      if(ThisTask >= firsttask && ThisTask <= lasttask)
	{
	  *filenr = firstfile;
	  *master = firsttask;
	  *last = lasttask;
	}
    }
}


#ifdef HAVE_HDF5
void read_header_attributes_in_hdf5(char *fname)
{
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;


  hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen(hdf5_file, "/Header");


  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Entropy_ICs");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_entropy_instead_u);
  H5Aclose(hdf5_attribute);

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);
}
#endif
