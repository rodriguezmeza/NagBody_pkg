
// DEJAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <mpi.h>

//#include <fcntl.h>
//#include <sys/stat.h>
//#include <sys/types.h>
//#include <sys/file.h>
//#include <unistd.h>
//#include <gsl/gsl_rng.h>

#include "globaldefs.h"
#include "protodefs.h"

//#include "../../../General_libs/mpi/mpi_proto.h"


static FILE *fd;

static void in(int *x, int modus);
static void byten(void *x, size_t n, int modus);


void restart(int modus)
{
  char buf[200], buf_bak[200], buf_mv[500];
  double save_PartAllocFactor, save_TreeAllocFactor;
  int i, nprocgroup, masterTask, groupTask, old_MaxPart, old_MaxNodes;
  global_data_all_processes all_task0;


  sprintf(buf, "%s%s.%d", gd.OutputDir, gd.RestartFile, ThisTask);
  sprintf(buf_bak, "%s%s.%d.bak", gd.OutputDir, gd.RestartFile, ThisTask);
  sprintf(buf_mv, "mv %s %s", buf, buf_bak);


  if((NTask < gd.NumFilesWrittenInParallel))
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun_mpi(0,2131);
    }

  nprocgroup = NTask / gd.NumFilesWrittenInParallel;

  if((NTask % gd.NumFilesWrittenInParallel))
    {
      nprocgroup++;
    }

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))
	{
	  if(modus)
	    {
	      if(!(fd = fopen(buf, "r")))
		{
		  printf("Restart file '%s' not found.\n", buf);
		  endrun_mpi(0,7870);
		}
	    }
	  else
	    {
	      system(buf_mv);

	      if(!(fd = fopen(buf, "w")))
		{
		  printf("Restart file '%s' cannot be opened.\n", buf);
		  endrun_mpi(0,7878);
		}
	    }


	  save_PartAllocFactor = gd.PartAllocFactor;
	  save_TreeAllocFactor = gd.TreeAllocFactor;

	  byten(&gd, sizeof(global_data_all_processes), modus);

	  if(ThisTask == 0 && modus > 0)
	    all_task0 = gd;

	  if(modus > 0 && groupTask == 0)
	    {
	      MPI_Bcast(&all_task0, sizeof(global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }

	  old_MaxPart = gd.MaxPart;
	  old_MaxNodes = gd.TreeAllocFactor * gd.MaxPart;

	  if(modus)
	    {
	      if(gd.PartAllocFactor != save_PartAllocFactor)
		{
		  gd.PartAllocFactor = save_PartAllocFactor;
		  gd.MaxPart = gd.PartAllocFactor * (gd.TotNumPart / NTask);
		  gd.MaxPartSph = gd.PartAllocFactor * (gd.TotN_gas / NTask);
		  save_PartAllocFactor = -1;
		}

	      if(gd.TreeAllocFactor != save_TreeAllocFactor)
		{
		  gd.TreeAllocFactor = save_TreeAllocFactor;
		  save_TreeAllocFactor = -1;
		}

	      if(all_task0.Time != gd.Time)
		{
		  printf("The restart file on task=%d is not consistent with the one on task=0\n", ThisTask);
		  fflush(stdout);
		  endrun_mpi(0,16);
		}

	      allocate_memory();
	    }

	  in(&NumPart, modus);

	  if(NumPart > gd.MaxPart)
	    {
	      printf
		("it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		 NumPart / (((double) gd.TotNumPart) / NTask));
	      printf("fatal error\n");
	      endrun_mpi(0,22);
	    }

	  byten(&P[0], NumPart * sizeof(particle_data), modus);

	  in(&N_gas, modus);

	  if(N_gas > 0)
	    {
	      if(N_gas > gd.MaxPartSph)
		{
		  printf
		    ("SPH: it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		     N_gas / (((double) gd.TotN_gas) / NTask));
		  printf("fatal error\n");
		  endrun_mpi(0,222);
		}
	      byten(&SphP[0], N_gas * sizeof(sph_particle_data), modus);
	    }

	  byten(gsl_rng_state(random_generator), gsl_rng_size(random_generator), modus);


	  if(modus)
	    {
	      ngb_treeallocate(MAX_NGB);

	      force_treeallocate(gd.TreeAllocFactor * gd.MaxPart, gd.MaxPart);
	    }


	  in(&Numnodestree, modus);
	  in(&NTopleaves, modus);

	  if(Numnodestree > MaxNodes)
	    {
	      printf
		("Tree storage: it seems you have reduced(!) 'PartAllocFactor' below the value needed to load the restart file (task=%d). "
		 "Numnodestree=%d  MaxNodes=%d\n", ThisTask, Numnodestree, MaxNodes);
	      endrun_mpi(0,221);
	    }

	  byten(Nodes_base, Numnodestree * sizeof(NODE), modus);
	  byten(Extnodes_base, Numnodestree * sizeof(extNODE), modus);

	  byten(Father, NumPart * sizeof(int), modus);

	  byten(Nextnode, NumPart * sizeof(int), modus);
	  byten(Nextnode + gd.MaxPart, MAXTOPNODES * sizeof(int), modus);

	  byten(DomainStartList, NTask * sizeof(int), modus);
	  byten(DomainEndList, NTask * sizeof(int), modus);
	  byten(DomainTask, MAXTOPNODES * sizeof(int), modus);
	  byten(DomainNodeIndex, MAXTOPNODES * sizeof(int), modus);
	  byten(DomainTreeNodeLen, MAXTOPNODES * sizeof(FLOAT), modus);
	  byten(DomainHmax, MAXTOPNODES * sizeof(FLOAT), modus);
	  byten(DomainMoment, MAXTOPNODES * sizeof(DomainNODE), modus);

	  byten(DomainCorner, 3 * sizeof(double), modus);
	  byten(DomainCenter, 3 * sizeof(double), modus);
	  byten(&DomainLen, sizeof(double), modus);
	  byten(&DomainFac, sizeof(double), modus);
	  byten(&DomainMyStart, sizeof(int), modus);
	  byten(&DomainMyLast, sizeof(int), modus);

	  if(modus)
	    if(gd.PartAllocFactor != save_PartAllocFactor || gd.TreeAllocFactor != save_TreeAllocFactor)
	      {
		for(i = 0; i < NumPart; i++)
		  Father[i] += (gd.MaxPart - old_MaxPart);

		for(i = 0; i < NumPart; i++)
		  if(Nextnode[i] >= old_MaxPart)
		    {
		      if(Nextnode[i] >= old_MaxPart + old_MaxNodes)
			Nextnode[i] += (gd.MaxPart - old_MaxPart) + (MaxNodes - old_MaxPart);
		      else
			Nextnode[i] += (gd.MaxPart - old_MaxPart);
		    }

		for(i = 0; i < Numnodestree; i++)
		  {
		    if(Nodes_base[i].u.d.sibling >= old_MaxPart)
		      {
			if(Nodes_base[i].u.d.sibling >= old_MaxPart + old_MaxNodes)
			  Nodes_base[i].u.d.sibling +=
			    (gd.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
			else
			  Nodes_base[i].u.d.sibling += (gd.MaxPart - old_MaxPart);
		      }

		    if(Nodes_base[i].u.d.father >= old_MaxPart)
		      {
			if(Nodes_base[i].u.d.father >= old_MaxPart + old_MaxNodes)
			  Nodes_base[i].u.d.father += (gd.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
			else
			  Nodes_base[i].u.d.father += (gd.MaxPart - old_MaxPart);
		      }

		    if(Nodes_base[i].u.d.nextnode >= old_MaxPart)
		      {
			if(Nodes_base[i].u.d.nextnode >= old_MaxPart + old_MaxNodes)
			  Nodes_base[i].u.d.nextnode +=
			    (gd.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
			else
			  Nodes_base[i].u.d.nextnode += (gd.MaxPart - old_MaxPart);
		      }
		  }

		for(i = 0; i < MAXTOPNODES; i++)
		  if(Nextnode[i + gd.MaxPart] >= old_MaxPart)
		    {
		      if(Nextnode[i + gd.MaxPart] >= old_MaxPart + old_MaxNodes)
			Nextnode[i + gd.MaxPart] += (gd.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
		      else
			Nextnode[i + gd.MaxPart] += (gd.MaxPart - old_MaxPart);
		    }

		for(i = 0; i < MAXTOPNODES; i++)
		  if(DomainNodeIndex[i] >= old_MaxPart)
		    {
		      if(DomainNodeIndex[i] >= old_MaxPart + old_MaxNodes)
			DomainNodeIndex[i] += (gd.MaxPart - old_MaxPart) + (MaxNodes - old_MaxNodes);
		      else
			DomainNodeIndex[i] += (gd.MaxPart - old_MaxPart);
		    }
	      }

	  fclose(fd);
	}
      else
	{
	  if(modus > 0 && groupTask == 0)
	    {
	      MPI_Bcast(&all_task0, sizeof(global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }
}



void byten(void *x, size_t n, int modus)
{
  if(modus)
    my_fread(x, n, 1, fd);
  else
    my_fwrite(x, n, 1, fd);
}


void in(int *x, int modus)
{
  if(modus)
    my_fread(x, 1, sizeof(int), fd);
  else
    my_fwrite(x, 1, sizeof(int), fd);
}
