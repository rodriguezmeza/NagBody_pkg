
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS


#include "globaldefs.h"
#include "protodefs.h"


/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over single
 *  timesteps. The loop terminates when the cpu-time limit is reached, when a
 *  `stop' file is found in the output directory, or when the simulation ends
 *  because we arrived at TimeMax.
 */
void run(void)
{
  FILE *fd;
  int stopflag = 0;
  char stopfname[200], contfname[200];
  double t0, t1;


  sprintf(stopfname, "%sstop", gd.OutputDir);
  sprintf(contfname, "%scont", gd.OutputDir);
  unlink(contfname);

  do				// main loop
    {
      t0 = second();

      find_next_sync_point_and_drift();	// find next synchronization point and drift particles to this time.
					 // If needed, this function will also write an output file
					 // at the desired time.


      every_timestep_stuff();	// write some info to log-files


      domain_Decomposition();	// do domain decomposition if needed


      compute_accelerations(0);	// compute accelerations for
				 // the particles that are to be advanced


      // check whether we want a full energy statistics
      if((gd.Time - gd.TimeLastStatistics) >= gd.TimeBetStatistics)
	{
#ifdef COMPUTE_POTENTIAL_ENERGY
	  compute_potential();
#endif
	  energy_statistics();	// compute and output energy statistics
	  gd.TimeLastStatistics += gd.TimeBetStatistics;
	}

      advance_and_find_timesteps();	// 'kick' active particles in
					 // momentum space and compute new
					 // timesteps for them

      gd.NumCurrentTiStep++;

      // Check whether we need to interrupt the run
      if(ThisTask == 0)
	{
	  // Is the stop-file present? If yes, interrupt the run.
	  if((fd = fopen(stopfname, "r")))
	    {
	      fclose(fd);
	      stopflag = 1;
	      unlink(stopfname);
	    }

	  // are we running out of CPU-time ? If yes, interrupt run.
	  if(CPUThisRun > 0.85 * gd.TimeLimitCPU)
	    {
	      printf("reaching time-limit. stopping.\n");
	      stopflag = 2;
	    }
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag)
	{
	  restart(0);		// write restart file
	  MPI_Barrier(MPI_COMM_WORLD);

	  if(stopflag == 2 && ThisTask == 0)
	    {
	      if((fd = fopen(contfname, "w")))
		fclose(fd);
	    }

	  if(stopflag == 2 && gd.ResubmitOn && ThisTask == 0)
	    {
	      close_outputfiles();
	      system(gd.ResubmitCommand);
	    }
	  return;
	}

      // is it time to write a regular restart-file? (for security)
      if(ThisTask == 0)
	{
	  if((CPUThisRun - gd.TimeLastRestartFile) >= gd.CpuTimeBetRestartFile)
	    {
	      gd.TimeLastRestartFile = CPUThisRun;
	      stopflag = 3;
	    }
	  else
	    stopflag = 0;
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag == 3)
	{
	  restart(0);		// write an occasional restart file
	  stopflag = 0;
	}

      t1 = second();

      gd.CPU_Total += timediff(t0, t1);
      CPUThisRun += timediff(t0, t1);
    }
  while(gd.Ti_Current < TIMEBASE && gd.Time <= gd.TimeMax);

  restart(0);

  savepositions(gd.SnapshotFileCount++);	// write a last snapshot
						 // file at final time (will
						 // be overwritten if
						 // gd.TimeMax is increased
						 // and the run is continued)

}


/*! This function finds the next synchronization point of the system (i.e. the
 *  earliest point of time any of the particles needs a force computation),
 *  and drifts the system to this point of time.  If the system drifts over
 *  the desired time of a snapshot file, the function will drift to this
 *  moment, generate an output, and then resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, min, min_glob, flag, *temp;
  double timeold;
  double t0, t1;

  t0 = second();

  timeold = gd.Time;

  for(n = 1, min = P[0].Ti_endstep; n < NumPart; n++)
    if(min > P[n].Ti_endstep)
      min = P[n].Ti_endstep;

  MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  // We check whether this is a full step where all particles are synchronized
  flag = 1;
  for(n = 0; n < NumPart; n++)
    if(P[n].Ti_endstep > min_glob)
      flag = 0;

  MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

#ifdef PMGRID
  if(min_glob >= gd.PM_Ti_endstep)
    {
      min_glob = gd.PM_Ti_endstep;
      Flag_FullStep = 1;
    }
#endif

  // Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active
  for(n = 0, NumForceUpdate = 0; n < NumPart; n++)
    {
      if(P[n].Ti_endstep == min_glob)
#ifdef SELECTIVE_NO_GRAVITY
        if(!((1 << P[n].Type) & (SELECTIVE_NO_GRAVITY)))
#endif
          NumForceUpdate++;
    }

  // note: NumForcesSinceLastDomainDecomp has type "long long"
  temp = malloc(NTask * sizeof(int));
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
  for(n = 0; n < NTask; n++)
    gd.NumForcesSinceLastDomainDecomp += temp[n];
  free(temp);



  t1 = second();

  gd.CPU_Predict += timediff(t0, t1);

  while(min_glob >= gd.Ti_nextoutput && gd.Ti_nextoutput >= 0)
    {
      move_particles(gd.Ti_Current, gd.Ti_nextoutput);

      gd.Ti_Current = gd.Ti_nextoutput;

      if(gd.ComovingIntegrationOn)
	gd.Time = gd.TimeBegin * exp(gd.Ti_Current * gd.Timebase_interval);
      else
	gd.Time = gd.TimeBegin + gd.Ti_Current * gd.Timebase_interval;

#ifdef OUTPUTPOTENTIAL
      gd.NumForcesSinceLastDomainDecomp = 1 + gd.TotNumPart * gd.TreeDomainUpdateFrequency;
      domain_Decomposition();
      compute_potential();
#endif
      savepositions(gd.SnapshotFileCount++);	// write snapshot file

      gd.Ti_nextoutput = find_next_outputtime(gd.Ti_nextoutput + 1);
    }

  move_particles(gd.Ti_Current, min_glob);

  gd.Ti_Current = min_glob;

  if(gd.ComovingIntegrationOn)
    gd.Time = gd.TimeBegin * exp(gd.Ti_Current * gd.Timebase_interval);
  else
    gd.Time = gd.TimeBegin + gd.Ti_Current * gd.Timebase_interval;

  gd.TimeStep = gd.Time - timeold;
}



/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
int find_next_outputtime(int ti_curr)
{
  int i, ti, ti_next, iter = 0;
  double next, time;

  ti_next = -1;

  if(gd.OutputListOn)
    {
      for(i = 0; i < gd.OutputListLength; i++)
	{
	  time = gd.OutputListTimes[i];

	  if(time >= gd.TimeBegin && time <= gd.TimeMax)
	    {
	      if(gd.ComovingIntegrationOn)
		ti = log(time / gd.TimeBegin) / gd.Timebase_interval;
	      else
		ti = (time - gd.TimeBegin) / gd.Timebase_interval;

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    ti_next = ti;

		  if(ti_next > ti)
		    ti_next = ti;
		}
	    }
	}
    }
  else
    {
      if(gd.ComovingIntegrationOn)
	{
	  if(gd.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun_mpi(0,13123);
	    }
	}
      else
	{
	  if(gd.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun_mpi(0,13123);
	    }
	}

      time = gd.TimeOfFirstSnapshot;

      iter = 0;

      while(time < gd.TimeBegin)
	{
	  if(gd.ComovingIntegrationOn)
	    time *= gd.TimeBetSnapshot;
	  else
	    time += gd.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun_mpi(0,110);
	    }
	}

      while(time <= gd.TimeMax)
	{
	  if(gd.ComovingIntegrationOn)
	    ti = log(time / gd.TimeBegin) / gd.Timebase_interval;
	  else
	    ti = (time - gd.TimeBegin) / gd.Timebase_interval;

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(gd.ComovingIntegrationOn)
	    time *= gd.TimeBetSnapshot;
	  else
	    time += gd.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun_mpi(0,111);
	    }
	}
    }

  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	// this will prevent any further output

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(gd.ComovingIntegrationOn)
	next = gd.TimeBegin * exp(ti_next * gd.Timebase_interval);
      else
	next = gd.TimeBegin + ti_next * gd.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g\n\n", next);
    }

  return ti_next;
}




/*! This routine writes one line for every timestep to two log-files.  In
 *  FdInfo, we just list the timesteps that have been done, while in FdCPU the
 *  cumulative cpu-time consumption in various parts of the code is stored.
 */
void every_timestep_stuff(void)
{
  double z;

  if(ThisTask == 0)
    {
      if(gd.ComovingIntegrationOn)
	{
	  z = 1.0 / (gd.Time) - 1;
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
		  gd.NumCurrentTiStep, gd.Time, z, gd.TimeStep,
		  log(gd.Time) - log(gd.Time - gd.TimeStep));
	  printf("\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", gd.NumCurrentTiStep,
		 gd.Time, z, gd.TimeStep, log(gd.Time) - log(gd.Time - gd.TimeStep));
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nBegin Step %d, Time: %g, Systemstep: %g\n", gd.NumCurrentTiStep, gd.Time,
		  gd.TimeStep);
	  printf("\nBegin Step %d, Time: %g, Systemstep: %g\n", gd.NumCurrentTiStep, gd.Time, gd.TimeStep);
	  fflush(FdInfo);
	}

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", gd.NumCurrentTiStep, gd.Time, NTask);

      fprintf(FdCPU,
	      "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	      gd.CPU_Total, gd.CPU_Gravity, gd.CPU_Hydro, gd.CPU_Domain, gd.CPU_Potential,
	      gd.CPU_Predict, gd.CPU_TimeLine, gd.CPU_Snapshot, gd.CPU_TreeWalk, gd.CPU_TreeConstruction,
	      gd.CPU_CommSum, gd.CPU_Imbalance, gd.CPU_HydCompWalk, gd.CPU_HydCommSumm,
	      gd.CPU_HydImbalance, gd.CPU_EnsureNgb, gd.CPU_PM, gd.CPU_Peano);
      fflush(FdCPU);
    }

  set_random_numbers();
}


/*! This routine first calls a computation of various global quantities of the
 *  particle distribution, and then writes some statistics about the energies
 *  in the various particle components to the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      gd.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

      fflush(FdEnergy);
    }
}
