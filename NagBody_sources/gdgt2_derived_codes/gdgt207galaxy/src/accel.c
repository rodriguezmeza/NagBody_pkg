
#include "globaldefs.h"
#include "protodefs.h"


void compute_accelerations(int mode)
{
  double tstart, tend;

  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

#ifdef PMGRID
  if(gd.PM_Ti_endstep == gd.Ti_Current)
    {
      tstart = second();
      long_range_force();
      tend = second();
      gd.CPU_PM += timediff(tstart, tend);
    }
#endif

  tstart = second();

  gravity_tree();

  if(gd.TypeOfOpeningCriterion == 1 && gd.Ti_Current == 0)
    gravity_tree();
  tend = second();
  gd.CPU_Gravity += timediff(tstart, tend);

#ifdef FORCETEST
  gravity_forcetest();
#endif

  if(gd.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}

      tstart = second();
      density();
      tend = second();
      gd.CPU_Hydro += timediff(tstart, tend);

      tstart = second();
      force_update_hmax();
      tend = second();
      gd.CPU_Predict += timediff(tstart, tend);


      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}

      tstart = second();
      hydro_force();
      tend = second();
      gd.CPU_Hydro += timediff(tstart, tend);
    }

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}
