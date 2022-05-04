
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS


#include "globaldefs.h"
#include "protodefs.h"


#ifdef PMGRID

void long_range_init(void)
{
/* #ifdef PERIODIC
  pm_init_periodic();
#ifdef PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif
#else */
  pm_init_nonperiodic();
//#endif
}


void long_range_init_regionsize(void)
{
/* #ifdef PERIODIC
#ifdef PLACEHIGHRESREGION
//    if(RestartFlag != 1)
    if(cmd.RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
#else */
//    if(RestartFlag != 1)
    if(cmd.RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
//#endif
}


void long_range_force(void)
{
  int i;

//#ifndef PERIODIC
  int j;
  double fac;
//#endif


  for(i = 0; i < NumPart; i++)
    P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

#ifdef NOGRAVITY
  return;
#endif


/* #ifdef PERIODIC
  pmforce_periodic();
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
  if(i == 1)
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1);
    }
  if(i == 1)
    endrun_mpi(0,68686);
#endif
#else */
  i = pmforce_nonperiodic(0);
  if(i == 1)
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(0);
    }
  if(i == 1)
    endrun_mpi(0,68687);
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
  if(i == 1)
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

      for(i = 0; i < NumPart; i++)
	P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

      i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);
    }
  if(i != 0)
    endrun_mpi(0,68688);
#endif
//#endif


//#ifndef PERIODIC
  if(gd.ComovingIntegrationOn)
    {
      fac = 0.5 * gd.Hubble * gd.Hubble * gd.Omega0;

      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }

  if(gd.ComovingIntegrationOn == 0)
    {
      fac = gd.OmegaLambda * gd.Hubble * gd.Hubble;

      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }
//#endif

}


#endif
