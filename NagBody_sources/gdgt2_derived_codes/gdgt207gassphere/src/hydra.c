
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS

#include "globaldefs.h"
#include "protodefs.h"


static double hubble_a, atime, hubble_a2, fac_mu, fac_vsic_fix, a3inv, fac_egy;

/* #ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif */


void hydro_force(void)
{
  long long ntot, ntotleft;
  int i, j, k, n, ngrp, maxfill, source, ndone;
  int *nbuffer, *noffset, *nsend_local, *nsend, *numlist, *ndonelist;
  int level, sendTask, recvTask, nexport, place;
  double soundspeed_i;
  double tstart, tend, sumt, sumcomm;
  double timecomp = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  MPI_Status status;

/* #ifdef PERIODIC
  boxSize = gd.BoxSize;
  boxHalf = 0.5 * gd.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif */

  if(gd.ComovingIntegrationOn)
    {
      hubble_a = gd.Omega0 / (gd.Time * gd.Time * gd.Time)
	+ (1 - gd.Omega0 - gd.OmegaLambda) / (gd.Time * gd.Time) + gd.OmegaLambda;

      hubble_a = gd.Hubble * sqrt(hubble_a);
      hubble_a2 = gd.Time * gd.Time * hubble_a;

      fac_mu = pow(gd.Time, 3 * (GAMMA - 1) / 2) / gd.Time;

      fac_egy = pow(gd.Time, 3 * (GAMMA - 1));

      fac_vsic_fix = hubble_a * pow(gd.Time, 3 * GAMMA_MINUS1);

      a3inv = 1 / (gd.Time * gd.Time * gd.Time);
      atime = gd.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = fac_vsic_fix = a3inv = fac_egy = 1.0;

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].Ti_endstep == gd.Ti_Current)
	NumSphUpdate++;
    }

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
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

      tstart = second();
      for(nexport = 0, ndone = 0; i < N_gas && nexport < gd.BunchSizeHydro - NTask; i++)
	if(P[i].Ti_endstep == gd.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;

	    hydro_evaluate(i, 0);

	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      {
			HydroDataIn[nexport].Pos[k] = P[i].Pos[k];
			HydroDataIn[nexport].Vel[k] = SphP[i].VelPred[k];
		      }
		    HydroDataIn[nexport].Hsml = SphP[i].Hsml;
		    HydroDataIn[nexport].Mass = P[i].Mass;
		    HydroDataIn[nexport].DhsmlDensityFactor = SphP[i].DhsmlDensityFactor;
		    HydroDataIn[nexport].Density = SphP[i].Density;
		    HydroDataIn[nexport].Pressure = SphP[i].Pressure;
		    HydroDataIn[nexport].Timestep = P[i].Ti_endstep - P[i].Ti_begstep;

		    soundspeed_i = sqrt(GAMMA * SphP[i].Pressure / SphP[i].Density);
		    HydroDataIn[nexport].F1 = fabs(SphP[i].DivVel) /
		      (fabs(SphP[i].DivVel) + SphP[i].CurlVel +
		       0.0001 * soundspeed_i / SphP[i].Hsml / fac_mu);

		    HydroDataIn[nexport].Index = i;
		    HydroDataIn[nexport].Task = j;
		    nexport++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timecomp += timediff(tstart, tend);

      qsort(HydroDataIn, nexport, sizeof(hydrodata_in), hydro_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);




      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
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
	      if(maxfill >= gd.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      MPI_Sendrecv(&HydroDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A,
				   &HydroDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(hydrodata_in), MPI_BYTE,
				   recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    hydro_evaluate(j, 1);
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  tstart = second();
	  MPI_Barrier(MPI_COMM_WORLD);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  tstart = second();
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
	      if(maxfill >= gd.BunchSizeHydro)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      MPI_Sendrecv(&HydroDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B,
				   &HydroDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(hydrodata_out),
				   MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, &status);

		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = HydroDataIn[source].Index;

			  for(k = 0; k < 3; k++)
			    SphP[place].HydroAccel[k] += HydroDataPartialResult[source].Acc[k];

			  SphP[place].DtEntropy += HydroDataPartialResult[source].DtEntropy;

			  if(SphP[place].MaxSignalVel < HydroDataPartialResult[source].MaxSignalVel)
			    SphP[place].MaxSignalVel = HydroDataPartialResult[source].MaxSignalVel;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

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



  tstart = second();

  for(i = 0; i < N_gas; i++)
    if(P[i].Ti_endstep == gd.Ti_Current)
      {
	SphP[i].DtEntropy *= GAMMA_MINUS1 / (hubble_a2 * pow(SphP[i].Density, GAMMA_MINUS1));
#ifdef SPH_BND_PARTICLES
	if(P[i].ID == 0)
	  {
	    SphP[i].DtEntropy = 0;
	    for(k = 0; k < 3; k++)
	      SphP[i].HydroAccel[k] = 0;
	  }
#endif
      }

  tend = second();
  timecomp += timediff(tstart, tend);

  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      gd.CPU_HydCompWalk += sumt / NTask;
      gd.CPU_HydCommSumm += sumcomm / NTask;
      gd.CPU_HydImbalance += sumimbalance / NTask;
    }
}


void hydro_evaluate(int target, int mode)
{
  int j, k, n, timestep, startnode, numngb;
  FLOAT *pos, *vel;
  FLOAT mass, h_i, dhsmlDensityFactor, rho, pressure, f1, f2;
  double acc[3], dtEntropy, maxSignalVel;
  double dx, dy, dz, dvx, dvy, dvz;
  double h_i2, hinv, hinv4;
  double p_over_rho2_i, p_over_rho2_j, soundspeed_i, soundspeed_j;
  double hfc, dwk_i, vdotr, vdotr2, visc, mu_ij, rho_ij, vsig;
  double h_j, dwk_j, r, r2, u, hfc_visc;

#ifndef NOVISCOSITYLIMITER
  double dt;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = SphP[target].VelPred;
      h_i = SphP[target].Hsml;
      mass = P[target].Mass;
      dhsmlDensityFactor = SphP[target].DhsmlDensityFactor;
      rho = SphP[target].Density;
      pressure = SphP[target].Pressure;
      timestep = P[target].Ti_endstep - P[target].Ti_begstep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
      f1 = fabs(SphP[target].DivVel) /
	(fabs(SphP[target].DivVel) + SphP[target].CurlVel +
	 0.0001 * soundspeed_i / SphP[target].Hsml / fac_mu);
    }
  else
    {
      pos = HydroDataGet[target].Pos;
      vel = HydroDataGet[target].Vel;
      h_i = HydroDataGet[target].Hsml;
      mass = HydroDataGet[target].Mass;
      dhsmlDensityFactor = HydroDataGet[target].DhsmlDensityFactor;
      rho = HydroDataGet[target].Density;
      pressure = HydroDataGet[target].Pressure;
      timestep = HydroDataGet[target].Timestep;
      soundspeed_i = sqrt(GAMMA * pressure / rho);
      f1 = HydroDataGet[target].F1;
    }


  acc[0] = acc[1] = acc[2] = dtEntropy = 0;
  maxSignalVel = 0;

  p_over_rho2_i = pressure / (rho * rho) * dhsmlDensityFactor;
  h_i2 = h_i * h_i;

  startnode = gd.MaxPart;
  do
    {
      numngb = ngb_treefind_pairs(&pos[0], h_i, &startnode);

      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

/* #ifdef PERIODIC
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif */
	  r2 = dx * dx + dy * dy + dz * dz;
	  h_j = SphP[j].Hsml;
	  if(r2 < h_i2 || r2 < h_j * h_j)
	    {
	      r = sqrt(r2);
	      if(r > 0)
		{
		  p_over_rho2_j = SphP[j].Pressure / (SphP[j].Density * SphP[j].Density);
		  soundspeed_j = sqrt(GAMMA * p_over_rho2_j * SphP[j].Density);
		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];
		  vdotr = dx * dvx + dy * dvy + dz * dvz;

		  if(gd.ComovingIntegrationOn)
		    vdotr2 = vdotr + hubble_a2 * r2;
		  else
		    vdotr2 = vdotr;

		  if(r2 < h_i2)
		    {
		      hinv = 1.0 / h_i;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }
		  else
		    {
		      dwk_i = 0;
		    }

		  if(r2 < h_j * h_j)
		    {
		      hinv = 1.0 / h_j;
#ifndef  TWODIMS
		      hinv4 = hinv * hinv * hinv * hinv;
#else
		      hinv4 = hinv * hinv * hinv / boxSize_Z;
#endif
		      u = r * hinv;
		      if(u < 0.5)
			dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		      else
			dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		    }
		  else
		    {
		      dwk_j = 0;
		    }

		  if(soundspeed_i + soundspeed_j > maxSignalVel)
		    maxSignalVel = soundspeed_i + soundspeed_j;

		  if(vdotr2 < 0)
		    {
		      mu_ij = fac_mu * vdotr2 / r;

                      vsig = soundspeed_i + soundspeed_j - 3 * mu_ij;

		      if(vsig > maxSignalVel)
			maxSignalVel = vsig;

		      rho_ij = 0.5 * (rho + SphP[j].Density);
		      f2 =
			fabs(SphP[j].DivVel) / (fabs(SphP[j].DivVel) + SphP[j].CurlVel +
						0.0001 * soundspeed_j / fac_mu / SphP[j].Hsml);

		      visc = 0.25 * gd.ArtBulkViscConst * vsig * (-mu_ij) / rho_ij * (f1 + f2);

#ifndef NOVISCOSITYLIMITER
		      dt = imax(timestep, (P[j].Ti_endstep - P[j].Ti_begstep)) * gd.Timebase_interval;
		      if(dt > 0 && (dwk_i + dwk_j) < 0)
			{
			  visc = dmin(visc, 0.5 * fac_vsic_fix * vdotr2 /
				      (0.5 * (mass + P[j].Mass) * (dwk_i + dwk_j) * r * dt));
			}
#endif
		    }
		  else
		    visc = 0;

		  p_over_rho2_j *= SphP[j].DhsmlDensityFactor;

		  hfc_visc = 0.5 * P[j].Mass * visc * (dwk_i + dwk_j) / r;

		  hfc = hfc_visc + P[j].Mass * (p_over_rho2_i * dwk_i + p_over_rho2_j * dwk_j) / r;

		  acc[0] -= hfc * dx;
		  acc[1] -= hfc * dy;
		  acc[2] -= hfc * dz;
		  dtEntropy += 0.5 * hfc_visc * vdotr2;
		}
	    }
	}
    }
  while(startnode >= 0);

  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
	SphP[target].HydroAccel[k] = acc[k];
      SphP[target].DtEntropy = dtEntropy;
      SphP[target].MaxSignalVel = maxSignalVel;
    }
  else
    {
      for(k = 0; k < 3; k++)
	HydroDataResult[target].Acc[k] = acc[k];
      HydroDataResult[target].DtEntropy = dtEntropy;
      HydroDataResult[target].MaxSignalVel = maxSignalVel;
    }
}


int hydro_compare_key(const void *a, const void *b)
{
  if(((hydrodata_in *) a)->Task < (((hydrodata_in *) b)->Task))
    return -1;
  if(((hydrodata_in *) a)->Task > (((hydrodata_in *) b)->Task))
    return +1;
  return 0;
}
