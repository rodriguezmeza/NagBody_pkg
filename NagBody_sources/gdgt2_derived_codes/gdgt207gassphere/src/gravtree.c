
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS


#include "globaldefs.h"
#include "protodefs.h"


void gravity_tree(void)
{
  long long ntot;
  int numnodes, nexportsum = 0;
  int i, j, iter = 0;
  int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double ewaldcount;
  double costtotal, ewaldtot, *costtreelist, *ewaldlist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm;

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  long long ntotleft;
  int ndone, maxfill, ngrp;
  int k, place;
  int level, sendTask, recvTask;
  double ax, ay, az;
  MPI_Status status;
#endif

  if(gd.ComovingIntegrationOn)
    set_softenings();


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

  costtotal = ewaldcount = 0;

  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


#ifndef NOGRAVITY
  if(ThisTask == 0)
    printf("Begin tree force.\n");


#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (SELECTIVE_NO_GRAVITY)))
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif


  noffset = malloc(sizeof(int) * NTask);
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;
  ntotleft = ntot;

  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      tstart = second();
      for(nexport = 0, ndone = 0; i < NumPart && nexport < gd.BunchSizeForce - NTask; i++)
	if(P[i].Ti_endstep == gd.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;
#ifndef PMGRID
	    costtotal += force_treeevaluate(i, 0, &ewaldcount);
#else
	    costtotal += force_treeevaluate_shortrange(i, 0);
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
		    nexportsum++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

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
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
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
	    {
#ifndef PMGRID
	      costtotal += force_treeevaluate(j, 1, &ewaldcount);
#else
	      costtotal += force_treeevaluate_shortrange(j, 1);
#endif
	    }
	  tend = second();
	  timetree += timediff(tstart, tend);

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
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];

			  P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
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


//#ifndef PERIODIC
#ifndef PMGRID
  if(gd.ComovingIntegrationOn)
    {
      fac = 0.5 * gd.Hubble * gd.Hubble * gd.Omega0 / gd.G;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == gd.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
//#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == gd.Ti_Current)
      {
#ifdef PMGRID
	ax = P[i].GravAccel[0] + P[i].GravPM[0] / gd.G;
	ay = P[i].GravAccel[1] + P[i].GravPM[1] / gd.G;
	az = P[i].GravAccel[2] + P[i].GravPM[2] / gd.G;
#else
	ax = P[i].GravAccel[0];
	ay = P[i].GravAccel[1];
	az = P[i].GravAccel[2];
#endif
	P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
      }


  if(gd.TypeOfOpeningCriterion == 1)
    gd.ErrTolTheta = 0;

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == gd.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] *= gd.G;


//#ifndef PERIODIC
#ifndef PMGRID
  if(gd.ComovingIntegrationOn == 0)
    {
      fac = gd.OmegaLambda * gd.Hubble * gd.Hubble;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == gd.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
//#endif

#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif

  if(ThisTask == 0)
    printf("tree is done.\n");

#else

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == gd.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#endif





  timetreelist = malloc(sizeof(double) * NTask);
  timecommlist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);
  numnodeslist = malloc(sizeof(int) * NTask);
  ewaldlist = malloc(sizeof(double) * NTask);
  nrecv = malloc(sizeof(int) * NTask);

  numnodes = Numnodestree;

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      gd.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", gd.NumCurrentTiStep, gd.Time, gd.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (gd.TotNumOfForces / 1000000000), (int) (gd.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);

      fac = NTask / ((double) gd.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	  ewaldtot += ewaldlist[i];
	}
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (gd.TreeAllocFactor * gd.MaxPart));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot, ((double) ewaldtot) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

      gd.CPU_TreeWalk += sumt / NTask;
      gd.CPU_Imbalance += sumimbalance / NTask;
      gd.CPU_CommSum += sumcomm / NTask;
    }

  free(nrecv);
  free(ewaldlist);
  free(numnodeslist);
  free(costtreelist);
  free(timecommlist);
  free(timetreelist);
}



void set_softenings(void)
{
  int i;

  if(gd.ComovingIntegrationOn)
    {
      if(gd.SofteningGas * gd.Time > gd.SofteningGasMaxPhys)
        gd.SofteningTable[0] = gd.SofteningGasMaxPhys / gd.Time;
      else
        gd.SofteningTable[0] = gd.SofteningGas;
      
      if(gd.SofteningHalo * gd.Time > gd.SofteningHaloMaxPhys)
        gd.SofteningTable[1] = gd.SofteningHaloMaxPhys / gd.Time;
      else
        gd.SofteningTable[1] = gd.SofteningHalo;
      
      if(gd.SofteningDisk * gd.Time > gd.SofteningDiskMaxPhys)
        gd.SofteningTable[2] = gd.SofteningDiskMaxPhys / gd.Time;
      else
        gd.SofteningTable[2] = gd.SofteningDisk;
      
      if(gd.SofteningBulge * gd.Time > gd.SofteningBulgeMaxPhys)
        gd.SofteningTable[3] = gd.SofteningBulgeMaxPhys / gd.Time;
      else
        gd.SofteningTable[3] = gd.SofteningBulge;
      
      if(gd.SofteningStars * gd.Time > gd.SofteningStarsMaxPhys)
        gd.SofteningTable[4] = gd.SofteningStarsMaxPhys / gd.Time;
      else
        gd.SofteningTable[4] = gd.SofteningStars;
      
      if(gd.SofteningBndry * gd.Time > gd.SofteningBndryMaxPhys)
        gd.SofteningTable[5] = gd.SofteningBndryMaxPhys / gd.Time;
      else
        gd.SofteningTable[5] = gd.SofteningBndry;
    }
  else
    {
      gd.SofteningTable[0] = gd.SofteningGas;
      gd.SofteningTable[1] = gd.SofteningHalo;
      gd.SofteningTable[2] = gd.SofteningDisk;
      gd.SofteningTable[3] = gd.SofteningBulge;
      gd.SofteningTable[4] = gd.SofteningStars;
      gd.SofteningTable[5] = gd.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    gd.ForceSoftening[i] = 2.8 * gd.SofteningTable[i];

  gd.MinGasHsml = gd.MinGasHsmlFractional * gd.ForceSoftening[0];
}


int grav_tree_compare_key(const void *a, const void *b)
{
  if(((gravdata_index *) a)->Task < (((gravdata_index *) b)->Task))
    return -1;

  if(((gravdata_index *) a)->Task > (((gravdata_index *) b)->Task))
    return +1;

  return 0;
}
