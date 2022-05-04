
// DEJAR COMOVING INTEGRATION ON; DEJAR PERIODIC; DEJAR MAKEGLASS
//
// DEJAR LOS IF´S DE ComovingIntegrationOn,
//
//
// DEJAR LAS VARIABLES Y CONSTANTES ASOCIADAS A ComovingIntegrationOn AND PeriodicBoundariesOn
// BoxSize, Hubble, Omegas ...
//

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <time.h>
//#include <mpi.h>

#include "globaldefs.h"
#include "protodefs.h"
//#include "../../../General_libs/mpi/mpi_proto.h"


static int last;

#define NTAB 1000
static float tabfac, shortrange_table[NTAB], shortrange_table_potential[NTAB];

static int first_flag = 0;


//#ifdef PERIODIC
#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
#define EN  64
static FLOAT fcorrx[EN + 1][EN + 1][EN + 1];
static FLOAT fcorry[EN + 1][EN + 1][EN + 1];
static FLOAT fcorrz[EN + 1][EN + 1][EN + 1];
static FLOAT potcorr[EN + 1][EN + 1][EN + 1];
static double fac_intp;
//#endif


int force_treebuild(int npart)
{
  Numnodestree = force_treebuild_single(npart);

  force_update_pseudoparticles();

  force_flag_localnodes();

  TimeOfLastTreeConstruction = gd.Time;

  return Numnodestree;
}


int force_treebuild_single(int npart)
{
  int i, j, subnode = 0, parent, numnodes;
  int nfree, th, nn, no;
  NODE *nfreep;
  double lenhalf, epsilon;
  peanokey key;


  nfree = gd.MaxPart;
  nfreep = &Nodes[nfree];

  nfreep->len = DomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];
  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;

  numnodes = 1;
  nfreep++;
  nfree++;

  force_create_empty_nodes(gd.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);

  nfreep = &Nodes[nfree];
  parent = -1;

  for(i = 0; i < npart; i++)
    {

      epsilon = gd.ForceSoftening[P[i].Type];

      key = peano_hilbert_key((P[i].Pos[0] - DomainCorner[0]) * DomainFac,
			      (P[i].Pos[1] - DomainCorner[1]) * DomainFac,
			      (P[i].Pos[2] - DomainCorner[2]) * DomainFac, BITS_PER_DIMENSION);

      no = 0;
      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;
      th = DomainNodeIndex[no];

      while(1)
	{
	  if(th >= gd.MaxPart)
	    {
	      subnode = 0;
	      if(P[i].Pos[0] > Nodes[th].center[0])
		subnode += 1;
	      if(P[i].Pos[1] > Nodes[th].center[1])
		subnode += 2;
	      if(P[i].Pos[2] > Nodes[th].center[2])
		subnode += 4;

	      nn = Nodes[th].u.suns[subnode];

	      if(nn >= 0)
		{
		  parent = th;
		  th = nn;
		}
	      else
		{
		  Nodes[th].u.suns[subnode] = i;
		  break;
		}
	    }
	  else
	    {
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

	      if(subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if(subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if(subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;

	      subnode = 0;
	      if(P[th].Pos[0] > nfreep->center[0])
		subnode += 1;
	      if(P[th].Pos[1] > nfreep->center[1])
		subnode += 2;
	      if(P[th].Pos[2] > nfreep->center[2])
		subnode += 4;
#ifndef NOTREERND
	      if(nfreep->len < 1.0e-3 * epsilon)
		{
		  subnode = (int) (8.0 * get_random_number((0xffff & P[i].ID) + P[i].GravCost));
		  P[i].GravCost += 1;
		  if(subnode >= 8)
		    subnode = 7;
		}
#endif
	      nfreep->u.suns[subnode] = th;

	      th = nfree;

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if((numnodes) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("for particle %d\n", i);
		  dump_particles();
		  endrun_mpi(0,1);
		}
	    }
	}
    }

  force_insert_pseudo_particles();

  last = -1;

  force_update_node_recursive(gd.MaxPart, -1, -1);

  if(last >= gd.MaxPart)
    {
      if(last >= gd.MaxPart + MaxNodes)
	Nextnode[last - MaxNodes] = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  return numnodes;
}


void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
			      int *nextfree)
{
  int i, j, k, n, sub, count;

  if(TopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
	for(j = 0; j < 2; j++)
	  for(k = 0; k < 2; k++)
	    {
	      sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

	      count = i + 2 * j + 4 * k;

	      Nodes[no].u.suns[count] = *nextfree;

	      Nodes[*nextfree].len = 0.5 * Nodes[no].len;
	      Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * 0.25 * Nodes[no].len;

	      for(n = 0; n < 8; n++)
		Nodes[*nextfree].u.suns[n] = -1;

	      if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
		DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;

	      *nextfree = *nextfree + 1;
	      *nodecount = *nodecount + 1;

	      if((*nodecount) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("in create empty nodes\n");
		  dump_particles();
		  endrun_mpi(0,11);
		}

	      force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
				       bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
	    }
    }
}


void force_insert_pseudo_particles(void)
{
  int i, index, subnode, nn, th;

  for(i = 0; i < NTopleaves; i++)
    {
      index = DomainNodeIndex[i];

      DomainMoment[i].mass = 0;
      DomainMoment[i].s[0] = Nodes[index].center[0];
      DomainMoment[i].s[1] = Nodes[index].center[1];
      DomainMoment[i].s[2] = Nodes[index].center[2];
    }

  for(i = 0; i < NTopleaves; i++)
    {
      if(i < DomainMyStart || i > DomainMyLast)
	{
	  th = gd.MaxPart;

	  while(1)
	    {
	      if(th >= gd.MaxPart)
		{
		  if(th >= gd.MaxPart + MaxNodes)
		    endrun_mpi(0,888);

		  subnode = 0;
		  if(DomainMoment[i].s[0] > Nodes[th].center[0])
		    subnode += 1;
		  if(DomainMoment[i].s[1] > Nodes[th].center[1])
		    subnode += 2;
		  if(DomainMoment[i].s[2] > Nodes[th].center[2])
		    subnode += 4;

		  nn = Nodes[th].u.suns[subnode];

		  if(nn >= 0)
		    {
		      th = nn;
		    }
		  else
		    {
		      Nodes[th].u.suns[subnode] = gd.MaxPart + MaxNodes + i;

		      break;
		    }
		}
	      else
		{
		  endrun_mpi(0,889);
		}
	    }
	}
    }
}


void force_update_node_recursive(int no, int sib, int father)
{
  int j, jj, p, pp, nextsib, suns[8];
  FLOAT hmax;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  FLOAT maxsoft;
#endif
#endif
  particle_data *pa;
  double s[3], vs[3], mass;

  if(no >= gd.MaxPart && no < gd.MaxPart + MaxNodes)
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];

      if(last >= 0)
	{
	  if(last >= gd.MaxPart)
	    {
	      if(last >= gd.MaxPart + MaxNodes)
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;
      vs[0] = 0;
      vs[1] = 0;
      vs[2] = 0;
      hmax = 0;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      maxsofttype = 7;
      diffsoftflag = 0;
#else
      maxsoft = 0;
#endif
#endif

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, no);


	      if(p >= gd.MaxPart)
		{
		  if(p >= gd.MaxPart + MaxNodes)
		    {
		      /* nothing to be done here because the mass of the
		       * pseudo-particle is still zero. This will be changed
		       * later.
		       */
		    }
		  else
		    {
		      mass += Nodes[p].u.d.mass;
		      s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
		      s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
		      s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
		      vs[0] += Nodes[p].u.d.mass * Extnodes[p].vs[0];
		      vs[1] += Nodes[p].u.d.mass * Extnodes[p].vs[1];
		      vs[2] += Nodes[p].u.d.mass * Extnodes[p].vs[2];

		      if(Extnodes[p].hmax > hmax)
			hmax = Extnodes[p].hmax;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		      diffsoftflag |= (Nodes[p].u.d.bitflags >> 5) & 1;

		      if(maxsofttype == 7)
			{
			  maxsofttype = (Nodes[p].u.d.bitflags >> 2) & 7;
			}
		      else
			{
			  if(((Nodes[p].u.d.bitflags >> 2) & 7) != 7)
			    {
			      if(gd.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] >
				 gd.ForceSoftening[maxsofttype])
				{
				  maxsofttype = ((Nodes[p].u.d.bitflags >> 2) & 7);
				  diffsoftflag = 1;
				}
			      else
				{
				  if(gd.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] <
				     gd.ForceSoftening[maxsofttype])
				    diffsoftflag = 1;
				}
			    }
			}
#else
		      if(Nodes[p].maxsoft > maxsoft)
			maxsoft = Nodes[p].maxsoft;
#endif
#endif
		    }
		}
	      else
		{
		  pa = &P[p];

		  mass += pa->Mass;
		  s[0] += pa->Mass * pa->Pos[0];
		  s[1] += pa->Mass * pa->Pos[1];
		  s[2] += pa->Mass * pa->Pos[2];
		  vs[0] += pa->Mass * pa->Vel[0];
		  vs[1] += pa->Mass * pa->Vel[1];
		  vs[2] += pa->Mass * pa->Vel[2];

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		  if(maxsofttype == 7)
		    {
		      maxsofttype = pa->Type;
		    }
		  else
		    {
		      if(gd.ForceSoftening[pa->Type] > gd.ForceSoftening[maxsofttype])
			{
			  maxsofttype = pa->Type;
			  diffsoftflag = 1;
			}
		      else
			{
			  if(gd.ForceSoftening[pa->Type] < gd.ForceSoftening[maxsofttype])
			    diffsoftflag = 1;
			}
		    }
#else
		  if(pa->Type == 0)
		    {
		      if(SphP[p].Hsml > maxsoft)
			maxsoft = SphP[p].Hsml;
		    }
		  else
		    {
		      if(gd.ForceSoftening[pa->Type] > maxsoft)
			maxsoft = gd.ForceSoftening[pa->Type];
		    }
#endif
#endif
		  if(pa->Type == 0)
		    if(SphP[p].Hsml > hmax)
		      hmax = SphP[p].Hsml;
		}
	    }
	}

      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	  vs[0] /= mass;
	  vs[1] /= mass;
	  vs[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      Nodes[no].u.d.bitflags = 4 * maxsofttype + 32 * diffsoftflag;
#else
      Nodes[no].u.d.bitflags = 0;
      Nodes[no].maxsoft = maxsoft;
#endif
#else
      Nodes[no].u.d.bitflags = 0;
#endif

      Extnodes[no].vs[0] = vs[0];
      Extnodes[no].vs[1] = vs[1];
      Extnodes[no].vs[2] = vs[2];
      Extnodes[no].hmax = hmax;

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;
    }
  else
    {
      if(last >= 0)
	{
	  if(last >= gd.MaxPart)
	    {
	      if(last >= gd.MaxPart + MaxNodes)
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if(no < gd.MaxPart)
	Father[no] = father;
    }

}


void force_update_pseudoparticles(void)
{
  force_exchange_pseudodata();

  force_treeupdate_pseudos();
}


void force_exchange_pseudodata(void)
{
  int i, no;
  MPI_Status status;
  int level, sendTask, recvTask;

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      DomainMoment[i].s[0] = Nodes[no].u.d.s[0];
      DomainMoment[i].s[1] = Nodes[no].u.d.s[1];
      DomainMoment[i].s[2] = Nodes[no].u.d.s[2];
      DomainMoment[i].vs[0] = Extnodes[no].vs[0];
      DomainMoment[i].vs[1] = Extnodes[no].vs[1];
      DomainMoment[i].vs[2] = Extnodes[no].vs[2];
      DomainMoment[i].mass = Nodes[no].u.d.mass;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
#else
      DomainMoment[i].maxsoft = Nodes[no].maxsoft;
#endif
#endif
    }

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainMoment[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(DomainNODE),
		     MPI_BYTE, recvTask, TAG_DMOM,
		     &DomainMoment[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(DomainNODE),
		     MPI_BYTE, recvTask, TAG_DMOM, MPI_COMM_WORLD, &status);
    }
}


void force_treeupdate_pseudos(void)
{
  int i, k, no;
  FLOAT sold[3], vsold[3], snew[3], vsnew[3], massold, massnew, mm;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  FLOAT maxsoft;
#endif
#endif

  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	for(k = 0; k < 3; k++)
	  {
	    sold[k] = Nodes[no].u.d.s[k];
	    vsold[k] = Extnodes[no].vs[k];
	  }
	massold = Nodes[no].u.d.mass;

	for(k = 0; k < 3; k++)
	  {
	    snew[k] = DomainMoment[i].s[k];
	    vsnew[k] = DomainMoment[i].vs[k];
	  }
	massnew = DomainMoment[i].mass;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	maxsofttype = (DomainMoment[i].bitflags >> 2) & 7;
	diffsoftflag = (DomainMoment[i].bitflags >> 5) & 1;
#else
	maxsoft = DomainMoment[i].maxsoft;
#endif
#endif
	do
	  {
	    mm = Nodes[no].u.d.mass + massnew - massold;
	    for(k = 0; k < 3; k++)
	      {
		if(mm > 0)
		  {
		    Nodes[no].u.d.s[k] =
		      (Nodes[no].u.d.mass * Nodes[no].u.d.s[k] + massnew * snew[k] - massold * sold[k]) / mm;
		    Extnodes[no].vs[k] =
		      (Nodes[no].u.d.mass * Extnodes[no].vs[k] + massnew * vsnew[k] -
		       massold * vsold[k]) / mm;
		  }
	      }
	    Nodes[no].u.d.mass = mm;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	    diffsoftflag |= (Nodes[no].u.d.bitflags >> 5) & 1;

	    if(maxsofttype == 7)
	      maxsofttype = (Nodes[no].u.d.bitflags >> 2) & 7;
	    else
	      {
		if(((Nodes[no].u.d.bitflags >> 2) & 7) != 7)
		  {
		    if(gd.ForceSoftening[((Nodes[no].u.d.bitflags >> 2) & 7)] >
		       gd.ForceSoftening[maxsofttype])
		      {
			maxsofttype = ((Nodes[no].u.d.bitflags >> 2) & 7);
			diffsoftflag = 1;
		      }
		    else
		      {
			if(gd.ForceSoftening[((Nodes[no].u.d.bitflags >> 2) & 7)] <
			   gd.ForceSoftening[maxsofttype])
			  diffsoftflag = 1;
		      }
		  }
	      }

	    Nodes[no].u.d.bitflags = (Nodes[no].u.d.bitflags & 3) + 4 * maxsofttype + 32 * diffsoftflag;
#else
	    if(Nodes[no].maxsoft < maxsoft)
	      Nodes[no].maxsoft = maxsoft;
	    maxsoft = Nodes[no].maxsoft;
#endif
#endif
	    no = Nodes[no].u.d.father;

	  }
	while(no >= 0);
      }
}


void force_flag_localnodes(void)
{
  int no, i;

  for(i = 0; i < NTopleaves; i++)
    {
      no = DomainNodeIndex[i];

      while(no >= 0)
	{
	  if((Nodes[no].u.d.bitflags & 1))
	    break;

	  Nodes[no].u.d.bitflags |= 1;

	  no = Nodes[no].u.d.father;
	}
    }

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      /*
         if(DomainMoment[i].mass > 0)
       */
      {
	no = DomainNodeIndex[i];

	while(no >= 0)
	  {
	    if((Nodes[no].u.d.bitflags & 2))
	      break;

	    Nodes[no].u.d.bitflags |= 2;

	    no = Nodes[no].u.d.father;
	  }
      }
    }
}


void force_update_len(void)
{
  int i, no;
  MPI_Status status;
  int level, sendTask, recvTask;

  force_update_node_len_local();

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      DomainTreeNodeLen[i] = Nodes[no].len;
    }

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainTreeNodeLen[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_NODELEN,
		     &DomainTreeNodeLen[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_NODELEN, MPI_COMM_WORLD, &status);
    }

  force_update_node_len_toptree();
}


void force_update_node_len_local(void)
{
  int i, p, k, no;
  FLOAT dist, distmax;

  for(i = 0; i < NumPart; i++)
    {
      no = Father[i];

      for(k = 0, distmax = 0; k < 3; k++)
	{
	  dist = P[i].Pos[k] - Nodes[no].center[k];
	  if(dist < 0)
	    dist = -dist;
	  if(dist > distmax)
	    distmax = dist;
	}

      if(distmax + distmax > Nodes[no].len)
	{
	  Nodes[no].len = distmax + distmax;
	  p = Nodes[no].u.d.father;

	  while(p >= 0)
	    {
	      distmax = Nodes[p].center[0] - Nodes[no].center[0];
	      if(distmax < 0)
		distmax = -distmax;
	      distmax = distmax + distmax + Nodes[no].len;

	      if(0.999999 * distmax > Nodes[p].len)
		{
		  Nodes[p].len = distmax;
		  no = p;
		  p = Nodes[p].u.d.father;
		}
	      else
		break;
	    }
	}
    }
}


void force_update_node_len_toptree(void)
{
  int i, no, p;
  FLOAT distmax;

  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	if(Nodes[no].len < DomainTreeNodeLen[i])
	  Nodes[no].len = DomainTreeNodeLen[i];

	p = Nodes[no].u.d.father;

	while(p >= 0)
	  {
	    distmax = Nodes[p].center[0] - Nodes[no].center[0];
	    if(distmax < 0)
	      distmax = -distmax;
	    distmax = distmax + distmax + Nodes[no].len;

	    if(0.999999 * distmax > Nodes[p].len)
	      {
		Nodes[p].len = distmax;
		no = p;
		p = Nodes[p].u.d.father;
	      }
	    else
	      break;
	  }
      }
}


void force_update_hmax(void)
{
  int i, no;
  MPI_Status status;
  int level, sendTask, recvTask;

  force_update_node_hmax_local();

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      DomainHmax[i] = Extnodes[no].hmax;
    }

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainHmax[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_HMAX,
		     &DomainHmax[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_HMAX, MPI_COMM_WORLD, &status);
    }

  force_update_node_hmax_toptree();
}


void force_update_node_hmax_local(void)
{
  int i, p, no;

  for(i = 0; i < N_gas; i++)
    {

      no = Father[i];

      if(SphP[i].Hsml > Extnodes[no].hmax)
	{

	  Extnodes[no].hmax = SphP[i].Hsml;
	  p = Nodes[no].u.d.father;

	  while(p >= 0)
	    {
	      if(Extnodes[no].hmax > Extnodes[p].hmax)
		{
		  Extnodes[p].hmax = Extnodes[no].hmax;
		  no = p;
		  p = Nodes[p].u.d.father;
		}
	      else
		break;
	    }
	}
    }
}


void force_update_node_hmax_toptree(void)
{

  int i, no, p;

  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	if(Extnodes[no].hmax < DomainHmax[i])
	  Extnodes[no].hmax = DomainHmax[i];

	p = Nodes[no].u.d.father;

	while(p >= 0)
	  {
	    if(Extnodes[no].hmax > Extnodes[p].hmax)
	      {
		Extnodes[p].hmax = Extnodes[no].hmax;
		no = p;
		p = Nodes[p].u.d.father;
	      }
	    else
	      break;
	  }
      }
}


int force_treeevaluate(int target, int mode, double *ewaldcountsum)
{
  NODE *nop = 0;
  int no, ninteractions, ptype;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z, aold;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
//#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = gd.BoxSize;
  boxhalf = 0.5 * gd.BoxSize;
//#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = gd.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
      aold = gd.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }

#ifndef UNEQUALSOFTENINGS
  h = gd.ForceSoftening[ptype];
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;
#endif
  no = gd.MaxPart;

  while(no >= 0)
    {
      if(no < gd.MaxPart)
	{

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;

	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= gd.MaxPart + MaxNodes)
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (gd.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }
	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;

	  mass = nop->u.d.mass;
	}
//#ifdef PERIODIC
      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);
//#endif
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < gd.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < gd.ForceSoftening[P[no].Type])
		h = gd.ForceSoftening[P[no].Type];
	    }
#else
	  h = gd.ForceSoftening[ptype];
	  if(h < gd.ForceSoftening[P[no].Type])
	    h = gd.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }


	  if(gd.ErrTolTheta)
	    {
	      if(nop->len * nop->len > r2 * gd.ErrTolTheta * gd.ErrTolTheta)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = gd.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7)
            {
              if(mass > 0)
                endrun_mpi(0,986);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < gd.ForceSoftening[maxsofttype])
                {
                  h = gd.ForceSoftening[maxsofttype];
                  if(r2 < h * h)
                    {
                      if(((nop->u.d.bitflags >> 5) & 1))
                        {
                          no = nop->u.d.nextnode;
                          continue;
                        }
                    }
                }
            }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	fac = mass / (r2 * r);
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
	  h3_inv = h_inv * h_inv * h_inv;
#endif
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	}

      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;

      ninteractions++;
    }

  if(mode == 0)
    {
      P[target].GravAccel[0] = acc_x;
      P[target].GravAccel[1] = acc_y;
      P[target].GravAccel[2] = acc_z;
      P[target].GravCost = ninteractions;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
      GravDataResult[target].w.Ninteractions = ninteractions;
    }

//#ifdef PERIODIC
  *ewaldcountsum += force_treeevaluate_ewald_correction(target, mode, pos_x, pos_y, pos_z, aold);
//#endif

  return ninteractions;
}


#ifdef PMGRID
int force_treeevaluate_shortrange(int target, int mode)
{
  NODE *nop = 0;
  int no, ptype, ninteractions, tabindex;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z, aold;
  double eff_dist;
  double rcut, asmth, asmthfac, rcut2, dist;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
//#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = gd.BoxSize;
  boxhalf = 0.5 * gd.BoxSize;
//#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = gd.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
      aold = gd.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }

  rcut = gd.Rcut[0];
  asmth = gd.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(((1 << ptype) & (PLACEHIGHRESREGION)))
    {
      rcut = gd.Rcut[1];
      asmth = gd.Asmth[1];
    }
#endif
  rcut2 = rcut * rcut;

  asmthfac = 0.5 / asmth * (NTAB / 3.0);

#ifndef UNEQUALSOFTENINGS
  h = gd.ForceSoftening[ptype];
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;
#endif
  no = gd.MaxPart;

  while(no >= 0)
    {
      if(no < gd.MaxPart)
	{
	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
//#ifdef PERIODIC
	  dx = NEAREST(dx);
	  dy = NEAREST(dy);
	  dz = NEAREST(dz);
//#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  mass = P[no].Mass;
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < gd.ForceSoftening[P[no].Type])
		h = gd.ForceSoftening[P[no].Type];
	    }
#else
	  h = gd.ForceSoftening[ptype];
	  if(h < gd.ForceSoftening[P[no].Type])
	    h = gd.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else
	{
	  if(no >= gd.MaxPart + MaxNodes)
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (gd.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  mass = nop->u.d.mass;

	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
//#ifdef PERIODIC
	  dx = NEAREST(dx);
	  dy = NEAREST(dy);
	  dz = NEAREST(dz);
//#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 > rcut2)
	    {
	      eff_dist = rcut + 0.5 * nop->len;
//#ifdef PERIODIC
	      dist = NEAREST(nop->center[0] - pos_x);
/*#else
	      dist = nop->center[0] - pos_x;
#endif */
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
//#ifdef PERIODIC
	      dist = NEAREST(nop->center[1] - pos_y);
/*#else
	      dist = nop->center[1] - pos_y;
#endif */
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
//#ifdef PERIODIC
	      dist = NEAREST(nop->center[2] - pos_z);
/* #else
	      dist = nop->center[2] - pos_z;
#endif */
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }


	  if(gd.ErrTolTheta)
	    {
	      if(nop->len * nop->len > r2 * gd.ErrTolTheta * gd.ErrTolTheta)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = gd.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7)
            {
              if(mass > 0)
                endrun_mpi(0,987);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < gd.ForceSoftening[maxsofttype])
                {
                  h = gd.ForceSoftening[maxsofttype];
                  if(r2 < h * h)
                    {
                      if(((nop->u.d.bitflags >> 5) & 1))
                        {
                          no = nop->u.d.nextnode;
                          
                          continue;
                        }
                    }
                }
            }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif
	  no = nop->u.d.sibling;

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	fac = mass / (r2 * r);
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
	  h3_inv = h_inv * h_inv * h_inv;
#endif
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	}

      tabindex = (int) (asmthfac * r);

      if(tabindex < NTAB)
	{
	  fac *= shortrange_table[tabindex];

	  acc_x += dx * fac;
	  acc_y += dy * fac;
	  acc_z += dz * fac;

	  ninteractions++;
	}
    }


  if(mode == 0)
    {
      P[target].GravAccel[0] = acc_x;
      P[target].GravAccel[1] = acc_y;
      P[target].GravAccel[2] = acc_z;
      P[target].GravCost = ninteractions;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
      GravDataResult[target].w.Ninteractions = ninteractions;
    }

  return ninteractions;
}

#endif


//#ifdef PERIODIC
int force_treeevaluate_ewald_correction(int target, int mode, double pos_x, double pos_y, double pos_z,
					double aold)
{
  NODE *nop = 0;
  int no, cost;
  double dx, dy, dz, mass, r2;
  int signx, signy, signz;
  int i, j, k, openflag;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double acc_x, acc_y, acc_z;
  double boxsize, boxhalf;

  boxsize = gd.BoxSize;
  boxhalf = 0.5 * gd.BoxSize;

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  cost = 0;

  no = gd.MaxPart;

  while(no >= 0)
    {
      if(no < gd.MaxPart)
	{

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= gd.MaxPart + MaxNodes)
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (gd.MaxPart + MaxNodes)]] = 1;
		}

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
	  mass = nop->u.d.mass;
	}

      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);

      if(no < gd.MaxPart)
	no = Nextnode[no];
      else
	{
	  openflag = 0;

	  r2 = dx * dx + dy * dy + dz * dz;

	  if(gd.ErrTolTheta)
	    {
	      if(nop->len * nop->len > r2 * gd.ErrTolTheta * gd.ErrTolTheta)
		{
		  openflag = 1;
		}
	    }
	  else
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  openflag = 1;
		}
	      else
		{
		  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
			{
			  if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			    {
			      openflag = 1;
			    }
			}
		    }
		}
	    }

	  if(openflag)
	    {

	      u = nop->center[0] - pos_x;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      u = nop->center[1] - pos_y;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      u = nop->center[2] - pos_z;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(nop->len > 0.20 * boxsize)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }

	  no = nop->u.d.sibling;

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 1))
		continue;
	    }
	}

      if(dx < 0)
	{
	  dx = -dx;
	  signx = +1;
	}
      else
	signx = -1;

      if(dy < 0)
	{
	  dy = -dy;
	  signy = +1;
	}
      else
	signy = -1;

      if(dz < 0)
	{
	  dz = -dz;
	  signz = +1;
	}
      else
	signz = -1;

      u = dx * fac_intp;
      i = (int) u;
      if(i >= EN)
	i = EN - 1;
      u -= i;
      v = dy * fac_intp;
      j = (int) v;
      if(j >= EN)
	j = EN - 1;
      v -= j;
      w = dz * fac_intp;
      k = (int) w;
      if(k >= EN)
	k = EN - 1;
      w -= k;

      f1 = (1 - u) * (1 - v) * (1 - w);
      f2 = (1 - u) * (1 - v) * (w);
      f3 = (1 - u) * (v) * (1 - w);
      f4 = (1 - u) * (v) * (w);
      f5 = (u) * (1 - v) * (1 - w);
      f6 = (u) * (1 - v) * (w);
      f7 = (u) * (v) * (1 - w);
      f8 = (u) * (v) * (w);

      acc_x += mass * signx * (fcorrx[i][j][k] * f1 +
			       fcorrx[i][j][k + 1] * f2 +
			       fcorrx[i][j + 1][k] * f3 +
			       fcorrx[i][j + 1][k + 1] * f4 +
			       fcorrx[i + 1][j][k] * f5 +
			       fcorrx[i + 1][j][k + 1] * f6 +
			       fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8);

      acc_y += mass * signy * (fcorry[i][j][k] * f1 +
			       fcorry[i][j][k + 1] * f2 +
			       fcorry[i][j + 1][k] * f3 +
			       fcorry[i][j + 1][k + 1] * f4 +
			       fcorry[i + 1][j][k] * f5 +
			       fcorry[i + 1][j][k + 1] * f6 +
			       fcorry[i + 1][j + 1][k] * f7 + fcorry[i + 1][j + 1][k + 1] * f8);

      acc_z += mass * signz * (fcorrz[i][j][k] * f1 +
			       fcorrz[i][j][k + 1] * f2 +
			       fcorrz[i][j + 1][k] * f3 +
			       fcorrz[i][j + 1][k + 1] * f4 +
			       fcorrz[i + 1][j][k] * f5 +
			       fcorrz[i + 1][j][k + 1] * f6 +
			       fcorrz[i + 1][j + 1][k] * f7 + fcorrz[i + 1][j + 1][k + 1] * f8);
      cost++;
    }

  if(mode == 0)
    {
      P[target].GravAccel[0] += acc_x;
      P[target].GravAccel[1] += acc_y;
      P[target].GravAccel[2] += acc_z;
      P[target].GravCost += cost;
    }
  else
    {
      GravDataResult[target].u.Acc[0] += acc_x;
      GravDataResult[target].u.Acc[1] += acc_y;
      GravDataResult[target].u.Acc[2] += acc_z;
      GravDataResult[target].w.Ninteractions += cost;
    }

  return cost;
}

//#endif


void force_treeevaluate_potential(int target, int mode)
{
  NODE *nop = 0;
  int no, ptype;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pot, pos_x, pos_y, pos_z, aold;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
//#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = gd.BoxSize;
  boxhalf = 0.5 * gd.BoxSize;
//#endif

  pot = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = gd.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
      aold = gd.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }

#ifndef UNEQUALSOFTENINGS
  h = gd.ForceSoftening[ptype];
  h_inv = 1.0 / h;
#endif
  no = gd.MaxPart;

  while(no >= 0)
    {
      if(no < gd.MaxPart)
	{

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= gd.MaxPart + MaxNodes)
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (gd.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
	  mass = nop->u.d.mass;
	}

//#ifdef PERIODIC
      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);
//#endif
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < gd.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < gd.ForceSoftening[P[no].Type])
		h = gd.ForceSoftening[P[no].Type];
	    }
#else
	  h = gd.ForceSoftening[ptype];
	  if(h < gd.ForceSoftening[P[no].Type])
	    h = gd.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  if(gd.ErrTolTheta)	
	    {
	      if(nop->len * nop->len > r2 * gd.ErrTolTheta * gd.ErrTolTheta)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = gd.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
            {
              if(mass > 0)
                endrun_mpi(0,988);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < gd.ForceSoftening[maxsofttype])
                {
                  h = gd.ForceSoftening[maxsofttype];
                  if(r2 < h * h)
                    {
                      if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
                        {
                          no = nop->u.d.nextnode;
                          continue;
                        }
                    }
                }
            }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;	/* node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	pot -= mass / r;
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
#endif
	  u = r * h_inv;

	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	  pot += mass * h_inv * wp;
	}
//#ifdef PERIODIC
      pot += mass * ewald_pot_corr(dx, dy, dz);
//#endif
    }

  if(mode == 0)
    P[target].Potential = pot;
  else
    GravDataResult[target].u.Potential = pot;
}


#ifdef PMGRID
void force_treeevaluate_potential_shortrange(int target, int mode)
{
  NODE *nop = 0;
  int no, ptype, tabindex;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pot, pos_x, pos_y, pos_z, aold;
  double eff_dist, fac, rcut, asmth, asmthfac;
  double dxx, dyy, dzz;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif

//#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = gd.BoxSize;
  boxhalf = 0.5 * gd.BoxSize;
//#endif

  pot = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = gd.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
      aold = gd.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }

  rcut = gd.Rcut[0];
  asmth = gd.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(((1 << ptype) & (PLACEHIGHRESREGION)))
    {
      rcut = gd.Rcut[1];
      asmth = gd.Asmth[1];
    }
#endif
  asmthfac = 0.5 / asmth * (NTAB / 3.0);

#ifndef UNEQUALSOFTENINGS
  h = gd.ForceSoftening[ptype];
  h_inv = 1.0 / h;
#endif

  no = gd.MaxPart;

  while(no >= 0)
    {
      if(no < gd.MaxPart)	/* single particle */
	{

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= gd.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (gd.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
	  mass = nop->u.d.mass;
	}

//#ifdef PERIODIC
      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);
//#endif
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < gd.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < gd.ForceSoftening[P[no].Type])
		h = gd.ForceSoftening[P[no].Type];
	    }
#else
	  h = gd.ForceSoftening[ptype];
	  if(h < gd.ForceSoftening[P[no].Type])
	    h = gd.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  if(no >= gd.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (gd.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node which does not contain local particles */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  eff_dist = rcut + 0.5 * nop->len;

	  dxx = nop->center[0] - pos_x;	/* observe the sign ! */
	  dyy = nop->center[1] - pos_y;	/* this vector is -y in my thesis notation */
	  dzz = nop->center[2] - pos_z;
//#ifdef PERIODIC
	  dxx = NEAREST(dxx);
	  dyy = NEAREST(dyy);
	  dzz = NEAREST(dzz);
//#endif
	  if(dxx < -eff_dist || dxx > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(dyy < -eff_dist || dyy > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(dzz < -eff_dist || dzz > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(gd.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * gd.ErrTolTheta * gd.ErrTolTheta)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = gd.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
            {
              if(mass > 0)
                endrun_mpi(0,989);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < gd.ForceSoftening[maxsofttype])
                {
                  h = gd.ForceSoftening[maxsofttype];
                  if(r2 < h * h)
                    {
                      /* bit-5 signals that there are particles of
                       * different softening in the node
                       */
                      if(((nop->u.d.bitflags >> 5) & 1))
                        {
                          no = nop->u.d.nextnode;
                          continue;
                        }
                    }
                }
	    }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = gd.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif
	  no = nop->u.d.sibling;	/* node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      tabindex = (int) (r * asmthfac);

      if(tabindex < NTAB)
	{
	  fac = shortrange_table_potential[tabindex];

	  if(r >= h)
	    pot -= fac * mass / r;
	  else
	    {
#ifdef UNEQUALSOFTENINGS
	      h_inv = 1.0 / h;
#endif
	      u = r * h_inv;

	      if(u < 0.5)
		wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	      else
		wp =
		  -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						       u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
	      pot += fac * mass * h_inv * wp;
	    }
	}
    }

  if(mode == 0)
    P[target].Potential = pot;
  else
    GravDataResult[target].u.Potential = pot;
}

#endif


void force_treeallocate(int maxnodes, int maxpart)
{
  int i;
  size_t bytes;
  double allbytes = 0;
  double u;

  MaxNodes = maxnodes;

  if(!(Nodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
      endrun_mpi(0,3);
    }
  allbytes += bytes;

  if(!(Extnodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(extNODE))))
    {
      printf("failed to allocate memory for %d tree-extnodes (%g MB).\n", MaxNodes,
	     bytes / (1024.0 * 1024.0));
      endrun_mpi(0,3);
    }
  allbytes += bytes;

  Nodes = Nodes_base - gd.MaxPart;
  Extnodes = Extnodes_base - gd.MaxPart;

  if(!(Nextnode = malloc(bytes = (maxpart + MAXTOPNODES) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n", maxpart + MAXTOPNODES,
	     bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(!(Father = malloc(bytes = (maxpart) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(first_flag == 0)
    {
      first_flag = 1;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for BH-tree. %d\n\n", allbytes / (1024.0 * 1024.0),
	       sizeof(NODE) + sizeof(extNODE));

      tabfac = NTAB / 3.0;

      for(i = 0; i < NTAB; i++)
	{
	  u = 3.0 / NTAB * (i + 0.5);
	  shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
	  shortrange_table_potential[i] = erfc(u);
	}
    }
}


void force_treefree(void)
{
  free(Father);
  free(Nextnode);
  free(Extnodes_base);
  free(Nodes_base);
}


#ifdef FORCETEST
int force_treeevaluate_direct(int target, int mode)
{
  double epsilon;
  double h, h_inv, dx, dy, dz, r, r2, u, r_inv, fac;
  int i, ptype;
  double pos_x, pos_y, pos_z;
  double acc_x, acc_y, acc_z;

//#ifdef PERIODIC
  double fcorr[3];
//#endif
//#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = gd.BoxSize;
  boxhalf = 0.5 * gd.BoxSize;
//#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
    }

  for(i = 0; i < NumPart; i++)
    {
      epsilon = dmax(gd.ForceSoftening[P[i].Type], gd.ForceSoftening[ptype]);

      h = epsilon;
      h_inv = 1 / h;

      dx = P[i].Pos[0] - pos_x;
      dy = P[i].Pos[1] - pos_y;
      dz = P[i].Pos[2] - pos_z;

//#ifdef PERIODIC
      while(dx > boxhalf)
	dx -= boxsize;
      while(dy > boxhalf)
	dy -= boxsize;
      while(dz > boxhalf)
	dz -= boxsize;
      while(dx < -boxhalf)
	dx += boxsize;
      while(dy < -boxhalf)
	dy += boxsize;
      while(dz < -boxhalf)
	dz += boxsize;
//#endif
      r2 = dx * dx + dy * dy + dz * dz;

      r = sqrt(r2);

      u = r * h_inv;

      if(u >= 1)
	{
	  r_inv = 1 / r;

	  fac = P[i].Mass * r_inv * r_inv * r_inv;
	}
      else
	{
	  if(u < 0.5)
	    fac = P[i].Mass * h_inv * h_inv * h_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      P[i].Mass * h_inv * h_inv * h_inv * (21.333333333333 -
						   48.0 * u + 38.4 * u * u -
						   10.666666666667 * u * u *
						   u - 0.066666666667 / (u * u * u));
	}

      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;

//#ifdef PERIODIC
      if(u > 1.0e-5)
	{
	  ewald_corr(dx, dy, dz, fcorr);

	  acc_x += P[i].Mass * fcorr[0];
	  acc_y += P[i].Mass * fcorr[1];
	  acc_z += P[i].Mass * fcorr[2];
	}
//#endif
    }

  if(mode == 0)
    {
      P[target].GravAccelDirect[0] = acc_x;
      P[target].GravAccelDirect[1] = acc_y;
      P[target].GravAccelDirect[2] = acc_z;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
    }

  return NumPart;
}
#endif


void dump_particles(void)
{
  FILE *fd;
  char buffer[200];
  int i;

  sprintf(buffer, "particles%d.dat", ThisTask);
  fd = fopen(buffer, "w");
  my_fwrite(&NumPart, 1, sizeof(int), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(FLOAT), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(FLOAT), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);

  fclose(fd);
}


//#ifdef PERIODIC

void ewald_init(void)
{
  int i, j, k, beg, len, size, n, task, count;
  double x[3], force[3];
  char buf[200];
  FILE *fd;

  if(ThisTask == 0)
    {
      printf("initialize Ewald correction...\n");
      fflush(stdout);
    }

#ifdef DOUBLEPRECISION
  sprintf(buf, "ewald_spc_table_%d_dbl.dat", EN);
#else
  sprintf(buf, "ewald_spc_table_%d.dat", EN);
#endif

  if((fd = fopen(buf, "r")))
    {
      if(ThisTask == 0)
	{
	  printf("\nreading Ewald tables from file `%s'\n", buf);
	  fflush(stdout);
	}

      my_fread(&fcorrx[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      my_fread(&fcorry[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      my_fread(&fcorrz[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      my_fread(&potcorr[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      fclose(fd);
    }
  else
    {
      if(ThisTask == 0)
	{
	  printf("\nNo Ewald tables in file `%s' found.\nRecomputing them...\n", buf);
	  fflush(stdout);
	}

      size = (EN + 1) * (EN + 1) * (EN + 1) / NTask;

      beg = ThisTask * size;
      len = size;
      if(ThisTask == (NTask - 1))
	len = (EN + 1) * (EN + 1) * (EN + 1) - beg;

      for(i = 0, count = 0; i <= EN; i++)
	for(j = 0; j <= EN; j++)
	  for(k = 0; k <= EN; k++)
	    {
	      n = (i * (EN + 1) + j) * (EN + 1) + k;
	      if(n >= beg && n < (beg + len))
		{
		  if(ThisTask == 0)
		    {
		      if((count % (len / 20)) == 0)
			{
			  printf("%4.1f percent done\n", count / (len / 100.0));
			  fflush(stdout);
			}
		    }

		  x[0] = 0.5 * ((double) i) / EN;
		  x[1] = 0.5 * ((double) j) / EN;
		  x[2] = 0.5 * ((double) k) / EN;

		  ewald_force(i, j, k, x, force);

		  fcorrx[i][j][k] = force[0];
		  fcorry[i][j][k] = force[1];
		  fcorrz[i][j][k] = force[2];

		  if(i + j + k == 0)
		    potcorr[i][j][k] = 2.8372975;
		  else
		    potcorr[i][j][k] = ewald_psi(x);

		  count++;
		}
	    }

      for(task = 0; task < NTask; task++)
	{
	  beg = task * size;
	  len = size;
	  if(task == (NTask - 1))
	    len = (EN + 1) * (EN + 1) * (EN + 1) - beg;

#ifdef DOUBLEPRECISION
	  MPI_Bcast(&fcorrx[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorry[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorrz[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	  MPI_Bcast(&potcorr[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
#else
	  MPI_Bcast(&fcorrx[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorry[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorrz[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&potcorr[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
#endif
	}

      if(ThisTask == 0)
	{
	  printf("\nwriting Ewald tables to file `%s'\n", buf);
	  fflush(stdout);

	  if((fd = fopen(buf, "w")))
	    {
	      my_fwrite(&fcorrx[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      my_fwrite(&fcorry[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      my_fwrite(&fcorrz[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      my_fwrite(&potcorr[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      fclose(fd);
	    }
	}
    }

  fac_intp = 2 * EN / gd.BoxSize;

  for(i = 0; i <= EN; i++)
    for(j = 0; j <= EN; j++)
      for(k = 0; k <= EN; k++)
	{
	  potcorr[i][j][k] /= gd.BoxSize;
	  fcorrx[i][j][k] /= gd.BoxSize * gd.BoxSize;
	  fcorry[i][j][k] /= gd.BoxSize * gd.BoxSize;
	  fcorrz[i][j][k] /= gd.BoxSize * gd.BoxSize;
	}

  if(ThisTask == 0)
    {
      printf("initialization of periodic boundaries finished.\n");
      fflush(stdout);
    }
}


#ifdef FORCETEST
void ewald_corr(double dx, double dy, double dz, double *fper)
{
  int signx, signy, signz;
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    {
      dx = -dx;
      signx = +1;
    }
  else
    signx = -1;

  if(dy < 0)
    {
      dy = -dy;
      signy = +1;
    }
  else
    signy = -1;

  if(dz < 0)
    {
      dz = -dz;
      signz = +1;
    }
  else
    signz = -1;

  u = dx * fac_intp;
  i = (int) u;
  if(i >= EN)
    i = EN - 1;
  u -= i;
  v = dy * fac_intp;
  j = (int) v;
  if(j >= EN)
    j = EN - 1;
  v -= j;
  w = dz * fac_intp;
  k = (int) w;
  if(k >= EN)
    k = EN - 1;
  w -= k;

  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);

  fper[0] = signx * (fcorrx[i][j][k] * f1 +
		     fcorrx[i][j][k + 1] * f2 +
		     fcorrx[i][j + 1][k] * f3 +
		     fcorrx[i][j + 1][k + 1] * f4 +
		     fcorrx[i + 1][j][k] * f5 +
		     fcorrx[i + 1][j][k + 1] * f6 +
		     fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8);

  fper[1] = signy * (fcorry[i][j][k] * f1 +
		     fcorry[i][j][k + 1] * f2 +
		     fcorry[i][j + 1][k] * f3 +
		     fcorry[i][j + 1][k + 1] * f4 +
		     fcorry[i + 1][j][k] * f5 +
		     fcorry[i + 1][j][k + 1] * f6 +
		     fcorry[i + 1][j + 1][k] * f7 + fcorry[i + 1][j + 1][k + 1] * f8);

  fper[2] = signz * (fcorrz[i][j][k] * f1 +
		     fcorrz[i][j][k + 1] * f2 +
		     fcorrz[i][j + 1][k] * f3 +
		     fcorrz[i][j + 1][k + 1] * f4 +
		     fcorrz[i + 1][j][k] * f5 +
		     fcorrz[i + 1][j][k + 1] * f6 +
		     fcorrz[i + 1][j + 1][k] * f7 + fcorrz[i + 1][j + 1][k + 1] * f8);
}
#endif


double ewald_pot_corr(double dx, double dy, double dz)
{
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    dx = -dx;

  if(dy < 0)
    dy = -dy;

  if(dz < 0)
    dz = -dz;

  u = dx * fac_intp;
  i = (int) u;
  if(i >= EN)
    i = EN - 1;
  u -= i;
  v = dy * fac_intp;
  j = (int) v;
  if(j >= EN)
    j = EN - 1;
  v -= j;
  w = dz * fac_intp;
  k = (int) w;
  if(k >= EN)
    k = EN - 1;
  w -= k;

  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);

  return potcorr[i][j][k] * f1 +
    potcorr[i][j][k + 1] * f2 +
    potcorr[i][j + 1][k] * f3 +
    potcorr[i][j + 1][k + 1] * f4 +
    potcorr[i + 1][j][k] * f5 +
    potcorr[i + 1][j][k + 1] * f6 + potcorr[i + 1][j + 1][k] * f7 + potcorr[i + 1][j + 1][k + 1] * f8;
}


double ewald_psi(double x[3])
{
  double alpha, psi;
  double r, sum1, sum2, hdotx;
  double dx[3];
  int i, n[3], h[3], h2;

  alpha = 2.0;

  for(n[0] = -4, sum1 = 0; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
	  sum1 += erfc(alpha * r) / r;
	}

  for(h[0] = -4, sum2 = 0; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
	  if(h2 > 0)
	    sum2 += 1 / (M_PI * h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * cos(2 * M_PI * hdotx);
	}

  r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

  psi = M_PI / (alpha * alpha) - sum1 - sum2 + 1 / r;

  return psi;
}


void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
  double alpha, r2;
  double r, val, hdotx, dx[3];
  int i, h[3], n[3], h2;

  alpha = 2.0;

  for(i = 0; i < 3; i++)
    force[i] = 0;

  if(iii == 0 && jjj == 0 && kkk == 0)
    return;

  r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

  for(i = 0; i < 3; i++)
    force[i] += x[i] / (r2 * sqrt(r2));

  for(n[0] = -4; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

	  val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);

	  for(i = 0; i < 3; i++)
	    force[i] -= dx[i] / (r * r * r) * val;
	}

  for(h[0] = -4; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];

	  if(h2 > 0)
	    {
	      val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);

	      for(i = 0; i < 3; i++)
		force[i] -= h[i] * val;
	    }
	}
}

//#endif
