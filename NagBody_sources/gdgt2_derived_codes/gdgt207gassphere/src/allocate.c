
// QUITAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS


#include "globaldefs.h"
#include "protodefs.h"


void allocate_commbuffers(void)
{
  size_t bytes;

  Exportflag = malloc(NTask * sizeof(char));
  DomainStartList = malloc(NTask * sizeof(int));
  DomainEndList = malloc(NTask * sizeof(int));

  TopNodes = malloc(MAXTOPNODES * sizeof(topnode_data));

  DomainWork = malloc(MAXTOPNODES * sizeof(double));
  DomainCount = malloc(MAXTOPNODES * sizeof(int));
  DomainCountSph = malloc(MAXTOPNODES * sizeof(int));
  DomainTask = malloc(MAXTOPNODES * sizeof(int));
  DomainNodeIndex = malloc(MAXTOPNODES * sizeof(int));
  DomainTreeNodeLen = malloc(MAXTOPNODES * sizeof(FLOAT));
  DomainHmax = malloc(MAXTOPNODES * sizeof(FLOAT));
  DomainMoment = malloc(MAXTOPNODES * sizeof(DomainNODE));

  if(!(CommBuffer = malloc(bytes = gd.BufferSize * 1024 * 1024)))
    {
      printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun_mpi(0,2);
    }

  gd.BunchSizeForce =
    (gd.BufferSize * 1024 * 1024) / (sizeof(gravdata_index) + 2 * sizeof(gravdata_in));

  if(gd.BunchSizeForce & 1)
    gd.BunchSizeForce -= 1;

  GravDataIndexTable = (struct gravdata_index *) CommBuffer;
  GravDataIn = (struct gravdata_in *) (GravDataIndexTable + gd.BunchSizeForce);
  GravDataGet = GravDataIn + gd.BunchSizeForce;
  GravDataOut = GravDataIn;
  GravDataResult = GravDataGet;


  gd.BunchSizeDensity =
    (gd.BufferSize * 1024 * 1024) / (2 * sizeof(densdata_in) + 2 * sizeof(densdata_out));

  DensDataIn = (densdata_in *) CommBuffer;
  DensDataGet = DensDataIn + gd.BunchSizeDensity;
  DensDataResult = (densdata_out *) (DensDataGet + gd.BunchSizeDensity);
  DensDataPartialResult = DensDataResult + gd.BunchSizeDensity;

  gd.BunchSizeHydro =
    (gd.BufferSize * 1024 * 1024) / (2 * sizeof(hydrodata_in) + 2 * sizeof(hydrodata_out));

  HydroDataIn = (hydrodata_in *) CommBuffer;
  HydroDataGet = HydroDataIn + gd.BunchSizeHydro;
  HydroDataResult = (hydrodata_out *) (HydroDataGet + gd.BunchSizeHydro);
  HydroDataPartialResult = HydroDataResult + gd.BunchSizeHydro;

  gd.BunchSizeDomain =
    (gd.BufferSize * 1024 * 1024) / (sizeof(particle_data) + sizeof(sph_particle_data) +
				      sizeof(peanokey));

  if(gd.BunchSizeDomain & 1)
    gd.BunchSizeDomain -= 1;

  DomainPartBuf = (particle_data *) CommBuffer;
  DomainSphBuf = (sph_particle_data *) (DomainPartBuf + gd.BunchSizeDomain);
  DomainKeyBuf = (peanokey *) (DomainSphBuf + gd.BunchSizeDomain);


  if(ThisTask == 0)
    {
      printf("\nAllocated %d MByte communication buffer per processor.\n\n", gd.BufferSize);
      printf("Communication buffer has room for %d particles in gravity computation\n", gd.BunchSizeForce);
      printf("Communication buffer has room for %d particles in density computation\n", gd.BunchSizeDensity);
      printf("Communication buffer has room for %d particles in hydro computation\n", gd.BunchSizeHydro);
      printf("Communication buffer has room for %d particles in domain decomposition\n", gd.BunchSizeDomain);
      printf("\n");
    }
}



void allocate_memory(void)
{
  size_t bytes;
  double bytes_tot = 0;

  if(gd.MaxPart > 0)
    {
      if(!(P = malloc(bytes = gd.MaxPart * sizeof(particle_data))))
	{
	  printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun_mpi(0,1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for particle storage. %d\n\n", bytes_tot / (1024.0 * 1024.0), sizeof(particle_data));
    }

  if(gd.MaxPartSph > 0)
    {
      bytes_tot = 0;

      if(!(SphP = malloc(bytes = gd.MaxPartSph * sizeof(sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphP' (%g MB) %d.\n", bytes / (1024.0 * 1024.0), sizeof(sph_particle_data));
	  endrun_mpi(0,1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of SPH data. %d\n\n", bytes_tot / (1024.0 * 1024.0), sizeof(sph_particle_data));
    }
}



void free_memory(void)
{
  if(gd.MaxPartSph > 0)
    free(SphP);

  if(gd.MaxPart > 0)
    free(P);
}

