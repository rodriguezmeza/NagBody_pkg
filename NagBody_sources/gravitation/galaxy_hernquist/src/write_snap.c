//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <errno.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_math.h>

//#include "../../../General_libs/general/stdinc.h"
#include "globaldefs.h"
#include "protodefs.h"


int *Id;

int      NumPart;
  double   Time, Redshift;
  int      ntot_type[6];
  double   mtot_type[6], ScaleTable[12];
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;

FILE *fd = 0;
#define NTAB 256
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);


void MainLoop(void)
{
  char path[200], input_fname[200], basename[200];
  int  type, snapshot_number, files, n, i, k, pc;
  double MassTable[6];
  float *block;
  int   *blockid;
  int    blksize, blockmaxlen, maxidlen;
  size_t bytes;


#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}
#define BUFFER 10


  files=1;                               // number of files per snapshot
  BoxSize=0.0;
  Omega0=0.0;
  OmegaLambda=0.0;
  HubbleParam=0.0;

  cmd.rmaxh=cmd.ah*sqrt(cmd.masscut)/(1.0-sqrt(cmd.masscut));

  cmd.rmaxb=cmd.ab*sqrt(cmd.masscut)/(1.0-sqrt(cmd.masscut));
//

  NumPart=0;
  ntot_type[0]=cmd.Ngas;
  ntot_type[1]=cmd.Nhalo;
  ntot_type[2]=cmd.Ndisk;
  ntot_type[3]=cmd.Nbulge;
  ntot_type[4]=0;
  ntot_type[5]=0;

  mtot_type[0]=cmd.Mgas;
  mtot_type[1]=cmd.Mhalo;
  mtot_type[2]=cmd.Mdisk;
  mtot_type[3]=cmd.Mbulge;
  mtot_type[4]=0.0;
  mtot_type[5]=0.0;

    printf("\n");

if (cmd.Ngas!=0)
  {
   MassTable[0]=mtot_type[0]/((double) ntot_type[0]);
   ScaleTable[0]=cmd.ag;
   ScaleTable[1]=cmd.zg;
   ScaleTable[2]=cmd.rmaxg;
   printf("Gas Cutoff Radius %g\n", cmd.rmaxg);
  }
if (cmd.Nhalo!=0)
  { 
   MassTable[1]=mtot_type[1]/((double) ntot_type[1]);
   ScaleTable[3]=cmd.ah;
   ScaleTable[4]=cmd.gammah;
   ScaleTable[5]=cmd.rmaxh;
   printf("Halo Cutoff Radius %g\n", cmd.rmaxh);
  }
if (cmd.Ndisk!=0)
  {
   MassTable[2]=mtot_type[2]/((double) ntot_type[2]);
   ScaleTable[6]=cmd.ad;
   ScaleTable[7]=cmd.zd;
   ScaleTable[8]=cmd.rmaxd;
   printf("Disk Cutoff Radius %g\n", cmd.rmaxd);
  }
if (cmd.Nbulge!=0)
  {
   MassTable[3]=mtot_type[3]/((double) ntot_type[3]);
   ScaleTable[9]=cmd.ab;
   ScaleTable[10]=cmd.gammab;
   ScaleTable[11]=cmd.rmaxb;
   printf("Bulge Cutoff Radius %g\n", cmd.rmaxb);
  }header1.npart;

   MassTable[4]=0.0;
   MassTable[5]=0.0;



  /* fill file header */

  for(n = 0; n < 6; n++)
    {
      header1.npart[n] = ntot_type[n];
      header1.npartTotal[n] = (unsigned int) ntot_type[n];
      header1.mass[n] = MassTable[n];
      NumPart += ntot_type[n];
    }

  for(n = 0; n < 6; n++)
      printf("Masses: %g\n",MassTable[n]);

    printf("Proportions individual masses: mb/md = %g : mb/mh = %g : md/mh = %g \n",
           MassTable[3]/MassTable[2],
           MassTable[3]/MassTable[1],
           MassTable[2]/MassTable[1]);
    
    printf("Proportions total masses: Mb/Md = %g : Mb/Mh = %g : Md/Mh = %g \n",
           (ntot_type[3]*MassTable[3])/(ntot_type[2]*MassTable[2]),
           (ntot_type[3]*MassTable[3])/(ntot_type[1]*MassTable[1]),
           (ntot_type[2]*MassTable[2])/(ntot_type[1]*MassTable[1]));
    
  header1.time = Time;

#ifdef ComovingIntegrationOn
  header1.redshift = 1.0 / Time - 1;
#endif
  header1.redshift = 0;

  header1.flag_sfr = 0;
  header1.flag_feedback = 0;
  header1.flag_cooling = 0;

#ifdef COOLING
  header1.flag_cooling = 1;
#endif
#ifdef SFR
  header1.flag_sfr = 1;
  header1.flag_feedback = 1;
#endif

  header1.num_files = files;
  header1.BoxSize = BoxSize;
  header1.Omega0 = Omega0;
  header1.OmegaLambda = OmegaLambda;
  header1.HubbleParam = HubbleParam;


  allocate_memory();

  /* open file and write header */


        if(!(fd = fopen(cmd.snapoutfile, "w")))
            {
              printf("can't open file `%s' for writing snapshot.\n", basename);
              exit(0);
            }

          blksize = sizeof(header1);
          SKIP;
          my_fwrite(&header1, sizeof(header1), 1, fd);
          SKIP;

     mtot_type[3]+=cmd.Mdisk;
     create_hernquist(mtot_type, ntot_type, ScaleTable);

 if (cmd.Ndisk!=0)
     mtot_type[3]=cmd.Mbulge;
     create_disk(mtot_type, ntot_type, ScaleTable);



/* ------------WRITE------------- */

  if(!(block = malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double) bytes);
      exit(0);
    }

       blockmaxlen = bytes / (3 * sizeof(float));
       maxidlen = bytes / (sizeof(int));
       blockid = (int *) block;

  /* write coordinates */

          blksize = NumPart * sizeof(float) * 3;
          SKIP;

  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
        {
          block[3 * pc + k] = P[i].Pos[k];
        }

      pc++;

      if(pc == blockmaxlen)
        {
          my_fwrite(block, sizeof(float), 3 * pc, fd);
          pc = 0;
        }
    }

  if(pc > 0)
         my_fwrite(block, sizeof(float), 3 * pc, fd);
         SKIP;


  /* write velocities */

          blksize = NumPart * sizeof(float) * 3;
          SKIP;

  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
        block[3 * pc + k] = P[i].Vel[k];

      pc++;

      if(pc == blockmaxlen)
        {
          my_fwrite(block, sizeof(float), 3 * pc, fd);
          pc = 0;
        }
    }
  if(pc > 0)
         my_fwrite(block, sizeof(float), 3 * pc, fd);
         SKIP;



  /* write particle ID */

       blksize = sizeof(int) * NumPart;
       SKIP;

  for(i = 0, pc = 0; i < NumPart; i++)
    {
      blockid[pc] = i+1;
      pc++;

      if(pc == maxidlen)
        {
          my_fwrite(blockid, sizeof(int), pc, fd);
          pc = 0;
        }
    }
  if(pc > 0)
    {
          my_fwrite(blockid, sizeof(int), pc, fd);
    }
         SKIP;


  /* write zero temperatures if needed */

  if (cmd.Ngas > 0)
  {
         blksize = sizeof(float) * NumPart;
         SKIP;

  for(i = 0, pc = 0; i < NumPart; i++)
    {
      block[pc] = 0;

      pc++;

      if(pc == blockmaxlen)
        {
          my_fwrite(block, sizeof(float), pc, fd);
          pc = 0;
        }
    }
  if(pc > 0)
         my_fwrite(block, sizeof(float), pc, fd);
         SKIP;
  }

  free(block);
  fclose(fd);


}


// this routine allocates the memory for the particle data.

int allocate_memory(void)
{
  printf("allocating memory...\n");

// size of particle_data data type
  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

//  P--;   /* start with offset 1 */


  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

//  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
}


/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  *  better stop.
 *   */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured: %s\n", 0, strerror(errno));
      fflush(stdout);
      exit(0);
    }
  return nwritten;
}

void EndRun(void)
{
    char   buf[200];
    FILE *fd;

    fclose(gd.outlog);
    printf("\nFinal CPU time : %lf\n\n", cputime() - gd.cpuinit);
}


