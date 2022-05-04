
// DEJAR COMOVING INTEGRATION ON; QUITAR PERIODIC; QUITAR MAKEGLASS


//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/general/getparam.h"

#define global

#include "globaldefs.h"
#include "protodefs.h"
#include "cmdline_defs.h"

int main(int argc, char **argv)
{
  double t0, t1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

	if(ThisTask==0) {
		stdoutstrm = openfile_mpi(ThisTask, STDOUTFILE, "w");
	}

  if(NTask <= 1)
    {
        fprintf_mpi_flush(STDOUT, ThisTask,
            "This is a massively parallel code.\n%s\n%s\n%s",
            "but you are running with 1 processor only.",
            "Compared to an equivalent serial code",
            "there is some unnecessary overhead...");
    }

  for(PTask = 0; NTask > (1 << PTask); PTask++);

  if (argc >1)
    strcpy(ParameterFile, argv[1]);

  if(argc >= 3)
      cmd.RestartFlag = atoi(argv[2]);
  else
      cmd.RestartFlag = 0;

  gd.CPU_TreeConstruction  =
  gd.CPU_TreeWalk          =
  gd.CPU_Gravity           =
  gd.CPU_Potential         =
  gd.CPU_Domain            =
  gd.CPU_Snapshot          =
  gd.CPU_Total             =
  gd.CPU_CommSum           =
  gd.CPU_Imbalance         =
  gd.CPU_Hydro             =
  gd.CPU_HydCompWalk       =
  gd.CPU_HydCommSumm       =
  gd.CPU_HydImbalance      =
  gd.CPU_EnsureNgb         =
  gd.CPU_Predict           =
  gd.CPU_TimeLine          =
  gd.CPU_PM                =
  gd.CPU_Peano             = 0;

  CPUThisRun = 0;

  t0 = second();

  if (ThisTask==0)
    InitParam(argv, defv);

  startrun(argv[0], HEAD1, HEAD2, HEAD3);

  t1 = second();
  CPUThisRun += timediff(t0, t1);
  gd.CPU_Total += timediff(t0, t1);

  run();

	if(ThisTask==0) {
		fclose(stdoutstrm);
	}

  MPI_Finalize();

  return 0;
}


