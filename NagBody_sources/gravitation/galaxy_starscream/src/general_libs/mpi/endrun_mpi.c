

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include <stdarg.h>

#include "mpi_proto.h"
#include "../general/switchs.h"

void endrun_mpi(int ThisTask, int ierr)
{
	if(ierr) {
		fprintf(stdout,
			"task %d: endrun called with an error level of %d\n\n\n", 
			ThisTask, ierr);
#ifdef DEBUG
        terminate_processes();
        raise(SIGABRT);
        sleep(60);
#else
		MPI_Abort(MPI_COMM_WORLD, ierr);
#endif
		exit(0);
	}

	MPI_Finalize();
	exit(0);
}


void *allocate_mpi(int ThisTask, long int nb, string fmt)
{
    void *mem;

    mem = calloc(nb, 1); 
    if (mem == NULL) {
		printf("%s (%ld bytes).\n",fmt, nb);
		endrun_mpi(ThisTask, 1);
	}
    return (mem);
}


void fprintf_mpi(FILE *fd, int ThisTask, string fmt, ...)
{
    va_list ap;

	if (ThisTask==0) {
		va_start(ap, fmt);
		vfprintf(fd, fmt, ap);	
		va_end(ap);
	}
}

void fprintf_mpi_flush(FILE *fd, int ThisTask, string fmt, ...)
{
    va_list ap;

	if (ThisTask==0) {
		va_start(ap, fmt);
		vfprintf(fd, fmt, ap);	
		va_end(ap);
		fflush(fd);
	}
}

FILE *openfile_mpi(int ThisTask, char *buf, char mode[2])
{
	FILE *fd;

	if(!(fd = fopen(buf, mode))) {
		fprintf(stdout, "error in opening file '%s'\n", buf);
		endrun_mpi(ThisTask, 1);
	}

	return (fd);
}

