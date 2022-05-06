#include "../general/stdinc.h"

void endrun_mpi(int, int);
void *allocate_mpi(int, long int, string);
void fprintf_mpi(FILE *, int, string, ...);
void fprintf_mpi_flush(FILE *, int, string, ...);
FILE *openfile_mpi(int, char *, char mode[2]);
