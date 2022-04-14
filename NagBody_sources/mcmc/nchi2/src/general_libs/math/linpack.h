
#ifndef _linpack_h
#define _linpack_h

#include <stdio.h>

//#include "../general/stdinc.h"

void lubksb(double  **a, int n, int *indx, double  b[]);
void ludcmp(double  **a, int n, int *indx, double  *d);

void choldc(double  **a, int n, double  p[]);
void cholsl(double  **a, int n, double  p[], double  b[], double  x[]);


#endif // ! _linpack_h


