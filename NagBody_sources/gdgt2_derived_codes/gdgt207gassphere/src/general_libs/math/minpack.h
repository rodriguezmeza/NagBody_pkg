
#ifndef _minpack_h
#define _minpack_h

#include <math.h>

// LevMarq method
void mrqmin(double x[], double y[], double sig[], int ndata, double a[],
            int ia[], int ma, double **covar, double **alpha, double *chisq,
            void (*funcs)(double, double [], double *, double [], int), double *alamda);
void covsrt(double **covar, int ma, int ia[], int mfit);
void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
            int ia[], int ma, double **alpha, double beta[], double *chisq,
            void (*funcs)(double, double [], double *, double [], int));

// BEGIN :: AMOEBA METHOD
void amoeba(double **p, double y[], int ndim, double ftol,
            double (*funk)(double []), int *iter);
double amotry(double **p, double y[], double psum[], int ndim,
              double (*funk)(double []), int ihi, double fac);
// END :: AMOEBA METHOD


// BEGIN :: POWELL METHOD
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
            double (*func)(double []));
void linmin(double p[], double xi[], int n, double *fret,
            double (*func)(double []));
double brent(double ax, double bx, double cx,
             double (*f)(double), double tol, double *xmin);
double f1dim(double x);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
            double *fc, double (*func)(double));
// END :: POWELL METHOD


// BEGIN :: FRPRMN METHOD
void frprmn(double p[], int n, double ftol, int *iter, double *fret,
            double (*func)(double []), void (*dfunc)(double [], double []));
// END :: FRPRMN METHOD



// BEGIN :: DFPMIN METHOD
void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
            double (*func)(double []), void (*dfunc)(double [], double []));
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
            double *f, double stpmax, int *check, double (*func)(double []));
// END :: DFPMIN METHOD


// BEGIN :: AMEBSA METHOD
void amebsa(double **p, double y[], int ndim, double pb[],    double *yb,
            double ftol, double (*funk)(double []), int *iter, double temptr);
double amotsa(double **p, double y[], double psum[], int ndim, double pb[],
              double *yb, double (*funk)(double []), int ihi, double *yhi, double fac);
// END :: AMEBSA METHOD

#endif // ! _minpack_h


