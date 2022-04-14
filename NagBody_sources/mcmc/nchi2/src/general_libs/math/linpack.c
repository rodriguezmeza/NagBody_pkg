

//#include "../general/stdinc.h"
//#include "mathfns.h"
//#include "../general/constant.h"
//#include "quads.h"
//#include "numrec.h"

//#include "nrutil.h"
#include "numrec.h"
#include "linpack.h"
#include <math.h>

void choldc(double **a, int n, double p[])
{
    void nrerror(char error_text[]);
    int i,j,k;
    double sum;

    for (i=1;i<=n;i++) {
        for (j=i;j<=n;j++) {
            for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
            if (i == j) {
                if (sum <= 0.0)
                    nrerror("choldc failed");
                p[i]=sqrt(sum);
            } else a[j][i]=sum/p[i];
        }
    }
}

void cholsl(double **a, int n, double p[], double b[], double x[])
{
    int i,k;
    double sum;

    for (i=1;i<=n;i++) {
        for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
        x[i]=sum/p[i];
    }
    for (i=n;i>=1;i--) {
        for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
        x[i]=sum/p[i];
    }
}

void lubksb(double **a, int n, int *indx, double b[])
{
    int i,ii=0,ip,j;
    double sum;

    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}

//#include <math.h>
//#define NRANSI
//#include "nrutil.h"
#define TINY 1.0e-20

void ludcmp(double **a, int n, int *indx, double *d)
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv;

    vv=dvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,n);
}
#undef TINY
//#undef NRANSI


