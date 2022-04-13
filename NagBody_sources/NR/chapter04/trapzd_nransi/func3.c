
/* Driver for routine trapzd */

#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr_nransi.h"
#include "nrutil_nransi.h"

// Funcion interpolada

//#define NP 10
#define PI 3.1415926

//#define NR_END 1
//#define FREE_ARG char*
/*
float *vector(long nl, long nh)
{
    float *v;

    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-NR_END));
}
*/

int NP;
float uno;


float func3(float x)
{
//    float ftmp;
    int i,nfunc;
    float f,y,yp1,ypn,*xa,*ya,*y2;

    xa=vector(1,NP);
    ya=vector(1,NP);
    y2=vector(1,NP);
    for (i=1;i<=NP;i++) {
        xa[i]=uno*i/NP;
        ya[i]=exp(xa[i]);
    }
    yp1=exp(xa[1]);
    ypn=exp(xa[NP]);
    spline(xa,ya,NP,yp1,ypn,y2);

    splint(xa,ya,y2,NP,x,&y);

    
    free_vector(y2,1,NP);
    free_vector(ya,1,NP);
    free_vector(xa,1,NP);

    
    return y;
}

//#undef NP
#undef PI

//



#undef NRANSI
