
/* Driver for routine trapzd */

#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr_nransi.h"
#include "nrutil_nransi.h"
#include "func3.h"

#define NMAX 14
#define PIO2 1.5707963

/* Test function */
float func(float x)
{
    return (x*x)*(x*x-2.0)*sin(x);
}

float func2(float x)
{
    return (x*x-2.0)*sin(x);
}

/* Integral of test function */
float fint(float x)
{
	return 4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*(x*x)+28.0)*cos(x);
}


int main(void)
{
	int i;
	float a=0.0,b=PIO2,s;

	printf("\nIntegral of func with 2^(n-1) points\n");
	printf("Actual value of integral is %10.6f\n",fint(b)-fint(a));
	printf("%6s %24s\n","n","approx. integral");
    for (i=1;i<=NMAX;i++) {
        s=trapzd(func,a,b,i);
        printf("%6d %20.6f\n",i,s);
    }
    printf("\n\n Func2 \n");
    for (i=1;i<=NMAX;i++) {
        s=trapzd(func2,a,b,i);
        printf("%6d %20.6f\n",i,s);
    }
    
   NP=5;
    uno=2.0;

    printf("\n\n Func3 \n");
    for (i=1;i<=NMAX;i++) {
        s=trapzd(func3,a,b,i);
        printf("%6d %20.6f\n",i,s);
    }
	return 0;
}
#undef NRANSI
