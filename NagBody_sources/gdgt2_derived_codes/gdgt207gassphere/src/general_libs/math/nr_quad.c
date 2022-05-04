
#include <math.h>
#include "numrec.h"
#include "nr_quad.h"

/* ########################## QGAUS ############################################ */

float qgaus(float (*func)(float), float a, float b)
{
	int j;
	float xr,xm,dx,s;
	static float x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static float w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
	}
	return s *= xr;
}

/* ########################################################################### */

/* ########################## QGAUS_01 ####################################### */

float qgaus_01(float (*func)(float), float a, float b)
{
	int j;
	float xr,xm,dx,s;
	static float x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static float w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
	}
	return s *= xr;
}

/* ########################################################################### */


/* ############################ QROMB ######################################## */

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

float qromb(float (*func)(float), float a, float b)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	float ss,dss;
	float s[JMAXP+1],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

/* ########################################################################### */


/* ############################ QROMO ######################################## */

#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

float qromo(float (*func)(float), float a, float b,
	float (*choose)(float(*)(float), float, float, int))
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void nrerror(char error_text[]);
	int j;
	float ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

/* ########################################################################### */

/* ############################ QROMO_01 ######################################## */

#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

float qromo_01(float (*func)(float), float a, float b,
	float (*choose)(float(*)(float), float, float, int))
{
	void polint_01(float xa[], float ya[], int n, float x, float *y, float *dy);
	void nrerror(char error_text[]);
	int j;
	float ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint_01(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

/* ########################################################################### */


/* ############################# QSIMP ####################################### */

#define EPS 1.0e-6
#define JMAX 20

float qsimp(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	int j;
	float s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	nrerror("Too many steps in routine qsimp");
	return 0.0;
}
#undef EPS
#undef JMAX

/* ########################################################################### */


/* ############################# QTRAP ####################################### */

#define EPS 1.0e-5
#define JMAX 20

float qtrap(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	int j;
	float s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	nrerror("Too many steps in routine qtrap");
	return 0.0;
}
#undef EPS
#undef JMAX

/* ########################################################################### */


/* ############################## QUAD3D ##################################### */

/* QUAD3D est'a desactiviado ya que falta ver como
se definen las funciones yy1, yy2, z1, y z2. */

/*
static float xsav,ysav;
static float (*nrfunc)(float,float,float);

float quad3d(float (*func)(float, float, float), float x1, float x2)
{
	float qgaus(float (*func)(float), float a, float b);
	float f1(float x);

	nrfunc=func;
	return qgaus(f1,x1,x2);
}

float f1(float x)
{
	float qgaus(float (*func)(float), float a, float b);
	float f2(float y);
	float yy1(float),yy2(float);

	xsav=x;
	return qgaus(f2,yy1(x),yy2(x));
}

float f2(float y)
{
	float qgaus(float (*func)(float), float a, float b);
	float f3(float z);
	float z1(float,float),z2(float,float);

	ysav=y;
	return qgaus(f3,z1(xsav,y),z2(xsav,y));
}

float f3(float z)
{
	return (*nrfunc)(xsav,ysav,z);
}
*/

/* ########################################################################### */


/* ############################ POLINT ####################################### */

#define NRANSI

void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=nr_vector(1,n);
	d=nr_vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}
#undef NRANSI

/* ########################################################################### */

/* ############################ POLINT_01 ####################################### */

#define NRANSI

void polint_01(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=nr_vector(1,n);
	d=nr_vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}
#undef NRANSI

/* ########################################################################### */


/* ########################## TRAPZD ######################################### */

#define FUNC(x) ((*func)(x))

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDEXP ######################################### */

#define FUNC(x) ((*funk)(-log(x))/(x))

float midexp(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=exp(-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDEXP_01 ######################################### */

#define FUNC(x) ((*funk)(-log(x))/(x))

float midexp_01(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=exp(-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDINF ######################################### */

#define FUNC(x) ((*funk)(1.0/(x))/((x)*(x)))

float midinf(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,b,a;
	static float s;
	int it,j;

	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}
#undef FUNC

/* ########################################################################### */

/* ########################## MIDINF_01 ######################################### */

#define FUNC(x) ((*funk)(1.0/(x))/((x)*(x)))

float midinf_01(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,b,a;
	static float s;
	int it,j;

	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDPNT ######################################### */

#define FUNC(x) ((*func)(x))

float midpnt(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del,ddel;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */

/* ########################## MIDPNT_01 ######################################### */

#define FUNC(x) ((*func)(x))

float midpnt_01(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del,ddel;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDSQL ######################################### */

#define FUNC(x) (2.0*(x)*(*funk)(aa+(x)*(x)))

float midsql(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDSQL_01 ######################################### */

#define FUNC(x) (2.0*(x)*(*funk)(aa+(x)*(x)))

float midsql_01(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDSQU ######################################### */

#define FUNC(x) (2.0*(x)*(*funk)(bb-(x)*(x)))

float midsqu(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */


/* ########################## MIDSQU_01 ######################################### */

#define FUNC(x) (2.0*(x)*(*funk)(bb-(x)*(x)))

float midsqu_01(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

/* ########################################################################### */
