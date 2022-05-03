
#include "../math/numrec.h"
#include "../math/minpack.h"

// LevMarq method
void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
            int ma, double **covar, double **alpha, double *chisq,
            void (*funcs)(double, double [], double *, double [], int), double *alamda)
{
    void covsrt(double **covar, int ma, int ia[], int mfit);
    void gaussj(double **a, int n, double **b, int m);
    void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
                int ia[], int ma, double **alpha, double beta[], double *chisq,
                void (*funcs)(double, double [], double *, double [], int));
    int j,k,l;
    static int mfit;
    static double ochisq,*atry,*beta,*da,**oneda;
    
    if (*alamda < 0.0) {
        atry=dvector(1,ma);
        beta=dvector(1,ma);
        da=dvector(1,ma);
        for (mfit=0,j=1;j<=ma;j++)
            if (ia[j]) mfit++;
        oneda=dmatrix(1,mfit,1,1);
        *alamda=0.001;
        mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
        ochisq=(*chisq);
        for (j=1;j<=ma;j++) atry[j]=a[j];
    }
    for (j=1;j<=mfit;j++) {
        for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
        covar[j][j]=alpha[j][j]*(1.0+(*alamda));
        oneda[j][1]=beta[j];
    }
    gaussj(covar,mfit,oneda,1);
    for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
    if (*alamda == 0.0) {
        covsrt(covar,ma,ia,mfit);
        covsrt(alpha,ma,ia,mfit);
        free_dmatrix(oneda,1,mfit,1,1);
        free_dvector(da,1,ma);
        free_dvector(beta,1,ma);
        free_dvector(atry,1,ma);
        return;
    }
    for (j=0,l=1;l<=ma;l++)
        if (ia[l]) atry[l]=a[l]+da[++j];
    mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
    if (*chisq < ochisq) {
        *alamda *= 0.1;
        ochisq=(*chisq);
        for (j=1;j<=mfit;j++) {
            for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
            beta[j]=da[j];
        }
        for (l=1;l<=ma;l++) a[l]=atry[l];
    } else {
        *alamda *= 10.0;
        *chisq=ochisq;
    }
}

void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
            int ma, double **alpha, double beta[], double *chisq,
            void (*funcs)(double, double [], double *, double [], int))
{
    int i,j,k,l,m,mfit=0;
    double ymod,wt,sig2i,dy,*dyda;
    
    dyda=dvector(1,ma);
    for (j=1;j<=ma;j++)
        if (ia[j]) mfit++;
    for (j=1;j<=mfit;j++) {
        for (k=1;k<=j;k++) alpha[j][k]=0.0;
        beta[j]=0.0;
    }
    *chisq=0.0;
    for (i=1;i<=ndata;i++) {
        (*funcs)(x[i],a,&ymod,dyda,ma);
        sig2i=1.0/(sig[i]*sig[i]);
        dy=y[i]-ymod;
        for (j=0,l=1;l<=ma;l++) {
            if (ia[l]) {
                wt=dyda[l]*sig2i;
                for (j++,k=0,m=1;m<=l;m++)
                    if (ia[m]) alpha[j][++k] += wt*dyda[m];
                beta[j] += dy*wt;
            }
        }
        *chisq += dy*dy*sig2i;
    }
    for (j=2;j<=mfit;j++)
        for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
    free_dvector(dyda,1,ma);
}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
void covsrt(double **covar, int ma, int ia[], int mfit)
{
    int i,j,k;
    double swap;
    
    for (i=mfit+1;i<=ma;i++)
        for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
    k=mfit;
    for (j=ma;j>=1;j--) {
        if (ia[j]) {
            for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
                for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
                    k--;
        }
    }
}
#undef SWAP



// BEGIN :: AMOEBA METHOD

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define TINY 1.0e-10
#define NMAX 5000
#define GET_PSUM \
for (j=1;j<=ndim;j++) {\
    for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
        psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

// minimize in N-dimensions by downhill simplex method
void amoeba(double **p, double y[], int ndim, double ftol,
            double (*funk)(double []), int *nfunk)
{
    double amotry(double **p, double y[], double psum[], int ndim,
                  double (*funk)(double []), int ihi, double fac);
    int i,ihi,ilo,inhi,j,mpts=ndim+1;
    double rtol,sum,swap,ysave,ytry,*psum;
    
    psum=dvector(1,ndim);
    *nfunk=0;
    GET_PSUM
    for (;;) {
        ilo=1;
        ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
        for (i=1;i<=mpts;i++) {
            if (y[i] <= y[ilo]) ilo=i;
            if (y[i] > y[ihi]) {
                inhi=ihi;
                ihi=i;
            } else if (y[i] > y[inhi] && i != ihi) inhi=i;
        }
        rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        if (rtol < ftol) {
            SWAP(y[1],y[ilo])
            for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
                break;
        }
        if (*nfunk >= NMAX) nrerror("NMAX exceeded");
        *nfunk += 2;
        ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
        if (ytry <= y[ilo])
            ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi]) {
            ysave=y[ihi];
            ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave) {
                for (i=1;i<=mpts;i++) {
                    if (i != ilo) {
                        for (j=1;j<=ndim;j++)
                            p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                        y[i]=(*funk)(psum);
                    }
                }
                *nfunk += ndim;
                GET_PSUM
            }
        } else --(*nfunk);
    }
    free_dvector(psum,1,ndim);
}
#undef TINY
#undef SWAP
#undef GET_PSUM
#undef NMAX
//#undef NRANSI

//#define NRANSI
//#include "nrutil_nransi.h"

double amotry(double **p, double y[], double psum[], int ndim,
              double (*funk)(double []), int ihi, double fac)
{
    int j;
    double fac1,fac2,ytry,*ptry;
    
    ptry=dvector(1,ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=(*funk)(ptry);
    if (ytry < y[ihi]) {
        y[ihi]=ytry;
        for (j=1;j<=ndim;j++) {
            psum[j] += ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
    free_dvector(ptry,1,ndim);
    return ytry;
}
//#undef NRANSI

// END :: AMOEBA METHOD


// BEGIN :: POWELL METHOD

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define TINY 1.0e-25
#define ITMAX 200

// minimize in N-dimensions by Powell's method
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
            double (*func)(double []))
{
    void linmin(double p[], double xi[], int n, double *fret,
                double (*func)(double []));
    int i,ibig,j;
    double del,fp,fptt,t,*pt,*ptt,*xit;
    
    pt=dvector(1,n);
    ptt=dvector(1,n);
    xit=dvector(1,n);
    *fret=(*func)(p);
    for (j=1;j<=n;j++) pt[j]=p[j];
    for (*iter=1;;++(*iter)) {
        fp=(*fret);
        ibig=0;
        del=0.0;
        for (i=1;i<=n;i++) {
            for (j=1;j<=n;j++) xit[j]=xi[j][i];
            fptt=(*fret);
            linmin(p,xit,n,fret,func);
            if (fptt-(*fret) > del) {
                del=fptt-(*fret);
                ibig=i;
            }
        }
        if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
            free_dvector(xit,1,n);
            free_dvector(ptt,1,n);
            free_dvector(pt,1,n);
            return;
        }
        if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
        for (j=1;j<=n;j++) {
            ptt[j]=2.0*p[j]-pt[j];
            xit[j]=p[j]-pt[j];
            pt[j]=p[j];
        }
        fptt=(*func)(ptt);
        if (fptt < fp) {
            t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
            if (t < 0.0) {
                linmin(p,xit,n,fret,func);
                for (j=1;j<=n;j++) {
                    xi[j][ibig]=xi[j][n];
                    xi[j][n]=xit[j];
                }
            }
        }
    }
}
#undef TINY
#undef ITMAX
//#undef NRANSI

//#define NRANSI
//#include "nrutil_nransi.h"
#define TOL 2.0e-4

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
    double brent(double ax, double bx, double cx,
                 double (*f)(double), double tol, double *xmin);
    double f1dim(double x);
    void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
                double *fc, double (*func)(double));
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;
    
    ncom=n;
    pcom=dvector(1,n);
    xicom=dvector(1,n);
    nrfunc=func;
    for (j=1;j<=n;j++) {
        pcom[j]=p[j];
        xicom[j]=xi[j];
    }
    ax=0.0;
    xx=1.0;
    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
    *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
    for (j=1;j<=n;j++) {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_dvector(xicom,1,n);
    free_dvector(pcom,1,n);
}
#undef TOL
//#undef NRANSI

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
             double *xmin)
{
    int iter;
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0.0;
    
    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(x);
    for (iter=1;iter<=ITMAX;iter++) {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
        if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1) {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            else {
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN(tol1,xm-x);
            }
        } else {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
        fu=(*f)(u);
        if (fu <= fx) {
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
        } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu <= fv || v == x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }
    nrerror("Too many iterations in brent");
    *xmin=x;
    return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
//#undef NRANSI

//#define NRANSI
//#include "nrutil_nransi.h"

extern int ncom;
extern double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(double x)
{
    int j;
    double f,*xt;
    
    xt=dvector(1,ncom);
    for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    f=(*nrfunc)(xt);
    free_dvector(xt,1,ncom);
    return f;
}
//#undef NRANSI


//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double))
{
    double ulim,u,r,q,fu,dum;
    
    *fa=(*func)(*ax);
    *fb=(*func)(*bx);
    if (*fb > *fa) {
        SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(*func)(*cx);
    while (*fb > *fc) {
        r=(*bx-*ax)*(*fb-*fc);
        q=(*bx-*cx)*(*fb-*fa);
        u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
        (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
        ulim=(*bx)+GLIMIT*(*cx-*bx);
        if ((*bx-u)*(u-*cx) > 0.0) {
            fu=(*func)(u);
            if (fu < *fc) {
                *ax=(*bx);
                *bx=u;
                *fa=(*fb);
                *fb=fu;
                return;
            } else if (fu > *fb) {
                *cx=u;
                *fc=fu;
                return;
            }
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(*func)(u);
        } else if ((*cx-u)*(u-ulim) > 0.0) {
            fu=(*func)(u);
            if (fu < *fc) {
                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                SHFT(*fb,*fc,fu,(*func)(u))
            }
        } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
            u=ulim;
            fu=(*func)(u);
        } else {
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(*func)(u);
        }
        SHFT(*ax,*bx,*cx,u)
        SHFT(*fa,*fb,*fc,fu)
    }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
//#undef NRANSI

// END :: POWELL METHOD


// BEGIN :: FRPRMN METHOD
//frprmn minimize in N-dimensions by conjugate gradient

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
            double (*func)(double []), void (*dfunc)(double [], double []))
{
    void linmin(double p[], double xi[], int n, double *fret,
                double (*func)(double []));
    int j,its;
    double gg,gam,fp,dgg;
    double *g,*h,*xi;
    
    g=dvector(1,n);
    h=dvector(1,n);
    xi=dvector(1,n);
    fp=(*func)(p);
    (*dfunc)(p,xi);
    for (j=1;j<=n;j++) {
        g[j] = -xi[j];
        xi[j]=h[j]=g[j];
    }
    for (its=1;its<=ITMAX;its++) {
        *iter=its;
        linmin(p,xi,n,fret,func);
        if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
            FREEALL
            return;
        }
        fp= *fret;
        (*dfunc)(p,xi);
        dgg=gg=0.0;
        for (j=1;j<=n;j++) {
            gg += g[j]*g[j];
            dgg += (xi[j]+g[j])*xi[j];
        }
        if (gg == 0.0) {
            FREEALL
            return;
        }
        gam=dgg/gg;
        for (j=1;j<=n;j++) {
            g[j] = -xi[j];
            xi[j]=h[j]=g[j]+gam*h[j];
        }
    }
    nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
//#undef NRANSI

// END :: FRPRMN METHOD


// BEGIN :: DFPMIN METHOD

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_dvector(xi,1,n);free_dvector(pnew,1,n); \
free_dmatrix(hessin,1,n,1,n);free_dvector(hdg,1,n);free_dvector(g,1,n); \
free_dvector(dg,1,n);

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
            double(*func)(double []), void (*dfunc)(double [], double []))
{
    void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
                double *f, double stpmax, int *check, double (*func)(double []));
    int check,i,its,j;
    double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
    double *dg,*g,*hdg,**hessin,*pnew,*xi;
    
    dg=dvector(1,n);
    g=dvector(1,n);
    hdg=dvector(1,n);
    hessin=dmatrix(1,n,1,n);
    pnew=dvector(1,n);
    xi=dvector(1,n);
    fp=(*func)(p);
    (*dfunc)(p,g);
    for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) hessin[i][j]=0.0;
        hessin[i][i]=1.0;
        xi[i] = -g[i];
        sum += p[i]*p[i];
    }
    stpmax=STPMX*FMAX(sqrt(sum),(double)n);
    for (its=1;its<=ITMAX;its++) {
        *iter=its;
        lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
        fp = *fret;
        for (i=1;i<=n;i++) {
            xi[i]=pnew[i]-p[i];
            p[i]=pnew[i];
        }
        test=0.0;
        for (i=1;i<=n;i++) {
            temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
            if (temp > test) test=temp;
        }
        if (test < TOLX) {
            FREEALL
            return;
        }
        for (i=1;i<=n;i++) dg[i]=g[i];
        (*dfunc)(p,g);
        test=0.0;
        den=FMAX(*fret,1.0);
        for (i=1;i<=n;i++) {
            temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
            if (temp > test) test=temp;
        }
        if (test < gtol) {
            FREEALL
            return;
        }
        for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
        for (i=1;i<=n;i++) {
            hdg[i]=0.0;
            for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
        }
        fac=fae=sumdg=sumxi=0.0;
        for (i=1;i<=n;i++) {
            fac += dg[i]*xi[i];
            fae += dg[i]*hdg[i];
            sumdg += SQR(dg[i]);
            sumxi += SQR(xi[i]);
        }
        if (fac > sqrt(EPS*sumdg*sumxi)) {
            fac=1.0/fac;
            fad=1.0/fae;
            for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
            for (i=1;i<=n;i++) {
                for (j=i;j<=n;j++) {
                    hessin[i][j] += fac*xi[i]*xi[j]
                    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
                    hessin[j][i]=hessin[i][j];
                }
            }
        }
        for (i=1;i<=n;i++) {
            xi[i]=0.0;
            for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
        }
    }
    nrerror("too many iterations in dfpmin");
    FREEALL
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
//#undef NRANSI

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
            double *f, double stpmax, int *check, double (*func)(double []))
{
    int i;
    double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;
    
    *check=0;
    for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
    sum=sqrt(sum);
    if (sum > stpmax)
        for (i=1;i<=n;i++) p[i] *= stpmax/sum;
    for (slope=0.0,i=1;i<=n;i++)
        slope += g[i]*p[i];
    if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
    test=0.0;
    for (i=1;i<=n;i++) {
        temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
        if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;) {
        for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
        *f=(*func)(x);
        if (alam < alamin) {
            for (i=1;i<=n;i++) x[i]=xold[i];
            *check=1;
            return;
        } else if (*f <= fold+ALF*alam*slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope/(2.0*(*f-fold-slope));
            else {
                rhs1 = *f-fold-alam*slope;
                rhs2=f2-fold-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else {
                    disc=b*b-3.0*a*slope;
                    if (disc < 0.0) tmplam=0.5*alam;
                    else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
                    else tmplam=-slope/(b+sqrt(disc));
                }
                if (tmplam > 0.5*alam)
                    tmplam=0.5*alam;
            }
        }
        alam2=alam;
        f2 = *f;
        alam=FMAX(tmplam,0.1*alam);
    }
}
#undef ALF
#undef TOLX
//#undef NRANSI

// END :: DFPMIN METHOD


// BEGIN :: AMEBSA METHOD

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define GET_PSUM \
for (n=1;n<=ndim;n++) {\
    for (sum=0.0,m=1;m<=mpts;m++) sum += p[m][n];\
        psum[n]=sum;}
extern long idum;
double tt;

// amebsa simulated annealing in continuous spaces
void amebsa(double **p, double y[], int ndim, double pb[], double *yb, double ftol,
            double (*funk)(double []), int *iter, double temptr)
{
    double amotsa(double **p, double y[], double psum[], int ndim, double pb[],
                  double *yb, double (*funk)(double []), int ihi, double *yhi, double fac);
    double ran1(long *idum);
    int i,ihi,ilo,j,m,n,mpts=ndim+1;
    double rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum;
    
    psum=dvector(1,ndim);
    tt = -temptr;
    GET_PSUM
    for (;;) {
        ilo=1;
        ihi=2;
        ynhi=ylo=y[1]+tt*log(ran1(&idum));
        yhi=y[2]+tt*log(ran1(&idum));
        if (ylo > yhi) {
            ihi=1;
            ilo=2;
            ynhi=yhi;
            yhi=ylo;
            ylo=ynhi;
        }
        for (i=3;i<=mpts;i++) {
            yt=y[i]+tt*log(ran1(&idum));
            if (yt <= ylo) {
                ilo=i;
                ylo=yt;
            }
            if (yt > yhi) {
                ynhi=yhi;
                ihi=i;
                yhi=yt;
            } else if (yt > ynhi) {
                ynhi=yt;
            }
        }
        rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
        if (rtol < ftol || *iter < 0) {
            swap=y[1];
            y[1]=y[ilo];
            y[ilo]=swap;
            for (n=1;n<=ndim;n++) {
                swap=p[1][n];
                p[1][n]=p[ilo][n];
                p[ilo][n]=swap;
            }
            break;
        }
        *iter -= 2;
        ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,-1.0);
        if (ytry <= ylo) {
            ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,2.0);
        } else if (ytry >= ynhi) {
            ysave=yhi;
            ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,0.5);
            if (ytry >= ysave) {
                for (i=1;i<=mpts;i++) {
                    if (i != ilo) {
                        for (j=1;j<=ndim;j++) {
                            psum[j]=0.5*(p[i][j]+p[ilo][j]);
                            p[i][j]=psum[j];
                        }
                        y[i]=(*funk)(psum);
                    }
                }
                *iter -= ndim;
                GET_PSUM
            }
        } else ++(*iter);
    }
    free_dvector(psum,1,ndim);
}
#undef GET_PSUM
//#undef NRANSI

//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"

extern long idum;
extern double tt;

double amotsa(double **p, double y[], double psum[], int ndim, double pb[],
              double *yb, double (*funk)(double []), int ihi, double *yhi, double fac)
{
    double ran1(long *idum);
    int j;
    double fac1,fac2,yflu,ytry,*ptry;
    
    ptry=dvector(1,ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=1;j<=ndim;j++)
        ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=(*funk)(ptry);
    if (ytry <= *yb) {
        for (j=1;j<=ndim;j++) pb[j]=ptry[j];
        *yb=ytry;
    }
    yflu=ytry-tt*log(ran1(&idum));
    if (yflu < *yhi) {
        y[ihi]=ytry;
        *yhi=yflu;
        for (j=1;j<=ndim;j++) {
            psum[j] += ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
    free_dvector(ptry,1,ndim);
    return yflu;
}
//#undef NRANSI

// END :: AMEBSA METHOD
