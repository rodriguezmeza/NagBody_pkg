/*==============================================================================
    MODULE: nminchi2.c				[nchi2]
    Written by: Mario A. Rodriguez-Meza
    Starting date:	January 2018
    Purpose:
    Language: C
    Use:
    Routines and functions:
    External modules, routines and headers:
    Comments and notes:
    Info: Mario A. Rodriguez-Meza
    Depto. de Fisica, ININ
        Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
        e-mail: marioalberto.rodriguez@inin.gob.mx
        https://github.com/rodriguezmeza

    Major revisions:
    Copyright: (c) 2005-2018 Mar.  All Rights Reserved
================================================================================
    Legal matters:
    The author does not warrant that the program and routines it contains
    listed below are free from error or suitable for particular applications,
    and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"


global void grid_Chi2(real (*ymodel)(real, real *))
{
    real Chi2tmp;
    real Chi2redtmp;
    int *N;
    real *p1, *p2, *Chi2red;
    int i, j;
    real *dp;
    int pN;
    real **mp;
    int *ic;

//    if (nObs(pObs[gd.filecount])<=pd.nparams) {
    if (pObsT.no<=pd.nparams) {
        error("\nNumber of observed points is equal or less than the number of parameters!\n\n");
    }

    N = nr_ivector(1, pd.nparams);
    pN = 1;
    for (i=1; i<=pd.nparams; i++) {
        N[i] = cmd.gridN;
        pN *= N[i];
    }
    
    dp = dvector(1,pd.nparams);
    ic = nr_ivector(1, pd.nparams);

    for (i=1; i<=pd.nparams; i++) {
        dp[i] = (pd.rparams[i][2]-pd.rparams[i][1])/((real)(N[i]-1));
        fprintf(stdout,"Gridding parameters space: %d %d %g\n",i, N[i], dp[i]);
    }

    fprintf(stdout,"Product of N's: %d\n",pN);
    mp = dmatrix(1, pN, 1, pd.nparams + 2);

    int *ppN;
    int *ib;
    ib = nr_ivector(1, pd.nparams);
    ppN = nr_ivector(1, pd.nparams);
    real chi2redmin, *pmin;
    pmin = nr_dvector(1, pd.nparams);
    
    chi2redmin = 1.0e10;

    ppN[pd.nparams] = N[pd.nparams];
    ib[pd.nparams] = 1.0;
    for (j=pd.nparams-1; j>=1; j--) {
        ib[j] = 1.0;
        ppN[j] *= N[j];
    }
    for (i=1; i<=pN; i++) {
        if (i > ib[pd.nparams]*N[pd.nparams]) {
            ++ib[pd.nparams];
            ic[pd.nparams] = 1;
            pd.params[pd.nparams] = pd.rparams[pd.nparams][1] + dp[pd.nparams]*((real)(ic[pd.nparams]-1));
        } else {
            ic[pd.nparams] = i - (ib[pd.nparams]-1)*N[pd.nparams];
            pd.params[pd.nparams] = pd.rparams[pd.nparams][1] + dp[pd.nparams]*((real)(ic[pd.nparams]-1));
        }
        for (j=pd.nparams-1; j>=1; j--){
            if (ib[j+1] > ib[j]*N[j]) {
                ++ib[j];
                ic[j] = 1;
                pd.params[j] = pd.rparams[j][1] + dp[j]*((real)(ic[j]-1));
            } else {
                ic[j] = ib[j+1] - (ib[j]-1)*N[j];
                pd.params[j] = pd.rparams[j][1] + dp[j]*((real)(ic[j]-1));
            }
        }
        Chi2tmp = gd.modelChi2(ymodel, pd.params);
//        Chi2redtmp = Chi2tmp/(nObs(pObs[gd.filecount])-pd.nparams);
        Chi2redtmp = Chi2tmp/(pObsT.no-pd.nparams);
        mp[i][1] = (real)i;
        mp[i][2] = Chi2redtmp;
        for (j=3; j<=pd.nparams+2; j++){
            mp[i][j] =  pd.params[j-2];
        }
        
        if (chi2redmin > Chi2redtmp) {
            chi2redmin = Chi2redtmp;
            for (j=1; j<=pd.nparams; j++){
                pmin[j] =  pd.params[j];
            }
//            for (j=1; j<=pd.nparams; j++) {
//                fprintf(stdout,"%e ",pd.params[j]);
//            }
//            fprintf(stdout," chi2redmin: %e\n",chi2redmin);
        }
//        for (j=1; j<=pd.nparams; j++) {
//            fprintf(stdout,"%e ",pd.params[j]);
//        }
//        fprintf(stdout," chi2redtmp: %e\n",Chi2redtmp);
    }

    fprintf(stdout,"\nSorting %d Chi2red values",pN);
    for (j=1; j<=pd.nparams; j++) {
        pd.params[j] = pmin[j];
    }
    fprintf(stdout,"\nMinimum Chi2red value and its parameters: %e ", chi2redmin);
    for (j=1; j<=pd.nparams; j++) {
        fprintf(stdout,"%e ",pd.params[j]);
    }
    fprintf(stdout,"\n\n");

    free_dmatrix(mp, 1, pN, 1, pd.nparams + 2);
    free_ivector(N, 1, pd.nparams);
    free_ivector(ic, 1, pd.nparams);
    free_dvector(dp,1, pd.nparams);
}

// Not working for combined observations (more than one file)
global void min_Chi2_LevMarq(void (*funcs)(double, double [], double *, double [], int))
{
    int i,*ia,iter,itst,j,k;
    double alamda,chisq,ochisq,**covar,**alpha;

    ia=nr_ivector(1,pd.nparams);
    covar=dmatrix(1,pd.nparams,1,pd.nparams);
    alpha=dmatrix(1,pd.nparams,1,pd.nparams);

    for (i=1;i<=pd.nparams;i++) ia[i]=1;           // None parameter is fixed...

    fprintf(stdout,"\nmin_Chi2_LevMarq:: Initial parameters:\n");
    for (i=1;i<=pd.nparams;i++) fprintf(stdout,"%g ",pd.params[i]);
    fprintf(stdout,"\n");
    
    for (iter=1;iter<=1;iter++) {
        alamda = -1;
        mrqmin(xObs(pObs[gd.filecount]),yObs(pObs[gd.filecount]),
               sigmaObs(pObs[gd.filecount]),nObs(pObs[gd.filecount]),
               pd.params,ia,pd.nparams,covar,alpha,&chisq,funcs,&alamda);

        k=1;
        itst=0;
        for (;;) {
            printf("\n%s %2d %17s %10.4f %10s %9.2e\n","Iteration #",k,
                   "chi-squared:",chisq,"alamda:",alamda);
            printf("%8s %8s\n",
                   "a[1]","a[2],...");
            for (i=1;i<=pd.nparams;i++) printf("%g ",pd.params[i]);
            printf("\n");
            k++;
            ochisq=chisq;
            mrqmin(xObs(pObs[gd.filecount]),yObs(pObs[gd.filecount]),
                   sigmaObs(pObs[gd.filecount]),nObs(pObs[gd.filecount]),
                   pd.params,ia,pd.nparams,covar,alpha,&chisq,funcs,&alamda);

            if (chisq > ochisq)
                itst=0;
            else if (fabs(ochisq-chisq) < 0.1)
                itst++;
            if (itst < 4) continue;
            alamda=0.0;
            mrqmin(xObs(pObs[gd.filecount]),yObs(pObs[gd.filecount]),
                   sigmaObs(pObs[gd.filecount]),nObs(pObs[gd.filecount]),
                   pd.params,ia,pd.nparams,covar,alpha,&chisq,funcs,&alamda);
            printf("\nUncertainties:\n");
            for (i=1;i<=pd.nparams;i++) printf("%g ",rsqrt(covar[i][i]));
            printf("\n");
            break;
        }
    }
    free_dmatrix(alpha,1,pd.nparams,1,pd.nparams);
    free_dmatrix(covar,1,pd.nparams,1,pd.nparams);
    free_ivector(ia,1,pd.nparams);
}


global real P_value(real Chi2, real nu)
{
    real Pvtmp;
    
    Pvtmp = gammp(nu/2.0, Chi2/2.0);
    
    return Pvtmp;
}

global real Q_value(real Chi2, real nu)
{
    real Qvtmp;
    
    Qvtmp = gammq(nu/2.0, Chi2/2.0);
    
    return Qvtmp;
}


// BEGIN :: AMOEBA METHOD

double func_amoeba(double x[])
{
    real functmp;

    functmp = gd.modelChi2(gd.modelfun, x);

    return functmp;
}

global void min_Chi2_amoeba(void (*funcs)(double, double [], double *, double [], int))
{
    int i,nfunc,j,ndim;
    double *x,*y,**p;
    int NP, MP;
    real *dp;
    int N;
    
    NP = pd.nparams;
    MP = NP+1;
    ndim = NP;
    
    printf("\namoeba method of minimization... \n");
    y=dvector(1,MP);
    p=dmatrix(1,MP,1,NP);
    
    dp = dvector(1,pd.nparams);
    N = cmd.gridNSimplex;

    for (i=1; i<=pd.nparams; i++) {
        dp[i] = (pd.rparams[i][2]-pd.rparams[i][1])/((real)(N-1));
        fprintf(stdout,"Amoeba:: gridding size: %d %g\n",i, dp[i]);
    }

    if (cmd.ucparams) {
        for (j=1;j<=NP;j++)
            pd.params[j]=(pd.rparams[j][1]+pd.rparams[j][2])/2.0;
    }

    for (i=1;i<=MP;i++) {
        for (j=1;j<=NP;j++)
        pd.params[j]=p[i][j]=(i == (j+1) ? pd.params[j]+dp[j] : pd.params[j]);
        y[i]=func_amoeba(pd.params);
    }

    printf("Vertices of initial 3-d simplex and\n");
    printf("function values at the vertices:\n\n");
    printf("%3s %10s %12s %12s\n\n",
           "i","param[1]","param[2]...","and Chi2");
    for (i=1;i<=MP;i++) {
        printf("%3d ",i);
        for (j=1;j<=NP;j++) printf("%e ",p[i][j]);
        printf("%e\n",y[i]);
    }
//
    amoeba(p,y,ndim,cmd.ftolmin,func_amoeba,&nfunc);
    printf("\nNumber of function evaluations: %3d\n",nfunc);
    printf("Vertices of final 3-d simplex and\n");
    printf("function values at the vertices:\n\n");
    printf("%3s %10s %12s %12s\n\n",
           "i","param[1]","param[2]...","and Chi2");
    for (i=1;i<=MP;i++) {
        printf("%3d ",i);
        for (j=1;j<=NP;j++) printf("%16.8e ",p[i][j]);
        printf("%16.8e\n",y[i]);
    }

    real ymin;
    ymin = y[1];
    for (j=1;j<=NP;j++)
        pd.params[j] = p[1][j];
    for (i=1;i<=MP;i++) {
        if (ymin > y[i]) {
            ymin = y[i];
            for (j=1;j<=NP;j++)
                pd.params[j] = p[i][j];
        }
    }
    printf("\nFinal: %10s %12s %12s\n\n",
           "param[1]","param[2]...","and Chi2");
        for (j=1;j<=NP;j++) printf("%16.8e ",pd.params[j]);
        printf("%16.8e\n",func_amoeba(pd.params));

    free_dvector(dp,1, pd.nparams);

    free_dmatrix(p,1,MP,1,NP);
    free_dvector(y,1,MP);
}

// END :: AMOEBA METHOD


// BEGIN :: POWELL METHOD

double func_powell(double x[])
{
    real functmp;
    functmp = gd.modelChi2(gd.modelfun, x);
    return functmp;
}

global void min_Chi2_powell(void (*funcs)(double, double [], double *, double [], int))
{
    int i,iter,j;
    double fret,**xi;

    printf("\nPowell method of minimization... \n");
    xi=dmatrix(1,pd.nparams,1,pd.nparams);
    for (i=1;i<=pd.nparams;i++)
        for (j=1;j<=pd.nparams;j++)
            xi[i][j]=(i == j ? 1.0 : 0.0);
    powell(pd.params,xi,pd.nparams,cmd.ftolmin,&iter,&fret,func_powell);
    printf("Iterations: %3d\n\n",iter);
    printf("Minimum found at: \n");
    for (i=1;i<=pd.nparams;i++) printf("%12.6f",pd.params[i]);
    printf("\n\nMinimum function value = %12.6f \n\n",fret);
    free_dmatrix(xi,1,pd.nparams,1,pd.nparams);
}

// END :: POWELL METHOD


// BEGIN :: FRPRMN METHOD

//frprmn minimize in N-dimensions by conjugate gradient

#define NDIM 3
#define PIO2 1.5707963

double func_frprmn(double x[])
{
    real functmp;
    functmp = gd.modelChi2(gd.modelfun, x);
    return functmp;
}

void dfunc_frprmn(double x[],double df[])
{
    real Chi2tmp;
    int i, j;
    real y;
    real *dfdp;

    dfdp=dvector(1,pd.nparams);
    set_model_Interp(gd.modelfun, pd.params);

    for (j=1; j<=pd.nparams; j++) {
        df[j] = 0.0;
        for (i=1; i<=nObs(pObs[gd.filecount]); i++) {
            if (sigmaObs(pObs[gd.filecount])[i]==0.)
                sigmaObs(pObs[gd.filecount])[i] = 1.0;
            gd.modelfunwd(xObs(pObs[gd.filecount])[i], pd.params, &y, dfdp, pd.nparams);
            y = model_Interp(xObs(pObs[gd.filecount])[i]);
            Chi2tmp += dfdp[j]*(yObs(pObs[gd.filecount])[i] - y)
                        /rsqr(sigmaObs(pObs[gd.filecount])[i]);
        }
        df[j] *= -2.0;
    }

    free_dvector(dfdp,1,pd.nparams);
}

global void min_Chi2_frprmn(void (*funcs)(double, double [], double *, double [], int))
{
    int iter,k;
    double angl,fret;
    int i;

    printf("\nConjugate gradient method of minimization... \n");
    frprmn(pd.params,pd.nparams,cmd.ftolmin,&iter,&fret,func_frprmn,dfunc_frprmn);
        printf("Iterations: %3d\n",iter);
    printf("Minimum found at: \n");
    for (i=1;i<=pd.nparams;i++) printf("%12.6f",pd.params[i]);
        printf("\nFunc. value at solution %14f\n",fret);
}

#undef NDIM
#undef PIO2

// END :: FRPRMN METHOD


// BEGIN :: DFPMIN METHOD

// dfpmin minimize in N-dimensions by variable metric method

static int nfunc,ndfunc;

double func_dfpmin(double x[])
{
    real functmp;
    nfunc++;
    functmp = gd.modelChi2(gd.modelfun, x);
    return functmp;
}

void dfunc_dfpmin(double x[],double df[])
{
    real Chi2tmp;
    int i, j;
    real y;
    real *dfdp;

    ndfunc++;

    dfdp=dvector(1,pd.nparams);
    set_model_Interp(gd.modelfun, pd.params);

    for (j=1; j<=pd.nparams; j++) {
        df[j] = 0.0;
        for (i=1; i<=nObs(pObs[gd.filecount]); i++) {
            if (sigmaObs(pObs[gd.filecount])[i]==0.)
                sigmaObs(pObs[gd.filecount])[i] = 1.0;
            gd.modelfunwd(xObs(pObs[gd.filecount])[i], pd.params, &y, dfdp, pd.nparams);
            y = model_Interp(xObs(pObs[gd.filecount])[i]);
            Chi2tmp += dfdp[j]*(yObs(pObs[gd.filecount])[i] - y)
                            /rsqr(sigmaObs(pObs[gd.filecount])[i]);
        }
        df[j] *= -2.0;
    }

    free_dvector(dfdp,1,pd.nparams);
}

#define NDIM 2
#define GTOL 1.0e-4

global void min_Chi2_dfpmin(void (*funcs)(double, double [], double *, double [], int))
{
    int iter;
    double fret;
    int i;

    printf("\nQuasi-Newton method of minimization... \n");
    nfunc=ndfunc=0;
    dfpmin(pd.params,pd.nparams,GTOL,&iter,&fret,func_dfpmin,dfunc_dfpmin);
    printf("Iterations: %3d\n",iter);
    printf("Func. evals: %3d\n",nfunc);
    printf("Deriv. evals: %3d\n",ndfunc);
    printf("Minimum found at: \n");
    for (i=1;i<=pd.nparams;i++) printf("%12.6f",pd.params[i]);
    printf("\nFunc. value at solution %14.6g\n",fret);
}

#undef NDIM
#undef GTOL

// END :: DFPMIN METHOD


// BEGIN :: AMEBSA METHOD

#define NP 4
#define MP 5
#define N 4
#define RAD 0.3
#define AUG 2.0

long idum=(-64);

double tfunk(double p[])
{
    int j;
    double q,r,sumd=0.0,sumr=0.0;
    static double wid[N+1]={0.0,1.0,3.0,10.0,30.0};
    
    for (j=1;j<=N;j++) {
        q=p[j]*wid[j];
        r=(double)(q >= 0 ? (int)(q+0.5) : -(int)(0.5-q));
        sumr += q*q;
        sumd += (q-r)*(q-r);
    }
    return 1+sumr*(1+(sumd > RAD*RAD ? AUG : AUG*sumd/(RAD*RAD)));
}

global void min_Chi2_amebsa(void (*funcs)(double, double [], double *, double [], int))
{
    int i,iiter,iter,j,jiter,ndim=NP,nit;
    double temptr,yb,ybb;
    double **p,*x,*y,*pb;
    static double xoff[NP+1]={0.0,10.0,10.0,10.0,10.0};

    printf("\nSimulated annealing + simplex method of minimization... \n");
    p=dmatrix(1,MP,1,NP);
    x=dvector(1,NP);
    y=dvector(1,MP);
    pb=dvector(1,NP);
    for (i=1;i<=MP;i++)
        for (j=1;j<=NP;j++) p[i][j]=0.0;
    for (;;) {
        for (j=2;j<=MP;j++) p[j][j-1]=1.0;
        for (i=1;i<=MP;i++) {
            for (j=1;j<=NP;j++) x[j]=(p[i][j] += xoff[j]);
            y[i]=tfunk(x);
        }
        yb=1.0e30;
        printf("Input t, iiter:\n");
        if (scanf("%f %d",&temptr,&iiter) == EOF) break;
        ybb=1.0e30;
        nit=0;
        for (jiter=1;jiter<=100;jiter++) {
            iter=iiter;
            temptr *= 0.8;
            amebsa(p,y,ndim,pb,&yb,cmd.ftolmin,tfunk,&iter,temptr);
            nit += iiter-iter;
            if (yb < ybb) {
                ybb=yb;
                printf("%6d %10.3e ",nit,temptr);
                for (j=1;j<=NP;j++) printf("%10.5f ",pb[j]);
                printf("%15.7e\n",yb);
            }
            if (iter > 0) break;
        }
        printf("Vertices of final 3-D simplex and\n");
        printf("float values at the vertices:\n");
        printf("%3s %10s %12s %12s %14s\n\n",
               "i","x[i]","y[i]","z[i]","function");
        for (i=1;i<=MP;i++) {
            printf("%3d ",i);
            for (j=1;j<=NP;j++) printf("%12.6f ",p[i][j]);
            printf("%15.7e\n",y[i]);
        }
        printf("%3d ",99);
        for (j=1;j<=NP;j++) printf("%12.6f ",pb[j]);
        printf("%15.7e\n",yb);
    }
    free_dvector(pb,1,NP);
    free_dvector(y,1,MP);
    free_dvector(x,1,NP);
    free_dmatrix(p,1,MP,1,NP);
    printf("Normal completion\n");
}

#undef NP
#undef MP

#undef N
#undef RAD
#undef AUG

// END :: AMEBSA METHOD

