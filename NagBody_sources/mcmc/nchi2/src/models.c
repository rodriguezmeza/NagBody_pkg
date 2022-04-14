/*==============================================================================
	MODULE: models.c			[nchi2]
	Written by: M.A. Rodriguez-Meza
	Starting date: January 2018
	Purpose: Rutines to create several types of models
	Language: C
	Use: 'testmodel();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		https://github.com/rodriguezmeza

	Mayor revisions: January 22, 2018
	Copyright: (c) 2005-2020 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#define global
#include "globaldefs.h"

local void set_model_Interp_init(void);

// TEMPLATE FOR THE USER
#include "model_user.h"

local void model_string_to_int(string, int *);
local void set_model_Interp_general(real (*ymodel)(real, real *), real params[]);

#include "models/model_rcbecho.h"
#include "models/model_rcpiso.h"
#include "models/model_sncpl.h"
#include "models/model_baocpl.h"
#include "models/model_gseir_ci.h"

void set_model(void)
{
    int model_int;

    model_string_to_int(cmd.model, &model_int);
    switch (model_int){

    case RCBECHO: Model_RCBECHO(); break;
    case RCPISO: Model_RCPISO(); break;
    case GSEIRINFCUM: Model_gSEIR_CI(); break;
    case SNCPL: Model_SNCPL(); break;
    case BAOCPL: Model_BAOCPL(); break;
    case USERMODEL: set_Model_USER(); break; // In models_user.h

    default: error("\nUnknown model type %s\n\n",cmd.model);
    }
    
    fprintf(stdout,"\nNumber of parameters to fit: %d\n",pd.nparams);
}

void test_model(void)
{
    int model_int;

    model_string_to_int(cmd.model, &model_int);
    switch (model_int){

    case RCBECHO: test_RCBECHO(); break;
    case RCPISO: test_RCPISO(); break;
    case GSEIRINFCUM: test_gSEIR_CI(); break;
    case SNCPL: test_SNCPL(); break;
    case BAOCPL: test_BAOCPL(); break;
    case USERMODEL: test_USER(); break; // In models_user.h

    default: error("\nUnknown model type to test %s\n\n",cmd.model);
    }
}

global void set_model_Interp(real (*ymodel)(real, real *), real params[])
{
    int model_int;

    model_string_to_int(cmd.model, &model_int);
    switch (model_int){

    case RCBECHO: set_model_Interp_general(ymodel, params); break;
    case RCPISO: set_model_Interp_general(ymodel, params); break;
    case GSEIRINFCUM: set_model_Interp_gSEIR_CI(ymodel, params); break;
    case SNCPL: set_model_Interp_general(ymodel, params); break;
    case BAOCPL: set_model_Interp_general(ymodel, params); break;
    case USERMODEL: set_model_Interp_general(ymodel, params); break;

    default: error("\nUnknown model type to set interpolation %s\n\n",cmd.model);
    }
}

local void model_string_to_int(string model_str,int *model_int)
{
	*model_int = -1;
    if (strcmp(model_str,"RCBECHO") == 0)           *model_int=RCBECHO;
    if (strcmp(model_str,"rcbecho") == 0)           *model_int=RCBECHO;
    if (strcmp(model_str,"RCPISO") == 0)            *model_int=RCPISO;
    if (strcmp(model_str,"rcpiso") == 0)            *model_int=RCPISO;
    if (strcmp(model_str,"GSEIRINFCUM") == 0)       *model_int=GSEIRINFCUM;
    if (strcmp(model_str,"gSEIRInfCum") == 0)       *model_int=GSEIRINFCUM;
    if (strcmp(model_str,"gseirinfcum") == 0)       *model_int=GSEIRINFCUM;
    if (strcmp(model_str,"SNCPL") == 0)             *model_int=SNCPL;
    if (strcmp(model_str,"sncpl") == 0)             *model_int=SNCPL;
    if (strcmp(model_str,"BAOCPL") == 0)             *model_int=BAOCPL;
    if (strcmp(model_str,"baocpl") == 0)             *model_int=BAOCPL;
    if (strcmp(model_str,"USER") == 0)              *model_int=USERMODEL;
    if (strcmp(model_str,"user") == 0)              *model_int=USERMODEL;
}

global void Model_end(void)
{
    gd.modelend();

    free_dmatrix(pd.rparams,1,pd.nparams,1,2);
    free_dvector(pd.dparams,1,pd.nparams);
    free_dvector(pd.params,1,pd.nparams);

    free_dvector(gd.yInterp2,1,gd.nInterp);
    free_dvector(gd.yInterp,1,gd.nInterp);
    free_dvector(gd.xInterp,1,gd.nInterp);
}

local void set_model_Interp_init(void)
{
    gd.nInterp = nObs(pObs[gd.filecount]);
    gd.xInterp = dvector(1,gd.nInterp);
    gd.yInterp = dvector(1,gd.nInterp);
    gd.yInterp2 = dvector(1,gd.nInterp);
}

local void set_model_Interp_general(real (*ymodel)(real, real *), real params[])
{
    int i;

    for (i=1; i<gd.nInterp; i++) {
        gd.xInterp[i] = xObs(pObs[gd.filecount])[i];
        gd.yInterp[i] = ymodel(xObs(pObs[gd.filecount])[i], params);
    }
    spline(gd.xInterp,gd.yInterp,gd.nInterp,1.0e30,1.0e30,gd.yInterp2);
}

global real model_Interp(real x)
{
    real yftmp;

    if ( x < gd.xInterp[1] || x > gd.xInterp[gd.nInterp] )
        fprintf(gd.outlog,
                "\n\nmodelInterpolation: warning! :: x is out of range (xInterp[1], xInterp[nInterp])... %g %g %g\n",
                x, gd.xInterp[1], gd.xInterp[gd.nInterp]);

    splint(gd.xInterp,gd.yInterp,gd.yInterp2,gd.nInterp,x,&yftmp);

    return (yftmp);
}

/*
local void set_pd_to_md(void)
{
    md.nparams = pd.nparams;
    md.p1cv = pd.params[1];
    md.p1min = pd.rparams[1][1];
    md.p1max = pd.rparams[1][2];
    md.dp1 = pd.dparams[1];
    md.p2cv = pd.params[2];
    md.p2min = pd.rparams[2][1];
    md.p2max = pd.rparams[2][2];
    md.dp2 = pd.dparams[2];
}
*/


// Additional auxiliary routines or general model independent

// First derivatives of the model respect parameters
// Numerical fourth order computation
// x : x independent coord; a : params; y : model;
// na : number of parameters; dyda : derivatives
global void model_wdyda(double x, double a[],
                         double *y, double dyda[], int na)
{
    int i, j, ic;
    double fac,ex,arg;
    real **ma;
    real *va;

    va = dvector(1,na);
    ma = dmatrix(1,na*5,1,na);

    for (i=1; i<=5*na; i++)
        for (j=1; j<=na; j++)
            ma[i][j] = a[j];
    
    for (j=1; j<=na; j++) {
        ic = 5*(j-1);
        for (i=1; i<=5; i++) {
            ma[ic+i][j] += pd.dparams[j]*((real)(i-3));
        }
    }

    *y = gd.modelfun(x,a);

    for (i=1; i<=na; i++) {
        ic = 5*(i-1)+1;
        dyda[i] = 0.0;
        for (j=1; j<=na; j++) {
            va[j] = ma[ic][j];
        }
        dyda[i] += 1.0*gd.modelfun(x,va);
        for (j=1; j<=na; j++) {
            va[j] = ma[ic+1][j];
        }
        dyda[i] -= 8.0*gd.modelfun(x,va);
        for (j=1; j<=na; j++) {
            va[j] = ma[ic+2][j];
        }
        dyda[i] += 0.0*gd.modelfun(x,va);
        for (j=1; j<=na; j++) {
            va[j] = ma[ic+3][j];
        }
        dyda[i] += 8.0*gd.modelfun(x,va);
        for (j=1; j<=na; j++) {
            va[j] = ma[ic+4][j];
        }
        dyda[i] -= 1.0*gd.modelfun(x,va);
        dyda[i] /= (12.0*pd.dparams[i]);
/*
        dyda[i] = (
                    1.0*gSEIR_CI(x,ma[ic])
                    -8.0*gSEIR_CI(x,ma[ic+1])
                    +0.0*gSEIR_CI(x,ma[ic+2])
                    +8.0*gSEIR_CI(x,ma[ic+3])
                    -1.0*gSEIR_CI(x,ma[ic+4])
                    )/(12.0*da[i]);
 */
    }

    free_dmatrix(ma,1,na*5,1,na);
    free_dvector(va,1,na);
}


// First derivatives of the model respect parameters
// Numerical fourth order computation
// x : x independent coord; a : params; y : model;
// na : number of parameters; dyda : derivatives
global void model_wd1yda(double x, double a[],
                         double *y, double d1yda[], int na)
{
    int i, j, k, l;
    real *fv;
    real *va;

    fv = dvector(1,na);
    va = dvector(1,na);

    fv[1] = +1.0;
    fv[2] = -8.0;
    fv[3] = 0.0;
    fv[4] = +8.0;
    fv[5] = -1.0;

    for (i=1; i<=na; i++) {
        for (j=1; j<=na; j++) {
            va[j] = a[j];
        }
        va[i] -= 2.0*pd.dparams[i];
        d1yda[i] = 0.0;
        for (j=1; k<=5; k++) {
            d1yda[i] += fv[j]*gd.modelfun(x,va);
            va[i] += pd.dparams[i];
        }
        d1yda[i] /= (12.0*pd.dparams[i]);
    }

    free_dvector(va,1,na);
    free_dvector(fv,1,na);
}

// Second derivatives of the model respect parameters
// Numerical fourth order computation
// x : x independent coord; a : params; y : model;
// na : number of parameters; dyda : derivatives
global void model_wd2yda(double x, double a[],
                         double **d2yda, int na)
{
    int i, j, k, l;
    real *fv;
    real *fvd;
    real *va;

    fv = dvector(1,na);
    fvd = dvector(1,na);
    va = dvector(1,na);

    fv[1] = +1.0;
    fv[2] = -8.0;
    fv[3] = 0.0;
    fv[4] = +8.0;
    fv[5] = -1.0;

    fvd[1] = -1.0;
    fvd[2] = +16.0;
    fvd[3] = -30.0;
    fvd[4] = +16.0;
    fvd[5] = -1.0;

    for (i=1; i<=na; i++) {
        for (j=1; j<=na; j++) {
            va[j] = a[j];
        }
        va[i] -= 2.0*pd.dparams[i];
        for (j=1; j<=na; j++) {
            va[j] -= 2.0*pd.dparams[j];
            d2yda[i][j] = 0.0;
            if (i==j) {
                for (k=1; k<=5; k++) {
                    d2yda[i][i] += fv[k]*gd.modelfun(x,va);
                    va[i] += pd.dparams[i];
                }
                d2yda[i][i] /= (12.0*pd.dparams[i]*pd.dparams[i]);
            } else {
                for (k=1; k<=5; k++) {
                    for (l=1; l<=5; l++){
                        d2yda[i][j] += fv[l]*gd.modelfun(x,va);
                        va[j] += pd.dparams[j];
                    }
                    d2yda[i][j] *= fv[k];
                    va[i] += pd.dparams[i];
                }
                d2yda[i][j] /= (144.0*pd.dparams[i]*pd.dparams[j]);
            }
        }
    }

    free_dvector(va,1,na);
    free_dvector(fvd,1,na);
    free_dvector(fv,1,na);
}

// Second derivatives of the Chi2 respect parameters
// Numerical fourth order computation
// x : x independent coord; a : params; y : model;
// na : number of parameters; d2yda : derivatives
global void model_wd2yda_chi2_old(double a[],
                         double **d2yda, int na)
{
    int i, j, k, l, n;
    real *fv;
//    real *fvd;
    real *va;
    real *vai;
    real *vaj;

    fv = dvector(1,na);
//    fvd = dvector(1,na);
    va = dvector(1,na);
    vai = dvector(1,na);
    vaj = dvector(1,na);

    fv[1] = +1.0;
    fv[2] = -8.0;
    fv[3] = 0.0;
    fv[4] = +8.0;
    fv[5] = -1.0;

//    fvd[1] = -1.0;
//    fvd[2] = +16.0;
//    fvd[3] = -30.0;
//    fvd[4] = +16.0;
//    fvd[5] = -1.0;

    for (i=1; i<=na; i++) {
        for (j=1; j<=na; j++) {
// BEGIN two main loops i, j
            vai[i] = a[i] - 2.0*pd.dparams[i];
            vaj[j] = a[j] - 2.0*pd.dparams[j];
            d2yda[i][j] = 0.0;
            if (i==j) {
                for (l=1; l<=5; l++) {
                    for (k=1; k<=5; k++) {
                        for (n=1; n<=na; n++) {
                            va[n] = a[n];
                        }
                        va[i] = vai[i];
                        va[j] = vaj[j];
                        d2yda[j][j] += fv[k]*gd.modelChi2(gd.modelfun, va);
                        vaj[j] += pd.dparams[j];
                    }
                    vai[i] += pd.dparams[i];
                }
                d2yda[j][j] /= (12.0*pd.dparams[j]*pd.dparams[j]);
            } else {
                for (l=1; l<=5; l++) {
                    for (k=1; k<=5; k++) {
                        for (n=1; n<=na; n++) {
                            va[n] = a[n];
                        }
                        va[i] = vai[i];
                        va[j] = vaj[j];
                        d2yda[i][j] += fv[k]*gd.modelChi2(gd.modelfun, va);
                        vaj[j] += pd.dparams[j];
                    }
                    d2yda[i][j] *= fv[l];
                    vai[i] += pd.dparams[i];
                }
                d2yda[i][j] /= (144.0*pd.dparams[i]*pd.dparams[j]);
            }
// END two main loops i, j
        }
    }

    free_dvector(vaj,1,na);
    free_dvector(vai,1,na);
    free_dvector(va,1,na);
//    free_dvector(fvd,1,na);
    free_dvector(fv,1,na);
}

// Second derivatives of the Chi2 respect parameters
// Numerical fourth order computation
// x : x independent coord; a : params; y : model;
// na : number of parameters; d2yda : derivatives
global void model_wd2yda_chi2(double a[],
                         double **d2yda, int na)
{
    int i, j, k;
    real y;
    real *dfdp;

    dfdp=dvector(1,na);
//    set_model_Interp(gd.modelfun, a);

    for (i=1; i<=na; i++) {
        for (j=1; j<=na; j++) {
            d2yda[i][j] = 0.0;
            for (k=1; k<=nObs(pObs[gd.filecount]); k++) {
                gd.modelfunwd(xObs(pObs[gd.filecount])[k], a, &y, dfdp, na);
                d2yda[i][j] += dfdp[i]*dfdp[j]/rsqr(sigmaObs(pObs[gd.filecount])[k]);
            }
        }
    }

    free_dvector(dfdp,1,na);
}

