/*==============================================================================
	MODULE: models_user.h			[nchi2]
==============================================================================*/

// ==========================================
// Begin: BAO CPL Model

#define BAOCPL                      4

// Begin:: Public interfaces:
global real baoCPL(real z, real params[]);
global real Chi2_BAOCPL(real (*ymodel)(real, real *), real params[]);
global void Model_BAOCPL_end(void);
// End :: Public interfaces
//

local void Model_BAOCPL(void);
local void test_BAOCPL(void);

// params[1] === w0   :: range: -10 < p1 < 10
// params[2] === w1   :: range: -10 < p2 < 10
typedef struct {
    // Local parameters and constant
    real c;
    real h;
    real omi;
    real obi;
    real ome;
    real obe;
    real om;
    real ob;
// Combined obs data
//    int nObsTotal;
// Fit:: Parameter's priors
    int nparams;
    real p1;
    real p1cv;
    real p1min;
    real p1max;
    real dp1;
    real p2;
    real p2cv;
    real p2min;
    real p2max;
    real dp2;
} model_data_baocpl, *model_data_baocpl_ptr;

local model_data_baocpl mdbaocpl;

local real rrBAO(real z);

local real Chi2_BAOCPL_dataset1(real (*ymodel)(real, real *), real params[],
real xobs[], real yobs[], real sigma[], int nobs);
local real Chi2_BAOCPL_WiggleZ(real (*ymodel)(real, real *), real params[],
real xobs[], real yobs[], real sigma[], int nobs);

#define md          mdbaocpl
local void Model_BAOCPL(void)
{
    strcpy(gd.model_comment, "BAOCPL Model");

    gd.filecount=0;
    InputObsDataTable(gd.filenames[gd.filecount], gd.filecount);
    if (gd.nfiles == 2) {
        InputObsDataTable(gd.filenames[gd.filecount+1], gd.filecount+1);
//        md.nObsTotal = nObs(pObs[gd.filecount]) + nObs(pObs[gd.filecount+1]);
        pObsT.no = nObs(pObs[gd.filecount]) + nObs(pObs[gd.filecount+1]);
        fprintf(stdout,"\nTotal observed data: %d\n",pObsT.no);
    } else
        error("\n\nMOdel_BAOCPL :: two files are needed\n");

    set_model_Interp_init();

    md.c = 2.99792458e8;
    md.h = 70.0;
    md.omi = 0.295;
    md.obi = 0.044898;
    md.ome = 0.034;
    md.obe = 0.00054;
    md.om = md.omi;
    md.ob = md.obi;

    if (!strnull(cmd.model_paramfile)) {
        pd.nparams = md.nparams;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
        pd.params[1] = md.p1cv;
        pd.rparams[1][1] = md.p1min;
        pd.rparams[1][2] = md.p1max;
        pd.dparams[1] = md.dp1;
        pd.params[2] = md.p2cv;
        pd.rparams[2][1] = md.p2min;
        pd.rparams[2][2] = md.p2max;
        pd.dparams[2] = md.dp2;
    } else {
#define model_parameter_null    "model_parameters_null"
        fprintf(stdout,
        "\nModel_BAOCPL: no parameter file was given...\nUsing a default set...'\n");
        pd.nparams = 2;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
// Priors:: In CosmoMC: central value, min, max, step and ...
        pd.params[1] = -1.0;
        pd.rparams[1][1] = -1.5;
        pd.rparams[1][2] = 0.;
        pd.dparams[1] = 0.05;
        pd.params[2] = 0.0;
        pd.rparams[2][1] = -3.0;
        pd.rparams[2][2] = 3.0;
        pd.dparams[2] = 0.05;

        md.nparams = pd.nparams;
        md.p1cv = pd.params[1];
        md.p1min = pd.rparams[1][1];
        md.p1max = pd.rparams[1][2];
        md.dp1 = pd.dparams[1];
        md.p2cv = pd.params[2];
        md.p2min = pd.rparams[2][1];
        md.p2max = pd.rparams[2][2];
        md.dp2 = pd.dparams[2];

        PrintModelParameterFile(model_parameter_null);
#undef model_parameter_null
    }

    gd.modelfun = &baoCPL;
    gd.modelfunwd = &model_wdyda;
    gd.modelChi2 = &Chi2_BAOCPL;
    gd.modelend = &Model_BAOCPL_end;
}

local void test_BAOCPL(void)
{
    real Chi2tmp;
    real Chi2red;

    fprintf(stdout,"\n\nTesting %s...\n",gd.model_comment);

    Chi2tmp = gd.modelChi2(gd.modelfun, pd.params);
    Chi2red = Chi2tmp/(pObsT.no-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g\n",pd.params[1],pd.params[2],Chi2tmp,Chi2red);

// Chi2red = 0.6151224925532289
    pd.params[1] = -0.7752023692413474;
    pd.params[2] = -2.5912950877901255;
    Chi2tmp = Chi2_BAOCPL(baoCPL, pd.params);
    Chi2red = Chi2tmp/(pObsT.no-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g\n",pd.params[1],pd.params[2],Chi2tmp,Chi2red);

    plotting_model_table(gd.modelfun,gd.fpfnamemodeltest);

    fprintf(stdout,"done.\n\n");
}

global real Chi2_BAOCPL(real (*ymodel)(real, real *), real params[])
{
    real Chi2tmp;

    Chi2tmp = Chi2_BAOCPL_dataset1(ymodel, params,
            xObs(pObs[gd.filecount]), yObs(pObs[gd.filecount]),
            sigmaObs(pObs[gd.filecount]), nObs(pObs[gd.filecount]));
    Chi2tmp += Chi2_BAOCPL_WiggleZ(ymodel, params,
            xObs(pObs[gd.filecount+1]), yObs(pObs[gd.filecount+1]),
            sigmaObs(pObs[gd.filecount+1]), nObs(pObs[gd.filecount+1]));

    return Chi2tmp;
}

local real Chi2_BAOCPL_dataset1(real (*ymodel)(real, real *),
                 real params[],
                 real xobs[], real yobs[], real sigma[], int nobs)
{
    real Chi2tmp;
    int i, j;
    real *diffs;
    real **percinvcov;
    real *matmul;
    
    diffs = dvector(1,nobs);
    percinvcov = dmatrix(1,nobs,1,nobs);
    matmul = dvector(1,nobs);

    percinvcov[1][1] = 4444.0;
    percinvcov[1][2] = 0.0;
    percinvcov[1][3] = 0.0;
    percinvcov[2][1] = 0.0;
    percinvcov[2][2] = 215156.0;
    percinvcov[2][3] = 0.0;
    percinvcov[3][1] = 0.0;
    percinvcov[3][2] = 0.0;
    percinvcov[3][3] = 721487.0;

    set_model_Interp(ymodel, params);

    for (i=1; i<=nobs; i++) {
//        diffs[i] = yobs[i] - model_Interp(xobs[i]);
        diffs[i] = (yobs[i] - ymodel(xobs[i], params));
    }

    for (i=1; i<=nobs; i++) {
        matmul[i] = 0.0;
        for (j=1; j<=nobs; j++) {
            matmul[i] = percinvcov[i][j] * diffs[j];
        }
    }
    Chi2tmp = 0.0;
    for (i=1; i<=nobs; i++)
        Chi2tmp += diffs[i] * matmul[i];

    free_dvector(matmul,1,nobs);
    free_dmatrix(percinvcov, 1, nobs, 1, nobs);
    free_dvector(diffs,1,nobs);

    return Chi2tmp;
}

local real Chi2_BAOCPL_WiggleZ(real (*ymodel)(real, real *),
                 real params[],
                 real xobs[], real yobs[], real sigma[], int nobs)
{
    real Chi2tmp;
    int i, j;
    real *diffs;
    real **percinvcov;
    real *matmul;
    
    diffs = dvector(1,nobs);
    percinvcov = dmatrix(1,nobs,1,nobs);
    matmul = dvector(1,nobs);

    percinvcov[1][1] = 1040.3;
    percinvcov[1][2] = -807.5;
    percinvcov[1][3] = 336.8;
    percinvcov[2][1] = -807.5;
    percinvcov[2][2] = 3720.3;
    percinvcov[2][3] = -1551.9;
    percinvcov[3][1] = 336.8;
    percinvcov[3][2] = -1551.9;
    percinvcov[3][3] = 2914.9;

    set_model_Interp(ymodel, params);

    for (i=1; i<=nobs; i++) {
//        diffs[i] = yobs[i] - model_Interp(xobs[i]);
        diffs[i] = (yobs[i] - ymodel(xobs[i], params));
    }

    for (i=1; i<=nobs; i++) {
        matmul[i] = 0.0;
        for (j=1; j<=nobs; j++) {
            matmul[i] = percinvcov[i][j] * diffs[j];
        }
    }
    Chi2tmp = 0.0;
    for (i=1; i<=nobs; i++)
        Chi2tmp += diffs[i] * matmul[i];

    free_dvector(matmul,1,nobs);
    free_dmatrix(percinvcov, 1, nobs, 1, nobs);
    free_dvector(diffs,1,nobs);

    return Chi2tmp;
}

#define mparamfinal     "modelparamfinal"
global void Model_BAOCPL_end(void)
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

    PrintModelParameterFile(mparamfinal);
}
#undef mparamfinal

global real baoCPL(real z, real params[])
{
#define w0    params[1]
#define w1      params[2]
    real res;
//    real dL;
//    real kk;
    real Dv;
    real expde;
    real Efun;
    real rszdrag;
    real baotest;
    
    expde = rpow(1.0+z,3.0*(1.+w0+w1))*rexp(-3.0*w1*z/(1.0+z));
    Efun = rsqrt( md.om*rpow(1.+z,3.0) + (1.-md.om)*expde);

    Dv = rpow(
              (md.c/(1000.*md.h)) * rsqr(rrBAO(z)) * z/Efun,
              1.0/3.0);
    
    rszdrag = 153.5* rpow(md.ob* rsqr(md.h/100.0)/0.02273,(-0.134))
    *rpow(md.om* rsqr(md.h/100.0)/0.1326,(-0.255));
    
    baotest = rszdrag/Dv;

//    kk = (md.c/md.h)/1000.0;
//    dL = kk*rrSN(z)*(1.+z);
//    res = 5.0*rlog10(dL) + 25.0;
    res = baotest;

    return res;
#undef w0
#undef w1
}

#define w0    pd.params[1]
#define w1    pd.params[2]
// Here it is defined CPL model:
local real Efun_int_bao(real z)
{
    real expde;
    real Efun;

    expde = rpow(1.0+z,3.0*(1.+w0+w1))*rexp(-3.0*w1*z/(1.0+z));
    Efun = rsqrt( md.om*rpow(1.+z,3.0) + (1.-md.om)*expde);

    return 1.0/Efun;
}
#undef w0
#undef w1

#define KK  5
// Integrate of 1/Efun wrt z, from z1, z2
local real rrBAO(real z)
{
    real result;
    real zmin;

    zmin = 0.0;

    result = (md.c/(1000.*md.h))
        *qromo(Efun_int_bao,zmin,z,midpnt,cmd.epsq,KK);

    return result;
}
#undef KK

// Read model parameters:
#define READPARAMS() \
    { \
        IPName(md.nparams,"nparams");   \
        RPName(md.p1cv,"p1cv");         \
        RPName(md.p1min,"p1min");       \
        RPName(md.p1max,"p1max");       \
        RPName(md.dp1,"dp1");           \
        RPName(md.p2cv,"p2cv");         \
        RPName(md.p2min,"p2min");       \
        RPName(md.p2max,"p2max");       \
        RPName(md.dp2,"dp2");           \
    }
//

// Write model parameters:
#define WRITEPARAMS() \
{ \
        fprintf(fdout,FMTI,"nparams",md.nparams);   \
        fprintf(fdout,FMTR,"p1cv",md.p1cv);         \
        fprintf(fdout,FMTR,"p1min",md.p1min);       \
        fprintf(fdout,FMTR,"p1max",md.p1max);       \
        fprintf(fdout,FMTR,"dp1",md.dp1);           \
        fprintf(fdout,FMTR,"p2cv",md.p2cv);         \
        fprintf(fdout,FMTR,"p2min",md.p2min);       \
        fprintf(fdout,FMTR,"p2max",md.p2max);       \
        fprintf(fdout,FMTR,"dp2",md.dp2);           \
}
//

//
#include "../models_io.h"

#undef READPARAMS
#undef WRITEPARAMS

#undef md

// End: BAO CPL Model
// ==========================================

