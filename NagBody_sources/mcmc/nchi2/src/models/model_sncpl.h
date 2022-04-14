/*==============================================================================
	MODULE: models_user.h			[nchi2]
==============================================================================*/

// ==========================================
// Begin: SN CPL Model

#define SNCPL                       3

// Begin:: Public interfaces:
//global real rcUSER(real x, real params[]);
//global real Chi2_USER(real (*ymodel)(real, real *), real params[]);
//global void Model_USER_end(void);

global real snCPL(real z, real params[]);
global real Chi2_SNCPL(real (*ymodel)(real, real *), real params[]);
global void Model_SNCPL_end(void);
// End :: Public interfaces
//

local void Model_SNCPL(void);
local void test_SNCPL(void);

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
} model_data_sncpl, *model_data_sncpl_ptr;

local model_data_sncpl mdsncpl;

local real rrSN(real z);

#define md          mdsncpl
local void Model_SNCPL(void)
{
    strcpy(gd.model_comment, "SNCPL Model");

    gd.filecount=0;
    InputObsDataTable(gd.filenames[gd.filecount], gd.filecount);
    pObsT.no = nObs(pObs[gd.filecount]);

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
        "Model_SNCPL: no parameter file was given...\nUsing a default set...'\n");
        pd.nparams = 2;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
// Priors:: In CosmoMC: central value, min, max, step and ...
        pd.params[1] = -1.0;
        pd.rparams[1][1] = -2.0;
        pd.rparams[1][2] = 0.1;
        pd.dparams[1] = 0.01;
        pd.params[2] = 0.0;
        pd.rparams[2][1] = -5.0;
        pd.rparams[2][2] = 5.0;
        pd.dparams[2] = 0.01;
        
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

    gd.modelfun = &snCPL;
    gd.modelfunwd = &model_wdyda;
    gd.modelChi2 = &Chi2_SNCPL;
    gd.modelend = &Model_SNCPL_end;
}

local void test_SNCPL(void)
{
//    strcpy(gd.model_comment, "SNCPL Model");
    real Chi2tmp;
    real Chi2red;

    fprintf(stdout,"\n\nTesting %s...\n",gd.model_comment);

//    pd.params[1] = -1.0;
//    pd.params[2] = 0.0;

    Chi2tmp = gd.modelChi2(gd.modelfun, pd.params);
    Chi2red = Chi2tmp/(nObs(pObs[gd.filecount])-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g\n",pd.params[1],pd.params[2],Chi2tmp,Chi2red);

    pd.params[1] = -0.99822;
    pd.params[2] = -0.331855;
    Chi2tmp = Chi2_SNCPL(gd.modelfun, pd.params);
    Chi2red = Chi2tmp/(nObs(pObs[gd.filecount])-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g\n",pd.params[1],pd.params[2],Chi2tmp,Chi2red);

    plotting_model_table(gd.modelfun,gd.fpfnamemodeltest);

    fprintf(stdout,"done.\n\n");
}

global real Chi2_SNCPL(real (*ymodel)(real, real *), real params[])
{
    real Chi2tmp;
    int i;
    real *diffs;
    real AG, BG, CG;

    diffs = dvector(1,nObs(pObs[gd.filecount]));

    set_model_Interp(ymodel, params);

    for (i=1; i<=nObs(pObs[gd.filecount]); i++) {
        diffs[i] = yObs(pObs[gd.filecount])[i]
                    - model_Interp(xObs(pObs[gd.filecount])[i]);
    }
//    diffs[i] = (yobs[i] - ymodel(xobs[i], params));

// Marginalized version
    AG = 0.;
    for (i=1; i<=nObs(pObs[gd.filecount]); i++)
        AG += rsqr(diffs[i])/sigmaObs(pObs[gd.filecount])[i];
//    AG += rsqr(yobs[i] - model_Interp(xobs[i]))/sigma[i];

    BG = 0.;
    for (i=1; i<=nObs(pObs[gd.filecount]); i++)
        BG += diffs[i]/sigmaObs(pObs[gd.filecount])[i];

    CG = 0.;
    for (i=1; i<=nObs(pObs[gd.filecount]); i++)
        CG += 1.0/sigmaObs(pObs[gd.filecount])[i];

    Chi2tmp = AG - rsqr(BG)/CG;

// Normal version:
/*
    Chi2tmp = 0.0;
    for (i=1; i<=nobs; i++) {
        if (sigma[i]==0.)
            sigma[i] = 1.0;
        Chi2tmp += rsqr(
                        (yobs[i] - model_Interp(xobs[i]))/sigma[i]
                        );
    }
//    (yobs[i] - ymodel(xobs[i], params))/sigma[i]
*/

    free_dvector(diffs,1,nObs(pObs[gd.filecount]));

    return Chi2tmp;
}

#define mparamfinal     "modelparamfinal"
global void Model_SNCPL_end(void)
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


global real snCPL(real z, real params[])
{
#define w0    params[1]
#define w1      params[2]
    real res;
    real dL;
    real kk;

//    kk = rlog(10.0)*(md.c/md.h)/1000.0;
//    kk = 1.0/rlog(10.0)*(md.c/md.h)/1000.0;
//    kk =(md.h/md.c)*1000.0;
//    kk =(md.h/md.c)*1000.0/rlog(10.0);
//    kk =md.h;
//    kk = 1.0;

    kk = (md.c/md.h)/1000.0;
    dL = kk*rrSN(z)*(1.+z);
    res = 5.0*rlog10(dL) + 25.0;

// These lines give CL results in the test routine:
//    dL = rrSN(z)*(1.+z);
//    res = 5.0*rlog10(dL);


    return res;
#undef w0
#undef w1
}

#define w0    pd.params[1]
#define w1    pd.params[2]
// Here it is defined CPL model:
local real Efun_int(real z)
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
local real rrSN(real z)
{
    real result;
    real zmin;

    zmin = 0.0;

    result=qromo(Efun_int,zmin,z,midpnt,cmd.epsq,KK);

    return result;
}
#undef KK
/*
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
*/
#undef md


// End: SN CPL Model
// ==========================================

