/*==============================================================================
	MODULE: models_rcpiso.h			[nchi2]
==============================================================================*/

// ==========================================
// Begin: RCPISO Model

#define RCPISO                      1

// Begin:: Public interfaces:
//global real rcUSER(real x, real params[]);
//global real Chi2_USER(real (*ymodel)(real, real *), real params[]);
//global void Model_USER_end(void);

global real rcPISO(real x, real params[]);
global void rcPISO_wdyda(double x, double a[], double *y, double dyda[], int na);
global real Chi2_RCPISO(real (*ymodel)(real, real *), real params[]);
global void Model_RCPISO_end(void);

// End :: Public interfaces
//

local void Model_RCPISO(void);
local void test_RCPISO(void);

local void Model_RCPISO(void);
local void test_RCPISO(void);

typedef struct {
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
} model_data_rcpiso, *model_data_rcpiso_ptr;

local model_data_rcpiso mdrcpiso;

#define md          mdrcpiso

// params[1] === rhos   :: range: 0 < p1 < 10000
// params[2] === rs     :: range: 0 < p2 < 10
local void Model_RCPISO(void)
{
    
    strcpy(gd.model_comment, "RCPISO Model");
    
    gd.filecount=0;
    InputObsDataTable(gd.filenames[gd.filecount], gd.filecount);
    pObsT.no = nObs(pObs[gd.filecount]);

    set_model_Interp_init();

    if (!strnull(cmd.model_paramfile)) {
        pd.nparams = md.nparams;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
        pd.params[1] = md.p1cv; // 4899.44;
        pd.rparams[1][1] = md.p1min; // 0.0001;
        pd.rparams[1][2] = md.p1max; // 10000.0;
        pd.dparams[1] = md.dp1;
        pd.params[2] = md.p2cv; // 3.65868;
        pd.rparams[2][1] = md.p2min; // 0.0001;
        pd.rparams[2][2] = md.p2max; // 10.0;
        pd.dparams[2] = md.dp2;
    } else {
#define model_parameter_null    "model_parameters_null"
        fprintf(stdout,
        "Model_RCPISO: no parameter file was given...\nUsing a default set...'\n");
        pd.nparams = 2;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
        pd.params[1] = 4899.44;
        pd.rparams[1][1] = 0.0001;
        pd.rparams[1][2] = 10000.0;
        pd.dparams[1] = 1.0;       // dummy
        pd.params[2] = 3.65868;
        pd.rparams[2][1] = 0.0001;
        pd.rparams[2][2] = 10.0;
        pd.dparams[2] = 1.0;       // dummy

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

    gd.modelfun = &rcPISO;
    gd.modelfunwd = &rcPISO_wdyda;
    gd.modelChi2 = &Chi2_RCPISO;
    gd.modelend = &Model_RCPISO_end;
}

local void test_RCPISO(void)
{
    
    strcpy(gd.model_comment, "RCPISO Model");

    real Chi2tmp;
    real Chi2red;

    fprintf(stdout,"\n\nTesting %s model...\n",gd.model_comment);

//    pd.params[1] = 4899.44;
//    pd.params[2] = 3.65868;

    Chi2tmp = gd.modelChi2(gd.modelfun, pd.params);
    Chi2red = Chi2tmp/(nObs(pObs[gd.filecount])-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g\n",pd.params[1],pd.params[2],Chi2tmp,Chi2red);
    plotting_model_table(gd.modelfun,gd.fpfnamemodeltest);

    fprintf(stdout,"done.\n\n");
}

global real Chi2_RCPISO(real (*ymodel)(real, real *), real params[])
{
    real Chi2tmp;
    int i;

    set_model_Interp(ymodel, params);

    Chi2tmp = 0.0;
    for (i=1; i<=nObs(pObs[gd.filecount]); i++) {
        if (sigmaObs(pObs[gd.filecount])[i]==0.)
            sigmaObs(pObs[gd.filecount])[i] = 1.0;
        Chi2tmp += rsqr(
                        (
                         yObs(pObs[gd.filecount])[i]
                         - model_Interp(xObs(pObs[gd.filecount])[i]))
                        /sigmaObs(pObs[gd.filecount])[i]
                        );
    }
//    (yobs[i] - ymodel(xobs[i], params))/sigma[i]

    return Chi2tmp;
}

#define mparamfinal     "modelparamfinal"
global void Model_RCPISO_end(void)
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

#define rhos    params[1]
#define rs      params[2]
global real rcPISO(real x, real params[])
{
    real vcTtmp;
    real mDM;
// PISO model density profile:
// rhos/(1 + (x/rs)^2)

    mDM = rpow(rs,3.0)*rhos* (x/rs - ratan(x/rs));
    
    vcTtmp = rsqrt( mDM/x );
    
    return vcTtmp;
}
#undef rhos
#undef rs

#define rhos    a[1]
#define rs      a[2]
global void rcPISO_wdyda(double x, double a[], double *y, double dyda[], int na)
{
    int i;
    double fac,ex,arg;
    real vcTotal, mDM;

    *y = rcPISO(x,a);
    vcTotal = rcPISO(x,a);
    mDM = rpow(rs,3.0)*rhos* (x/rs - ratan(x/rs));

    dyda[1] = mDM /
    (2.0*rhos*x*vcTotal);

    dyda[2] = (1.0 /
            (2.0*x*vcTotal))
            *(
              (3.0/rs)*mDM
              + (rpow(rs,3.0)*rhos)*( -(x/rsqr(rs)) + x/(rsqr(rs)*(1. + rsqr(x)/rsqr(rs))) )
            );
}
#undef rhos
#undef rs
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



// End: USER Model
// ==========================================

