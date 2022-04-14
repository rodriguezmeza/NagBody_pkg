/*==============================================================================
	MODULE: models_user.h			[nchi2]
==============================================================================*/

// ==========================================
// Begin: USER Model

#define RCBECHO                     0

// Begin:: Public interfaces:
//global real rcUSER(real x, real params[]);
//global real Chi2_USER(real (*ymodel)(real, real *), real params[]);
//global void Model_USER_end(void);

global real vcTotalBECHO(real x, real params[]);
global void vcTBECHO(double x, double a[], double *y, double dyda[], int na);
global real Chi2_RCBECHO(real (*ymodel)(real, real *), real params[]);
global void Model_RCBECHO_end(void);
// End :: Public interfaces
//

local void Model_RCBECHO(void);
local void test_RCBECHO(void);

local void Model_RCBECHO(void);
local void test_RCBECHO(void);

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
} model_data_rcbecho, *model_data_rcbecho_ptr;

local model_data_rcbecho mdrcbecho;

#define md          mdrcbecho

// params[1] === rhos   :: range: 0 < p1 < 10000
// params[2] === rs     :: range: 0 < p2 < 10

local void Model_RCBECHO(void)
{
    strcpy(gd.model_comment, "RCBECHO Model");

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
        "Model_RCBECHO: no parameter file was given...\nUsing a default set...'\n");
        pd.nparams = 2;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
        pd.params[1] = 4899.44;
        pd.rparams[1][1] = 0.0001;
        pd.rparams[1][2] = 10000.0;
        pd.dparams[1] = 10.0;       // dummy
        pd.params[2] = 3.65868;
        pd.rparams[2][1] = 0.0001;
        pd.rparams[2][2] = 10.0;
        pd.dparams[2] = 0.01;       // dummy
        
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

    gd.modelfun = &vcTotalBECHO;
    gd.modelfunwd = &vcTBECHO;
    gd.modelChi2 = &Chi2_RCBECHO;
    gd.modelend = &Model_RCBECHO_end;
}

local void test_RCBECHO(void)
{
    strcpy(gd.model_comment, "RCBECHO Model");
    
    real Chi2tmp;
    real Chi2red;

    fprintf(stdout,"\n\nTesting %s...\n",gd.model_comment);

//    pd.params[1] = 4899.44;
//    pd.params[2] = 3.65868;
// Chi2     = 1.72175                       // For galaxy F563-V2_rotmod.txt
// Chi2red  = 0.215219

    Chi2tmp = gd.modelChi2(vcTotalBECHO, pd.params);
    Chi2red = Chi2tmp/(nObs(pObs[gd.filecount])-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g\n",pd.params[1],pd.params[2],Chi2tmp,Chi2red);
    plotting_model_table(gd.modelfun,gd.fpfnamemodeltest);

    fprintf(stdout,"done.\n\n");

}

global real Chi2_RCBECHO(real (*ymodel)(real, real *), real params[])
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
                         -model_Interp(xObs(pObs[gd.filecount])[i]))
                        /sigmaObs(pObs[gd.filecount])[i]
                        );
    }
//    (yobs[i] - ymodel(xobs[i], params))/sigma[i]

    return Chi2tmp;
}

#define mparamfinal     "modelparamfinal"
global void Model_RCBECHO_end(void)
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

local real IrhoDM(real x, real params[])
{
    real Irhotmp;
    
    Irhotmp = (1.0/4.0)
    *(
      -2.0*rexp(-rsqr(x)/rsqr(params[2]))*(x/params[2])
      + rsqrt(PI)*erff(x/params[2])
      );
    
    return (Irhotmp);
}

local real massDM(real x, real params[])
{
    real masstmp;
    
    masstmp = params[1]*rpow(params[2],3.0)*(1.0/4.0)
            *(
              -2.0*rexp(-rsqr(x)/rsqr(params[2]))*(x/params[2])
              + rsqrt(PI)*erff(x/params[2])
              );

    return (masstmp);
}

local real vcDM(real x, real params[])
{
    real vcDMtmp;
    
    vcDMtmp = rsqrt(massDM(x,params)/x);
    
    return vcDMtmp;
}

global real vcTotalBECHO(real x, real params[])
{
    real vcTtmp;
    
    vcTtmp = rsqrt( rsqr(vcDM(x,params)) );
    
    return vcTtmp;
}


global void vcTBECHO(double x, double a[], double *y, double dyda[], int na)
{
    int i;
    double fac,ex,arg;
    real vcTotal;

    *y = vcTotalBECHO(x,a);
    vcTotal = vcTotalBECHO(x,a);

    dyda[1] = massDM(x,a) /
    (2.0*a[1]*x*vcTotal);

    dyda[2] = (massDM(x,a) /
            (2.0*a[2]*x*vcTotal))
            *(
                3.0-(x/a[2])*(1.0/IrhoDM(x,a))
                *rexp(-rsqr(x/a[2]))*rsqr(x/a[2])
            );
}
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

