/*==============================================================================
	MODULE: models_user.h			[nchi2]
==============================================================================*/

// ==========================================
// Begin: gSEIR Model (cumulative infected data only)

#define GSEIRINFCUM                 2

// Begin:: Public interfaces:
global real gSEIR_CI(real x, real params[]);
global real Chi2_gSEIR_CI(real (*ymodel)(real, real *), real params[]);
global void Model_gSEIR_CI_end(void);
// End :: Public interfaces
//

local void Model_gSEIR_CI(void);
local void test_gSEIR_CI(void);
local void set_model_Interp_gSEIR_CI(real (*ymodel)(real, real *), real params[]);


// Variables: S (susceptible,1), P (isolated,2), E (exposed,3),
// I (infected,4), Q (quarantine,5)
// Parameters:
// p(1) = alpha :: protection rate
// p(2) = beta :: infection rate
// p(3) = gamma :: inverse of the average latent time
// p(4) = delta :: inverse of the average quarantine time
typedef struct {
// Local parameters and constant
    real Ntotal;
    real I0;
    int method_int;
    real xnow;
    real xstop;
    real eps;
    real dx;
    real dxmin;
    int neqdiffs;
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
    real p3;
    real p3cv;
    real p3min;
    real p3max;
    real dp3;
    real p4;
    real p4cv;
    real p4min;
    real p4max;
    real dp4;
} model_data, *model_data_ptr;

local model_data md_gseir_ci;

#define md      md_gseir_ci

local void integration(double ystart[], int nvar, double x1, double x2, double eps, double h1,
double hmin, int *nok, int *nbad, int maxnsteps,
void (*derivsin)(double, double [], double []));
local void gSEIR_cumulative_infected(double x,double y[],double dydx[]);
void integration_method_string_to_int(string method_str,int *method_int);
local real gSEIR_cumulative_infected_function(real t, double params[]);


local void Model_gSEIR_CI(void)
{
    char integration_method[10];
    strcpy(gd.model_comment, "Cumulative infected gSEIR model");

//    strcpy(integration_method, "rkqs");
    strcpy(integration_method, "bsstep");

    integration_method_string_to_int(integration_method, &md.method_int);

    gd.filecount=0;
    InputObsDataTable(gd.filenames[gd.filecount], gd.filecount);
    pObsT.no = nObs(pObs[gd.filecount]);

    set_model_Interp_init();

// Total population :: example for Mexico
    md.Ntotal = 126200000.0;
// Initial infected population :: example for Mexico at day 55 (csv data)
    md.I0 = 53.0;
    md.neqdiffs = 5;

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
        pd.params[3] = md.p3cv;
        pd.rparams[3][1] = md.p3min;
        pd.rparams[3][2] = md.p3max;
        pd.dparams[3] = md.dp3;
        pd.params[4] = md.p4cv;
        pd.rparams[4][1] = md.p4min;
        pd.rparams[4][2] = md.p4max;
        pd.dparams[4] = md.dp4;
    } else {
#define model_parameter_null    "model_parameters_null"
        fprintf(stdout,
        "Model_gSEIR_CI: no parameter file was given...\nUsing a default set...'\n");
        pd.nparams = 4;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
// Prior: min, max of parameters
        pd.params[1] = 0.03558266;
        pd.rparams[1][1] = 0.00001;
        pd.rparams[1][2] = 0.1;
        pd.dparams[1] = 0.00001;
        pd.params[2] = 0.433269;
        pd.rparams[2][1] = 0.00001;
        pd.rparams[2][2] = 1.0;
        pd.dparams[2] = 0.00001;
        pd.params[3] = 2.41942;
        pd.rparams[3][1] = 0.00001;
        pd.rparams[3][2] = 5.0;
        pd.dparams[3] = 0.00001;
        pd.params[4] = 0.00258743;
        pd.rparams[4][1] = 0.00001;
        pd.rparams[4][2] = 0.01;
        pd.dparams[4] = 0.000001;
        
        md.nparams = pd.nparams;
        md.p1cv = pd.params[1];
        md.p1min = pd.rparams[1][1];
        md.p1max = pd.rparams[1][2];
        md.dp1 = pd.dparams[1];
        md.p2cv = pd.params[2];
        md.p2min = pd.rparams[2][1];
        md.p2max = pd.rparams[2][2];
        md.dp2 = pd.dparams[2];
        md.p3cv = pd.params[3];
        md.p3min = pd.rparams[3][1];
        md.p3max = pd.rparams[3][2];
        md.dp3 = pd.dparams[3];
        md.p4cv = pd.params[4];
        md.p4min = pd.rparams[4][1];
        md.p4max = pd.rparams[4][2];
        md.dp4 = pd.dparams[4];

        PrintModelParameterFile(model_parameter_null);
#undef model_parameter_null
    }

    gd.modelfun = &gSEIR_CI;
    gd.modelfunwd = &model_wdyda;
    gd.modelChi2 = &Chi2_gSEIR_CI;
    gd.modelend = &Model_gSEIR_CI_end;
}

local void test_gSEIR_CI(void)
{
//    strcpy(gd.model_comment, "Cumulative infected gSEIR model");
    real Chi2tmp;
    real Chi2red;

    fprintf(stdout,"\n\nTesting %s...\n",gd.model_comment);

//    pd.params[1] = 0.03558266;
//    pd.params[2] = 0.433269;
//    pd.params[3] = 2.41942;
//    pd.params[4] = 0.00258743;

    Chi2tmp = gd.modelChi2(gd.modelfun, pd.params);
    Chi2red = Chi2tmp/(nObs(pObs[gd.filecount])-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g %g %g\n",
            pd.params[1],pd.params[2],pd.params[3],pd.params[4],Chi2tmp,Chi2red);
    plotting_model_table(gd.modelfun,gd.fpfnamemodeltest);

    fprintf(stdout,"done.\n\n");
}

global real Chi2_gSEIR_CI(real (*ymodel)(real, real *), real params[])
{
    real Chi2tmp;
    int i;

    set_model_Interp(ymodel, params);

    Chi2tmp = 0.0;
    for (i=1; i<=nObs(pObs[gd.filecount]); i++) {
        if (sigmaObs(pObs[gd.filecount])[i]==0.)
            sigmaObs(pObs[gd.filecount])[i] = 1.0;
        Chi2tmp += rsqr(
                        (yObs(pObs[gd.filecount])[i] - model_Interp(xObs(pObs[gd.filecount])[i]))/sigmaObs(pObs[gd.filecount])[i]
                        );
    }

    return Chi2tmp;
}

#define mparamfinal     "modelparamfinal"
global void Model_gSEIR_CI_end(void)
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
    md.p3cv = pd.params[3];
    md.p3min = pd.rparams[3][1];
    md.p3max = pd.rparams[3][2];
    md.dp3 = pd.dparams[3];
    md.p4cv = pd.params[4];
    md.p4min = pd.rparams[4][1];
    md.p4max = pd.rparams[4][2];
    md.dp4 = pd.dparams[4];

    PrintModelParameterFile(mparamfinal);
}
#undef mparamfinal

global real gSEIR_CI(real x, real params[])
{
    real rtmp;
    
    rtmp = gSEIR_cumulative_infected_function(x, params);
    
    return rtmp;
}

// Variables: S (susceptible,1), P (isolated,2), E (exposed,3),
// I (infected,4), Q (quarantine,5)
#define alpha   pd.params[1]
#define beta   pd.params[2]
#define gamma   pd.params[3]
#define delta   pd.params[4]
local void gSEIR_cumulative_infected(double x,double y[],double dydx[])
{
    nrhs++;

//    Ntotal = y[1] + y[2] + y[3] + y[4] + y[5];

// S :: susceptible
    dydx[1] = -alpha*y[1] - beta*y[1]*y[4]/md.Ntotal;
// P :: isolated
    dydx[2]= alpha*y[1];
// E :: exposed
    dydx[3] = -gamma*y[3] + beta*y[1]*y[4]/md.Ntotal;
// I :: infected
    dydx[4] = gamma*y[3] - delta*y[4];
// Q :: quarantine
    dydx[5] = delta*y[4];
}
#undef alpha
#undef beta
#undef gamma
#undef delta

local void set_model_Interp_gSEIR_CI(real (*ymodel)(real, real *), real params[])
{
    int nbad,nok;
    double *ystart;
    int maxnsteps;
//
    int i;
    real xstop=500.0;
//
    ystart=dvector(1,md.neqdiffs);
    xp=dvector(1,200);
    yp=dmatrix(1,md.neqdiffs,1,200);
    md.xnow = 0.0;
    md.xstop = xstop;
    md.eps=1.0e-4;
    maxnsteps=10000;
    md.dxmin = 0.0;
    md.dx = 2.0/5.0;

    ystart[1]=md.Ntotal - md.I0;
    ystart[2]=0.0;
    ystart[3]=0.0;
    ystart[4]=md.I0;
    ystart[5]=0.0;
    
    nrhs=0;
    kmax=100; // maximum number of steps that can be stored
    dxsav=(md.xstop-md.xnow)/20.0;

    integration(ystart,md.neqdiffs,md.xnow,md.xstop,md.eps,md.dx,md.dxmin,
                &nok,&nbad,maxnsteps,gSEIR_cumulative_infected);
//
    gd.nInterp = kount;
    gd.xInterp = dvector(1,gd.nInterp);
    gd.yInterp = dvector(1,gd.nInterp);
    gd.yInterp2 = dvector(1,gd.nInterp);

    for (i=1; i<gd.nInterp; i++) {
        gd.xInterp[i] = xp[i];
        gd.yInterp[i] = yp[5][i];
    }
    spline(gd.xInterp,gd.yInterp,gd.nInterp,1.0e30,1.0e30,gd.yInterp2);
//
    free_dmatrix(yp,1,md.neqdiffs,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,md.neqdiffs);
}

local real gSEIR_cumulative_infected_function(real t, double params[])
{
    int nbad,nok;
    double *ystart;
    real ptmp;
    int maxnsteps;
//
    ystart=dvector(1,md.neqdiffs);
    xp=dvector(1,200);
    yp=dmatrix(1,md.neqdiffs,1,200);
    md.xnow = 0.0;
    md.xstop = t;
    md.eps=1.0e-4;
    maxnsteps=10000;
    md.dxmin = 0.0;
    md.dx = 2.0/5.0;

    ystart[1]=md.Ntotal - md.I0;
    ystart[2]=0.0;
    ystart[3]=0.0;
    ystart[4]=md.I0;
    ystart[5]=0.0;
    
    nrhs=0;
    kmax=100; // maximum number of steps that can be stored
    dxsav=(md.xstop-md.xnow)/20.0;

    integration(ystart,md.neqdiffs,md.xnow,md.xstop,md.eps,md.dx,md.dxmin,
                &nok,&nbad,maxnsteps,gSEIR_cumulative_infected);

    ptmp = yp[5][kount];
//
    free_dmatrix(yp,1,md.neqdiffs,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,md.neqdiffs);
    
    return ptmp;
}


#define BSSTEP          0
#define NULLMETHOD      1
#define RKQS            2
local void integration(double ystart[], int nvar, double x1, double x2, double eps, double h1,
                       double hmin, int *nok, int *nbad, int maxnsteps,
                       void (*derivsin)(double, double [], double []))
{
    switch (md.method_int) {
        case BSSTEP:
//            odeint(ystart,nvar,xnow,xstop,eps,dx,dxmin,nok,nbad,maxnsteps,derivsin,bsstep);
            odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,maxnsteps,derivsin,bsstep);
            break;
//
        case RKQS:
            odeint(ystart,nvar,md.xnow,md.xstop,md.eps,md.dx,md.dxmin,nok,nbad,maxnsteps,derivsin,rkqs);
            break;
//
        case NULLMETHOD:
            odeint(ystart,nvar,md.xnow,md.xstop,md.eps,md.dx,md.dxmin,nok,nbad,maxnsteps,derivsin,bsstep);
            break;
//
        default:
            odeint(ystart,nvar,md.xnow,md.xstop,md.eps,md.dx,md.dxmin,nok,nbad,maxnsteps,derivsin,bsstep);
            break;
    }

//    exit(1);
}

void integration_method_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"bsstep") == 0) {
        *method_int = BSSTEP;
//        strcpy(gd.integration_method_comment, "bsstep integration method");
    }
//
    if (strcmp(method_str,"rkqs") == 0) {
        *method_int = RKQS;
//        strcpy(gd.integration_method_comment, "rkqs integration method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
//        strcpy(gd.integration_method_comment,
//               "null integration method ... running deafult (bsstep)");
//        fprintf(stdout,"\n\tintegration: default integration method (bsstep)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = BSSTEP;
//        strcpy(gd.integration_method_comment,
//               "Unknown integration method ... running deafult (bsstep)");
//        fprintf(stdout,"\n\tintegration: Unknown method... %s ",integration_method);
//        fprintf(stdout,
//                "\n\trunning default integration method (bsstep)...\n");
    }
}
#undef BSSTEP
#undef RKQS
#undef NULLMETHOD

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
        RPName(md.p3cv,"p3cv");         \
        RPName(md.p3min,"p3min");       \
        RPName(md.p3max,"p3max");       \
        RPName(md.dp3,"dp3");           \
        RPName(md.p4cv,"p4cv");         \
        RPName(md.p4min,"p4min");       \
        RPName(md.p4max,"p4max");       \
        RPName(md.dp4,"dp4");           \
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
        fprintf(fdout,FMTR,"p3cv",md.p3cv);         \
        fprintf(fdout,FMTR,"p3min",md.p3min);       \
        fprintf(fdout,FMTR,"p3max",md.p3max);       \
        fprintf(fdout,FMTR,"dp3",md.dp3);           \
        fprintf(fdout,FMTR,"p4cv",md.p4cv);         \
        fprintf(fdout,FMTR,"p4min",md.p4min);       \
        fprintf(fdout,FMTR,"p4max",md.p4max);       \
        fprintf(fdout,FMTR,"dp4",md.dp4);           \
}
//

//
#include "../models_io.h"

#undef READPARAMS
#undef WRITEPARAMS
*/

#undef md


// End: gSEIR Model (cumulative infected data only)
// ==========================================

