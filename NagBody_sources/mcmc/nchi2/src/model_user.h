/*==============================================================================
	MODULE: models_user.h			[nchi2]
==============================================================================*/

// ==========================================
// Begin: USER Model

// AS IS HERE IS EQUIVALENT TO RCPISO MODEL... results must be (almost...) the same.

#define USERMODEL 100

// Begin:: Public interfaces:
global real rcUSER(real x, real params[]);
global real Chi2_USER(real (*ymodel)(real, real *), real params[]);
global void Model_USER_end(void);
// End :: Public interfaces
//
// User model :: reading/printing parameters
// Public interfaces but exclusive of user model:
//global void ReadModelParameterFile_user(char *);
//global void PrintModelParameterFile_user(char *);
//global void CheckParameters_model_user(void);
//

local void test_USER(void);

//=================================================================
// Begin: Definitions of reading tables and making interpolation defs
//
#define SIXPI2  59.21762640653615
//int nallfunT;
//char fpfnameallfun[100];
//local real sigma2;
local real M1sigma2v;
local real M1sigma2Psi;
local real M1f0;
//
// -----------------------------------------
// B :: PSL table structure
local real *kpslT;
local real *pslT;
local real *pslT2;
local int npslfunT;
// E :: PSL table structure
// -----------------------------------------
//
// -----------------------------------------
// B :: FK table structure
local real *kfkT;
local real *fkT;
local real *fkT2;
local int nfkfunT;
// E :: FK table structure
// -----------------------------------------
//
// -----------------------------------------
// B :: monopole table structure 2
local real *kM1MnbT;
local real *M1MnbT;
local real *M1MnbT2;
local int nM1MnbT;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 3
local real *kM1Mb1T;
local real *M1Mb1T;
local real *M1Mb1T2;
local int nM1Mb1T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 4
local real *kM1Mb12T;
local real *M1Mb12T;
local real *M1Mb12T2;
local int nM1Mb12T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 5
local real *kM1Mb2T;
local real *M1Mb2T;
local real *M1Mb2T2;
local int nM1Mb2T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 6
local real *kM1Mb22T;
local real *M1Mb22T;
local real *M1Mb22T2;
local int nM1Mb22T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 7
local real *kM1Mb2b1T;
local real *M1Mb2b1T;
local real *M1Mb2b1T2;
local int nM1Mb2b1T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 8
local real *kM1Mb2bs2T;
local real *M1Mb2bs2T;
local real *M1Mb2bs2T2;
local int nM1Mb2bs2T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 9
local real *kM1Mbs2T;
local real *M1Mbs2T;
local real *M1Mbs2T2;
local int nM1Mbs2T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 10
local real *kM1Mbs22T;
local real *M1Mbs22T;
local real *M1Mbs22T2;
local int nM1Mbs22T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 11
local real *kM1Mbs2b1T;
local real *M1Mbs2b1T;
local real *M1Mbs2b1T2;
local int nM1Mbs2b1T;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 12
local real *kM1Mb1b3nlT;
local real *M1Mb1b3nlT;
local real *M1Mb1b3nlT2;
local int nM1Mb1b3nlT;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 13
local real *kM1Mb3nlT;
local real *M1Mb3nlT;
local real *M1Mb3nlT2;
local int nM1Mb3nlT;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: monopole table structure 14
local real *kM1MeftT;
local real *M1MeftT;
local real *M1MeftT2;
local int nM1MeftT;
// E :: monopole table structure
// -----------------------------------------
// -----------------------------------------
// B :: quadrupole table structure 2
local real *kM1QnbT;
local real *M1QnbT;
local real *M1QnbT2;
local int nM1QnbT;
// E :: quadrupole table structure
// -----------------------------------------
// -----------------------------------------
// B :: quadrupole table structure 3
local real *kM1Qb1T;
local real *M1Qb1T;
local real *M1Qb1T2;
local int nM1Qb1T;
// E :: quadrupole table structure
// -----------------------------------------
// -----------------------------------------
// B :: quadrupole table structure 4
local real *kM1Qb12T;
local real *M1Qb12T;
local real *M1Qb12T2;
local int nM1Qb12T;
// E :: quadrupole table structure
// -----------------------------------------
// -----------------------------------------
// B :: quadrupole table structure 5
local real *kM1Qb2T;
local real *M1Qb2T;
local real *M1Qb2T2;
local int nM1Qb2T;
// E :: quadrupole table structure
// -----------------------------------------
// -----------------------------------------
// B :: quadrupole table structure 6
local real *kM1Qbs2T;
local real *M1Qbs2T;
local real *M1Qbs2T2;
local int nM1Qbs2T;
// E :: quadrupole table structure
// -----------------------------------------
// -----------------------------------------
// B :: quadrupole table structure 7
local real *kM1Qb3nlT;
local real *M1Qb3nlT;
local real *M1Qb3nlT2;
local int nM1Qb3nlT;
// E :: quadrupole table structure
// -----------------------------------------
// -----------------------------------------
// B :: quadrupole table structure 8
local real *kM1QeftT;
local real *M1QeftT;
local real *M1QeftT2;
local int nM1QeftT;
// E :: quadrupole table structure
// -----------------------------------------
//
// Hexapoles:
// -----------------------------------------
// B :: hexapoles table structure 2
local real *kM1HnbT;
local real *M1HnbT;
local real *M1HnbT2;
local int nM1HnbT;
// E :: hexapoles table structure
// -----------------------------------------
// -----------------------------------------
// B :: hexapoles table structure 3
local real *kM1Hb1T;
local real *M1Hb1T;
local real *M1Hb1T2;
local int nM1Hb1T;
// E :: hexapoles table structure
// -----------------------------------------
// -----------------------------------------
// B :: hexapoles table structure 4
local real *kM1Hb12T;
local real *M1Hb12T;
local real *M1Hb12T2;
local int nM1Hb12T;
// E :: hexapoles table structure
// -----------------------------------------
// -----------------------------------------
// B :: hexapoles table structure 5
local real *kM1HeftT;
local real *M1HeftT;
local real *M1HeftT2;
local int nM1HeftT;
// E :: hexapoles table structure
// -----------------------------------------

local void read_allfunctions(void);
local real sigma2v_function_int(real y);
local real sigma2v_function(void);
local real sigma2Psi_function_int(real y);
local real sigma2Psi_function(void);
local real M1pk(real kv);
local real M1fk(real kv);
local real M1Mnb(real kv);
local real M1Mb1(real kv);
local real M1Mb12(real kv);
local real M1Mb2(real kv);
local real M1Mb22(real kv);
local real M1Mb2b1(real kv);
local real M1Mb2bs2(real kv);
local real M1Mbs2(real kv);
local real M1Mbs22(real kv);
local real M1Mbs2b1(real kv);
local real M1Mb1b3nl(real kv);
local real M1Mb3nl(real kv);
local real M1Meft(real kv);
// Quadrupoles:
local real M1Qnb(real kv);
local real M1Qb1(real kv);
local real M1Qb12(real kv);
local real M1Qb2(real kv);
local real M1Qbs2(real kv);
local real M1Qb3nl(real kv);
local real M1Qeft(real kv);
// Hexapoles:
local real M1Hnb(real kv);
local real M1Hb1(real kv);
local real M1Hb12(real kv);
local real M1Heft(real kv);

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);
local  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS);
local  real fkInterpolation_nr(real k, double kfk[], double pfk[], int nfk);

// (*Construct Multipoles*)
local real factor(real k);
local real M1HexadecapoleME(real k);
local real M1QuadrupoleME(real k);
local real M1MonopoleME(real k);
local real DNLOell4(real k, real b1);
local real DNLOell2(real k, real b1);
local real DNLOell0(real k, real b1);
local real PM1MEell4(real k, real b1, real b2, real bs2, real b3nl, real cell4);
local real PM1MEell2(real k, real b1, real b2, real bs2, real b3nl, real cell2);
local real PM1MEell0(real k, real b1, real b2, real bs2, real b3nl, real cell0);

// End: Definitions of reading tables and making interpolation defs
//=================================================================

//(*Parameters, BIAS and EFT*)
// b1 = 1.45;    (*Prior  1, 2.5*)
//
// b2 = -0.2;    (*Prior  -1, 1*)
//
// bs2 = -4/7 (b1 - 1);   (*PriorFixed. This is the coevolution value. But it can \
// change a little. Maybe add a small shift*)
//
// b3nl = 32/315 (b1 - 1);    (*PriorFixed. This is the coevolution value. But it can \
// change a little. Maybe add a small shift*)
//
// cell0 = -12;   (*Prior  -40, 0. It is strictly negative*)
//
// cell2 = -32;  (*Prior  -60, 0. It is strictly negative*)
//
// cell4 = -1;   (*Prior  -50, 0. It is strictly negative BUT WE DO NOT \
//WANT TO FIX TO HEXADECAPOLE. TOO NOISY*)
//
// ctilde = -0.25;  (*Prior  -1, 1   but small steps ~ 0.1*)

typedef struct {
    int nparams;
    real p1;        // b1
    real p1cv;
    real p1min;
    real p1max;
    real dp1;

    real p2;        // b2
    real p2cv;
    real p2min;
    real p2max;
    real dp2;

    real p3;        // cell0
    real p3cv;
    real p3min;
    real p3max;
    real dp3;

    real p4;        // cell2
    real p4cv;
    real p4min;
    real p4max;
    real dp4;

    real p5;        // cell4
    real p5cv;
    real p5min;
    real p5max;
    real dp5;

    real p6;        // ctilde
    real p6cv;
    real p6min;
    real p6max;
    real dp6;
    
    real bs2;
    real brnl;

} model_data_user, *model_data_user_ptr;

local model_data_user mduser;

#define md          mduser

// params[1] === b1     :: (1.45) range: 1 < p1 < 2.5
// params[2] === b2     :: (-0.2) range: -1 < p2 < 1
// params[3] === cell0  :: (-12) range: -40 < p3 < 0
// params[4] === cell2  :: (-32) range: -60 < p4 < 0
// params[5] === cell4  :: (-1) range: -50 < p5 < 0
// params[6] === ctilde :: (-0.25) range: -1 < p6 < 1, dp=0.1

//#define b1      params[1]
//#define b2      params[2]
//#define cell0   params[3]
//#define cell2   params[4]
//#define cell4   params[5]
//#define ctilde  params[6]

local real b1;
local real b2;
local real cell0;
local real cell2;
local real cell4;
local real ctilde;

local real bs2;
local real b3nl;


local void set_Model_USER(void)
{
    strcpy(gd.model_comment, "USER Model");
    
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
#define usrmodel_parameter_null    "usrmodel_parameters_null"
        fprintf(stdout,
        "set_Model_USER: no parameter file was given...\nUsing a default set...'\n");
        pd.nparams = 2;
        pd.params = dvector(1,pd.nparams);
        pd.dparams = dvector(1,pd.nparams);
        pd.rparams = dmatrix(1,pd.nparams,1,2);
        pd.params[1] = 4899.44;
        pd.rparams[1][1] = 0.0001;
        pd.rparams[1][2] = 10000.0;
        pd.dparams[1] = 10.0; // 10.0
        pd.params[2] = 3.65868;
        pd.rparams[2][1] = 0.0001;
        pd.rparams[2][2] = 10.0;
        pd.dparams[2] = 0.01;  // 0.01

        md.nparams = pd.nparams;
        md.p1cv = pd.params[1];
        md.p1min = pd.rparams[1][1];
        md.p1max = pd.rparams[1][2];
        md.dp1 = pd.dparams[1];
        md.p2cv = pd.params[2];
        md.p2min = pd.rparams[2][1];
        md.p2max = pd.rparams[2][2];
        md.dp2 = pd.dparams[2];

        PrintModelParameterFile_user(usrmodel_parameter_null);
#undef usrmodel_parameter_null
    }

    gd.modelfun = &rcUSER;
    gd.modelfunwd = &model_wdyda;
//    gd.modelfunwd = &model_wd1yda;
    gd.modelChi2 = &Chi2_USER;
    gd.modelend = &Model_USER_end;
    
    read_allfunctions();
}

local void test_USER(void)
{
    real Chi2tmp;
    real Chi2red;

    fprintf(stdout,"\n\nTesting %s...\n",gd.model_comment);

    Chi2tmp = gd.modelChi2(gd.modelfun, pd.params);
    Chi2red = Chi2tmp/(nObs(pObs[gd.filecount])-pd.nparams);
    fprintf(stdout,"Chi2: %g %g %g %g\n",pd.params[1],pd.params[2],Chi2tmp,Chi2red);
    plotting_model_table(gd.modelfun,gd.fpfnamemodeltest);

    fprintf(stdout,"done.\n\n");
}

global real Chi2_USER(real (*ymodel)(real, real *), real params[])
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
                         -model_Interp(xObs(pObs[gd.filecount])[i])
                         )
                        /sigmaObs(pObs[gd.filecount])[i]
                        );
    }
//    (yobs[i] - ymodel(xobs[i], params))/sigma[i]

    return Chi2tmp;
}

#define mparamfinal     "modelparamfinal"
global void Model_USER_end(void)
{
    real **covm;
    int i, j, k;
    int *indx;
    real d;
    real *col;
    real **alpha, **atmp;
    real sum;
    real **ident;
    real *b, *x;

    alpha = dmatrix(1,pd.nparams,1,pd.nparams);
    atmp = dmatrix(1,pd.nparams,1,pd.nparams);
    covm = dmatrix(1,pd.nparams,1,pd.nparams);
    ident = dmatrix(1,pd.nparams,1,pd.nparams);
    col = dvector(1,pd.nparams);
    indx = ivector(1,pd.nparams);
    b = dvector(1,pd.nparams);
    x = dvector(1,pd.nparams);

    md.nparams = pd.nparams;
    md.p1cv = pd.params[1];
    md.p1min = pd.rparams[1][1];
    md.p1max = pd.rparams[1][2];
    md.dp1 = pd.dparams[1];
    md.p2cv = pd.params[2];
    md.p2min = pd.rparams[2][1];
    md.p2max = pd.rparams[2][2];
    md.dp2 = pd.dparams[2];
    
    PrintModelParameterFile_user(mparamfinal);


    model_wd2yda_chi2(pd.params, alpha, pd.nparams);
    
    fprintf(stdout,"\nalpha matrix:\n");
    for (i=1; i<=pd.nparams; i++) {
        for (j=1; j<=pd.nparams; j++) {
            fprintf(stdout,"%e ", alpha[i][j]);
        }
        fprintf(stdout,"\n");
    }
    for (i=1; i<=pd.nparams; i++)
        for (j=1; j<=pd.nparams; j++)
            atmp[i][j]=alpha[i][j];

    ludcmp(atmp, pd.nparams, indx, &d);
    for (j=1; j<=pd.nparams; j++) {
        for (i=1; i<=pd.nparams; i++) col[i]=0.;
        col[j] = 1.0;
        lubksb(atmp, pd.nparams, indx, col);
        for (i=1; i<=pd.nparams; i++) covm[i][j]=col[i];
    }

    fprintf(stdout,"\ncovm matrix:\n");
    for (i=1; i<=pd.nparams; i++) {
        for (j=1; j<=pd.nparams; j++) {
            fprintf(stdout,"%e ", covm[i][j]);
        }
        fprintf(stdout,"\n");
    }

// Check inverse:
    for (i=1; i<=pd.nparams; i++) {
        for (j=1; j<=pd.nparams; j++) {
            sum=0.0;
            for (k=1; k<=pd.nparams; k++) sum += alpha[i][k]*covm[k][j];
            ident[i][j] = sum;
        }
    }
    fprintf(stdout,"\nproduct matrix:\n");
    for (i=1; i<=pd.nparams; i++) {
        for (j=1; j<=pd.nparams; j++) {
            fprintf(stdout,"%e ", ident[i][j]);
        }
        fprintf(stdout,"\n");
    }


// Cholesky decomposition
    for (i=1; i<=pd.nparams; i++)
        for (j=1; j<=pd.nparams; j++)
            atmp[i][j]=alpha[i][j];
    choldc(atmp, pd.nparams, col);
    for (j=1; j<=pd.nparams; j++) {
        for (i=1; i<=pd.nparams; i++) b[i]=0.;
        b[j] = 1.0;
        cholsl(atmp, pd.nparams, col, b, x);
        for (i=1; i<=pd.nparams; i++) covm[i][j]=x[i];
    }

/*
// LU decomposition: A = L . U
// Cholesky: A = L . L^(T), where [L^(T)]_ij = L_ji
// Obtain the lower triangle of L^-1
    for (i=1; i<=pd.nparams; i++) {
        covm[i][i]=1.0/col[i];
        for (j=i+1; j<=pd.nparams; j++) {
            sum=0.0;
            for (k=i; k<j; k++) sum -= atmp[j][k]*atmp[k][i];
            covm[j][i] = sum/col[j];
        }
    }
*/

    fprintf(stdout,"\ncovm matrix Cholesky decomposition:\n");
    for (i=1; i<=pd.nparams; i++) {
        for (j=1; j<=pd.nparams; j++) {
            fprintf(stdout,"%e ", covm[i][j]);
        }
        fprintf(stdout,"\n");
    }

// Check inverse:
    for (i=1; i<=pd.nparams; i++) {
        for (j=1; j<=pd.nparams; j++) {
            sum=0.0;
            for (k=1; k<=pd.nparams; k++) sum += alpha[i][k]*covm[k][j];
            ident[i][j] = sum;
        }
    }
    fprintf(stdout,"\nproduct matrix:\n");
    for (i=1; i<=pd.nparams; i++) {
        for (j=1; j<=pd.nparams; j++) {
            fprintf(stdout,"%e ", ident[i][j]);
        }
        fprintf(stdout,"\n");
    }

// Errors:
    fprintf(stdout,"\nErrors:\n");
    for (i=1; i<=pd.nparams; i++)
        fprintf(stdout,"%e\n", rsqrt(covm[i][i]) );

    free_dvector(x,1,pd.nparams);
    free_dvector(b,1,pd.nparams);
    free_ivector(indx,1,pd.nparams);
    free_dvector(col,1,pd.nparams);
    free_dmatrix(ident,1,pd.nparams,1,pd.nparams);
    free_dmatrix(covm,1,pd.nparams,1,pd.nparams);
    free_dmatrix(atmp,1,pd.nparams,1,pd.nparams);
    free_dmatrix(alpha,1,pd.nparams,1,pd.nparams);
}
#undef mparamfinal


#define rhos    params[1]
#define rs      params[2]
global real rcUSER(real x, real params[])
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
#include "model_user_io.h"

#undef READPARAMS
#undef WRITEPARAMS

#undef md

//#define suffixModel     "GRz05"
//#define suffixModel    "parameters_null-mgpt"

//int nPSLT;
//real *kPS;
//real *pPS;
//real *pPS2;

//int nfk;
//real *kfk;
//real *pfk;
//real *pfk2;

local void read_allfunctions(void)
{
    int i;
    char fpfnameallfun[100];
    char suffixModel[100];
    strcpy(suffixModel,"GRz05");

    fprintf(stdout,"\n\nReading all functions...\n");

// ------------------------ 1 and 2 --------------------------------------------
    sprintf(fpfnameallfun,"Input/TheoryFiles/linear/PSL_%s.dat",suffixModel);
    fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 2, &npslfunT);

    if (npslfunT <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", npslfunT);
    fprintf(gd.outlog,"\nRead %d points...",npslfunT);

    kpslT = dvector(1,npslfunT);
    pslT = dvector(1,npslfunT);
    pslT2 = dvector(1,npslfunT);

    for (i=0; i<npslfunT; i++) {
        kpslT[i+1] = inout_xval[i];
        pslT[i+1] = inout_yval[i];
    }
    spline(kpslT,pslT,npslfunT,1.0e30,1.0e30,pslT2);
// -----------------------------------------------------------------------------

// ------------------------ 1 and 2 --------------------------------------------
    sprintf(fpfnameallfun,"Input/TheoryFiles/linear/fk_%s.dat",suffixModel);
    fprintf(gd.outlog,"\n\nread_allfunctions: Reading fkT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 2, &nfkfunT);

    if (nfkfunT <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nfkfunT);
    fprintf(gd.outlog,"\nRead %d points...",nfkfunT);

    kfkT = dvector(1,nfkfunT);
    fkT = dvector(1,nfkfunT);
    fkT2 = dvector(1,nfkfunT);

    for (i=0; i<nfkfunT; i++) {
//        fprintf(stdout,"\ninout values: %g %g", inout_xval[i], inout_yval[i]);
        kfkT[i+1] = inout_xval[i];
        fkT[i+1] = inout_yval[i];
    }
    spline(kfkT,fkT,nfkfunT,1.0e30,1.0e30,fkT2);
// -----------------------------------------------------------------------------

/*
    kPS = dvector(1,nPSLT);
    pPS = dvector(1,nPSLT);
    pPS2 = dvector(1,nPSLT);
    kfk = dvector(1,nfk);
    pfk = dvector(1,nfk);
    pfk2 = dvector(1,nfk);
    spline(kPS,pPS,nPSLT,1.0e30,1.0e30,pPS2);
    spline(kPS,pfk,nPSLT,1.0e30,1.0e30,pfk2);
*/

    M1f0 = fkT[1];
    M1sigma2v = sigma2v_function();
    M1sigma2Psi = sigma2Psi_function();
    fprintf(stdout,"\nf0, sigma2v, sigma2Psi: %g %g %g",M1f0, M1sigma2v, M1sigma2Psi);
    fprintf(stdout,"\n1/sqrt(sigma2v), 1/sqrt(sigma2Psi): %g %g\n",
            1/rsqrt(M1sigma2v), 1/rsqrt(M1sigma2Psi));

//==============================================================================
// M1monopole:

// ------------------------ 1 and 2 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
    fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 2, &nM1MnbT);

    if (nM1MnbT <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1MnbT);
    fprintf(gd.outlog,"\nRead %d points...",nM1MnbT);

    kM1MnbT = dvector(1,nM1MnbT);
    M1MnbT = dvector(1,nM1MnbT);
    M1MnbT2 = dvector(1,nM1MnbT);

    for (i=0; i<nM1MnbT; i++) {
        kM1MnbT[i+1] = inout_xval[i];
        M1MnbT[i+1] = inout_yval[i];
    }
    spline(kM1MnbT,M1MnbT,nM1MnbT,1.0e30,1.0e30,M1MnbT2);
// -----------------------------------------------------------------------------

// ------------------------ 1 and 3 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 3, &nM1Mb1T);

    if (nM1Mb1T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb1T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb1T);

    kM1Mb1T = dvector(1,nM1Mb1T);
    M1Mb1T = dvector(1,nM1Mb1T);
    M1Mb1T2 = dvector(1,nM1Mb1T);

    for (i=0; i<nM1Mb1T; i++) {
        kM1Mb1T[i+1] = inout_xval[i];
        M1Mb1T[i+1] = inout_yval[i];
    }
    spline(kM1Mb1T,M1Mb1T,nM1Mb1T,1.0e30,1.0e30,M1Mb1T2);
// -----------------------------------------------------------------------------

// ------------------------ 1 and 4 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 4, &nM1Mb12T);

    if (nM1Mb12T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb12T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb12T);

    kM1Mb12T = dvector(1,nM1Mb12T);
    M1Mb12T = dvector(1,nM1Mb12T);
    M1Mb12T2 = dvector(1,nM1Mb12T);

    for (i=0; i<nM1Mb12T; i++) {
        kM1Mb12T[i+1] = inout_xval[i];
        M1Mb12T[i+1] = inout_yval[i];
    }
    spline(kM1Mb12T,M1Mb12T,nM1Mb12T,1.0e30,1.0e30,M1Mb12T2);

// -----------------------------------------------------------------------------


// ------------------------ 1 and 5 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 5, &nM1Mb2T);

    if (nM1Mb2T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb2T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb2T);

    kM1Mb2T = dvector(1,nM1Mb2T);
    M1Mb2T = dvector(1,nM1Mb2T);
    M1Mb2T2 = dvector(1,nM1Mb2T);

    for (i=0; i<nM1Mb2T; i++) {
        kM1Mb2T[i+1] = inout_xval[i];
        M1Mb2T[i+1] = inout_yval[i];
    }
    spline(kM1Mb2T,M1Mb2T,nM1Mb2T,1.0e30,1.0e30,M1Mb2T2);


// -----------------------------------------------------------------------------

    
// ------------------------ 1 and 6 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 6, &nM1Mb22T);

    if (nM1Mb22T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb22T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb22T);

    kM1Mb22T = dvector(1,nM1Mb22T);
    M1Mb22T = dvector(1,nM1Mb22T);
    M1Mb22T2 = dvector(1,nM1Mb22T);

    for (i=0; i<nM1Mb22T; i++) {
        kM1Mb22T[i+1] = inout_xval[i];
        M1Mb22T[i+1] = inout_yval[i];
    }
    spline(kM1Mb22T,M1Mb22T,nM1Mb22T,1.0e30,1.0e30,M1Mb22T2);


// -----------------------------------------------------------------------------

// ------------------------ 1 and 7 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 7, &nM1Mb2b1T);

    if (nM1Mb2b1T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb2b1T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb2b1T);

    kM1Mb2b1T = dvector(1,nM1Mb2b1T);
    M1Mb2b1T = dvector(1,nM1Mb2b1T);
    M1Mb2b1T2 = dvector(1,nM1Mb2b1T);

    for (i=0; i<nM1Mb2b1T; i++) {
        kM1Mb2b1T[i+1] = inout_xval[i];
        M1Mb2b1T[i+1] = inout_yval[i];
    }
    spline(kM1Mb2b1T,M1Mb2b1T,nM1Mb2b1T,1.0e30,1.0e30,M1Mb2b1T2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 8 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 8, &nM1Mb2bs2T);

    if (nM1Mb2bs2T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb2bs2T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb2bs2T);

    kM1Mb2bs2T = dvector(1,nM1Mb2bs2T);
    M1Mb2bs2T = dvector(1,nM1Mb2bs2T);
    M1Mb2bs2T2 = dvector(1,nM1Mb2bs2T);

    for (i=0; i<nM1Mb2bs2T; i++) {
        kM1Mb2bs2T[i+1] = inout_xval[i];
        M1Mb2bs2T[i+1] = inout_yval[i];
    }
    spline(kM1Mb2bs2T,M1Mb2bs2T,nM1Mb2bs2T,1.0e30,1.0e30,M1Mb2bs2T2);


// -----------------------------------------------------------------------------

// ------------------------ 1 and 9 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 9, &nM1Mbs2T);

    if (nM1Mbs2T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mbs2T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mbs2T);

    kM1Mbs2T = dvector(1,nM1Mbs2T);
    M1Mbs2T = dvector(1,nM1Mbs2T);
    M1Mbs2T2 = dvector(1,nM1Mbs2T);

    for (i=0; i<nM1Mbs2T; i++) {
        kM1Mbs2T[i+1] = inout_xval[i];
        M1Mbs2T[i+1] = inout_yval[i];
    }
    spline(kM1Mbs2T,M1Mbs2T,nM1Mbs2T,1.0e30,1.0e30,M1Mbs2T2);


// -----------------------------------------------------------------------------

// ------------------------ 1 and 10 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 10, &nM1Mbs22T);

    if (nM1Mbs22T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mbs22T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mbs22T);

    kM1Mbs22T = dvector(1,nM1Mbs22T);
    M1Mbs22T = dvector(1,nM1Mbs22T);
    M1Mbs22T2 = dvector(1,nM1Mbs22T);

    for (i=0; i<nM1Mbs22T; i++) {
        kM1Mbs22T[i+1] = inout_xval[i];
        M1Mbs22T[i+1] = inout_yval[i];
    }
    spline(kM1Mbs22T,M1Mbs22T,nM1Mbs22T,1.0e30,1.0e30,M1Mbs22T2);


// -----------------------------------------------------------------------------

// ------------------------ 1 and 11 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 11, &nM1Mbs2b1T);

    if (nM1Mbs2b1T <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mbs2b1T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mbs2b1T);

    kM1Mbs2b1T = dvector(1,nM1Mbs2b1T);
    M1Mbs2b1T = dvector(1,nM1Mbs2b1T);
    M1Mbs2b1T2 = dvector(1,nM1Mbs2b1T);

    for (i=0; i<nM1Mbs2b1T; i++) {
        kM1Mbs2b1T[i+1] = inout_xval[i];
        M1Mbs2b1T[i+1] = inout_yval[i];
    }
    spline(kM1Mbs2b1T,M1Mbs2b1T,nM1Mbs2b1T,1.0e30,1.0e30,M1Mbs2b1T2);


// -----------------------------------------------------------------------------

// ------------------------ 1 and 12 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 12, &nM1Mb1b3nlT);

    if (nM1Mb1b3nlT <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb1b3nlT);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb1b3nlT);

    kM1Mb1b3nlT = dvector(1,nM1Mb1b3nlT);
    M1Mb1b3nlT = dvector(1,nM1Mb1b3nlT);
    M1Mb1b3nlT2 = dvector(1,nM1Mb1b3nlT);

    for (i=0; i<nM1Mb1b3nlT; i++) {
        kM1Mb1b3nlT[i+1] = inout_xval[i];
        M1Mb1b3nlT[i+1] = inout_yval[i];
    }
    spline(kM1Mb1b3nlT,M1Mb1b3nlT,nM1Mb1b3nlT,1.0e30,1.0e30,M1Mb1b3nlT2);


// -----------------------------------------------------------------------------

// ------------------------ 1 and 13 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 13, &nM1Mb3nlT);

    if (nM1Mb3nlT <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1Mb3nlT);
    fprintf(gd.outlog,"\nRead %d points...",nM1Mb3nlT);

    kM1Mb3nlT = dvector(1,nM1Mb3nlT);
    M1Mb3nlT = dvector(1,nM1Mb3nlT);
    M1Mb3nlT2 = dvector(1,nM1Mb3nlT);

    for (i=0; i<nM1Mb3nlT; i++) {
        kM1Mb3nlT[i+1] = inout_xval[i];
        M1Mb3nlT[i+1] = inout_yval[i];
    }
    spline(kM1Mb3nlT,M1Mb3nlT,nM1Mb3nlT,1.0e30,1.0e30,M1Mb3nlT2);


// -----------------------------------------------------------------------------

// ------------------------ 1 and 14 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell0_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading pslT from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 14, &nM1MeftT);

    if (nM1MeftT <= 1)
        error("\nInputPSLfunctionTable: nallfunT = %d is absurd\n\n", nM1MeftT);
    fprintf(gd.outlog,"\nRead %d points...",nM1MeftT);

    kM1MeftT = dvector(1,nM1MeftT);
    M1MeftT = dvector(1,nM1MeftT);
    M1MeftT2 = dvector(1,nM1MeftT);

    for (i=0; i<nM1MeftT; i++) {
        kM1MeftT[i+1] = inout_xval[i];
        M1MeftT[i+1] = inout_yval[i];
    }
    spline(kM1MeftT,M1MeftT,nM1MeftT,1.0e30,1.0e30,M1MeftT2);
// -----------------------------------------------------------------------------
//==============================================================================

//==============================================================================
// Quadrupoles:

// ------------------------ 1 and 2 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell2_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Qnb from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 2, &nM1QnbT);

    if (nM1QnbT <= 1)
        error("\nInput M1Qnb Table: nM1QnbT = %d is absurd\n\n", nM1QnbT);
    fprintf(gd.outlog,"\nRead %d points...",nM1QnbT);

    kM1QnbT = dvector(1,nM1QnbT);
    M1QnbT = dvector(1,nM1QnbT);
    M1QnbT2 = dvector(1,nM1QnbT);

    for (i=0; i<nM1QnbT; i++) {
        kM1QnbT[i+1] = inout_xval[i];
        M1QnbT[i+1] = inout_yval[i];
    }
    spline(kM1QnbT,M1QnbT,nM1QnbT,1.0e30,1.0e30,M1QnbT2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 3 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell2_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Qb1 from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 3, &nM1Qb1T);

    if (nM1Qb1T <= 1)
        error("\nInput M1Qb1 Table: nM1Qb1T = %d is absurd\n\n", nM1Qb1T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Qb1T);

    kM1Qb1T = dvector(1,nM1Qb1T);
    M1Qb1T = dvector(1,nM1Qb1T);
    M1Qb1T2 = dvector(1,nM1Qb1T);

    for (i=0; i<nM1Qb1T; i++) {
        kM1Qb1T[i+1] = inout_xval[i];
        M1Qb1T[i+1] = inout_yval[i];
    }
    spline(kM1Qb1T,M1Qb1T,nM1Qb1T,1.0e30,1.0e30,M1Qb1T2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 4 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell2_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Qb12 from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 4, &nM1Qb12T);

    if (nM1Qb12T <= 1)
        error("\nInput M1Qb12 Table: nM1Qb12T = %d is absurd\n\n", nM1Qb12T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Qb12T);

    kM1Qb12T = dvector(1,nM1Qb12T);
    M1Qb12T = dvector(1,nM1Qb12T);
    M1Qb12T2 = dvector(1,nM1Qb12T);

    for (i=0; i<nM1Qb12T; i++) {
        kM1Qb12T[i+1] = inout_xval[i];
        M1Qb12T[i+1] = inout_yval[i];
    }
    spline(kM1Qb12T,M1Qb12T,nM1Qb12T,1.0e30,1.0e30,M1Qb12T2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 5 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell2_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Qb2 from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 5, &nM1Qb2T);

    if (nM1Qb2T <= 1)
        error("\nInput M1Qb2 Table: nM1Qb2T = %d is absurd\n\n", nM1Qb2T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Qb2T);

    kM1Qb2T = dvector(1,nM1Qb2T);
    M1Qb2T = dvector(1,nM1Qb2T);
    M1Qb2T2 = dvector(1,nM1Qb2T);

    for (i=0; i<nM1Qb2T; i++) {
        kM1Qb2T[i+1] = inout_xval[i];
        M1Qb2T[i+1] = inout_yval[i];
    }
    spline(kM1Qb2T,M1Qb2T,nM1Qb2T,1.0e30,1.0e30,M1Qb2T2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 6 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell2_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Qbs2 from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 6, &nM1Qbs2T);

    if (nM1Qbs2T <= 1)
        error("\nInput M1Qbs2 Table: nM1Qbs2T = %d is absurd\n\n", nM1Qbs2T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Qbs2T);

    kM1Qbs2T = dvector(1,nM1Qbs2T);
    M1Qbs2T = dvector(1,nM1Qbs2T);
    M1Qbs2T2 = dvector(1,nM1Qbs2T);

    for (i=0; i<nM1Qbs2T; i++) {
        kM1Qbs2T[i+1] = inout_xval[i];
        M1Qbs2T[i+1] = inout_yval[i];
    }
    spline(kM1Qbs2T,M1Qbs2T,nM1Qbs2T,1.0e30,1.0e30,M1Qbs2T2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 7 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell2_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Qb3nl from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 7, &nM1Qb3nlT);

    if (nM1Qb3nlT <= 1)
        error("\nInput M1Qb3nl Table: nM1Qb3nlT = %d is absurd\n\n", nM1Qb3nlT);
    fprintf(gd.outlog,"\nRead %d points...",nM1Qb3nlT);

    kM1Qb3nlT = dvector(1,nM1Qb3nlT);
    M1Qb3nlT = dvector(1,nM1Qb3nlT);
    M1Qb3nlT2 = dvector(1,nM1Qb3nlT);

    for (i=0; i<nM1Qb3nlT; i++) {
        kM1Qb3nlT[i+1] = inout_xval[i];
        M1Qb3nlT[i+1] = inout_yval[i];
    }
    spline(kM1Qb3nlT,M1Qb3nlT,nM1Qb3nlT,1.0e30,1.0e30,M1Qb3nlT2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 8 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell2_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Qeft from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 8, &nM1QeftT);

    if (nM1QeftT <= 1)
        error("\nInput M1Qeft Table: nM1QeftT = %d is absurd\n\n", nM1QeftT);
    fprintf(gd.outlog,"\nRead %d points...",nM1QeftT);

    kM1QeftT = dvector(1,nM1QeftT);
    M1QeftT = dvector(1,nM1QeftT);
    M1QeftT2 = dvector(1,nM1QeftT);

    for (i=0; i<nM1QeftT; i++) {
        kM1QeftT[i+1] = inout_xval[i];
        M1QeftT[i+1] = inout_yval[i];
    }
    spline(kM1QeftT,M1QeftT,nM1QeftT,1.0e30,1.0e30,M1QeftT2);

// -----------------------------------------------------------------------------
//==============================================================================

//==============================================================================
// Hexapoles:

// ------------------------ 1 and 2 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell4_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Hnb from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 2, &nM1HnbT);

    if (nM1HnbT <= 1)
        error("\nInput M1Hnb Table: nM1HnbT = %d is absurd\n\n", nM1HnbT);
    fprintf(gd.outlog,"\nRead %d points...",nM1HnbT);

    kM1HnbT = dvector(1,nM1HnbT);
    M1HnbT = dvector(1,nM1HnbT);
    M1HnbT2 = dvector(1,nM1HnbT);

    for (i=0; i<nM1HnbT; i++) {
        kM1HnbT[i+1] = inout_xval[i];
        M1HnbT[i+1] = inout_yval[i];
    }
    spline(kM1HnbT,M1HnbT,nM1HnbT,1.0e30,1.0e30,M1HnbT2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 3 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell4_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Hb1 from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 3, &nM1Hb1T);

    if (nM1Hb1T <= 1)
        error("\nInput M1Hb1 Table: nM1Hb1T = %d is absurd\n\n", nM1Hb1T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Hb1T);

    kM1Hb1T = dvector(1,nM1Hb1T);
    M1Hb1T = dvector(1,nM1Hb1T);
    M1Hb1T2 = dvector(1,nM1Hb1T);

    for (i=0; i<nM1Hb1T; i++) {
        kM1Hb1T[i+1] = inout_xval[i];
        M1Hb1T[i+1] = inout_yval[i];
    }
    spline(kM1Hb1T,M1Hb1T,nM1Hb1T,1.0e30,1.0e30,M1Hb1T2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 4 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell4_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Hb12 from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 4, &nM1Hb12T);

    if (nM1Hb12T <= 1)
        error("\nInput M1Hb12 Table: nM1Hb12T = %d is absurd\n\n", nM1Hb12T);
    fprintf(gd.outlog,"\nRead %d points...",nM1Hb12T);

    kM1Hb12T = dvector(1,nM1Hb12T);
    M1Hb12T = dvector(1,nM1Hb12T);
    M1Hb12T2 = dvector(1,nM1Hb12T);

    for (i=0; i<nM1Hb12T; i++) {
        kM1Hb12T[i+1] = inout_xval[i];
        M1Hb12T[i+1] = inout_yval[i];
    }
    spline(kM1Hb12T,M1Hb12T,nM1Hb12T,1.0e30,1.0e30,M1Hb12T2);

// -----------------------------------------------------------------------------

// ------------------------ 1 and 5 --------------------------------------------
    sprintf(fpfnameallfun,
            "Input/TheoryFiles/Ps_MomentumExpansion_ell4_%s.dat",suffixModel);
        fprintf(gd.outlog,"\n\nread_allfunctions: Reading M1Heft from file %s...",
        fpfnameallfun);
    inout_InputData(fpfnameallfun, 1, 5, &nM1HeftT);

    if (nM1HeftT <= 1)
        error("\nInput M1Heft Table: nM1HeftT = %d is absurd\n\n", nM1HeftT);
    fprintf(gd.outlog,"\nRead %d points...",nM1HeftT);

    kM1HeftT = dvector(1,nM1HeftT);
    M1HeftT = dvector(1,nM1HeftT);
    M1HeftT2 = dvector(1,nM1HeftT);

    for (i=0; i<nM1HeftT; i++) {
        kM1HeftT[i+1] = inout_xval[i];
        M1HeftT[i+1] = inout_yval[i];
    }
    spline(kM1HeftT,M1HeftT,nM1HeftT,1.0e30,1.0e30,M1HeftT2);

// -----------------------------------------------------------------------------
//==============================================================================
}

//==============================================================================
// (*Construct Multipoles*)

local real PM1MEell0(real k, real b1, real b2, real bs2, real b3nl, real cell0)
{
    real Stmp;
    Stmp =   M1Mnb(k) + b1* M1Mb1(k) + rsqr(b1)* M1Mb12(k) + b2* M1Mb2(k)
    + rsqr(b2)* M1Mb22(k) + b2*b1* M1Mb2b1(k) + b2*bs2* M1Mb2bs2(k)
    + bs2* M1Mbs2(k) + rsqr(bs2)* M1Mbs22(k) + bs2*b1* M1Mbs2b1(k)
    + b1*b3nl* M1Mb1b3nl(k) + b3nl* M1Mb3nl(k) + cell0* M1Meft(k);

    return (Stmp);
}

local real PM1MEell2(real k, real b1, real b2, real bs2, real b3nl, real cell2)
{
    real Stmp;
    Stmp =  b1* M1Qb1(k) + rsqr(b1)* M1Qb12(k) + b2* M1Qb2(k) + bs2* M1Qbs2(k) +
   b3nl* M1Qb3nl(k) + M1Qnb(k) + cell2* M1Qeft(k);
    return (Stmp);
}

local real PM1MEell4(real k, real b1, real b2, real bs2, real b3nl, real cell4)
{
    real Stmp;

    Stmp = b1* M1Hb1(k) + rsqr(b1)* M1Hb12(k) + M1Hnb(k) + cell4* M1Heft(k);
    
    return (Stmp);
}

local real DNLOell0(real k, real b1)
{
    real Stmp;

    Stmp = rsqr(b1)/5.0 + (2.0* b1* M1fk(k))/7.0 + rsqr(M1fk(k))/9.0;
    
    return (Stmp);
}

local real DNLOell2(real k, real b1)
{
    real Stmp;

    Stmp = (4.0* rsqr(b1))/7.0 + (20.0* b1* M1fk(k))/21.0 + (40.0* rsqr(M1fk(k)))/99.0;
    
    return (Stmp);
}

local real DNLOell4(real k, real b1)
{
    real Stmp;

    Stmp = (8.0* rsqr(b1))/35.0 + (48.0* b1* M1fk(k))/77.0
    + (48.0* rsqr(M1fk(k)))/143.0;
    
    return (Stmp);
}

local real M1MonopoleME(real k)
{
    real Stmp;

    Stmp = PM1MEell0(k, b1, b2, bs2, b3nl, cell0)
    + ctilde* DNLOell0(k, b1) * rpow(M1f0,4.0)*  rpow(k,4.0)* rsqr(M1sigma2v)* M1pk(k);

    return (Stmp);
}

local real M1QuadrupoleME(real k)
{
    real Stmp;

    Stmp = PM1MEell2(k, b1, b2, bs2, b3nl, cell2)
    + ctilde* DNLOell2(k, b1)* rpow(M1f0,4.0)* rpow(k,4.0)* rsqr(M1sigma2v)* M1pk(k);

    return (Stmp);
}

local real M1HexadecapoleME(real k)
{
    real Stmp;

    Stmp = PM1MEell4(k, b1, b2, bs2, b3nl, cell4) +
     ctilde* DNLOell4(k, b1)* rpow(M1f0,4.0)* rpow(k,4.0)* rsqr(M1sigma2v)* M1pk(k);

    return (Stmp);
}

local real factor(real k)
{
    real Stmp;

    Stmp = rpow(k,1.5);

    return (Stmp);
}

//==============================================================================


local real M1pk(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kpslT, pslT, npslfunT, pslT2);
    return (Stmp);
}

local real M1fk(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kfkT, fkT, nfkfunT, fkT2);
    return (Stmp);
}

local real M1Mnb(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1MnbT, M1MnbT, nM1MnbT, M1MnbT2);
    return (Stmp);
}

local real M1Mb1(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb1T, M1Mb1T, nM1Mb1T, M1Mb1T2);
    return (Stmp);
}

local real M1Mb12(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb12T, M1Mb12T, nM1Mb12T, M1Mb12T2);
    return (Stmp);
}

local real M1Mb2(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb2T, M1Mb2T, nM1Mb2T, M1Mb2T2);
    return (Stmp);
}

local real M1Mb22(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb22T, M1Mb22T, nM1Mb22T, M1Mb22T2);
    return (Stmp);
}

local real M1Mb2b1(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb2b1T, M1Mb2b1T, nM1Mb2b1T, M1Mb2b1T2);
    return (Stmp);
}

local real M1Mb2bs2(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb2bs2T, M1Mb2bs2T, nM1Mb2bs2T, M1Mb2bs2T2);
    return (Stmp);
}

local real M1Mbs2(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mbs2T, M1Mbs2T, nM1Mbs2T, M1Mbs2T2);
    return (Stmp);
}

local real M1Mbs22(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mbs22T, M1Mbs22T, nM1Mbs22T, M1Mbs22T2);
    return (Stmp);
}

local real M1Mbs2b1(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mbs2b1T, M1Mbs2b1T, nM1Mbs2b1T, M1Mbs2b1T2);
    return (Stmp);
}

local real M1Mb1b3nl(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb1b3nlT, M1Mb1b3nlT, nM1Mb1b3nlT, M1Mb1b3nlT2);
    return (Stmp);
}

local real M1Mb3nl(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Mb3nlT, M1Mb3nlT, nM1Mb3nlT, M1Mb3nlT2);
    return (Stmp);
}

local real M1Meft(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1MeftT, M1MeftT, nM1MeftT, M1MeftT2);
    return (Stmp);
}


//==============================================================================
// Quadrupoles:

local real M1Qnb(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1QnbT, M1QnbT, nM1QnbT, M1QnbT2);
    return (Stmp);
}

local real M1Qb1(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Qb1T, M1Qb1T, nM1Qb1T, M1Qb1T2);
    return (Stmp);
}

local real M1Qb12(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Qb12T, M1Qb12T, nM1Qb12T, M1Qb12T2);
    return (Stmp);
}

local real M1Qb2(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Qb2T, M1Qb2T, nM1Qb2T, M1Qb2T2);
    return (Stmp);
}

local real M1Qbs2(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Qbs2T, M1Qbs2T, nM1Qbs2T, M1Qbs2T2);
    return (Stmp);
}

local real M1Qb3nl(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Qb3nlT, M1Qb3nlT, nM1Qb3nlT, M1Qb3nlT2);
    return (Stmp);
}

local real M1Qeft(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1QeftT, M1QeftT, nM1QeftT, M1QeftT2);
    return (Stmp);
}
//==============================================================================

//==============================================================================
// Hexapoles:

local real M1Hnb(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1HnbT, M1HnbT, nM1HnbT, M1HnbT2);
    return (Stmp);
}

local real M1Hb1(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Hb1T, M1Hb1T, nM1Hb1T, M1Hb1T2);
    return (Stmp);
}

local real M1Hb12(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1Hb12T, M1Hb12T, nM1Hb12T, M1Hb12T2);
    return (Stmp);
}

local real M1Heft(real kv){
    real Stmp;
    Stmp = Interpolation_nr(kv, kM1HeftT, M1HeftT, nM1HeftT, M1HeftT2);
    return (Stmp);
}

//==============================================================================


#define KK  5
#define EPSQ 1.0e-6

// sigma2v :: integration over the extended MG power spectrum
local real sigma2v_function_int(real y)
{
    real p;
    real fkv, PSL;
    real Stmp;

    p = rpow(10.0,y);
//    p = y;
    fkv = Interpolation_nr(p, kfkT, fkT, nfkfunT, fkT2);
    PSL = Interpolation_nr(p, kpslT, pslT, npslfunT, pslT2);
//    fkv = fkInterpolation_nr(p, gd.kfk, gd.pfk, gd.nfk);
//    PSL = psInterpolation_nr(p, gd.kPS, gd.pPS, gd.nPSLT);

    Stmp = rsqr(fkv/M1f0)*PSL;
    
//    printf("\nStmp: %g", Stmp);

    return p*Stmp;
}

local real sigma2v_function(void)
{
    real result;
    real kmin, kmax;
    real ymin, ymax;

    kmin = kpslT[1];
    kmax = kpslT[npslfunT];
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma2v_function_int,ymin,ymax,midpnt,EPSQ,KK);
//    result= (1.0/SIXPI2)
//    *qromo(sigma2v_function_int,kmin,kmax,midpnt,EPSQ,KK);

    return result;

}

// sigma2Psi :: integration over the extended MG power spectrum
local real sigma2Psi_function_int(real y)
{
    real p;
    real fkv, PSL;
    real Stmp;

    p = rpow(10.0,y);
//    p = y;
//    fkv = Interpolation_nr(p, kfkT, fkT, nfkfunT, fkT2);
    PSL = Interpolation_nr(p, kpslT, pslT, npslfunT, pslT2);
//    fkv = fkInterpolation_nr(p, gd.kfk, gd.pfk, gd.nfk);
//    PSL = psInterpolation_nr(p, gd.kPS, gd.pPS, gd.nPSLT);

    Stmp = PSL;
    
//    printf("\nStmp: %g", Stmp);

    return p*Stmp;
}

local real sigma2Psi_function(void)
{
    real result;
    real kmin, kmax;
    real ymin, ymax;

    kmin = kpslT[1];
    kmax = kpslT[npslfunT];
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma2Psi_function_int,ymin,ymax,midpnt,EPSQ,KK);
//    result= (1.0/SIXPI2)
//    *qromo(sigma2v_function_int,kmin,kmax,midpnt,EPSQ,KK);

    return result;

}

#undef EPSQ
#undef KK

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;
    
    if ( k < kPS[1] || k > kPS[nPS] )
        fprintf(gd.outlog,"\n\nInterpolation_nr: warning! :: k is out of range (kPS[1], kPS[nPS])... %g %g %g\n",k, kPS[1], kPS[nPS]);

    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    
    return (psftmp);
}
/*
local  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS)
{
    real psftmp;

    splint(kPS,pPS,gd.pPS2,nPS,k,&psftmp);
    
    return (psftmp);
}

local  real fkInterpolation_nr(real k, double kfk[], double pfk[], int nfk)
{
    real fkftmp;

    splint(kfk,pfk,gd.pfk2,nfk,k,&fkftmp);
    
    return (fkftmp);
}
*/
// End: USER Model
// ==========================================

