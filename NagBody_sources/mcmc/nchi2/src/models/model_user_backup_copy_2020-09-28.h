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
} model_data_user, *model_data_user_ptr;

local model_data_user mduser;

#define md          mduser

// params[1] === rhos   :: range: 0 < p1 < 10000
// params[2] === rs     :: range: 0 < p2 < 10

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

// End: USER Model
// ==========================================

