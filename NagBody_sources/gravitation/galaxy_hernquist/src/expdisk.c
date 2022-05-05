/*
 * MKEXPDISK:   set up a (bare) exponential disk. (based on mkbaredisk)
 *		has table output for debugging
 *		more mode creation methods
 *	original	created		Josh Barnes
 *	15-nov-90	helpvec		PJT
 *	24-mar-94	ansi fix
 *	24-mar-97	proto fixes	pjt
 *	29-mar-97	SINGLEPREC fixed, ndisk= now nbody=	pjt
 *      29-may-01       Add time      PJT
 *       8-sep-01       gsl/xrandom
 */

//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_bessel.h>

//#include "../../../General_libs/general/stdinc.h"
#include "globaldefs.h"
#include "protodefs.h"
#include "spline.h"
#include "xrandom.h"


/* DEFAULT INPUT PARAMETERS 
    "nbody=1024\n	  number of particles ",
    "alpha=4.0\n	  inverse exponential scale length ",
    "rcut=1.25\n	  outer cutoff radius ",
    "mdisk=1.0\n	  disk mass ",
    "Qtoomre=1.0\n	  Toomre's Q parameter ",
    "gamma=1.0\n	  fudge factor for v_eff(r) if > 0 ",
    "z0=0.025\n	          vertical scaleheight (softening) ",
    "mode=1\n             creation mode: 1=josh 2=kruit/searle ",
    "tab=f\n		  table output also? ",
    "zerocm=t\n           center the snapshot?",
    "VERSION=1.2d\n	  7-nov-05 PJT",
*/


double alpha, rcut, mdisk, Qtoomre, gammas, z0, mcut;
double ah, ab, mhalo, mbulge;

int  Qtab;
local int  ndisk, nstart;
int  cmode;

double gdisk(double);
double grandom(double,double);
void inittables(void);
void makedisk(void);
void centersnap(void);

int  nhalo, nstart_halo, nbulge, nstart_bulge;
void centersnap_wrspct_halo_and_bulge(void);


 FILE *ff;
 FILE *fh;


double sqr(double x)
{
    return x*x;
}

create_disk(double mtot_type[], int ntot_type[], double ScaleTable[])
{

    printf("set up a (bare) exponential disk\n");

    ah = ScaleTable[3];
    ab = ScaleTable[9];
    alpha = ScaleTable[6];
    z0 = ScaleTable[7];
    rcut = ScaleTable[8];
    ndisk = ntot_type[2];
    mhalo = mtot_type[1];
    mdisk = mtot_type[2];
    mbulge = mtot_type[3];
    Qtoomre = 2.0;
    gammas = 0.0;
    cmode = 2;
    Qtab = 1;
    nstart=ntot_type[0]+ntot_type[1];

    nstart_halo=ntot_type[0];
    nhalo=ntot_type[1];
    nstart_bulge=ntot_type[0]+ntot_type[1]+ntot_type[2];
    nbulge=ntot_type[3];

     ff = fopen("qtab.dat", "w");
     fh = fopen("RC.dat", "w");


    inittables();
    makedisk();
} // END

#define NTAB 256

double rcir[NTAB];
double vcir[4*NTAB];
double vcirs[4*NTAB];       // definido en globaldefs.h para transferirlo a hernquist.c
double mdsk[NTAB];
double rdsk[4*NTAB];




void inittables()
{

    int i;

    rcir[0] = vcir[0] = vcirs[0] = 0.0;
    for (i = 1; i < NTAB; i++) {
	rcir[i] = rcut * ((double) i) / (NTAB - 1);
        vcirs[i] = mhalo*rcir[i]/pow(rcir[i]+ah, 2.0) + mbulge*rcir[i]/pow(rcir[i]+ab, 2.0);
	vcir[i] = sqrt(- gdisk(rcir[i]) * rcir[i] + vcirs[i]);
    fprintf (fh,"%g %g %g %g \n", rcir[i], sqrt(vcirs[i]), sqrt(-gdisk(rcir[i])*rcir[i]), sqrt(vcir[i]));
    }

    spline(&vcir[NTAB], &rcir[0], &vcir[0], NTAB);

    for (i = 0; i < NTAB; i++) {
	rdsk[i] = rcut * pow(((double) i) / (NTAB - 1), 1.5);
	mdsk[i] = 1 - (1 + alpha * rdsk[i]) * exp(- alpha * rdsk[i]);
    }

    spline(&rdsk[NTAB], &mdsk[0], &rdsk[0], NTAB);
}

double gdisk(double rad)
{
    double x;

    x = 0.5 * alpha * rad;
//    return - alpha*alpha * mdisk * x *
//	      (bessi0(x) * bessk0(x) - bessi1(x) * bessk1(x));
    return - alpha*alpha * mdisk * x *
       (gsl_sf_bessel_I0(x) * gsl_sf_bessel_K0(x) - gsl_sf_bessel_I1(x) * gsl_sf_bessel_K1(x));
}

void makedisk()
{
    int i, k, nzero=0;
    double mdsk_i, rad_i, theta_i, vcir_i, omega, Aoort, kappa;
    double mu, sig_r, sig_t, sig_z, vrad_i, veff_i, vorb_i;
    k=0;

   if (Qtab)
     fprintf (ff, "# rad_i,mdsk_i,vcir_i,omega,kappa,Aoort,mu,sig_r,sig_t,sig_z,veff_i, Qtoomre, sig_t/sig_r,sig_z/sig_r, 1.5-(sqr(sig_t) + 0.5*sqr(sig_z))/sqr(sig_r)\n");

    for (i = nstart; i < ndisk + nstart; i++) {

	mcut = mdsk[NTAB-1] / ndisk;
	mdsk_i = mdsk[NTAB-1] * ((double) k + 1.0) / ndisk;
	rad_i = seval(mdsk_i, &mdsk[0], &rdsk[0], &rdsk[NTAB], NTAB);
	theta_i = xrandom(0.0, 2*M_PI);

	P[i].Pos[0] = (float) (rad_i * sin(theta_i));	/* assign positions */
	P[i].Pos[1] = (float) (rad_i * cos(theta_i));
	P[i].Pos[2] = (float) grandom(0.0, 0.5 * z0);

	vcir_i = seval(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB);
	omega = vcir_i / rad_i;

	Aoort = - 0.5 *
	    (spldif(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB) - omega);
	if (omega - Aoort < 0.0)
	    printf("rad_i, omega, Aoort = %f %f %f\n", rad_i, omega, Aoort);

	kappa = 2 * sqrt(omega*omega - Aoort * omega);
	mu = alpha*alpha * mdisk * exp(- alpha * rad_i) / (2*M_PI);

	if (cmode==1) {                 /* Straight from Josh - mkbaredisk*/
	   sig_r = 3.358 * Qtoomre * mu / kappa;
	   sig_t = 0.5 * sig_r * kappa / omega;
	   sig_z = 0.5 * sig_r;
	} else if (cmode==2) {
	   sig_z = sqrt(M_PI * mu * z0);          /* isothermal sech sheet */
           sig_r = 2.0 * sig_z;                 /* with constant scaleheight */
           Qtoomre = sig_r * kappa / (3.358 * mu);  /* See vdKruit/Searle */
	   sig_t = 0.5 * sig_r * kappa / omega;
        } else {
	    printf("illegal mode=%d",cmode);
            exit(0);
            }

	vrad_i = grandom(0.0, sig_r);

// Asymmetric drift
	if (gammas > 0.0) 		/* Josh' method: averaged */
	   veff_i = (vcir_i*vcir_i +
			(gammas - 3*alpha*rad_i) * sig_r*sig_r);
	else				/* a little more accurate */
	   veff_i = sqr(vcir_i) - sqr(sig_r) * 
	      (sqr(sig_t/sig_r) - 1.5 + 0.5*sqr(sig_z/sig_r) + 2*alpha*rad_i);
	if (veff_i < 0.0) {
            nzero++;
            veff_i = 0.0;
        } else
            veff_i = sqrt(veff_i);

	vorb_i = veff_i + grandom(0.0, sig_t);

	P[i].Vel[0] = (float)			/* assign velocities */
	  ((vrad_i * P[i].Pos[0] + vorb_i * P[i].Pos[1]) / rad_i);
	P[i].Vel[1] = (float)
	  ((vrad_i * P[i].Pos[1] - vorb_i * P[i].Pos[0]) / rad_i);
	P[i].Vel[2] = (float) grandom(0.0, sig_z);

	if (Qtab) {
            fprintf (ff,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
            rad_i,mdsk_i,vcir_i,omega,kappa,Aoort,mu,sig_r,sig_t,sig_z,veff_i,
            Qtoomre,
            sig_t/sig_r,sig_z/sig_r,
            1.5-(sqr(sig_t) + 0.5*sqr(sig_z))/sqr(sig_r) );
        }
        k++;
    }

    if (nzero)
        printf("Warning: %d stars with too little orbital motion\n",nzero);

    centersnap();
//    centersnap_wrspct_halo();
}


void centersnap()
{
    
    int i;
    float xcm, ycm, zcm;
    float vxcm, vycm, vzcm;
    float mtot;
    
    xcm = ycm = zcm = vxcm = vycm = vzcm = 0.0;
    mtot = 0.0;
    
    for(i=nstart; i<ndisk+nstart; i++) {
        xcm += mcut*P[i].Pos[0];
        ycm += mcut*P[i].Pos[1];
        zcm += mcut*P[i].Pos[2];
        vxcm += mcut*P[i].Vel[0];
        vycm += mcut*P[i].Vel[1];
        vzcm += mcut*P[i].Vel[2];
        mtot += mcut;
    }
    xcm /= mtot;
    ycm /= mtot;
    zcm /= mtot;
    vxcm /= mtot;
    vycm /= mtot;
    vzcm /= mtot;
    fprintf(stdout,"Centering snapshot (%g,%g,%g) (%g,%g,%g)\n",
            xcm,ycm,zcm,vxcm,vycm,vzcm);
    
    
    for(i=nstart; i<ndisk+nstart; i++) {
        P[i].Pos[0] -= xcm;
        P[i].Pos[1] -= ycm;
        P[i].Pos[2] -= zcm;
        P[i].Vel[0] -= vxcm;
        P[i].Vel[1] -= vycm;
        P[i].Vel[2] -= vzcm;
    }
}

void centersnap_wrspct_halo()
{
    
    int i;
    float xcm, ycm, zcm;
    float vxcm, vycm, vzcm;
    float mtot;
    
    xcm = ycm = zcm = vxcm = vycm = vzcm = 0.0;
    mtot = 0.0;
    
    for(i=nstart_halo; i<nhalo+nstart_halo; i++) {
        xcm += mcut*P[i].Pos[0];
        ycm += mcut*P[i].Pos[1];
        zcm += mcut*P[i].Pos[2];
        vxcm += mcut*P[i].Vel[0];
        vycm += mcut*P[i].Vel[1];
        vzcm += mcut*P[i].Vel[2];
        mtot += mcut;
    }
    xcm /= mtot;
    ycm /= mtot;
    zcm /= mtot;
    vxcm /= mtot;
    vycm /= mtot;
    vzcm /= mtot;
    fprintf(stdout,"Centering bulge and disk with respect to halo (%g,%g,%g) (%g,%g,%g)\n",
            xcm,ycm,zcm,vxcm,vycm,vzcm);
    
    for(i=nstart_bulge; i<nbulge+nstart_bulge; i++) {
        P[i].Pos[0] -= xcm;
        P[i].Pos[1] -= ycm;
        P[i].Pos[2] -= zcm;
        P[i].Vel[0] -= vxcm;
        P[i].Vel[1] -= vycm;
        P[i].Vel[2] -= vzcm;
    }

    for(i=nstart; i<ndisk+nstart; i++) {
        P[i].Pos[0] -= xcm;
        P[i].Pos[1] -= ycm;
        P[i].Pos[2] -= zcm;
        P[i].Vel[0] -= vxcm;
        P[i].Vel[1] -= vycm;
        P[i].Vel[2] -= vzcm;
    }
}


/* Gaussian RND generator */
double grandom(double mean, double sdev)
{
    static double v1, v2, s;
    static int gcount = 0;

    if (gcount) {
        gcount = 0;
        return mean + sdev * v2 * sqrt(-2.0 * log(s) / s);
    }

    do {                                /* loop until */
        v1 = xrandom(-1.0, 1.0);        /* two points within */
        v2 = xrandom(-1.0, 1.0);        /* the unit square */
        s = v1*v1 + v2*v2;              /* have radius such that */
    } while (s >= 1.0);                 /* point fall inside unit circle */

    gcount = 1;
    return mean + sdev * v1 * sqrt(-2.0 * log(s) / s);

}

