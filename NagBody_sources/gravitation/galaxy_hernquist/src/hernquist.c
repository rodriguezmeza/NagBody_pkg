/*
 * MKHERNQUIST: 
 *
 * Generate a two-component, isotropic, spherical galaxy/halo model using a
 * superposition of two Hernquist models.  The galaxy is assumed to be a
 * unit Hernquist model with M=a=1 while the halo mass and scale radius are
 * free parameters specified on the command line.  The distribution
 * function is calculated on the fly using Eddington's formula and then
 * sampled to generate a N-body realization of a galaxy.  The galaxy and
 * halo particles are dumped in succession -- first half are galaxy
 * particles.  
 *
 *                  original version                John Dubinski - 1999
 *
 *      9-mar-99    Adopted for NEMO                Peter Teuben
 *                   - removed sorting and fcut stuff (snapsort,unbind)
 *                   - float->real (including the NR stuff; no query/ran1)
 *      9-sep-01    gsl/xrandom
 *     27-mar-03    fixed double/float usage or qromb (it never worked before)
 *	10-ene-12   for Gadget RGF
 */
//#include <time.h>
//#include <sys/types.h>
//#include <sys/times.h>
//#include <unistd.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_math.h>

//#include "../../../General_libs/general/stdinc.h"
#include "globaldefs.h"
#include "protodefs.h"


#define CONST (1.0/sqrt(32.)/(M_PI*M_PI*M_PI))
#define NTAB 256

 static gsl_rng *my_r = NULL;
 static const gsl_rng_type *my_T;
 static char env_type[] = "GSL_RNG_TYPE";    /* environment variables    */
 static char env_seed[] = "GSL_RNG_SEED";    /* used by GSL_RNG routines */

int galaxy;
local int nobj1, nobj2, nobj, ngas, ndisk, shiftblk;

int n, idum1;
double dE, mg, mh, mp, rmax1, rmax2;
double a1, a2, b, mu1, mu2, E, rcut1; 

/* local functions */

double rad(double eta), f(double E2), drho2_d(double u), pot(double r),
       radpsi(double p), rho1(double r), rho2(double r);
void   cofm(void), radmass(double);
float  drho2(float u);

// random generator defined below
double xrandom(double,double);
int init_xrandom(int);
int set_xrandom(int);


void create_hernquist(double mtot_type[], int ntot_type[], double ScaleTable[])
{

	int i, seed;
	double vx, vy, vz;
	double fmax, f0, f1, v2, vmax, vmax2;
        double lim1, lim2;
        float xcm, ycm, zcm;
        float vxcm, vycm, vzcm;
        float mtot;

// halo:
        mu1 = mtot_type[1];
        a1 =  ScaleTable[3];
        rmax1 = ScaleTable[5];
        nobj1 = ntot_type[1];
     if (rmax1 != 0.0)
        {
        lim1=rmax1*rmax1/pow(rmax1+a1,2.0);
	} else {
        lim1=1.0;
        }

// bulge:
        mu2 = mtot_type[3];
        a2 =  ScaleTable[9];
        rmax2 = ScaleTable[11];
        nobj2 = ntot_type[3];
     if (rmax2 != 0.0)
        {
        lim2=rmax2*rmax2/pow(rmax2+a2,2.0);
	} else {
        lim2=1.0;
        }
        
        ngas=ntot_type[0];
        ndisk=ntot_type[2];
        rcut1 = ScaleTable[8];
        nobj  = nobj1 + nobj2;

        mh = mu1/nobj1;
        mg = mu2/nobj2;

//        seed = init_xrandom(3276495886);
        seed = init_xrandom(cmd.seed);
	b = a1 - a2;

	xcm = ycm = zcm = vxcm = vycm = vzcm = mtot = 0;

	for(i=ngas; i<ngas+nobj; i++) {
		double eta, radius, cth, sth, phi;
		double psi0;


	if( i<nobj1+ngas ) {
                if( mu1 == 0.0 ) break;
		galaxy = 0;
		eta = (double) xrandom(0.0,lim1);
		radius = rad(eta);
		radius *= a1;
                shiftblk=0;
                mp=mh;
	} else {

                if( mu2 == 0.0 ) break;
		galaxy = 1;
                eta = (double) xrandom(0.0,lim2);
                radius = rad(eta);
		radius *= a2;
                shiftblk=ndisk;
                mp=mg;
        }
		phi = xrandom(0.0,2*M_PI);
		cth = xrandom(-1.0,1.0);
		sth = sqrt(1.0 - cth*cth);

		P[i+shiftblk].Pos[0] = (float) (radius*sth*cos(phi));
		P[i+shiftblk].Pos[1] = (float) (radius*sth*sin(phi));
		P[i+shiftblk].Pos[2] = (float) (radius*cth);

		psi0 = -pot(radius);
		vmax2 = 2.0*psi0;
		vmax = sqrt(vmax2);
		fmax = f(psi0);
		f0 = fmax; f1 = 1.1*fmax;        /* just some initial values */
		while( f1 > f0 ) {

                        /* pick a velocity */
			v2 = 2.0*vmax2;
			while( v2 > vmax2 ) {    /* rejection technique */
				vx = vmax*xrandom(-1.0,1.0);
				vy = vmax*xrandom(-1.0,1.0);
				vz = vmax*xrandom(-1.0,1.0);
				v2 = vx*vx + vy*vy + vz*vz;
			}

			f0 = f((psi0 - 0.5*v2));
			f1 = fmax*xrandom(0.0,1.0);

		}
		P[i+shiftblk].Vel[0] = vx;
		P[i+shiftblk].Vel[1] = vy;
		P[i+shiftblk].Vel[2] = vz;

                xcm += mp*P[i+shiftblk].Pos[0];
                ycm += mp*P[i+shiftblk].Pos[1];
                zcm += mp*P[i+shiftblk].Pos[2];
                vxcm += mp*vx;
                vycm += mp*vy;
                vzcm += mp*vz;
                mtot += mp;


        } 

        xcm /= mtot;
        ycm /= mtot;
        zcm /= mtot;
        vxcm /= mtot;
        vycm /= mtot;
        vzcm /= mtot;
        fprintf(stdout,"Centering Halo+Bulge (%g,%g,%g) (%g,%g,%g)\n",
                xcm,ycm,zcm,vxcm,vycm,vzcm);

       if (nobj1 > 0)
        {
        for(i=ngas; i<nobj1; i++) {
                P[i].Pos[0] -= xcm;
                P[i].Pos[1] -= ycm;
                P[i].Pos[2] -= zcm;
                P[i].Vel[0] -= vxcm;
                P[i].Vel[1] -= vycm;
                P[i].Vel[2] -= vzcm;
                }
        }

       if (nobj2 > 0)
        {
        for(i=ngas+nobj1+ndisk; i<ngas+nobj1+ndisk+nobj2; i++) {
                P[i].Pos[0] -= xcm;
                P[i].Pos[1] -= ycm;
                P[i].Pos[2] -= zcm;
                P[i].Vel[0] -= vxcm;
                P[i].Vel[1] -= vycm;
                P[i].Vel[2] -= vzcm;
                }
        }

    //    radmass(rcut1);

}
// END


	
double rad(double eta)
{
	double sqeta;
	sqeta = sqrt(eta);
	return(sqeta/(1.0-sqeta));
}


/* Distribution function of the 2-component Hernquist model */

double f(double E2)
{
	double ans;

	E = E2;
	ans = qromb(drho2,0.0,2.0*sqrt(E))*CONST;

	return ans;
}

float drho2(float u)
{
  return (float)drho2_d( (double) u);
}

double drho2_d(double u)
{
	double ans;
	double r;
	double c0, c1, c02, c03, c12, c13, c04, c14;
	double r2, r3;
	double p;
	double dp1, dp2, drho1r, drho2r;

	p = E - 0.25*u*u;
	if( p <= 0 )
		return 0;
	else
		r = radpsi(p);

	r2 = r*r; r3 = r2*r;

	c0 = a1 + r;  c1 = a2 + r;
	c02 = c0*c0; c12 = c1*c1;
	c03 = c02*c0; c13 = c12*c1; c04 = c02*c02; c14 = c12*c12;
	dp1 = -mu1/c02 - mu2/c12;
	dp2 = 2*mu1/c03 + 2*mu2/c13;

	if( !galaxy ) { /* its the halo */
		drho1r = - mu1*a1*(a1 + 4*r)/(r2*c04);
		drho2r = 2*mu1*a1*(10*r2 + 5*a1*r + a1*a1)/(r3*c04*c0);
	}
	else { /* its the galaxy */
		drho1r = - mu2*a2*(a2 + 4*r)/(r2*c14);
		drho2r = 2*mu2*a2*(10*r2 + 5*a2*r + a2*a2)/(r3*c14*c1);
	}

	ans = drho2r/(dp1*dp1) - drho1r/(dp1*dp1*dp1)*dp2;

	return ans;
}

double pot(double r)
{
	double p1, p2;

	p1 = -mu1/(a1 + r);
	p2 = -mu2/(a2 + r);

	return p1+p2;
}


double radpsi(double p)
{
	double b0, c0;
	double p1;
	double ans;

	p1 = 1/p;

//	b0 = 0.5*((1+a) - (1+mu)*p1);
//	c0 = a*(1-p1) - mu*p1;

	b0 = 0.5*((a1+a2) - (mu1+mu2)*p1);
	c0 = a2*(a1-mu1*p1) - mu2*a1*p1;

	ans = -b0 + sqrt(b0*b0 - c0);
	return ans;
}

void radmass(double rcut1)
{
        int i;
        double rcirs[NTAB];

    rcirs[0] = vcirs[0] = 0.0;
    for (i = 1; i < NTAB; i++) {
        rcirs[i] = rcut1 * ((double) i) / (NTAB - 1);
        vcirs[i] = mu1*rcirs[i]/pow(rcirs[i]+a1, 2.0) + mu2*rcirs[i]/pow(rcirs[i]+a2, 2.0);
    }


}
		
double rho1(double r)
{
	return mu1*a1/r/pow(a1+r,3.0);
}

double rho2(double r)
{
	return mu2*a2/r/pow(a2+r,3.0);
}


// produce a random number in between [xl,xh]
double xrandom(double xl, double xh)
{
    double retval;

     if (my_r == NULL) printf("GSL init_xrandom was never called\n");

    for(;;) {
            retval = gsl_rng_uniform(my_r);
            if (retval<0.0 || retval>1.0) {
            fprintf(stderr,"xrandom: spinning again, out of bounds [%g]",retval);
            continue;
            }
            break;
    }
    return xl + retval*(xh-xl);

}




// initialize the GSL RNG
int init_xrandom(int iseed)
{
    char my_gsl_type[64], my_gsl_seed[64];

          if (iseed > 0) {
          sprintf(my_gsl_seed,"%s=%u",env_seed,iseed);
          } else {
          iseed = set_xrandom(iseed);

          sprintf(my_gsl_seed,"%s=%u",env_seed,iseed); /* Produce a Seed */
          putenv(my_gsl_seed);

          sprintf(my_gsl_type,"%s=%s",env_type,"taus"); /* Type of RNG */
          putenv(my_gsl_type);



    gsl_rng_env_setup();                          /* initialize the rng (name/seed) setup */
    my_T = gsl_rng_default;
    my_r = gsl_rng_alloc(my_T);

    printf("GSL generator type: %s\n",gsl_rng_name(my_r));
    printf("GSL seed = %lu\n",gsl_rng_default_seed);
    printf("GSL first value = %lu\n",gsl_rng_get(my_r));

    return (int) gsl_rng_default_seed;

    }
    return set_xrandom(iseed);
}

//   generate seed  
int set_xrandom(int dum)
{
    int retval;
    struct tms buffer;

    if (dum <= 0) {
        if (dum == -1)
            retval = idum1 = (int) times(&buffer);   /* clock cycles */
        else if (dum == -2)
            retval = idum1 = (int) getpid();         /* process id */
        else            /* normally if dum==0 */
            retval = idum1 = (int) time(0);          /* seconds since 1970 */
    } else
        retval = idum1 = dum;               /* use supplied seed in argument */

    return retval;
}

