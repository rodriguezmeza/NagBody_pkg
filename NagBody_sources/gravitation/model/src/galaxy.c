/*==============================================================================
	MODULE: galaxy.c			[model]
	Written by: Mario A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Galaxy model
	Language: C
	Use:
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: July 23, 2007
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#define global
#include "globaldefs.h"

local real ab,sb,xr;


local real ad,xr,z0,xz,sr,v0;


local real Rmaxb, Rmaxd, Rmaxh, Rmaxg, Rmaxdm;

#define outputfile "galaxy.dat"

local void bulge_distribution(int, int, int, vector, vector, real);
static float fx_bulge_v(float);
local void disk_distribution(int, int, int, vector, vector, real);

static float fx_disk_r(float);
static float fx_disk_z(float);
static float fx_disk_vr(float);
static float fx_disk_vphi(float);
local void halo_distribution(int, int, int, vector, vector, real);
local real m2_rhob(real, real);
local real m2_rhoh(real, real);
local real m2_rhod(real, real, real);
local void galaxy_phasesp_out(void);
local void BDH_avg_disk(int, int);

#define GALAXY_1 "galaxy1"
#define GALAXY_2 "galaxy2"
#define GALAXY_DAT "galaxy.dat"
#define BDHGGALAXY_DAT "bdhggal.dat"
#define BDHGDMGALAXY_DAT "bdhgdmg.dat"
#define BDHGALAXY_IN "bdhgal.in"
#define BDHGGALAXY_IN "bdhggal.in"
#define BDHGDMGALAXY_IN "bdhgdmg.in"
#define TWOBDHGALAXY_PAR_IN "2bdhpar.in"


void Two_BDH_Galaxies_Coll(void)
{
	int nbody1,nbody2;
	int nbod1,nbod2;
	int ndim1,ndim2;
	real tnow1,tnow2,tmpx,tmpy,tmpz;
	vector Rcm1,Rcm2,Vcm1,Vcm2;
	stream in_galaxy_1, in_galaxy_2, galaxy_dat, galaxy_in;
    struct stat buf;
	bodyptr p,q;

	real alpha1, beta1, gamma1;
	real alpha2, beta2, gamma2;
	real Rmaxb1, Rmaxd1, Rmaxh1;
	real Rmaxb2, Rmaxd2, Rmaxh2;

	strcpy(gd.model_comment, "Two BDH Galaxies Collision Model");

    if (stat(GALAXY_1, &buf) != 0)               
        error("The file galaxy1 does not exist!");         
    else 
	    in_galaxy_1 = stropen(GALAXY_1, "r");

    if (stat(GALAXY_2, &buf) != 0)               
        error("The file galaxy2 does not exist!");         
    else 
	    in_galaxy_2 = stropen(GALAXY_2, "r");

    if (stat(GALAXY_DAT, &buf) != 0)               
        error("The file galaxy.dat does not exist!");         
    else 
        galaxy_dat = stropen(GALAXY_DAT, "r");         

    if (stat(TWOBDHGALAXY_PAR_IN, &buf) != 0)               
        error("The file 2bdhpar.in does not exist!");         
    else 
        galaxy_in = stropen(TWOBDHGALAXY_PAR_IN, "r");         

	fscanf(galaxy_dat,"%d %d %d %d %d %d %d %lf %lf %lf",
		&Nbi1,&Nbf1,&Ndi1,&Ndf1,&Nhi1,&Nhf1,&nbod1,&Rmaxb1,&Rmaxd1,&Rmaxh1);
	printf("Galaxy (1) data: %d %d %d %d %d %d %d %lf %lf %lf\n",
		Nbi1,Nbf1,Ndi1,Ndf1,Nhi1,Nhf1,nbod1,Rmaxb1,Rmaxd1,Rmaxh1);
	fscanf(galaxy_dat,"%d %d %d %d %d %d %d %lf %lf %lf",
		&Nbi2,&Nbf2,&Ndi2,&Ndf2,&Nhi2,&Nhf2,&nbod2,&Rmaxb2,&Rmaxd2,&Rmaxh2);
	printf("Galaxy (2) data: %d %d %d %d %d %d %d %lf %lf %lf\n",
		Nbi2,Nbf2,Ndi2,Ndf2,Nhi2,Nhf2,nbod2,Rmaxb2,Rmaxd2,Rmaxh2);
	fclose(galaxy_dat);

	CLRV(Rcm1);
	CLRV(Vcm1);
	CLRV(Rcm2);
	CLRV(Vcm2);
	printf("Rcm and Vcm of galaxy (1) is set to zero.\n");
	fscanf(galaxy_in,"%lf%lf%lf", &alpha1, &beta1, &gamma1);
    in_vector(galaxy_in, Rcm1);
    in_vector(galaxy_in, Vcm1);
	fscanf(galaxy_in,"%lf%lf%lf", &alpha2, &beta2, &gamma2);
    in_vector(galaxy_in, Rcm2);
    in_vector(galaxy_in, Vcm2);
	printf("galaxy (1):\n");
	printf("%f %f %f\n",alpha1,beta1,gamma1);
	printf("%f %f %f\n",Rcm1[0],Rcm1[1],Rcm1[2]);
	printf("%f %f %f\n",Vcm1[0],Vcm1[1],Vcm1[2]);
	printf("galaxy (2):\n");
	printf("%f %f %f\n",alpha2,beta2,gamma2);
	printf("%f %f %f\n",Rcm2[0],Rcm2[1],Rcm2[2]);
	printf("%f %f %f\n",Vcm2[0],Vcm2[1],Vcm2[2]);

	fscanf(in_galaxy_1,"%d",&nbody1);
	fscanf(in_galaxy_1,"%d",&ndim1);
	fscanf(in_galaxy_1,"%f",&tnow1);
	fscanf(in_galaxy_2,"%d",&nbody2);
	fscanf(in_galaxy_2,"%d",&ndim2);
	fscanf(in_galaxy_2,"%f",&tnow2);

	if (NDIM != ndim1 || NDIM != ndim2) 
		error("Dimensions are different!\n");
	if (nbod1 != nbody1) 
		error("Number of bodies of galaxy (1) does not match!\n");
	if (nbod2 != nbody2) 
		error("Number of bodies of galaxy (2) does not match!\n");

	printf("First data: %d %d %d %d %d %d\n",Nbi2,Nbf2,Ndi2,Ndf2,Nhi2,Nhf2);
	Nbi2 += nbody1;
	Nbf2 += nbody1;
	Ndi2 += nbody1;
	Ndf2 += nbody1;
	Nhi2 += nbody1;
	Nhf2 += nbody1;
	printf("Second data: %d %d %d %d %d %d\n",Nbi2,Nbf2,Ndi2,Ndf2,Nhi2,Nhf2);
    if (stat(GALAXY_DAT, &buf) != 0)               
        error("The file galaxy.dat does not exist!");         
    else 
        galaxy_dat = stropen(GALAXY_DAT, "a");         
	fprintf(galaxy_dat,"\n%d %d %d %d %d %d %d",
		Nbi2,Nbf2,Ndi2,Ndf2,Nhi2,Nhf2,nbody2);
	fclose(galaxy_dat);

    error("Creating bdh galaxies finished\n");


	cmd.nbody = nbody1 + nbody2;
	gd.tnow = 0.0;
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
    for (p = bodytab; p < bodytab+nbody1; p++)   
        in_real(in_galaxy_1, &Mass(p));               
    for (p = bodytab+nbody1; p < bodytab+cmd.nbody; p++)   
        in_real(in_galaxy_2, &Mass(p));

	printf("nbody1,nbody2: %d %d\n",nbody1,nbody2);
	printf("Rcm1: %14.7f %14.7f %14.7f\n",Rcm1[0],Rcm1[1],Rcm1[2]);
	printf("Vcm1: %14.7f %14.7f %14.7f\n",Vcm1[0],Vcm1[1],Vcm1[2]);
	printf("Rcm2: %14.7f %14.7f %14.7f\n",Rcm2[0],Rcm2[1],Rcm2[2]);
	printf("Vcm2: %14.7f %14.7f %14.7f\n",Vcm2[0],Vcm2[1],Vcm2[2]);

	for (p = bodytab; p < bodytab+nbody1; p++) {
        in_vector(in_galaxy_1, Pos(p));
		ADDV(Pos(p),Pos(p),Rcm1);
	}
    for (p = bodytab+nbody1; p < bodytab+cmd.nbody; p++) {
        in_vector(in_galaxy_2, Pos(p));
		ADDV(Pos(p),Pos(p),Rcm2);
	}
    for (p = bodytab; p < bodytab+nbody1; p++) {
        in_vector(in_galaxy_1,Vel(p));
		ADDV(Vel(p),Vel(p),Vcm1);
	}
    for (p = bodytab+nbody1; p < bodytab+cmd.nbody; p++) {
        in_vector(in_galaxy_2, Vel(p));
		ADDV(Vel(p),Vel(p),Vcm2);
	}
    for (p = bodytab; p < bodytab+cmd.nbody; p++)   
        Type(p) = BODY;                         
    for (p = bodytab; p < bodytab+cmd.nbody-1; p++) {
		for (q = p+1; q < bodytab+cmd.nbody; q++) {
			if ( (Pos(p)[0] == Pos(q)[0]) && (Pos(p)[1] == Pos(q)[1]) &&
				(Pos(p)[2] == Pos(q)[2]) ) {
				printf("Pos(%d)= %14.7f %14.7f %14.7f\n",
					p-bodytab,Pos(p)[0],Pos(p)[1],Pos(p)[2]);
				printf("Pos(%d)= %14.7f %14.7f %14.7f\n",
					q-bodytab,Pos(q)[0],Pos(q)[1],Pos(q)[2]);
			}
		}
	}
	printf("finished two galaxy coll initial conditions\n");
	fclose(in_galaxy_1);
	fclose(in_galaxy_2);
}

void BDH_Galaxy(void)
{
	vector Rcmb,Rcmd,Rcmh,Vcmb,Vcmd,Vcmh;
	real Mb,Md,Mh;
	int np;
	stream galaxy_dat, galaxy_in;
    struct stat buf;

	strcpy(gd.model_comment, "BDH Galaxy Model");

    if (stat(BDHGALAXY_IN, &buf) != 0)               
        error("The file bdhgal.in does not exist!");         
    else 
	    galaxy_in = stropen(BDHGALAXY_IN, "r");

	fscanf(galaxy_in,"%d%d%d", &Nb, &Nd, &Nh);

	cmd.nbody=Nb+Nd+Nh;
	printf("Nb,Nd,Nh,nbody= %d %d %d %d\n",Nb,Nd,Nh,cmd.nbody);
	Nbi = 0;
	Nbf = Nb-1;
	Ndi = Nb;
	Ndf = Nb+Nd-1;
	Nhi = Nb+Nd;
	Nhf = Nb+Nd+Nh-1;
    if (stat(GALAXY_DAT, &buf) != 0)               
        galaxy_dat = stropen(GALAXY_DAT, "w");         
    else 
        galaxy_dat = stropen(GALAXY_DAT, "a");         
	fprintf(galaxy_dat,"\n%d %d %d %d %d %d %d",
		Nbi,Nbf,Ndi,Ndf,Nhi,Nhf,cmd.nbody);
	CLRV(Rcmb);
	CLRV(Rcmd);
	CLRV(Rcmh);
	CLRV(Vcmb);
	CLRV(Vcmd);
	CLRV(Vcmh);
	Mb = 0.0625;
	Md = 0.1875;
	Mh = 1.0;
	Rmaxb = 1.0;
	Rmaxd = 2.0;
	Rmaxh = 5.0;

	fprintf(galaxy_dat,"%14.7f %14.7f %14.7f",Rmaxb,Rmaxd,Rmaxh);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
	printf("bulge\n");
	bulge_distribution(Nbi,Nbf,Nb,Rcmb,Vcmb,Mb);
	printf("disk\n");
	disk_distribution(Ndi,Ndf,Nd,Rcmd,Vcmd,Md);
	printf("halo\n");
	halo_distribution(Nhi,Nhf,Nh,Rcmh,Vcmh,Mh);
	printf("Nbi,Nbf,Ndi,Ndf,Nhi,Nhf: %d %d %d %d %d %d\n",
		Nbi,Nbf,Ndi,Ndf,Nhi,Nhf);
	printf("finished initial conditions\n");
	fclose(galaxy_dat);
}


#define NZ 50
#define NBMAX 20
#define X1 0.0
#define X2 10.0
#define AB	0.04168

void bulge_distribution(int Nbi, int Nbf, int Nb, 
	vector Rcmb, vector Vcmb, real Mb)
{
	real fsen,mass_i,fv,rob,frb,fb;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real rtheta,rphi,fr;
	bodyptr p;
	float tol,xb1[NBMAX],xb2[NBMAX];
	int nbr=NBMAX;

	fb = rsqr( AB/Rmaxb + ((real) 1.0) );
	mass_i = Mb / ((real) Nb);
	for (p = bodytab+Nbi; p < bodytab+Nbf+1; p++) {
        Type(p) = BODY;
		rtheta	= racos( ((real) 1.0) - ((real) 2.0) * xrandom(0.0,1.0) );
		rphi	= ((real) 2.0) * PI * xrandom(0.0,1.0);
		fr = ((real) 1.5) * Rmaxb;
		while (fr > Rmaxb) {
			xr = xrandom(0.0,1.0);
			while( (xr >= 1.0) || (xr <= 0.0) ) {
				xr = xrandom(0.0,1.0);
			}
			fr	= ( xr/fb + rsqrt(xr/fb) ) * AB / ( ((real) 1.0) - xr/fb );
		}
		Mass(p)	= mass_i;
		fsen	= fr*rsin(rtheta);
		Pos(p)[0]	= fsen*rcos(rphi);
		Pos(p)[1]	= fsen*rsin(rphi);
		Pos(p)[2]	= fr*rcos(rtheta);

		rtheta	= racos( ((real) 1.0) - ((real) 2.0) * xrandom(0.0,1.0) );
		rphi	= ((real) 2.0) * PI * xrandom(0.0,1.0);
		xr = xrandom(0.0,1.0);
		frb=fr+AB;
		rob=fr/AB;
		if ( (fr/AB) >= ((real) 2000.0) ) {
			sb = ((real) 1.0) / (((real) 5.0) * fr);
		} else {
			sb = ( ((real) 12.0)*fr*rpow(frb,3)/rpow(AB,4) )*rlog((frb)/fr)
				-(fr/frb)*( ((real) 25.0)+ ((real) 52.0)*rob+				
				((real) 42.0)*rsqr(rob)+ ((real) 12.0)*rpow(rob,3) );
		}
		if (sb <= 0.0) 
				error("bulge: sb < 0; fr/ab=, fr/ab\n");
		sb = rsqrt( ( ((real) 1.0) / ( ((real) 12.0) * AB ) ) * sb );
		zbrak(fx_bulge_v,X1,X2,NZ,xb1,xb2,&nbr);
		if (nbr == 0) nrerror("bulge_v: no roots founds");
		tol=(1.0e-6)*(xb1[1]+xb2[1])/2.0;
		fv  =zbrent(fx_bulge_v,xb1[1],xb2[1],tol);
		fsen	= fv*rsin(rtheta);
		Vel(p)[0]	= fsen*rcos(rphi);
		Vel(p)[1]	= fsen*rsin(rphi);
		Vel(p)[2]	= fv*rcos(rtheta);
	}
	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for (p = bodytab+Nbi; p < bodytab+Nbf+1; p++){
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0]*Mass(p);
		sum_1y += Pos(p)[1]*Mass(p);
		sum_1z += Pos(p)[2]*Mass(p);
	}
	sum_0  = sum_0/Mb;
	sum_1x = sum_1x/Mb;
	sum_1y = sum_1y/Mb;
	sum_1z = sum_1z/Mb;
	printf("Bulge distribution\n");
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for (p = bodytab+Nbi; p < bodytab+Nbf+1; p++){
		sum_2x += rsqr(Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += rsqr(Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += rsqr(Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x/Mb;
	sum_2y = sum_2y/Mb;
	sum_2z = sum_2z/Mb;
	printf("2do. momento: %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z,m2_rhob(AB,Rmaxb));
	for (p = bodytab+Nbi; p < bodytab+Nbf; p++){
		ADDV(Pos(p),Pos(p),Rcmb);
		ADDV(Vel(p),Vel(p),Vcmb);
	}
	printf("Bulge finished!\n");
}

static float fx_bulge_v(float v)
{
	real arg;
	float out;
	arg = v/(SQRT2 * sb);
	out = merff(arg) - ( ((real) 2.0)/SQRTPI)*arg*rexp(-arg*arg) - xr;
	return(out);
}

real pr_b(real r, real dr, real Rmax, real ab)
{
	real Nb;
	Nb = ((real) 2.0) * ab * rsqr(ab/Rmax + ((real) 1.0)) * r * dr;
	return( Nb / rpow(r + ab,3.0) );
}

local real m2_rhob(real a, real rmax)
{
	real out,f1;

	f1 = rsqr(a/rmax + ((real) 1.0));
	out = f1*a*(rmax*(6.0*rpow(a,2.0) + 9.0*a*rmax + 2.0*rpow(rmax,2.0)) + 
			6.0*a*rpow(a + rmax,2.0)*rlog(a) - 
			6.0*a*rpow(a + rmax,2.0)*rlog(a + rmax))/
			(rpow(a + rmax,2.0));
	return(out);
}

#undef NZ
#undef NBMAX
#undef X1
#undef X2


#define NZ 5
#define NBMAX 50
#define X1 0.0
#define X2 10.0

#define NZY 5
#define NBYMAX 50
#define Y1 -100.0
#define Y2 100.0

#define NZXX 100
#define NBXMAX 50
#define XX1 -10.0
#define XX2 10.0

#define NZZ 5
#define NBZMAX 50
#define Z1 -5.0
#define Z2  5.0

#define NZW 5
#define NBWMAX 100
#define W1 -2.0
#define W2  5.0

#define AD	0.083333333333333333333333333333333
#define Z0	0.007

local void disk_distribution(int Ndi, int Ndf, int Nd, 
	vector Rcmd, vector Vcmd, real Md)
{
	real Ort_A,Ort_B,disk_Q,kappa,sphi,sigma,sz;
	real mass_i,vphi,vz,vr;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real rtheta,r,z;
	bodyptr p;
	real fx1,fx2;
	float tol;
	float xb1[NBMAX],xb2[NBMAX];
	float yb1[NBYMAX],yb2[NBYMAX];
	float xxb1[NBXMAX],xxb2[NBXMAX];
	float zb1[NBZMAX],zb2[NBZMAX];
	float wb1[NBWMAX],wb2[NBWMAX];
	int nb=NBMAX;
	int nby=NBYMAX;
	int nbxx=NBXMAX;
	int nbz=NBZMAX;
	int nbw=NBWMAX;

	mass_i = Md / ((real) Nd);
	for ( p = bodytab+Ndi; p < bodytab+Ndf+1; p++) {
        Type(p) = BODY;
		rtheta	= ((real) 2.0) * PI_D * xrandom(0.0,1.0);
		r = ((real) 1.5) * Rmaxd;
		while (r > Rmaxd) {
			xr = xrandom(0.0,1.0);
			zbrak(fx_disk_r,X1,X2,NZ,xb1,xb2,&nb);
			if (nb == 0) nrerror("disk_r: no roots founds");
			tol=(1.0e-6)*(rabs(xb1[1]+xb2[1]))/2.0;
			r=zbrent(fx_disk_r,xb1[1],xb2[1],tol);
		}

		xz = xrandom(0.0,1.0);
		zbrak(fx_disk_z,Y1,Y2,NZY,yb1,yb2,&nby);
		if (nby == 0) nrerror("disk_z: no roots founds");
		tol=(1.0e-6)*(rabs(yb1[1]+yb2[1]))/2.0;
		z=zbrent(fx_disk_z,yb1[1],yb2[1],tol);

		Mass(p)	  = mass_i;
		Pos(p)[0] = r*rcos(rtheta);
		Pos(p)[1] = r*rsin(rtheta);
		Pos(p)[2] = z;

		v0 = 1.0;
		sigma = Md * rexp(-r/AD)/(TWOPI_D * rsqr(AD));
		sz = rsqrt(PI_D * sigma * Z0);
		Ort_A = 14.5;
		Ort_B = -12.0;
		sr = ((real) 2.0) * sz;
		sphi = rsqrt(-Ort_B/(Ort_A-Ort_B)) * sr;
		kappa = 1.0;
		disk_Q = sr * kappa / ( ((real) 3.36) * sigma ); 

		xr=xrandom(0.0,1.0);
		zbrak(fx_disk_vr,XX1,XX2,NZXX,xxb1,xxb2,&nbxx);
		if (nbxx == 0) {
			printf("xr,sr: %14.7f %14.7f\n",xr,sr);
			printf("nbxx,nzxx: %14.7f %14.7f\n",nbxx,NZXX);
			printf("xx1,xx2: %14.7f %14.7f\n",XX1,XX2);
			fx1 = fx_disk_vr(XX1);
			fx2 = fx_disk_vr(XX2);
			printf("fx1,fx2: %14.7f %14.7f\n",fx1,fx2);
			nrerror("disk_vr: no roots founds");
		}
		tol=(1.0e-6)*(rabs(xxb1[1]+xxb2[1]))/2.0;
		vr=zbrent(fx_disk_vr,xxb1[1],xxb2[1],tol);

		xr=xrandom(0.0,1.0);
		sr=sz;
		zbrak(fx_disk_vr,Z1,Z2,NZZ,zb1,zb2,&nbz);
		if (nbz == 0) nrerror("disk_vz: no roots founds");
		tol=(1.0e-6)*(rabs(zb1[1]+zb2[1]))/2.0;
		vz=zbrent(fx_disk_vr,zb1[1],zb2[1],tol);

		xr=xrandom(0.0,1.0);
		sr=sphi;
		zbrak(fx_disk_vphi,W1,W2,NZW,wb1,wb2,&nbw);
		if (nbw == 0) nrerror("disk_vphi: no roots founds");
		tol=(1.0e-6)*(rabs(wb1[1]+wb2[1]))/2.0;
		vphi = zbrent(fx_disk_vphi,wb1[1],wb2[1],tol);

		Vel(p)[0] = vr*rcos(rtheta) - vphi*rsin(rtheta);
		Vel(p)[1] = vr*rsin(rtheta) + vphi*rcos(rtheta);
		Vel(p)[2] = vz;
	} 
	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for ( p = bodytab+Ndi; p < bodytab+Ndf+1; p++) {
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0] * Mass(p);
		sum_1y += Pos(p)[1] * Mass(p);
		sum_1z += Pos(p)[2] * Mass(p);
	}
	sum_0=sum_0/Md;
	sum_1x = sum_1x/Md;
	sum_1y = sum_1y/Md;
	sum_1z = sum_1z/Md;
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for ( p = bodytab+Ndi; p < bodytab+Ndf+1; p++) {
		sum_2x += (Pos(p)[0]-sum_1x) * (Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += (Pos(p)[1]-sum_1y) * (Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += (Pos(p)[2]-sum_1z) * (Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x/Md;
	sum_2y = sum_2y/Md;
	sum_2z = sum_2z/Md;
	printf("2do. momento: %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z,m2_rhod(AD,Z0,Rmaxd));
	for ( p = bodytab+Ndi; p < bodytab+Ndf+1; p++) {
		ADDV(Pos(p),Pos(p),Rcmd);
		ADDV(Vel(p),Vel(p),Vcmd);
	}
	printf("disk finished!\n");
}

static float fx_disk_r(float r)
{
	float out,fd;
	fd = ((real) 1.0) / 
		( ((real) 1.0) - rexp(-Rmaxd/AD) * ( ((real) 1.0) + Rmaxd/AD ) );
	out = rexp(-r/AD)*( ((real) 1.0) + r/AD) + xr/fd - ((real) 1.0);
	return(out);
}

static float fx_disk_z(float z)
{
	float out;
	out = rtanh(z/Z0) + ((real) 1.0) - ((real) 2.0) * xz;
	return(out);
}

static float fx_disk_vr(float vr_1)
{
	float out;
	out = merff(vr_1/(SQRT2*sr)) + ((real) 1.0) - ((real) 2.0) * xr;
	return(out);
}



static float fx_disk_vphi(float v)
{
	float out;
	out = merff( (v-v0)/(SQRT2*sr)) + ((real) 1.0) - ((real) 2.0) * xr;
	return(out);
}

real pr_d(real r, real dr, real Rmax, real ad)
{
	real Nd, fd;

	fd = 1 / ( 1 - rexp(-Rmax/ad) * ( 1 + Rmax/ad ) );
	Nd = fd * r * dr / rsqr(ad);
	return( Nd * rexp(-r/ad) );
}

local real m2_rhod(real a, real z0, real rmax)
{
	real out, fd, era;
	era = rexp(-rmax/a);
	fd = 1 / ( 1 - era * ( 1 + rmax/a ) );
	out = fd * ( 6*rsqr(a) - era * ( 6*a + 6*rmax + 3 * rpow(rmax,2)/a
			+ rpow(rmax,3)/rsqr(a) ) + rsqr(PI_D * z0) * fd / 12 );
	return(out);
}
#undef NZ
#undef NBMAX
#undef X1
#undef X2

#undef NZY
#undef NBYMAX
#undef Y1
#undef Y2

#undef NZXX
#undef NBXMAX
#undef XX1
#undef XX2

#undef NZZ
#undef NBZMAX
#undef Z1
#undef Z2

#undef NZW
#undef NBWMAX
#undef W1
#undef W2

#undef AD
#undef Z0

#define NZ 5
#define NBMAX 50
#define X1 0.0
#define X2 10.0

#define AH	0.1

void halo_distribution(int Nhi, int Nhf, int Nh, 
	vector Rcmh, vector Vcmh, real Mh)
{
	real rsen,fsen,mass_i,xrc,fh; 
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real rtheta,rphi,r,fv;
	bodyptr p;
	float tol,xb1[NBMAX],xb2[NBMAX];
	int nbr=NBMAX;

	fh = rpow( AH/Rmaxh + ((real) 1.0), 3 );
	mass_i = Mh / ((real) Nh);
	for (p = bodytab+Nhi; p < bodytab+Nhf+1; p++) {
        Type(p) = BODY;
		rtheta	= racos( ((real) 1.0) - ((real) 2.0) * xrandom(0.0,1.0) );
		rphi	= ((real) 2.0) * PI_D * xrandom(0.0,1.0);
		r = ((real) 1.5) * Rmaxh;
		while (r > Rmaxh) {
			xr = xrandom(0.0,1.0);
			while (xr >= ((real) 1.0)) {
				xr = xrandom(0.0,1.0);
			}
			xrc = rpow(xr/fh,1.0/3.0);
			r	= AH * xrc / ( 1 - xrc );
		}
		Mass(p)	= mass_i;
		rsen	= r * rsin(rtheta);
		Pos(p)[0]	= rsen*rcos(rphi);
		Pos(p)[1]	= rsen*rsin(rphi);
		Pos(p)[2]	= r*rcos(rtheta);

		rtheta	= racos( ((real) 1.0) - ((real) 2.0) * xrandom(0.0,1.0) );
		rphi	= ((real) 2.0) * PI_D * xrandom(0.0,1.0);
		xr = xrandom(0.0,1.0);
		sb = rsqrt( ( AH + ((real) 6.0)*r ) / ((real) 30.0) )/( r + AH );
		zbrak(fx_bulge_v,X1,X2,NZ,xb1,xb2,&nbr);
		if (nbr == 0) nrerror("halo_v: no roots founds");
		tol=(1.0e-6)*(rabs(xb1[1]+xb2[1]))/2.0;
		fv=zbrent(fx_bulge_v,xb1[1],xb2[1],tol);
		fsen	= fv*rsin(rtheta);
		Vel(p)[0]	= fsen*rcos(rphi);
		Vel(p)[1]	= fsen*rsin(rphi);
		Vel(p)[2]	= fv*rcos(rtheta);
	}
	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for (p = bodytab+Nhi; p < bodytab+Nhf+1; p++) {
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0] * Mass(p);
		sum_1y += Pos(p)[1] * Mass(p);
		sum_1z += Pos(p)[2] * Mass(p);
	}
	sum_0 = sum_0 / Mh;
	sum_1x = sum_1x / Mh;
	sum_1y = sum_1y / Mh;
	sum_1z = sum_1z / Mh;
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for (p = bodytab+Nhi; p < bodytab+Nhf+1; p++) {
		sum_2x += (Pos(p)[0]-sum_1x) * (Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += (Pos(p)[1]-sum_1y) * (Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += (Pos(p)[2]-sum_1z) * (Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x / Mh;
	sum_2y = sum_2y / Mh;
	sum_2z = sum_2z / Mh;
	printf("2do. momento: %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z,m2_rhoh(AH,Rmaxh));
	for (p = bodytab+Nhi; p < bodytab+Nhf+1; p++) {
		ADDV(Pos(p),Pos(p),Rcmh);
		ADDV(Vel(p),Vel(p),Vcmh);
	}
	printf("halo finished!\n");
}

real pr_h(real r, real dr, real Rmax, real ah)
{
	real Nh, fh;
	fh = rpow( ah/Rmax + ((real) 1.0), 3.0 );
	Nh = ((real) 3.0) * ah * fh * rsqr(r) * dr;
	return( Nh / rpow(r+ah,4.0) );
}

local real m2_rhoh(real a, real rmax)
{
	real out, fh;
	
	fh = rpow( a/rmax + ((real) 1.0), 3.0 );
	out = fh*a*(rmax*(12*rpow(a,3) + 30*rpow(a,2)*rmax + 
			22*a*rpow(rmax,2) + 3*rpow(rmax,3)) + 
			12*a*rpow(a + rmax,3)*rlog(a) - 
			12*a*rpow(a + rmax,3)*rlog(a + rmax))/
			rpow(a + rmax,3);
	return(out);
}
#undef NZ
#undef NBMAX
#undef X1
#undef X2


#define BULGE	"bulg%03d"
#define DISK	"disk%03d"
#define HALO	"halo%03d"

void outputdata_BDH_Galaxy(void)
{
	outputdata_body(BULGE,Nbi,Nbf);
	outputdata_body(DISK,Ndi,Ndf);
	outputdata_body(HALO,Nhi,Nhf);
}

#undef BULGE
#undef DISK
#undef HALO

#define BULGE1	"blg1%03d"
#define DISK1	"dsk1%03d"
#define HALO1	"hlo1%03d"
#define BULGE2	"blg2%03d"
#define DISK2	"dsk2%03d"
#define HALO2	"hlo2%03d"

void outputdata_Two_BDH_Galaxies_Coll(void)
{
	outputdata_body(BULGE1,Nbi1,Nbf1);
	outputdata_body(DISK1,Ndi1,Ndf1);
	outputdata_body(HALO1,Nhi1,Nhf1);
	outputdata_body(BULGE2,Nbi2,Nbf2);
	outputdata_body(DISK2,Ndi2,Ndf2);
	outputdata_body(HALO2,Nhi2,Nhf2);
}
#undef BULGE1
#undef DISK1
#undef HALO1
#undef BULGE2
#undef DISK2
#undef HALO2

#define outphasefile	"phs%03d"

local void galaxy_phasesp_out(void)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;

	bodyptr p;
	real prad,pphi,vphi,vr;

    sprintf(namebuf, outphasefile, gd.nstep);  
    if (stat(namebuf, &buf) != 0)        
        outstr = stropen(namebuf, "w");  
    else                                 
        outstr = stropen(namebuf, "a");  

	for ( p = bodytab+Ndi; p < bodytab+Ndf+1; p++) {
		prad = rsqrt( Pos(p)[0]*Pos(p)[0]+Pos(p)[1]*Pos(p)[1] );
		pphi=angle(0.0,0.0,Pos(p)[0],Pos(p)[1]);
		vphi = ( Pos(p)[0]*Vel(p)[1] - Pos(p)[1]*Vel(p)[0] )/prad;
		vr = ( Pos(p)[0]*Vel(p)[0] + Pos(p)[1]*Vel(p)[1] )/prad;
		fprintf(outstr,"%14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n",
			prad,pphi,Pos(p)[3],vr,vphi,Vel(p)[2],rabs(vphi));
	}
}

#define outavgdisk	"avgdisk"

local void BDH_avg_disk(int ni, int nf)
{
    char namebuf[256];
    struct stat buf;
    stream outstr;

	bodyptr p;
	real prad, vphi, vr;
	real vravg, vphiavg, vzavg, Norm, sigmar, sigmaphi, sigmaz;

    real amabs;
    vector tmpv;
#if !defined(THREEDIM)
	real tmps;
#endif
#if defined(THREEDIM)
	vector amvec;
#else
	real amvec;
#endif

    sprintf(namebuf, outavgdisk, gd.nstep);        
    if (stat(namebuf, &buf) != 0)               
        outstr = stropen(namebuf, "w");         
    else                                        
        outstr = stropen(namebuf, "a");         

#if defined(THREEDIM)
    CLRV(amvec); 
#else
	amvec=0.0;
#endif    

	Norm = ((real) nf - ni + 1);
	vravg = 0.0;
	vphiavg = 0.0;
	vzavg = 0.0;
	sigmar = 0.0;
	sigmaphi = 0.0;
	sigmaz = 0.0;

	for ( p = bodytab+ni; p < bodytab + nf; p++) {
		prad = rsqrt( Pos(p)[0]*Pos(p)[0]+Pos(p)[1]*Pos(p)[1] );
		vphi = ( Pos(p)[0]*Vel(p)[1] - Pos(p)[1]*Vel(p)[0] )/prad;
		vr = ( Pos(p)[0]*Vel(p)[0] + Pos(p)[1]*Vel(p)[1] )/prad;
		if (prad == 0.0) {
			vphi=0.0;
		}else{
			vphi = ( Pos(p)[0]*Vel(p)[1] - Pos(p)[1]*Vel(p)[0] )/prad;
		}
		vravg += vr;
		vphiavg += vphi;
		vzavg += Vel(p)[2];
#if defined(THREEDIM)
        CROSSVP(tmpv, Vel(p), Pos(p));          
        MULVS(tmpv, tmpv, Mass(p));
        ADDV(amvec, amvec, tmpv);
#endif
#if defined(TWODIM)
        CROSSVP(tmps, Vel(p), Pos(p));          
        amvec += tmps * Mass(p);
#endif
#if defined(ONEDIM)
		amvec = 0.0;
#endif
	}

#if defined(THREEDIM)
    ABSV(amabs, amvec);                         
#else
	amabs=rabs(amvec);
#endif

	vravg = vravg / Norm;
	vphiavg = vphiavg / Norm;
	vzavg = vzavg / Norm;
	for ( p = bodytab+ni; p < bodytab + nf; p++) {
		prad = rsqrt( Pos(p)[0]*Pos(p)[0]+Pos(p)[1]*Pos(p)[1] );
		vphi = ( Pos(p)[0]*Vel(p)[1] - Pos(p)[1]*Vel(p)[0] )/prad;
		vr = ( Pos(p)[0]*Vel(p)[0] + Pos(p)[1]*Vel(p)[1] )/prad;
		if (prad == 0.0) {
			vphi=0.0;
		}else{
			vphi = ( Pos(p)[0]*Vel(p)[1] - Pos(p)[1]*Vel(p)[0] )/prad;
		}
		sigmar += rsqr( vr - vravg );
		sigmaphi += rsqr( vphi - vphiavg );
		sigmaz += rsqr( Vel(p)[2] - vzavg );
	}
	sigmar = rsqrt(sigmar)/Norm;
	sigmaphi = rsqrt(sigmaphi)/Norm;
	sigmaz = rsqrt(sigmaz)/Norm;

	out_real_mar(outstr, gd.tnow);
	out_real_mar(outstr, sigmar);
	out_real_mar(outstr, sigmaphi);
	out_real_mar(outstr, sigmaz);
	out_real(outstr, amabs);
	fclose(outstr);
}

