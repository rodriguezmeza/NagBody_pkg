/*==============================================================================
	MODULE: models.c			[models]
	Written by: M.A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Rutines to create several types of models
	Language: C
	Use: 'testdata();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: July 23, 2007
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#define global
#include "globaldefs.h"


local void model_string_to_int(string, int *);

local void Model3D_Plummer(void);
local void Model3D_Plummer_Finite(void);
local void Freefall_Plummer_model(void);
local void UniformSphere_model(void);

local void Gaussian_Spheroid_model(void);
local void PertGaussianSphere_model(void);
local void PertGaussianSphere_d_model(void);
local void K_PertGaussianSphere_model(void);
local void PertIsothermalSphere_model(void);

local void Disk_Section(void);
local void Disk_Complete(void);
static float fx_disk_r(float);
static float fx_disk_z(float);
static float fx_disk_vr(float);
static float fx_disk_vphi(float);
local real m2_rhod(real, real, real);
local real xz, sr, v0;
local real Rmaxd;

local void UniformCubicBox_model(void);

//----------------- PASAR ESTAS RUTINAS A DATANALY_GRAV ------------------------
local void two_dist_phasesp_out(void);
local void two_dist_avg();

local void two_dist_cofm(void);				// Rutina para recalcular el CM


static float fx_pgs_phi(float);
static float fx_pgs_r(float);
static float fx_pgs_r_d(float);

local real Rmaxdist1, Rmaxdist2;
local real mtot;							// Usados en two_dist_cofm
local vector cmpos;                             
local vector cmvel;                             

local real xr,xphi;								
local real a_pgs,m_pgs,R_pgs,b_pgs;				

#define PLUMMER						0
#define PLUMMERFINITE				1
#define UNIFORMSPHERE				4
#define FREEFALLPLUMMER				5

#define PERTGAUSSIANSPHERE			302
#define PERTGAUSSIANSPHERED			303
#define K_PERTGAUSSIANSPHERE		304
#define GAUSSIAN_SPHEROID			305
#define PERTISOTHERMALSPHERE		308

#define UNIFORMCUBICBOX				400

#define GALAXY						500
#define TWO_GALAXIES_COLL			501

#define DISK_SECTION				601
#define DISK_COMPLETE				602

void testdata(void)
{
	int model_int;

	model_string_to_int(cmd.model, &model_int);
	switch (model_int){

	case PLUMMER: Model3D_Plummer(); break;
	case PLUMMERFINITE: Model3D_Plummer_Finite(); break;
	case FREEFALLPLUMMER: Freefall_Plummer_model(); break;

	case UNIFORMSPHERE: UniformSphere_model(); break;

	case PERTGAUSSIANSPHERE: PertGaussianSphere_model(); break;
	case PERTGAUSSIANSPHERED: PertGaussianSphere_d_model(); break;
	case K_PERTGAUSSIANSPHERE: K_PertGaussianSphere_model(); break;
	case GAUSSIAN_SPHEROID: Gaussian_Spheroid_model(); break;
	case PERTISOTHERMALSPHERE: PertIsothermalSphere_model(); break;

	case UNIFORMCUBICBOX: UniformCubicBox_model(); break;

	case GALAXY: BDH_Galaxy(); break;
	case TWO_GALAXIES_COLL: Two_BDH_Galaxies_Coll(); break;

	case DISK_SECTION: Disk_Section(); break;
	case DISK_COMPLETE: Disk_Complete(); break;

	default: error("\nUnknown model type %s\n\n",cmd.model);
	}
}

local void model_string_to_int(string model_str,int *model_int)
{
	*model_int = -1;
	if (strcmp(model_str,"plummer") == 0)				*model_int=PLUMMER;
	if (strcmp(model_str,"plummer-finite") == 0)		*model_int=1;
	if (strcmp(model_str,"uniformsphere") == 0)			*model_int=4;
	if (strcmp(model_str,"freefall-plummer") == 0)		*model_int=5;

	if (strcmp(model_str,"pertgaussiansphere") == 0)	*model_int=302;
	if (strcmp(model_str,"pertgaussiansphere_d") == 0)	*model_int=303;
	if (strcmp(model_str,"k-pertgaussiansphere") == 0)	*model_int=304;
	if (strcmp(model_str,"gaussianspheroid") == 0)		*model_int=GAUSSIAN_SPHEROID;
	if (strcmp(model_str,"pertisothermalsphere") == 0)	*model_int=308;

	if (strcmp(model_str,"uniformcubicbox") == 0)		*model_int=400;

	if (strcmp(model_str,"galaxy") == 0)				*model_int=GALAXY;
	if (strcmp(model_str,"two_galaxies_coll") == 0)		*model_int=501;
	if (strcmp(model_str,"bdhg_galaxy") == 0)			*model_int=502;
	if (strcmp(model_str,"bdhgdm_galaxy") == 0)			*model_int=503;
	if (strcmp(model_str,"disk-section") == 0)			*model_int=DISK_SECTION;
	if (strcmp(model_str,"disk-complete") == 0)			*model_int=DISK_COMPLETE;

}

#undef PLUMMER
#undef PLUMMERFINITE
#undef UNIFORMSPHERE
#undef FREEFALLPLUMMER

#undef PERTGAUSSIANSPHERE
#undef PERTGAUSSIANSPHERED
#undef K_PERTGAUSSIANSPHERE
#undef GAUSSIAN_SPHEROID
#undef PERTISOTHERMALSPHERE

#undef UNIFORMCUBICBOX

#undef GALAXY
#undef TWO_GALAXIES_COLL


// COMIENZAN MODELOS DE PLUMMER (Unificarlos en uno solo!!!) -------------------

#define MFRAC  0.999  

local void Model3D_Plummer(void)
{
    real rsc, vsc, r, v, x, y;
    vector rcm, vcm;
    bodyptr p;

    strcpy(gd.model_comment, "3D Plummer Model");
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
                                                
    rsc = (3 * PI) / 16;                        
    vsc = rsqrt(1.0 / rsc);                     
    CLRV(rcm);                                  
    CLRV(vcm);                                  
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        Type(p) = BODY;                         
        Mass(p) = 1.0 / cmd.nbody;                  
        x = xrandom(0.0, MFRAC);                
        r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1); 
        pickshell(Pos(p), NDIM, rsc * r);       
        do {                                    
            x = xrandom(0.0, 1.0);              
            y = xrandom(0.0, 0.1);              
        } while (y > x*x * rpow(1 - x*x, 3.5)); 
        v = x * rsqrt(2.0 / rsqrt(1 + r*r));    
        pickshell(Vel(p), NDIM, vsc * v);       
        ADDMULVS(rcm, Pos(p), 1.0 / cmd.nbody);     
        ADDMULVS(vcm, Vel(p), 1.0 / cmd.nbody);     
    }
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        SUBV(Pos(p), Pos(p), rcm);              
        SUBV(Vel(p), Vel(p), vcm);              
    }

    for (p = bodytab; p < bodytab+cmd.nbody; p++) {	// Adding angular velocity
		Pos(p)[0] += cmd.cmx;						// cmpos and cmvel...
		Pos(p)[1] += cmd.cmy;
		Pos(p)[2] += cmd.cmz;
		Vel(p)[0] += cmd.vcmx - cmd.omega0 * Pos(p)[1];
		Vel(p)[1] += cmd.vcmy + cmd.omega0 * Pos(p)[0];
		Vel(p)[2] += cmd.vcmz;
    }
}

local void Model3D_Plummer_Finite(void)
{
    real rsc, vsc, r, v, x, y;
    vector rcm, vcm;
    bodyptr p;
	real Rmax, s;

	strcpy(gd.model_comment, "3D Plummer Model Finite R");

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
                                         
    printf("\n Total mass: %g\n\n",cmd.Mtotal);

//	Rmax = 2.0;
    rsc = (3 * PI) / 16;                        
    vsc = rsqrt(1.0 / rsc);                     
    CLRV(rcm);                                  
    CLRV(vcm);                                  
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        Type(p) = BODY;
//        Mass(p) = 1.0 / cmd.nbody;
        Mass(p) = cmd.Mtotal / cmd.nbody;
		do {
			x = xrandom(0.0, MFRAC);                
			r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1); 
			pickshell(Pos(p), NDIM, rsc * r);       
			ABSV(s,Pos(p));
		} while ( s > cmd.Rmax );
        do {                                    
            x = xrandom(0.0, 1.0);              
            y = xrandom(0.0, 0.1);              
        } while (y > x*x * rpow(1 - x*x, 3.5)); 
        v = x * rsqrt(2.0 / rsqrt(1 + r*r));    
        pickshell(Vel(p), NDIM, vsc * v);       
        ADDMULVS(rcm, Pos(p), 1.0 / cmd.nbody);     
        ADDMULVS(vcm, Vel(p), 1.0 / cmd.nbody);     
    }
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        SUBV(Pos(p), Pos(p), rcm);              
        SUBV(Vel(p), Vel(p), vcm);              
    }

    for (p = bodytab; p < bodytab+cmd.nbody; p++) {	// Adding angular velocity
		Pos(p)[0] += cmd.cmx;						// cmpos and cmvel...
		Pos(p)[1] += cmd.cmy;
		Pos(p)[2] += cmd.cmz;
		Vel(p)[0] += cmd.vcmx - cmd.omega0 * Pos(p)[1];
		Vel(p)[1] += cmd.vcmy + cmd.omega0 * Pos(p)[0];
		Vel(p)[2] += cmd.vcmz;
    }
}

local void Freefall_Plummer_model(void)
{
    real rsc, r, x, y;
    vector rcm;
    bodyptr p;

    strcpy(gd.model_comment, "Free falling Plummer Model");

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

    rsc = (3 * PI) / 16;
    CLRV(rcm);
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        Type(p) = BODY;                         
        Mass(p) = 1.0 / cmd.nbody;                  
        x = xrandom(0.0, MFRAC);                
        r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1); 
        pickshell(Pos(p), NDIM, rsc * r);       
        do {                                    
            x = xrandom(0.0, 1.0);
            y = xrandom(0.0, 0.1);
        } while (y > x*x * rpow(1 - x*x, 3.5));
		CLRV(Vel(p));
        ADDMULVS(rcm, Pos(p), 1.0 / cmd.nbody);     
    }
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        SUBV(Pos(p), Pos(p), rcm);              
}

// TERMINAN MODELOS DE PLUMMER -------------------------------------------------


// COMIENZAN DISTRIBUCIONES ESFERICAS (Unificar estas rutinas!!!!) -------------

// Esfera Uniforme
// Uso: models model_type=uniformsphere
// Units: M=R=G=1

// With this routine we may obtain the standard model: a sphere with uniform
// density rotating with angular velocity omega0=1. We have to do absvel=0
// vx=vy=vz=0. We may also obtain the standard model with random velocity,
// just do omega0=0, vx=vy=vz=0.

//Modelo para probar el colapso adiabático unidimensional 
//de una esfera de gas. Ver: Hernquist & Katz; Sprigel.

//Notas: El archivo de parámetros debe tener T (=10) 
//diferente de cero. Si T=0 aparecen problemas con la
//búsqueda de vecinos: ngb_treefind. (Averiguar
//por que sucede esto, T no debe de influir en la
//búsqueda de vecinos).
//La información que se necesita en el archivo de
//parámetros requiere de las unidades de longitud,
//de masa, velocidad y el valor de G.

// GENERA UNA ESFERA UNIFORME. El formato en que se guarda la
// define como esfera de gas o esfera de n-cuerpos no colisionales.
// The output format defines the system as a gas or as a collisionless N-body. 

local void UniformSphere_model(void)
{
	real rsen,mass_i;
	real rphi,rtheta,r;
	int np;
	bodyptr p;
	vector cmpos, cmvel, tmpv;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z,the_sum_2;

	strcpy(gd.model_comment, "Uniform sphere model");

	mass_i = cmd.Mtotal/((real) cmd.nbody);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

	CLRV(cmpos);
    CLRV(cmvel);

	for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
		Type(p) = BODY;
		Mass(p) = mass_i;
		rphi	= 2.0 * PI_D * xrandom(0.0, 1.0);
		rtheta	= racos(1.0 - 2.0 * xrandom(0.0, 1.0));
//		r = Rmax*rsqrt(xrandom(0.0,1.0));					// Checar cual es cual...
		r		= cmd.Rmax*rpow(xrandom(0.0, 1.0),1.0/3.0);
		rsen	= r*rsin(rtheta);
		Pos(p)[0] = rsen*rcos(rphi);
		Pos(p)[1] = rsen*rsin(rphi);
		Pos(p)[2] = r*rcos(rtheta);
		Vel(p)[0] = cmd.absvel * xrandom(-1.0, 1.0);
		Vel(p)[1] = cmd.absvel * xrandom(-1.0, 1.0);
		Vel(p)[2] = cmd.absvel * xrandom(-1.0, 1.0);
        MULVS(tmpv, Pos(p), Mass(p));
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));
        ADDV(cmvel, cmvel, tmpv);
	}
    DIVVS(cmpos, cmpos, cmd.Mtotal);
    DIVVS(cmvel, cmvel, cmd.Mtotal);
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {	// Correcting cmpos and cmvel...
		SUBV(Pos(p),Pos(p),cmpos)
		SUBV(Vel(p),Vel(p),cmvel)
    }
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {	// Adding angular velocity
		Pos(p)[0] += cmd.cmx;						// cmpos and cmvel...
		Pos(p)[1] += cmd.cmy;
		Pos(p)[2] += cmd.cmz;
		Vel(p)[0] += cmd.vcmx - cmd.omega0 * Pos(p)[1];
		Vel(p)[1] += cmd.vcmy + cmd.omega0 * Pos(p)[0];
		Vel(p)[2] += cmd.vcmz;
    }

	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0] * Mass(p);
		sum_1y += Pos(p)[1] * Mass(p);
		sum_1z += Pos(p)[2] * Mass(p);
	}
	sum_0=sum_0/cmd.Mtotal;
	sum_1x = sum_1x/cmd.Mtotal;
	sum_1y = sum_1y/cmd.Mtotal;
	sum_1z = sum_1z/cmd.Mtotal;
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_2x += (Pos(p)[0]-sum_1x) * (Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += (Pos(p)[1]-sum_1y) * (Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += (Pos(p)[2]-sum_1z) * (Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x/cmd.Mtotal;
	sum_2y = sum_2y/cmd.Mtotal;
	sum_2z = sum_2z/cmd.Mtotal;
	the_sum_2 = ((real) 0.5)*rsqr(cmd.Rmax);	// Checar este valor analitico...
	printf("2do. momento (numerical and analytical): %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z, the_sum_2);
}

local void Gaussian_Spheroid_model(void)
{
	real mass_i;
	bodyptr p;
	vector cmpos, cmvel, tmpv;

	strcpy(gd.model_comment, "Gaussian-Spheroid Model");

	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
	mass_i = cmd.Mtotal/((real) cmd.nbody);

	CLRV(cmpos);
    CLRV(cmvel);

	for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        Type(p) = BODY;
        Mass(p) = mass_i;
		Pos(p)[0]	= grandom(0.0,cmd.a_spheroid);
		Pos(p)[1]	= grandom(0.0,cmd.b_spheroid);
		Pos(p)[2]	= grandom(0.0,cmd.c_spheroid);
		Vel(p)[0] = cmd.absvel * xrandom(-1.0, 1.0);
		Vel(p)[1] = cmd.absvel * xrandom(-1.0, 1.0);
		Vel(p)[2] = cmd.absvel * xrandom(-1.0, 1.0);
        MULVS(tmpv, Pos(p), Mass(p));
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));
        ADDV(cmvel, cmvel, tmpv);
	}
    DIVVS(cmpos, cmpos, cmd.Mtotal);
    DIVVS(cmvel, cmvel, cmd.Mtotal);
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {	// Correcting cmpos and cmvel...
		SUBV(Pos(p),Pos(p),cmpos)
		SUBV(Vel(p),Vel(p),cmvel)
    }
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {	// Adding angular velocity
		Pos(p)[0] += cmd.cmx;						// cmpos and cmvel...
		Pos(p)[1] += cmd.cmy;
		Pos(p)[2] += cmd.cmz;
		Vel(p)[0] += cmd.vcmx - cmd.omega0 * Pos(p)[1];
		Vel(p)[1] += cmd.vcmy + cmd.omega0 * Pos(p)[0];
		Vel(p)[2] += cmd.vcmz;
    }
    for (p = bodytab; p < bodytab+cmd.nbody; p++)	// Adding Id
		Id(p) = p-bodytab+1;
}

#define NZ 50
#define NBMAX 20
#define X1 0.0
#define X2 7.0

#define NZY 50
#define NBYMAX 50
#define Y1 0.0
#define Y2 100.0

local void	PertGaussianSphere_model(void)
{
	real rsen,Mtotal,mass_i,omega_0,rho0,b,R_cgs;
	real rphi,rtheta,r;
	int np;
	bodyptr p;
	float tol,xb1[NBMAX],xb2[NBMAX];
	float yb1[NBYMAX],yb2[NBYMAX];
	int nbr=NBMAX;
	int nby=NBYMAX;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real the_sum_2,Rob,Rob2,numerador,denominador;
    stream outstr;
	real InternalEnergy;

	strcpy(gd.model_comment, "Perturbed Gaussian Sphere Model");

	b_pgs = 0.57776137;
	Rob = ((real) 1.0)/b_pgs;
	Rob2 = Rob*Rob;
	Mtotal = 1.0;
	R_pgs = 1.0;
	omega_0 = 0.9705314;				
	cmd.SoundSpeed = 0.369331351;
	a_pgs = 0.1;
	m_pgs = 2.0;
	printf("input num. of particles of sphere [power of 2]:");
	scanf("%d", &np);
	cmd.nbody=rpow(2, np);
	mass_i = Mtotal/((real) cmd.nbody);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
	for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        Type(p) = BODY;
        Mass(p) = mass_i;

		do {
		xr = xrandom(0.0,1.0);
		zbrak(fx_pgs_r,Y1,Y2,NZY,yb1,yb2,&nby);
		tol=(1.0e-6)*(yb1[1]+yb2[1])/2.0;
		r  = b_pgs * zbrent(fx_pgs_r,yb1[1],yb2[1],tol);
		} while (r > R_pgs);

		do{											
			xphi = xrandom(0.0,1.0);
			zbrak(fx_pgs_phi,X1,X2,NZ,xb1,xb2,&nbr);
		} while(nbr == 0);
		if (nbr == 0) nrerror("PertGaussianSphere: no root found");
		tol=(1.0e-6)*(xb1[1]+xb2[1])/2.0;
		rphi  =zbrent(fx_pgs_phi,xb1[1],xb2[1],tol);

		rtheta	= racos(1.0 - 2.0 * xrandom(0.0, 1.0));

		rsen	= r*rsin(rtheta);

		Pos(p)[0]	= rsen*rcos(rphi);
		Pos(p)[1]	= rsen*rsin(rphi);
		Pos(p)[2]	= r*rcos(rtheta);
		Vel(p)[0]	= - omega_0 * Pos(p)[1];
		Vel(p)[1]	=   omega_0 * Pos(p)[0];
		Vel(p)[2]	= 0.0;
	}
	
	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0] * Mass(p);
		sum_1y += Pos(p)[1] * Mass(p);
		sum_1z += Pos(p)[2] * Mass(p);
	}
	sum_1x = sum_1x/sum_0;
	sum_1y = sum_1y/sum_0;
	sum_1z = sum_1z/sum_0;
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_2x += (Pos(p)[0]-sum_1x) * (Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += (Pos(p)[1]-sum_1y) * (Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += (Pos(p)[2]-sum_1z) * (Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x/sum_0;
	sum_2y = sum_2y/sum_0;
	sum_2z = sum_2z/sum_0;
	Rob = R_pgs/b_pgs;
	Rob2 = rsqr(Rob);
	numerador = -((real) 0.5)*Rob*(Rob2+((real) 1.5))*rexp(-Rob2)+
				((real) 3.0/8.0)*SQRTPI*merff(Rob);
	denominador = -((real) 0.5)*Rob*rexp(-Rob2)+
				((real) 0.25)*SQRTPI*merff(Rob);
	the_sum_2 = rsqr(b_pgs)*numerador/denominador;
	printf("2do. momento (numerical and analytical): %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z, the_sum_2);
	printf("Perturbed Gaussian Sphere model finished\n");
}

local void	PertIsothermalSphere_model(void)
{
/*
Rutina para generar la condicion inicial de una esfera
con densidad uniforme con una perturbacion a*cos(phi).
La esfera est'a en rotaci'on.

Desarrollado: Julio, 2001.
Fechas de actualizacion: Julio, 2001;

Uso: model model_type=pertisothermalsphere ofmt=gadget11-ascii

Los parametros para tener en cuenta en la generaci'on de
condiciones iniciales y la simulaciones correspondientes son:
SoundSpeed, a_pgs, omega_0. Tambien recordar que las simulaciones
se deben hacer considerando que es un colapso isot'ermico.

Considerar el buen uso de SoundSpeed ya que es necesario
para calcular la energia interna. El archivo debe salvarse
en formato gadget-sph-ascii.

Other posible values for the parameters are:

(1)	SoundSpeed = 0.31547;
	a_pgs = 0.5; 
	omega_0 = 0.775183;
	m_pgs = 2.0;

(2)	SoundSpeed = 0.394338;
	a_pgs = 0.1;
	omega_0 = 0.698783;				
	m_pgs = 2.0;

	Mtotal = 1.0;
	R_pgs = 1.0;
*/

	real rsen,Mtotal,mass_i,omega_0; 
	real rphi,rtheta,r;
	int np;
	bodyptr p;
	float tol,xb1[NBMAX],xb2[NBMAX];
	int nbr=NBMAX;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real the_sum_2,Rob,Rob2,numerador,denominador;

	strcpy(gd.model_comment, "Perturbed Isothermal Uniform Sphere Model");

	R_pgs = cmd.Rmax;					// Given in the command line
	a_pgs = cmd.a_p;					// Given in the command line
	omega_0 = cmd.omega0;				// Given in the command line
	m_pgs = cmd.m_p;					// Given in the command line

	mass_i = Mtotal/((real) cmd.nbody);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
	for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        Type(p) = BODY;
        Mass(p) = mass_i;

		xr = xrandom(0.0,1.0);
		r = rpow(xr,1.0/3.0);

		do{											
			xphi = xrandom(0.0,1.0);
			zbrak(fx_pgs_phi,X1,X2,NZ,xb1,xb2,&nbr);
		} while(nbr == 0);
		if (nbr == 0) nrerror("PertIsothermalSphere: no root found");
		tol=(1.0e-6)*(xb1[1]+xb2[1])/2.0;
		rphi  =zbrent(fx_pgs_phi,xb1[1],xb2[1],tol);

		rtheta	= racos(1.0 - 2.0 * xrandom(0.0, 1.0));

		rsen	= r*rsin(rtheta);

		Pos(p)[0]	= rsen*rcos(rphi);
		Pos(p)[1]	= rsen*rsin(rphi);
		Pos(p)[2]	= r*rcos(rtheta);
		Vel(p)[0]	= - omega_0 * Pos(p)[1];
		Vel(p)[1]	=   omega_0 * Pos(p)[0];
		Vel(p)[2]	= 0.0;
	}
	
	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0] * Mass(p);
		sum_1y += Pos(p)[1] * Mass(p);
		sum_1z += Pos(p)[2] * Mass(p);
	}
	sum_1x = sum_1x/sum_0;
	sum_1y = sum_1y/sum_0;
	sum_1z = sum_1z/sum_0;
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_2x += (Pos(p)[0]-sum_1x) * (Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += (Pos(p)[1]-sum_1y) * (Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += (Pos(p)[2]-sum_1z) * (Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x/sum_0;
	sum_2y = sum_2y/sum_0;
	sum_2z = sum_2z/sum_0;
	the_sum_2 = (3.0/5.0);
	printf("2do. momento (numerical and analytical): %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z, the_sum_2);
}

local void	PertGaussianSphere_d_model(void)
{
	real rsen,Mtotal,mass_i,omega_0,rho0,b;
	real rphi,rtheta,r,deltar,deltatheta,deltaphi,x1,x2,xtheta;
	int nphi,ntheta,nr,iphi,itheta,ir,ip;
	bodyptr p;
	float tol,xb1[NBMAX],xb2[NBMAX];
	float yb1[NBYMAX],yb2[NBYMAX];
	int nbr=NBMAX;
	int nby=NBYMAX;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real the_sum_2,Rob,Rob2,numerador,denominador;

	strcpy(gd.model_comment, "Perturbed Gaussian Sphere deterministic Model");

    printf("%s\n", gd.model_comment);

	Mtotal = 1.0;
	R_pgs = 1.0;
	omega_0 = 0.9705314;				
	rho0 = 1.0;
	a_pgs = 0.1;
	m_pgs = 2.0;
	b_pgs = 0.577;
	printf("input num. of particles of sphere [nphi,ntheta,nr]:");
	scanf("%d %d %d", &nphi, &ntheta, &nr);
	cmd.nbody=nphi*ntheta*nr;
	mass_i = Mtotal/((real) cmd.nbody);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

	x1=0.00001;
	x2=0.99999;
	deltar = (x2-x1)/((real) nr-1);
	deltatheta = (x2-x1)/((real) ntheta-1);
	deltaphi = (x2-x1)/((real) nphi-1);

	p = bodytab;
	for (ir = 1; ir <= nr; ir++) { 
	  xr = x1 + ((real) ir-1)*deltar;
	  zbrak(fx_pgs_r_d,Y1,Y2,NZY,yb1,yb2,&nby);
	  if (nby == 0) nrerror("PertGaussianSphere: no root found for r");
	  tol=(1.0e-6)*(yb1[1]+yb2[1])/2.0;
	  r  = b_pgs * zbrent(fx_pgs_r_d,yb1[1],yb2[1],tol);
	  for (iphi = 1; iphi <= nphi; iphi++) { 
		xphi = x1 + ((real) iphi-1)*deltaphi;
		zbrak(fx_pgs_phi,X1,X2,NZ,xb1,xb2,&nbr);
		if (nbr == 0) nrerror("PertGaussianSphere: no root found for phi");
		tol=(1.0e-6)*(xb1[1]+xb2[1])/2.0;
		rphi  =zbrent(fx_pgs_phi,xb1[1],xb2[1],tol);
	    for (itheta = 1; itheta <= ntheta; itheta++) { 
		  xtheta = x1 + ((real) itheta-1)*deltatheta;
		  rtheta	= racos(1.0 - 2.0 * xtheta);
		  rsen	= r*rsin(rtheta);
          Type(p) = BODY;
          Mass(p) = mass_i;
		  Pos(p)[0]	= rsen*rcos(rphi);
		  Pos(p)[1]	= rsen*rsin(rphi);
		  Pos(p)[2]	= r*rcos(rtheta);
		  Vel(p)[0]	= - omega_0 * Pos(p)[1];
		  Vel(p)[1]	=   omega_0 * Pos(p)[0];
		  Vel(p)[2]	= 0.0;
	      p++;
		}
	  }
	}
	
	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0] * Mass(p);
		sum_1y += Pos(p)[1] * Mass(p);
		sum_1z += Pos(p)[2] * Mass(p);
	}
	sum_0=sum_0/Mtotal;
	sum_1x = sum_1x/Mtotal;
	sum_1y = sum_1y/Mtotal;
	sum_1z = sum_1z/Mtotal;
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_2x += (Pos(p)[0]-sum_1x) * (Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += (Pos(p)[1]-sum_1y) * (Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += (Pos(p)[2]-sum_1z) * (Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x/Mtotal;
	sum_2y = sum_2y/Mtotal;
	sum_2z = sum_2z/Mtotal;
	Rob = R_pgs/b_pgs;
	Rob2 = rsqr(Rob);
	numerador = -((real) 0.5)*Rob*(Rob2+((real) 1.5))*rexp(-Rob2)+
				((real) 3.0/8.0)*SQRTPI*merff(Rob);
	denominador = -((real) 0.5)*Rob*rexp(-Rob2)+
				((real) 0.25)*SQRTPI*merff(Rob);
	the_sum_2 = rsqr(b_pgs)*numerador/denominador;
	printf("2do. momento (numerical and analytical): %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z, the_sum_2);
	printf("Perturbed Gaussian Sphere deterministic model finished\n");
}

static float fx_pgs_phi(float phi)
{
	real arg;
	float out;
	arg = m_pgs*phi;
	out = (a_pgs/m_pgs) * rsin(arg) + phi - TWOPI*xphi;
	return(out);
}

static float fx_pgs_r(float t)
{
	float out,Rob;
	Rob = R_pgs/b_pgs;
/*	out = (merff(t)-xr)*SQRTPI/((real) 4.0) - (t/((real) 2.0))*rexp(-t*t); */
	out = ((real) 0.25) * SQRTPI * ( merff(t) - xr*merff(Rob) )	
			+ ((real) 0.5) * ( xr*Rob*rexp(-Rob*Rob) - t*rexp(-t*t) );
	return(out);
}

static float fx_pgs_r_d(float t)
{
	float out;
	out = ( merff(t)-xr*merff(R_pgs/b_pgs) )*SQRTPI/((real) 4.0) 
		- ( t/((real) 2.0) )*rexp(-t*t) 
		+ ( R_pgs/(b_pgs*((real) 2.0)) )*rexp(-rsqr(R_pgs/b_pgs))*xr;
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

local void	K_PertGaussianSphere_model(void)
{
	real Rmax,rsen,Mtotal,mass_i,omega_0,rho0,b;
	real rphi,rtheta,r,pert,x0,y0;
	int np;
	bodyptr p,q;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real the_sum_2,Rob,Rob2,numerador,denominador;
    stream instr;
	int i,j,k,pp,jp;
	double xran[50000];
	double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10;
	double x1,x,y1,y,z1,z,delta,delta3,fran;

    real rsc, vsc, v; 
	vector rcm,vcm;

	strcpy(gd.model_comment, "K Perturbed Gaussian Sphere Model");

    printf("%s\n", gd.model_comment);

	Mtotal = 1.0;
	Rmax = 1.0;
	omega_0 = 0.9705314;				
	a_pgs = 0.1;
	m_pgs = 2.0;
	b_pgs = 0.577*Rmax;
	rho0 = Mtotal/rpow(b_pgs*SQRTPI,3.0);
	printf("input num. of divisions of the cube:");
	scanf("%d", &np);

	x1 = y1 = z1 = -Rmax;
	delta = ((real) 2.0)*Rmax / ((real) np-1);
	cmd.nbody=0;

	for (i=1;i<=np;i++) {
		x = x1 + ( ((real) i) - ((real) 0.5) ) * delta;
		for (j=1;j<=np;j++){
			y = y1 + ( ((real) j) - ((real) 0.5) ) * delta;
			for (k=1;k<=np;k++){
				z = z1 + ( ((real) k) - ((real) 0.5) ) * delta;
				r = rsqrt( x*x + y*y + z*z );
				if (r <= Rmax) cmd.nbody++;
			}
		}
	}

	printf("number of bodies: %d\n",cmd.nbody);

	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

	pp=0;
	p=bodytab;
	fran=1.0e-1;								

	delta3 = rpow(delta,3.0);

	for (i=1;i<=np;i++) {
		x = x1 + ( ((real) i) - ((real) 0.5) ) * delta;
		for (j=1;j<=np;j++){
			y = y1 + ( ((real) j) - ((real) 0.5) ) * delta;
			for (k=1;k<=np;k++){
				z = z1 + ( ((real) k) - ((real) 0.5) ) * delta;
				r = rsqrt( x*x + y*y + z*z );
				if (r <= Rmax) {
					Type(p) = BODY;                         

					Pos(p)[0]=x*(1.0+fran*xrandom(0.0,1.0));
					Pos(p)[1]=y*(1.0+fran*xrandom(0.0,1.0));
					Pos(p)[2]=z*(1.0+fran*xrandom(0.0,1.0));

					pp++;
					Vel(p)[0]=-omega_0*y;
					Vel(p)[1]= omega_0*x;
					Vel(p)[2]=0.0;

					x0 = Pos(p)[0];
					y0 = Pos(p)[1];
					rphi = angle(0.0,0.0,x0,y0);
					pert = ((real) 1.0) + a_pgs*rcos(m_pgs*rphi);

					Mass(p)=rho0*rexp(-rsqr(r/b_pgs))*delta3*pert;
					p++;
				}
			}
		}
	}

	printf("pp: %d\n",pp);
	if ( pp != cmd.nbody ) error("nbody not equal pp");

	sum_0=0.0;
	sum_1x=0.0;
	sum_1y=0.0;
	sum_1z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_0 += Mass(p);
		sum_1x += Pos(p)[0] * Mass(p);
		sum_1y += Pos(p)[1] * Mass(p);
		sum_1z += Pos(p)[2] * Mass(p);
	}
	sum_0=sum_0/Mtotal;
	sum_1x = sum_1x/Mtotal;
	sum_1y = sum_1y/Mtotal;
	sum_1z = sum_1z/Mtotal;
	printf("momento de orden cero: %14.7f\n",sum_0);
	printf("1er. momento: %14.7f %14.7f %14.7f\n",sum_1x,sum_1y,sum_1z);
	sum_2x=0.0;
	sum_2y=0.0;
	sum_2z=0.0;
	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		sum_2x += (Pos(p)[0]-sum_1x) * (Pos(p)[0]-sum_1x) * Mass(p);
		sum_2y += (Pos(p)[1]-sum_1y) * (Pos(p)[1]-sum_1y) * Mass(p);
		sum_2z += (Pos(p)[2]-sum_1z) * (Pos(p)[2]-sum_1z) * Mass(p);
	}
	sum_2x = sum_2x/Mtotal;
	sum_2y = sum_2y/Mtotal;
	sum_2z = sum_2z/Mtotal;
	Rob = Rmax/b_pgs;
	Rob2 = rsqr(Rob);
	numerador = -((real) 0.5)*Rob*(Rob2+((real) 1.5))*rexp(-Rob2)+
				((real) 3.0/8.0)*SQRTPI*merff(Rob);
	denominador = -((real) 0.5)*Rob*rexp(-Rob2)+
				((real) 0.25)*SQRTPI*merff(Rob);
	the_sum_2 = rsqr(b_pgs)*numerador/denominador;
	printf("2do. momento (numerical and analytical): %14.7f %14.7f\n",
		sum_2x+sum_2y+sum_2z, the_sum_2);
	printf("K_Perturbed Gaussian Sphere model finished\n");
}

// TERMINAN DISTRIBUCIONES ESFERICAS -------------------------------------------


/* 
Caja Cubica Uniforme
Posiciones y velocidades distribuidas aleatoriamente siguiendo
una distribucion uniforme.

Uso: model model_type=uniformcubicbox
*/

local void UniformCubicBox_model(void)
{
	real rsen,mass_i;
	real rphi,rtheta,r;
	int np;
	bodyptr p;
        real absvel;

	strcpy(gd.model_comment, "Uniform cubic box model");

	mass_i = cmd.Mtotal/((real) cmd.nbody);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
	for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
            Type(p) = BODY;
            Mass(p) = mass_i;
            Pos(p)[0]	= cmd.Rmax * xrandom(-1.0, 1.0);
            Pos(p)[1]	= cmd.Rmax * xrandom(-1.0, 1.0);
            Pos(p)[2]	= cmd.Rmax * xrandom(-1.0, 1.0);
            Vel(p)[0]	= absvel * xrandom(-1.0, 1.0);
            Vel(p)[1]	= absvel * xrandom(-1.0, 1.0);
            Vel(p)[2]	= absvel * xrandom(-1.0, 1.0);
	}
}


// COMIENZA RUTINA DISK_SECTION ------------------------------------------------

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
//#define W1 -2.0
//#define W2  5.0
#define W1 -5.0					// Funciona mejor con estos parametros ...
#define W2  10.0

#define AD	0.083333333333333333333333333333333		// alpha^-1 = 1/12
#define Z0	0.007

local void	Disk_Section(void)
{
	int Ndi, Ndf, Nd;
	vector Rcmd, Vcmd;
	real Md;

	real Ort_A,Ort_B,disk_Q,kappa,sphi,sigma,sz;
	real mass_i,vphi,vz,vr;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real rtheta,r,z;
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
	real Rmax, Rmin, thetamax, thetamin, zmax, zmin;
	bodyptr p, q, btabtmp;
	real xbox, ybox, zbox;
	int nbodytmp;
	vector cmpos, cmvel, tmpv;
	real tmass;
	int k;

//	Md = cmd.Mtotal;
	Md = 0.1875*cmd.Mtotal;
	Ndi = 0;
	Ndf = cmd.nbody-1;
	Nd = cmd.nbody;
	CLRV(Rcmd);
	CLRV(Vcmd);
	Md = 0.1875;
	Rmaxd = 1.0;
	Rmax = 0.226;					// 0.2+0.025 (1 kpc/ 40 kpc = 0.025)
	Rmin = 0.2;						//  8/40=0.2 que es 8 kpc
	thetamin = -0.0626;				// dl = rmin (thetmax - thetamin) = 1 kpc
	thetamax = 0.0626;				//    = 8 kpc ( 1/8 ) = 1 kpc
	zmin = -0.0125;
	zmax = 0.0125;

	xbox = 0.025;					// Simulation box of 1 kpc size 
	ybox = 0.025;
	zbox = 0.025;

	mass_i = Md / ((real) Nd);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

	for ( p = bodytab+Ndi; p < bodytab+Ndf+1; p++) {
        Type(p) = BODY;
		rtheta = xrandom(thetamin,thetamax);

		r = xrandom(Rmin,Rmax);

		xz = xrandom(0.0,1.0);
		zbrak(fx_disk_z,Y1,Y2,NZY,yb1,yb2,&nby);
		if (nby == 0) nrerror("disk_z: no roots founds");
		tol=(1.0e-6)*(rabs(yb1[1]+yb2[1]))/2.0;
		z=zbrent(fx_disk_z,yb1[1],yb2[1],tol);
//		z = xrandom(zmin,zmax);								// Uniform z

		Mass(p)	  = mass_i;
		Pos(p)[0] = r*rcos(rtheta);
		Pos(p)[1] = r*rsin(rtheta);
		Pos(p)[2] = z;

//		v0 = 1.0;
		v0 = 0.0;
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

        ADDMULVS(Rcmd, Pos(p), Mass(p));     
        ADDMULVS(Vcmd, Vel(p), Mass(p));     
	}

    DIVVS(Rcmd, Rcmd, Md);                  
    DIVVS(Vcmd, Vcmd, Md);

    for (p = bodytab+Ndi; p < bodytab+Ndf+1; p++) { 
        SUBV(Pos(p), Pos(p), Rcmd);              
        SUBV(Vel(p), Vel(p), Vcmd);              
    }

	fprintf(stdout,"CM position: %g %g %g\n",Rcmd[0],Rcmd[1],Rcmd[2]);
	fprintf(stdout,"CM velocity: %g %g %g\n",Vcmd[0],Vcmd[1],Vcmd[2]);

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
		sum_2x+sum_2y+sum_2z, m2_rhod(AD,Z0,Rmaxd));

    btabtmp = (bodyptr) allocate(cmd.nbody * sizeof(body));
	q = btabtmp;
	nbodytmp = 0;
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
		if ( (Pos(p)[0] < 0.5*xbox && Pos(p)[0] > - 0.5*xbox) &&
			 (Pos(p)[1] < 0.5*ybox && Pos(p)[1] > - 0.5*ybox) &&
			 (Pos(p)[2] < 0.5*zbox && Pos(p)[2] > - 0.5*zbox) ) {

			Mass(q) = Mass(p);
			Id(q) = Id(p);
			Type(q) = Type(p);
			SETV(Pos(q), Pos(p));
			SETV(Vel(q), Vel(p));
			++q; ++nbodytmp;

		}
	}

	CLRV(cmpos);
	CLRV(cmvel);
	tmass=0.0;
	DO_BODY(p, btabtmp, btabtmp+nbodytmp) {
		tmass += Mass(p);
		MULVS(tmpv, Pos(p), Mass(p));           
		ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           
        ADDV(cmvel, cmvel, tmpv);
	}
    DIVVS(cmpos, cmpos, tmass);
    DIVVS(cmvel, cmvel, tmass);
	fprintf(stdout,"\n\nNumber of bodies in the simulation box : %d %d\n",nbodytmp,q-btabtmp);
	fprintf(stdout,"CM Pos: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
	fprintf(stdout,"CM Vel: %g %g %g\n",cmvel[0],cmvel[1],cmvel[2]);

	free(bodytab);
	cmd.nbody=nbodytmp;
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
	q=bodytab;
	DO_BODY(p, btabtmp, btabtmp+nbodytmp) {
		Mass(q) = Mass(p);
		Id(q) = Id(p);
		Type(q) = Type(p);
		SUBV(Pos(q), Pos(p), cmpos);
		SUBV(Vel(q), Vel(p), cmvel);
		++q;
	}
	DO_BODY(q, bodytab, bodytab+cmd.nbody) {
		for (k=0; k<3; k++) {
			Pos(q)[k] += 0.0125;
		}
	}

	printf("disk-section finished!\n");
}

local void	Disk_Complete(void)
{
	int Ndi, Ndf, Nd;
	vector Rcmd, Vcmd;
	real Md;

	real Ort_A,Ort_B,disk_Q,kappa,sphi,sigma,sz;
	real mass_i,vphi,vz,vr;
	real sum_0,sum_1x,sum_1y,sum_1z,sum_2x,sum_2y,sum_2z;
	real rtheta,r,z;
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
	bodyptr p;

	Md = 0.1875*cmd.Mtotal;
	Ndi = 0;
	Ndf = cmd.nbody-1;
	Nd = cmd.nbody;
	CLRV(Rcmd);
	CLRV(Vcmd);
	Md = 0.1875;
	Rmaxd = 1.0;

	mass_i = Md / ((real) Nd);
	bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

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
//		z = xrandom(zmin,zmax);								// Uniform z

		Mass(p)	  = mass_i;
		Pos(p)[0] = r*rcos(rtheta);
		Pos(p)[1] = r*rsin(rtheta);
		Pos(p)[2] = z;

		v0 = 1.0;
//		v0 = 0.0;
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
		sum_2x+sum_2y+sum_2z, m2_rhod(AD,Z0,Rmaxd));
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

// TERMINA RUTINA DISK_SECTION -------------------------------------------------


// Rutina para recalcular el CM de dos distribuciones...

local void	two_dist_cofm(void)
{
    bodyptr p;
	vector tmpv;
	
    mtot = 0.0;
	CLRV(cmpos);
    CLRV(cmvel);
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
        mtot += Mass(p);
        MULVS(tmpv, Pos(p), Mass(p));
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));
        ADDV(cmvel, cmvel, tmpv);
    }
    DIVVS(cmpos, cmpos, mtot);
    DIVVS(cmvel, cmvel, mtot);
    for (p = bodytab; p < bodytab+cmd.nbody; p++) { 
		SUBV(Pos(p),Pos(p),cmpos)
		SUBV(Vel(p),Vel(p),cmvel)
    }

}


//----------------- PASAR ESTAS RUTINAS A ANALYSIS_GRAV ------------------------

#define outphasefile	"phs%04d"

local void two_dist_phasesp_out(void)
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

	for ( p = bodytab; p < bodytab+cmd.nbody; p++) {
		prad = rsqrt( Pos(p)[0]*Pos(p)[0]+Pos(p)[1]*Pos(p)[1] );
		pphi=angle(0.0,0.0,Pos(p)[0],Pos(p)[1]);
		vphi = ( Pos(p)[0]*Vel(p)[1] - Pos(p)[1]*Vel(p)[0] )/prad;
		vr = ( Pos(p)[0]*Vel(p)[0] + Pos(p)[1]*Vel(p)[1] )/prad;
		fprintf(outstr,"%14.7E %14.7E %14.7E %14.7E %14.7E %14.7E %14.7E\n",
			prad,pphi,Pos(p)[2],vr,vphi,Vel(p)[2],rabs(vphi));
	}
}

#define outavg2dist	"avg2dist"

local void two_dist_avg()
{
    char namebuf[256];
    struct stat buf;
    stream outstr;

	bodyptr p;
	real prad, vphi, vr;
	real vravg, vphiavg, vzavg, Norm, sigmar, sigmaphi, sigmaz;

    real amabs,amabs1,amabs2,mtot1,mtot2;
    vector tmpv;
#if !defined(THREEDIM)
	real tmps;
#endif
#if defined(THREEDIM)
	vector amvec,amvec1,amvec2,cmpos1,cmvel1,cmpos2,cmvel2;
#else
	real amvec,amvec1,amvec2,cmpos1,cmvel1,cmpos2,cmvel2;
#endif

    sprintf(namebuf, outavg2dist, gd.nstep);        
    if (stat(namebuf, &buf) != 0)               
        outstr = stropen(namebuf, "w");         
    else                                        
        outstr = stropen(namebuf, "a");         

#if defined(THREEDIM)
    CLRV(amvec); 
    CLRV(amvec1); 
    CLRV(amvec2); 
#else
	amvec=0.0;
	amvec1=0.0;
	amvec2=0.0;
#endif    

	Norm = ((real) cmd.nbody);
	vravg = 0.0;
	vphiavg = 0.0;
	vzavg = 0.0;
	sigmar = 0.0;
	sigmaphi = 0.0;
	sigmaz = 0.0;

	for ( p = bodytab; p < bodytab + cmd.nbody; p++) {
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
	for ( p = bodytab; p < bodytab + cmd.nbody; p++) {
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

    mtot1 = 0.0;
	CLRV(cmpos1);
    CLRV(cmvel1);
    for (p = bodytab+Ndisti1; p < bodytab+Ndistf1+1; p++) { 
        mtot1 += Mass(p);
        MULVS(tmpv, Pos(p), Mass(p));
        ADDV(cmpos1, cmpos1, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));
        ADDV(cmvel1, cmvel1, tmpv);
    }
    DIVVS(cmpos1, cmpos1, mtot1);
    DIVVS(cmvel1, cmvel1, mtot1);
    for (p = bodytab+Ndisti1; p < bodytab+Ndistf1+1; p++) { 
		SUBV(Pos(p),Pos(p),cmpos1)
		SUBV(Vel(p),Vel(p),cmvel1)
    }
	for ( p = bodytab+Ndisti1; p < bodytab + Ndistf1+1; p++) {
        CROSSVP(tmpv, Vel(p), Pos(p));          
        MULVS(tmpv, tmpv, Mass(p));
        ADDV(amvec1, amvec1, tmpv);
	}
    ABSV(amabs1, amvec1);                       
	
    mtot2 = 0.0;
	CLRV(cmpos2);
    CLRV(cmvel2);
    for (p = bodytab+Ndisti2; p < bodytab+Ndistf2+1; p++) { 
        mtot2 += Mass(p);
        MULVS(tmpv, Pos(p), Mass(p));
        ADDV(cmpos2, cmpos2, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));
        ADDV(cmvel2, cmvel2, tmpv);
    }
    DIVVS(cmpos2, cmpos2, mtot2);
    DIVVS(cmvel2, cmvel2, mtot2);
    for (p = bodytab+Ndisti2; p < bodytab+Ndistf2+1; p++) { 
		SUBV(Pos(p),Pos(p),cmpos2)
		SUBV(Vel(p),Vel(p),cmvel2)
    }
	for ( p = bodytab+Ndisti2; p < bodytab + Ndistf2+1; p++) {
        CROSSVP(tmpv, Vel(p), Pos(p));          
        MULVS(tmpv, tmpv, Mass(p));
        ADDV(amvec2, amvec2, tmpv);
	}
    ABSV(amabs2, amvec2);                       
	
	out_real_mar(outstr, gd.tnow);
	out_real_mar(outstr, amabs1/mtot1);
	out_real_mar(outstr, amabs2/mtot2);
	out_real_mar(outstr, amabs/(mtot1+mtot2) );
	out_real_mar(outstr, sigmar);
	out_real_mar(outstr, sigmaphi);
	out_real(outstr, sigmaz);
	fclose(outstr);
}

void mass_equatorial_slab(void)
{
	real factor,size,Mtot;
    real dmax,d,delta;
    bodyptr p;
    int i,k,l,ix,iy;
    vector rcm;
    stream outstr;
	real cell[100],super[100][100];
 
	Mtot=0.0;
    CLRV(rcm);
    for (p = bodytab; p < bodytab+cmd.nbody; p++) {
        ADDMULVS(rcm, Pos(p), Mass(p));
		Mtot += Mass(p);
	}
	DIVVS(rcm, rcm, Mtot);

    dmax = 0.0;
	size = 1.0;
    for (p = bodytab; p < bodytab+cmd.nbody; p++)
        for (k = 0; k < NDIM; k++) {
            d = rabs(Pos(p)[k] - rcm[k]);
            if (d > dmax)
                dmax = d;
        }
    while (size < 2 * dmax)
        size = 2 * size;
	factor = 0.1;

	l=30;											

	delta=2*dmax/((real) l-1);
	for (i=0;i<=l;i++)
	    cell[i] = -dmax + ((real) i-1)*delta-0.001*delta;

	printf("dmax, size = %lf %lf\n",dmax,size);

    outstr = stropen("eqslab", "w");         

	for (iy=0; iy<l; iy++)
		for (ix=0; ix<l; ix++)
			super[ix][iy] = 0.0;

	for (p = bodytab; p < bodytab+cmd.nbody; p++) {
		if ( rabs(Pos(p)[2]-rcm[2]) < factor*size ) {
			fprintf(outstr," %lf %lf %lf %lf\n",
				Pos(p)[0],Pos(p)[1],Pos(p)[2],Mass(p));
			ix=0;
			while ( cell[ix+1] <= Pos(p)[0] && ix < l)
				ix++;
			iy=0;
			while ( cell[iy+1] <= Pos(p)[1] && ix < l )
				iy++;

			super[ix][iy] += Mass(p);

		}
	}
	fclose(outstr);

    outstr = stropen("super", "w");         

	for (iy=0; iy<l; iy++){
		for (ix=0; ix<l-1; ix++)
			fprintf(outstr," %14.6e ",super[ix][iy]);
		fprintf(outstr," %14.6e\n",super[l-1][iy]);
	}
	fclose(outstr);
}

// TERMINAN RUTINAS PARA PASAR A ANALYSIS_GRAV ---------------------------------
