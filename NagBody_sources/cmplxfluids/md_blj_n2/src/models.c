/*==============================================================================
	MODULE: models.c				[md_blj_n2]
	Written by: M.A. Rodriguez-Meza
	Starting date: Enero 2005
	Purpose: Initialize md_lj_tree
	Language: C
	Use: 'startrun();'
	Routines and functions: testdata, ReadParameterFile, PrintParameterFile
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: December 2006
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#include "globaldefs.h"

local void testdata_sp_01(void);
local void testdata_sp_02(void);			// To compare with chatos' papers ...
local void testdata_sp_03(void);
local void testdata_sp_04(void);
local void testdata_sp_05(void);
local void testdata_sp_06(void);			// FCC structure default ...
local void MinimumDistanceBBodies(real);

#define NULLMETHOD	0
#define TESTDATA1	1
#define CHATO_IC	2
#define CELLS_IC	3
#define TESTDATA4	4
#define TESTDATA5	5
#define FCC_IC		6						// FCC structure default ...

void testdata(void)
{
    switch(cmd.icModel) {
        case TESTDATA1:
            testdata_sp_01(); break;
        case CHATO_IC:
            testdata_sp_02(); break;		// Chato IC ...
        case CELLS_IC:						// CELLS_IC MODEL ...
            testdata_sp_03(); break;
        case TESTDATA4:
            testdata_sp_04(); break;
        case TESTDATA5:
            testdata_sp_05(); break;
        case FCC_IC:						// FCC structure default ...
            testdata_sp_06(); break;
        case NULLMETHOD:
            printf("\n\tdefault IC ...\n");
            testdata_sp_02(); break;
        default:
            printf("\n\tUnknown IC method...");
            printf("\n\tdefault IC method...\n"); 
            testdata_sp_02(); break;
    }
}

#undef NULLMETHOD
#undef TESTDATA1
#undef CHATO_IC
#undef CELLS_IC
#undef TESTDATA4
#undef TESTDATA5
#undef FCC_IC

local void testdata_sp_01(void)							// CHECK 2D --- OK!!!
{
    vector vcm, b, l;
	bodyptr p;
	real ls, Ekin, tmass, velsq;
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real lmin, lmax;
	vector dr;
    double cpustart;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_01: Parallelepiped : smaller parallelepipeds\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [1] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2,	gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

#if (NDIM==3)
	ls=rpow(gd.nbody,1./3.);
#else
	ls=rsqrt(gd.nbody);
#endif
	DIVVS(l,gdforce.Box,ls);

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g %g",
		l[0],l[1],l[2],ls);
#else
	fprintf(stdout,"\nSimulation box : %g %g",	gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g",
		l[0],l[1],ls);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

#if (NDIM==3)
	b[2]=0.;
	while(b[2]<gdforce.Box[2]) {
#endif
		b[1]=0.;
		while(b[1]<gdforce.Box[1]) {
			b[0]=0.;
			while(b[0]<gdforce.Box[0]) {
				++ibox;
				++i; ++i1;
				if (i1<=gd.nbody1) {
					DO_COORD(k)
						Pos(p)[k] = b[k]+0.001*l[k];
					Id(p) = p-bodytab+1;
					Type(p) = BODY1;
					Mass(p) = gd.mass1;
					DO_COORD(k)
						Vel(p)[k] = grandom(0.0,1.0);
					tmass += Mass(p);
					ADDMULVS(vcm, Vel(p), Mass(p));
					++p;
				} else { 
					--i1; 
					--i;
					++i; 
					++i2;
					if (i2<=gd.nbody2) {
						DO_COORD(k)
							Pos(p)[k] = b[k]+0.001*l[k];
						Id(p) = p-bodytab+1;
						Type(p) = BODY2;
						Mass(p) = gd.mass2;
						DO_COORD(k)
							Vel(p)[k] = grandom(0.0,1.0);
						tmass += Mass(p);
						ADDMULVS(vcm, Vel(p), Mass(p));
						++p;
					} else { 
						--i2; 
						--i; 
					}
				}

				b[0] += l[0];
			}
			b[1] += l[1];
		}
#if (NDIM==3)
		b[2] += l[2];
	}
#endif

	fprintf(gd.outlog,
		"\nTotal number of particles created: %d %d in %d invoqued boxes\n",
		p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin
	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	DO_BODY(p, bodytab, bodytab+gd.nbody)
		VVSAdd(Pos(p), -0.5, gdforce.Box);		// Box center is at (0,0,0)

	Diagnose();

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 01

local void testdata_sp_02(void)						// Chato IC ...
{
    vector vcm, b, l;
	bodyptr p;
	real ls, Ekin, tmass, velsq;
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real lmin, lmax;
	real vol;
    double cpustart;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_02: Parallelepiped : smaller cubic boxes\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [2] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

#if (NDIM==3)
	vol = gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]/gd.nbody;
	ls=rpow(vol,1./3.);
	l[0]=l[1]=l[2]=ls;
#else
	vol = gdforce.Box[0]*gdforce.Box[1]/gd.nbody;
	ls=rsqrt(vol);
	l[0]=l[1]=ls;
#endif

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g %g",
		l[0],l[1],l[2],ls);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g",
		l[0],l[1],ls);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

#if (NDIM==3)
	b[2]=0.;
	while(b[2]<gdforce.Box[2]) {
#endif
		b[1]=0.;
		while(b[1]<gdforce.Box[1]) {
			b[0]=0.;
			while(b[0]<gdforce.Box[0]) {
				++ibox;
				++i; ++i1;
				if (i1<=gd.nbody1) {
					DO_COORD(k) {
						Pos(p)[k] = b[k]+0.5*l[k];
						if(Pos(p)[k]>gdforce.Box[k])
							Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
					}
					Id(p) = p-bodytab+1;
					Type(p) = BODY1;
					Mass(p) = gd.mass1;
					DO_COORD(k)
						Vel(p)[k] = grandom(0.0,1.0);
					tmass += Mass(p);
					ADDMULVS(vcm, Vel(p), Mass(p));
					++p;
				} else { 
					--i1; 
					--i;
					++i; 
					++i2;
					if (i2<=gd.nbody2) {
						DO_COORD(k) {
							Pos(p)[k] = b[k]+0.5*l[k];
							if(Pos(p)[k]>gdforce.Box[k])
								Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
						}
						Id(p) = p-bodytab+1;
						Type(p) = BODY2;
						Mass(p) = gd.mass2;
						DO_COORD(k)
							Vel(p)[k] = grandom(0.0,1.0);
						tmass += Mass(p);
						ADDMULVS(vcm, Vel(p), Mass(p));
						++p;
					} else { 
						--i2; 
						--i; 
					}
				}

				b[0] += l[0];
			}
			b[1] += l[1];
		}
#if (NDIM==3)
		b[2] += l[2];
	}
#endif

	fprintf(gd.outlog,
		"\nTotal number of particles created: %d %d in %d invoqued boxes\n",
		p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	DO_BODY(p, bodytab, bodytab+gd.nbody)
		VVSAdd(Pos(p), -0.5, gdforce.Box);		// Box center is at (0,0,0)

	Diagnose();

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nCHATO_IC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 02

local void testdata_sp_05(void)							// CHECK 2D --- OK!!!
{
    vector vcm, b, l;
	bodyptr p, q;
	real ls, Ekin, tmass, velsq; 
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real dist, distMin, lmin, lmax, drpq, drpq2;
	vector dr;
	int ido;
    double cpustart;
	int iran, ranmax=200;

    cpustart = cputime();                       
	fprintf(gd.outlog,
		"\ntestdata_sp_05: Parallelepiped : Spatially random positions\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [3] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

#if (NDIM==3)
	ls=rpow(gd.nbody,1./3.);
#else
	ls=rsqrt(gd.nbody);
#endif
	DIVVS(l,gdforce.Box,ls);

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
#endif
	lmin *= 1.0;   // Set appropriately. Smaller the scale faster the routine...
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

	for (i=1; i<=gd.nbody; i++) {
		++i1;
		if (i1<=gd.nbody1) {
			iran=0;
			do {
				DO_COORD(k)
					Pos(p)[k] = xrandom(-0.5*gdforce.Box[k],0.5*gdforce.Box[k]);
				distMin = lmax;
				for (q=bodytab; q<p; q++) {
					DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
					DO_COORD(k)
						dr[k]=dr[k] - 
							((real)(nint(dr[k]/gdforce.Box[k])))*gdforce.Box[k];
					DOTVP(drpq2, dr, dr);
					drpq = rsqrt(drpq2);
					distMin = MIN(distMin, drpq);
				}
				++iran; 
				if (iran>ranmax) 
					error("\ntestdata_sp_05: [1] body %d : %s",
						p-bodytab,"Maximum number of iter reached\n");
			} while (distMin < lmin);
			Id(p) = p-bodytab+1;
			Type(p) = BODY1;
			Mass(p) = gd.mass1;
			DO_COORD(k)
				Vel(p)[k] = grandom(0.0,1.0);
			tmass += Mass(p);
			ADDMULVS(vcm, Vel(p), Mass(p));
			++p;
		} else {
			--i1; ++i2;
			if (i2<=gd.nbody2) {
				iran=0;
				do {
					DO_COORD(k)
						Pos(p)[k] = 
							xrandom(-0.5*gdforce.Box[k], 0.5*gdforce.Box[k]);
					distMin = lmax;
					for (q=bodytab; q<p; q++) {
						DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
						DO_COORD(k)
							dr[k]=dr[k] - 
							((real)(nint(dr[k]/gdforce.Box[k])))*gdforce.Box[k];
						DOTVP(drpq2, dr, dr);
						drpq = rsqrt(drpq2);
						distMin = MIN(distMin, drpq);
					}
					++iran;
					if (iran>ranmax)
						error("\ntestdata_sp_05: [2] body %d : %s",
							p-bodytab,"Maximum number of iter reached\n");
				} while (distMin < lmin);
				Id(p) = p-bodytab+1;
				Type(p) = BODY2;
				Mass(p) = gd.mass2;
				DO_COORD(k)
					Vel(p)[k] = grandom(0.0,1.0);
				tmass += Mass(p);
				ADDMULVS(vcm, Vel(p), Mass(p));
				++p;
			} else {
				--i2;
			}
		}
	}

	fprintf(gd.outlog,
		"\nTotal number of particles created: %d %d in %d invoqued boxes\n",
		p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	Diagnose();

#if (NDIM==3)
		fprintf(gd.outlog,"\ndensity = %g %g\n",
			tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
		fprintf(gd.outlog,"\ndensity = %g %g\n",
			tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 05

local void testdata_sp_04(void)							// CHECK 2D --- OK!!!
{
    vector vcm, b, l;
	bodyptr p;
	real ls, Ekin, tmass, velsq, L;
	int k, i, i1, i2;
	int numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, ibox;
	real lmin, lmax;
	real vol;
    double cpustart;

    cpustart = cputime();                       

	fprintf(gd.outlog,"\ntestdata_sp_04: Cubic box : smaller cubic boxes\n");

    strcpy(gd.model_comment, "Simple Cubic Box Model [4] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

	vol = (((real)gd.nbody1)*gd.mass1+((real)gd.nbody2)*gd.mass2)/cmd.density;
#if (NDIM==3)
	L=rpow(vol,1./3.);
	gdforce.Box[0]=gdforce.Box[1]=gdforce.Box[2]=L;
	gd.Lx=gd.Ly=gd.Lz=L;
	ls=L/rpow(gd.nbody,1./3.);
	l[0]=l[1]=l[2]=ls;
#else
	L=rsqrt(vol);
	gdforce.Box[0]=gdforce.Box[1]=L;
	gd.Lx=gd.Ly=L;
	ls=L/rsqrt(gd.nbody);
	l[0]=l[1]=ls;
#endif

#if (NDIM==3)
	fprintf(gd.outlog,"\nBox : %g %g %g",
		gdforce.Box[0],gdforce.Box[1],gdforce.Box[2]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g %g",l[0],l[1],l[2],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	numSmallBoxesZ = gdforce.Box[2]/l[2];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesZ, 
		numSmallBoxesX*numSmallBoxesY*numSmallBoxesZ, 
		(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2])/(l[0]*l[1]*l[2]) );
#else
	fprintf(gd.outlog,"\nBox : %g %g",gdforce.Box[0],gdforce.Box[1]);
	fprintf(gd.outlog,"\nl and ls : %g %g %g",l[0],l[1],ls);
	numSmallBoxesX = gdforce.Box[0]/l[0];
	numSmallBoxesY = gdforce.Box[1]/l[1];
	fprintf(gd.outlog,"\n\nTotal number of small boxes: %d %d %d %g\n", 
		numSmallBoxesX, numSmallBoxesY, numSmallBoxesX*numSmallBoxesY,
				(gdforce.Box[0]*gdforce.Box[1])/(l[0]*l[1]) );
#endif

	Ekin=0.0; tmass=0.0;
    CLRV(vcm);                                  
	CLRV(b);

	fprintf(gd.outlog,"\nCreating bodies: %d %d\n",gd.nbody1, gd.nbody2);
	p = bodytab; i=0; i1=i2=0;
	ibox=0;

	lmin = MIN(l[0],l[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(l[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g %g",
		l[0],l[1],l[2],ls);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nBox per body & scaling length : l and ls : %g %g %g",
		l[0],l[1],ls);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

#if (NDIM==3)
	b[2]=0.;
	while(b[2]<gdforce.Box[2]) {
#endif
		b[1]=0.;
		while(b[1]<gdforce.Box[1]) {
			b[0]=0.;
			while(b[0]<gdforce.Box[0]) {
				++ibox;
				++i; ++i1;
				if (i1<=gd.nbody1) {
					DO_COORD(k) {
						Pos(p)[k] = b[k]+0.5*l[k];
						if(Pos(p)[k]>gdforce.Box[k])
							Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
					}
					Id(p) = p-bodytab+1;
					Type(p) = BODY1;
					Mass(p) = gd.mass1;
					DO_COORD(k)
						Vel(p)[k] = grandom(0.0,1.0);
					tmass += Mass(p);
					ADDMULVS(vcm, Vel(p), Mass(p));
					++p;
				} else { 
					--i1; 
					--i;
					++i; 
					++i2;
					if (i2<=gd.nbody2) {
						DO_COORD(k) {
							Pos(p)[k] = b[k]+0.5*l[k];
							if(Pos(p)[k]>gdforce.Box[k])
								Pos(p)[k] = 0.5*(gdforce.Box[k]+b[k]);
						}
						Id(p) = p-bodytab+1;
						Type(p) = BODY2;
						Mass(p) = gd.mass2;
						DO_COORD(k)
							Vel(p)[k] = grandom(0.0,1.0);
						tmass += Mass(p);
						ADDMULVS(vcm, Vel(p), Mass(p));
						++p;
					} else { 
						--i2; 
						--i; 
					}
				}

				b[0] += l[0];
			}
			b[1] += l[1];
		}
#if (NDIM==3)
		b[2] += l[2];
	}
#endif

	fprintf(gd.outlog,"\n%s: %d %d in %d invoqued boxes\n",
		"Total number of particles created",p-bodytab, i, ibox);
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

	DIVVS(vcm,vcm,tmass);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	DO_BODY(p, bodytab, bodytab+gd.nbody)
		VVSAdd(Pos(p), -0.5, gdforce.Box);		// Box center is at (0,0,0)

	Diagnose();

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

//	MinimumDistanceBBodies(lmax);
    fprintf(stdout,"\nIC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 04

local void testdata_sp_03(void)						// CELLS_IC MODEL ...
{
	vector c, gap, vSum, l, sc, sgap, ctmp;
	vectorI scell;
	int nx, ny, nz, i, i1, i2, si1, si2, j, nxs, nys, nzs;
    double cpustart;
	real tmass, lmin, lmax, ls;
	bodyptr p;
	real velMag;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_03: Simple parallelepiped box : smaller cubic boxes\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [3] created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

	VDiv (gap, gdforce.Box, gd.numUcell);

#if (NDIM==3)
	ls=rpow((gd.nbc1*gd.mass1+gd.nbc2*gd.mass2)/cmd.density, 1./3.);
	l[0]=l[1]=l[2]=ls;
#else							// TENGO DUDA SOBRE EL PAPEL DE GD.NUMUCELL...
	ls=gd.numUcell[0]*rsqrt((gd.nbc1*gd.mass1+gd.nbc2*gd.mass2)/cmd.density);
	l[0]=l[1]=ls;
#endif

	j=1;
	while (gd.nbc1+gd.nbc2 > Cube(j)) ++j;			// ES LO MISMO EN DOS? ...
#if (NDIM==3)
	scell[0]=scell[1]=scell[2]=j;
#else
	scell[0]=scell[1]=j;
#endif

	VDiv (sgap, l, scell);
#if (NDIM==3)
	fprintf(stdout,"\nsGap : %g %g %g", sgap[0], sgap[1], sgap[2]);
#else
	fprintf(stdout,"\nsGap : %g %g",sgap[0],sgap[1]);
#endif

	lmin = MIN(sgap[0],sgap[1]);
	lmax = MAX(gdforce.Box[0],gdforce.Box[1]);
#if (NDIM==3)
	lmin = MIN(sgap[2],lmin);
	lmax = MAX(gdforce.Box[2],lmax);
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nUnit cell : l : %g %g %g",l[0],l[1],l[2]);
	fprintf(stdout,"\nGap : %g %g %g", gap[0], gap[1], gap[2]);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nUnit cell : l : %g %g",l[0],l[1]);
	fprintf(stdout,"\nGap : %g %g",gap[0],gap[1]);
#endif
	fprintf(stdout,"\nlmin, lmax : %g %g", lmin, lmax);

	fprintf(stdout,"\nNumber of bodies in %d smaller cubes : %d %d",
			Cube(j), gd.nbc1, gd.nbc2);

	p = bodytab; i=i1=i2=0;
	tmass = 0.;

#if (NDIM==3)
	for (nz = 0; nz < gd.numUcell[2]; nz++) {
#endif
		for (ny = 0; ny < gd.numUcell[1]; ny++) {
			for (nx = 0; nx < gd.numUcell[0]; nx++) {
#if (NDIM==3)
				VSet (c, nx, ny, nz);
#else
				VSet (c, nx, ny);
#endif
				VMul (c, c, gap);
				si1=si2=0;
#if (NDIM==3)
				for (nzs = 0; nzs < scell[2]; nzs++) {
#endif
					for (nys = 0; nys < scell[1]; nys++) {
						for (nxs = 0; nxs < scell[0]; nxs++) {
#if (NDIM==3)
							VSet (sc, nxs + 0.5, nys + 0.5, nzs + 0.5);
#else
							VSet (sc, nxs + 0.5, nys + 0.5);
#endif
							VMul (sc, sc, sgap);
							VAdd(ctmp, c, sc);
							VVSAdd(ctmp, -0.5, gdforce.Box);
							++i; ++si1;
							if (si1<=gd.nbc1) {
								SETV(Pos(p), ctmp);
								Id(p) = p-bodytab+1;
								Type(p) = BODY1;
								Mass(p) = gd.mass1;
								tmass += Mass(p);
								++p; ++i1;
							} else { 
								--si1; 
								--i;
								++i; 
								++si2;
								if (si2<=gd.nbc2) {
									SETV(Pos(p), ctmp);
									Id(p) = p-bodytab+1;
									Type(p) = BODY2;
									Mass(p) = gd.mass2;
									tmass += Mass(p);
									++p; ++i2;
								} else { 
									--si2; 
									--i; 
								}
							}
						}
					}
#if (NDIM==3)
				}
#endif
			}
		}
#if (NDIM==3)
	}
#endif

	velMag = sqrt (NDIM * (1. - 1. / gd.nbody) * cmd.temperature);
	VZero(vSum);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		VRand(Vel(p));
		VScale(Vel(p), velMag);
		VVAdd(vSum, Vel(p));
	}
	DO_BODY(p, bodytab, bodytab+gd.nbody) 
		VVSAdd (Vel(p), - 1. / gd.nbody, vSum);

	Diagnose();

	fprintf(gd.outlog,"\n%s: %d in %d invoqued boxes\n",
		"Total number of particles created",p-bodytab, 
		VProd(gd.numUcell));
	fprintf(gd.outlog,"bodies types: %d %d\n", i1, i2);

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

    fprintf(stdout,"\nCELLS_IC - cpu time : %g\n\n", cputime() - cpustart);
}  // End icModel 03

local void testdata_sp_06(void)						// FCC structure default ...
{
	vector c, gap, vSum, pos;
	int nx, ny, nz, i, i1, i2, j, i1tot, i2tot;
    double cpustart;
	real tmass, Ekin, velsq;
	bodyptr p;
//	real velMag;

    cpustart = cputime();                       

	fprintf(gd.outlog,
		"\ntestdata_sp_06: Simple parallelepiped box : smaller cubic boxes (FCC)\n");

    strcpy(gd.model_comment, "Simple Parallelepiped Box Model [6] (FCC) created");

	fprintf(stdout,"\nCreating bodies: %d %d %d",
		gd.nbody1, gd.nbody2, gd.nbody);
	fprintf(stdout,"\nBodies per cell: %d %d",gd.nbc1, gd.nbc2);

    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));

	VDiv (gap, gdforce.Box, gd.numUcell);

#if (NDIM==3)
	fprintf(stdout,"\nSimulation box : %g %g %g", 
		gdforce.Box[0], gdforce.Box[1], gdforce.Box[2]);
	fprintf(stdout,"\nUnit cell (Gap) : %g %g %g", gap[0], gap[1], gap[2]);
#else
	fprintf(stdout,"\nSimulation box : %g %g", gdforce.Box[0], gdforce.Box[1]);
	fprintf(stdout,"\nUnit cel (Gap) : %g %g",gap[0],gap[1]);
#endif

	p = bodytab; i1tot=i2tot=0;
	tmass = 0.;

#if (NDIM==3)
	for (nz = 0; nz < gd.numUcell[2]; nz++) {
#endif
		for (ny = 0; ny < gd.numUcell[1]; ny++) {
			for (nx = 0; nx < gd.numUcell[0]; nx++) {
#if (NDIM==3)
				VSet (c, nx+0.25, ny+0.25, nz+0.25);
#else
				VSet (c, nx+0.25, ny+0.25);
#endif
				VMul (c, c, gap);
				VVSAdd (c, -0.5, gdforce.Box);
				i1=i2=1;
#if (NDIM==3)
				for (j = 0; j < 4; j ++) {
#else
				for (j = 0; j < 2; j ++) {
#endif
						SETV(Pos(p), c);
					if (i1<=gd.nbc1) {
#if (NDIM==3)
						if (j != 3) {
							if (j != 0) Pos(p)[0] += 0.5 * gap[0];
							if (j != 1) Pos(p)[1] += 0.5 * gap[1];
							if (j != 2) Pos(p)[2] += 0.5 * gap[2];
						}
#else
						if (j != 1) {
							Pos(p)[0] += 0.5 * gap[0];
							Pos(p)[1] += 0.5 * gap[1];
						}
#endif
						Id(p) = p-bodytab+1;
						Type(p) = BODY1;
						Mass(p) = gd.mass1;
						tmass += Mass(p);
						++p; ++i1; ++i1tot;
					} else {
						if (i2<=gd.nbc2) {
#if (NDIM==3)
							if (j != 3) {
								if (j != 0) Pos(p)[0] += 0.5 * gap[0];
								if (j != 1) Pos(p)[1] += 0.5 * gap[1];
								if (j != 2) Pos(p)[2] += 0.5 * gap[2];
							}
#else
							if (j != 1) {
								Pos(p)[0] += 0.5 * gap[0];
								Pos(p)[1] += 0.5 * gap[1];
							}
#endif
							Id(p) = p-bodytab+1;
							Type(p) = BODY2;
							Mass(p) = gd.mass2;
							tmass += Mass(p);
							++p; ++i2;  ++i2tot;
						}
					}
#if (NDIM==3)
				}
#else
				}
#endif
			}
		}
#if (NDIM==3)
	}
#endif

	fprintf(gd.outlog,"\n%s: %d in %d invoqued boxes\n",
		"Total number of particles created",p-bodytab, 
		VProd(gd.numUcell));
	fprintf(gd.outlog,"bodies types: %d %d\n", i1tot, i2tot);

//	velMag = rsqrt(NDIM*(1. - 1./(real)gd.nbody)*gd.kB*cmd.temperature/gd.mass);
	VZero(vSum);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		VRand(Vel(p));
//		VScale(Vel(p), velMag);
		VVAdd(vSum, Vel(p));
	}
	DO_BODY(p, bodytab, bodytab+gd.nbody) 
		VVSAdd (Vel(p), - 1. / gd.nbody, vSum);

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
//        SUBV(Vel(p), Vel(p), vcm);
        DOTVP(velsq, Vel(p), Vel(p));           
		Ekin+=Mass(p)*velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin before and total mass=%g %g",Ekin,tmass);

	gd.vMag = rsqrt( NDIM*((real)gd.nbody-1.0)*gd.kB*cmd.temperature/2.0 );
	AdjustTemp(Ekin);			// Temperature sets Ekin

	Ekin=0.0;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += Mass(p) * velsq;
    }
	Ekin *= 0.5;
	fprintf(gd.outlog,"\nEkin after=%g\n",Ekin);

	Diagnose();

#if (NDIM==3)
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]*gdforce.Box[2]), cmd.density);
#else
	fprintf(gd.outlog,"\ndensity = %g %g\n",
		tmass/(gdforce.Box[0]*gdforce.Box[1]), cmd.density);
#endif

    fprintf(stdout,"\nFCC-IC - cpu time : %g\n\n", cputime() - cpustart);
} // End icModel 6

local void MinimumDistanceBBodies(real lmax)			// CHECK 2D --- OK!!!
{
	bodyptr p, q;
	int k;
	real distMin, drpq, drpq2;
	vector dr;

	fprintf(gd.outlog,"\nFinding minimum distances between bodies ... ");
	distMin = lmax;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		DO_BODY(q, bodytab, p) {
			DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
			DO_COORD(k)
				dr[k]=dr[k]-((real)(nint(dr[k]/gdforce.Box[k])))*gdforce.Box[k];
			DOTVP(drpq2, dr, dr);
			drpq = rsqrt(drpq2);
			distMin = MIN(distMin, drpq);
		}
	}
	fprintf(gd.outlog,"distMin=%g\n",distMin);
	fflush(gd.outlog);
}

