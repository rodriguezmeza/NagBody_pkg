/*==============================================================================
	MODULE: timestep.c				[md_blj]
	Written by: Mario A. Rodriguez-Meza
	Starting date:	February 2005
	Purpose: compute a timestep evolution of the system
	Language: C
	Use: tree_ljforce(); stepsystem()
	Routines and functions:	tree_ljforce, stepsystem
	External modules, routines and headers:	stdinc, vectmath, mathfns,
		global_defs, proto_defs
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions:
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

local void forcecalc(bodyptr, int);  

local void stepsystem_0(void);
local void stepsystem_1(void);
local void stepsystem_2(void);
local void AdvanceVel(real);
local void BackVel(real);
local void AdvanceVel_KinEnergy(real, realptr);
local void AdvancePos(real);
local void EstimateTimeMin(real);
local void BoundaryCondition(void);
local void AdjustCenterOfMass(void);

void MainLoop(void)										// CHECK 2D --- OK!!!
{
    if (gd.nstep == 0) {                           
        tree_ljforce();
        output();
    }
    if (gd.dtime != 0.0)
        while (cmd.tstop - gd.tnow > 0.01*gd.dtime) {
            stepsystem();
            output();
			checkstop();
			if (gd.stopflag) break;
        }
}

void tree_ljforce(void)									// CHECK 2D --- OK!!!
{
    bodyptr p;

	DO_BODY(p,bodytab,bodytab+gd.nbody)
        Update(p) = TRUE;
	if (gd.forcemethod_int!=4 && gd.forcemethod_int!=5 && gd.forcemethod_int!=6)
		maketree(bodytab, gd.nbody, cmd.options, &gdtree, &gdforce);
	forcecalc(bodytab, gd.nbody);
}

#define FIXEDDT		0
#define FIXEDDT1		1
#define VARIABLEDT	2

void stepsystem(void)									// CHECK 2D --- OK!!!
{
    double cpustart;

    cpustart = cputime();                       

    switch(cmd.intMethod) {
        case FIXEDDT:
            stepsystem_0(); break;
        case FIXEDDT1:
            stepsystem_1(); break;
        case VARIABLEDT:
            stepsystem_2(); break;
        default:
            stepsystem_1(); break;
    }
// Estudiar otra forma de dar esta misma informacion ...
//	fprintf(stdout,"stepsystem CPU additional time: %g\n",
//			cputime()-cpustart-gdforce.cpuforce);
}

#undef FIXEDDT
#undef FIXEDDT1
#undef VARIABLEDT

local void stepsystem_0(void)							// CHECK 2D --- OK!!!
{
    bodyptr p;
    real velsq, Ekin, fEkin, Scale;

	DO_BODY(p,bodytab,bodytab+gd.nbody) {
        ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);     
        ADDMULVS(Pos(p), Vel(p), gd.dtime);
    }

	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		VWrapAll (Pos(p));
	}

	Diagnose();
    tree_ljforce(); 

	Ekin=0.0;
	DO_BODY(p,bodytab,bodytab+gd.nbody) {
        ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += 0.5 * Mass(p) * velsq;
    }

	if (cmd.adjustTemperature) 
		AdjustTemp(Ekin);
	else
		if (gd.nstep < cmd.stepEquil)
			AdjustTemp(Ekin);

	if (cmd.adjustCenterOfMass) AdjustCenterOfMass();

    gd.nstep++;
    gd.nstepNew++;					// In order to restore work ...
    gd.tnow = gd.tnow + gd.dtime;
}

void stepsystem_1(void)									// CHECK 2D --- OK!!!
{
    bodyptr p;
    real velsq, Ekin;
	int k;

	DO_BODY(p,bodytab,bodytab+gd.nbody)
        ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);     

	gd.VelMax=0.;
	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		DOTVP(velsq, Vel(p), Vel(p));
		if (gd.VelMax<velsq) gd.VelMax=velsq;
	}
	gd.VelMax=rsqrt(gd.VelMax);
	gd.TimeMin = gd.LBoxMin/gd.VelMax;
	if (gd.TimeMin < gd.dtime) 
		fprintf(stdout,
			"\nstepsystem : Warning : TimeMin=%g is less than dtime=%g\n",
			gd.TimeMin, gd.dtime);

	DO_BODY(p,bodytab,bodytab+gd.nbody)
        ADDMULVS(Pos(p), Vel(p), gd.dtime);

	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		VWrapAll (Pos(p));
	}

	Diagnose();
    tree_ljforce(); 

	Ekin=0.0;
	DO_BODY(p,bodytab,bodytab+gd.nbody) {
        ADDMULVS(Vel(p), Acc(p), 0.5 * gd.dtime);
        DOTVP(velsq, Vel(p), Vel(p));
        Ekin += 0.5 * Mass(p) * velsq;
    }

//	if (cmd.adjustTemperature) AdjustTemp(Ekin);

	if (cmd.adjustTemperature) 
		AdjustTemp(Ekin);
	else
		if (gd.nstep < cmd.stepEquil)
			AdjustTemp(Ekin);


	if (cmd.adjustCenterOfMass) AdjustCenterOfMass();

    gd.nstep++;
    gd.nstepNew++;					// In order to restore work ...
    gd.tnow = gd.tnow + gd.dtime;
}

local void stepsystem_2(void)							// CHECK 2D --- OK!!!
{
	real dt, dtime, Ekin;
	int i, j=1, iMax=200, jMax=100;
	short flag;
	real accint_time;

	dtime=gd.dtime;
	accint_time=0.0;
fbegin:
	dt=dtime;
	flag=0; i=1;
sbegin:
	AdvanceVel(0.5*dt);
	EstimateTimeMin(dtime);
// Estudiar otra forma de dar esta misma informacion ...
//	fprintf(gd.outlog,"internal timestep iteration %d\n",i);
	if (gd.TimeMin < dt) {
		BackVel(0.5*dt);
		dt /= 2.0;
		i++; flag=1;
		if (i>iMax) error("\n\nstepsystem_2 [1]: %s\n\t\t  %s",
						"Unable to advance the system.","Try a smaller dt");
		goto sbegin;
	}
	AdvancePos(dt);
	BoundaryCondition();
	Diagnose();
	tree_ljforce();
	AdvanceVel_KinEnergy(0.5*dt,&Ekin);

//	if (cmd.adjustTemperature) AdjustTemp(Ekin);
	if (cmd.adjustTemperature) 
		AdjustTemp(Ekin);
	else
		if (gd.nstep < cmd.stepEquil)
			AdjustTemp(Ekin);

	if (cmd.adjustCenterOfMass) AdjustCenterOfMass();

	accint_time+=dt;
// ¡Esta consumiendo mucha espacio en disco para el archivo asociado a "stdout"! 
// Estudiar otra estrategia para dar informacion al respecto pero mas barata ...
//	fprintf(stdout,"accumulated internal time step %g\n",accint_time);

	if (flag) {
		dtime=dtime-dt;
		j++;
		if (j>jMax) error("\n\nstepsystem_2 [2]: %s\n\t\t  %s",
			"Unable to advance the system.","Try a smaller dt\n");
		goto fbegin;
	}
    gd.nstep++;
    gd.nstepNew++;					// In order to restore work ...
    gd.tnow += gd.dtime;

}

local void AdvanceVel(real dt)							// CHECK 2D --- OK!!!
{
    bodyptr p;

	DO_BODY(p,bodytab,bodytab+gd.nbody)
        ADDMULVS(Vel(p), Acc(p), dt);     
}

local void BackVel(real dt)								// CHECK 2D --- OK!!!
{
    bodyptr p;

	DO_BODY(p,bodytab,bodytab+gd.nbody)
        ADDMULVS(Vel(p), Acc(p), -dt);     
}

local void AdvanceVel_KinEnergy(real dt, realptr Ekin)	// CHECK 2D --- OK!!!
{
    bodyptr p;
	real velsq;

	*Ekin=0.0;
	DO_BODY(p,bodytab,bodytab+gd.nbody) {
        ADDMULVS(Vel(p), Acc(p), dt);
        DOTVP(velsq, Vel(p), Vel(p));
        *Ekin += 0.5 * Mass(p) * velsq;
    }
}

local void AdvancePos(real dt)							// CHECK 2D --- OK!!!
{
    bodyptr p;

	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		ADDMULVS(Pos(p), Vel(p), dt);
	}
}

void AdjustTemp(real KinEnergy)							// CHECK 2D --- OK!!!
{
    bodyptr p;
	real vFac;

	gd.kinEnergySave += KinEnergy;
	if (gd.nstep % cmd.stepAdjustTemperature == 0) {
		gd.kinEnergySave /= cmd.stepAdjustTemperature;
//	vFac = gd.vMag/rsqrt(KinEnergy);
		vFac = gd.vMag/rsqrt(gd.kinEnergySave);
		DO_BODY(p, bodytab, bodytab+gd.nbody)
			MULVS(Vel(p), Vel(p), vFac);
		gd.kinEnergySave = 0.0;
	}
}

local void AdjustCenterOfMass(void)						// CHECK 2D --- OK!!!
{
    bodyptr p;
	real mtot;
	vector cmpos, cmvel, tmpv;

printf("\nAdjusting Centor Of Mass ...\n");

    CLRV(cmpos);
    CLRV(cmvel);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
        mtot += Mass(p);                        
        MULVS(tmpv, Pos(p), Mass(p));           
        ADDV(cmpos, cmpos, tmpv);
        MULVS(tmpv, Vel(p), Mass(p));           
        ADDV(cmvel, cmvel, tmpv);
    }
    DIVVS(cmpos, cmpos, mtot);                  
    DIVVS(cmvel, cmvel, mtot);
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		SUBV(Pos(p),Pos(p),cmpos);
//		SUBV(Vel(p),Vel(p),cmvel);
	}
}

local void EstimateTimeMin(real dtime)					// CHECK 2D --- OK!!!
{
    bodyptr p;
	real velsq;

	gd.VelMax=0.;
	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		DOTVP(velsq, Vel(p), Vel(p));
		if (gd.VelMax<velsq) gd.VelMax=velsq;
	}
	gd.VelMax=rsqrt(gd.VelMax);
	gd.TimeMin = gd.LBoxMin/gd.VelMax;
	if (gd.TimeMin < dtime) 
		fprintf(gd.outlog,
			"\nstepsystem : Warning : TimeMin=%g is less than dtime=%g\n",
			gd.TimeMin, dtime);
}

local void BoundaryCondition(void)						// CHECK 2D --- OK!!!
{
    bodyptr p;

	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		VWrapAll (Pos(p));
	}
}

void Diagnose()											// CHECK 2D --- OK!!!
{
    bodyptr p;
	int k;

	DO_BODY(p,bodytab,bodytab+gd.nbody) {
		DO_COORD(k)
			if(Pos(p)[k]<-0.5*gdforce.Box[k] || Pos(p)[k]>0.5*gdforce.Box[k])
				error("\nDiagnose : body=%d is out boundaries\n\n",
					p-bodytab+1);
	}
}


#define BARNES			0
#define NULLMETHOD		1
#define NORMAL			2
#define NBLIST			3
#define DIRECT			4
#define CELLSMETHOD		5
#define DIRECT2			6
#define NORMAL2			7
#define BARNES2			8
#define NORMAL3			9
#define CELLSMETHOD3	10


local void forcecalc(bodyptr btab, int nbody)			// CHECK 2D --- OK!!!
{
    switch(gd.forcemethod_int) {
        case BARNES:
            ljforcecalc_barnes(btab, nbody, &gdforce, &gdtree); break;
        case BARNES2:
            ljforcecalc_barnes2(btab, nbody, &gdforce, &gdtree); break;
        case NULLMETHOD:
            printf("\n\trunning default method (Normal)...\n");
            ljforcecalc_normal(btab, nbody, &gdforce, &gdtree); break;
        case NORMAL:
            ljforcecalc_normal(btab, nbody, &gdforce, &gdtree); break;
        case NORMAL2:
            ljforcecalc_normal2(btab, nbody, &gdforce, &gdtree); break;
        case NORMAL3:
            ljforcecalc_normal3(btab, nbody, &gdforce, &gdtree); break;
        case NBLIST:
            ljforcecalc_nblist(btab, nbody, &gdforce, &gdtree); break;
        case DIRECT:
            ljforcecalc_direct(btab, nbody, &gdforce); break;
        case DIRECT2:
            ljforcecalc_direct2(btab, nbody, &gdforce); break;
        case CELLSMETHOD:
            ljforcecalc_cellsmethod(btab, nbody, &gdforce); break;
        case CELLSMETHOD3:
            ljforcecalc_cellsmethod3(btab, nbody, &gdforce); break;
        default:
            fprintf(stdout,"\n\tforcecalc_method: Unknown method...");
            fprintf(stdout,
				"\n\trunning default force calculation method (Normal)...\n"); 
            ljforcecalc_normal(btab, nbody, &gdforce, &gdtree); break;
    }
}

#undef BARNES
#undef NULLMETHOD
#undef NORMAL
#undef NBLIST
#undef DIRECT 
#undef CELLSMETHOD
#undef DIRECT2
#undef NORMAL2
#undef BARNES2
#undef NORMAL3
#undef CELLSMETHOD3


// -COMIENZAN RUTINAS NOSE-HOOVER THERMOSTAT - FEEDBACK PT ---------------------

/*

void SingleStep ()
{
  ++ stepCount;
  timeNow = stepCount * deltaT;
  PredictorStep ();
  PredictorStepPT ();
  ApplyBoundaryCond ();
  UpdateCellSize ();
  UnscaleCoords ();
  ComputeForces ();
  ComputeDerivsPT ();
  CorrectorStep ();
  CorrectorStepPT ();
  ApplyBoundaryCond ();
  EvalProps ();
  AccumProps (1);
  if (stepCount % stepAvg == 0) {
    AccumProps (2);
    PrintSummary (stdout);
    AccumProps (0);
  }
}

#define PCR4(r, ro, v, a, a1, a2, t)                        \
   r.t = ro.t + deltaT * v.t +                              \
   wr * (cr[0] * a.t + cr[1] * a1.t + cr[2] * a2.t)
#define PCV4(r, ro, v, a, a1, a2, t)                        \
   v.t = (r.t - ro.t) / deltaT +                            \
   wv * (cv[0] * a.t + cv[1] * a1.t + cv[2] * a2.t)

#define PR(t)                                               \
   PCR4 (mol[n].r, mol[n].r, mol[n].rv,                     \
   mol[n].ra, mol[n].ra1, mol[n].ra2, t)
#define PRV(t)                                              \
   PCV4 (mol[n].r, mol[n].ro, mol[n].rv,                    \
   mol[n].ra, mol[n].ra1, mol[n].ra2, t)
#define CR(t)                                               \
   PCR4 (mol[n].r, mol[n].ro, mol[n].rvo,                   \
   mol[n].ra, mol[n].ra1, mol[n].ra2, t)
#define CRV(t)                                              \
   PCV4 (mol[n].r, mol[n].ro, mol[n].rv,                    \
   mol[n].ra, mol[n].ra1, mol[n].ra2, t)

void PredictorStep ()
{
  real cr[] = {19.,-10.,3.}, cv[] = {27.,-22.,7.}, div = 24., wr, wv;
  int n;

  wr = Sqr (deltaT) / div;
  wv = deltaT / div;
  DO_MOL {
    mol[n].ro = mol[n].r;
    mol[n].rvo = mol[n].rv;
    PR (x);
    PRV (x);
    PR (y);
    PRV (y);
    PR (z);
    PRV (z);
    mol[n].ra2 = mol[n].ra1;
    mol[n].ra1 = mol[n].ra;
  }
}

void CorrectorStep ()
{
  real cr[] = {3.,10.,-1.}, cv[] = {7.,6.,-1.}, div = 24., wr, wv;
  int n;

  wr = Sqr (deltaT) / div;
  wv = deltaT / div;
  DO_MOL {
    CR (x);
    CRV (x);
    CR (y);
    CRV (y);
    CR (z);
    CRV (z);
  }
}


#undef PCR4
#undef PCV4

#define PCR4(r, ro, v, a, a1, a2)                           \
   r = ro + deltaT * v +                                    \
   wr * (cr[0] * a + cr[1] * a1 + cr[2] * a2)
#define PCV4(r, ro, v, a, a1, a2)                           \
   v = (r - ro) / deltaT +                                  \
   wv * (cv[0] * a + cv[1] * a1 + cv[2] * a2)

void PredictorStepPT ()
{
  real cr[] = {19.,-10.,3.}, cv[] = {27.,-22.,7.}, div = 24., e, wr, wv;

  wr = Sqr (deltaT) / div;
  wv = deltaT / div;
  varSo = varS;
  varSvo = varSv;
  varVo = varV;
  varVvo = varVv;
  PCR4 (varS, varS, varSv, varSa, varSa1, varSa2);
  PCV4 (varS, varSo, varSv, varSa, varSa1, varSa2);
  PCR4 (varV, varV, varVv, varVa, varVa1, varVa2);
  PCV4 (varV, varVo, varVv, varVa, varVa1, varVa2);
  varSa2 = varSa1;
  varVa2 = varVa1;
  varSa1 = varSa;
  varVa1 = varVa;
  e = pow (varV, 1. / NDIM);
  VSetAll (region, e);
}

void CorrectorStepPT ()
{
  real cr[] = {3.,10.,-1.}, cv[] = {7.,6.,-1.}, div = 24., e, wr, wv;

  wr = Sqr (deltaT) / div;
  wv = deltaT / div;
  PCR4 (varS, varSo, varSvo, varSa, varSa1, varSa2);
  PCV4 (varS, varSo, varSvo, varSa, varSa1, varSa2);
  PCR4 (varV, varVo, varVvo, varVa, varVa1, varVa2);
  PCV4 (varV, varVo, varVvo, varVa, varVa1, varVa2);
  e = pow (varV, 1. / NDIM);
  VSetAll (region, e);
}


#undef VWrap

#define VWrap(v, t)                                         \
   if (v.t >= 0.5) v.t -= 1.;                               \
   else if (v.t < -0.5) v.t += 1.

void ApplyBoundaryCond ()
{
  int n;

  DO_MOL VWrapAll (mol[n].r);
}

void ComputeDerivsPT ()
{
  real aFac, vFac;
  int n;

  vvSum = 0.;
  DO_MOL vvSum += VLenSq (mol[n].rv);
  vvSum *= pow (varV, 2./3.);
  g1Sum = vvSum - 3. * nMol * temperature;
  g2Sum = vvSum + virSum - 3. * extPressure * varV;
  aFac = pow (varV, -1./3.);
  vFac = - varSv / varS - 2. * varVv / (3. * varV);
  DO_MOL VSSAdd (mol[n].ra, aFac, mol[n].ra, vFac, mol[n].rv);
  varSa = Sqr (varSv) / varS + g1Sum * varS / massS;
  varVa = varSv * varVv / varS +
     g2Sum * Sqr (varS) / (3. * massV * varV);
}

void InitFeedbackVars ()
{
  varS = 1.;
  varSv = 0.;
  varSa = varSa1 = varSa2 = 0.;
  varV = Cube (region.x);
  varVv = 0.;
  varVa = varVa1 = varVa2 = 0.;
}

void ScaleCoords ()
{
  real fac;
  int n;

  fac = pow (varV, -1. / 3.);
  DO_MOL VScale (mol[n].r, fac);
}

void UnscaleCoords ()
{
  real fac;
  int n;

  fac = pow (varV, 1. / 3.);
  DO_MOL VScale (mol[n].r, fac);
}

void ScaleVels ()
{
  real fac;
  int n;

  fac = pow (varV, -1. / 3.);
  DO_MOL VScale (mol[n].rv, fac);
}

void UpdateCellSize ()
{
  VSCopy (cells, 1. / rCut, region);
  cells.x = Min (cells.x, maxEdgeCells);
  cells.y = Min (cells.y, maxEdgeCells);
  cells.z = Min (cells.z, maxEdgeCells);
}

*/

// -TERMINAN RUTINAS NOSE-HOOVER THERMOSTAT - FEEDBACK PT ----------------------
