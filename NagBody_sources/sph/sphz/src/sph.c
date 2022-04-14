//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)


#include "../../../General_libs/zeno/clib/stdinc.h"
#include "../../../General_libs/zeno/clib/mathfns.h"
#include "../../../General_libs/zeno/clib/vectmath.h"
#include "../../../General_libs/zeno/clib/getparam.h"
#include "../../../General_libs/zeno/clib/datatypes.h"
#include "../../../General_libs/zeno/clib/gsp.h"

#include "sphcode.h"
#include "kdtree.h"
#include "smooth.h"
#include "../../../General_libs/zeno/libnbody/fixbody.h"

//#include "globaldefs.h"
#include "switches.h"
#include "protodefs.h"

local void setupbody(void);
local void startnewrun(void);
local void startoldrun(void);
local void microstep(real *);
local void update(int, real);
local real coolfunc(bodyptr);
local void evaluate(bool);
local void gravcalc(bool);
local void sphcalc(bool);
local void sum_density(smxptr, int, int);
local void est_tau(smxptr, int, int);
local void sum_derivs(smxptr, int, int);
local void setlevels(smxptr);
local void starform(real);


void startrun(void)
{
    infile = getparam("in");
    outfile = getparam("out");
    savefile = getparam("save");
    restfile = getparam("restore");
    setupbody();
    if (! strnull(infile))
	startnewrun();
    else if (! strnull(restfile))
	startoldrun();
    else
	error("%s: must supply file to input or restore\n", getargv0());
    if (scanopt(options, "new-tout"))
	tout = tnow + dtout;
}


local void setupbody(void)
{
    define_body(sizeof(body), Precision, NDIM);
    define_body_offset(TypeTag, BodyOffset(Type));
    define_body_offset(PosTag, BodyOffset(Pos));
    define_body_offset(VelTag, BodyOffset(Vel));
    define_body_offset(MassTag, BodyOffset(Mass));
    define_body_offset(SmoothTag, BodyOffset(Smooth));
    define_body_offset(PhiTag, BodyOffset(Phi));
    define_body_offset(AccTag, BodyOffset(Acc));
    define_body_offset(RhoTag, BodyOffset(Rho));
#if defined(ENTROPY)
    define_body_offset(EntFuncTag, BodyOffset(EntFunc));
#else
    define_body_offset(UinternTag, BodyOffset(Uintern));
#endif
    define_body_offset(UdotIntTag, BodyOffset(UdotInt));
#if defined(RADIATING)
    define_body_offset(UdotRadTag, BodyOffset(UdotRad));
#endif
#if defined(COMPVISC)
    define_body_offset(UdotVisTag, BodyOffset(UdotVis));
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    define_body_offset(TauTag, BodyOffset(Tau));
#endif
#if defined(STARFORM)
    define_body_offset(BirthTag, BodyOffset(Birth));
#endif
}


local void startnewrun(void)
{
    real dt1, dt2;

    gamma0 = getdparam("gamma");
#if defined(RADIATING)
    uradmax = getdparam("uradmax");
    lambmax = getdparam("lambmax");
#endif
#if defined(DIFFUSING)
    sigmastar = getdparam("sigmastar");
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    opacity = getdparam("opacity");
#endif
#if defined(CONDUCTING)
    conduct = getdparam("conduct");
#endif
#if defined(STARFORM)
    starprob = getdparam("starprob");
    rhoindx = getdparam("rhoindx");
    udotindx = getdparam("udotindx");
    (void) initstate(getiparam("starseed"), randstate, RANDSIZE);
#endif
    alpha = getdparam("alpha");
    beta = getdparam("beta");
    nsmooth = getiparam("nsmooth");
    nbucket = getiparam("nbucket");
    slope0 = getdparam("slope");
    courant = getdparam("courant");
    dtime = (sscanf(getparam("dtime"), "%f/%f", &dt1, &dt2) == 2 ?
	     dt1 / dt2 : getdparam("dtime"));
    fdrag = getdparam("fdrag");
#if defined(GRAVITY)
    eps = getdparam("eps");
    usequad = getbparam("usequad");
#if !defined(QUICKSCAN)
    theta = getdparam("theta");
#endif
#elif defined(EXTGRAV)
    gspfile = getparam("gravgsp");
#endif
    options = getparam("options");
    outputs = getparam("outputs");
    tstop = getdparam("tstop");
    dtout = (sscanf(getparam("dtout"), "%f/%f", &dt1, &dt2) == 2 ?
	     dt1 / dt2 : getdparam("dtout"));
    inputdata();
    tout = tnow;
    rsize = 1.0;
    nstep = 0;
    eradiate = 0.0;
}


local void startoldrun(void)
{
    real dt1, dt2;

    restorestate(restfile);
    if (getparamstat("gamma") & ARGPARAM)
	gamma0 = getdparam("gamma");
#if defined(RADIATING)
    if (getparamstat("uradmax") & ARGPARAM)
	uradmax = getdparam("uradmax");
    if (getparamstat("lambmax") & ARGPARAM)
	lambmax = getdparam("lambmax");
#endif
#if defined(DIFFUSING)
    if (getparamstat("sigmastar") & ARGPARAM)
	sigmastar = getdparam("sigmastar");
#endif    
#if defined(DIFFUSING) || defined(OPAQUE)
    if (getparamstat("opacity") & ARGPARAM)
	opacity = getdparam("opacity");
#endif
#if defined(CONDUCTING)
    if (getparamstat("conduct") & ARGPARAM)
	conduct = getdparam("conduct");
#endif
#if defined(STARFORM)
    if (getparamstat("starprob") & ARGPARAM)
	starprob = getdparam("starprob");
    if (getparamstat("rhoindx") & ARGPARAM)
	rhoindx = getdparam("rhoindx");
    if (getparamstat("udotindx") & ARGPARAM)
	udotindx = getdparam("udotindx");
    if (getparamstat("starseed") & ARGPARAM)
        (void) initstate(getiparam("starseed"), randstate, RANDSIZE);
#endif
    if (getparamstat("alpha") & ARGPARAM)
	alpha = getdparam("alpha");
    if (getparamstat("beta") & ARGPARAM)
	beta = getdparam("beta");
    if (getparamstat("nsmooth") & ARGPARAM)
	nsmooth = getiparam("nsmooth");
    if (getparamstat("nbucket") & ARGPARAM)
	nbucket = getiparam("nbucket");
    if (getparamstat("slope") & ARGPARAM)
	slope0 = getdparam("slope");
    if (getparamstat("courant") & ARGPARAM)
	courant = getdparam("courant");
    if (getparamstat("fdrag") & ARGPARAM)
	fdrag = getdparam("fdrag");
#if defined(GRAVITY)
    if (getparamstat("eps") & ARGPARAM)
	eps = getdparam("eps");
    if (getparamstat("usequad") & ARGPARAM)
	usequad = getbparam("usequad");
#if !defined(QUICKSCAN)
    if (getparamstat("theta") & ARGPARAM)
	theta = getdparam("theta");
#endif
#endif

    if (getparamstat("options") & ARGPARAM)
	options = getparam("options");
    if (getparamstat("outputs") & ARGPARAM)
	outputs = getparam("outputs");
    if (getparamstat("tstop") & ARGPARAM)
	tstop = getdparam("tstop");
    if (getparamstat("dtout") & ARGPARAM)
	dtout = (sscanf(getparam("dtout"), "%f/%f", &dt1, &dt2) == 2 ?
		 dt1 / dt2 : getdparam("dtout"));
}


void initstep(void)
{
    real ftmp, dt;
    bodyptr p;

    for (p = btab; p < btab+nbody; p++)	{
        Frequency(p) = 1.0 / dtime;
	Flags(p) = INCLUDE | UPDATE;
    }
    evaluate(TRUE);
    levmax = 0;
    for (p = btab; p < btab+nbody; p++) {
	CurLevel(p) = NewLevel(p);
	levmax = MAX(levmax, CurLevel(p));
	dt = dtime / (2 << CurLevel(p));
	SETV(Vmid(p), Vel(p));
	ADDMULVS(Vmid(p), Acc(p), dt);
    }
}


void macrostep(void)
{
    real fstep;
    int lev, n, m;

#if defined(TRACESTEP)
    fprintf(stdout, "\n    %12s%8s%8s\n", "fstep", "level", "update");
#endif
    fstep = 0.0;
    microstep(&fstep);
    while (fstep < 1.0) {
	lev = levmax;
	n = (1 << levmax) * fstep;
	for (m = n ^ (n-1); m > 1; m = m >> 1)
	    lev--;
	update(lev, fstep);
	microstep(&fstep);
    }
    if (fstep > 1.0)
	error("%s: fstep = %f\n", getargv0(), fstep);
    nstep++;
    update(0, fstep);
#if defined(STARFORM)
    if (!scanopt(options, "micro-star"))
        starform(dtime);
#endif
}


local void microstep(real *fstep)
{
    real dt, edotrad, udotrad;
    bodyptr p;
    vector vtmp;

    dt = dtime / (1 << levmax);
    edotrad = 0.0;
    for (p = btab; p < btab+nbody; p++) {
	SETV(vtmp, Vel(p));
	ADDMULVS(vtmp, Acc(p), 0.5 * dt);
	ADDMULVS(Pos(p), vtmp, dt);
	ADDMULVS(Vel(p), Acc(p), dt);
	if (Gas(p)) {
#if defined(ENTROPY)
#  if defined(ADIABATIC)
	    udotrad = - UdotInt(p);
#  elif defined(RADIATING)
	    udotrad = UdotRad(p);
#  else
	    udotrad = 0.0;
#  endif
	    EntFunc(p) += dt * (UdotInt(p) + udotrad) *
			    (gamma0 - 1) / rpow(Rho(p), gamma0 - 1);
	    if (EntFunc(p) < 0.0)
		error("%s: negative entropy function\n", getargv0());
#else
#  if defined(ISOTHERMAL)
	    udotrad = - UdotInt(p);
#  elif defined(RADIATING)
	    udotrad = UdotRad(p);
#  else
	    udotrad = 0.0;
#  endif
	    Uintern(p) += dt * (UdotInt(p) + udotrad);
	    if (Uintern(p) < 0.0)
		error("%s: negative internal energy\n", getargv0());
#endif
	    edotrad -= Mass(p) * udotrad;
	}
    }
    eradiate += dt * edotrad;
    tnow += dt;
#if defined(STARFORM)
    if (scanopt(options, "micro-star"))
        starform(dt);
#endif
    *fstep += 1.0 / (1 << levmax);
}


local void update(int lev, real fstep)
{
    int nupdate, oldlev, newlev;
    bodyptr p;
    real dt1, dt2, x;
    vector vdif, aold;

    nupdate = 0;
    for (p = btab; p < btab+nbody; p++) {
	Flags(p) = INCLUDE | (CurLevel(p) < lev ? 0 : UPDATE);

	if (Update(p))
	    nupdate++;
    }
#if defined(TRACESTEP)
    fprintf(stdout, "    %12.5f%8d%8d\n", fstep, lev, nupdate);
#endif
    fflush(stdout);
    evaluate(fstep == 1.0);
    levmax = lev;
    for (p = btab; p < btab+nbody; p++)
	if (Update(p)) {
	    oldlev = CurLevel(p);
	    newlev = MAX(NewLevel(p), MAX(oldlev - 1, lev));

	    CurLevel(p) = newlev;
	    levmax = MAX(levmax, CurLevel(p));
	    dt1 = dtime / (2 << oldlev);
	    dt2 = dtime / (2 << newlev);
	    SUBV(vdif, Vel(p), Vmid(p));
	    DIVVS(aold, vdif, dt1);
	    SETV(Vel(p), Vmid(p));
	    ADDMULVS(Vel(p), Acc(p), 0.75 * dt1);
	    ADDMULVS(Vel(p), aold, 0.25 * dt1);
	    x = (3 + dt2 / dt1) / 4;
	    ADDMULVS(Vmid(p), Acc(p), x * (dt1 + dt2));
	    ADDMULVS(Vmid(p), aold, (1 - x) * (dt1 + dt2));
	}
}


local void evaluate(bool report)
{
    bodyptr p;

    if (report)
	outputhead(stdout);
    gravcalc(report);
    sphcalc(report);
    if (fdrag != 0.0)
        for (p = btab; p < btab+nbody; p++)
            if (Gas(p) && Update(p)) {
		ADDMULVS(Acc(p), Vel(p), - fdrag);

	}
}


local void gravcalc(bool report)
{
    bodyptr p;
    real r, mrinv3;

#if defined(GRAVITY)
    maketree(btab, nbody);
    gravforce();
    if (report)
	report_force(stdout, nbody);
#elif defined(EXTGRAV)
    for (p = btab; p < btab+nbody; p++)
	if (Update(p) && gravgsp != NULL) {
	    r = absv(Pos(p));
	    mrinv3 = - mass_gsp(gravgsp, r) / (r*r*r);
	    MULVS(Acc(p), Pos(p), mrinv3);
	    Phi(p) = phi_gsp(gravgsp, r);
	} else if (Update(p)) {
	    CLRV(Acc(p));
	    Phi(p) = 0.0;
	}
#elif defined(NOACCEL)
    for (p = btab; p < btab+nbody; p++)
	if (Update(p)) {
	    CLRV(Acc(p));
	    Phi(p) = 0.0;
	}
#endif
}


local void sphcalc(bool report)
{
    kdxptr kd;
    smxptr sm;
    bodyptr p;
    real freqrad;

    kd = init_kdtree(btab, nbody, ngas);
    build_kdtree(kd, nbucket);
    sm = init_smooth(kd, nsmooth, slope0);
    for (p = btab; p < btab+nbody; p++)
        if (Gas(p)) {
	    Rho(p) = 0.0;
	    UdotInt(p) = 0.0;
#if defined(RADIATING)
	    UdotRad(p) = 0.0;
#endif
#if defined(COMPVISC)
	    UdotVis(p) = 0.0;
#endif
	}
    smooth(sm, sum_density);
#if defined(DIFFUSING) || defined(OPAQUE)
    resmooth(sm, est_tau);
#endif
    for (p = btab; p < btab+nbody; p++)
	if (Gas(p)) {
#if defined(ENTROPY)
	    Press(p) = EntFunc(p) * rpow(Rho(p), gamma0);

#else
	    Press(p) = (gamma0 - 1) * Uintern(p) * Rho(p);

#endif
	}
    resmooth(sm, sum_derivs);
#if defined(RADIATING)
    for (p = btab; p < btab+nbody; p++)
	if (Gas(p)) {
#if defined(DIFFUSING) || defined(OPAQUE)
	    UdotRad(p) -= (Tau(p) < 1.0 ? coolfunc(p) : 0.0);
#if defined(DIFFUSING)
	    if (Surface(p))
		UdotRad(p) -= sigmastar * rsqr(rsqr(Uintern(p))) /
			        (2 * Smooth(p) * Rho(p));
#endif
#else
	    UdotRad(p) -= coolfunc(p);
#endif
	    freqrad = (ABS(UdotRad(p)) / Uintern(p)) / courant;
	    Frequency(p) = MAX(Frequency(p), freqrad);

	}
#endif
    setlevels(sm);
    if (report)
	report_smooth(sm, stdout, options);
    finish_smooth(sm);
    finish_kdtree(kd);
}


local void sum_density(smxptr sm, int pi, int nball)
{
    bodyptr bi = sm->kd->bptr[pi], bj;
    real hinv2, wsc, rhinv2, wsm;
    int j;

    hinv2 = 4 / sm->r2ball[pi];
    wsc = 0.5 * rsqrt(hinv2) * hinv2 / PI;
    for (j = 0; j < nball; ++j) {
        bj = sm->kd->bptr[sm->inlist[j]];
        rhinv2 = sm->r2list[j] * hinv2;
        WSmooth(wsm, wsc, rhinv2, sm->coefs);
        Rho(bi) += wsm * Mass(bj);
        Rho(bj) += wsm * Mass(bi);
    }
}

#if defined(DIFFUSING) || defined(OPAQUE)


local void est_tau(smxptr sm, int pi, int nball)
{
    bodyptr bi = sm->kd->bptr[pi], bj;
    real hinv2, dwsc, rhinv2, dwsm, dwsmrinv, fij;
    vector gradrho, rij, gradW;
    int j;

    hinv2 = 4 / sm->r2ball[pi];
    dwsc = hinv2 * hinv2 / PI;
    CLRV(gradrho);
    for (j = 0; j < nball; ++j) {
        bj = sm->kd->bptr[sm->inlist[j]];
	if (bi != bj) {
            rhinv2 = sm->r2list[j] * hinv2;
	    dWSmooth(dwsm, dwsc, rhinv2, sm->coefs);
	    dwsmrinv = dwsm / sqrt(sm->r2list[j]);
	    SUBV(rij, Pos(bi), Pos(bj));
	    fij = Rho(bj) / Rho(bi) - 1.0;
	    ADDMULVS(gradrho, rij, Mass(bj) * fij * dwsmrinv);
	}
    }
    Tau(bi) = opacity * rsqr(Rho(bi)) / absv(gradrho);
}

#endif


local void sum_derivs(smxptr sm, int pi, int nball)
{
    bodyptr bi = sm->kd->bptr[pi], bj;
    real hinv2, dwsc, epsi, hsi, csi, divv, mumax, rhinv2, dwsm, dwsmrinv;
    real epsj, rv, hsavg, csavg, rhoavg, muvis, pivis, ascale, freqcou;
    real chii, chij, chiavg, *r2list = sm->r2list;
    int *inlist = sm->inlist, j, lev;
    vector rij, vij;

    hinv2 = 4 / sm->r2ball[pi];
    dwsc = 0.5 * hinv2 * hinv2 / PI;
    epsi = Press(bi) / (Rho(bi) * Rho(bi));
    hsi = rsqrt(sm->r2ball[pi]) / 2;
    csi = rsqrt(gamma0 * Press(bi) / Rho(bi));
    divv = mumax = 0.0;
    for (j = 0; j < nball; j++) {
	bj = sm->kd->bptr[inlist[j]];
        if (bi != bj) {
            rhinv2 = r2list[j] * hinv2;
            dWSmooth(dwsm, dwsc, rhinv2, sm->coefs);
            dwsmrinv = dwsm / rsqrt(r2list[j]);
            epsj = Press(bj) / (Rho(bj) * Rho(bj));
            SUBV(rij, Pos(bi), Pos(bj));
            SUBV(vij, Vel(bi), Vel(bj));
            DOTVP(rv, rij, vij);
	    hsavg = (hsi + rsqrt(sm->r2ball[inlist[j]]) / 2) / 2;
	    muvis = rv / (r2list[j] / hsavg + ETA2 * hsavg);
	    mumax = MAX(mumax, ABS(muvis));
            if (rv < 0) {
                csavg = (csi + rsqrt(gamma0 * Press(bj) / Rho(bj))) / 2;
                rhoavg = (Rho(bi) + Rho(bj)) / 2;
                pivis = muvis * (beta * muvis - alpha * csavg) / rhoavg;
#if defined(ENTROPY)
		UdotInt(bi) += Mass(bj) * 0.5 * pivis * rv * dwsmrinv;
		UdotInt(bj) += Mass(bi) * 0.5 * pivis * rv * dwsmrinv;
#endif
            } else
                pivis = 0.0;
#if !defined(NOACCEL)
            ascale = (epsi + epsj + pivis) * dwsmrinv;
	    if (Update(bi)) {
		ADDMULVS(Acc(bi), rij, - Mass(bj) * ascale);
	    }
	    if (Update(bj)) {
		ADDMULVS(Acc(bj), rij, Mass(bi) * ascale);
	    }
#endif

#if !defined(ENTROPY)
	    UdotInt(bi) += Mass(bj) * (epsi + 0.5 * pivis) * rv * dwsmrinv;
	    UdotInt(bj) += Mass(bi) * (epsj + 0.5 * pivis) * rv * dwsmrinv;
#  if defined(COMPVISC)
	    UdotVis(bi) += Mass(bj) * 0.5 * pivis * rv * dwsmrinv;
	    UdotVis(bj) += Mass(bi) * 0.5 * pivis * rv * dwsmrinv;
#  endif
#  if defined(RADIATING) && defined(DIFFUSING)
	    chii = (16.0/3.0) * (sigmastar / opacity) *
		     (Tau(bi) > 1.0 ? rqbe(Uintern(bi)) / Rho(bi) : 0.0);
	    chij = (16.0/3.0) * (sigmastar / opacity) *
		     (Tau(bj) > 1.0 ? rqbe(Uintern(bj)) / Rho(bj) : 0.0);
	    chiavg = 0.5 * (chii + chij);
	    UdotInt(bi) += Mass(bj) * (chiavg / (Rho(bi) * Rho(bj))) *
			     (Uintern(bi) - Uintern(bj)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
	    UdotInt(bj) += Mass(bi) * (chiavg / (Rho(bi) * Rho(bj))) *
			     (Uintern(bj) - Uintern(bi)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
	    if ((Tau(bi) > 1.0) != (Tau(bj) > 1.0))
		SetFlag(bi, SURFACE);
#  endif
#  if defined(CONDUCTING)
	    UdotInt(bi) += Mass(bj) * (conduct / (Rho(bi) * Rho(bj))) *
			     (Uintern(bi) - Uintern(bj)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
	    UdotInt(bj) += Mass(bi) * (conduct / (Rho(bi) * Rho(bj))) *
			     (Uintern(bj) - Uintern(bi)) * dwsm *
			       rsqrt(r2list[j]) / (r2list[j] + ETA2 * hsavg);
#  endif
#endif
            divv -= 2 * Mass(bj) * rv * dwsmrinv / Rho(bi);
        }
    }
    if (Update(bi))
	Frequency(bi) = ((csi + 1.2 * (alpha*csi + beta*mumax)) / hsi +
			   ABS(divv)) / courant;
}

#if defined(RADIATING)


local real coolfunc(bodyptr p)
{
    real logx, f1, f2, f3, udotrad;

    logx = log10(Uintern(p) / uradmax);
    f1 = -0.17 * rsqr(logx + 1.23) - 3.14;
    f2 = -1.88 * rsqr(rsqr(logx));
    f3 = (logx > 0.97 ? 0.0 : -2.0 * rsqr(rsqr(logx - 0.97)))  - 1.60;
    udotrad = lambmax * Rho(p) * (rdex(f1) + rdex(f2) + rdex(f3));
    return (udotrad);
}    

#endif


local void setlevels(smxptr sm)
{
    bodyptr p;
    real freqmin;
    int lev, nmax;

    sm->freqmax = sm->freqsum = 0.0;
    for (p = btab; p < btab+nbody; p++) {
        if (Gas(p) && Update(p)) {
#if defined(CONDUCTING)
	    freqmin = (ABS(UdotInt(p)) / Uintern(p)) / courant;
	    if (Frequency(p) < freqmin)
	        Frequency(p) = freqmin;
#endif
	    sm->freqmax = MAX(sm->freqmax, Frequency(p));
	    sm->freqsum += Frequency(p);
	    lev = rceil(rlog2(MAX(dtime * Frequency(p), 1.0)));
	    if (lev > MAXLEV) {
	        eprintf("[%s: courant condition violated; level = %d]\n",
		        getargv0(), lev);
		lev = MAXLEV;
	    }
	    NewLevel(p) = lev;
	}
	if (! Gas(p))
	    NewLevel(p) = 0;
    }
    if (! scanopt(options, "nolimit")) {
	lev = 0;
	for (p = btab; p < btab+nbody; p++)
	    if (NewLevel(p) > lev) {
		lev = NewLevel(p);
		nmax = 1;
	    } else if (NewLevel(p) == lev)
		nmax++;
	if (nmax < nsmooth/4)
	    for (p = btab; p < btab+nbody; p++)
		if (NewLevel(p) == lev)
		    NewLevel(p)--;
    }
    if (scanopt(options, "fixstep")) {
	if (scanopt(options, "lockstep"))
	    error("%s: fixstep and lockstep incompatable\n", getargv0());
	for (p = btab; p < btab+nbody; p++)
	    NewLevel(p) = 0;
    }
    if (scanopt(options, "lockstep")) {
	lev = 0;
	for (p = btab; p < btab+nbody; p++)
	    lev = MAX(NewLevel(p), lev);
	for (p = btab; p < btab+nbody; p++)
	    NewLevel(p) = lev;
    }
}

#if defined(STARFORM)


local void starform(real dt)
{
    bodyptr p;
    real frac;

    for (p = btab; p < btab+nbody; p++)
	if (Gas(p)) {
#if defined(COMPVISC)
	    frac = starprob * rpow(Rho(p), rhoindx) *
		     rpow(UdotVis(p), udotindx) * dt;
#else	
	    frac = starprob * rpow(Rho(p), rhoindx) *
		     rpow(MAX(UdotInt(p), 0.0), udotindx) * dt;
#endif
	    if (xrandom(0.0, 1.0) < frac) {
	        Type(p) = BODY | STAR;
		eradiate += Mass(p) * Uintern(p);
		Birth(p) = tnow;
		ngas--;
	    }
	}
}

#endif
