//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)


#include "../../../General_libs/zeno/clib/stdinc.h"
#include "../../../General_libs/zeno/clib/mathfns.h"
#include "../../../General_libs/zeno/clib/vectmath.h"
#include "../../../General_libs/zeno/clib/getparam.h"
#include "../../../General_libs/zeno/clib/strset.h"
#include "../../../General_libs/zeno/clib/filestruct.h"
#include "../../../General_libs/zeno/clib/gsp.h"
#include "sphcode.h"
#include "../../../General_libs/zeno/libnbody/fixbody.h"

#include <sys/types.h>
#include <sys/stat.h>


local void diagnostics(void);


local real mtot;
local real etot[4];
local matrix keten;
local matrix peten;
local vector cmpos;
local vector cmvel;
local vector amvec;


void inputdata(void)
{
    stream istr, gstr;
    string intags[MaxBodyFields];
    bodyptr p;

    btab = NULL;
    istr = stropen(infile, "r");
    get_history(istr);
    if (! get_snap(istr, &btab, &nbody, &tnow, intags, FALSE))
	error("inputdata: no data in input file\n");
    strclose(istr);
    if (! set_member(intags, MassTag))
	error("inputdata: mass data missing\n");
    if (! set_member(intags, PosTag))
	error("inputdata: position data missing\n");
    if (! set_member(intags, VelTag))
	error("inputdata: velocity data missing\n");
#if defined(ENTROPY)
    if (! set_member(intags, EntFuncTag))
	error("inputdata: entropy function data missing\n");
#else
    if (! set_member(intags, UinternTag))
	error("inputdata: internal energy data missing\n");
#endif
    ngas = 0;
    for (p = btab; p < btab+nbody; p++) {
#if defined(ENTROPY)
	Type(p) = BODY | (EntFunc(p) >= 0 ? GAS : 0);

#else
	Type(p) = BODY | (Uintern(p) >= 0 ? GAS : 0);

#endif
	if (Gas(p))
	    ngas++;
    }
#if defined(EXTGRAV)
    if (! strnull(gspfile)) {
	gstr = stropen(gspfile, "r");
	get_history(gstr);
	gravgsp = get_gsprof(gstr);
	strclose(gstr);
    } else
	gravgsp = NULL;
#endif
    if (scanopt(options, "reset-time"))
	tnow = 0.0;
}


void startoutput(stream ostr, string defv[])
{
    int i;

    fprintf(ostr, "\n%s  [v%s]\n", headline, getversion());
    for (i = 1; defv[i][0] == ';'; i++)
        fprintf(ostr, "    %s\n", defv[i]+1);
    fprintf(ostr, "\n %9s %9s %9s %9s %9s %9s %9s\n", "nbody", "ngas",
	    "nsmooth", "courant", "dtime", "dtout", "tstop");
    fprintf(ostr, " %9d %9d %9d %9.4f %9.6f %9.6f %9.4f\n\n", nbody, ngas,
	    nsmooth, courant, dtime, dtout, tstop);
    fprintf(ostr, " %9s", "gamma");
#if defined(RADIATING)
    fprintf(ostr, " %9s %9s", "uradmax", "lambmax");
#endif
#if defined(DIFFUSING)
    fprintf(ostr, " %9s", "sigmastar");
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    fprintf(ostr, " %9s", "opacity");
#endif
#if defined(CONDUCTING)
    fprintf(ostr, " %9s", "conduct");
#endif
#if defined(STARFORM)
    fprintf(ostr, " %9s %9s %9s", "starprob", "rhoindx", "udotindx");
#endif
    fprintf(ostr, "\n");
    fprintf(ostr, " %9.4f", gamma0);
#if defined(RADIATING)
    fprintf(ostr, " %9.4g %9.4g", uradmax, lambmax);
#endif
#if defined(DIFFUSING)
    fprintf(ostr, " %9.4g", sigmastar);
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    fprintf(ostr, " %9.4g", opacity);
#endif
#if defined(CONDUCTING)
    fprintf(ostr, " %9.4g", conduct);
#endif
#if defined(STARFORM)
    fprintf(ostr, " %9.4g %9.2f %9.2f", starprob, rhoindx, udotindx);
#endif
    fprintf(ostr, "\n\n");
#if defined(GRAVITY)
    fprintf(ostr, " %9s %9s", "eps", "usequad");
#if !defined(QUICKSCAN)
    fprintf(ostr, " %9s", "theta");
#endif
#endif
    fprintf(ostr, " %9s %9s %9s %9s\n", "alpha", "beta", "slope", "fdrag");
#if defined(GRAVITY)
    fprintf(ostr, " %9.4f %9s", eps, usequad ? "true" : "false");
#if !defined(QUICKSCAN)
    fprintf(ostr, " %9.2f", theta);
#endif
#endif
    fprintf(ostr, " %9.2f %9.2f %9.4f %9.4f\n", alpha, beta, slope0, fdrag);
    if (! strnull(options))
        fprintf(ostr, "\n\toptions: %s\n", options);
    fflush(ostr);
    if (! strnull(savefile))
	savestate(savefile);
}


void outputhead(stream ostr)
{
    fprintf(ostr, "\n");
    fprintf(ostr, "########## %9s %9s    ################"
	    "############################\n", "Time", "Nstep");
    fprintf(ostr, "########## %9.4f %9d    ################"
	    "############################\n", tnow, nstep);
    fflush(ostr);
}


void outputdata(stream ostr)
{
    real teff;
    int n;
    string reqtags[MaxBodyFields], *outtags, *alltags;
    char namebuf[256];
    struct stat buf;
    stream outstr;

    diagnostics();
    fprintf(ostr, "\n %9s %9s %9s %9s %9s %9s %9s %7s\n",
	   "Etot", "Eint", "Ekin", "Epot", "Erad",
	    "|Jtot|", "|Vcom|", "CPUtot");
    fprintf(ostr, " %9.5f %9.5f %9.5f %9.5f %9.5f %9.6f %9.6f %7.1f\n",
	   etot[0], etot[1], etot[2], etot[3], eradiate,
	    absv(amvec), absv(cmvel), cputime());
    teff = tnow + (dtime > 0 ? dtime/8 : 0);
    if (! strnull(outfile) && teff >= tout) {
	n = 0;
	if (scanopt(outputs, MassTag) || (nstep == 0))
	    reqtags[n++] = MassTag;
	reqtags[n++] = PosTag;
	reqtags[n++] = VelTag;
#if defined(ENTROPY)
#  if defined(ADIABATIC)
	if (scanopt(outputs, UinternTag) || (nstep == 0))
	    reqtags[n++] = EntFuncTag;
#  else
	reqtags[n++] = EntFuncTag;
#  endif
#else
#  if defined(ISOTHERMAL)
	if (scanopt(outputs, UinternTag) || (nstep == 0))
	    reqtags[n++] = UinternTag;
#  else
	reqtags[n++] = UinternTag;
#  endif
#endif
	reqtags[n] = NULL;
	outtags = burststring(outputs, ", ");
	alltags = set_union(reqtags, outtags);
	sprintf(namebuf, outfile, nstep);
	if (stat(namebuf, &buf) != 0) {
	    outstr = stropen(namebuf, "w");
	    put_history(outstr);
	} else
	    outstr = stropen(namebuf, "a");
	put_snap(outstr, &btab, &nbody, &tnow, alltags);
	strclose(outstr);
	free(outtags);
	free(alltags);
	fprintf(ostr, "\n\tdata output to file %s at time %f\n",
		namebuf, tnow);
	tout += dtout;
    }
    fflush(ostr);
    if (! strnull(savefile))
	savestate(savefile);
}


local void diagnostics(void)
{
    register bodyptr p;
    real uint, velsq;
    vector tmpv;
    matrix tmpt;

    mtot = 0;
    etot[1] = etot[2] = etot[3] = 0;
    CLRM(keten);
    CLRM(peten);
    CLRV(amvec);
    CLRV(cmpos);
    CLRV(cmvel);
    for (p = btab; p < btab+nbody; p++) {
	mtot += Mass(p);
	if (Gas(p)) {
#if defined(ENTROPY)
	    uint = EntFunc(p) * rpow(Rho(p), gamma0 - 1) / (gamma0 - 1);
#else
	    uint = Uintern(p);
#endif
	    etot[1] += Mass(p) * uint;
	}
	DOTVP(velsq, Vel(p), Vel(p));
	etot[2] += 0.5 * Mass(p) * velsq;
#if defined(GRAVITY)
	etot[3] += 0.5 * Mass(p) * Phi(p);
#elif defined(EXTGRAV)
	etot[3] += Mass(p) * Phi(p);
#endif
	MULVS(tmpv, Vel(p), 0.5 * Mass(p));
	OUTVP(tmpt, tmpv, Vel(p));
	ADDM(keten, keten, tmpt);
	MULVS(tmpv, Pos(p), Mass(p));
	OUTVP(tmpt, tmpv, Acc(p));
	ADDM(peten, peten, tmpt);
	CROSSVP(tmpv, Vel(p), Pos(p));
	MULVS(tmpv, tmpv, Mass(p));
	ADDV(amvec, amvec, tmpv);
	MULVS(tmpv, Pos(p), Mass(p));
	ADDV(cmpos, cmpos, tmpv);
	MULVS(tmpv, Vel(p), Mass(p));
	ADDV(cmvel, cmvel, tmpv);
    }
    etot[0] = etot[1] + etot[2] + etot[3] + eradiate;

    DIVVS(cmpos, cmpos, mtot);
    DIVVS(cmvel, cmvel, mtot);
}


void savestate(string pattern)
{
    char namebuf[256];
    stream str;

    sprintf(namebuf, pattern, nstep & 1);
    str = stropen(namebuf, "w!");
    put_string(str, "program", getargv0());
    put_string(str, "version", getversion());
    put_string(str, "headline", headline);
    put_data(str, "gamma", RealType, &gamma0, 0);
#if defined(RADIATING)
    put_data(str, "uradmax", RealType, &uradmax, 0);
    put_data(str, "lambmax", RealType, &lambmax, 0);
#endif
#if defined(CONDUCTING)
    put_data(str, "conduct", RealType, &conduct, 0);
#endif
#if defined(DIFFUSING)
    put_data(str, "sigmastar", RealType, &sigmastar, 0);
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    put_data(str, "opacity", RealType, &opacity, 0);
#endif
#if defined(STARFORM)
    put_data(str, "starprob", RealType, &starprob, 0);
    put_data(str, "rhoindx", RealType, &rhoindx, 0);
    put_data(str, "udotindx", RealType, &udotindx, 0);
#endif
    put_data(str, "alpha", RealType, &alpha, 0);
    put_data(str, "beta", RealType, &beta, 0);
    put_data(str, "nsmooth", IntType, &nsmooth, 0);
    put_data(str, "nbucket", IntType, &nbucket, 0);
    put_data(str, "slope", RealType, &slope0, 0);
    put_data(str, "courant", RealType, &courant, 0);
    put_data(str, "dtime", RealType, &dtime, 0);
    put_data(str, "fdrag", RealType, &fdrag, 0);
#if defined(GRAVITY)
    put_data(str, "eps", RealType, &eps, 0);
    put_data(str, "usequad", BoolType, &usequad, 0);
#if !defined(QUICKSCAN)
    put_data(str, "theta", RealType, &theta, 0);
#endif
#elif defined(EXTGRAV)
    put_string(str, "gspfile", gspfile);
    if (! strnull(gspfile))
	put_gsprof(str, gravgsp);
#endif
    put_string(str, "options", options);
    put_string(str, "outputs", outputs);
    put_data(str, "tstop", RealType, &tstop, 0);
    put_data(str, "dtout", RealType, &dtout, 0);
    put_data(str, "nstep", IntType, &nstep, 0);
    put_data(str, "levmax", IntType, &levmax, 0);
    put_data(str, "rsize", RealType, &rsize, 0);
    put_data(str, "eradiate", RealType, &eradiate, 0);
    put_data(str, "tout", RealType, &tout, 0);
    put_data(str, "tnow", RealType, &tnow, 0);
#if defined(STARFORM)
    put_data(str, "randstate", ByteType, randstate, RANDSIZE, 0);
#endif
    put_data(str, "nbody", IntType, &nbody, 0);
    put_data(str, "ngas", IntType, &ngas, 0);
    put_data(str, "btab", AnyType, btab, nbody, sizeof(body), 0);
    strclose(str);
}


void restorestate(string file)
{
    stream str, gstr;
    string program, version;

    str = stropen(file, "r");
    program = get_string(str, "program");
    version = get_string(str, "version");
    if (! streq(program, getargv0()) ||
	  ! streq(version, getversion()))
	eprintf("[%s: input state file may be outdated]\n", getargv0());
    headline = get_string(str, "headline");
    get_data(str, "gamma", RealType, &gamma0, 0);
#if defined(RADIATING)
    put_data(str, "uradmax", RealType, &uradmax, 0);
    put_data(str, "lambmax", RealType, &lambmax, 0);
#endif
#if defined(CONDUCTING)
    get_data(str, "conduct", RealType, &conduct, 0);
#endif
#if defined(DIFFUSING)
    get_data(str, "sigmastar", RealType, &sigmastar, 0);
#endif
#if defined(DIFFUSING) || defined(OPAQUE)
    get_data(str, "opacity", RealType, &opacity, 0);
#endif
#if defined(STARFORM)
    get_data(str, "starprob", RealType, &starprob, 0);
    get_data(str, "rhoindx", RealType, &rhoindx, 0);
    get_data(str, "udotindx", RealType, &udotindx, 0);
#endif
    get_data(str, "alpha", RealType, &alpha, 0);
    get_data(str, "beta", RealType, &beta, 0);
    get_data(str, "nsmooth", IntType, &nsmooth, 0);
    get_data(str, "nbucket", IntType, &nbucket, 0);
    get_data(str, "slope", RealType, &slope0, 0);
    get_data(str, "courant", RealType, &courant, 0);
    get_data(str, "dtime", RealType, &dtime, 0);
    get_data(str, "fdrag", RealType, &fdrag, 0);
#if defined(GRAVITY)
    get_data(str, "eps", RealType, &eps, 0);
    get_data(str, "usequad", BoolType, &usequad, 0);
#if !defined(QUICKSCAN)
    get_data(str, "theta", RealType, &theta, 0);
#endif
#elif defined(EXTGRAV)
    gspfile = get_string(str, "gspfile");
    if (! strnull(gspfile))
	gravgsp = get_gsprof(str);
    else
	gravgsp = NULL;
#endif
    options = get_string(str, "options");
    outputs = get_string(str, "outputs");
    get_data(str, "tstop", RealType, &tstop, 0);
    get_data(str, "dtout", RealType, &dtout, 0);
    get_data(str, "nstep", IntType, &nstep, 0);
    get_data(str, "levmax", IntType, &levmax, 0);
    get_data(str, "rsize", RealType, &rsize, 0);
    get_data(str, "eradiate", RealType, &eradiate, 0);
    get_data(str, "tout", RealType, &tout, 0);
    get_data(str, "tnow", RealType, &tnow, 0);
#if defined(STARFORM)
    get_data(str, "randstate", ByteType, randstate, RANDSIZE, 0);
    (void) setstate(randstate);
#endif
    get_data(str, "nbody", IntType, &nbody, 0);
    get_data(str, "ngas", IntType, &ngas, 0);
    btab = (bodyptr) allocate(nbody * sizeof(body));
    get_data(str, "btab", AnyType, btab, nbody, sizeof(body), 0);
    strclose(str);
}
