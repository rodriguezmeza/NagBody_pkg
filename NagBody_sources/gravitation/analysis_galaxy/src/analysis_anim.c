/*==============================================================================
	MODULE: analysis_anim.c			[analysis_galaxy]
	Written by: M.A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Initialize datanaly_md
	Language: C
	Use: 'snap_anim_gnuplot();', 'rdf_anim_gnuplot();', 'vel_anim_gnuplot();',
		'snap_conversion();', 'thermo_avg();', 'snap_less_nbody();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: mar@nuclear.inin.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: January 2007;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

//#include "../../../General_libs/general/stdinc.h"
//#include "../../../General_libs/math/vectmath.h"
#include "globaldefs.h"
#include "protodefs.h"

//#include <sys/stat.h>

void snap_conversion(void)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			if (exist_snap) {
				++snapcount;
				outputdata(cmd.out, cmd.outfmt, cmd.infmt, snapcount, gd.nbody,
					gd.tnow, &hdr, cmd.options);
				free_memory();
				free(bodytab);
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
			++step;
	}
}

void snap_less_nbody(void)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;
	real reductionFac;							// PONERLO EN LA LINEA DE COMANDOS...

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			if (exist_snap) {
				++snapcount;
				SnapLessNBody(&gd.nbody, snapcount, reductionFac, cmd.options, 
							&gdtree, &gdforce);
				outputdata(cmd.out, cmd.outfmt, cmd.infmt, snapcount, gd.nbody,
					gd.tnow, &hdr, cmd.options);
				free_memory();
				free(bodytab);
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
			++step;
	}
}

void extract_sets(void)
{
	int i, step=0, ndim, nbodytmp, snapcount=0;
	bool exist_snap;
	bodyptr p, q, btabtmp;
	real tnow=0.0, tmass;
	vector cmpos, cmvel, tmpv;

	inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
	if (!exist_snap)
		error("\n\nextract_sets: file %s doesnot exist\n",cmd.in);

	for (i=0; i<gd.nbodiesSets; i++)
		if (gd.bodyIDMin[i]>gd.nbody || gd.bodyIDMax[i]>gd.nbody)
			error("\n\nPrintBodiesSets: IDMin or IDMax out of range\n");

	nbodytmp = gd.bodyIDMax[0]-gd.bodyIDMin[0]+1
			  +gd.bodyIDMax[1]-gd.bodyIDMin[1]+1
			  +gd.bodyIDMax[2]-gd.bodyIDMin[2]+1;

    btabtmp = (bodyptr) allocate(nbodytmp * sizeof(body));

	q = btabtmp;
	DO_BODY(p, bodytab, bodytab+gd.nbody) {
		for (i=0; i<gd.nbodiesSets; i++) {
			if (Id(p)>=gd.bodyIDMin[i] && Id(p)<=gd.bodyIDMax[i]) {
				Mass(q) = Mass(p);
				Id(q) = Id(p);
				Type(q) = Type(p);
				SETV(Pos(q), Pos(p));
				SETV(Vel(q), Vel(p));
				++q;
			}
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
	fprintf(stdout,"\n\nNumber of bodies in sets : %d %d\n",nbodytmp,q-btabtmp);
	fprintf(stdout,"CM Pos: %g %g %g\n",cmpos[0],cmpos[1],cmpos[2]);
	fprintf(stdout,"CM Vel: %g %g %g\n",cmvel[0],cmvel[1],cmvel[2]);

	free(bodytab);
	gd.nbody=nbodytmp;
    bodytab = (bodyptr) allocate(gd.nbody * sizeof(body));
	q=bodytab;
	DO_BODY(p, btabtmp, btabtmp+nbodytmp) {
		Mass(q) = Mass(p);
		Id(q) = Id(p);
		Type(q) = Type(p);
		SUBV(Pos(q), Pos(p), cmpos);
		SUBV(Vel(q), Vel(p), cmvel);
		++q;
	}

	outputdata(cmd.out, cmd.outfmt, cmd.infmt, snapcount, gd.nbody, tnow, &hdr,
		cmd.options);
}

void acf(string fname, string fnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
    char namebuf[200], namebuf1[200];
	int ndim;

	gd.Acf_flag=0;
	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			if (exist_snap) {
				++snapcount;
				EvalAcf(fname, fnametmp);
//				if (gd.Acf_flag) {
//					animation_acf(fname, snapcount);
//					free(xval); free(yval);
//				}
				free_memory();
				free(bodytab);
				gd.Acf_flag=0;
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
		++step;
	}
}

void locate_bodiesid(string fname, string fnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			if (exist_snap) {
				++snapcount;
				LocateBodiesID(fname, fnametmp);
				free_memory();
				free(bodytab);
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
		++step;
	}
}

void two_bdh_galaxies(string fname, string fnametmp)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim, nbody1, nbody2, nbody;
	bodyptr p, q, btab1, btab2;
	string options;
	real tnow;

	inputdata(gd.filenames[0], gd.filenamesfmt[0], gd.headerfmt, step, &nbody1, &ndim, &tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
	if (!exist_snap)
		error("\n\ntwo_bdh_galaxies: file %s doesnot exist\n",gd.filenames[0]);
    btab1 = (bodyptr) allocate(nbody1 * sizeof(body));
	q=btab1;
	DO_BODY(p, bodytab, bodytab+nbody1) {
		Mass(q) = Mass(p);
		Id(q) = Id(p);
		Type(q) = Type(p);
		SETV(Pos(q), Pos(p));
		SETV(Vel(q), Vel(p));
		++q;
	}
	free(bodytab);

	inputdata(gd.filenames[1], gd.filenamesfmt[1], gd.headerfmt, step, &nbody2, &ndim, &tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
	if (!exist_snap)
		error("\n\ntwo_bdh_galaxies: file %s doesnot exist\n",gd.filenames[1]);
    btab2 = (bodyptr) allocate(nbody2 * sizeof(body));
	q=btab2;
	DO_BODY(p, bodytab, bodytab+nbody2) {
		Mass(q) = Mass(p);
		Id(q) = Id(p);
		Type(q) = Type(p);
		SETV(Pos(q), Pos(p));
		SETV(Vel(q), Vel(p));
		++q;
	}
	free(bodytab);

	nbody=nbody1+nbody2;
	tnow=0.0;
    bodytab = (bodyptr) allocate(nbody * sizeof(body));

	q=bodytab;
	DO_BODY(p, btab1, btab1+nbody1) {
		Mass(q) = Mass(p);
		Id(q) = Id(p);
		Type(q) = Type(p);
		ADDV(Pos(q), Pos(p), gd.cmpos1);
		ADDV(Vel(q), Vel(p), gd.cmvel1);
		++q;
	}
	DO_BODY(p, btab2, btab2+nbody2) {
		Mass(q) = Mass(p);
		Id(q) = Id(p);
		Type(q) = Type(p);
		ADDV(Pos(q), Pos(p), gd.cmpos2);
		ADDV(Vel(q), Vel(p), gd.cmvel2);
		++q;
	}

	outputdata(cmd.out,cmd.outfmt,cmd.infmt,snapcount,nbody,tnow,&hdr,options);
	free(btab1);
	free(btab2);
	free(bodytab);
}

