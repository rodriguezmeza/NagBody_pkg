/*==============================================================================
	MODULE: analysis_anim.c			[analysis_grav]
	Written by: Mario A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Initialize analysis_grav
	Language: C
	Use: 'snap_anim_gnuplot();', 'rdf_anim_gnuplot();', 'vel_anim_gnuplot();',
		'snap_conversion();', 'thermo_avg();', 'snap_less_nbody();'
	Routines and functions:
	Modules, routines and external headers:
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: January 2007;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"

/*
#define tmpcgmfile	"tmp.cgm"	// Nombre que se usa en los scripts de gnuplot

void snap_anim_gnuplot(string snapname, string snapnametmp)
{
	int step, moreSteps, stepLimit;
	bool exist_snap;
    char   buf[200], namebuf[200], namebuf1[200];
	int ndim;

	stepLimit = 10000;
	step = 0;
	moreSteps = 1;
	while (moreSteps) {
//		readin_snap(step, &exist_snap);
		inputdata(cmd.in, cmd.infmt, step, &gd.nbody, &ndim, &gd.tnow, 
			&exist_snap, &hdr, cmd.options, gd.model_comment);
		if (exist_snap) {
//			PrintSnap(stdout, gd.nbody, cmd.options);
			PrintSnap(snapname, snapnametmp, gd.nbody, cmd.options);
			if (scanopt(cmd.options, "snap_cgm")) {
				sprintf(buf,"gnuplot snap_cgm.gnu");
				printf("\nsystem: %s",buf);
				system(buf);
				sprintf(namebuf1, cmd.out, step);
				sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
				sprintf(buf,"mv %s %s",tmpcgmfile,namebuf);
				printf("\nsystem: %s\n",buf);
				system(buf);
			}
		}
		++step;
		if (step >= stepLimit) moreSteps = 0;
	}
}

void rdf_anim_gnuplot(string fname, string fnametmp)
{
	int step=0, moreSteps=1, stepLimit=10000;
	bool exist_snap;
    char   buf[200], namebuf[200], namebuf1[200];
	int ndim;

	gd.Rdf_flag=0;
	while (moreSteps) {
//		readin_snap(step, &exist_snap);
		inputdata(cmd.in, cmd.infmt, step, &gd.nbody, &ndim, &gd.tnow, 
			&exist_snap, &hdr, cmd.options, gd.model_comment);
		if (exist_snap) {
			EvalRdf(fname, fnametmp);

			if (scanopt(cmd.options, "rdf_cgm") && gd.Rdf_flag) {
				sprintf(buf,"gnuplot rdf_cgm.gnu");
				printf("\nsystem: %s",buf);
				system(buf);
				sprintf(namebuf1, cmd.out, step);
				sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
				sprintf(buf,"mv %s %s",tmpcgmfile,namebuf);
				printf("\nsystem: %s\n",buf);
				system(buf);
				gd.Rdf_flag=0;
			}
		}
		
		++step;
		if (step >= stepLimit) moreSteps = 0;
	}
}

void vel_anim_gnuplot(string fname, string fnametmp)
{
	int step=0, moreSteps=1, stepLimit=10000;
	bool exist_snap;
    char   buf[200], namebuf[200], namebuf1[200];
	int ndim;

	gd.Vel_flag=0;
	while (moreSteps) {
//		readin_snap(step, &exist_snap);
		inputdata(cmd.in, cmd.infmt, step, &gd.nbody, &ndim, &gd.tnow, 
			&exist_snap, &hdr, cmd.options, gd.model_comment);
		if (exist_snap) {
			EvalVelDist(fname, fnametmp);

			if (scanopt(cmd.options, "vel_cgm") && gd.Vel_flag) {
				sprintf(buf,"gnuplot vel_cgm.gnu");
				printf("\nsystem: %s",buf);
				system(buf);
				sprintf(namebuf1, cmd.out, step);
				sprintf(namebuf, "%s/%s.%s", cmd.basedir, namebuf1, cmd.outfmt);
				sprintf(buf,"mv %s %s",tmpcgmfile,namebuf);
				printf("\nsystem: %s\n",buf);
				system(buf);
				gd.Vel_flag=0;
			}
		}
		++step;
		if (step >= stepLimit) moreSteps = 0;
	}
}

#undef tmpcgmfile
*/

void snap_conversion(void)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
//			Header_to_Global();
			if (exist_snap) {
				++snapcount;
//				Global_to_Header();
				outputdata(cmd.out, cmd.outfmt, cmd.infmt, snapcount, 
					gd.nbody, gd.tnow, &hdr, cmd.options);
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
//	real reductionFac;							// PONERLO EN LA LINEA DE COMANDOS...

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
//			Header_to_Global();
			if (exist_snap) {
				++snapcount;
//				Global_to_Header();
				SnapLessNBody(&gd.nbody, snapcount, cmd.reductionFac, cmd.options, 
							&gdtree, &gdforce);
				outputdata(cmd.out, cmd.outfmt, cmd.infmt, snapcount, 
					gd.nbody, gd.tnow, &hdr, cmd.options);
				free_memory();
				free(bodytab);
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
			++step;
	}
}

void groups_catalog(void)
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;

	cmd.fsnap=0;						// Solo un snap es procesado...

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
//			inputdata(step, &exist_snap);
//			inputdata(in, infmt, step, &nbody, &ndim, &exist_snap);
//			inputdata(in, infmt, step, &nbody, &ndim, &tnow, &exist_snap, 
//				&hdr, options, model_comment);
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, &gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
//printf("\n\nAqui voy [1]\n");
			if (exist_snap) {
				++snapcount;
				GroupsCatalog(snapcount);
				free_memory();
				free(bodytab);
			}
		} else
			if (step >= cmd.fsnap) moreSteps = 0;
			++step;
	}
}

/*
void thermo_avg(void)
{
    stream outstr;
    struct stat buf;
	int i,npoints;
	real sumx, sumy, sumerry;

	if (scanopt(cmd.options, "errorbar"))
		if (!readin_pl_file_3c(cmd.in,gd.column1,gd.column2,gd.column3,
			gd.row1,gd.row2,&npoints))
			error("\nthermo_avg : error reading file %s\n\n",cmd.in);
	else
		if (!readin_pl_file(cmd.in,gd.column1,gd.column2,gd.row1,gd.row2,&npoints))
			error("\nthermo_avg : error reading file %s\n\n",cmd.in);

	sumx = sumy = 0.;
	if (scanopt(cmd.options, "errorbar")) sumerry = 0.;
	for (i=0; i<npoints; i++) {
		sumx += xval[i];
		sumy += yval[i];
		if (scanopt(cmd.options, "errorbar")) sumerry += zval[i];
	}
	sumx /= npoints; sumy /= npoints;
	if (scanopt(cmd.options, "errorbar")) sumerry /= npoints;

    if (scanopt(cmd.options, "overwrite"))			// Append or overwrite ...
        outstr = stropen(cmd.out, "w!");
	else
		if (stat(cmd.out, &buf) != 0)
			error("\nthermo_avg : Error : file %s doesnot exist to append!\n\n",
				cmd.out);
		else
			outstr = stropen(cmd.out, "a");
	if (scanopt(cmd.options, "errorbar"))
		fprintf(outstr,"%g %g %g %g\n", sumx, sumy, sumy-sumerry, sumy+sumerry);
	else
		fprintf(outstr,"%g %g\n",sumx,sumy);
	fclose(outstr);
}
*/

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
