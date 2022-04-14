/*==============================================================================
	MODULE: analysis_anim.c			[analysis_md]
	Written by: Mario A. Rodriguez-Meza
	Starting date: January 2005
	Purpose: Initialize analysis_md
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
	Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#include "globaldefs.h"


typedef struct {
	real length, time, mass, energy, temperature, density, pressure;
} units_data;

local void PrintUnits(stream, units_data);

void snap_conversion(void)								// CHECK 2D --- OK!!!
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
				&gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			Header_to_Global();
			if (exist_snap) {
				++snapcount;
				Global_to_Header();
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

void snap_less_nbody(void)								// CHECK 2D --- OK!!!
{
	int step=0, moreSteps=1, snapcount=0;
	bool exist_snap;
	int ndim;
//	real reductionFac;							// PONERLO EN LA LINEA DE COMANDOS...

	while (moreSteps) {
		if (step >= cmd.isnap && step <= cmd.fsnap) {
			inputdata(cmd.in, cmd.infmt, gd.headerfmt, step, 
				&gd.nbody, &ndim, &gd.tnow, 
				&exist_snap, &hdr, cmd.options, gd.model_comment);
			Header_to_Global();
			if (exist_snap) {
				++snapcount;
				Global_to_Header();
				SnapLessNBody(&gd.nbody, snapcount, cmd.reductionFac, cmd.options, 
							&gdtree, &gdforce);
				Global_to_Header();
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

void thermo_avg(void)									// CHECK 2D --- OK!!!
{
    stream outstr;
    struct stat buf;
	int i,npoints;
	real sumx, sumy, sumerry;

	if (scanopt(cmd.options, "errorbar")) {
		if (!readin_pl_file_3c(cmd.in,gd.column1,gd.column2,gd.column3,
			gd.row1,gd.row2,&npoints))
			error("\nthermo_avg : error reading file %s\n\n",cmd.in);
	} else {
		if (!readin_pl_file(cmd.in,gd.column1,gd.column2,gd.row1,gd.row2,&npoints))
			error("\nthermo_avg : error reading file %s\n\n",cmd.in);
	}

	sumx = sumy = 0.;
	if (scanopt(cmd.options, "errorbar")) sumerry = 0.;
	for (i=0; i<npoints; i++) {
		sumx += npltd.xval[i];
		sumy += npltd.yval[i];
		if (scanopt(cmd.options, "errorbar")) sumerry += npltd.zval[i];
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
	if (scanopt(cmd.options, "errorbar")) {
		fprintf(outstr,"%g %g %g %g\n", sumx, sumy, sumy-sumerry, sumy+sumerry);
	} else {
		fprintf(outstr,"%g %g\n",sumx,sumy);
	}
	fclose(outstr);
}

#define unitsfile	"units.dat"

void units_conversion(void)								// CHECK 2D --- OK!!!
{
    stream outstr;
//	real length, time, mass, energy, temperature, density, pressure;
	units_data unit;

	outstr = stropen(unitsfile, "w!");

	unit.length = gd.unitLength;
	unit.energy = gd.unitEnergy;
	unit.mass = gd.unitMass;

	unit.time = rsqrt( gd.unitMass * rsqr(gd.unitLength) / gd.unitEnergy );
	unit.temperature = gd.unitEnergy / BOLTZMANNCONSTANT;
#if (NDIM==3)
	unit.density = gd.unitMass/Cube(gd.unitLength);
	unit.pressure = gd.unitEnergy/Cube(gd.unitLength);
#else
	unit.density = gd.unitMass/rsqr(gd.unitLength);
	unit.pressure = gd.unitEnergy/rsqr(gd.unitLength);
#endif

	PrintUnits(outstr, unit);
	PrintUnits(stdout, unit);

	fclose(outstr);
}

local void PrintUnits(stream outstr, units_data unit)	// CHECK 2D --- OK!!!
{
	fprintf(outstr,"\n\n");
	fprintf(outstr,"Quantities in MKS system of units: \n");
	fprintf(outstr,"length = %g meters\n", unit.length);
	fprintf(outstr,"energy = %g Joules\n", unit.energy);
	fprintf(outstr,"mass = %g Kilograms\n", unit.mass);
	fprintf(outstr,"time = %g sec\n", unit.time);
	fprintf(outstr,"temperature = %g Kelvin\n", unit.temperature);
#if (NDIM==3)
	fprintf(outstr,"density = %g Kilograms/m^3\n", unit.density);
	fprintf(outstr,"pressure = %g Joules/m^3\n", unit.pressure);
#else
	fprintf(outstr,"density = %g Kilograms/m^2\n", unit.density);
	fprintf(outstr,"pressure = %g Joules/m^2\n", unit.pressure);
#endif

	fprintf(outstr,"\nQuantities in other useful units: \n");
	fprintf(outstr,"length = %g Angstrom\n", unit.length/ANGSTROM);
	fprintf(outstr,"energy = %g Kelvin\n", unit.energy/BOLTZMANNCONSTANT);
	fprintf(outstr,"mass = %g grams\n", unit.mass/KILO);
	fprintf(outstr,"time = %g ps\n", unit.time/NANOSEC);
	fprintf(outstr,"temperature = %g Kelvin\n", unit.temperature);
	fprintf(outstr,"density = %g grams/cm^3\n", unit.density*KILO/1.0e6);
	fprintf(outstr,"pressure = %g Pa\n", unit.pressure);
}

#undef unitsfile
