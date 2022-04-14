//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)


#ifndef _sphcode_h
#define _sphcode_h

#include "sphdefs.h"


#define RANDSIZE 32
#define ETA2     0.01

global string infile;
global string outfile;
global string savefile;
global string restfile;
global string gspfile;
global real gamma0;
global real uradmax;
global real lambmax;
global real sigmastar;
global real opacity;
global real conduct;
global real starprob;
global real rhoindx;
global real udotindx;
global real alpha;
global real beta;
global int nsmooth;
global int nbucket;
global real slope0;
global real dtime;
global real courant;
global real fdrag;
global gsprof *gravgsp;
global string outputs;
global real tstop;
global real dtout;
global string headline;
global int nstep;
global int levmax;
global real eradiate;
global real tnow;
global real tout;
global char randstate[RANDSIZE];
global int nbody;
global int ngas;
global bodyptr btab;


void inputdata(void);
void startoutput(stream, string []);
void outputhead(stream);
void outputdata(stream);
void savestate(string);
void restorestate(string);

#endif // ! _sphcode_h
