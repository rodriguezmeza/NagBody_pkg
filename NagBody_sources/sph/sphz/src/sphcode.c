//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)


#define global

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
#include "cmdline_defs.h"
#include "protodefs.h"


int main(int argc, string argv[])
{
    initparam(argv, defv);
    headline = defv[0] + 1;
    startrun();
    startoutput(stdout, defv);
    if (strnull(restfile))
	initstep();
    outputdata(stdout);
    while (tstop - tnow > 0.01 * dtime) {
	macrostep();
	outputdata(stdout);
    }
    return (0);
}

