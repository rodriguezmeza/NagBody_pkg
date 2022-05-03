/*
 * SNAPCENTER.C: routines to compute center-of-mass coordinates.
 */

#include "../clib/stdinc.h"
#include "../clib/vectmath.h"
#include "phatbody.h"
#include "snapcenter.h"

/*
 * WEIGHT: macro to select weight field for c. of m. calculations.
 */

#define Weight(bp, woff)  SelectReal(bp, woff)

/*
 * SNAPCMPOS, SNAPCMVEL: compute center of mass position and velocity.
 */

void snapcmpos(vector cmpos, bodyptr btab, int nbody, int woff)
{
    double mtot, cmtmp[NDIM];
    bodyptr bp;

    mtot = 0.0;
    CLRV(cmtmp);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	mtot = mtot + Weight(bp, woff);
	ADDMULVS(cmtmp, Pos(bp), Weight(bp, woff));
    }
    DIVVS(cmpos, cmtmp, mtot);
}

void snapcmvel(vector cmvel, bodyptr btab, int nbody, int woff)
{
    double mtot, cmtmp[NDIM];
    bodyptr bp;

    mtot = 0.0;
    CLRV(cmtmp);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	mtot = mtot + Weight(bp, woff);
	ADDMULVS(cmtmp, Vel(bp), Weight(bp, woff));
    }
    DIVVS(cmvel, cmtmp, mtot);
}

/*
 * SNAPCENTER: transform to the weighted center of mass.
 */

void snapcenter(bodyptr btab, int nbody, int woff)
{
    vector cmpos, cmvel;
    bodyptr bp;

    snapcmpos(cmpos, btab, nbody, woff);
    snapcmvel(cmvel, btab, nbody, woff);
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
	SUBV(Pos(bp), Pos(bp), cmpos);
	SUBV(Vel(bp), Vel(bp), cmvel);
    }
}
