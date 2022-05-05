/*
 * BODYIO.C: Input/Output code for phat body arrays.
 */

#include "../clib/stdinc.h"
#include "../clib/filestruct.h"
#include "../clib/vectdefs.h"
#include "phatbody.h"
#include <string.h>

local void put_parameters(stream, int *, real *);
local void put_particles(stream, bodyptr *, int *, string *);

local void get_parameters(stream, int *, real *);
local void add_known_fields(stream);
local void get_particles(stream, bodyptr *, int *, string *);

local void set_mask(int *, int, int, int);

/*
 * PUT_SNAP: write snapshot to structured output stream.
 */

void put_snap(stream ostr, bodyptr *btab, int *nbody, real *tnow,
	      string *tags)
{
    put_set(ostr, SnapShotTag);
    put_parameters(ostr, nbody, tnow);
    put_particles(ostr, btab, nbody, tags);
    put_tes(ostr, SnapShotTag);
}

/*
 * PUT_PARAMETERS: write snapshot parameters.
 */

local void put_parameters(stream ostr, int *nbody, real *tnow)
{
    put_set(ostr, ParametersTag);
    put_data(ostr, NBodyTag, IntType, nbody, 0);
    put_data(ostr, TimeTag, RealType, tnow, 0);
    put_tes(ostr, ParametersTag);
}

/*
 * PUT_PARTICLES: write particle data to output file.
 */

local void put_particles(stream ostr, bodyptr *btab, int *nbody,
			 string *tags)
{
    string *tp, type;
    ps_field *pbf;
    int mask[32];

    put_set(ostr, ParticlesTag);
    for (tp = tags; *tp != NULL; tp++) {	/* loop over list of tags   */
	for (pbf = phatbody; pbf->name != NULL; pbf++)
	    if (streq(pbf->name, *tp))		/* look for name in struct  */
		break;
	if (pbf->name == NULL)
	    error("put_particles: field %s unknown\n", *tp);
	if (pbf->offset == BadOffset)
	    error("put_particles: field %s undefined\n", *tp);
	set_mask(mask, SizeofBody, pbf->offset, type_length(pbf->type));
	type = type_base(pbf->type);		/* get base type of datum   */
	if (strlen(pbf->type) == 1)
	    put_data_masked(ostr, *tp, type, *btab, *nbody, 0, mask);
	else
	    put_data_masked(ostr, *tp, type, *btab, *nbody, NDIM, 0, mask);
    }
    put_tes(ostr, ParticlesTag);
}

/*
 * GET_SNAP: read snapshot from structured input stream.
 */

bool get_snap(stream istr, bodyptr *btab, int *nbody, real *tnow,
	      string *tags, bool expand)
{
    if (get_tag_ok(istr, SnapShotTag)) {
	get_set(istr, SnapShotTag);
	get_parameters(istr, nbody, tnow);
	if (expand && *btab == NULL)		/* don't expand after alloc */
	    add_known_fields(istr);
	get_particles(istr, btab, nbody, tags);
	get_tes(istr, SnapShotTag);
	return (TRUE);
    } else
	return (FALSE);    
}

#ifndef TimeFuzz
#  define TimeFuzz  0.001               /* uncertainty in time comparison   */
#endif

bool get_snap_t(stream istr, bodyptr *btab, int *nbody, real *tnow,
		string *tags, bool expand, string times)
{
    bool success;
    int ntmp;
    real ttmp;

    success = FALSE;
    while (! success && get_tag_ok(istr, SnapShotTag)) {
	get_set(istr, SnapShotTag);
	get_parameters(istr, &ntmp, &ttmp);
	if (streq(times, "all") || within(ttmp, times, TimeFuzz)) {
	    success = TRUE;
	    *nbody = ntmp;
	    *tnow = ttmp;
	    if (expand && *btab == NULL)	/* don't expand after alloc */
		add_known_fields(istr);
	    get_particles(istr, btab, nbody, tags);
	}
	get_tes(istr, SnapShotTag);
    }
    return (success);
}

/*
 * GET_PARAMETERS: read snapshot parameters.
 */

local void get_parameters(stream istr, int *nbody, real *tnow)
{
    get_set(istr, ParametersTag);
    if (get_tag_ok(istr, NBodyTag))
        get_data(istr, NBodyTag, IntType, nbody, 0);
    else if (get_tag_ok(istr, NobjTag))
        get_data(istr, NobjTag, IntType, nbody, 0);
    else
        error("get_snap: need %s or %s in %s\n",
	      NBodyTag, NobjTag, ParametersTag);
    if (get_tag_ok(istr, TimeTag))
	get_data(istr, TimeTag, RealType, tnow, 0);
    else {
	*tnow = 0.0;
	eprintf("[get_snap: time defaults to zero]\n");
    }
    get_tes(istr, ParametersTag);
}

/*
 * ADD_KNOWN_FIELDS: add all known body fields present in input set.
 */

local void add_known_fields(stream istr)
{
    string *tp, tags[MaxBodyFields];
    ps_field *pbf;

    get_set(istr, ParticlesTag);
    tp = tags;
    for (pbf = phatbody; pbf->name != NULL; pbf++)
    	if (get_tag_ok(istr, pbf->name))
	    *tp++ = pbf->name;
    *tp = NULL;
    get_tes(istr, ParticlesTag);
    layout_body(tags, Precision, NDIM);
}

/*
 * GET_PARTICLES: read particle data from input file.
 */

local void get_particles(stream istr, bodyptr *btab, int *nbody,
			 string *tags)
{
    string *tp, type;
    int len, mask[32];
    ps_field *pbf;

    if (*btab == NULL) {
	*btab = (bodyptr) allocate(*nbody * SizeofBody);
	eprintf("[get_snap: allocating %d bodies (size %d bytes)]\n",
		*nbody, SizeofBody);
    }
    get_set(istr, ParticlesTag);
    tp = tags;
    if (get_tag_ok(istr, PhaseTag) &&
	  PosField.offset != BadOffset && VelField.offset != BadOffset) {
	len = type_length(PosField.type);
	if (VelField.offset - PosField.offset != len)
	    error("get_particles: %s and %s not contiguous", PosTag, VelTag);
	set_mask(mask, SizeofBody, PosField.offset, 2 * len);
	type = type_base(PosField.type);
	get_data_masked(istr, PhaseTag, type, *btab, *nbody, 2, NDIM, 0, mask);
	*tp++ = PosTag;
	*tp++ = VelTag;
    }
    for (pbf = phatbody; pbf->name != NULL; pbf++)
	if (get_tag_ok(istr, pbf->name) && pbf->offset != BadOffset) {
	    set_mask(mask, SizeofBody, pbf->offset, type_length(pbf->type));
	    type = type_base(pbf->type);
	    if (strlen(pbf->type) == 1)
		get_data_masked(istr, pbf->name, type, *btab, *nbody, 0, mask);
	    else
		get_data_masked(istr, pbf->name, type, *btab,
				*nbody, NDIM, 0, mask);
	    *tp++ = pbf->name;
	}
    *tp = NULL;
    get_tes(istr, ParticlesTag);
}

/*
 * SET_MASK: initialize byte mask for field.
 */

local void set_mask(int *mask, int size, int offset, int length)
{
    int i;

    i = 0;
    if (offset > 0)
	mask[i++] = -offset;			/* skip preceding stuff     */
    mask[i++] = length;				/* copy field itself        */
    mask[i++] = -(size - offset - length);	/* skip rest of structure   */
    mask[i] = 0;
}
