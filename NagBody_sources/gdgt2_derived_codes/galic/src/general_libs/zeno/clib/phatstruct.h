/*
 * PHATSTRUCT.H: definitions for "phat" structures.
 */

#ifndef _phatstruct_h
#define _phatstruct_h

/* 
 * PS_FIELD: structure describing one element of a phat structure.  A
 * phat structure is described by an array of ps_fields.  The zeroth
 * field refers to the entire structure and its length is the length
 * of the whole structure.  The array is terminated by a ps_field with
 * a NULL name.  Compound data types are supported, but they must be
 * homogenious.
 */

typedef struct {
    string type;		/* type code following "datatypes.h"        */
    string name;		/* text for name of field                   */
    int offset;			/* offset in bytes from start of struct     */
    int length;			/* length in bytes of this field            */
} ps_field;

/*
 * BADOFFSET: offset value for structure fields which aren't actualized.
 */

#define BadOffset  -1

/*
 * LAYOUT_STRUCT: compute offsets and length for structure with named fields.
 */

void layout_struct(ps_field *, string *);

/*
 * NEW_FIELD: define a new field of given type and name.
 */

void new_field(ps_field *, string, string);

/*
 * DEFINE_STRUCT, DEFINE_OFFSET: describe total length and individual field
 * offsets of static structure.
 */

void define_struct(ps_field *, string, int);
void define_offset(ps_field *, string, int);

/*
 * Generic selector macros.
 */

#define SelectByte(p,o)    (*((byte *) ((byte *)(p) + (o))))
#define SelectChar(p,o)    (*((char *) ((byte *)(p) + (o))))
#define SelectShort(p,o)   (*((short *) ((byte *)(p) + (o))))
#define SelectBool(p,o)    (*((bool *) ((byte *)(p) + (o))))
#define SelectInt(p,o)     (*((int *) ((byte *)(p) + (o))))
#define SelectLong(p,o)    (*((long *) ((byte *)(p) + (o))))
#define SelectFloat(p,o)   (*((float *) ((byte *)(p) + (o))))
#define SelectDouble(p,o)  (*((double *) ((byte *)(p) + (o))))
#define SelectReal(p,o)    (*((real *) ((byte *)(p) + (o))))
#define SelectVect(p,o)    ((real *) ((byte *)(p) + (o)))

#endif  /* ! _phatstruct_h */
