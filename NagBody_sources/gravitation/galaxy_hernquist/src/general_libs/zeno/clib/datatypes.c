/*
 * DATATYPES.C: table and utility routines for data type package.
 */

#include "stdinc.h"
#include "getparam.h"
#include "datatypes.h"
#include <string.h>

/*
 * TYPEINFO: struct for info about data types.
 */

typedef struct {
    string type;				/* name from datatypes.h    */
    string name;				/* C-style decleration      */
    string fmt1;				/* low-precision format     */
    string fmt2;				/* high-precision format    */
    int size;					/* size of type in bytes    */
} typeinfo;

/*
 * TYPETABLE: table of basic types and associated data.
 */

local typeinfo typetable[] = {
    { AnyType,    "any",    "%#o",  "%#o",     sizeof(byte),   },
    { CharType,   "char",   "%c",   "%c",      sizeof(char),   },
    { ByteType,   "byte",   "%#o",  "%#o",     sizeof(byte),   },
    { ShortType,  "short",  "%#o",  "%#o",     sizeof(short),  },
    { IntType,    "int",    "%#o",  "%#o",     sizeof(int),    },
    { LongType,   "long",   "%#lo", "%#lo",    sizeof(long),   },
    { FloatType,  "float",  "%#g",  "%15.8e",  sizeof(float),  },
    { DoubleType, "double", "%#g",  "%23.16e", sizeof(double), },
    { NULL,       NULL,     NULL,   NULL,      0,              },
};

local typeinfo *find_type(string);		/* lookup by 1st element    */
local typeinfo *find_name(string);		/* lookup by C-style name   */

/*
 * TYPE_LENGTH: compute total size of basic or compound type.
 */

int type_length(string type)
{
    int count;
    char *tp;

    count = 0;
    for (tp = type; *tp != (char) NULL; tp++)	/* loop over all elements   */
	count += find_type(tp)->size;		/* add size of each element */
    return (count);
}

/*
 * TYPE_NAME: return C-style name of type.
 */

string type_name(string type)
{
    return (find_type(type)->name);
}

/*
 * NAME_TYPE: return basic type of C-style name.
 */

string name_type(string name)
{
    return (find_name(name)->type);
}

/*
 * TYPE_FMT: return format string for type.
 */

string type_fmt(string type, bool lowprec)
{
    return (lowprec ? find_type(type)->fmt1 : find_type(type)->fmt2);
}

/*
 * TYPE_BASE: return basic type of homogenious compound type.
 */

string type_base(string type)
{
    char *tp;

    for (tp = type; *tp != (char) NULL; tp++)
	if (*tp != type[0])
	    error("base_type in %s: type %s inhomogenious\n",
		  getargv0(), type);
    return (tp - 1);
}

/*
 * FIND_TYPE: return typeinfo for given type.  Note: for compound
 * types, the typeinfo for the first element is returned.
 */

local typeinfo *find_type(string type)
{
    typeinfo *tp;

    for (tp = typetable; tp->type != NULL; tp++)
	if (tp->type[0] == type[0])
	    return (tp);
    error("find_type in %s: type %c unknown\n", getargv0(), type[0]);
    return (NULL);
}

/*
 * FIND_NAME: return typeinfo for named type.
 */

local typeinfo *find_name(string name)
{
    typeinfo *tp;

    for (tp = typetable; tp->type != NULL; tp++)
	if (streq(tp->name, name))
	    return (tp);
    error("find_name in %s: type %s unknown\n", getargv0(), name);
    return (NULL);
}

#ifdef TESTBED

string defv[] = {		";Test datatypes routines",
    "type=" IntType,		";Data type string",
    "name=int",			";C-style type name",
    "VERSION=1",		";Josh Barnes  15 October 1997",
    NULL,
};

void main(int argc, string argv[])
{
    string type, name;

    initparam(argv, defv);
    type = getparam("type");
    name = getparam("name");
    printf("type_length(\"%s\") = %d\n", type, type_length(type));
    printf("type_name(\"%s\") = %s\n", type, type_name(type));
    printf("type_base(\"%s\") = %s\n", type, type_base(type));
    printf("type_fmt(\"%s\", TRUE) = %s\n", type, type_fmt(type, TRUE));
    printf("type_fmt(\"%s\", FALSE) = %s\n", type, type_fmt(type, FALSE));
    printf("name_type(\"%s\") = %s\n", name, name_type(name));
}

#endif
