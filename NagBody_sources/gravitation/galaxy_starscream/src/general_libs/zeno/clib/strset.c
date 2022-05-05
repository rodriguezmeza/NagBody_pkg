/*
 * STRSET.C: code for simple string set package.
 */

#include "stdinc.h"
#include "strset.h"

#include <stdarg.h>

/*
 * Sets of strings are represented by null-terminated vectors of string
 * pointers.  This strategy keeps the code simple, but also implies
 * that many operations (construction, union, intersection, difference)
 * are O(N^2) in the number of elements.  The maximum number of elements
 * is therefore limited to a relatively small value.
 */

#define MaxElements  64

/*
 * SET_CONS: construct a set from a null-terminated list of elements.
 */

string *set_cons(string first, ...)
{
    int n;
    string name, strs[MaxElements+1];
    va_list ap;

    n = 0;
    strs[n++] = first;
    strs[n] = NULL;
    va_start(ap, first);
    while ((name = va_arg(ap, string)) != NULL)
        if (set_member(strs, name))
	    eprintf("[set_cons: element %s duplicated]\n", name);
        else {
	    if (n == MaxElements)
	        error("set_cons: too many elements\n");
	    strs[n++] = name;
	    strs[n] = NULL;
	}
    va_end(ap);
    return (set_copy(strs));
}

/*
 * SET_COPY: copy set (but not member strings).
 */

string *set_copy(string *set)
{
    return ((string *) copxstr(set, sizeof(string)));
}

/*
 * SET_LENGTH: return count of set members.
 */

int set_length(string *set)
{
    return (xstrlen(set, sizeof(string)) - 1);
}

/*
 * SET_MEMBER: test if element is member of set.
 */

bool set_member(string *set, string element)
{
    string *sp;

    for (sp = set; *sp != NULL; sp++)
        if (streq(*sp, element))
	    return (TRUE);
    return (FALSE);
}

/*
 * SET_SUBSET: return true if every member of set2 is in set1.
 */

bool set_subset(string *set1, string *set2)
{
    string *sp;

    for (sp = set2; *sp != NULL; sp++)
        if (! set_member(set1, *sp))
	    return (FALSE);
    return (TRUE);
}

/*
 * SET_EQUAL: return true if sets have same members.
 */

bool set_equal(string *set1, string *set2)
{
    return (set_subset(set1, set2) && set_subset(set2, set1));
}

/*
 * SET_UNION: construct union of two sets.
 */

string *set_union(string *set1, string *set2)
{
    int n;
    string *sp, strs[MaxElements+1];

    n = 0;
    for (sp = set1; *sp != NULL; sp++)
        strs[n++] = *sp;
    for (sp = set2; *sp != NULL; sp++)
        if (! set_member(set1, *sp)) {
	    if (n > MaxElements)
	        error("set_union: too many elements\n");
	    strs[n++] = *sp;
	}
    strs[n] = NULL;
    return (set_copy(strs));
}

/*
 * SET_INTER: construct intersection of two sets.
 */

string *set_inter(string *set1, string *set2)
{
    int n;
    string *sp, strs[MaxElements+1];

    n = 0;
    for (sp = set2; *sp != NULL; sp++)
        if (set_member(set1, *sp))
	    strs[n++] = *sp;
    strs[n] = NULL;
    return (set_copy(strs));
}

/*
 * SET_DIFF: construct difference of two sets.
 */

string *set_diff(string *set1, string *set2)
{
    int n;
    string *sp, strs[MaxElements+1];

    n = 0;
    for (sp = set1; *sp != NULL; sp++)
        if (! set_member(set2, *sp))
	    strs[n++] = *sp;
    strs[n] = NULL;
    return (set_copy(strs));
}

#ifdef TESTBED

#include "getparam.h"

string defv[] = {
    "set1=foo,bar,fum,fie",
    "set2=fie,bar,waldo,frodo",
    "VERSION=1.0",
    NULL,
};

void print_set(string, string *);

void main(int argc, string argv[])
{
    string *set1, *set2;

    initparam(argv, defv);
    set1 = burststring(getparam("set1"), ", ");
    set2 = burststring(getparam("set2"), ", ");
    printf("set_length(set1) = %d\n", set_length(set1));
    printf("set_subset(set1,set2) = %s\n",
	   set_subset(set1, set2) ? "TRUE" : "FALSE");
    printf("set_equal(set1,set2) = %s\n",
	   set_equal(set1, set2) ? "TRUE" : "FALSE");
    print_set("set_copy(set1):", set_copy(set1));
    print_set("set_union(set1,set2):", set_union(set1, set2));
    print_set("set_inter(set1,set2):", set_inter(set1, set2));
    print_set("set_diff(set1,set2):", set_diff(set1, set2));
}

void print_set(string label, string *set)
{
    string *sp;

    printf("%s", label);
    for (sp = set; *sp != NULL; sp++)
        printf(" %s", *sp);
    printf("\n");
}

#endif /* TESTBED */
