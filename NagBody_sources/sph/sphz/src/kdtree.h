//
// Edited and modified by M.A. Rodriguez-Meza (2010)
//
// Adapted from zeno 
// (see http://www.ifa.hawaii.edu/faculty/barnes/barnes.html)

#ifndef _kdtree_h
#define _kdtree_h

typedef struct {
    real minb[3];
    real maxb[3];
} bound;


typedef struct {
    bound bnd;
    int dim;
    real split;
    int first;
    int last;
} kdnode;


#define KDROOT		1
#define Lower(i)	((i)<<1)
#define Upper(i)	(((i)<<1)+1)
#define Parent(i)	((i)>>1)
#define Sibling(i) 	(((i)&1)?(i)-1:(i)+1)

#define SetNext(i) {				\
    while (i&1)					\
        i=i>>1;					\
    ++i;					\
}


typedef struct {
    int ngas;
    bodyptr *bptr;
    bound bnd;
    int nnode;
    int nsplit;
    kdnode *ntab;
} kdcontext, *kdxptr;

kdxptr init_kdtree(bodyptr, int, int);
void build_kdtree(kdxptr, int);
void finish_kdtree(kdxptr);

#endif  // ! _kdtree_h
