/*
 * BODY.H: Structure and accessor macro definitions for a vanilla-flavored
 * body structure to be used in conjunction with SnapShot binary files.
 */
#ifndef _body_h
#define _body_h


#define _body_h_dens

typedef struct {
    real   bodymass;			/* mass of body			    */
    vector bodyphase[2];		/* phase-space coordinates	    */
    real   bodyphi;			/* gravitational potential	    */
    vector bodyacc;			/* gravitational acceleration	    */
    real   bodyaux;			/* misc. real value assoc. w. body  */
    int    bodykey;			/* misc. int. value assoc. w. body  */
#ifdef _body_h_dens
    real   bodydens;			/* density associated w. body       */
    real   bodyeps;                     /* softening length w. body         */
#endif
} NBbody;

typedef int  (*btiproc)(NBbody *, real, int);
typedef real (*btrproc)(NBbody *, real, int);


#define Body     NBbody

#define NBMass(b)  ((b)->bodymass)
#define NBPhase(b) ((b)->bodyphase)
#define NBPos(b)   ((b)->bodyphase[0])
#define NBVel(b)   ((b)->bodyphase[1])
#define NBPhi(b)   ((b)->bodyphi)
#define NBAcc(b)   ((b)->bodyacc)
#define NBAux(b)   ((b)->bodyaux)
#define NBKey(b)   ((b)->bodykey)
#ifdef _body_h_dens
#define NBDens(b)  ((b)->bodydens)
#define NBEps(b)   ((b)->bodyeps)
#endif

#endif
