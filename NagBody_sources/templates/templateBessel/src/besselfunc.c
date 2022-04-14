/*==============================================================================
 MODULE: besselfunc.c				[templateBessel]
 Written by: M.A. Rodriguez-Meza
 Starting date: January 2005
 Purpose: Compute bessel functions
 Language: C
 Use: ...
 Routines and functions: ...
 Modules, routines and external headers:
 Comments and notes: Reference: Abramowitz & Stegun, chapter 9
 Info: M.A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.astro.inin.mx/mar
 
 Major revisions: July 2007; January 2018;
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.																		!
 ==============================================================================*/

#include "../../../General_libs/general/stdinc.h"
#include "../../../General_libs/general/getparam.h"
#include "../../../General_libs/math/mathfns.h"

#include "protodefs.h"

// BESSI0: modified Bessel function I0(x).

double bessi0(double x)
{
    double t, tt, ti, u;

    t = ABS(x) / 3.75;
    tt = t * t;
    if (tt < 1.0) {
	  u = 1.0 + tt * (3.5156229 + tt * (3.0899424 +
			tt * (1.2067492 + tt * (0.2659732 + tt * (0.0360768 +
					  tt * 0.0045813)))));
	  return (u);
    } else {
	  ti = 1.0 / t;
	  u = 0.39894228 + ti * (0.01328592 + ti * (0.00225319 +
			ti * (-0.00157565 + ti * (0.00916281 + ti * (-0.02057706 +
				ti * (0.02635537 + ti * (-0.01647633 +
						      ti * 0.00392377)))))));
	  return (u * exp(ABS(x)) / sqrt(ABS(x)));
    }
}


// BESSI1: modified Bessel function I1(x). 

double bessi1(double x)
{
    double t, tt, ti, u;

    t = x / 3.75;
    tt = t * t;
    if (tt < 1.0) {
	  u = 0.5 + tt * (0.87890594 + tt * (0.51498869 +
			tt * (0.15084934 + tt * (0.02658733 + tt * (0.00301532 +
					  tt * 0.00032411)))));
	  return (u * x);
    } else {
	  if (t < 0.0)
	    error("bessi1: invalid for x < -3.75\n");
	  ti = 1.0 / t;
	  u = 0.39894228 + ti * (-0.03988024 + ti * (-0.00362018 + 
			ti * (0.00163801 + ti * (-0.01031555 + ti * (0.02282967 +
					ti * (-0.02895312 + ti * (0.01787654 +
						      ti * -0.00420059)))))));
	  return (u * exp(ABS(x)) / sqrt(ABS(x)));
    }
}

// BESSK0: modified Bessel function K0(x). 

double bessk0(double x)
{
    double t, tt, ti, u;

    if (x < 0.0)
	error("bessk0: negative argument\n");
    t = x / 2.0;
    if (t < 1.0) {
	  tt = t * t;
	  u = -0.57721566 + tt * (0.42278420 + tt * (0.23069756 +
			tt * (0.03488590 + tt * (0.00262698 + tt * (0.00010750 +
					  tt * 0.00000740)))));
	  return (u - log(t) * bessi0(x));
    } else {
	  ti = 1.0 / t;
	  u = 1.25331414 + ti * (-0.07832358 + ti * (0.02189568 +
			ti * (-0.01062446 + ti * (0.00587872 + ti * (-0.00251540 +
					  ti * 0.00053208)))));
	  return (u * exp(- x) / sqrt(x));
    }
}


// BESSK1: modified Bessel function K1(x). 

double bessk1(double x)
{
    double t, tt, ti, u;

    if (x < 0.0)
	  error("bessk1: negative argument\n");
    t = x / 2.0;
    if (t < 1.0) {
	  tt = t * t;
	  u = 1.0 +
	    tt * (0.15443144 + tt * (-0.67278579 + tt * (-0.18156897 +
			      tt * (-0.01919402 + tt * (-0.00110404 +
					  tt * -0.00004686)))));
	  return (u / x + log(t) * bessi1(x));
    } else {
	  ti = 1.0 / t;
	  u = 1.25331414 + ti * (0.23498619 + ti * (-0.03655620 +
			ti * (0.01504268 + ti * (-0.00780353 + ti * (0.00325614 +
			ti * -0.00068245)))));
	  return (u * exp(- x) / sqrt(x));
    }
}


