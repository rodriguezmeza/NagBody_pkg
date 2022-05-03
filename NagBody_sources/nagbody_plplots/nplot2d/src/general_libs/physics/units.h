

#define  CM_PER_MPC			3.085678e24
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7



/* Module to convert from arbitrary system units in which
	G=1, M* = 1 and R* = 1; in galaxy simulations.
	G = 6.672x10^(-8) cm^3 /(g s^2).
	Unit of mass (M*):
		gram= 2.2x10^11 Mo (Solar mass) = 2.2x10^11 x 1.989x10^33 g
			= 4.3758 x 10^44 g.
	Unit of length (R*):
		centimeter	= 40 kpc = 40 x 3.086 x 10^21 cm
					= 1.2344 x 10^23 cm.
	Unit of time (s^2 = (R*)^3/(G M*)):
		(R*)^3 = 1.880908803584 x 10^69 cm^3
		G M*	= 2.91953376 x 10^37 cm^3 / s^2
		(R*)^3/(G M*)	= 0.6442497186893 x 10^32 s^2
		s	= 0.802651679553 x 10^16 sec
		second	= 250 Myr = 250 x 10^6 x 365 x 24 x 60 x 60 s
				= 7.884 x 10^15 s.
*/
/*
#define GRAM		4.97592E44
#define CENTIMETER	2.2817E22
#define SEC			8.02651679553E15
*/

/* Module to convert from arbitrary system units in which
	G=1, M* = 1 and R* = 1; in colapse simulations.
	G = 6.672x10^(-8) cm^3 /(g s^2).
	Unit of mass (M*):
		gram= 1 Mo (Solar mass) = 1.989x10^33 g
	Unit of length (R*):
		centimeter	= 5 x 10^16 cm = 0.016 pc.
	Unit of time (s^2 = (R*)^3/(G M*)): 
		s	= 9.705314 x 10^11 sec
*/

#define GRAM		1.989E33
#define CENTIMETER	5.0E16
#define SEC			9.705314E11


