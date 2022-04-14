/*==============================================================================
	HEADER: models.h			[model]
	Written by: M.A. Rodriguez-Meza
	Starting date: October 1999
	Purpose: proto definitios of some model routines
	Language: C
	Use: '#include "...."
	Use in routines and functions: datanaly_md (main)
	External headers: None
	Comments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Major revisions: July 23, 2007
	Copyright: (c) 2005-2010 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their	use.
==============================================================================*/

#ifndef _models_h
#define _models_h

/*
int Nb, Nd, Nh, Nbi, Nbf, Ndi, Ndf, Nhi, Nhf;
int    Nb1,Nd1,Nh1,Nbi1,Nbf1,Ndi1,Ndf1,Nhi1,Nhf1;
int Nb2,Nd2,Nh2,Nbi2,Nbf2,Ndi2,Ndf2,Nhi2,Nhf2;
int Ng, Ngi, Ngf;
int Nsph, Nsphi, Nsphf;
int Ndm, Ndmi, Ndmf;
*/
//  gcc11 :: To avoid Error :: duplicate symbol '_Ndm' in:
local int Nb, Nd, Nh, Nbi, Nbf, Ndi, Ndf, Nhi, Nhf;
local int    Nb1,Nd1,Nh1,Nbi1,Nbf1,Ndi1,Ndf1,Nhi1,Nhf1;
local int Nb2,Nd2,Nh2,Nbi2,Nbf2,Ndi2,Ndf2,Nhi2,Nhf2;
local int Ng, Ngi, Ngf;
local int Nsph, Nsphi, Nsphf;
local int Ndm, Ndmi, Ndmf;


void BDH_Galaxy(void);
void Two_BDH_Galaxies_Coll(void);
void outputdata_BDH_Galaxy(void);
void outputdata_Two_BDH_Galaxies_Coll(void);
void dist_rho_bdh(void);
real pr_b(real, real, real, real);
real pr_d(real, real, real, real);
real pr_h(real, real, real, real);


//int Ndist1, Ndisti1, Ndistf1, Ndist2, Ndisti2, Ndistf2;
//  gcc11 :: To avoid Error :: duplicate symbol '_Ndist2' in:
local int Ndist1, Ndisti1, Ndistf1, Ndist2, Ndisti2, Ndistf2;

//----------------- PASAR ESTA RUTINA A DATANALY_GRAV --------------------------
void mass_equatorial_slab(void);

#endif // !_cmodels_h
