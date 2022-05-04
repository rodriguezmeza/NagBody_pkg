/*==============================================================================
	MODULE: lic.c					[General_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date: May 2006
	Purpose: routines to license the code
	Language: C
	Use: 'MakeLicItem();'
	Routines and functions:
	Modules, routines and external headers: stdinc.h, machines.h, <string.h>,
		<sys/types.h>, <sys/stat.h>
	Coments and notes:
	Info: M.A. Rodriguez-Meza
		Depto. de Fisica, ININ
		Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
		e-mail: marioalberto.rodriguez@inin.gob.mx
		http://www.astro.inin.mx/mar

	Mayor revisions: January 2007;
	Copyright: (c) 2005-2011 Mar.  All Rights Reserved
================================================================================
	Legal matters:
	The author does not warrant that the program and routines it contains
	listed below are free from error or suitable for particular applications,
	and he disclaims all liability from any consequences arising from their use.
==============================================================================*/

#include "machines.h"
#include "lic.h"

#include "stdinc.h"
#include "../io/inout.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

local void lic_root(void);
local void lic(void);
local void ReadInLicItem(stream, license *);
local void GetSystemDate(int *, int *, int *);
local int GetSystemMonth(char *);
local bool AcceptLic(license *, char *, date *, date *);
local void PrintOutLicItems(license *);

#define NAGUSERFILE		".nagusertmp"
#define NUSERS			100

#define ANABELLE		"/home/anabelle"
#define MARVEGA			"/home/mar"
#define MARDRACO		"/home/mar"
#define MARKANBALAM		"/global/home/ext/mardgz"
#define RUSLANKANBALAM	"/global/home/ext/ruslan"
#define MARLAP			"/Users/mar"
#define GUEST			"/home/guest"
#define GLOBALUSER		"/home/mar"
#define JAMR			"/Users/jamr"
#define NACHO			"/home/Nacho"
#define GIACOMO			"/home/giacomo"
#define JERONIMO		"/home/ikari"
#define CESAR			"/home/cesar"
#define JERONIMO2		"/home/jeronimo"
#define YAYA			"/home/yaya"

#define LMGRCODE		"lmgr"

//char usersv[7][200] = {
string usersv[] = {
ANABELLE, MARVEGA, MARDRACO, MARKANBALAM, RUSLANKANBALAM, MARLAP, JAMR,
NACHO, JERONIMO, GIACOMO, CESAR, JERONIMO2, YAYA,
NULL
};

void LicDriver(char *rcode)
{
//	fprintf(stdout,"\nRunning code : %s\n",rcode);
	if (strcmp(rcode,LMGRCODE)==0)
		lic_root();
	else
		lic();
}

string rootv[] = {MARLAP, NULL };

local void lic_root(void)
{
    char filebuf[200], sysbuf[200], pathname[200];
    struct stat buf;
    stream instr;
//	char username[200];
	char namebuf[200];
	int i;
	license userlic;

//	fprintf(stdout,"\nInside lic_root\n");

	sprintf(sysbuf,"echo $HOME > %s",NAGUSERFILE);
	system(sysbuf);
    sprintf(filebuf, "%s",NAGUSERFILE);
    if (stat(filebuf, &buf) != 0)
		error("\nUnable to open file %s\n",filebuf);
	else
        instr = stropen(filebuf, "r");
	ReadInString(instr, pathname);
    fclose(instr);

    sprintf(filebuf, "%s/.nagbody",pathname);
    if (stat(filebuf, &buf) != 0)
		error("\nerror opening file %s\n",filebuf);
	else
        instr = stropen(filebuf, "r");
	ReadInLicItem(instr, &userlic);
    fclose(instr);

	for (i=0; i< NUSERS; i++)
		if (rootv[i]==NULL)
			error("\nlic : [1] : invalid root ... quitting\n");
		else {
			sprintf(namebuf,"%s",rootv[i]);
			if(strcmp(userlic.pathusername,namebuf)==0)
				return;
		}

	error("\nlic : [2] : invalid root ... quitting\n");
}

local void lic(void)
{
    char gfilebuf[200], filebuf[200], sysbuf[200], pathname[200];
    struct stat buf;
    stream instr;
	int i, flag=1;
	license *userlic;
	date *ld, *sd;

//	fprintf(stdout,"\nInside lic\n");

	userlic = (license *)allocate(sizeof(license));
	ld = (date *)allocate(sizeof(date));
	sd = (date *)allocate(sizeof(date));

	sprintf(gfilebuf,"%s/.nagbody",GLOBALUSER);
	if (stat(gfilebuf, &buf) !=0) {
		sprintf(sysbuf,"echo $HOME > .nagusertmp");
		system(sysbuf);
		sprintf(filebuf, "%s",NAGUSERFILE);
		if (stat(filebuf, &buf) != 0)
			exit(0);
		else
			instr = stropen(filebuf, "r");

		ReadInString(instr, pathname);
		fclose(instr);

		sprintf(filebuf, "%s/.nagbody",pathname);
		if (stat(filebuf, &buf) != 0)
			error("\nlic [1]: invalid user or machine : %s ... quitting\n\n",pathname);
		else
			instr = stropen(filebuf, "r");
	} else
		instr = stropen(gfilebuf, "r");

	ReadInLicItem(instr, userlic);
    fclose(instr);

//	PrintOutLicItems(userlic);

//	if (!(sscanf(userlic->date, "%ld-%ld-%ld", &ld->year, &ld->month, &ld->day) == 3))
	if (!(sscanf(userlic->date, "%d-%d-%d", &ld->year, &ld->month, &ld->day) == 3))
		error("\nlic : date must be in the form of yyyy-mm-dd\n\n");
	GetSystemDate(&sd->year, &sd->month, &sd->day);

	for (i=0; i< NUSERS; i++)
		if (usersv[i]==NULL) {
			flag=0; break;
		} else
			if (AcceptLic(userlic, usersv[i], ld, sd))
				return;

	if (flag==0 && userlic->id==0)
		if (AcceptLic(userlic, GUEST, ld, sd))
			return;
		else
			error("\nlic [2]: invalid user or machine ... quitting\n\n");
	else
		error("\nlic [3]: invalid user or machine ... quitting\n\n");
}

#undef NAGUSERFILE
#undef NUSERS

#undef ANABELLE
#undef MARVEGA
#undef MARDRACO
#undef MARKANBALAM
#undef MARLAP
#undef GUEST
#undef GLOBALUSER
#undef JAMR
#undef JERONIMO
#undef CESAR
#undef JERONIMO2
#undef YAYA

#undef LMGRCODE

local bool AcceptLic(license *lic, char *user, date *ld, date *sd)
{
	bool flag=FALSE;
												// Puede ser necesario en
												// algunas maquinas agregar
												// un espacio en blanco al
												// final de cada linea en el 
												// archivo '.nagbody'.
												// Este fue el caso de las
												// Mac's del Chato...
	if(strcmp(lic->pathusername,user)==0 )
		if (sd->year <= ld->year) {					// HAY QUE CAMBIAR LAS 
			if (sd->month <= ld->month)				// PREGUNTAS A NEGACION...
				if (sd->day <= ld->day)
					flag=TRUE;
		}

	return flag;
}

local void ReadInLicItem(stream instr, license *userlic)
{
	char word[200];

	ReadInString(instr, word);
	userlic->id = atoi(word);

	ReadInString(instr, word);
	strcpy(userlic->date, word);

	ReadInString(instr, word);
	strcpy(userlic->licid, word);

	ReadInString(instr, word);
	strcpy(userlic->machineid, word);

	ReadInString(instr, word);
	strcpy(userlic->passwd, word);

	ReadInString(instr, word);
	strcpy(userlic->pathusername, word);
}

local void PrintOutLicItems(license *userlic)
{
	fprintf(stdout,"\n\nLic.ID : %d\n",userlic->id);
	fprintf(stdout,"Lic.DATE : %s\n",userlic->date);
	fprintf(stdout,"Lic.LICID : %s\n",userlic->licid);
	fprintf(stdout,"Lic.MACHINEID : %s\n",userlic->machineid);
	fprintf(stdout,"Lic.PASSWD : %s\n",userlic->passwd);
	fprintf(stdout,"Lic.PATHUSERNAME : %s\n\n",userlic->pathusername);
}

#define NAGDATEFILE		".datetmp"

#define JAN				"Jan"
#define FEB				"Feb"
#define MAR				"Mar"
#define APR				"Apr"
#define MAY				"May"
#define JUN				"Jun"
#define JUL				"Jul"
#define AUG				"Aug"
#define SEP				"Sep"
#define OCT				"Oct"
#define NOV				"Nov"
#define DEC				"Dec"

string monthv[] = {
JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC, NULL
};

local void GetSystemDate(int *year, int *month, int *day)
{
    char filebuf[256], sysbuf[200];
    struct stat buf;
    stream instr;
	char date[200], m[4], hour[9], dw[4], hc[4];

	sprintf(sysbuf,"date > %s",NAGDATEFILE);
	system(sysbuf);
    sprintf(filebuf, "%s",NAGDATEFILE);
    if (stat(filebuf, &buf) != 0) {
		fprintf(stdout,"\nUnable to open file %s\n",filebuf);
        exit(0);
    } else {
        instr = stropen(filebuf, "r");
    }
	ReadInLineString(instr, date);
    fclose(instr);

//	if (!(sscanf(date, "%3s %3s %ld %8s %3s %ld", dw, m, day, hour, hc, year) == 6))
	if (!(sscanf(date, "%3s %3s %d %8s %3s %d", dw, m, day, hour, hc, year) == 6))
		error("\nGetSystemDate : date format is incorrect\n\n");

	*month = GetSystemMonth(m);
}

local int GetSystemMonth(char *m)
{
	int i;

	for (i=0; i< 12; i++)
		if(strcmp(m,monthv[i])==0)
			break;
	return ++i;
}

#undef NAGDATEFILE

#undef JAN
#undef FEB
#undef MAR
#undef APR
#undef MAY
#undef JUN
#undef JUL
#undef AUG
#undef SEP
#undef OCT
#undef NOV
#undef DEC
