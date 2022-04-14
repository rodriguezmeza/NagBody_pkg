't" t
.TH ngenic 1 "January 2019" UNIX "NagBody PROJECT"
.na
.nh   

.SH NAME
ngenic - A code for cosmological initial conditions.
.SH SYNOPSIS
\fBngenic\fR
.sp

.SH DESCRIPTION
\fBngenic\fR - A code for cosmological initial conditions by Volker Springel.

This code can be used to create initial conditions for N-body
simulations of cosmological structure formation, within the standard
LCDM paradigm. Optionally the code can also create additional SPH
particles for mixed dark matter/gas simulations.  The code is an MPI
code and can be used to setup fairly large initial conditions if
desired. The output format is that of the GADGET code, in its basic
file format 1.

.SH OPTIONS
See documentation in doc/cosmo_background_1lpt_2lpt/ngenic directory.
.sp

.SH EXAMPLES
See documentation in doc/cosmo_background_1lpt_2lpt/ngenic directory.

.SH SEE ALSO
gdgt207lcdm(1), cola(1), lpicola(1) and mglpicola(1)

.SH COPYRIGHT
Copyright (C) 1999-2019
.br
Mario A. Rodriguez-Meza
.br
