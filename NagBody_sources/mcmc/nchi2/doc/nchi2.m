't" t
.TH nagbdmcmc 1 "January 2018" UNIX "NagBody PROJECT"
.na
.nh   

.SH NAME
nagbdmcmc - TEMPLATE for developing ODE solvers codes
.SH SYNOPSIS
\fBnagbdmcmc\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ] 
.sp

.SH DESCRIPTION
\fBnagbdmcmc\fR - MCMC to analyse data versus models.

.SH OPTIONS
All the options have the structure
.sp
\fIoption_name\fR=<option_value>
.sp
except option \fB-help\fR.
.sp

Options and their possible values are:

.IP "\fB-help\fR" 12
By writting

.sp
nagbdmcmc -help
.sp

you will get the list of all parameters and their default values. An option may have an alias which is a short name of the option. If an option has an alias in the list above it comes after its description surrounded by brackets tagged with 'a:'. For example,

.sp
option_name=<option_value>	... Description ... [a: opt]
.sp
here 'opt' is the short name of the option. In the description of the options below, when an option has an alias it will be noted in the same way but before its description.

.IP "\fBparamfile\fR" 12
is the name file with the values of the input parameters. Overwrite parameters
values below. You may input this filename by only writing:
.sp
nagbdmcmc paramfile=parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. Comment lines start with an "%". Follow each name option with a blank space and the option value. The order of the option lines does not matter.  Also you may create an example input file by executing
.sp
nagbdmcmc
.sp
This will run the \fBnagbdmcmc\fR code with default values and when it finish you will have in your running directory the file "parameters_null-nagbdmcmc-usedvalues". Now you may edit this file to adapt to your own plotting parameters. This file can be overwritten so it may be helpful to change this file name to whatever apropriate.

.IP "\fBx\fR" 12
is the value of x to start the computation.

.IP "\fBdx\fR" 12
[a: dx] is the x step. Can be given as a ration of two numbers. Can not be zero.

.IP "\fBxstop\fR" 12
is the x value to stop the computation.

.IP "\fBmaxnsteps\fR" 12
[a: maxn]is the maximum number of steps.

.IP "\fBout\fR" 12
[a: o] You give here the name of the file where it will be written the computation results. The output is written in column form.

.IP "\fBdxout\fR" 12
[a: dxo] You give here the x step to generate output of computation.

.IP "\fBdxoutinfo\fR" 12
[a: dxoinfo] You give here the x step to generate output of computation to the standard output.

.IP "\fBoptions\fR" 12
[a: opt] you may give here various code behavior options.

.SH EXAMPLES
nagbdmcmc x=1.5 dx=1/8

.sp
You will have computation written in a file named with option \fBout\fR. You may plot the solutions using \fBnplot2d\fR code.

.SH SEE ALSO
nplot2d(1)

.SH COPYRIGHT
Copyright (C) 1999-2018
.br
Mario A. Rodriguez-Meza
.br
