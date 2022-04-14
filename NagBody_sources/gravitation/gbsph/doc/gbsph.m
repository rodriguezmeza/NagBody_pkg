't" t
.TH gbsph 1 "January 2005" UNIX "NagBody PROJECT"
.na
.nh   

.SH NAME
\fBgbsph\fR - Hierarchical (tree force calculation) simulation of a selfgravitating N-general body system

.SH SYNOPSIS
\fBgbsph\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ] 
.sp

.SH DESCRIPTION
\fBgbsph\fR - Simulates the evolution of an N-body system interacting with a potential
proportional to 1/r and/or a Yukawa potential.
The force is computed by walking a hierarchical tree. Code complexity is NlogN.

.SH OPTIONS
The options have the structure
.sp
\fIoption_name\fR = <option_value>

.sp
except option \fB-help\fR.
.sp
Options and their possible values are:

.IP "\fB-help\fR" 12
By writting

.sp
gbsph -help
.sp

you will get the list of all parameters and their default values.
An option may have an alias which is a short name of the option. If an option
has an alias in the list above it comes after its description
surrounded by brackets tagged with 'a:'. For example,

.sp
option_name=<option_value>	... Description ... [a: opt]

.sp
here 'opt' is the short name of the option. In the description of the options
below, when an option has an alias it will be noted in the same way but before
its description.

.IP "\fBparamfile\fR" 12
is the name file with the values of the input parameters. Overwrite parameters
values below. You may also input this filename by only writing:
.sp 
gbsph parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. Comment lines start
with an "%". Follow each name option with a blank space and the option value.
The order of the option lines does not matter.
Also you may create an example input file
by executing
.sp
gbsph
.sp
This will run the \fBgbsph\fR code with default values generating
a test data set internally,
and when it finish you will have in your
running directory the file "parameters_null-gbsph-usedvalues".
Now you may edit this file to adapt it
to your own simulation parameters.
It may be helpful to change this file name to whatever apropriate.

.IP "\fBforcecalc_method\fR" 12
[a: fmeth] force calculation method. Could be <barnes> (default), an interaction list is
generated; <normal>, the force is computed walking the tree in the stadard way;
<direct>, the force is computed directly (N^2 complexity);

.IP "\fBforce_models\fR" 12
[a: fmod] force models to include. Could be <gravity_plummer> (default);
<scalar_field_potential_eps>. You may give several models separated by a "coma".

.IP "\fBeps\fR" 12
it is the smoothing force length. The spline smoothing parameters is 'h=2.8*eps'.
Therefore, fixing the \fBeps\fR value it fixes the 'h' value.

.IP "\fBtheta\fR" 12
Cell opening angle parameter that defines the accuracy of force computation. 
Direct force calculation is done if \fBtheta\fR = <0>.
Ussually a value of <0.8> is apropriate to get good accuracy of force computation.
Default value is <1.0>.

.IP "\fBusequad\fR" 12
if true compute quadrupole moment contributions to the gravitational and Yukawa forces.

.IP "\fBnbody\fR" 12
it is the number of bodies to simulate.
.sp
Note: it is recommended that this number is given as power of 2.

.IP "\fBdtime\fR" 12
[a: dt] It is the integration time step. Could be given as a ratio, <dt1/dt2>.
Default value is <1/32>.

.IP "\fBtstop\fR" 12
is the time to stop the simulation. Default value is <2>.

.IP "\fBseed\fR" 12
it is the Random number seed for the test run. It can be any integer number, its 
default value is <-1>.

.IP "\fBicfile\fR" 12
[a: ic] you give here the name of the file with the N-body initial data. The data structure
of this file follows the standard snap format.
If this parameter is <null> then a Plummer sphere will be constructed internally
by the code.

.IP "\fBicfilefmt\fR" 12
[a: icfmt] is the format of the initial condition file name given with parameter \fBicfile\fR.
Formats are:  <snap-bin> (Binary), <snap-ascii> (ASCII), or <snap-pv> (ASCII-columns).
Initial condition info for bodies needed are: mass, position, velocity, ID. Also
needed are number of simulation bodies, dimension of the system, time to start
the simultions. All quantities must be given in code units.

.IP "\fBsnapout\fR" 12
[a: o] you give here the pattern of the snapshot names. The format follows
as the ones used in C-language for integers ("%0#d").

.IP "\fBsnapoutfmt\fR" 12
[a: ofmt] you tell the code the format of the snap output. 
See \fBicfilefmt\fR for the snap formats.

.IP "\fBdtout\fR" 12
this defines the time steps for output info. 
I can be given as an output frequency, just type in <dt1/dt2>.
Default value is <8/32>. If \fBdtime\fR = <1/32> then every
8 time step info will be written to the stdout.

.IP "\fBstatefile\fR" 12
[a: state] you give here the name of a file where the run state will be saved. If it is null no run
state will be saved.

.IP "\fBrestorefile\fR" 12
[a: restore] if it is not null a run will be restarted from the data stored in this file. 
Note also that if you restart from a file state, the code path saved in the state file
must be the same as the path you are restarting. You may find the path used in the simulation
by seeing the file "parameters-usedvalues".

.IP "\fBoptions\fR" 12
[a: opt] you may give here various code behavior options.
They are, <reset-time> (inputdata); <out-phi> (outputdata); 
<out-acc> (outputdata). Also, you may give the opening cell 
criterion, <bh86> or <sw94>.
If you save a state file, this parameter is saved.

.IP "\fBdm_lambda\fR" 12
parameter lambda of Yukawa potential.

.IP "\fBdm_alpha\fR" 12
parameter alpha of Yukawa potential.

.IP "\fBG\fR" 12
Local value of the constant of gravity. Its value sets
the system of units used by the code. Natural code units system
must be defined by G=1, the default value of this parameter.

.IP "\fBdm_a\fR" 12
Scalar field slope coeficient that characterizes the rate of temporal change
of the SF conected adiabatically.

.IP "\fBdm_time\fR" 12
Scalar field perturbation time. Adiabatic perturbation
with a scalar field can be done setting \fBdm_time\fR
to non-zero value. Scalar field strength parameter
is a function of time: 

.sp
alpha(t) = (Exp(a t)-1) alpha0 / (Exp(a tp) - 1)

.sp
where a is the slope parameter (\fBdm_a\fR parameter) and tp is the time for adiabatic
perturbation, i.e., \fBdm_time\fR here, and alpha0 is the
final value of alpha, that is set by parameter \fBdm_alpha\fR. 

.IP "\fBeps_pot\fR" 12
energy scale of the external potential if included.

.IP "\fBsigma_pot\fR" 12
length scale of the external potential if included.

.IP "\fBx_pot\fR" 12
x-position of an external potential.

.IP "\fBy_pot\fR" 12
y-position of an external potential.

.IP "\fBz_pot\fR" 12
z-position of an external potential.


.SH UNITS
Units are such that gravitational constant is G=1. Units of 
mass, length, and velocity that are given in an initial condition
file must be consistent with this value of G. If no initial 
condition file name is given in the parameter \fBicfile\fR then
a Plummer model will be constructed its units are such that,
G=M=-4E=1, where M and E are total mass and total energy of the sphere
(see Aarseth, SJ, Henon, M, and Wielen, R (1974) Astron. and Astrophys. 37, 183).


.SH STOPPING THE CODE
Once the \fBgbsph\fR is running you may allways stop it by executing the command

echo > stop

You must be in the same directory were the process were lunched.
A file named 'stop-state' is 
saved containing the state of the simulation run. This option permits us to
stop the simulation and change some physical, numerical, or control
parameters such as
\fBtstop\fR, \fBeps\fR or \fBoptions\fR.

.SH EXAMPLES
The following command will run the code untill time is 2, with a time
step of 1/64, and every 8 time steps will write a snap,
.sp
gbsph nbody=4096 dtime=1/64 out=snap%03d
.sp

If you run
gbsph nbody=4096 dtime=1/64 out=snap%03d statefile=save
.sp
and if before it terminates normally it is interrupted you may resume the run with
.sp
gbsph restorefile=save

.SH SEE ALSO
nbody_n2(1)

.SH COPYRIGHT
Copyright (C) 1999-2007
.br
M.A. Rodriguez-Meza
.br
