't" t
.TH nbody_n2 1 "January 2009" UNIX "NagBody PROJECT"
.na
.nh   

.SH NAME
\fBnbody_n2\fR - (Direct force calculation) simulation of a 
selfgravitating N-body system

.SH SYNOPSIS
\fBnbody_n2\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ] 
.sp

.SH DESCRIPTION
\fBnbody_n2\fR - Simulates the evolution of a N-body system interacting with a potential
proportional to 1/r.
The force is computed directly. Code complexity is N*N.

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
nbody_n2 -help
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
is the name file with the values of the input parameters.
You may also input this filename by only writing:
.sp 
nbody_n2 parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice.
Comment lines start
with an "%". Follow each name option with a blank space and the option value.
The order of the option lines does not matter.
Also you may create an example input file
by executing
.sp
nbody_n2
.sp
This will run the \fBnbody_n2\fR code with default values generating
a test data set internally,
and when it finish you will have in your
running directory the file "parameters_null-nbody_n2-usedvalues".
Now you may edit this file to adapt it
to your own simulation parameters.
It may be helpful to change this file name to whatever apropriate.

.sp
Note: You may modify values of some parameters in the parameter input file. 
Just, in the command line,
after the parameter input file name write the name of the parameter, "=", 
and its value.
See description below of parameters and its possible values.

.IP "\fBeps\fR" 12
It is the smoothing force length. The interaction potential is proportional
to 1/Sqrt(r^2 + eps^2). In this way the singularity is avoided.

.IP "\fBnbody\fR" 12
It is the number of bodies to simulate.
.sp
Note: it is recommended that this number is given as power of 2.

.IP "\fBdtime\fR" 12
[a: dt] is the time step integration. It can be given as an integration
frequency just type in <dt1/dt2>.

.IP "\fBtstop\fR" 12
is the time to stop the simulation.

.IP "\fBseed\fR" 12
It is the Random number seed for the test run.

.IP "\fBicfile\fR" 12
[a: ic] you give here the name of the file with the N-body initial data.

.IP "\fBicfilefmt\fR" 12
[a: icfmt] is the format of the initial condition file. Values are:
<snap-bin> (Binary) or <snap-ascii> (ASCII).

.IP "\fBsnapout\fR" 12
[a: o] you give here the name structure for the output of N-body snaps.
The format follows
as the ones used in C-language for integers ("%0#d").

.IP "\fBsnapoutfmt\fR" 12
[a: ofmt] you tell the code the format of the snaps output. 
There are two options <snap-ascii>,
<snap-bin>, or <snap-pv>.
Brakets are not written.
The n-body standard format binary or ASCII is a file with 
n-body data written as follows:

nbody
.br
NDIM
.br
time
.br
mass (of all the particles)
.br
x y z (position for all the particles)
.br
vx vy vz (velocity for all the particles)
.br

And in the colummns format the particle data is written in column form as

# nbody NDIM time
.br
# nbody value NDIM value time value
.br
Id mass x y z vx vy vz (for all the particles)

.IP "\fBdtout\fR" 12
This is the time step for a snap output. It can be given as an output frequency, just
type in <dt1/dt2>.

.IP "\fBdtoutinfo\fR" 12
This is the time step output for other info such energy statistics.
It can be given as an output frequency, just
type in <dt1/dt2>.

.IP "\fBstatefile\fR" 12
[a: state] you give here the name of a file where the run state will be saved.
If it is null no run
state will be saved.

.IP "\fBstepState\fR" 12
you give here the number of time step integrations to jump to save a state
of the run.

.IP "\fBrestorefile\fR" 12
[a: restore] if it is not null a run will be restarted from the data stored in this file.

.IP "\fBoptions\fR" 12
[a: opt] you may give here various code behavior options.
They are, <reset-time> (inputdata); <out-phi> (outputdata); 
<out-acc> (outputdata).
If you save a state file, this parameter is saved.

.SH UNITS
If an initial condition is given units are defined according to the position,
velocity and mass of the bodies in the file given with \fBicfile\fR.
If \fBicfile\fR is NULL then position, velocity, and mass of each body is
genereted according to a Plummer model with units such that, G=-E=M=1, G the
gravitational constant, M the total mass of the Plummer sphere.

.SH STOPPING THE CODE
Once the \fBnbody_n2\fR is running you may always stop it by executing the command

echo > stop

You must be in the same directory were the process were lunched.
A file named 'stop-state' is 
saved containing the state of the simulation run. This option permits us to
stop the simulation and change some numerical parameters such as
\fBeps\fR or \fBoptions\fR.

.SH EXAMPLES
By executing the command,

.br
nbody_n2 nbody=4096 dtime=1/32 out=snap%03d

.br
will run the code, generating internally a Plummer sphere
sampled with 4096 particles and then evolving it up to t=2
with a time step of 1/32.
Snaps of the evolution will be saved as 
snap000, snap004, ..., snap064, with format <snap-ascii>.

.SH ANIMATIONS
You may run the \fBanalysis_grav\fR code to see animation plots.
Run a simulation using, for example, the command:

nbody_n2 tstop=10 o=snap > output &
.br
then, use the command:

analysis_grav in=snap analysis_type=snap-anim

to see the animation of the simulation.
The x and y ranges may change according to 
particle positions during simulation, then, can be useful to set
\fBxrange\fR=<xmin:xmax> and \fBxrange\fR=<ymin:ymax>.

.SH SEE ALSO
gbsph(1), analysis_grav(1), nplot2d(1)

.SH COPYRIGHT
Copyright (C) 1999-2010
.br
M.A. Rodriguez-Meza
.br
