't" t
.TH model 1 "January 2005" UNIX "NagBody PROJECT"
.na
.nh   

.SH NAME
\fBmodel\fR - Models generator for the simulation of a selfgravitating 
N-general body system

.SH SYNOPSIS
\fBmodel\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ] 
.sp

.SH DESCRIPTION
\fBmodel\fR - Generates the initial condition for N-Body system and/or N-sph body
which is giving to NagBody codes to simulate the evolution of these systems.

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
models -help
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
models parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. Comment lines start
with an "%". Follow each name option with a blank space and the option value.
The order of the option lines does not matter. Also you may create an example input file
by executing
.sp
models
.sp
This will run the \fBmodel\fR code with default values and when it finish you will have in your
running directory the file "parameters_null-usedvalues". Now you may edit this file to adapt
to your own simulation parameters. It may be helpful to change this file name to whatever apropriate.

.IP "\fBin\fR" 12
you give here the name of a file or files with N-body data. The data structure
of these files follows the snap formats included in the NagBody project.
With this option you can combine or
manipulate the IC data giving in these files. 
Or you may change the format of a single file. 

.IP "\fBinfmt\fR" 12
is the format of the files in the option \fBin\fR. They are <snap-ascii>
(the standard format) <snap-bin> (Binary), <snap-pv>
(in columns), <gadget11-bin>, <gadget11-ascii>, <gadget11-bin-double>, 
<heitmann-ascii>. For more than 10^6 bodies we have <gadget11-ascii-long> and
<heitmann-ascii-long>. 
Default value is <snap-ascii>.
If in option \fBin\fR you give more than one input files you must give in
option \fBinfmt\fR a snap format for each one.

.IP "\fBout\fR" 12
[a: o] you give here the name file for the output of N-body initial data. Default value is <snapout.dat>.

.IP "\fBoutfmt\fR" 12
[a: ofmt] you tell the code the format of the snap output. 
There are several options <snap-ascii>, <snap-bin>, <snap-pv>, <gadget11-sph-ascii>,
<gadget11-normal-body-ascii>,
<gadget11-ascii-reducido>,
<typsy-bin>, <nemo>, <powmes-ascii>, <powmes-ascii-long>.
The last format added is <snap-gadget>, the standard gadget format, particles are written as dark matter particles, gadget type 1.
Angular brackets are not written.  Default value is <snap-ascii>.
.sp
Note that when changing snap format
you may need to supply additional header information by hand or chose an
appropriate output format. For example, to change snap format of a file from
<snap-ascii> to a gadget format you may choose <gadget-sph-ascii>
or <snap-normal-body-ascii> as output format but no 
<gadget11-bin>, <gadget11-ascii>, <gadget11-bin-double> gadget formats.

.IP "\fBmodel-type\fR" 12
you give here the model to generate. Possible values are: 
.br
<plummer> generates a Plummer sphere. Units are G=M=-4E=1 
(see: S.J. Aarseth, M. Henon, and R. Wielen, Astro. & Astrophys. 37, 183 (1974)). 
This is the default value.
.sp
<plummer-finite> generates a Plummer sphere with a maximum radius (Rmax=2). Units are G=M=-4E=1.
.sp
<freefall-plummer> generates a free fall (zero particles velocities) Plummer sphere with a 
maximum radius (Rmax=2). Units are G=M=-4E=1.
.sp
<uniformsphere> generates a sphere with uniform density of matter. Units are G=M=R=1, M total
mass, R radius of the sphere. The particles have a uniform random velocity. The velocity of the 
center of mass of the sphere  may be given and also the sphere may be rotating with an angular velocity
omega0 around z-axis. If saved with outfmt=gadget-sph-ascii it will be useful to test the standard
model of a gas collapsing isothermically with gadget ver. 1.1.
The format in which the file is saved defines it as a gas or as a collisionless n-body system.
.sp
<pertisothermalsphere> generates a sphere with uniform mass density with a azimuthal
perturbation of amplitud a_p=0.1 and second harmonic m_p=2.
The format in which the file is saved defines it as a gas or as a collisionless n-body system.
If saved as a gas, this model is useful for studying the isothermal colapse of a perturbed
uniform sphere.
.sp
<uniformcubicbox> generates a cubic box with a uniform density of bodies. Units are such
that G=M=R=1. Other needed parameters are Rmax and absvel which determines the size of the 
system and the maximum random speed of the particles. The center of the box is at (0,0,0).
Again, the format in which the file is saved defines it as a gas or as a collisionless n-body system.
.sp
<gaussianspheroid> generates a gaussian spheroid with axes scales a_spheroid, 
b_spheroid, and c_spheroid. Other parameteres may be given such cmpos, cmvel, omega0, absvel.
The Bonnort-Ebert model is obtained with a_spheroid=b_spheroid=2, c_spheroid=1, Mtotal=40,
absvel=0, and omega0=1.

.IP "\fBoptions\fR" 12
you may give here various code behavior options.
For example, you may convert the mass data in a snap file.
Just give \fBoptions\fR=<mass-transform> and \fBfactor\fR=<value>.
This will multiply bodies mass by <value>. Other possibilities
are: \fBoptions\fR=<pos-transform> and \fBoptions\fR=<vel-transform>,
that transform position and velocities of particles respectively.

.IP "\fBnbody\fR" 12
it is the number of bodies to generate.
.sp
Note: it is recommended that this number is given as power of 2.

.IP "\fBMtotal\fR" 12
is the total mass of the system.

.IP "\fBRmax\fR" 12
is the maximum radius of the system.

.IP "\fBvcmx\fR" 12
is the x-component of the velocity of center of mass of the system.

.IP "\fBvcmy\fR" 12
is the y-component of the velocity of center of mass of the system.

.IP "\fBvcmz\fR" 12
is the z-component of the velocity of center of mass of the system.

.IP "\fBcmx\fR" 12
is the x-component of center of mass of the system.

.IP "\fBcmy\fR" 12
is the y-component of center of mass of the system.

.IP "\fBcmz\fR" 12
is the z-component of center of mass of the system.

.IP "\fBabsvel\fR" 12
is the maximum absolute velocity for the random motion of the bodies.

.IP "\fBomega0\fR" 12
is the angular velocity of the system around the z-axis
(omega0=0.698783 for isothermal collapse of a sphere).

.IP "\fBa_p\fR" 12
is the amplitud of the azimuthal perturbation (Useful isothermal collapse of a sphere).

.IP "\fBm_p\fR" 12
is the mode of the azimuthal perturbation (Useful isothermal collapse of a sphere).

.IP "\fBSoundSpeed\fR" 12
is the sound speed (Useful for isothermal collapse of a sphere).

.IP "\fBa_spheroid\fR" 12
is the spheroid x-axis (Bonnort-Ebert case default, also give Mtotal=40).
.IP "\fBb_spheroid\fR" 12
is the spheroid y-axis (Bonnort-Ebert case default).
.IP "\fBc_spheroid\fR" 12
is the spheroid z-axis (Bonnort-Ebert case default).

.IP "\fBfactor\fR" 12
is a factor to transform snap file data.

.IP "\fBseed\fR" 12
it is the Random number seed for the random generator needed to generate some of the
models. It can be any integer number, its 
default value is '-1'.

.SH EXAMPLES
The following command will run the code taking the file "snap0010" which is
in the standard snap format (ASCII) and will transform it to the gadget normal body format
.sp
model in=snap0010 out=snap0010.gdt outfmt=gadget11-normal-body-ascii
.sp
The command
.sp
model model_type=plummer
.sp
will generate a Plummer sphere with 4096 bodies.
.sp
The command
.sp
model model_type=uniformsphere out=isothermgassphere outfmt=gadget-sph-ascii
.sp
will generate an isothermal gas sphere to study the adiabatical collapse of a gas sphere.
The format of the output snap is apropriate for gadget (ASCII).
.sp
A variation of the previous command
.sp
model model_type=uniformsphere out=uniformsphere outfmt=gadget11-normal-body-ascii
.sp
will produce a 4096 bodies sphere of uniform mass density free-fall collapsing. The output
is apropriate for runing this initial condition with gadget 1.1.

.sp
The following command
.sp
model model_type=gaussianspheroid out=bonnort Mtotal=40 omega0=1
.sp
will generate the Bonnort-Ebert model.
The format of the output is standard snap.

.sp
The command
.sp
model model_type=pertisothermalsphere out=pertisosphere outfmt=gadget-sph-ascii omega0=0.698783
.sp
will generate the perturbed isothermal gas sphere model apropriate for studying
the collapse problem in stellar formation.
The format of the output is for running with pgdgt_sph.

.sp
The command
.sp
model in=snapone,snaptwo,snaptree infmt=snap-ascii,snap-bin,gadget-sph-ascii out=snapmixture outfmt=snap-ascii
.sp
will combine in the file snapmixture the three snaps, snapone, snaptwo and snapthree with
formats snap-ascii, snap-bin, and gadget-sph-ascii, respectively.
The format of the output is the standard snap-ascii format. 
The snaps snapone, snaptwo, and snapthree may be created with \fBmodel\fR
using the option 'model_type', and options to change cmpos or cmvel, among other
options.

.SH SEE ALSO
nbody_n2(1)

.SH VERSION
Actual version is 1.2.

.SH COPYRIGHT
Copyright (C) 1999-2007
.br
M.A. Rodriguez-Meza
.br
