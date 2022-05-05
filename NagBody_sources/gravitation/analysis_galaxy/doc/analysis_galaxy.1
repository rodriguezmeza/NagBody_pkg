't" t
.TH analysis_galaxy 1 "January 2006" UNIX "NagBody PROJECT"
.na
.nh   

.SH NAME
analysis_md - Code for data analysis of simulations of a MD system
.SH SYNOPSIS
\fBanalysis_galaxy\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ] 
.sp


.SH DESCRIPTION
analysis_galaxy - Analyses data from simulation of the evolution of a N-body MD system 
interacting with a potential like Lennard-Jones or other similar potentials.
It tooks a snap or
a sequence of snaps and process them according with the options below.


.SH OPTIONS
The options have the structure
.sp
\fIoption_name\fR = <option_value>
.sp
Options and their possible values are:

.IP "\fBparamfile\fR" 12
is the name file with the values of the input parameters. Overwrite parameters
values below. You may also input this filename by only writing:
.sp 
analysis_md parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. 
Comment lines start
with an "%". Follow each name option with a blank space and the option value.
The order of the option lines does not matter. Also you may create an example 
input file
by executing
.sp
analysis_galaxy
.sp
This will run the \fBanalysis_galaxy\fR code with default values and when it 
finish you will have in your
running directory the file "parameters_null-usedvalues". Now you may edit this 
file to adapt
to your own simulation parameters. It may be helpful to change this file name 
to whatever apropriate.

.IP "\fB-help\fR" 12
with this option you get the list of all parameters and their default values.

.IP "\fBanalysis_type\fR" 12
[a: at] you give here the type of analysis you want to do. 
First, we have the animations we can do with the internal visualization package.
.br 

.br
<snap-anim>
will give us an x-y animation of the snaps. With the option \fBusingcolums\fR below
you may see other projections. In the case of binary mixtures you may visualize
with diferent colours the individual species, you just use
\fBoptions\fR=<type-anim>.
.br

.br
<snap-animtrajectory>
will give us an x-y animation of the particle trajectories.
With the option \fBusingcolums\fR below
you may see other projections. 
.br

.br
<general-fx-anim>
This option will give you an animation of a general function tabulated in a
file or a system of files.

.br
<rhotheta-anim>
Gives an azimuthal profile animation of several quantities, like,
disc density, halo density, radial aceleration, tangential aceleration,
and tangential velocity.

.br
<vcr-anim>
Gives you an animation of the circular velocity vs distance to the center
of mass.

.br
<acf-anim>
Show you snaps of the autocorrelation function.

.br
<bodies-anim>
Gives you snaps of selected bodies. Bodies can be selected individually
or by groups. Use \fBbodiesID\fR and \fBbodiesSets\fR options.

.br
<bodies-animtrajectory> 
work in the same you as <bodies-anim>.

.br
<locate-bodiesid>
Locate bodies IDs given the bodiesSets, RhoDeltaZ, RhoR, and RhoRDeltaR.

.br
Secondly, we have the snap conversion. 

.br
<snap-conversion>
If you type in this option 
followed by the apropriate I/O options: \fBin\fR, \fBinfmt\fR, \fBout\fR,
\fBoutfmt\fR, (see below) you may convert the snap to 
other several formats.
.br

.br
<snap-lessnbody> with this option you may reduce the number of bodies that
represent a N-body system. This may be useful to process a smaller version
of a snap prior to make production runs or further processing.
.br

.br
<two-bdh-galaxies> with this option you may combine two galaxy models. Give 
the two file names and the two file formats if necesary in options \fBin\fR
and \fBinfmt\fR.

.br

.IP "\fBin\fR" 12
you give here the name structure for the input files of N-body snaps 
(they are the output of N-Body codes). The format follows
the one used in C-language for integers ("%0#d"). You may give here the whole
trajectory. For example, you may type in "in=snap%03d". When used in
combination with \fBoptions\fr=<two-bdh-galaxies> you should give two
file names for this option separated by a ",".

.IP "\fBinfmt\fR" 12
you tell the code the format of the snaps input. There are several options. 
First, we have the standard formats: <snap-ascii>, <snap-bin>, or <snap-pv>. 
Second, we have the gadget format: <gadget11-bin> (binary) and <gadget11-ascii> 
which are the standard snap format of Gadget 1.1. Other option included is
<gadget11-bin-double> which is the option <gadget11-bin> but in double precision.
<gadget11-bin-double-reducido> and <gadget11-ascii-reducido> read snaps in gadget
format but without the cosmological parameters. Third, we have other useful
formats like <heitmann-ascii> or the two components
Lennard-Jones snap formats:
<snap-blj-ascii>, <snap-blj-pv>, and <snap-blj-bin>.
As with optin \fBin\fR, when used in
combination with \fBoptions\fr=<two-bdh-galaxies> you should give two
file formats for this option separated by a ",".

.IP "\fBout\fR" 12
you give here the name structure for the output files of N-body snaps. 
The format follows the option \fBin\fR above.

.IP "\fBoutfmt\fR" 12
you tell the code the format of the snaps output. There are several options. 
First, we have the standard formats (see \fBinfmt\fR):
<snap-ascii>, <snap-bin>, or <snap-pv>. 
Second, we have the gadget formats: <gadget11-ascii>, <gadget11-bin-double>,
<gadget-body-ascii>, and 
<gadget-sph-acii>. <gadget11-bin-double-reducido> and <gadget11-ascii-reducido>
write snaps in gadget format but without the cosmological parameters. In the case
of formats <gadget11-ascii> and <gadget11-bin-double> you have to provide the
distribution of body types and may be the mass tables, i.e., N_gas, N_halo, and
so on.
.br
The combination of I/O options are very useful for snap format conversion.

.IP "\fBbasedir\fR" 12
Base directory for output files.

.IP "\fBisnap\fR" 12
Initial snap number of the sequence snaps you want to procces.

.IP "\fBfsnap\fR" 12
Final snap number of the sequence of snaps.

.IP "\fBxrange\fR" 12
When plotting this option stablishes the range in x-axis.

.IP "\fByrange\fR" 12
When plotting this option stablishes the range in y-axis.

.IP "\fBzrange\fR" 12
When plotting in 3D this option stablishes the range in z-axis.

.IP "\fBusingcolumns\fR" 12
You tell the code what columns of the data file use to analyse.

.IP "\fBusingrows\fR" 12
You tell the code what rows of the data file use to analyse. This option
should be used in combination with <analysis_type=thermo_avg>.

.IP "\fBsizeHistRho\fR" 12
This option set the size of the rho historgram. This option must be use in
combination with <datanaly_type=rho_anim> or <datanaly_type=vc_anim>.

.IP "\fBrangeR\fR" 12
This option stablishes the range in r values.
This option must be use in combination with <datanaly_type=rho_anim>
or <datanaly_type=vc_anim>.

.IP "\fBdev\fR" 12
When plotting with the internal visualization package this option set the
output device.

.IP "\fBori\fR" 12
When plotting with the internal visualization package this option set the
the orientation of the display. When sending the output to a postscriipt
device, ori=1, produce an output in portrait form.

.IP "\fBa\fR" 12
When plotting with the internal visualization package this option set the
the aspect ratio of the display.

.IP "\fBgeo\fR" 12
When plotting with the internal visualization package this option set the
the window size in pixels.

.IP "\fBbg\fR" 12
When plotting with the internal visualization package this option set the
background color.

.IP "\fBncol0\fR" 12
When plotting with the internal visualization package this option set the
number of colors to allocate in cmap 0 (upper bound).

.IP "\fBncol1\fR" 12
When plotting with the internal visualization package this option set the
number of colors to allocate in cmap 1 (upper bound).


.IP "\fBoptions\fR" 12
you may give here various code behavior options.
.br

.br
<save>
In the case of animations, you may type in for this option <save> to
tell the code that it must save the images it produces to files instead of
to sended to the graphics device.
.br

.br
<rho-cgm>
If you type in <rho-cgm> in combination of \fBanalysis_type\fR=<rho-anim>,
the code will transform each snap file into a graphics file (CGM format). Gnuplot
is needed and also a gnuplot script file named "rho_cgm.gnu".


.SH EXAMPLES
The following command will produce a single file named "rhoaxes.dat" with
the density distribution of matter along the three axes for every
snap file created in the simulation (in this case with the snaps are
written in the format <snap-blj-ascii>,
.sp
analysis_md in=snap%03d infmt=snap-blj-ascii analysis_type=rhoaxes-anim
.sp
The command:
.sp
analysis_md in=snap%03d analysis_type=rho-gnuanim options=rho-cgm
.sp
will produce a graphic file (CGM format) of the distribution of particles
along the axes for every
snap of the simulation.

.sp

.sp
The following command will reduce the snap "indata" in "heitmann-ascii" format
to an initial factor of 1/8 and will produce the file "outdata" with "snap-pv"
format
.sp
analysis_md analysis_type=snap-lessnbody in=indata infmt=heitmann-ascii 
out=outdata outfmt=snap-pv fsnap=0 reductionFac=0.125

.SH ANIMATIONS WITH GNUPLOT
You may run gnuplot to see animation plots. Make the following script:
.sp
a=a+1 
.br
plot "snap.dat" w points pointtype 6 
.br
pause 1 
.br
if(a<500000) reread
.sp
Save the script in a file named "snap_anim.gnu". 
At the same time that datanaly_grav code is running with the command
.sp
analysis_md in=snap%03d analysis_type=snap-gnuanim
.sp
start 
gnuplot and execute the commands:
.sp
a=0
.br
load "snap_anim.gnu"
.sp
Then you will be able to see an animation plot of the snaps.
You may also write the following script to see the animation of the
radial distribution function,
.sp
a=a+1
.br
plot "rdf.dat" w l
.br
pause 1
.br
if(a<50000) reread
.br

.br
Note: in case of Mac OS X you need to run in X-Window session.


.SH SEE ALSO
md_lj_n2(1), md_lj_tree(1)
.SH COPYRIGHT
Copyright (C) 1999-2006
.br
M.A. Rodriguez-Meza
.br
