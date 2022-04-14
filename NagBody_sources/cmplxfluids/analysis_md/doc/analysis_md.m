't" t
.TH analysis_md 1 "January 2012" UNIX "NagBody PROJECT"
.na
.nh   

.SH NAME
\fBanalysis_md\fR - Code for data analysis and visualization of simulations of a 
MD binary system
.SH SYNOPSIS
\fBanalysis_md\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ] 
.sp


.SH DESCRIPTION
\fBanalysis_md\fR - Analyses data from simulation of the evolution of a 
N-body MD binary system 
interacting with a potential like Lennard-Jones or other similar potentials.
It tooks a snap or
a sequence of snaps and process them according with the options below.


.SH OPTIONS
The options have the structure
.sp
\fIoption_name\fR = <option_value>

.sp
except option \fB-help\fR.

.IP "\fB-help\fR" 12
By writting

analysis_md -help

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
analysis_md parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. Comment lines start
with an "%". Follow each name option with a blank space and the option value.
The order of the option lines does not matter. Also you may create an example input file
by executing
.sp
analysis_md analysis_type=snap-anim
.sp
This will run the \fBanalysis_md\fR code with default values, display nothing,
you just hit 'return' and
you will have in your running directory the file "parameters_null-analysis_md".
Now you may edit this file to adapt to your own analysis parameters. 
It may be helpful to change this file name to whatever apropriate.

.sp
Note: You may modify values of some parameters in the parameter input file. 
Just, in the command line,
after the parameter input file name write the name of the parameter, "=", 
and its value.
See description below of parameters and its possible values.

.IP "\fBanalysis_type\fR" 12
[a: at] you give here the type of analysis you want to do. 
First, we have the animations we can do with the internal visualization package.
.br 

.br
<snap-anim>
will give us an x-y animation of the snaps. With the option \fBusingcolums\fR below
you may see other projections (for the meaning of columns see file 'snap-analysis.dat').
In the case of binary mixtures you may visualize
with diferent colours the individual species, you just use
\fBoptions\fR=<type>. Also, it could be useful to change the size of the symbol,
using \fBsymbolsize\fR.
.br

.br
<snap-animtrajectory>
will give us an x-y animation of the particle trajectories.
With the option \fBusingcolums\fR below
you may see other projections. With \fBoptions\fR=<type> it will use different
colors for each species.  With \fBoptions\fR=<save-to-file> will save
into the file 'snap-analysis.dat'
all snaps info used to show a trajectory animation.
With \fBoptions\fR=<save> will save the last frame to a file whose name
is in \fBout\fR and with extension given in \fBoutfmt\fR.
With \fBoptions\fR=<save-animation> will save all frames in a trajectory
animation, their filenames are given in the \fBout\fR pattern and with
extension given in \fBoutfmt\fR, you may also use the \fBbasedir\fR option
to give a directory base. Note that you may combine in \fBoptions\fR
more than one value or case.
.br

.br
<snap-anim3d>
will give us a 3D view of the snap animations. 

.br
Note: If the system does not form
a cube, it will be necesary to set \fBxrange\fR, \fByrange\fR, and \fBzrange\fR
in such a way that xyz-axis form a cube containing the system. See an example
in the section EXAMPLES below.
.br

.br
<snap-animtrajectory3d>
will give us a 3D animation of the particle trajectories. See note above.
.br

.br
<slab-anim>
will give us a slab animation of the particles. Use this animation type in combination with
\fBxmin\fR, \fBxmax\fR, \fBymin\fR, \fBymax\fR, \fBzmin\fR, and \fBzmax\fR options, to define the slab geometry.
See also <snap-anim> animation type for more details.
Note that this animation type is only defined in 3D molecular dynamics animations.

.br

.br
<slab-animtrajectory>
will give us a slab animation of the particle trajectories. See <slab-anim> above.
Note that this animation type is only defined in 3D molecular dynamics animations.
.br

.br
<bodies-anim>
will give us a snap animations of selected bodies using the option
\fBbodiesID\fR.
Also, note that you must set \fBxrange\fR and \fByrange\fR to
see properly the animation.
.br

.br
<bodies-animtrajectory>
will give us a bodies animation of the particle trajectories of selected
bodies using the option \fBbodiesID\fR. As above, 
note that you must set \fBxrange\fR and \fByrange\fR to
see properly the animation. 
.br

.br
<rhoaxes-anim>
will give you an animation of the density distribution function along the axes.
This option will
produce a file named "snap-analysis.dat" with rho vs x, y or z;
rho1 vs x, y or z; and rho2 vs x, y or z. This file can be used to
improve your plot with \fBnplot2d\fR. 
Use options \fBsizeHist\fR and \fBstepAvg\fR described below
to adapt to your needs.

.br
<nfrecaxes-anim>
will give you an animation of the frequency of the density distribution
function along the axes.
This option will
produce a file named "rhoaxes.dat" with rho vs x, y or z.
Use options \fBsizeHist\fR and \fBstepAvg\fR described below
to adapt to your needs.

.br
<vel-anim>
This option will give you an animation of the distribution of particles speed.
Use options \fBsizeHist, \fBrangeVal\fR, and \fBstepAvg\fR
described below to adapt to your needs.

.br
<rdf-anim>
This option will give you an animation of the radial distribution function.
Use options \fBsizeHist, \fBrangeVal\fR, and \fBstepAvg\fR
described below to adapt to your needs. Included here are experimental data
of neutron dispersion by an Argon system. If you wish to include experimental
data give \fBoptions\fR=<experiments> and note that you should have in
the working directory a file with experimental data named "gdr_exp.txt".

.br
<general-fx-anim>
This option will give you an animation of a general function tabulated in a
file or a system of files.

.br
Secondly, we have the snap conversion. 

.br
<snap-conversion>
If you type in this option 
followed by the apropriate I/O options: \fBin\fR, \fBinfmt\fR, \fBout\fR,
\fBoutfmt\fR, (see below) you may convert snaps to 
other several formats. Note that the default snap format for \fBin\fR
is <snap-blj-ascii>. Notice also that you may allways convert snaps in
the \fBmd_blj_tree\fR code natural format that could
be binary, ascii or pv (see man pages of \fBmd_blj_tree\fR)
to other formats but not from an arbitrary format to the natural 
\fBmd_blj_tree\fR format.
.br

.br
<snap-lessnbody> with this option you may reduce the number of bodies that
represent an N-body system. This may be useful to process a smaller version
of a snap prior to make production runs or further processing. Note that the
format of the output is <snap-ascii> given that we have lost the header
of the input file. So, in order to visualize it or further processing it, you
should use other analysis or visualization code like, for example,
\fBanalysis_grav\fR, if it is available. Or, you may by hand edit the snap
output file and create an apropriate header in such a way that you end with
a snap in the \fBmd_blj_tree\fR ascii format. Therefore, you should save the snap
using \fBofmt\fR=<snap-blj-ascii>, then edit the file fixing values
of nbody and masses values for types '1' and '2'.
.br

.br
And thirdly, we have the option to compute averages of several observables,
and
conversion of units:

.br
<thermo-avg> will give you the average of the observables in the file
thermo.dat (the name may have been changed and have to be given with
the option \fBin\fR). With options \fBusingcolumns\fR
and \fBusingrows\fR you change the particular observables and the rows
you want to average. The option \fBoptions\fR=<overwrite> creates a new
file (given by the option \fBout\fR)
where the averages are saved, without this value the code appends
the averages to an existing file. Also, the option \fBoptions\fR=<errorbar>
will create two additional columns, with the error limits of the observable.
.br

.br
<units-conversion> will give you the several tipical scale physical quantities,
like, the tipical or unit of length, energy, time, mass, etc. The option
\fBunitsset\fR establishes the units set used to do the conversion 
in the form of <unit-length:unit-mass:unit-energy>. Default values are for
the atom Ar in MKS.
.br

.IP "\fBreductionFac\fR" 12
Reduction factor that must be used in combination with
\fBanalysis_type\fR=<snap-lessnbody>.

.IP "\fBbodiesID\fR" 12
[a: bid]
Bodies IDs to use with \fBanalysis_type\fR=<bodies-anim> and
\fBanalysis_type\fR=<bodies-animtrajectory>. The format
of this option is <b1,b2,...>. The limit in the number of bodies to
show is 20.

.IP "\fBin\fR" 12
you give here the name structure for the input files of N-body snaps 
(they are the output of N-Body codes). The format follows
the one used in C-language for integers ("%0#d"). You may give here the whole
trajectory. For example, you may type in "in=snap%03d".

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

.IP "\fBout\fR" 12
[a: o] you give here the name structure for the output files of N-body snaps. 
The format follows the option \fBin\fR above.

.IP "\fBoutfmt\fR" 12
[a: ofmt] you tell the code the format of the snaps output. There are several options. 
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
[a: xr] When plotting this option stablishes the range in x-axis.

.IP "\fByrange\fR" 12
[a: yr] When plotting this option stablishes the range in y-axis.

.IP "\fBzrange\fR" 12
When plotting in 3D this option stablishes the range in z-axis.

.IP "\fBxmin\fR" 12
When using <slab-anim> animation type option, \fBxmin\fR, \fBxmax\fR, \fBymin\fR, \fBymax\fR, 
\fBzmin\fR, and \fBzmax\fR, define the slab geometry. 
Note that this option is only defined in 3D molecular dynamics animations.

.IP "\fBxmax\fR" 12
When using <slab-anim> animation type option, \fBxmin\fR, \fBxmax\fR, \fBymin\fR, \fBymax\fR, 
\fBzmin\fR, and \fBzmax\fR, define the slab geometry.
Note that this option is only defined in 3D molecular dynamics animations.

.IP "\fBymin\fR" 12
When using <slab-anim> animation type option, \fBxmin\fR, \fBxmax\fR, \fBymin\fR, \fBymax\fR, 
\fBzmin\fR, and \fBzmax\fR, define the slab geometry.
Note that this option is only defined in 3D molecular dynamics animations.

.IP "\fBymax\fR" 12
When using <slab-anim> animation type option, \fBxmin\fR, \fBxmax\fR, \fBymin\fR, \fBymax\fR, 
\fBzmin\fR, and \fBzmax\fR, define the slab geometry.
Note that this option is only defined in 3D molecular dynamics animations.

.IP "\fBzmin\fR" 12
When using <slab-anim> animation type option, \fBxmin\fR, \fBxmax\fR, \fBymin\fR, \fBymax\fR, 
\fBzmin\fR, and \fBzmax\fR, define the slab geometry.
Note that this option is only defined in 3D molecular dynamics animations.

.IP "\fBzmax\fR" 12
When using <slab-anim> animation type option, \fBxmin\fR, \fBxmax\fR, \fBymin\fR, \fBymax\fR, 
\fBzmin\fR, and \fBzmax\fR, define the slab geometry.
Note that this option is only defined in 3D molecular dynamics animations.

.IP "\fBusingcolumns\fR" 12
[a: uc] You tell the code what columns of the data file use to analyse.

.IP "\fBusingrows\fR" 12
[a: ur] You tell the code what rows of the data file use to analyse. This option
should be used in combination with <analysis_type=thermo_avg>.

.IP "\fBxlabel\fR" 12
[a: xl]
X-axis label.

.IP "\fBylabel\fR" 12
[a: yl]
Y-axis label.

.IP "\fBplotlabel\fR" 12
[a: pl]
Plot label.

.IP "\fBplotjoined\fR" 12
[a: pj]
Plot with points joined with a line.

.IP "\fBwithdots\fR" 12
[a: wd]
Plot with dots. Set \fBplotjoined\fR=<false>.

.IP "\fBwithsymbols\fR" 12
[a: ws]

.IP "\fBstepAvg\fR" 12
Sets the number of steps jumps to save an histogram of
RhoAxes, NFrecAxes, Vel, and RDF.
This option must be use in combination with \fBanalysis_type\fR equal to
<rhoaxes-anim>, <nfrecaxes-anim> <rdf-anim> or <=vel-anim>.

.IP "\fBsizeHist\fR" 12
This option set the size of the historgram for RhoAxes, NFrecAxes, Vel, and RDF
analysis.
This option must be use in combination with \fBanalysis_type\fR equal to
<rhoaxes-anim>, <nfrecaxes-anim> <rdf-anim> or <=vel-anim>.

.IP "\fBrangeVal\fR" 12
This option stablishes the range values for
Vel, and RDF analysis.
This option must be use in combination with \fBanalysis_type\fR equal to
<rdf-anim> or <=vel-anim>.

.IP "\fBlabelfontsize\fR" 12
Font size for labels. It is a real value.

.IP "\fBlabelfontweight\fR" 12
Font weight for labels. It is an integer value.

.IP "\fBnlsize\fR" 12
Font size for tick labels. It is a real value.

.IP "\fBlinewidth\fR" 12
Line width. It is an integer value.

.IP "\fBaxeswidth\fR" 12
Axes width. It is an integer value.

.IP "\fBsymbolweight\fR" 12
Symbol weight. It is an integer value.

.IP "\fBsymbolcolor\fR" 12
Symbol color. It is an integer value.

.IP "\fBsymbolsize\fR" 12
Symbol size. It is a real value.

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
[a: opt] you may give here various code behavior options. Above with show
several cases. Here we just remaind you two cases.
.br

.br
<save>
In the case of animations, you may type in for this option <save> to
tell the code that it must save the images it produces to files instead of
to sended to the graphics device.
.br

.br
<type>
In the case of <snap-anim> you may use this value to show plots with
bodies colored according to their species.

.SH I/O FILES
\fBanalysis_md\fR almost always produces a file named "snap-analysis.dat".
This file has in column form the data needed to produced the corresponding
plot. For example, if you are using \fBanalysis_type\fR=<snap-anim>, the data
saved in this file is bodies spatial positions (3 values in 3D), bodies
velocities (3 values in 3D), body ID, and body type.

.SH EXAMPLES
The following command will produce a single file named "rhoaxes.dat" with
the density distribution of matter along the three axes for every
snap file created in the simulation,
.sp
analysis_md in=snap%03d analysis_type=rhoaxes-anim
.sp

.sp
Let us assume that you lauched a simulation with the following command:
.sp
md_blj icModel=2 tstop=100 d=0.8 o=snap > output &
.sp
Therefore, the simulation is running in the background. The simulation box, with
Lx=8.4653, Ly=8.4653 and Lz=17.8618 is a parallepiped.
Now, type in the
command:
analysis_md in=snap xr=-9:9 yr=-9:9 zr=-9:9 geo=800x600 xl=x yl=y zl=z
fsnap=20 at=snap-anim3d ss=0.5
.sp
that will display a animation of the system in 3D. Try the same command without
options \fBxrange\fR, \fByrange\fR,  \fBzrange\fR set as before, instead use
default values.

.SH MOVIES FROM SNAPS
Once you are satisfied how frames are displayed you may produce eps files
that latter can be unified into a single movie file with a third party software
(In case of Mac OS X use "GraphicConverter"). Just use the same command line
instruction you gave adding \fBoptions\fR=<save>, and setting properlly 
options \fBdev\fR, \fBout\fR, \fBoutfmt\fR, and \fBori\fR.

.br
Note: in case of Mac OS X you need to run in X-Window session.


.SH SEE ALSO
md_blj(1), md_blj_n2(1), nplot2d(1)
.SH COPYRIGHT
Copyright (C) 1999-2012
.br
M.A. Rodriguez-Meza
.br
