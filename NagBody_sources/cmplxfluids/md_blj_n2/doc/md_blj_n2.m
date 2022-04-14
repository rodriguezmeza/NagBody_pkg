't" t
.TH md_blj_n2 1 "January 2011" UNIX "NagBody PROJECT"
.na
.nh   
.SH NAME
\fBmd_blj_n2\fR - simulation of a 
binary Lennard-Jones mixture liquid

.SH SYNOPSIS
\fBmd_blj_n2\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ] 
.sp

.SH DESCRIPTION
\fBmd_blj_n2\fR - Simulates the evolution of a binary N-body system interacting with 
a Lennard-Jones potential or a general potential.
The force is computed using the direct scheme. 
Code complexity in this case 
is N*N. 
In general periodic boundary condition is used. However, the code permits
to simulate cavities or other finite systems.
Two ensambles are used, one is the NVT ensamble (canonical) other is the
NVE (microcanonical).

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
md_blj_n2 -help
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
md_blj_n2 parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. 
Comment lines start
with an "%". Follow each name option with a blank space and the option value.
The order of the option lines does not matter. Also you may create an example 
input file
by executing
.sp
md_blj_n2
.sp
This will run the \fBmd_blj_n2\fR code with default values and when it finish 
you will have in your
running directory the file "parameters_null-md-usedvalues". 
Now you may edit this file to adapt it
to your own simulation parameters. It may be helpful to change this file name 
to whatever apropriate.
.sp
Note: You may modify values of some parameters in the parameter input file. 
Just, in the command line,
after the parameter input file name write the name of the parameter, "=", 
and its value.
See description below of parameters and its possible values.

.IP "\fBforcecalc_method\fR" 12
[a: fm] is the force calculation method to use. We have implemented four 
methods: direct <direct>, 
nearest neighborh list <nblist>, cells method <cells>, 
and tree walk method <normal>. The default value is <cells>.

.IP "\fBpotentialType\fR" 12
[a: pT] This option is to tell the code the force model. Values are <0> for
a Lennard-Jones force, <1> a shifted Lennard-Jones force, or <2> for a force
model given as a table in a file whose name is introduced with option 
\fBfnamePot\fR, see below.

.IP "\fBfnamePot\fR" 12
When use option \fBpotentialType\fR=<2> you should give with this option
a file name where
potential and force values table are. This file must have a header 
consisting of four lines, each one starting with "#". The first line should
be "# nforcepot ri rf dr". Second line should be 
"# sigma11 sigma12 sigma22 eps11 eps12 eps22 fphi11 fphi12 fphi22 fa11 fa12 fa22".
The third lines must have the values of first line parameters starting with "#".
Fourth line must have the values of second line parameters starting with "#".
.sp
After the header the force and potential info is in a table with 8 columns.
First column is the position id, second column is the magnitude of the distance
between two particles. Third and fourth columns are the values of the potential
and magnitude of the acceleration as if the particles were of type "1". Fifth
and sixth is for the case as if particles were type "1" and "2". Seventh and
eighth as if particles were of type "2".

.IP "\fBeps11\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the depth of the 
potential when interacting particles types are '1' and '1'. 
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBeps12\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the depth of the 
potential when interacting particles types are '1' and '2'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBeps22\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the depth of the 
potential when interacting particles types are '2' and '2'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBsigma11\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the scale length of the 
potential when interacting particles types are '1' and '1'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBsigma12\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the scale length of the 
potential when interacting particles types are '1' and '2'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBsigma22\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the scale length of the 
potential when interacting particles types are '2' and '2'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBRcut11\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the cut radius of the 
potential when interacting particles types are '1' and '1'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBRcut12\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the cut radius of the 
potential when interacting particles types are '1' and '2'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBRcut22\fR" 12
If \fBpotentialType\fR=<0> or <1> this parameter gives the cut radius of the 
potential when interacting particles types are '2' and '2'.
See section 'UNITS' below for the apropriate values that can be given.

.IP "\fBtheta\fR" 12
is a force acuracy parameter. Not in use by now. It will be useful for long
range forces.

.IP "\fBusequad\fR" 12
is the option to include quadrupoles in the force calculation. 
Not in use by now. It will be useful for long
range forces.

.IP "\fBtemperature\fR" 12
[a: t] is the temperature of the simulation in units such that eps/kB is the 
unit of temperature (see below section on "UNITS").

.IP "\fBdensity\fR" 12
[a: d] is the density of the liquid in code units.

.IP "\fBadjustTemperature\fR" 12
Canonical ensamble will be used if true. If false microcanical ensamble will 
be used.

.IP "\fBstepAdjustTemperature\fR" 12
Number of steps to adjut temperature. Only active if adjustTemperature is true.

.IP "\fBadjustCenterOfMass\fR" 12
If true the center of mass of the system is adjusted. This is useful when
computing the bulk viscosity.

.IP "\fBstepEquil\fR" 12
is the number of step to begin equilibrium computations.
Its default value is 100 time steps.

.IP "\fBstepAvg\fR" 12
is the number of steps to average properties.
Can be set to a new value when a run is restored.

.IP "\fBcomputeRhoAxes\fR" 12
Compute density profile.

.IP "\fBstepRhoAxes\fR" 12
Number of step jumps to save a rho_axes histogram.

.IP "\fBstepAvgRhoAxes\fR" 12
Number of rho_axes histograms to average.

.IP "\fBsizeHistRhoAxes\fR" 12
array size for rho_axes histogram.

.IP "\fBcomputeNFrecAxes\fR" 12
Compute frequency distribution of density profile.

.IP "\fBstepNFrecAxes\fR" 12
Number of step jumps to save a nfrec_axes histogram.

.IP "\fBstepAvgNFrecAxes\fR" 12
Number of nfrec_axes histograms to average.

.IP "\fBsizeHistNFrecAxes\fR" 12
array size for nfrec_axes histogram.

.IP "\fBcomputeVelDist\fR" 12
If <true> the code compute the histogram of the magnitud of velocity of the
particles and it is saved in the file 'vel.dat'. Default value is <false>.

.IP "\fBstepVel\fR" 12
is the number of steps jumps to save a velocity histogram.

.IP "\fBstepAvgVel\fR" 12
is the number of velocity histograms to average. Its default value is 4.

.IP "\fBsizeHistVel\fR" 12
is the array size for the velocity histogram. Its default value is 50.

.IP "\fBrangeVel\fR" 12
is the range of velocities for the histogram computation.

.IP "\fBcomputeRdf\fR" 12
If <true> the code compute the radial distribution function
of the
particles and it is saved in the file 'rdf.dat'. Default value is <false>.

.IP "\fBstepRdf\fR" 12
is the number of steps jumps to save a radial distribution function (RDF) histogram.

.IP "\fBstepAvgRdf\fR" 12
is the number of RDF histograms to average. Its default value is 20.

.IP "\fBsizeHistRdf\fR" 12
is the array size for the RDF histogram. Its default value is 200.

.IP "\fBrangeRdf\fR" 12
is the range of positions for the RDF histogram computation.


.IP "\fBcomputePressAxes\fR" 12
Compute pressure profile.

.IP "\fBstepPressAxes\fR" 12
Number of step jumps to save a pressure profile histogram.

.IP "\fBstepAvgPressAxes\fR" 12
Number of pressure profile measurements to average.

.IP "\fBsizeHistPressAxes\fR" 12
array size for pressure profile histogram.


.IP "\fBcomputeChemPot\fR" 12
Compute pressure profile.

.IP "\fBstepChemPot\fR" 12
Number of step jumps to save a pressure profile histogram.

.IP "\fBstepAvgChemPot\fR" 12
Number of pressure profile measurements to average.

.IP "\fBsizeHistChemPot\fR" 12
array size for pressure profile histogram.

.IP "\fBnumTestBodies\fR" 12
number of test bodies to make a ChemPot measurement.


.IP "\fBcomputeDiffusion\fR" 12
Compute diffusion coefficient.

.IP "\fBstepDiffuse\fR" 12
Number of step jumps to save a diffusion measurement.

.IP "\fBstepAvgDiffuse\fR" 12
Number of diffusion coefficient measurements to average.

.IP "\fBnBuffDiffuse\fR" 12
size of buffer to save a diffusion coefficient measurements.

.IP "\fBnValDiffuse\fR" 12
number of values to save of diffusion coefficient measurement.


.IP "\fBcomputeVelAcf\fR" 12
Compute velocity autocorrelation function.

.IP "\fBcomputeBulkViscosity\fR" 12
Compute bulk viscosity.

.IP "\fBcomputeTransport\fR" 12
Compute transport properties.


.IP "\fBstepAcf\fR" 12
Number of step jumps to save an Acf measurement.

.IP "\fBstepAvgAcf\fR" 12
Number of Acf measurements to average.

.IP "\fBnBuffAcf\fR" 12
size of buffer to save a Acf measurements.

.IP "\fBnValAcf\fR" 12
number of values to save of Acf measurement.


.IP "\fBcomputeSTCorr\fR" 12
Compute space-time correlations.

.IP "\fBstepCorr\fR" 12
Number of step jumps to save a STCorr measurement.

.IP "\fBstepAvgCorr\fR" 12
Number of STCorr measurements to average.

.IP "\fBnBuffCorr\fR" 12
size of buffer to save a STCorr measurements.

.IP "\fBnValCorr\fR" 12
number of values to save of STCorr measurement.

.IP "\fBnFunCorr\fR" 12
number of values to save of STCorr measurement.


.IP "\fBnbody\fR" 12
it is the number of bodies to simulate. If icfile option is null the code will generate
a initial condition (a test run) were all the particles are distributed uniformly in a cubic box
with gaussian random velocities.
Therefore, nbody will be of the form n^3.

.IP "\fBdtime\fR" 12
[a: dt]
is the time step integration in code units. Can also be given in the form of p/q,
where q is expressed as powers of 2 so we can think in terms of integration frequency.

.IP "\fBtstop\fR" 12
is the time to stop the simulation.

.IP "\fBintMethod\fR" 12
Integration method. It can be <0> or <1> for fixed time step given by \fBdtime\fR.
Or can be <2> for an adaptive time step method, the value of \fBdtime\fR is
the maximum allowed value.

.IP "\fBicModel\fR" 12
gives the initial condition model to construct. Possible values are: 
<1>, <2>, <3>, <4>, <5>, and <6>. Default is <2>. Option value <6> 
generate an FCC structure, see EXAMPLES section for its use.

.IP "\fBseed\fR" 12
it is the random number seed for generating the initial condition model.

.IP "\fBunitCells\fR" 12
[a: uC]
Number of unit cells along axes apropriate for \fBicModel\fR=<3>.

.IP "\fBnbodyprop\fR" 12
gives the number of bodies of the mixture for initial condition model
constructed using the option \fBicModel\fR. Must be given as <nbody1/nbody2>
format.

.IP "\fBmassprop\fR" 12
Masses of the spieces for generating initial condition data using the 
option \fBicModel\fR, in the format <mass1/mass2>.

.IP "\fBLxLyprop\fR" 12
Base sides of the parallelepiped initial condition generated using
\fBicModel\fR=<2>, in the format <Lx/Ly>.

.IP "\fBicfile\fR" 12
[a: ic]
you give here the name of the file with the N-body initial data.

.IP "\fBicfilefmt\fR" 12
[a: icfmt]
is the format of the initial condition file: 'snap-blj-bin' (binary)
or 'snap-blj-ascii' (ASCII) or 'snap-blj-pv'. 
See \fBsnapoutfmt\fR option below for a description of the file formats.

.IP "\fBsnapout\fR" 12
[a: o]
you give here the name structure for the output of N-body snaps. The format follows
as the ones used in C-language for integers ("%0#d"). snaps will be written at
time steps according to \fBdtout\fR.

.IP "\fBsnapoutfmt\fR" 12
[a: ofmt]
you tell the code the format of the snaps output. There are three options: 
the ASCII
snap n-body format (snap-blj-ascii); the binary n-body format 
(snap-blj-bin);
and the columns format
(snap-blj-pv). Default is <snap-blj-ascii>.
The n-body format binary or ASCII is a file with 
n-body data written as follows:

nbody
.br
nbody1
.br
nbody2
.br
NDIM
.br
time
.br
temperature
.br
density
.br
mass1
.br
mass2
.br
Lx
.br
Ly
.br
Lz
.br
eps11
.br
eps12
.br
eps22
.br
sigma11
.br
sigma12
.br
sigma22
.br
Rcut11
.br
Rcut12
.br
Rcut22
.br
ID (of all the particles)
.br
Type (of all the particles)
.br
mass (of all the particles)
.br
x y z (position for all the particles)
.br
vx vy vz (velocity for all the particles)
.br
phi (of all the particles, optionally)
.br
ax ay az (of all the particles, optionally)
.br

And in the colummns format the particle data is written in column form as

# nbody nbody1 nbody2 NDIM time temperature density masses Ls epss sigmas Rcuts
.br
# nbody nbody1 nbody2 NDIM time temperature density masses Ls epss sigmas Rcuts (the values)
.br
Id Type mass x y z vx vy vz phi (optionally) ax ay az (optionally) (for all the particles)


.IP "\fBdtout\fR" 12
this is the time for output a snap file. The out files will be written every dtout time.
Can also be given in the form of p/q,
where q is expressed as powers of 2.

.IP "\fBdtoutinfo\fR" 12
this is the time for output in the stadar output (stdout). The output to the stdout
will be written every dtoutinfo time. Can also be given in the form of p/q,
where q is expressed as powers of 2. This option is useful for controlling the cpu time
consumed by output processing.

.IP "\fBstatefile\fR" 12
[state] you give here the name of a file where the run state will be saved. If it is null no run
state will be saved.

.IP "\fBstepState\fR" 12
you give here the number of time step integrations to jump to save a state
of the run.

.IP "\fBrestorefile\fR" 12
[restore] if it is not null a run will be restarted from the data stored in this file.

.IP "\fBoptions\fR" 12
[opt] you may give here various code behavior options. 
They are, <reset-time> (inputdata); 
<out-phi> (outputdata); <out-acc> (outputdata). Other values are:
<bh86> and <sw94> useful to test opening cells when using the tree scheme
force calcultion method.
If you save a state file, this parameter is saved.
Several values can be given separated by coma, \fBoptions\fR=<option1,option2,...>.
Can be set to a new values when a run is restored.

.SH UNITS
Units are such that eps=sigma=m=1, where eps is the potential depth of the Lennard-Jones potential,
sigma the cut radius, and m is the mass of each particle. 
We also have that the Boltzmann
constant kB=1 such that eps/kB is the unit of temperature.

.SH OUTPUT AND THERMODYNAMICS
\fBmd_blj_n2\fR code produces several output files by default: 
'md_blj_n2.log', 'thermo.dat', and 'parameters_null-md-usedvalues'.
If you instruct other additional files are produced, such as those instructed 
by snapout options. 
The file 'md_blj_n2.log' is the log file were some simulation parameters are 
save such as size of the simulation box, number of particle created if a 
initial condition file were
not used, and so on. In the file 'parameters_null-md-usedvalues' the parameters
values used in the simulation are saved.

Thermodynamics parameters such as pressure and its fluctuation are saved 
in file 'thermo.dat' at time steps according to the value of \fBstepAvg\fR
option.

.SH SPATIAL DIMENSIONS
\fBmd_blj_n2\fR code may be run in two or three spatial dimensions. 
To choose two dimensions 
edit the file
"vectdefs.h" in directory "NagBody_pkg/General_libs" and choose TWODIM. 
Recompile again the code.

.SH STOPPING THE CODE
Once the \fBmd_blj_n2\fR is running you may always stop it by executing the command

echo > stop

You must be in the same directory were the process were lunched.
A file named 'stop-state' is 
saved containing the state of the simulation run. This option permits us to
stop the simulation and change some numerical and/or physical parameters
such as temperature or
density.

.SH EXAMPLES
Executing

md_blj_n2

will run the default simulation which is consistent with the experimental data reported
in: J.L. Yarnell, et al, "Structure factor and radial distribution function for liquid
Argon at 85 K", Phys. Rev. A7 (1972) 2130. Parameter eps in the Lennard-Jones potential
that fit the experimental data is eps/kB = 119 K.


With the command:
md_blj_n2 nbody=4096 out=snap%03d dtout=2/256

we run a simulation with 4096 bodies and save snap data every other time step.

With the command:
md_blj_n2 statefile=state

we run a simulation were a state file is saved every dtout/dtime steps. If for some
reasons the simulation run is interrupted you may always restart the run ejecuting
the command:

md_blj_n2 restorefile=state

Or if you stop it the run using the command 'echo > stop' you may restart the run with
the command:

md_blj_n2 restorefile=stop-state
.sp
.sp
Ejecuting the command:
.sp
md_blj_n2 icModel=6 nbodyprop=2/2 tstop=100 d=0.8 o=snap stepvAvg=30 > output &
.sp
a run is launched in the background using an FCC initial condition whose basic
cell has 2/2 bodies proportion with a density=0.8, snaps will be saved every 5
time steps and properties will be computed every 30 steps.
.sp
Now, we may stop the simulation using 'echo > stop', then we re-launch the
simulation using the command:
.sp
md_blj_n2 restore=stop-state stepAvg=70 > output01 &
Therefore, properties are now computed every 70 time steps and the output of 
the run is save in output01.

.SH ANIMATIONS
You may run the \fBanalysis_md\fR code to see animation plots.
Run a simulation using, for example, the command:

md_blj_n2 icModel=6 nbodyprop=2/2 tstop=100 d=0.8 o=snap > output &
.br
then, use the command:

analysis_md in=snap analysis_type=snap-anim

to see the animation of the simulation. Instead of <snap-anim> you may
use <snap-animtrajectory> to see the trajectories of the particles.
If you want to distiguish the type of particles use also \fBoptions\fR=<types>.

The x and y ranges may change according to 
particle positions during simulation, then, can be useful to set
\fBxrange\fR=<xmin:xmax> and \fBxrange\fR=<ymin:ymax>.

Then you may execute the command:

analysis_md in=snap analysis_type=snap-animtrajectory options=type withsymbols=false
wd=true xr=-4.5:4.5 yr=-4.5:4.5 fsnap=20

to see trajectory animation distinguishing particles types with a x and y ranges
fix and showning up to 20 snaps in total.

.SH SEE ALSO
md_blj(1), analysis_md(1), nplot2d(1).

.SH COPYRIGHT
Copyright (C) 1999-2011
.br
M.A. Rodriguez-Meza
.br
