# The binary file structure of *P_SpaceChargeLight*

The binary data files, written by P_SpaceChargeLight contain integer and floating point values, as (currently) defined in the following:

integer: 4 byte, little endian

single: 4 byte, little endian


Note that all calculations are internally performed with double-precision (8 byte variables), however to achive compatibility to gnuplot's binary file reader, the binary files are stored in single precision, only.

Binary files always represent sectional planes through the full 3D grid. Thus, the final solution will be stored as a set of NodeCount_x-1 data files that represent y-z planes. For the interim solutions, the user can specify a sectional plane of interest, that is stored as binary file every k-th iteration step.

**Filestructure**


**num_col** Number of columns (integer)

**y_1,y_2,...,y_NodeCount_y-1** 	Coordinates of mesh points in y direction (single)

**z_1,D(x,y_1,z_1),D(x,y_2,z_1),...,D(x,y_NodeCount_y-1,z_1)**	Coordinate of first mesh point in z direction, followed by values of the specific 3D data array (single)

**z_2,D(x,y_1,z_2),D(x,y_2,z_2),...,D(x,y_NodeCount_y-1,z_2)**	Coordinate of second point in z direction, followed by values of the specific 3D data array (single)

**...**

**First example of a Gnuplot file**

The following example demonstrates, how to generate a 2D color/contour plot from the binary files with the help of Gnuplot.


> reset
> 
> set term png
> 
> set title "Electrostatic potential of a probe tip facing a semiconductor surface"
> 
> set xrange [ -20: 20] noreverse nowriteback
> 
> set yrange [ -20: 20] noreverse nowriteback
> 
> set cblabel "Phi [V]"
> 
> set pm3d interpolate 2,2
> 
> set pm3d map
> 
> set size square
> 
> set contour
> 
> set cntrparam levels 30
> 
> aFile = "phi_V=2.00V_x=32.bin"
> 
> aOutput = "phi.png"
> 
> set output aOutput
> 
> splot aFile binary matrix with pm3d
> 
> exit gnuplot;



**Second example of a Gnuplot file**

In this example, the interim solutions (every 100th k-step) are used to generate a movie that indicates the evolution of the solution. Such movies can be generated for the electrostatic potential and the carrier concentrations. This evolution of the solution must not be confused with a time evolution! P_SpaceChargeLight finds a steady-state solution and not a time-depended one. However, the visualization of the evolution of the solution may be helpful for error determination purposes.

> reset
> 
> set term gif animate
> 
> set view map
> 
> set xrange [ -20: 20] noreverse nowriteback
> 
> set yrange [ -20: 20] noreverse nowriteback
> 
> set cblabel "Phi [V]"
> 
> #set cbrange [0.0 : 0.7] noreverse nowriteback
> 
> set pm3d interpolate 2,2
> 
> set pm3d map
> 
> set hidden3d
> 
> set size square
> 
> j = 1000
> 
> aOutput = "phi_interim_solution.gif"
> 
> set output aOutput
> 
> do for [i=0:j] {
> 
> aFile = sprintf("phi_k=%d.bin",i*100)
> 
> set title aFile
> 
> splot aFile binary matrix with pm3d
> 
> }
> 
> exit gnuplot;

