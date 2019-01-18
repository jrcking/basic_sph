# basic_sph
A basic SPH code - weakly compressible formulation

It is about as basic possible whilst still being capable of simulating a dam break. It is not designed to give nice results, or to run fast, or to be adaptable.

Features:
    * Weakly Compressible formulation, with Tait EoS for pressure, and low sound speed;
    * Leap-frog time integration - as in Price 2012 JCP paper;
    * Cubic spline kernel function;
    * Slip wall BCs applied using mirror particles (but without solid boundary particles - lazy);
    * Density is periodically re-initialised by Shepard filtering;
    * Viscosity is included by artificial viscosity .

What is doesn't contain:
    * Tensile correction;
    * Kernel gradient correction;
    * Shifting.

The file simulation_parameters does as it says.
    * rlow and rhigh set the coordinates of the rectangular domain.
    * pblow and pbhigh set the coordinates of a rectangular box of particles for initialisation.

To compile, type "make" or try "gfortran ./src/*.f90 -o ksph"  <----- you might need to do this a few times to build the mod files.

To run: "./ksph". It will set up conditions for a dam break, and run until "tmax" (as set in parameter file), outputting every "dmp_period".

Output files will be stored in the "output" directory, which also contains a Matlab/Octave script to plot particle positions (ppos.m) and a script to display an animation (anim.m).

