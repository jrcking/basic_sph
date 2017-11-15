# basic_sph
A basic SPH code - weakly compressible formulation

It is about as basic possible whilst still providing a reasonable looking result for the dam break test case.

This uses the Weakly Compressible formulation. Leap-frog time integration. Cubic spline kernel function. Tait EoS for pressure, with artificially low sound speed.

Boundary conditions are applied using a basic mirror particle technique. Boundaries are slip-walls.

Density is re-initialised periodically, using a first-order corrected kernel.

THERE IS NO TENSILE CORRECTION.
THERE IS NO KERNEL GRADIENT CORRECTION.

The file simulation_parameters does as it says.
    * rlow and rhigh set the coordinates of the rectangular domain.
    * pblow and pbhigh set the coordinates of a rectangular box of particles for initialisation.


