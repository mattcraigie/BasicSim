======================
BasicSim Documentation
======================

Introduction
------------

BasicSim is a python package that provides a simple and easy to use 3D N-body simulation for cosmological-scale interactions. It uses a particle-mesh
technique which allows efficient calculation for large numbers of particles, allowing users to watch as cosmological strcture grows.
BasicSim provides one main function that
calculates gravitational interaction by solving the Poisson equation for gravitational potential and then calculating
the gradient to provide acceleration. With a preference for simplicity (and allowing me to actually finish it), it
uses a simple Euler integration method. The Fast Fourier Transform (FFT) technique that underlies the particle-mesh method scales with
O(N + G log[G]) where N is the number of particles and G is the number of grid points, in contrast with O(N^2) of the
direct point-point comparison method.

The particle-mesh technique used in BasicSim is developed from Andrey Kravtsov's incredibly useful presentation on producing a
particle-mesh cosmology simulation (available at https://astro.uchicago.edu/~andrey/talks/PM/pmpotsdam.pdf). It uses
a system of normalised comoving coordinates for position and velocity, and steps in the scale factor :math:`a` as the time step. It produces
an x-y projection of the 3D space in normalised comoving coordinates :math:`x` that may be recovered to a separation :math:`r`
by :math:`r=a x /r_0`. To enable BasicSim to run smoothly, it always begins at :math:`a = 0.4` and steps forward in steps of
:math:`0.001` per frame. The initial conditions are always a random postion distribution at this time. This is somewhat
unphysical, but an approximation that must be made for BasicSim. One artifact of the method used (and a major limitation of this code)
is the treatment of particles at :math:`a > 1`. Particles tend to 'freeze out' as the separation becomes too far for
the gravitational field to overcome. As a result, longer simulation times are not recommended (i.e. past about 1000 frames),
and if you have time to run a large simulation it may be worth favouring higher N for more interesting results.

Additionally, BasicSim provides the means to calculate and plot the power spectrum and correlation function at each step
throughout the simulation. To maintain the calculation speed, the resulting power spectrum and correlation functions are
quick and dirty again using FFTs. Changing to more precise pairwise methods requires O(N^2) making calculations far
less practical for the large number of particles intended to be used with BasicSim.

Check out the different sections of this documentation for help:

.. toctree::
   :maxdepth: 2

   install
   simfuncs
   documentation
   simulation
   modules


If you wish to contribute to BasicSim please raise an issue on the BasicSim `Github <https://github.com/mattcraigie/BasicSim>`_.
I would like to acknowledge the assistance of Dragan Huterer's notes on constructing density, as well as a significant contribution
from Cullan Howlett in understanding intricacies of particle-mesh and the resulting power spectrum.