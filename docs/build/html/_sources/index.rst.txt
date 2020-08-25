======================
BasicSim Documentation
======================

Introduction
------------

BasicSim is a python package that provides a simple and easy to use N-body simulation for cosmological-scale interactions. It uses a particle-mesh
technique which allows efficient calculation for large numbers of particles. BasicSim provides one main function that
calculates gravitational interaction by solving the Poisson equation for gravitational potential and then calculating
the gradient to provide acceleration. With a preference for simplicity (and allowing me to actually finish it), it
uses a simple Euler integration method. The Fast Fourier Transform (FFT) technique that underlies the particle-mesh method scales with
O(N + G log[G]) where N is the number of particles and G is the number of grid points, in contrast with O(N^2) of the
direct point-point comparison method.

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

If you wish to contribute to BasicSim please raise an issue on the BasicSim `Github <https://github.com/mattcraigie/BasicSim>`_.
