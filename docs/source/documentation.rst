==================
Full Documentation
==================

.. py:class:: simulation.pm.SimPM(filename, n_particles=100, n_grid=30, n_frames=1000, show_ps=False, ps_bins=50, gravity_factor=1)

   Creates and runs a simulation using a particle-mesh N-body implementation. An mp4 file of the simulation is saved at
   the location specified by filename. The particle-mesh method scales with O(N + G*log[G]), and operations are done in
   place so large number of particles can be computed in reasonable times. Setting show_ps=True will slow computation
   but is also calculated with O(G*log[G]).

   :param filename: name of the file to save the simulation video
   :param n_particles: the number of particles to use in the simulation.
   :param n_grid: number of grid cells along a single dimension, total grid cells will be n_grid**3
   :param n_frames: number of frames to run the simulation for
   :param show_ps: set to true to plot the power spectrum (k * P[k]) and correlation function (xi[r]) alongside the simulation
   :param ps_bins: the number of bins to use for the power spectrum and correlation function if show_ps=True
   :param gravity_factor: a scale factor used to increase or decrease the strength of gravity, for experimenting in simulations
   :return: returns None


