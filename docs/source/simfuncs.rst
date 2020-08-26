====================
Running a Simulation
====================

This is a quick guide to running your first custom simulation. After completing the installation and correcting any FFMpeg
issues, try running with the default settings

.. code-block:: python

   from simulation.pm import SimPM
   SimPM('default.mp4')

This will run a simulation with 100 particles on a 30^3 grid for 1000 frames, and save the video to the file 'default.mp4'.


Next, try running the same code with the correlation function and power spectrum plotting enabled,

.. code-block:: python

   SimPM('default_ps.mp4', show_ps=True)

Now, all that's left is to experiment! Remember that the time complexity scales with number of particles N and number of
grid points G as O(N + G log[G]), so pushing N to high values of N ~ 4000 can still be completed in reasonable times.
For example, you might try

.. code-block:: python

   SimPM('highn.mp4', n_particles=4000)
   SimPM('longtime.mp4', n_particles=100, n_frames=3000)
   SimPM('finergrid.mp4', n_particles=1000, n_grid=50, n_frames=500)

Using BasicSim with these higher number of particles results in structure formation in the simulation, which is
detectable in the power spectrum and correlation functions.

Another parameter that can be changed is the gravity factor. This scales the strength of gravity, and modifying it will
change structure formation.

.. code-block:: python

   # For gravity twice as strong
   SimPM('stronggrav.mp4', n_particles=4000, gravity_factor=2)

.. toctree::
   :hidden:

   contents
