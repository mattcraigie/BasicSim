====================
Running a Simulation
====================

This is a quick guide to running your first custom simulation. After completing the installation and correcting any FFMpeg
issues, try running with the default settings

.. code-block:: python

   import simulation as sim
   sim.run_sim_pm('default.mp4')

This will run a simulation with 100 particles on a 30^3 grid for 1000 frames, and save the video to the file 'default.mp4'.


Next, try running the same code with the correlation function and power spectrum plotting enabled,

.. code-block:: python

   sim.run_sim_pm('default_ps.mp4', show_ps=True)

Now, all that's left is to experiment! Remember that the time complexity scales with number of particles N and number of
grid points G as O(N + G log[G]), so pushing N to high values of N ~ 4000 can still be completed in reasonable times.
For example, you might try

.. code-block:: python

   sim.run_sim_pm('highn.mp4', n_particles=16**3)
   sim.run_sim_pm('longtime.mp4', n_particles=100, n_frames=3000)
   sim.run_sim_pm('finergrid.mp4', n_particles=1000, n_grid=50, n_frames=500)