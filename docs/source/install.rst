============
Installation
============

The most straightforward way to use BasicSim is:

1. Clone the github repository at https://github.com/mattcraigie/BasicSim using your favourite IDE or using::

    git clone https://github.com/mattcraigie/BasicSim

2. Activate your preferred python environment.
3. Navigate to the BasicSim toplevel directory and use::

    pip install -r requirements.txt

4. Install FFMpeg on your machine (https://ffmpeg.org/)
5. You're good to go! Create a script or run python interactively, and import the simulation module using

.. code-block:: python

   import simulation
   simulation.pm.run_sim_pm(filename)


You may also need to specify the FFMpeg install location by adding the lines

.. code-block:: python

   import matplotlib as mpl
   mpl.rcParams['animation.ffmpeg_path'] = r'path\\to\\FFMpeg\\bin\\ffmpeg.exe'

before running any BasicSim functions.


