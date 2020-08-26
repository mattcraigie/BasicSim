Welcome to BasicSim, the premier software for basic cosmology N-body simulations!

-------------
Documentation
-------------

Documentation and help with getting started is available at https://basicsim.readthedocs.io/


------------
Installation
------------

To install the requirements, activate your working environment and use a pip install:

    >> pip install -r requirements.txt

FFMpeg must also be installed (https://ffmpeg.org/). Unless you are using the default
FFMpeg install location, you will need to add the lines

    >> import matplotlib as mpl
    >> mpl.rcParams['animation.ffmpeg_path'] = r'path\\to\\FFMpeg\\bin\\ffmpeg.exe'

before calling the simulation functions.

