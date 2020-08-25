

BasicSim
====================================

#.. toctree::
#   :maxdepth: 2


Installation
------------
BasicSim is most easily used by cloning the github repository at https://github.com/mattcraigie/BasicSim.

BasicSim has the following dependencies::

    numpy
    matplotlib
    ffmpeg

To install these quickly, activate your python environment, navigate to the BasicSim directory and use::

    pip install -r requirements.txt

FFMpeg must also be installed (https://ffmpeg.org/). You may need to specify the FFMpeg install location, you will need
to add the lines.

To install these quickly, activate your python environment, navigate to the BasicSim directory and use::

    import matplotlib as mpl
    mpl.rcParams['animation.ffmpeg_path'] = r'path\\to\\FFMpeg\\bin\\ffmpeg.exe'

before running any BasicSim functions.


Contributing
------------

If you wish to contribute to BasicSim please raise an issue
via `Github <https://github.com/mattcraigie/BasicSim>`_.
