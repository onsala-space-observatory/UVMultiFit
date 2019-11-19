# Installation

Steps to install the `UVMULTIFIT` package:

## Install dependencies (i.e., GNU Scientific Library, so far).

### For Linux (Ubuntu-like), in a terminal:

    sudo apt-get install libgsl-dev libfftw3-dev

### For Mac OS, in a terminal:

    sudo port install gsl fft3
    export LIBRARY_PATH="/opt/local/lib"
    export LD_LIBRARY_PATH="/opt/local/lib"

(and either run the setup script from that terminal or update your
`LIBRARY_PATH`s in your configuration).

### Compile the C++ module. Just run:

    python setup.py build_ext --inplace

After this step, the file `_uvmultimodel.so` should have been created.

 * Copy the `uvmultifit.py` and `_uvmultimodel.so` files into a
   destination directory. The name of this directory could be, e.g.

        ~/.casa/Nordic_Tools/uvfit

 * Import the module into CASA. You can do it from the init.py file or
   from inside CASA. To import the module from the init file, add the
   following line at the end of your $HOME/.casa/init.py:

        UVMULTIFIT_PATH = "$HOME/.casa/Nordic_Tools/uvfit"
        import imp
        uvm = imp.load_source('uvmultifit', UVMULTIFIT_PATH+'/uvmultifit.py')

   To import the module from CASA, just execute the same line in a script,
   or at the CASA prompt.

 * ENJOY!

Any feedback and bug report should be sent either to the ARC Nordic
Node (contact@nordic-alma.se) or to the source maintainer
(michael.olberg@chalmers.se).
