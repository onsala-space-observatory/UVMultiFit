# Installation

Steps to install the `UVMULTIFIT` package:

## Install dependencies

### For Linux (Ubuntu-like), in a terminal:

    sudo apt-get install libgsl-dev libfftw3-dev

### For Mac OS, in a terminal:
_This may need an update!_

    sudo port install gsl fft3
    export LIBRARY_PATH="/opt/local/lib"
    export LD_LIBRARY_PATH="/opt/local/lib"

(and either run the setup script from that terminal or update your
`LIBRARY_PATH`s in your configuration).

## Clone the repository or update your installation

Pick a directory where you want to install `UVMultiFit`, e.g.

	cd ~
    mkdir -p .casa/Nordic_tools
	cd .casa/Nordic_tools
    git clone https://github.com/onsala-space-observatory/UVMultiFit.git
	cd UVMultiFit

If you already had a `git` based version installed, pull in the latest changes:

	git pull

## Compile the C++ module

    python setup.py build_ext --inplace

After this step, the file `_uvmultimodel.so` should have been
created. You may check that you can load the module and execute some
of its support functions:

    python test/test_uvmultimodel.py

which should produce output like this

> ....
> ----------------------------------------------------------------------
> Ran 4 tests in 0.000s
>
> OK

## Running the test suite

From the `UVMultiFit` directory created above, do

    cd test
	python TEST_ALL.py

## Use inside CASA

* Import the module into CASA. You can do it from the init.py file or
   from inside CASA. To import the module from the init file, add the
   following line at the end of your $HOME/.casa/init.py:

        UVMULTIFIT_PATH = "$HOME/.casa/Nordic_Tools/UVMultiFit"
        import imp
        uvm = imp.load_source('uvmultifit', UVMULTIFIT_PATH+'/uvmultifit.py')

   To import the module from CASA, just execute the same line in a script,
   or at the CASA prompt.

 * ENJOY!

Any feedback and bug report should be sent either to the ARC Nordic
Node (contact@nordic-alma.se) or to the source maintainer
(michael.olberg@chalmers.se).
