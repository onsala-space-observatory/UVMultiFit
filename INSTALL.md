# Installation

Steps to install the `UVMultiFit` package:

## Install dependencies

### CASA

These instructions assume you have a `casa` version 5.x installed,
either the normal release or with the ALMA and/or VLA pipeline
included. The instructions below use an installation of
`casa-pipeline-release-5.6.1-8.el7` (ALMA pipeline) as an
example. Earlier versions (4.x) are neither tested nor supported.

**Find out, where `casa` is installed** on your computer as you will need
this information during the installation of `UVMultiFit`.

### For Linux (Ubuntu-like), in a terminal:

    sudo apt-get install libgsl-dev libfftw3-dev

### For Mac OS, in a terminal:

If using Mac Ports:

    sudo port install gsl fftw-3

If using Homebrew:

    brew install gsl fftw

## Clone the repository or update your installation

Pick a directory where you want to install `UVMultiFit`, e.g.

    cd ~
    mkdir -p casa/NordicARC
    cd casa/NordicARC
    git clone https://github.com/onsala-space-observatory/UVMultiFit.git
    cd UVMultiFit

If you already had a `git` based version installed, pull in the latest changes:

    git pull

## Compile the C++ module

Open the `Makefile` and edit the first line, where it says

    CASADIR=/home/olberg/Python/casa-pipeline-release-5.6.1-8.el7

and replace with the path of your own `casa` installation. If you have
write access to the directories underneath `$CASADIR`, you may then
run

    make install

This should make `uvmultifit` available also to other users, running
from the same casa installation.

Alternatively, you can install to your local user `PYTHONPATH` by running

    make user

You may check that you can load the module and execute some of its
support functions. Move into the `test` folder and start up `casa`


    $ cd test
	$ <CASADIR>/bin/casa   # replace `CASADIR` with the correct path from above!

and run the short unit test in that folder:

    CASA <1>: %run test_uvmultimodel.py

which should produce output like this

>
> C++ shared library loaded successfully
>
> ....
> ----------------------------------------------------------------------
> Ran 4 tests in 0.000s
>
> OK

## Running the test suite

Still in the `test` directory, and with `casa` running, start the test suite via

    CASA <2>: execfile("TEST_ALL.py")

## Running your own model:

Now, this can be done from any directory.

* Start up `casa` and run

        CASA <1>: from NordicARC import uvmultifit as uvm
        CASA <2>: help(uvm.uvmultifit)                 # to get the help text
        CASA <3>: myfit = uvm.uvmultifit(vis=..., ...) # fit your model ...

 * ENJOY!

Any feedback and bug report should be sent either to the ARC Nordic
Node (contact@nordic-alma.se) or to the source maintainer
(michael.olberg@chalmers.se).
