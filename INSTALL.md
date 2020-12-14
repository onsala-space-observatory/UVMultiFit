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

Still in the `test` directory, and with `casa` running, run any of the tests via

    CASA <2>: execfile("<test>.py")

where `<test>` is one of:

   * `PROFILES_ELLIPTICITY.py`
   * `MULTICOMP.py`
   * `WIDEFIELD_MOSAIC.py`
   * `GAINS_EXPERIMENTAL.py`
   * `IMMULTIFIT.py`
   * `ONLYFLUX.py`
   * `FRINGE_FITTING.py`

### TEST 1+2: Fitting source ellipticity with different radial intensity profiles

The source appears unresolved in the image plane. The axis ratio and
position angle are, within uncertainties, compatible with the ones
used in the simulation.

All sources appear unresolved in the image plane.

The residuals are all very small, with the exception of the "power-2"
source:

Take into account that the "power-2" source has an *infinite* flux
density, since the source extends to infinity and its 2d integral
does not converge. Hence, the shortest baselines (i.e. the ones
corresponding to the spatial scales of the order of the image size
used in the simulation) have to be removed before the fitting (they
are biased toward lower flux densities).

Furthermore, this "window effect" from the image used by simobserve
introduces artifacts in the residual image. But such bad residuals are
*not* indicative of a wrong fit (these are due to the fact that we
have not used an exact "power-2" model for the simulation, since the
model image needs to be a truncated representation of the source).

A similar thing happens with the "expo" model, which also presents
some residual artifacts (notice the hint of a "missing-flux" effect in
the image).

### TEST 3: Fitting multiple sources simultaneously.

The source is a disc with a hole (the whole is modelled using an inner
disc with a negative flux density) and a gaussian with an offset
position. Notice that for the hole model to work, the surface
brightness of the two discs must be the same (in absolute value). All
sources appear unresolved in the image plane.

### TEST 4: Widefield mosaic.

Simulations of an l-band vla mosaic with three pointings, covering an
area of about 1 degree squared. the sources are distributed in a grid
accross the mosaic area. This dataset is useful to check the
primary-beam correction and the mosaic re-projection algorithm.

### TEST 5: A Gaussian source observed with PdB.

A time-varying gain amplitude is added to one of the antennas. We use
the gain solver feature in UVMULTIFIT, to solve for the source
parameters and the antenna gain simultaneously, Approximating the gain
with a "piece-wise" time function.

This gain-solver feature is still quite experimental. In realistic
situations, a lot of fitting parameters may be needed, to solve for
all the antenna gains with enough time resolution. This may slow down
the code significantly (and there may be risks of falling in a local
minimum of the chi2), although this method should be superior to the
conventional scheme, where the gaincal has no memory of the gain
values close by in time (and the calibration is solved independently
of the source model parametrization).

### TEST 6: IMMULTIFIT.

A simple model (two Gaussians) fitted to the gridded visibilities in
the image plane.

### TEST 7: A test of the "only_flux" mode,

which can be used to fit many (tens or hundreds!) source components
whose shapes and coordinates are fixed in the model (i.e., only the
flux densities are fitted). This mode can be used as an interface to
implement modelling based on compressed-sensing techniques.

Notice that simobserve introduces a small sub-pixel shift in RA, which
introduces small systematic residuals when the sources are fixed to
their nominal positions (look at the plot of residuals).

### TEST 8: An example of global fringe fitting, a-la UVMULTIFIT.

A synthetic C-band dataset is generated with realistic noise and 128
channels (covering 32mhz). The gains are estimated first by fitting
them to the peaks in delay-rate space (using the new helper rutine
"QuinnFF"). These estimates are then used as a-priori for the ordinary
global fringe fitting.

Notice that calibration tables are not generated (that feature is
beyond the purpose of UVMULTIFIT), but the calibrated data can still
be generated, either by dividing the data by the (gain-corrected)
model (which UVMULTIFIT writes in the model column when you set
write='model') or by setting write='calibrated'. the master script in
the test directory shows how to do this.

A fixed point source is used for the gff (centered at the nominal
source position), but of course more complicated structures (and even
*.model images generated with clean/tclean) can be used as models and
you can even fit the source parameters together with the gains!

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
