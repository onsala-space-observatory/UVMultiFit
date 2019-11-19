# SET OF TEST EXAMPLES OF UVMULTIFIT.

## HOW TO RUN:

Execute casa from this directory and then:

    execfile('TEST_ALL.py')

You may first need to configure this script, by setting these parameters:


    DoSimObs = True   # Re-generate all the simulated datasets
                      # (only needed for the first time you run the tests)

    DoFit = True      # Re-do all the example fits.

    casaexe='casa --nologger'  # execution command for CASA
                               # (may depend on your installation)

The script will create several diagnostic files in two directories:

- FIT_SCRIPTS: all the scripts used in the fits (with the model
               definitions, initial parameters, bounds, etc.). These
               scripts are automatically generated by the `MASTER*.py`
               scripts of the test suite (see below).

- ALLTESTS.RESULTS: plots of dirty images, post-fit residual images,
                    parameter plots (for spectral-line fits), and the
                    file `FIT_RESULTS.DAT`

- In `FIT_RESULTS.DAT`, the user will find the fitted parameter
  values, together with the true values used in the simulations, and
  the execution times of all fittings.


If you want to run one particular test, go to its directory
(`TEST*_*`) and execute the `MASTER*.py` script from `casa`. You will
need to manually define the `DoSimObs`, `DoFit` and `casaexe`
variables before running it.

The `MASTER*.py` scripts generate the simulated data and cfreate the
fitting scripts, which are then executed in independent casa sessions,
started by MASTER*.py through the `os.system` approach.

The fitting scripts created by `MASTER*.py` are based on the `test*.py`
scripts in the same directory. if you want to change any detail in the
simulation and/or fitting strategy used in the tests, you may need to
change both, the `MASTER*.py` and the `test*.py`.


## Description of the tests:

### TEST 1 - Fitting a source ellipticity.

The source appears unresolved in the image plane. The axis ratio and
position angle are, within uncertainties, compatible with the ones
used in the simulation.



TEST 2: - Fitting sources with different radial intensity profiles.

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



TEST 3: - Fitting multiple sources simultaneously.

The source is a disc with a hole (the whole is modelled using an inner
disc with a negative flux density) and a gaussian with an offset
position. Notice that for the hole model to work, the surface
brightness of the two discs must be the same (in absolute value). All
sources appear unresolved in the image plane.



TEST 4: - Widefield mosaic.

Simulations of an l-band vla mosaic with three pointings, covering an
area of about 1 degree squared. the sources are distributed in a grid
accross the mosaic area. This dataset is useful to check the
primary-beam correction and the mosaic re-projection algorithm.



TEST 5: - A Gaussian source observed with PdB.

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



TEST 6: - IMMULTIFIT.

A simple model (two Gaussians) fitted to the gridded visibilities in
the image plane.


TEST 7: - A test of the "only_flux" mode,

which can be used to fit many (tens or hundreds!) source components
whose shapes and coordinates are fixed in the model (i.e., only the
flux densities are fitted). This mode can be used as an interface to
implement modelling based on compressed-sensing techniques.

Notice that simobserve introduces a small sub-pixel shift in RA, which
introduces small systematic residuals when the sources are fixed to
their nominal positions (look at the plot of residuals).


TEST 8: - An example of global fringe fitting, a-la UVMULTIFIT.

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