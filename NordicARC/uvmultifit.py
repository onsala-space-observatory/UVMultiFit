# NordicARC/uvmultifit
"""A wrapper to perform a fit visibility data.

The function `uvmultifit` will create an instance of a `Modeler` class
as well as a `MeasurementSet` class and perform a fit.

Fit results will be written to an output file using method `save_results`.

"""
import sys
import numpy as np
import logging

from typing import Union

from .measurementset import MeasurementSet
from .modeler import Modeler

def save_results(outfile: str, results: dict, mdl: Modeler, ms: MeasurementSet) -> None:
    """Save fit results to file.

    Given the results of the fit, and the instances of the `Modeler` and
    `MeasurementSet` class, write the results to an output file.

    The method does not return anything.

    Args:
        outfile: name of the outputfile.
        results: results returned by `mdl.fit`
        mdl: instance of the `Modeler` class, which performed the fit.
        ms: instance of a `MeasurementSet` class, on which the fit was performed.
    """

    Nvis = results["Degrees of Freedom"]
    if isinstance(Nvis, list):
        DOG = int(np.average(Nvis))
    else:
        DOG = int(Nvis)
    ChiSq = results["Reduced Chi squared"]
    fitparams = results["Parameters"]
    fiterrors = results["Uncertainties"]
    if not mdl.OneFitPerChannel:
        prtpars = []
        for pp, ppar in enumerate(fitparams):
            prtpars.append(ppar)
            prtpars.append(fiterrors[pp])
    else:
        nspwtot = ms.spwlist[-1][3] + len(ms.spwlist[-1][2])
        prtpars = [[[] for spi in range(len(fitparams[sp]))] for sp in range(nspwtot)]
        for sp in range(nspwtot):
            for pp, ppar in enumerate(fitparams[sp]):
                for ppi, ppari in enumerate(ppar):
                    prtpars[sp][pp].append(ppari)
                    prtpars[sp][pp].append(fiterrors[sp][pp][ppi])

    with open(outfile, 'w') as outf:
        outf.write("# MODELFITTING RESULTS FOR MS: " + ','.join(ms.vis) + "\n")
        outf.write("# NOTICE THAT THE REDUCED CHI SQUARE SHOWN HERE IS THE VALUE\n")
        outf.write("# *BEFORE* THE RE-SCALING OF THE VISIBILITY WEIGHTS.\n")
        outf.write("# HOWEVER, THE PARAMETER UNCERTAINTIES ARE *ALREADY* GIVEN\n")
        outf.write("# FOR A REDUCED CHI SQUARE OF 1.\n")
        outf.write(f"# AVG. NUMBER OF DEGREES OF FREEDOM: {DOG}\n")
        if ms.pbeam:
            if isinstance(ms.dish_diameter, float):
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. "
                           "USING A DISH DIAMETER OF: {ms.dish_diameter:.3f} METERS\n")
            else:
                outf.write("# PRIMARY-BEAM CORRECTION HAS BEEN APPLIED. USING: \n")
                for i, name in enumerate(ms.antnames):
                    outf.write(f"#    FOR ANTENNA {name} A DIAMETER OF {ms.userDiameters[i]:.2f}m \n")

        else:
            outf.write("# PRIMARY-BEAM CORRECTION HAS NOT BEEN APPLIED.\n")

        outf.write("###########################################\n")
        outf.write("#\n# MODEL CONSISTS OF:\n")
        for m, mod in enumerate(mdl.model):
            outf.write("# '" + mod + "' with variables: " + mdl.var[m] + "\n")
        if len(mdl.fixed) > 0:
            outf.write("#\n#\n# FIXED MODEL CONSISTS OF:\n")
            if mdl.takeModel:
                outf.write("# THE MODEL COLUMN FOUND IN THE DATA\n")
            else:
                for m, mod in enumerate(mdl.fixed):
                    var2print = list(map(float, mdl.fixedvar[m].split(',')))
                    outf.write(
                        "# '" + mod + "' with variables: " + ' '.join(['%.3e'] * len(var2print)) % tuple(var2print))
                outf.write(f"#\n#  - AND SCALING FACTOR: {mdl.scalefix}\n")

        outf.write("#\n# INITIAL PARAMETER VALUES:\n")
        for p0i, p0 in enumerate(mdl.p_ini):
            if mdl.bounds is not None:
                outf.write(f"#  p[{p0i}] = {p0:.5e} with bounds: {str(mdl.bounds[p0i])}\n")
            else:
                outf.write(f"#  p[{p0i}] = {p0:.5e} with no bounds\n")
        outf.write("#\n##########################################\n")

        N = len(mdl.p_ini)
        parshead = [n for n in range(N) for i in range(2)]

        headstr = ("# Frequency (Hz)     " + "  p[%i]       err[%i]   " * N + "   red. ChiSq\n") \
            % tuple(parshead)
        outf.write(headstr)

        formatting = "%.12e   " + "%.4e " * (2 * N) + "   %.4e \n"
        if not mdl.OneFitPerChannel:
            toprint = tuple([np.average(ms.averfreqs)] + prtpars + [ChiSq])
            outf.write(formatting % toprint)

        else:
            for spwvi in ms.spwlist:
                for r, rr in enumerate(spwvi[2]):
                    k = spwvi[3] + r
                    # rang = rr[1]
                    for nu, freq in enumerate(ms.averfreqs[k]):
                        toprint = tuple([freq] + prtpars[k][nu] + [ChiSq[k][nu]])
                        outf.write(formatting % toprint)

def uvmultifit(vis: str, spw: int = 0, field: Union[int, str] = 0, scans: list = [],
               column: str = 'data', uniform=False, uvtaper=0.0, chanwidth=1, timewidth=1, stokes='I',
               MJDrange=[-1.0, -1.0], ldfac=1.22,
               phase_center='', pbeam=False, wgt_power=1.0, dish_diameter=0.0,
               model=['delta'], var=['p[0],p[1],p[2]'], p_ini=[0.0, 0.0, 1.0], bounds=None,
               OneFitPerChannel=True, only_flux=False, fixed=[], fixedvar=[], scalefix='1.0',
               phase_gains={}, amp_gains={}, method='levenberg', HankelOrder=80,
               LMtune=[1.0e-3, 10.0, 1.0e-5, 200, 1.0e-3], SMPtune=[1.0e-4, 1.0e-1, 200],
               NCPU=4, proper_motion=0.0,
               write='', outfile='modelfit.dat', cov_return=False, finetune=False):
    """Create instances of `Modeler` and `MeasurementSet` and perform fit.

    The first group of arguments (`vis` to `dish_diameter`) are used to instantiate a
    `MeasurementSet` object. The second group (`model` to `proper_motion`) are used to create
    a `Modeler` object.

    The final four arguments are used for performing the actual fit and saving the data to a file.

    Args:
       vis: Name of the directory containing the measurement set.
            Or list of measurement sets.
       spw: Integer defining the index of the spectral window.
            Or select many channel ranges using strings in the usual CASA way.
       field: The id (or name) of the target source in the measurement set.
       scans: The id numbers of the scans to load.
              Default means to load all the scans of the selected source.
              If multiple measurement sets are selected, this should be a list of lists
              (i.e., one list of scans per measurement set).
              For instance, if `vis = ["ms1.ms","ms2.ms"]`, then `scans = [[1, 2, 3],[]]`
              would select scans 1, 2, and 3 from the `ms1.ms` data set
              and all the scans from the `ms2.ms` dataset.
       column: The data columns, it should be either `data` or `corrected`.
    """

    logger = logging.getLogger("uvmultifit")
    logging.basicConfig(level=logging.INFO,
                        format='%(name)s - %(levelname)s - %(message)s')
    logger.info("started")
    mdl = Modeler(model, var, p_ini, bounds, OneFitPerChannel, only_flux, fixed, fixedvar,
                  scalefix, phase_gains, amp_gains, method, HankelOrder,
                  LMtune, SMPtune, NCPU, proper_motion)

    corrected = column == "corrected"
    ms = MeasurementSet(vis, spw, field, scans, corrected, uniform,
                        uvtaper, chanwidth, timewidth, stokes,
                        MJDrange, ldfac, phase_center, pbeam, wgt_power, dish_diameter)
    if not ms.check_measurementset():
        logger.error("failed to check measurement set")
        sys.exit(1)

    if not ms.read_data():
        logger.error("failed to read data")
        sys.exit(1)

    nspw = ms.Nspw
    if not mdl.init_data(nspw, ms):
        logger.error("failed to init data")
        sys.exit(1)

    if not mdl.init_model(nspw, ms.tArr, ms.averfreqs, ms.refpos):
        logger.error("failed to init model")
        sys.exit(1)

    write_model = {'': 0, 'model': 1, 'residuals': 2, 'calibrated': 3}
    write_model_index = write_model.get(write, -1)
    if write_model_index == -1:
        logger.error("keyword 'write' should be set to either '', 'model', 'residuals' or 'calibrated'")

    results = mdl.fit(ms, write_model_index, cov_return)
    if not results:
        logger.error("failed to fit model")
        sys.exit(1)

    if not finetune:
        #       if not self.fit():
        #             logger.error("failed fit!")
        #             return False
        if write_model_index in [1, 2, 3]:
            if timewidth == 1 and stokes not in ['Q', 'U', 'V']:
                logger.info("writing into measurement set(s)")
                # writeModel()
            else:
                msg = "cannot write into measurement set!\n"
                msg += "If you want to fill-in the model (or corrected) column:\n"
                msg += "    1.- 'timewidth' and 'chanwidth' should be set to 1\n"
                msg += "    2.- 'stokes' should be set to either 'I', 'PI', or a corr. product."
                logger.error(msg)

        logger.info("fit done!")

    column = "MODEL_DATA" if write_model_index == 1 else "CORRECTED_DATA"
    ms.write_model(column, mdl)
    save_results(outfile, results, mdl, ms)

    return results
