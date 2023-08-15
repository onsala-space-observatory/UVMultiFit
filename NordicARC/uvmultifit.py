"""A wrapper to perform a fit to visibility data.

The function `uvmultifit` will create an instance of a `Modeler` class
as well as a `MeasurementSet` class and perform a fit.

Fit results will be written to an output file using method `save_results`.

"""
import sys
import logging

from typing import List, Dict

import numpy as np

from .measurementset import MeasurementSet
from .modeler import Modeler

def save_results(outfile: str, results: Dict, mdl: Modeler, ms: MeasurementSet) -> None:
    """Save fit results to file.

    Given the results of the fit, and the instances of the `Modeler` and
    `MeasurementSet` class, write the results to an output file.

    The method does not return anything.

    Args:
        outfile (str): name of the outputfile.
        results (dict): results returned by `mdl.fit`
        mdl (Modeler): instance of the `Modeler` class, which performed the fit.
        ms (MeasurementSet): instance of a `MeasurementSet` class, on which the fit was performed.

    Returns:
        None
    """

    Nvis = results["Degrees of Freedom"]
    if isinstance(Nvis, list):
        DOG = int(np.average(Nvis))
    else:
        DOG = int(Nvis)
    ChiSq = results["Reduced Chi squared"]
    fitparams = results["Parameters"]
    fiterrors = results["Uncertainties"]
    if not mdl.spectral_mode:
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
        outf.write("#\n")

        N = len(mdl.p_ini)
        parshead = [n for n in range(N) for i in range(2)]

        #                    1         2         3         4         5         6
        #          01234567890123456789012345678901234567890123456789012345678901234567890123456789
        headstr = "#          GHz" + "        p[%i]       err[%i]" * N + "   red.ChiSq\n"
        outf.write(headstr % tuple(parshead))

        formatting = "%14.9f" + "%12.4e %12.4e" * N + "%12.4e\n"
        if not mdl.spectral_mode:
            toprint = tuple([np.average(ms.averfreqs)/1.0e9] + prtpars + [ChiSq])
            outf.write(formatting % toprint)

        else:
            for spwvi in ms.spwlist:
                for r, _ in enumerate(spwvi[2]):
                    k = spwvi[3] + r
                    # rang = rr[1]
                    for nu, freq in enumerate(ms.averfreqs[k]):
                        toprint = tuple([freq/1.0e9] + prtpars[k][nu] + [ChiSq[k][nu]])
                        outf.write(formatting % toprint)

def uvmultifit(vis: str, spw: str = '0', field: int = 0, scans: List = [],
               column: str = 'data', uniform: bool = False,
               uvtaper: float = 0.0, chanwidth: int = 1, timewidth: int = 1,
               stokes: str = 'I', MJDrange: List[float] = [-1.0, -1.0], ldfac: float = 1.22,
               phase_center: str = '', pbeam: bool = False, wgt_power: float = 1.0,
               dish_diameter: float = 0.0, model: List[str] = ['delta'],
               var: List[str] = ['p[0],p[1],p[2]'], p_ini: List[float] = [0.0, 0.0, 1.0],
               bounds: List = None, OneFitPerChannel: bool = True, only_flux: bool = False,
               fixed: List = [], fixedvar: List = [], scalefix: str = '1.0',
               amp_gains: Dict = {}, phase_gains: Dict = {}, method: str = 'levenberg',
               HankelOrder: int = 80, LMtune: List[float] = [1.0e-3, 10.0, 1.0e-5, 200, 1.0e-3],
               SMPtune: List[float] = [1.0e-4, 1.0e-1, 200], NCPU: int = 4,
               proper_motion: float = 0.0, write: str = '', outfile: str = 'modelfit.dat',
               cov_return: bool = False, finetune: bool = False) -> Dict:
    """Create instances of `Modeler` and `MeasurementSet` and perform a fit.

    The first group of arguments (`vis` to `dish_diameter`) are used to instantiate a
    `MeasurementSet` object.

    The second group (`model` to `proper_motion`) are used to create a `Modeler` object.

    The final four arguments are used when performing the actual fit and saving the data to a file.

    Args:
       vis (str): Name of the directory containing the measurement set.
            (Or list of measurement sets.)
       spw (int): Integer defining the index of the spectral window.
            It can also be a string (or list of strings). These are the spectral window(s) and
            channel selection(s) to be fitted. For each spw, the user can select many channel
            ranges in the usual CASA way.
            If one string is given, all measurement sets listed in 'vis' are supposed to have
            the same spectral configuration.
            If a list of strings is given, one per measurement set, this restriction doesn't apply.
       field (int or str):
            The id (or name) of the target source in the measurement set.
       scans (list): The id numbers of the scans to load.
            Default means to load all the scans of the selected source.
            If multiple measurement sets are selected, this should be a list of lists
            (i.e., one list of scans per measurement set).
            For instance, if `vis = ["ms1.ms","ms2.ms"]`, then `scans = [[1, 2, 3],[]]`
            would select scans 1, 2, and 3 from the `ms1.ms` data set
            and all the scans from the `ms2.ms` dataset.
       column (str): The data column, it should be either `data` or `corrected`.
       uniform (bool): Default is False. If True, the weights of all data are made equal.
            Notice that a UV taper can still be applied (i.e., by setting the `uvtaper` parameter).
       uvtaper (float): Default is 0.0. If not 0.0, the weights of the visibilities are multiplied by a
            Gaussian in Fourier space, whose HWHM is the value of `uvtaper`, *in meters*.
       chanwidth (int):
            Number of spectral channels to average in each chunk of data **before the fit**. BEWARE of
            the bandwidth-smearing effects if your baselines are long (and/or your source is large).
       timewidth (int)
            Number of time channels (i.e., integration times) to average in each chunk of data
            **before the fit**. The averaging is always cut at the scan boundaries (so a very large
            timewidth means to average down to one visibility per scan). BEWARE of the
            time-smearing effects!
       stokes (str):
            Polarization product. Can be any of `PI, I, Q, U, V`. It also accepts individual
            correlation products: `XX, YY, XY, YX, RR, LL, LR`, or `LR`. Default is `I`.
            If `PI` (which stands for *Polarization Independent*) is given, the program will
            compute `I` whenever possible and use either `XX, YY, RR` or `LL` otherwise.
            This way, the amount of data used in polarization-independent fits is maximized.
       MJDrange (list):
            List of two floats. These are the initial and final Modified Julian Dates of the data
            to be used in the fitting. Notice that all the scans asked by the user are loaded
            a-priori, and the MJDrange condition is applied *afterwards*. This way, variability
            studies can be performed efficiently, by loading all data at once and then setting
            different MJDranges iteratively. Default (i.e., <= 0.0) means not to select data based
            on JD time range.
       ldfac (float):
            Proportionality constant between the ratio 'lambda/Diameter' and the FWHM of the
            primary beam (assumed to be a Gaussian!). I.e.: FWHM = ldfac*lambda/Diameter.
            Normally, `ldfac = 1.22` should be fine, although 1.00 works better with data coming
            from simobserve.
       phase_center (str):
            The sky position where all components are referenced to. If an empty string is given,
            the phase center will be that of the first field id that is being read (i.e., if a
            mosaic is being read, the first pointing will be set as the phase center).
            If the string is not empty, the program expects to find a coordinate in CASA format
            (i.e., `J2000 RA Dec`, where **RA** is in format `00h00m00.0s` and **Dec** is in
            format `00d00m00.0s`).
       pbeam (bool):
            If false, the primary-beam correction is not applied. This is *very important* for
            fitting mosaic data.
       wgt_power (float):
            Default is 1. Power index of the visibility weights in the computation of the Chi
            square. `wgt_power = 1` would be the *statistically justified* value, although other
            values may be tested if the user suspects that the visibility uncertainties are not
            well estimated in his/her dataset.
       dish_diameter (float):
            In case that the antenna diameters cannot be read from the datasets, the user must
            provide the antenna diameters (in meters). This can be given as a single float (so
            the array is assumed to be homogeneous) or as a dictionary, whose keys are
            antenna names (or *regular expressions*, that match antenna-name patterns) and whose
            elements are the antenna diameters (in meters).

            *Note*: For primary beam correction of heterogeneous data sets, please use either
                    one concatenated Measurement Set (MS) for all your data or several MSs
                    with *the same* antenna table.
       model (list of str):
            List of strings (i.e., model components to fit). Each component is given as a string.
            Possible models are: `delta, disc, Gaussian, ring, sphere, bubble, expo, power-2,
            power-3`, and `GaussianRing`.
            If only one model component is being fitted, the `model` keyword can also be a string.
       var (list of str):
            List of strings (or just one string, if only one model component is being fitted).
            These are the variables for each model. The variables can be set to *any* algebraic
            expression involving the fitting parameters (being the ith parameter represented by
            `p[i]`) and the observing frequency in Hz (represented by `nu`). Any numpy function
            can also be called, using the prefix `np` (e.g., `p[0]*np.power(nu, p[1])`).
            See some examples below.
       p_ini (list):
            List of the initial values of the fitting parameters. This is expected to be a list of floats.
       bounds (list):
            List of boundaries (i.e., minimum and maximum allowed values) for the fitting
            parameters. 'None' means that no bound is defined. If the list is empty, no bounds
            are assumed for any parameter (see examples below).
       OneFitPerChannel (bool):
            If True, independent fits are performed to the different frequency channels, one by
            one. If False, one common fit is performed to all data. In this case, the user may
            want to fit for the spectral indices of the components, if the fractional bandwidth
            is wide.
       only_flux (bool):
       fixed (list):
            Like **model**, but defines model components with completely fixed variables (i.e.,
            whose variables are defined only by numbers; not fitting parameters). This model will
            be computed only once (i.e., just before the fit), hence making the code execution
            faster. The user can load the model column of the measurement set(s) as a fixed model,
            by setting `fixed='model_column'`.
       fixedvar (list):
            Like **var**, but refers to the **fixed** model. Hence, it is expected to be either a
            list of numbers or a list of strings representing numbers. This is not needed if
            `fixed = 'model_column'` (since, in that case, the model column is read *as is* from
            the measurement set).
       scalefix (str):
            String representing a function that sets a scaling factor for the fixed-model's total
            flux density. It *can be* a function of the fitting parameters (e.g., `scalefix='p[0]'`
            will multiply the overall flux density of the fixed model by `p[0]`) and can also be
            a function of the observing frequency `nu`.
       amp_gains (dict):
            Dictionary (default empty). Controls whether to solve for antenna amplitude gains
            and source-model parameters *simultaneously*. The keys of the dictionary are the
            indices of the antennas to self-calibrate. The values of the dictionary are strings,
            whose elements represent functions of the fitting parameters. See examples below.
       phase_gains (dict):
            Same as for amp_gains, but dealing with phase gains.

            *Note*: For the gain solver to work, please use either one concatenated Measurement
                      Set (MS) for all your data or ensure that all your MSs have the same
                      antenna table.
       method (str):
       HankelOrder (int):
            Only used for models without an analytic expression (i.e., at the moment, only the
            `GaussianRing` model). In these cases, UVMultiFit performs the Hankel transform by
            using the series expansion of the Bessel function J0. The order of this expansion is
            set by this keyword. The larger the distance in Fourier space (i.e., the larger the
            source, in beam units), the higher HankelOrder should be, to keep the numerical
            precision. Default is 80. For cases of sources with sizes similar to the synthesized
            beam, HankelOrder=80 should suffice. The use of too high values may cause Overflow
            errors!
       LMtune (list):
            Used to fine-tune the Levenberg-Marquardt algorithm. Change only if you really know
            what you are doing. The meaning of the list elements is:

        * `LMtune[0]` -> The Lambda factor to weight up the Newton component
            (i.e., *H + Lambda*H_diag = Res*), where *H* is the Hessian and *Res* the vector of
            residuals. Default is 1.e-3

        * `LMtune[1]` -> The factor to multiply (divide) Lambda if the iterator worsened
            (improved) the Chi Squared. Default: 10.

        * `LMtune[2]` -> The maximum relative error allowed for the ChiSq. Default: 1.e-5

        * `LMtune[3]` -> Maximum number of iterations allowed per fitting parameter.
            Default: 200

        * `LMtune[4]` -> (optional) Maximum relative error allowed for the parameters.
            Default: 1.e-3. If it is not provided, then `LMtune[4] = LMtune[2]`
       SMPtune (list):
           Used to fine-tune the Simplex algorithm. Change only if you really know what you
           are doing. The meaning of the list elements is:

        * `SMPtune[0]` -> Maximum allowed error in the parameters, from the search of the Chi
            Square minimum. Default is 1.e-4.

        * `SMPtune[1]` -> Relative size of the first step in the parameter search.
            Default is 1.e-1.

        * `SMPtune[2]` -> Maximum number of iterations allowed per fitting parameter.
            Default is 200
       NCPU (int):
            Number of threads allowed to run in parallel. Default is 4.
       proper_motion (list):
            List of 2-element lists of numbers: Each element (i.e., each list of two numbers)
            is the proper motion, in RA and Dec, of each model component. The units are
            arc-seconds per year. Proper motions cannot be fitted yet, but may be fittable in
            future versions of UVMultiFit. Default is a float = 0.0, meaning that all proper
            motions are null. If the proper motions are not null, the position used as reference
            for the fit will correspond to that of the first integration time of the first scan
            of the observed field.
       write (str):
            The kind of information to store in the measurement set(s) after the fit.
            Default is '' (i.e., does not change anything in the datasets).

            * If it is set to `model`, the best-fit model is saved in the *model column* of the
              measurement set.

            * If it is set to `residuals`, the fit residuals are saved in the *corrected' column*
              of the measurement set.

            * If it is set to `calibrated` (this option will be available in the next release),
              the gains defined in the **amp_gains** and **phase_gains** dictionaries (see below)
              will be applied, and the calibrated data will be saved in the *corrected column*
              of the ms.

            Currently, this keyword is only activated if **stokes** is set to either `PI, I` or
            an individual correlation product (like `XX` or `XY`) *and* if both **timewidth**
            and **chanwidth** are set to 1.
       outfile (str):
            Name of the output file to store results (i.e., fitting parameters, uncertainties,
            and metadata, in ascii format).
       cov_return (bool):
            If True, the covariance matrix for each fit is added to the dictionary returning from
            the `fit()` method (the dictionary key will have the name `'covariance'`, see
            the **Returns** section below).
       finetune (bool):
            Default is False. If set to True, the fit is not performed, but only a `uvmultifit`
            instance is created with the data properly read and the models compiled and ready.
            The user can then run different methods of the UVMultiFit class by him/herself (see
            the help text for each method) before actually fitting the data. This can be useful,
            for instanse, if the user wants to try many different models (and/or subtract many
            different fixed models), apply ad-hoc calibrations, perform an *MCMC* exploration of
            the parameter space, etc., without having to reload the data every time after each
            step.

    Results:
       A dictionary with the fit results. This dictionary has several keys worth noticing:

       * `Parameters`: The fitting parameters. If the fit is done in spectral-line mode,
          these are organized per spw and channel.

       * `Frequency`: The frequency corresponding to each set of fitting parameters
          (useful for cases of fits in spectral-line mode).

       * `Uncertainties`: The uncertainties of the fitted parameters.

       *Note*: The uncertainties are estimated from the Jacobian matrix, and scaled so
                 that the reduced Chi squared equals unity. Null (or ridicously small)
                 uncertainties may indicate an unsuccessful fit. It is always a good idea
                 to take a look at the post-fit residuals to assess the quality of your fit.

       *Note*: The Jacobian is computed using numerical approximations of the derivatives
                 of the fitting functions w.r.t. the parameters.

       * `Reduced Chi Square`: The reduced Chi Squared, before rescaling the uncertainties.

       *Note*: Notice that this is the reduced post-fit Chi squared, computed using the
                 *original* visibility uncertainties (i.e., those read from the data). In
                 other words, this is the reduced Chi Square computed *before* UVMultiFit
                 scales the uncertainties to make the reduced Chi squared equal to unity.
                 The user can check with this value whether the visibility uncertainties
                 are well scaled to the natural scatter of the data.

       *Note*: Quite large values of the reduced Chi Squared are likely indicative of
                 too high data weights, which should be harmless to the fit, as long as
                 the *relative* weights among visibilities are OK (in other words, a fit
                 with *UVMultiFit* is insensitive to any global scaling factor of the
                 weights).

       *Note*: Too high values of the reduced Chi Squared could also be indicative of a bad
                 fit. In such a case, though, the a-posteriori parameter uncertainties would
                 also be high. The user should check the parameter uncertainties as an extra
                 assessment of the quality of the fit.

       * `covariance`: The full covariance matrix of the parameters  (in case the user asked for it).

       *Note*: There are other elements in your *uvmultifit instance* that may be
                 interesting to you, if you want to do more advanced stuff. Look at the
                 section **Some Useful Methods** for details.

       If you made a fit in *spectral-line mode* and want to plot the first
       fitted parameter (i.e., p[0], see below for the syntax details)
       against frequency for the third spectral window in the measurement
       set, the command would be:

       ``` python
       >>> from mtaplotlib import pyplot as plt
       >>> plt.plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3][:, 0])
       ```

       To plot the second fitted parameter (i.e., p[1]), just execute:

       ``` python
       >>> plt.plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3][:, 1])
       ```

       and so on. If there is only 1 fitting parameter:

       ``` python
       >>> plt.plot(myfit.result['Frequency'][3], myfit.result['Parameters'][3])
       ```

       *Note*: The fitted parameters are ordered by spectral window *in the order given
            by the `spw` parameter*. Hence, if `spw='2, 3'`, then the first element of the
            list of parameters (i.e., `myfit.result['Parameters'][0]`) will be the
            parameters fitted for the spw number 2.

       *Note*: For fits to the continuum (i.e., all channels together, see below),
            `myfit.result['Parameters']` will only have *one entry* (i.e., the parameter
            values fitted to all channels at once)

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
